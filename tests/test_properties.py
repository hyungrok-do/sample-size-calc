"""Properties that every procedure must satisfy, checked over the whole catalog.

The reference-value tests (core/tests) compare selected calculations with
independent computations.  These tests instead sweep all procedures and all
solve-for modes with generic invariants, driven only by the procedures'
schemas:

* solving for N or for an effect and then computing the power at the solved
  value gives back the reported power, and N is the smallest one that works;
* precision procedures (no power) round-trip between N and the width;
* power never decreases when alpha grows;
* results do not depend on the number of threads;
* a seeded sweep of perturbed inputs never crashes, hangs, produces non-finite
  numbers or power outside [0, 1], and every error is a clear user message.
"""
from __future__ import annotations

import itertools
import math
import os
import random
import re
import time
from typing import Dict, List, Optional

import pytest

import ssc

CATALOG = ssc.catalog()["procedures"]
IDS = [p["id"] for p in CATALOG]
SCHEMAS = {pid: ssc.schema(pid) for pid in IDS}

# messages that betray an engine bug rather than a problem with the inputs
INTERNAL = ("missing input", "missing selection", "missing list input", "not bracketed", "internal",
            "std::", "nan", "unknown option", "out of range")


# ----------------------------------------------------------------------------- schema helpers
def visible(rule: list, lookup) -> bool:
    if not rule:
        return True
    for conj in rule:
        if all((lookup(c["key"]) in c.get("values", [])) != bool(c.get("negate", False)) for c in conj):
            return True
    return False


def params(pid: str) -> Dict[str, dict]:
    return {p["id"]: p for p in SCHEMAS[pid]["params"]}


def modes(pid: str) -> List[str]:
    return [o["value"] for o in SCHEMAS[pid]["solve_for"]]


def first_value(text: str) -> str:
    """'0.8 0.9' -> '0.8'; '10 to 100 by 10' -> '10' (vectors are kept whole)."""
    text = str(text).strip()
    return text.split()[0] if text else text


def base_inputs(pid: str) -> Dict[str, str]:
    """Defaults updated with the example, every list reduced to its first value."""
    schema = SCHEMAS[pid]
    inputs = {p["id"]: str(p["default"]) for p in schema["params"]}
    inputs.update({k: str(v) for k, v in (schema.get("example") or {}).get("inputs", {}).items()})
    ps = params(pid)
    for k, v in list(inputs.items()):
        if ps.get(k, {}).get("type") == "number":
            inputs[k] = first_value(v)
    return inputs


def is_visible(pid: str, key: str, inputs: Dict[str, str], solve: str) -> bool:
    p = params(pid)[key]
    return visible(p.get("visible", []), lambda k: solve if k == "$solve" else inputs.get(k, ""))


def run(pid: str, solve: str, inputs: Dict[str, str], threads: int = 0) -> dict:
    return ssc.run(pid, solve, inputs, threads=threads, power_curve=False)


def first_row(res: dict) -> Optional[dict]:
    if not res.get("ok") or not res.get("rows"):
        return None
    row = res["rows"][0]
    return None if "_error" in row else row


def fmt(v) -> str:
    return repr(float(v)) if not float(v).is_integer() else str(int(v))


POWER_KEYS = ("power", "pa", "pb", "pab")  # the factorial ANOVA reports one power per effect


def powers(row: dict) -> Dict[str, float]:
    return {k: row[k] for k in POWER_KEYS if isinstance(row.get(k), (int, float))}


def make_visible(pid: str, key: str, inputs: Dict[str, str], solve: str) -> Dict[str, str]:
    """Change one choice so that `key` becomes an input (e.g. effect size given as f)."""
    if is_visible(pid, key, inputs, solve):
        return inputs
    for other, p in params(pid).items():
        if p["type"] != "select":
            continue
        for o in p["options"]:
            trial = {**inputs, other: o["value"]}
            if is_visible(pid, key, trial, solve):
                return trial
    return inputs


def mapped_inputs(pid: str, row: dict, inputs: Dict[str, str], src: str, dst: str):
    """Inputs for solving `dst`, with the numbers solved under `src` copied from the result row.

    Returns (inputs, mapped parameter ids)."""
    out = dict(inputs)
    if src in params(pid) and src in row:
        out = make_visible(pid, src, out, dst)
    mapped = []
    for key, p in params(pid).items():
        if p["type"] != "number" or not isinstance(row.get(key), (int, float)):
            continue
        if is_visible(pid, key, out, dst) and not is_visible(pid, key, inputs, src):
            out[key] = fmt(row[key])
            mapped.append(key)
    return out, mapped


# ----------------------------------------------------------------------------- round trips
ROUND_TRIPS = [(pid, m) for pid in IDS if "power" in modes(pid) for m in modes(pid) if m != "power"]


@pytest.mark.parametrize("pid,mode", ROUND_TRIPS, ids=[f"{p}:{m}" for p, m in ROUND_TRIPS])
def test_solved_value_gives_back_the_power(pid, mode):
    inputs = base_inputs(pid)
    row = first_row(run(pid, mode, inputs))
    assert row is not None, f"{pid} could not solve {mode} for its example"
    back, mapped = mapped_inputs(pid, row, inputs, mode, "power")
    assert mapped, f"nothing in the {mode} result maps to the power inputs"
    check = first_row(run(pid, "power", back))
    assert check is not None, f"power at the solved {mapped} failed"
    got, want = powers(check), powers(row)
    assert got and set(got) == set(want)
    for key in want:
        assert got[key] == pytest.approx(want[key], rel=1e-7, abs=1e-9), (key, mapped, back)
    target = float(inputs["power"])
    # the powers that must reach the target (one per effect for the factorial ANOVA's "all" goal)
    goal = {"A": ["pa"], "B": ["pb"], "AB": ["pab"]}.get(str(row.get("goal")), list(want))
    if "target_power" in row:  # solved for a sample size: the power reaches the target ...
        assert min(want[k] for k in goal) >= target - 1e-12
        # ... and one fewer subject (per group) does not
        if len(mapped) == 1:
            key = mapped[0]
            p = params(pid)[key]
            for step in (1, 2):
                smaller = dict(back)
                value = float(back[key]) - step
                lo = p.get("min")
                if lo is not None and (value < lo or (p.get("min_open") and value <= lo)):
                    break
                smaller[key] = fmt(value)
                below = first_row(run(pid, "power", smaller))
                if below is not None:
                    lower = powers(below)
                    assert min(lower[k] for k in goal) < target, f"{key} = {value} already reaches the target"
                    break
    else:  # solved for an effect: the power equals the target
        assert want["power"] == pytest.approx(target, abs=2e-6)


PRECISION = [pid for pid in IDS if "power" not in modes(pid) and "n" in modes(pid) and len(modes(pid)) == 2]


@pytest.mark.parametrize("pid", PRECISION)
def test_precision_procedures_round_trip(pid):
    other = [m for m in modes(pid) if m != "n"][0]
    inputs = base_inputs(pid)
    row = first_row(run(pid, other, inputs))
    assert row is not None
    back, mapped = mapped_inputs(pid, row, inputs, other, "n")
    assert mapped
    again = first_row(run(pid, "n", back))
    assert again is not None
    # the original sample size achieves the width, so the solved one is no larger
    sizes = [k for k in params(pid) if is_visible(pid, k, inputs, other) and not is_visible(pid, k, inputs, "n")]
    for key in sizes:
        if key in again:
            assert again[key] <= float(inputs[key]) + 1e-9


ALPHA = [pid for pid in IDS if "power" in modes(pid) and "alpha" in params(pid)]


@pytest.mark.parametrize("pid", ALPHA)
def test_power_does_not_decrease_with_alpha(pid):
    inputs = base_inputs(pid)
    seen = []
    for alpha in ("0.01", "0.05", "0.1"):
        row = first_row(run(pid, "power", {**inputs, "alpha": alpha}))
        assert row is not None, alpha
        seen.append(powers(row))
    for key in seen[0]:
        values = [p[key] for p in seen]
        assert values == sorted(values), (key, values)


@pytest.mark.parametrize("pid", IDS)
def test_results_do_not_depend_on_threads(pid):
    ex = SCHEMAS[pid]["example"]
    one = ssc.run(pid, ex["solve_for"], ex["inputs"], threads=1)
    many = ssc.run(pid, ex["solve_for"], ex["inputs"], threads=4)
    assert one["rows"] == many["rows"]
    assert one.get("plots") == many.get("plots")


# ----------------------------------------------------------------------------- fuzz sweep
# cases per procedure and seed: SSC_FUZZ_CASES=300 SSC_FUZZ_SEED=7 runs a longer sweep
FUZZ_CASES = int(os.environ.get("SSC_FUZZ_CASES", "40"))
FUZZ_SEED = int(os.environ.get("SSC_FUZZ_SEED", "2026"))
FACTORS = (0.1, 0.5, 0.9, 1.1, 2.0, 5.0)


def perturb(p: dict, value: str, rng: random.Random) -> str:
    """A nearby or boundary value inside the parameter's declared range."""
    lo, hi = p.get("min"), p.get("max")
    try:
        v = float(value)
    except ValueError:
        return value
    if lo is not None and hi is not None and hi - lo <= 1.0001:  # probabilities, proportions, correlations
        eps = (hi - lo) * 0.02
        choice = rng.choice(["near-lo", "near-hi", "uniform", "uniform"])
        if choice == "near-lo":
            new = lo + eps * rng.uniform(0.05, 1)
        elif choice == "near-hi":
            new = hi - eps * rng.uniform(0.05, 1)
        else:
            new = rng.uniform(lo + eps, hi - eps)
    else:
        new = v * rng.choice(FACTORS) if v != 0 else rng.choice((-2.0, -0.5, 0.5, 2.0))
        if lo is not None:
            new = max(new, lo + (1 if p.get("integer") else 1e-6) * (1 if p.get("min_open") else 0))
        if hi is not None:
            new = min(new, hi)
    if p.get("integer"):
        new = max(round(new), math.ceil(lo) if lo is not None else round(new))
        if hi is not None:
            new = min(new, hi)
        return str(int(new))
    return f"{new:.6g}"


def fuzz_cases(pid: str, n: int, seed: int) -> List[tuple]:
    rng = random.Random(f"{pid}:{seed}")
    cases = []
    for _ in range(n):
        solve = rng.choice(modes(pid))
        inputs = base_inputs(pid)
        ps = params(pid)
        for key, p in ps.items():  # random choices, respecting the dependencies between them
            if p["type"] == "select" and rng.random() < 0.5:
                options = [o["value"] for o in p["options"]
                           if visible(o.get("visible", []), lambda k: solve if k == "$solve" else inputs.get(k, ""))]
                if options:
                    inputs[key] = rng.choice(options)
            elif p["type"] == "bool" and rng.random() < 0.3:
                inputs[key] = rng.choice(("true", "false"))
        numeric = [k for k, p in ps.items() if p["type"] == "number" and is_visible(pid, k, inputs, solve)]
        for key in rng.sample(numeric, k=min(len(numeric), rng.randint(1, 3))):
            inputs[key] = perturb(ps[key], inputs[key], rng)
        cases.append((solve, inputs))
    return cases


def check_result(pid: str, solve: str, inputs: dict, res: dict, seconds: float) -> List[str]:
    problems = []
    where = f"{pid} [{solve}] {inputs}"
    # simulations refit a model to every simulated study, so extreme designs (thousands of clusters)
    # legitimately take longer than closed-form calculations: a win-ratio search that ends at 4,117 patients per
    # group compares every pair in each of 500 studies, 18 s on 4 cores and 71 s on one.  The limit is for hangs.
    if seconds > (120 if SCHEMAS[pid].get("simulation") else 10):
        problems.append(f"slow ({seconds:.1f} s): {where}")
    if not res.get("ok"):
        msg = (res.get("error") or "") + " " + " ".join((res.get("input_errors") or {}).values())
        if any(s in msg.lower() for s in INTERNAL):
            problems.append(f"internal error '{msg.strip()}': {where}")
        return problems
    for row in res.get("rows", []):
        if "_error" in row:
            # a numerical routine's own words stay behind the generic message (in _detail)
            text = row["_error"] + " " + row.get("_detail", "")
            if any(s in text.lower() for s in INTERNAL):
                problems.append(f"internal error '{text.strip()}': {where}")
            continue
        for key, value in row.items():
            if isinstance(value, float) and not math.isfinite(value) and not math.isnan(value):
                problems.append(f"{key} = {value}: {where}")
            if isinstance(value, float) and math.isnan(value) and key in ("power", "n", "n1", "n2"):
                problems.append(f"{key} is NaN: {where}")
        p = row.get("power")
        if isinstance(p, float) and not (-1e-12 <= p <= 1 + 1e-12):
            problems.append(f"power {p} outside [0, 1]: {where}")
    return problems


@pytest.mark.parametrize("pid", IDS)
def test_fuzzed_inputs_are_handled(pid):
    problems = []
    for solve, inputs in fuzz_cases(pid, FUZZ_CASES, seed=FUZZ_SEED):
        start = time.perf_counter()
        res = run(pid, solve, inputs)
        problems += check_result(pid, solve, inputs, res, time.perf_counter() - start)
    assert not problems, "\n".join(problems[:10])


EFFECT_SOLVES = [(pid, m) for pid, m in ROUND_TRIPS if "power" in modes(pid) and m not in ("n", "k", "g", "m")]


@pytest.mark.parametrize("pid,mode", EFFECT_SOLVES, ids=[f"{p}:{m}" for p, m in EFFECT_SOLVES])
def test_unreachable_targets_give_clear_messages(pid, mode):
    """Targets at or below alpha, or effects too small to detect, fail with a plain message."""
    inputs = base_inputs(pid)
    for power in ("0.001", "0.999999"):
        start = time.perf_counter()
        res = run(pid, mode, {**inputs, "power": power})
        assert not check_result(pid, mode, inputs, res, time.perf_counter() - start)


EXTREMES = [
    {"alpha": "0.000001"},
    {"power": "0.9999"},
    {"alpha": "0.5", "power": "0.51"},
    {"dropout": "0.9"},
]


@pytest.mark.parametrize("pid", IDS)
def test_extreme_but_valid_inputs(pid):
    """Very small alpha, very high power or heavy dropout: a result or a plain message, quickly."""
    problems = []
    for extra in EXTREMES:
        for solve in modes(pid):
            inputs = {**base_inputs(pid), **{k: v for k, v in extra.items() if k in params(pid)}}
            start = time.perf_counter()
            res = run(pid, solve, inputs)
            problems += check_result(pid, solve, inputs, res, time.perf_counter() - start)
    assert not problems, "\n".join(problems[:10])


# Inputs from the long fuzzing sweep that once took from several seconds to
# minutes: tiny effects that need very large designs, an exact search far
# beyond what enumeration can afford, a slope solve that runs into the far tail
# of the noncentral t, and the exact correlation distribution next to ρ0 = −1.
SLOW_BEFORE = [
    ("prop.one.simon", "n", {"criterion": "optimal", "power": "0.8", "alpha": "0.05", "p0": "0.267458", "p1": "0.3"},
     {"n": 1410, "n1": 506, "r1": 140}),
    ("prop.one.simon", "n", {"criterion": "optimal", "power": "0.270909", "alpha": "0.05", "p0": "0.985303",
                             "p1": "0.989282"}, {"n": 1362, "n1": 181, "r1": 179}),
    ("prop.two.equiv", "n", {"method": "exact", "power": "0.8", "alpha": "0.05", "alloc": "equal", "dl": "-0.00641354",
                             "du": "0.995227", "es": "props", "p1": "0.6", "p2": "0.6"}, {}),
    ("prop.two", "n", {"method": "exact", "test": "fisher", "power": "0.8", "alpha": "0.05", "alloc": "equal",
                       "es": "props", "p1": "0.51", "p2": "0.5"}, {}),
    ("reg.slope", "b1", {"alternative": "two", "power": "0.8", "alpha": "0.05", "n": "250", "b0": "0", "sx": "5",
                         "sdtype": "sy", "sigma": "0.1"}, {}),
    ("diag.sens", "n", {"goal": "both", "power": "0.8", "alpha": "0.45", "prev": "0.2", "se0": "0.7",
                        "se1": "0.85", "sp0": "0.8", "sp1": "0.800195"}, {"n": 4919953}),
    ("corr.one", "n", {"alternative": "two", "method": "exact", "power": "0.8", "alpha": "0.05", "rho0": "-0.999999",
                       "rho1": "0.18"}, {"n": 4}),
    ("prop.paired.mcnemar", "n", {"alternative": "two", "method": "exact", "power": "0.983872", "alpha": "0.019839",
                                  "es": "props", "p10": "0.162728", "p01": "0.1"}, {"n": 1345}),
]


@pytest.mark.parametrize("pid,solve,inputs,expect", SLOW_BEFORE, ids=[f"{c[0]}:{c[1]}" for c in SLOW_BEFORE])
def test_formerly_slow_inputs_finish_quickly(pid, solve, inputs, expect):
    start = time.perf_counter()
    res = run(pid, solve, {**base_inputs(pid), **inputs})
    seconds = time.perf_counter() - start
    row = first_row(res)
    assert row is not None and not row.get("_error"), res
    assert seconds < 5, f"{pid} took {seconds:.1f} s"
    assert 0 <= row["power"] <= 1
    for key, value in expect.items():
        assert row[key] == value, (key, row[key], value)
    if pid in ("prop.two.equiv", "prop.two"):  # too large to search exactly: the normal approximation, said so
        assert any("normal approximation" in n for n in row.get("_notes", [])), row.get("_notes")


# ----------------------------------------------------------------------------- relations between procedures
def power_of(pid: str, inputs: dict) -> Optional[float]:
    row = first_row(run(pid, "power", {k: str(v) for k, v in inputs.items()}, threads=1))
    return None if row is None else row["power"]


def test_procedures_agree_where_they_coincide():
    """Special cases in which two procedures compute the same power: a paired t-test is a one-sample t-test on
    the differences, Welch's test with equal SDs and equal groups is the pooled t-test, a one-way ANOVA of two
    groups is the two-sided t-test (F = t²), a non-inferiority test with a vanishing margin is a one-sided
    superiority test, a cluster design with one subject per cluster and no ICC is the t-test, and a
    group-sequential design with one look is the fixed design."""
    rng = random.Random(2026)
    for _ in range(25):
        alt = rng.choice(["two", "upper", "lower"])
        test = rng.choice(["t", "z"])
        alpha = rng.choice([0.01, 0.025, 0.05, 0.1])
        n = rng.randint(4, 120)
        sign = -1 if alt == "lower" else 1
        d, sd = sign * rng.uniform(0.1, 1.2), rng.uniform(0.3, 3)
        t = power_of("means.two.t", dict(alternative=alt, test=test, alpha=alpha, n1=n, mu1=d, mu2=0, sigma=sd))
        pairs = [
            (power_of("means.paired.t", dict(alternative=alt, test=test, alpha=alpha, n=n, mud=d, sdtype="sdd", sdd=sd)),
             power_of("means.one.t", dict(alternative=alt, test=test, alpha=alpha, n=n, mu0=0, mu1=d, sigma=sd)), 1e-12),
            (power_of("means.two.welch", dict(alternative=alt, test=test, alpha=alpha, n1=n, mu1=d, mu2=0, sd1=sd, sd2=sd)),
             t, 1e-9),
            (power_of("means.two.cluster", dict(alternative=alt, test=test, alpha=alpha, k1=n, kr=1, m=1, icc=0, mu1=d,
                                                mu2=0, sigma=sd)), t, 1e-9),
            (power_of("means.two.t", dict(alternative={"two": "two", "upper": "lower", "lower": "upper"}[alt], test=test,
                                          alpha=alpha, n1=n, mu1=0, mu2=d, sigma=sd)), t, 1e-12),
            (power_of("means.anova", dict(alpha=alpha, n=n, estype="means", means=f"0 {abs(d)}", sigma=sd, alloc="equal")),
             power_of("means.two.t", dict(alternative="two", test="t", alpha=alpha, n1=n, mu1=abs(d), mu2=0, sigma=sd)),
             1e-9),
            (power_of("means.two.ni", dict(testtype="ni", direction="higher", test=test, alpha=alpha, n1=n, margin=1e-9,
                                           delta=abs(d), sigma=sd)),
             power_of("means.two.t", dict(alternative="upper", test=test, alpha=alpha, n1=n, mu1=abs(d), mu2=0,
                                          sigma=sd)), 1e-6),
        ]
        if alt != "lower":
            pairs.append((power_of("gs.means", dict(alternative=alt, looks=1, tfrac="1", alpha=alpha, n1=n, mu1=d, mu2=0,
                                                    sigma=sd)),
                          power_of("means.two.t", dict(alternative=alt, test="z", alpha=alpha, n1=n, mu1=d, mu2=0,
                                                       sigma=sd)), 1e-7))
        for i, (a, b, tol) in enumerate(pairs):
            assert a is not None and b is not None and a == pytest.approx(b, abs=tol), (i, alt, test, alpha, n, d, sd)
    for _ in range(10):
        hr, acc, fol, med = rng.uniform(0.4, 0.9), rng.uniform(3, 24), rng.uniform(3, 24), rng.uniform(4, 30)
        n = rng.randint(20, 400)
        common = dict(alpha=0.05, n1=n, hr=hr, accrual=acc, follow=fol)
        fixed = power_of("surv.logrank", dict(alternative="two", method="schoenfeld", survin="median", med2=med, eta=0,
                                              **common))
        assert power_of("gs.logrank", dict(alternative="two", looks=1, tfrac="1", med2=med, **common)) == \
            pytest.approx(fixed, abs=1e-6)


ONE_SIDED = [pid for pid in IDS if not SCHEMAS[pid].get("simulation") and "alternative" in params(pid)
             and {"upper", "lower"} <= {o["value"] for o in params(pid)["alternative"]["options"]}]


@pytest.mark.parametrize("pid", ONE_SIDED)
def test_an_effect_on_the_wrong_side_is_named(pid):
    """Solving for a sample size with the effect on the side opposite to a one-sided alternative says so, rather
    than searching in vain and reporting that no sample size reaches the target."""
    for solve in [m for m in modes(pid) if m in ("n", "k", "g")]:
        messages = []
        for side in ("upper", "lower"):
            res = run(pid, solve, {**base_inputs(pid), "alternative": side}, threads=1)
            row = first_row(res)
            messages.append(None if row is not None else res.get("error") or res["rows"][0]["_error"])
        assert messages.count(None) == 1, messages      # one side works, the other does not
        wrong = next(m for m in messages if m)
        assert "opposite" in wrong and "cannot be reached" not in wrong, wrong


DROPOUT = [pid for pid in IDS if not SCHEMAS[pid].get("simulation") and "dropout" in params(pid) and "n" in modes(pid)]


@pytest.mark.parametrize("pid", DROPOUT)
def test_dropout_inflates_each_group_and_sums_them(pid):
    for d in (0.05, 0.13, 0.3):
        res = run(pid, "n", {**base_inputs(pid), "dropout": str(d)}, threads=1)
        row = first_row(res)
        inflated = next((x for x in (res.get("dropout") or {}).get("rows", []) if x.get("_row") == 0), None)
        assert row is not None and inflated, (pid, d)
        groups = [k for k in ("n1", "n2") if inflated.get(f"{k}_inflated") is not None]
        sizes = groups or [k for k in ("n", "N") if inflated.get(f"{k}_inflated") is not None]  # one size, or the total
        assert sizes, inflated
        even = 2 if pid in ("means.crossover", "means.be.ratio") else 1  # two sequences of equal size
        for key in sizes:
            assert inflated[f"{key}_inflated"] == math.ceil(row[key] / (1 - d) / even - 1e-9) * even, (key, d)
        if groups and inflated.get("n_inflated") is not None:   # the total is the sum of the enrolled groups
            assert inflated["n_inflated"] == sum(inflated[f"{k}_inflated"] for k in groups), (d, inflated)


def test_a_comma_inside_a_number_is_not_read_as_a_list():
    """'2,500' could be a thousands separator, a decimal comma or a list of 2 and 500: rather than quietly make
    two scenarios, the engine asks for the number written plainly; a comma followed by a space separates values."""
    for text in ("2,500", "1,5", "1,000,000"):
        res = run("means.two.t", "power", {**base_inputs("means.two.t"), "n1": text})
        assert not res["ok"] and "without commas" in res["input_errors"]["n1"] and text in res["input_errors"]["n1"]
    res = run("means.two.t", "power", {**base_inputs("means.two.t"), "n1": "20, 40"})
    assert [r["n1"] for r in res["rows"]] == [20, 40]


# ----------------------------------------------------------------------------- found by the correctness sweep
def test_ratio_allocation_finds_the_smallest_n():
    """With N2 = ⌈r·N1⌉ and r < 1, the smallest N1 whose N2 is large enough lies below ⌈min/r⌉ (rounding up helps):
    the search starts there, so no smaller design reaches the target."""
    rng = random.Random(7)
    cases = [("means.two.t", {"mu1": "3", "mu2": "0", "sigma": "1"}), ("var.two", {"alternative": "lower", "ratio": "0.1"}),
             ("means.two.equiv", {"test": "z", "e": "2", "delta": "0"})]
    for pid, extra in cases:
        for _ in range(6):
            r = rng.choice(["0.1", "0.25", "0.3", "0.6"])
            inputs = {**base_inputs(pid), **extra, "alloc": "ratio", "r": r,
                      "power": rng.choice(["0.6", "0.7", "0.8"]), "alpha": rng.choice(["0.1", "0.2"])}
            row = first_row(run(pid, "n", inputs, threads=1))
            assert row is not None, (pid, inputs)
            for n1 in range(2, int(row["n1"])):   # every smaller allowed design falls short
                below = first_row(run(pid, "power", {**inputs, "n1": str(n1)}, threads=1))
                if below is not None:
                    assert below["power"] < float(inputs["power"]), (pid, inputs, n1, below["power"], row["n1"])


def test_fishers_exact_power_counts_the_most_extreme_tables():
    """With a large effect and many subjects, nearly all the power lies in tables far from the null mean."""
    for n1, p1, p2 in ((1000, "0.99", "0.01"), (2000, "0.9", "0.1"), (5000, "0.99", "0.01")):
        row = first_row(run("prop.two", "power", {**base_inputs("prop.two"), "test": "fisher", "method": "exact",
                                                  "n1": str(n1), "es": "props", "p1": p1, "p2": p2}, threads=1))
        assert row["power"] == pytest.approx(1, abs=1e-9), (n1, p1, p2)
    row = first_row(run("prop.two", "p1", {**base_inputs("prop.two"), "test": "fisher", "method": "exact", "n1": "1000",
                                           "es": "props", "p2": "0.05", "power": "0.8"}, threads=1))
    assert row is not None and row["power"] == pytest.approx(0.8, abs=1e-6)


@pytest.mark.parametrize("pid,inputs,first", [
    ("prop.one", {"test": "exact", "alternative": "two", "alpha": "0.05", "power": "0.8", "p0": "0.0085", "p1": "0.0425"}, 129),
    ("prop.one", {"test": "exact", "alternative": "upper", "alpha": "0.1", "power": "0.8", "p0": "0.00252", "p1": "0.0101"},
     423),
    ("rates.one.poisson", {"test": "exact", "alternative": "lower", "alpha": "0.1", "power": "0.8", "time": "1",
                           "lam0": "0.005", "lam1": "0.001"}, 778),
    ("rates.one.poisson", {"test": "z", "method": "exact", "alternative": "two", "alpha": "0.1", "power": "0.5",
                           "lam0": "0.01", "lam1": "0.002"}, 271),
])
def test_saw_toothed_power_with_rare_events_gives_the_first_n(pid, inputs, first):
    """Rare events make long saw-teeth (about 1/p0 subjects): the first sample size reaching the target can lie more
    than a few dozen below the one bisection finds (values from R's exact binomial and Poisson tests)."""
    row = first_row(run(pid, "n", {**base_inputs(pid), **inputs}, threads=1))
    assert row["n"] == first and row["power"] >= float(inputs["power"])


@pytest.mark.parametrize("pid,mode", EFFECT_SOLVES, ids=[f"{p}:{m}" for p, m in EFFECT_SOLVES])
def test_an_effect_solve_needs_a_target_above_alpha(pid, mode):
    """At no effect the power is the significance level: a target at or below it has no answer (not 'no effect')."""
    res = run(pid, mode, {**base_inputs(pid), "power": "0.001"}, threads=1)
    assert first_row(res) is None, res["rows"][0] if res.get("rows") else res


@pytest.mark.parametrize("pid", ["means.one.equiv", "means.two.equiv", "means.be.ratio", "prop.two.equiv"])
def test_equivalence_tests_need_alpha_below_one_half(pid):
    """Two one-sided tests at α ≥ 0.5 have overlapping regions that reject every outcome."""
    res = run(pid, "power", {**base_inputs(pid), "alpha": "0.6"})
    assert not res["ok"] and "less than 0.5" in res["input_errors"]["alpha"]


def test_impossible_null_hypotheses_are_refused():
    res = run("means.rm.anova", "power", {**base_inputs("means.rm.anova"), "means": "0 0.2 0.4", "rho": "-0.9", "n": "20"})
    assert first_row(res) is None and "−1/(m − 1)" in res["rows"][0]["_error"]
    for extra in ({"scale": "diff", "d0": "-0.7", "es": "props", "p1": "0.5", "p2": "0.5"},
                  {"scale": "ratio", "r0": "1.1", "h1": "upper", "es": "props", "p1": "0.916", "p2": "0.927"}):
        res = run("prop.two.ni", "power", {**base_inputs("prop.two.ni"), **extra, "method": "exact", "n1": "50"})
        assert first_row(res) is None and "boundary of the null hypothesis" in res["rows"][0]["_error"], extra
    res = run("means.one.ci", "n", {**base_inputs("means.one.ci"), "e": "1e-6"})
    assert "target half-width cannot be reached" in res["rows"][0]["_error"]


LEVEL_CASES = [pid for pid in IDS if "n" in modes(pid) and {"power", "alpha"} <= set(params(pid))]


@pytest.mark.parametrize("pid", LEVEL_CASES)
def test_a_target_power_at_or_below_alpha_is_refused(pid):
    """With no effect at all a test rejects with probability α, so a target power at or below α needs no study:
    the n-solve says so rather than returning the smallest design."""
    inputs = base_inputs(pid)
    if not (is_visible(pid, "power", inputs, "n") and is_visible(pid, "alpha", inputs, "n")):
        pytest.skip("power or alpha is not an input of this design")
    alpha = float(inputs["alpha"])
    res = run(pid, "n", {**inputs, "power": fmt(alpha / 2)}, threads=1)
    assert first_row(res) is None, (pid, res["rows"][0] if res.get("rows") else res)


@pytest.mark.parametrize("pid", [pid for pid in IDS if "alternative" in params(pid) and "alpha" in params(pid)])
def test_a_one_sided_test_needs_alpha_below_one_half(pid):
    """A one-sided test at α ≥ 0.5 rejects at least as often as not when there is no effect."""
    sides = [o["value"] for o in params(pid)["alternative"]["options"] if o["value"] != "two"]
    inputs = {**base_inputs(pid), "alternative": sides[0], "alpha": "0.6"}
    solve = "power" if "power" in modes(pid) else modes(pid)[0]
    if not (is_visible(pid, "alpha", inputs, solve) and is_visible(pid, "alternative", inputs, solve)):
        pytest.skip("alpha or the alternative is not an input of this design")
    res = run(pid, solve, inputs, threads=1)
    errors = [r.get("_error", "") for r in res.get("rows", [])] + list((res.get("input_errors") or {}).values())
    assert not res["ok"] or first_row(res) is None, pid
    assert any("0.5" in e for e in errors), (pid, errors)


@pytest.mark.parametrize("inputs,first", [
    ({"goal": "se", "prev": "0.224", "se0": "0.9", "se1": "0.995", "alpha": "0.025", "power": "0.8"}, 161),
    ({"goal": "sp", "prev": "0.583", "sp0": "0.86", "sp1": "0.99", "alpha": "0.05", "power": "0.8"}, 48),
])
def test_exact_sensitivity_power_gives_the_first_n(inputs, first):
    """Power is 0 until all-correct becomes significant, jumps above the target there, then decays below it for
    two teeth before the next jump: the first N must not be skipped (found by brute force over every N)."""
    row = first_row(run("diag.sens", "n", {**base_inputs("diag.sens"), **inputs}, threads=1))
    assert row["n"] == first and row["power"] >= 0.8
    res = run("diag.sens", "power", {**base_inputs("diag.sens"), **inputs, "n": f"2 to {first - 1} by 1"})
    assert res["ok"] and len(res["rows"]) == first - 2
    assert all(r["power"] < 0.8 for r in res["rows"])


def test_designs_that_cannot_estimate_the_effect_are_refused():
    # one step: every cluster is under control in period 1 and under intervention in period 2
    for pid in ("cluster.sw", "sim.sw"):
        res = run(pid, "power", {**base_inputs(pid), "steps": "1", "reps": "100"})
        assert not res["ok"] and "at least 2" in res["input_errors"]["steps"]
    res = run("cluster.sw", "power", {**base_inputs("cluster.sw"), "m": "0.5"})
    assert not res["ok"] and "whole" in res["input_errors"]["m"]
    # the unstructured MMRM needs T + 2 patients at the last visit
    res = run("sim.mmrm", "power", {**base_inputs("sim.mmrm"), "n1": "3", "reps": "100"})
    assert first_row(res) is None and "at the last visit" in res["rows"][0]["_error"]
    row = first_row(run("sim.mmrm", "n", {**base_inputs("sim.mmrm"), "visits": "12", "drop1": "0.5", "drop2": "0.5",
                                          "delta": "3", "reps": "100"}, threads=1))
    assert row["n1"] * 0.5 + row["n2"] * 0.5 >= 14
    res = run("sim.ordinal", "power", {**base_inputs("sim.ordinal"), "p2": "0 0 1 0 0", "reps": "100"})
    assert "same category" in res["rows"][0]["_error"]
    res = run("sim.winratio", "n", {**base_inputs("sim.winratio"), "hr_d": "1", "hr_h": "1", "tier3": "none", "reps": "100"})
    assert "equally likely" in res["rows"][0]["_error"]


def test_a_robust_simulation_does_not_depend_on_the_units():
    """The t, Welch and Mann–Whitney tests are unchanged by a common change of location and scale, so SDs near the
    largest allowed (σ2 = 1e30) give the power of SDs of 1 (SDs of 1e300 once overflowed to power 0)."""
    def row(scale):
        return first_row(run("sim.means.robust", "power", {**base_inputs("sim.means.robust"), "sd1": f"{scale:g}",
                                                           "sd2": f"{2 * scale:g}", "mu1": f"{scale:g}", "mu2": "0",
                                                           "n1": "40", "reps": "200"}, threads=1))
    one, huge = row(1), row(5e29)
    assert huge["power"] == one["power"] and huge["approx_power"] == pytest.approx(one["approx_power"], rel=1e-12)


# ----------------------------------------------------------------------------- found by the second round
CURVE_CASES = [(pid, mode) for pid in IDS for mode in modes(pid)]


@pytest.mark.parametrize("pid,mode", CURVE_CASES, ids=[f"{p}:{m}" for p, m in CURVE_CASES])
def test_every_power_curve_has_points_the_engine_accepts(pid, mode):
    """A curve drawn for a result has points (a choice hidden by the solved effect, or a row with several powers,
    once left every point blank), and starts no lower than its input allows."""
    res = ssc.run(pid, mode, base_inputs(pid), threads=1, power_curve=True, approximate=True)
    if first_row(res) is None:
        pytest.skip("no result to draw a curve for")
    for plot in [p for p in res["plots"] if p["kind"] == "curve"]:
        assert any(len(s["x"]) > 0 for s in plot["series"]), (pid, mode, plot["title"])
        spec = params(pid).get(plot["x_key"], {})
        lowest = min(x for s in plot["series"] for x in s["x"])
        if spec.get("min") is not None:
            assert lowest >= spec["min"] and not (spec.get("min_open") and lowest == spec["min"]), (pid, lowest)


def test_a_curve_names_the_input_it_varies():
    """With N1 given and N2 solved the curve runs along N2 (and a click must choose separate sizes); with a
    percentage allocation it runs along the total."""
    for inputs, key, with_ in (({"alloc": "n1", "n1": "66", "mu1": "0.5"}, "n2", {"alloc": "sep"}),
                               ({"alloc": "n2", "n2": "200", "mu1": "0.5"}, "n1", {"alloc": "sep"}),
                               ({"alloc": "pct", "pct1": "30", "mu1": "0.5"}, "ntotal", {}),
                               ({"mu1": "0.5"}, "n1", {})):
        res = ssc.run("means.two.t", "n", {**base_inputs("means.two.t"), **inputs}, power_curve=True)
        curve = next(p for p in res["plots"] if p["kind"] == "curve")
        assert curve["x_key"] == key and curve["x_with"] == with_, inputs


def test_an_exact_curve_takes_one_method_judged_at_its_largest_point():
    """With N2 given apart from N1 the curve varies N1 alone.  Each point was judged by scaling both groups by the
    curve's largest N1 over its own, which sent the small sizes (and their huge scaled N2) to the normal
    approximation while the large ones stayed exact; each point now matches its own exact power."""
    base = {**base_inputs("prop.two"), "method": "exact", "test": "zpool", "alloc": "sep", "n1": "200",
            "n2": "400000", "es": "props", "p1": "0.5", "p2": "0.45"}
    res = ssc.run("prop.two", "power", base, threads=1, power_curve=True)
    series = next(p for p in res["plots"] if p["kind"] == "curve")["series"][0]
    for x, y in list(zip(series["x"], series["y"]))[::10]:
        row = ssc.run("prop.two", "power", {**base, "n1": str(int(x))}, threads=1, power_curve=False)["rows"][0]
        assert y == pytest.approx(row["power"], abs=1e-12), x


ANALYTIC_CURVE_CASES = [(pid, mode) for pid in IDS if not ssc.schema(pid).get("simulation") for mode in modes(pid)]


@pytest.mark.parametrize("pid,mode", ANALYTIC_CURVE_CASES, ids=[f"{p}:{m}" for p, m in ANALYTIC_CURVE_CASES])
def test_a_curve_passes_through_its_answer(pid, mode):
    """A curve computed as its answer is gives the answer's power at the answer's size.  (Each answer was written
    into its curve, which hid a curve drawn for the wrong effect and spiked one drawn another way.)"""
    res = ssc.run(pid, mode, base_inputs(pid), threads=1, power_curve=True)
    if not res.get("ok") or any("normal approximation; the answers marked" in t for t in res.get("notes", [])):
        pytest.skip("no curve computed as the answer is")
    for plot in [p for p in res["plots"] if p["kind"] == "curve"]:
        for series in plot["series"]:
            for x, y in series.get("markers", []):
                if x in series["x"]:
                    assert series["y"][series["x"].index(x)] == pytest.approx(y, abs=1e-7), (pid, mode, x)


def test_a_solved_f_is_the_effect_of_its_curve():
    """means.anova solved for f drew its curve with the means it was solved next to (0.913 at the answer's n,
    against 0.800, with unequal groups)."""
    res = ssc.run("means.anova", "f", {**base_inputs("means.anova"), "alloc": "ratio", "n": "20"}, power_curve=True)
    series = res["plots"][0]["series"][0]
    assert series["y"][series["x"].index(20)] == pytest.approx(res["rows"][0]["power"], abs=1e-9)
    ys = series["y"]
    assert all(a < b for a, b in zip(ys, ys[1:])), "the curve rises with n"


def test_a_simulation_curve_says_it_is_the_approximation():
    """Simulated answers sit near the curve, which is the analytic approximation; answers from the approximation
    itself lie on it, and the note is not given."""
    note = "analytic approximation; the simulated answers"
    inputs = {**base_inputs("sim.ordinal"), "reps": "100"}
    res = ssc.run("sim.ordinal", "power", inputs, power_curve=True)
    assert any(note in t for t in res["notes"])
    res = ssc.run("sim.ordinal", "power", inputs, power_curve=True, approximate=True)
    assert not any(note in t for t in res["notes"])


@pytest.mark.parametrize("pid,solve,inputs", [
    ("means.two.ni", "power", {"alloc": "ratio", "r": "1e30"}),
    ("prop.two.ni", "power", {"alloc": "ratio", "r": "1e30"}),
    ("rates.two.poisson", "power", {"alloc": "ratio", "r": "1e30"}),
    ("gs.logrank", "power", {"alloc": "ratio", "r": "1e30"}),
    ("prop.cc", "power", {"m": "1e30"}),
    ("means.two.ni", "n", {"alloc": "ratio", "r": "1e15"}),
    ("prop.cc", "n", {"m": "9007199254740993"}),
])
def test_a_count_past_any_study_is_said(pid, solve, inputs):
    """Group sizes and controls of infinity, or of more than 10^15, were given as answers ("20 and inf", at a ratio of
    1e308; the largest ratio allowed, 1e30, still gives counts past any study)."""
    row = run(pid, solve, {**base_inputs(pid), **inputs})["rows"][0]
    assert "a count would be infinite or beyond 10^15" in row.get("_error", ""), row


def test_a_median_of_about_zero_has_every_subject_with_the_event():
    """A subnormal median made the hazard infinite and the event probability inf/inf: "-nan events".  Such a median
    is refused; the smallest allowed, 1e-30, has every subject with the event."""
    for solve in ("power", "n"):
        rows = [first_row(run("gs.logrank", solve, {**base_inputs("gs.logrank"), "med2": m})) for m in ("1e-30", "1e-20")]
        assert not rows[0].get("_error") and "nan" not in rows[0]["_statement"], rows[0]
        assert rows[0]["events"] == rows[1]["events"] and rows[0]["power"] == rows[1]["power"]
    row = first_row(run("surv.logrank", "n", {**base_inputs("surv.logrank"), "med2": "1e-30"}))
    assert not row.get("_error") and "nan" not in row["_statement"], row


def test_maxcombo_at_a_small_alpha():
    """MaxCombo's critical value and power came from 1 − P(every statistic below c), which the lattice computes
    to about 1e-5: from an alpha of about 1e-3 the critical value fell short, its sample size dropped below that of
    the best of its four tests on its own (6182 against 6188 at 1e-30), and a study's power at 1e-30 was 0."""
    base = {**base_inputs("sim.surv.nph"), "alternative": "two"}

    def approx(solve, inputs):
        return first_row(ssc.run("sim.surv.nph", solve, {**base, **inputs}, power_curve=False, approximate=True))

    for alpha in ("1e-4", "1e-30"):
        best = min(approx("n", {"alpha": alpha, "test": "fh", "rho": r, "gamma": g})["n"]
                   for r, g in (("0", "0"), ("0", "1"), ("1", "0"), ("1", "1")))
        assert approx("n", {"alpha": alpha, "test": "maxcombo"})["n"] >= best, alpha
    assert approx("power", {"alpha": "1e-30", "test": "maxcombo"})["power"] > 0
    # with no effect, the power is the level
    row = approx("power", {"alpha": "1e-30", "test": "maxcombo", "hr1": "1", "hr2": "1"})
    assert row["power"] == pytest.approx(1e-30, rel=1e-6), row


def test_maxcombo_where_the_tail_computation_takes_over():
    """The two ways of computing MaxCombo's probability of leaving the box are blended between 0.005 and 0.02.  When
    the level, the power and the p-value each switched by their own value at 0.01, the power fell as alpha rose
    (0.33087 at 0.009999, 0.33054 at 0.01), and with no effect it was not alpha just below 0.01."""
    base = {**base_inputs("sim.surv.nph"), "alternative": "two", "test": "maxcombo"}

    def power(inputs):
        return first_row(ssc.run("sim.surv.nph", "power", {**base, **inputs}, power_curve=False,
                                 approximate=True))["power"]

    alphas = ("0.0049", "0.005", "0.0051", "0.0099", "0.009999", "0.01", "0.0101", "0.0199", "0.02", "0.0201")
    powers = [power({"alpha": a}) for a in alphas]
    assert powers == sorted(powers), list(zip(alphas, powers))
    for a in ("0.004", "0.009999", "0.012", "0.019", "0.03"):
        assert power({"alpha": a, "hr1": "1", "hr2": "1"}) == pytest.approx(float(a), rel=1e-6), a


def test_hazards_whose_sum_overflows_keep_their_share_of_events():
    """With event and dropout hazards of 1e308 each, their sum overflowed and the event probability came out 1
    (half the answer).  Such hazards are refused; up to the largest allowed, 1e30, the share of events is that of
    any equal pair of hazards."""
    rows = [first_row(run("surv.logrank", "n", {**base_inputs("surv.logrank"), "survin": "hazard", "lam2": h,
                                                "eta": h, "hr": "0.5"})) for h in ("1e3", "1e29", "1e30")]
    assert all(r is not None for r in rows), rows
    assert len({r["n"] for r in rows}) == 1, [r["n"] for r in rows]


@pytest.mark.parametrize("dropout", ["0.9999999999", "0.999999", "0.9999985"])
def test_a_dropout_table_never_shows_a_count_past_any_study(dropout):
    """Each group's inflated size, and the total's (which two groups just under 10^15 take past it)."""
    res = run("means.two.t", "power", {**base_inputs("means.two.t"), "n1": "1000000000", "dropout": dropout})
    assert first_row(res) is not None
    for row in (res.get("dropout") or {}).get("rows", []):
        assert all(abs(v) <= 1e15 for k, v in row.items() if k.endswith(("_inflated", "_dropouts"))), row


def test_hazards_whose_sum_overflows_follow_the_time_at_risk():
    """No time at risk refuses every pair of hazards, the largest too; a tiny one gives nearly the probability of
    nearly equal pairs.  (The pairs of 8.9e307 and 1e308, whose sum overflowed, are refused.)"""
    base = {**base_inputs("surv.logrank"), "survin": "hazard", "hr": "0.5"}
    for h in ("8.9e29", "1e30"):
        row = run("surv.logrank", "n", {**base, "lam2": h, "eta": h, "accrual": "0", "follow": "0"})["rows"][0]
        assert row.get("_error"), row
    rows = [first_row(run("surv.logrank", "n", {**base, "lam2": h, "eta": h, "accrual": "0", "follow": "1e-30"}))
            for h in ("8.9e29", "1e30")]
    assert abs(rows[0]["n"] - rows[1]["n"]) <= 0.15 * rows[0]["n"], rows


@pytest.mark.parametrize("pid,inputs", [
    ("prop.two", {}), ("prop.two", {"test": "fisher"}), ("prop.cc", {"method": "exact"}),
    ("prop.two.ni", {"method": "exact"}), ("prop.two.equiv", {"method": "exact"}),
])
def test_an_exact_search_for_a_power_its_sums_cannot_reach_ends(pid, inputs):
    """A target power of 1 - 1e-14 is past what the exact sums resolve: the search climbed toward 10^7 a group and
    ran for minutes (prop.two over 300 s); past its budget it now uses the normal approximation and says so."""
    start = time.perf_counter()
    row = first_row(run(pid, "n", {**base_inputs(pid), **inputs, "power": "0.99999999999999"}))
    assert time.perf_counter() - start < 30
    assert row is not None and any("normal approximation was used" in t for t in row.get("_notes", [])), row


def test_a_tiny_percentage_is_refused_at_once():
    start = time.perf_counter()
    row = run("means.two.welch", "n", {**base_inputs("means.two.welch"), "alloc": "pct", "pct1": "1e-30"})["rows"][0]
    assert row.get("_error") == "the percentage leaves one group empty" and time.perf_counter() - start < 0.5
    row = first_row(run("means.two.welch", "n", {**base_inputs("means.two.welch"), "alloc": "pct", "pct1": "0.01"}))
    assert row["n1"] >= 2 and row["n2"] >= 2


def test_a_zero_width_space_is_read_as_nothing():
    """"0<ZWSP>.4" looks like 0.4 in the field; it was read as the two values 0 and .4."""
    res = run("agree.icc", "n", {**base_inputs("agree.icc"), "rho0": "0\u200b.4"})
    assert [r["rho0"] for r in res["rows"]] == [0.4]


def test_a_robust_simulation_whose_values_overflow_ends():
    """An outlier SD multiple of 5e307 or more made the draws inf / inf: the midranks of a NaN never ended, and the
    run could not be cancelled natively (CLI, MCP time limits). Past about 1.3e154 the outliers' variance overflows
    and every value would be 0.  A multiple past 1e30, the largest an input may have, is refused; 1e30 is
    simulated."""
    for solve in ("power", "n"):
        start = time.perf_counter()
        res = run("sim.means.robust", solve, {**base_inputs("sim.means.robust"), "dist": "contaminated",
                                              "kout": "1e308", "reps": "100"})
        assert time.perf_counter() - start < 10
        assert res["input_errors"]["kout"] == "values must be at most 1e30", res
    row = first_row(run("sim.means.robust", "power", {**base_inputs("sim.means.robust"), "dist": "contaminated",
                                                      "kout": "1e30", "reps": "100"}))
    assert row is not None and row["exkurt"] is not None and row["exkurt"] > 0, row

@pytest.mark.parametrize("pid,solve", [("reg.cox", "n"), ("rates.two.poisson", "n"), ("means.two.t", "n"),
                                       ("corr.one", "n"), ("prop.paired.mcnemar", "n"), ("cluster.sw", "g"),
                                       ("rates.one.poisson", "n")])
def test_a_two_sided_alpha_whose_half_rounds_to_zero(pid, solve):
    """Half of an alpha of 5e-324 rounded to 0, whose critical value is infinite: the target "cannot be reached
    with a sample size up to 100000000".  Such an alpha is refused; the smallest allowed, 1e-30, needs more
    subjects than a larger one, and two-sided more than one-sided."""
    res = run(pid, solve, {**base_inputs(pid), "alpha": "5e-324"})
    assert res["input_errors"] == {"alpha": "values must be at least 1e-30"}, res
    base = {**base_inputs(pid), "alpha": "1e-30"}
    two = first_row(run(pid, solve, base))
    larger = first_row(run(pid, solve, {**base, "alpha": "1e-29"}))
    assert two and larger, (two, larger)
    key = next(k for k in ("n", "g") if k in two)
    assert two[key] > larger[key], (two, larger)
    if "alternative" in params(pid):
        one = first_row(run(pid, solve, {**base, "alternative": "upper"}))
        assert one and one[key] < two[key], (one, two)


@pytest.mark.parametrize("pid,solve,key", [("sim.lmm.cluster", "k", "k1"), ("sim.lmm.longitudinal", "n", "n1")])
def test_a_simulation_at_a_two_sided_alpha_whose_half_rounds_to_zero(pid, solve, key):
    res = run(pid, solve, {**base_inputs(pid), "alpha": "5e-324", "reps": "100"})
    assert res["input_errors"] == {"alpha": "values must be at least 1e-30"}, res
    rows = [first_row(run(pid, solve, {**base_inputs(pid), "alpha": a, "reps": "100"})) for a in ("1e-30", "1e-20")]
    assert rows[0] and rows[1] and rows[0][key] > rows[1][key], rows


def test_the_exact_test_of_a_correlation_at_a_subnormal_alpha():
    """t² overflowed beyond 1e154: a critical correlation of 0, and 4 subjects reaching any power.  Such alphas are
    refused; at the smallest allowed the test needs more subjects than at a larger one."""
    base = base_inputs("corr.one")
    for alpha in ("1e-320", "5e-324"):
        assert run("corr.one", "n", {**base, "alpha": alpha})["input_errors"] == {"alpha": "values must be at least 1e-30"}
    row = first_row(run("corr.one", "n", {**base, "alpha": "1e-30"}))
    assert row["n"] > first_row(run("corr.one", "n", {**base, "alpha": "1e-29"}))["n"], row


def test_a_cv_whose_square_underflows():
    """√ln(1 + CV²) underflowed to 0 below a CV of about 1e-162: no power at any size, and no sample size.  Such a
    CV is refused, but 1 + CV² is 1 in floating point down to the smallest allowed, 1e-30."""
    base = base_inputs("means.be.ratio")
    small = first_row(run("means.be.ratio", "n", {**base, "cv": "1e-15"}))
    assert small["n"] == 4
    for cv in ("1e-20", "1e-30"):
        row = first_row(run("means.be.ratio", "n", {**base, "cv": cv}))
        assert row["n"] == 4 and row["sw"] == float(cv), row
        row = first_row(run("means.be.ratio", "power", {**base, "cv": cv, "n": "24"}))
        assert row["power"] > 0.999999, row


@pytest.mark.parametrize("pid,key", [("prop.two.ci", "n1"), ("prop.one.ci", "n")])
@pytest.mark.parametrize("conf", ["1e-15", "1e-30"])
def test_an_interval_at_a_confidence_near_zero(pid, key, conf):
    """Two-sided limits close on the estimate as the confidence goes to 0, to a width of 0 or of a rounding error;
    the search took that for a degenerate size and said the target width "cannot be reached"."""
    for method in [o["value"] for o in params(pid)["method"]["options"]]:
        inputs = {**base_inputs(pid), "conf": conf, "method": method}
        row = first_row(run(pid, "n", inputs))
        assert row and 0 <= row["width"] <= float(inputs["width"]), (method, row)
        if method in ("newcombe", "wald", "wilson", "ac", "jeffreys"):
            assert row[key] == 1, (method, row)
    row = first_row(run("prop.two.ci", "width", {**base_inputs("prop.two.ci"), "conf": "1e-30", "n1": "100"}))
    assert row["width"] == 0


def test_squares_that_underflow_in_a_variance():
    """σ² and D² underflowed to 0: a stepped-wedge variance of 0/0, a search for a mean difference from 0 to 0, a
    slope variance σe²/Sxx of infinity.  Such σ and D (below 1e-30) are refused; near the smallest and largest
    allowed, the answers scale with the inputs."""
    base = base_inputs("cluster.sw")
    for solve, key in (("power", "power"), ("g", "g"), ("m", "m")):
        rows = [first_row(run("cluster.sw", solve, {**base, "sigma": repr(s), "delta": repr(0.3 * s)}))
                for s in (1, 1e-29, 1e29)]
        assert all(rows) and rows[1][key] == pytest.approx(rows[0][key], rel=1e-9), rows
        assert rows[2][key] == pytest.approx(rows[0][key], rel=1e-9), rows
    base = base_inputs("means.two.t")
    one = first_row(run("means.two.t", "mu1", {**base, "sigma": "1"}))
    for sigma in (1e-8, 1e-15, 1e-30):
        row = first_row(run("means.two.t", "mu1", {**base, "sigma": repr(sigma)}))
        assert row["mu1"] / sigma == pytest.approx(one["mu1"], rel=1e-9) and row["power"] == pytest.approx(0.8)
    base = {**base_inputs("means.two.slopes"), "sdb": "0"}  # the slopes' SD is then σe/√Sxx, which scales with 1/D
    one = first_row(run("means.two.slopes", "delta", {**base, "duration": "1"}))
    row = first_row(run("means.two.slopes", "delta", {**base, "duration": "1e-30"}))
    assert row["delta"] * 1e-30 == pytest.approx(one["delta"], rel=1e-9), row
    base = {**base, "timing": "custom"}
    one = first_row(run("means.two.slopes", "delta", {**base, "times": "0 1 2"}))
    row = first_row(run("means.two.slopes", "delta", {**base, "times": "0 1e-30 2e-30"}))
    assert row["delta"] * 1e-30 == pytest.approx(one["delta"], rel=1e-9), row


def test_an_input_too_small_to_compute_with_is_named():
    row = run("prop.cc", "power", {**base_inputs("prop.cc"), "m": "1e-15"})["rows"][0]
    assert "M is too small: there would be no controls" in row.get("_error", ""), row


def test_the_rule_for_the_number_of_controls_is_said_in_every_mode():
    """Controls are M × cases rounded up, a product within 1e-9 above a whole number counting as that number, so that
    one of 1e-9 or less leaves none: the message said "(M × cases, rounded up)", which makes such a product 1, and the
    n solve spoke of an allocation ratio and group 2."""
    base = {**base_inputs("prop.cc"), "method": "normal"}
    says = "M is too small: there would be no controls (M × cases is 1e-9 or less)"
    for m in ("1e-11", "1e-30"):
        for solve in ("power", "or"):
            row = run("prop.cc", solve, {**base, "m": m, "ncase": "100"})["rows"][0]
            assert row.get("_error") == says, (solve, m, row)
        row = run("prop.cc", "n", {**base, "m": m})["rows"][0]
        assert row.get("_error") == "M is too small: even 10000000 cases would have fewer than 2 controls", (m, row)
    assert first_row(run("prop.cc", "power", {**base, "m": "1e-11", "ncase": "101"}))["nctrl"] == 1


@pytest.mark.parametrize("pid,solve,inputs,key,says", [
    # values that a method once computed with as 0 or infinity: a power of 1 at proportions near 0 ("8 cases and 8
    # controls"), a power of 1 from ∞/∞ ("group sample sizes of 1 and 1"), boundaries that spend nothing, an
    # infinite hazard, a slope variance of infinity, fits that all failed; each later said it was too small
    ("prop.cc", "n", {"p0": "5e-324"}, "p0", "values must be at least 1e-30"),
    ("prop.two", "n", {"es": "or", "p2": "5e-324"}, "p2", "values must be at least 1e-30"),
    ("rates.two.poisson", "n", {"time": "5e-324"}, "time", "values must be between 1e-30 and 1e30"),
    ("rates.two.poisson", "power", {"lam2": "5e-324"}, "lam2", "values must be between 1e-30 and 1e30"),
    ("rates.two.negbin", "n", {"r0": "5e-324"}, "r0", "values must be between 1e-30 and 1e30"),
    ("gs.means", "n", {"alpha": "1e-300"}, "alpha", "values must be at least 1e-30"),
    ("gs.props", "power", {"alpha": "5e-324", "boundary": "wt_obf"}, "alpha", "values must be at least 1e-30"),
    ("sim.winratio", "power", {"med_d": "5e-324"}, "med_d", "values must be between 1e-30 and 1e30"),
    ("sim.winratio", "power", {"med_h": "5e-324"}, "med_h", "values must be between 1e-30 and 1e30"),
    ("means.two.slopes", "delta", {"duration": "1e-160"}, "duration", "values must be between 1e-30 and 1e30"),
    ("means.two.slopes", "delta", {"timing": "custom", "times": "0 1e-160 2e-160"}, "times",
     "values must be 0 or between 1e-30 and 1e30 in size"),
    ("sim.lmm.longitudinal", "power", {"sde": "1e-100"}, "sde", "values must be between 1e-30 and 1e30"),
    ("sim.means.robust", "power", {"sd1": "1e300"}, "sd1", "values must be between 1e-30 and 1e30"),
    # each kind of limit: 0 allowed or not, either sign, one end or both; lists, ranges and lists of lists
    ("means.two.t", "power", {"mu1": "-1e31"}, "mu1", "values must be 0 or between 1e-30 and 1e30 in size"),
    ("means.two.t", "n", {"dropout": "0 1e-31"}, "dropout", "values must be 0 or at least 1e-30"),
    ("means.two.t", "power", {"sigma": "1e-31 to 1 by 0.5"}, "sigma", "values must be between 1e-30 and 1e30"),
    # each end of a range as typed, as it would be alone: one its values do not reach, and a start kept as typed
    # (rounded to 12 digits, it would be 1e-30)
    ("means.two.t", "power", {"mu1": "0 to 1.5e30 by 1e30"}, "mu1", "values must be 0 or between 1e-30 and 1e30 in size"),
    ("means.two.t", "power", {"mu1": "9.9999999999999e-31 to 9.9999999999999e-31"}, "mu1",
     "values must be 0 or between 1e-30 and 1e30 in size"),
    ("corr.one", "n", {"rho1": "-1e-40"}, "rho1", "values must be 0 or at least 1e-30 in size"),
    ("prop.two.equiv", "n", {"dl": "-1e-31"}, "dl", "values must be at least 1e-30 in size"),
    ("prop.two.equiv", "n", {"scale": "or", "oru": "1e31"}, "oru", "values must be at most 1e30"),
    ("means.anova", "n", {"means": "0 0.5 1; 0 1e31 2"}, "means", "values must be 0 or between 1e-30 and 1e30 in size"),
])
def test_a_value_below_1e_30_or_above_1e30_in_size_is_refused(pid, solve, inputs, key, says):
    """Every number that is not a count is 0 or between 1e-30 and 1e30 in size: no study's quantity lies beyond in
    any sensible unit, and within it the methods' products, squares and ratios stay far from overflow and underflow.
    The refusal names the input's limits."""
    res = run(pid, solve, {**base_inputs(pid), **inputs})
    assert not res["ok"] and res["input_errors"] == {key: says}, res


def test_values_at_the_size_limits_are_allowed():
    res = run("means.two.t", "power", {**base_inputs("means.two.t"), "sigma": "1e-30 1e30", "mu1": "0 -1e30 1e-30 -0"})
    assert res["ok"] and len(res["rows"]) == 8 and all("_error" not in r for r in res["rows"]), res
    res = run("means.anova", "n", {**base_inputs("means.anova"), "means": "0 -1e-30 1e30"})
    assert res["ok"] and first_row(res), res
    # counts keep their own limit of a billion
    res = run("corr.one", "power", {**base_inputs("corr.one"), "n": "1000000000"})
    assert res["ok"] and first_row(res), res


@pytest.mark.parametrize("pid,inputs,key", [
    ("means.two.t", {"mu1": "-0.3 to 0.3 by 0.1"}, "mu1"),
    ("means.two.t", {"mu1": "-3e-15 to 3e-15 by 1e-15", "sigma": "2e-15"}, "mu1"),
    ("means.two.t", {"mu1": "0.3e-29 to -0.3e-29 by 1e-30"}, "mu1"),
    ("means.crossover", {"delta": "-3e-15 to 3e-15 by 1e-15"}, "delta"),
])
def test_a_range_that_crosses_0_has_0_among_its_values(pid, inputs, key):
    """Where a range crosses 0, a + k * step leaves rounding noise: the label 5.55111512313e-17 for '-0.3 to 0.3 by
    0.1', and 3.9e-31 for '-3e-15 to 3e-15 by 1e-15', which was refused as below 1e-30 though the user typed no such
    value.  A value within a billionth of a step of 0 is 0."""
    res = run(pid, "power", {**base_inputs(pid), **inputs})
    assert res["ok"] and len(res["rows"]) == 7 and res["rows"][3][key] == 0, res
    assert all("_error" not in r for r in res["rows"]), res


@pytest.mark.parametrize("pid,solve,inputs,says", [
    # the start of a range, and its end where the steps reach it, as typed: checked as each would be alone (rounded to
    # 12 digits, 5.0000000000001 was 5 and 0.99999999999999 was 1)
    ("sim.lmm.cluster", "power", {"cv": "5.0000000000001 to 5.0000000000001"}, {"cv": "values must be at most 5"}),
    ("means.two.t", "power", {"n1": "10.0000000000001 to 12"}, {"n1": "values must be whole numbers"}),
    ("means.two.t", "power", {"mu1": "0.5 to 0.9999999999999 by 0.4999999999999", "mu2": "1"}, [0.5, 0.9999999999999]),
    ("means.two.t", "n", {"power": "0.99999999999999 to 0.99999999999999"}, [0.99999999999999]),
    ("means.two.t", "n", {"power": "0.5 to 1 by 0.2"}, [0.5, 0.7, 0.9]),  # an end the steps do not reach is no value
    # a range that ends at 0 has 0 among its values, whichever side the rounding noise of the sum fell on (0.9 - 3 *
    # 0.3 is 1.1e-16, 0.3 - 3 * 0.1 is -5.6e-17)
    ("means.two.t", "n", {"alpha": "0.9 to 0 by 0.3"}, {"alpha": "values must be greater than 0"}),
    ("means.two.t", "n", {"alpha": "0.3 to 0 by 0.1"}, {"alpha": "values must be greater than 0"}),
    ("means.two.t", "n", {"dropout": "0.9 to 0 by 0.3"}, [0.9, 0.6, 0.3, 0]),
])
def test_the_ends_of_a_range_are_values_as_typed(pid, solve, inputs, says):
    res = run(pid, solve, {**base_inputs(pid), **inputs})
    if isinstance(says, dict):
        assert not res["ok"] and res["input_errors"] == says, res
    else:
        key = next(iter(inputs))
        assert res["ok"] and [r.get("target_" + key, r.get(key)) for r in res["rows"]] == says, res


def test_a_list_of_group_means_that_crosses_0_is_computed():
    res = run("means.anova", "n", {**base_inputs("means.anova"), "means": "-3e-15 to 3e-15 by 1e-15", "sigma": "1e-15"})
    assert res["ok"] and first_row(res), res


@pytest.mark.parametrize("pid,solve,alpha", [("gs.logrank", "power", "1e-10"), ("gs.means", "power", "1e-12"),
                                             ("gs.means", "n", "1e-15"), ("gs.props", "n", "1e-15")])
def test_a_group_sequential_design_at_a_small_alpha(pid, solve, alpha):
    """The other rejection region adds less than a rounding error to the fixed design's power at a small alpha,
    which left its drift unbracketed ("the answer cannot be computed")."""
    row = first_row(run(pid, solve, {**base_inputs(pid), "alpha": alpha}))
    assert row and row["inflation"] == pytest.approx(1, abs=1e-3), row


@pytest.mark.parametrize("pid", ["gs.means", "gs.props", "gs.logrank"])
def test_power_family_spending_with_a_tiny_exponent(pid):
    """Below a ρ of about 1e-16, t^ρ rounds to 1, and α·t^ρ − α·s^ρ to 0, as if all of alpha were spent at the first
    analysis; the later analyses spend about α·ρ·ln(t/s), at boundaries far out that carry most of the power when the
    first analysis has little information.  References: the two-analysis design integrated with scipy (quad, relative
    error 1e-12).  Spent as rounded, '1e-4 1' with ρ = 1e-30 gave an inflation of 10000, and '0.01 1' with 1e-17 one
    of 100."""
    base = {**base_inputs(pid), "boundary": "power", "power": "0.8", "alpha": "0.05", "alternative": "two"}
    rows = [first_row(run(pid, "n", {**base, "looks": "3", "spacing": "equal", "rho": rho}))
            for rho in ("1e-15", "1e-16", "1e-30")]
    assert all(r and r["n"] == rows[0]["n"] for r in rows), rows
    assert all(r["inflation"] == pytest.approx(3, rel=1e-4) for r in rows), rows  # the first analysis's own power
    for tfrac, rho, inflation in (("0.01 1", "1e-17", 11.414589), ("1e-4 1", "1e-30", 19.556845),
                                  ("0.05 1", "1e-15", 9.489303)):
        row = first_row(run(pid, "n", {**base, "looks": "2", "spacing": "custom", "tfrac": tfrac, "rho": rho}))
        assert row and row["inflation"] == pytest.approx(inflation, rel=1e-5), (tfrac, rho, row)


def test_flat_spending_functions_with_many_analyses():
    """Twenty analyses of a nearly flat spending function: the power family with ρ = 1e-30 and Hwang–Shih–DeCani with
    γ = 40.  References: recursive integration with Simpson's rule on 2000 uniform intervals per analysis and the
    increments taken without cancellation (inflation 16.887305 and 7.494518).  With rounded increments they were 20
    (N 2512 for 2122) and 7.5079; with the grid of an alpha of 0.05 (r = 32), 16.867 (N 2120)."""
    base = {**base_inputs("gs.means"), "power": "0.8", "alpha": "0.05", "alternative": "two", "looks": "20",
            "spacing": "equal"}
    for extra, inflation in (({"boundary": "power", "rho": "1e-30"}, 16.887305),
                             ({"boundary": "hsd", "gamma": "40"}, 7.494518)):
        row = first_row(run("gs.means", "n", {**base, **extra}))
        assert row and row["inflation"] == pytest.approx(inflation, rel=1e-5), (extra, row)


def test_group_sequential_analyses_close_together_at_a_small_alpha():
    """Analyses closer than 20 equally spaced ones narrow the step from one to the next: below an alpha of 1e-8 such
    designs keep the grid of r = 128, as before the grid was graded (with 64, N was 2 more than with r = 256 here)."""
    base = {"stat": "pooled", "alpha": "1e-9", "spacing": "custom"}
    for extra, n in (({"looks": "3", "tfrac": "0.98 0.99 1", "boundary": "wt_obf"}, 2100),
                     ({"looks": "5", "tfrac": "0.2 0.4 0.6 0.999 1", "boundary": "hsd", "gamma": "-4"}, 2114)):
        row = first_row(run("gs.props", "n", {**base, **extra}))
        assert row and row["n"] == n, (extra, row)


@pytest.mark.parametrize("pid,solve,scale,fixed", [
    ("means.crossover", "delta", "sigma", {}), ("means.one.t", "mu1", "sigma", {"mu0": "0"}),
    ("means.paired.t", "mud", "sdd", {"sdtype": "sdd", "mud0": "0"}),
    ("means.two.repeated", "mu1", "sigma", {"mu2": "0"}), ("means.two.ancova", "mu1", "sigma", {"mu2": "0"}),
    ("means.two.cluster", "mu1", "sigma", {"mu2": "0"}), ("reg.slope", "b1", "sigma", {"b0": "0"})])
def test_an_effect_solved_at_a_small_sd_scales_with_it(pid, solve, scale, fixed):
    """The search's tolerance was 1e-12 absolute below 1: at an SD of 1e-15 the effect was the first value tried
    (the SD itself), and at 1e-8 it was off in the fourth digit."""
    base = {**base_inputs(pid), **fixed, scale: "1"}
    unit = first_row(run(pid, solve, base))[solve]
    for sd in (1e-8, 1e-15, 1e-30):
        row = first_row(run(pid, solve, {**base, scale: repr(sd)}))
        assert row and row[solve] / sd == pytest.approx(unit, rel=1e-6), (sd, row)


@pytest.mark.parametrize("pid,inputs,solve", [
    ("means.two.welch", {"sd1": 1, "sd2": 1.5, "mu1": 1, "mu2": 0}, "n"),
    ("means.two.welch", {"sd1": 1, "sd2": 1.5, "mu2": 0}, "mu1"),
    ("means.paired.t", {"sdtype": "corr", "sd1": 1, "sd2": 1, "mud": 0.5, "mud0": 0}, "n"),
    ("means.two.slopes", {"sdb": 0.5, "sde": 1, "delta": 0.5}, "n"),
    ("means.two.slopes", {"sdb": 0.5, "sde": 1}, "delta")])
def test_sds_whose_squares_underflow_together(pid, inputs, solve):
    """Every SD and mean difference scaled by 1e-80 or 1e-300: the Welch degrees of freedom (from the variances'
    squares, which underflow below SDs of 1e-77), the paired SD and the slopes' SD gave other answers, or none.
    Such SDs are refused; scaled by the smallest and largest factors that keep every input within 1e-30 to 1e30,
    the answers are those of the inputs as they are."""
    base = base_inputs(pid)
    one = first_row(run(pid, solve, {**base, **{k: str(v) for k, v in inputs.items()}}))
    for f in (2e-30, 5e29):
        scaled = {k: repr(v * f) if isinstance(v, (int, float)) else v for k, v in inputs.items()}
        row = first_row(run(pid, solve, {**base, **scaled}))
        assert row, (f, row)
        if solve == "n":
            assert row["n"] == one["n"] and row["power"] == pytest.approx(one["power"], rel=1e-9), (f, row, one)
        else:
            assert row[solve] / f == pytest.approx(one[solve], rel=1e-9), (f, row, one)


@pytest.mark.parametrize("pid", ["sim.lmm.cluster", "sim.glmm.cluster"])
def test_a_cluster_size_cv_whose_square_underflows(pid):
    """1/CV² was infinite below a CV of about 1e-154 and every cluster had 10⁹ subjects (a simulated power of
    100% at the defaults of sim.lmm.cluster, against 86% at a CV of 1e-8).  Such a CV is refused; at the smallest
    allowed, 1e-30, the sizes all round to the mean, as at a CV of 1e-8 (whose draws take the same random numbers)
    and with the approximation of a CV of 0."""
    base = {**base_inputs(pid), "reps": "100"}
    zero, small, smallest = (first_row(run(pid, "power", {**base, "cv": cv})) for cv in ("0", "1e-8", "1e-30"))
    assert smallest["power"] == small["power"], (smallest, small)
    assert smallest["approx_power"] == pytest.approx(zero["approx_power"], rel=1e-12), (smallest, zero)


@pytest.mark.parametrize("key", ["med_d", "med_h"])
def test_a_win_ratio_median_of_about_zero(key):
    """An infinite hazard (a median of 5e-324, now refused) made every simulated time no number: a simulated power
    of 0 with no approximation.  The smallest median allowed, 1e-30, is simulated and approximated."""
    row = first_row(run("sim.winratio", "power", {**base_inputs("sim.winratio"), key: "1e-30", "reps": "100"}))
    assert row and row["failures"] == 0 and row["power"] == pytest.approx(row["approx_power"], abs=0.15), row


def test_an_implied_p1_that_rounds_to_one_is_named():
    row = run("prop.two.ni", "n", {**base_inputs("prop.two.ni"), "es": "or", "or1": "1e30"})["rows"][0]
    assert row["_error"] == "the implied P1 = 1 is not strictly between 0 and 1", row


@pytest.mark.parametrize("method", ["normal", "exact"])
@pytest.mark.parametrize("p1,p2", [(0.6, 0.6), (0.3, 0.2), (0.95, 0.9), (0.05, 0.1), (0.999, 0.998)])
def test_widening_an_odds_ratio_margin_never_lowers_the_power(method, p1, p2):
    """Restricted estimates of the wrong sign gave 0% power once an upper limit was some 10^13 from 1, a 1 − P2
    that lost its digits lowered the exact power from a lower limit of 10^-20 on, and a non-inferiority margin whose
    boundary P1 rounds to 1 was refused."""
    odds = p1 / (1 - p1) / (p2 / (1 - p2))
    common = {"scale": "or", "p1": str(p1), "p2": str(p2), "n1": "60", "method": method}
    steps = [10.0 ** e for e in (0.2, 1, 4, 8, 12, 13, 14, 16, 20, 30)]
    lo, hi = min(odds, 1) / 1.5, max(odds, 1) * 1.5  # the limits lie on either side of 1
    for side in ("orl", "oru"):
        powers = []
        for f in steps:
            limit = lo / f if side == "orl" else hi * f
            if limit < 1e-30 or limit > 1e30:
                continue
            limits = {"orl": repr(lo), "oru": repr(hi), side: repr(limit)}
            row = first_row(run("prop.two.equiv", "power", {**common, **limits}))
            assert row, (side, limit)
            powers.append(row["power"])
        assert all(b >= a - 1e-12 for a, b in zip(powers, powers[1:])), (side, powers)
    for h1, sign in (("upper", -1), ("lower", 1)):
        powers = []
        for f in steps:
            margin = odds * f ** sign
            if margin < 1e-30 or margin > 1e30:
                continue
            row = first_row(run("prop.two.ni", "power", {**common, "h1": h1, "or0": repr(margin)}))
            assert row, (h1, margin)
            powers.append(row["power"])
        assert all(b >= a - 1e-12 for a, b in zip(powers, powers[1:])), (h1, powers)


def test_an_exact_odds_ratio_test_far_from_one_rejects_the_tables_it_should():
    """With P1 = 0.999 and P2 = 0.998 almost every table has a group whose subjects all responded.  For H1: OR >
    10^-13 the test rejects every table but those where all of group 2 responded (power 1 − 0.998^25), and for
    H1: OR > 10^13 none: the numerator ph1 − pt1, a rounding error off 0 with no information to weigh it, gave
    99.97% and 100%."""
    base = {"scale": "or", "method": "exact", "p1": "0.999", "p2": "0.998", "n1": "25", "h1": "upper", "alpha": "0.025"}
    row = first_row(run("prop.two.ni", "power", {**base, "or0": "1e-13"}))
    assert row["power"] == pytest.approx(1 - 0.998 ** 25, abs=1e-12), row
    row = first_row(run("prop.two.ni", "power", {**base, "or0": "1e13"}))
    assert row["power"] < 1e-12, row


def test_an_odds_ratio_margin_with_proportions_near_one_counts_no_table_where_everyone_responded():
    """The table where every subject of both groups responded rejected on an infinite statistic before it was made
    one without information: at P1 = P2 = 0.95, 20 per group and OR0 = 0.302 the exact power was 0.25488, where a
    60-digit enumeration gives 0.12637."""
    row = first_row(run("prop.two.ni", "power", {"scale": "or", "method": "exact", "p1": "0.95", "p2": "0.95",
                                                 "n1": "20", "or0": "0.302", "h1": "upper", "alpha": "0.025"}))
    assert row["power"] == pytest.approx(0.1263705114339569, abs=1e-12), row


@pytest.mark.parametrize("h1,d0", [("upper", "-5e-17"), ("upper", "-1e-20"), ("upper", "-1e-30"), ("lower", "1e-20")])
def test_a_difference_margin_next_to_zero_leaves_the_table_where_everyone_responded_unrejected(h1, d0):
    """At the boundary the table where every subject responded has P1 = 1 + D0 and P2 = 1: 1 − P1 = −D0 was lost to
    rounding below about 5.6e-17, so that the table had no variance and its numerator −D0 an infinite statistic.  The
    exact power and size at P1 = P2 = 0.95 and 20 per group were 0.1343 (0.0058 + 0.95^40), not the 0.0058 of a
    margin of −1e-16, and an equivalence test with limits of ±1e-20 had a power of 0.95^40."""
    base = {"scale": "diff", "method": "exact", "p1": "0.95", "p2": "0.95", "n1": "20", "alpha": "0.025"}
    near = first_row(run("prop.two.ni", "power", {**base, "h1": "upper", "d0": "-1e-16"}))
    row = first_row(run("prop.two.ni", "power", {**base, "h1": h1, "d0": d0}))
    assert row["power"] == pytest.approx(near["power"], abs=1e-12), (row, near)
    assert row["actual_alpha"] == pytest.approx(near["actual_alpha"], abs=1e-12), (row, near)
    equiv = first_row(run("prop.two.equiv", "power", {**base, "alpha": "0.05", "dl": "-1e-20", "du": "1e-20"}))
    assert equiv["power"] < 1e-12, equiv


def test_a_given_difference_is_judged_and_reported_as_given():
    """P1 = P2 + D1 loses a D1 of 1e-20 to rounding at P2 = 0.5: the n solve judged the actual effect from P1 − P2 = 0
    and refused it as lying outside H1: P1 − P2 > 1e-25 (the power cannot be reached, which it now says), and the
    row's D1 was 0."""
    base = {"scale": "diff", "es": "diff", "p2": "0.5", "d1": "1e-20", "d0": "1e-25", "h1": "upper", "alpha": "0.025"}
    for method in ("normal", "exact"):
        row = run("prop.two.ni", "n", {**base, "method": method, "power": "0.8"})["rows"][0]
        assert "the target power cannot be reached" in row.get("_error", ""), (method, row)
    assert first_row(run("prop.two.ni", "power", {**base, "method": "normal", "n1": "100"}))["d1"] == 1e-20


def test_a_superiority_test_with_a_negligible_difference_margin_needs_the_same_groups():
    """D0 = −1e-20 gave 2 per group at an actual alpha of 0.96 (the table where every subject responded rejected)."""
    base = {"scale": "diff", "method": "exact", "p1": "0.999", "p2": "0.99", "h1": "upper", "alpha": "0.025",
            "power": "0.8"}
    a, b = (first_row(run("prop.two.ni", "n", {**base, "d0": d0})) for d0 in ("-1e-15", "-1e-20"))
    assert a["n1"] == b["n1"] == 910 and b["actual_alpha"] < 0.025, (a, b)


@pytest.mark.parametrize("scale,inputs,power", [
    ("diff", {"p1": "1e-20", "p2": "1e-20", "d0": "-1e-25"}, 0.024999999998568537),
    ("diff", {"p1": "0.9999999999", "p2": "0.9999999999", "d0": "-1e-25"}, 0.025),
    ("diff", {"p1": "0.999999", "p2": "0.9999995", "d0": "-1e-7", "n1": "10000000"}, 0.0014652489849946286),
    ("diff", {"p1": "0.9999", "p2": "0.9999", "d0": "0", "n1": "1000"}, 0.025),
    ("diff", {"p1": "0.99986", "p2": "0.999876", "d0": "-1.23e-6", "n1": "100000", "h1": "lower"},
     0.047239565537847697),
    ("ratio", {"p1": "1e-20", "p2": "1e-20", "r0": "0.8"}, 0.026439090706363813)])
def test_score_tests_at_proportions_near_0_or_1_by_the_normal_approximation(scale, inputs, power):
    """The restricted estimates lost their digits near 0 and 1 (the closed-form root of the difference's cubic is
    off by 1e-16 or more whatever its size, and the smaller root of the ratio's quadratic cancelled to 0): at
    P1 = P2 = 1e-20 the power was 50% where it is about α.  Estimates 1e-4 to 1e-3 from 0 or 1 kept seven or eight
    digits: at P1 = P2 = 0.9999 and D0 = 0 the power was 0.024999998 where it is α.  The references evaluate the same
    formula to 60 or 120 digits."""
    row = first_row(run("prop.two.ni", "power", {"scale": scale, "method": "normal", "n1": "100", "h1": "upper",
                                                 "alpha": "0.025", **inputs}))
    assert row["power"] == pytest.approx(power, rel=1e-12, abs=0), row


@pytest.mark.parametrize("method,inputs", [
    ("exact", {"p1": "0.001", "p2": "0.2", "n1": "50", "r0": "1e-17", "h1": "lower"}),
    ("exact", {"p1": "1e-25", "p2": "0.2", "n1": "50", "r0": "1e-17", "h1": "lower"}),
    ("exact", {"p1": "0.5", "p2": "1e-19", "n1": "30", "r0": "1e18", "h1": "upper"}),
    ("normal", {"p1": "1e-25", "p2": "0.2", "n1": "50", "r0": "1e-17", "h1": "lower"})])
def test_a_ratio_margin_far_from_one_rejects_no_table_it_should_not(method, inputs):
    """(−B − √D)/(2A), the smaller root of the restricted estimates' quadratic, cancelled to 0 for a margin below about
    1e-16 (or above 1e16 with a small P2): the table had no variance and an infinite statistic.  The exact power and
    size were 95% and 99.999%, and the normal approximation 50%, where a 120-digit enumeration gives 0 (the actual
    ratio lies far inside H0)."""
    row = first_row(run("prop.two.ni", "power", {"scale": "ratio", "method": method, "alpha": "0.025", **inputs}))
    assert row["power"] < 1e-12 and row.get("actual_alpha", 0) < 1e-12, row


@pytest.mark.parametrize("inputs,size", [
    ({"scale": "or", "p1": "0.6", "p2": "0.6", "n1": "4", "or0": "1e30", "h1": "lower"}, 2.5984000000000001736e-30),
    ({"scale": "or", "p1": "0.999", "p2": "0.998", "n1": "25", "or0": "1e13", "h1": "lower"}, 5.01002004008015173e-15),
    ({"scale": "or", "p1": "0.6", "p2": "0.6", "n1": "20", "or0": "1e-30", "h1": "upper"}, 2.9998903152467980833e-29),
    ({"scale": "ratio", "p1": "0.3", "p2": "0.2", "n1": "30", "r0": "1e-20", "h1": "upper"}, 6e-20),
    ({"scale": "diff", "alloc": "sep", "n1": "6", "n2": "7", "p1": "0.0001", "p2": "1e-14", "d0": "0.9999999999999",
      "h1": "lower"}, 6.1018656711219812871e-13),
    ({"scale": "ratio", "alloc": "sep", "n1": "6", "n2": "7", "p1": "0.0001", "p2": "0.7", "r0": "1.4285714285714",
      "h1": "lower"}, 9.9279087084619924668e-15)])
def test_the_exact_size_is_taken_at_the_boundary_itself(inputs, size):
    """P1 was kept within 1e-9 of 0 and 1, so that the size at OR0 = 1e30 was 3.9e-9 where a 150-digit enumeration
    gives 2.6e-30.  On the difference and ratio scales 1 − P1 there carried the rounding of 1 − P2 (and of R0 × P2):
    sizes of 6.1023e-13 and 9.9295e-15 where 130-digit enumerations give 6.1019e-13 and 9.9279e-15."""
    row = first_row(run("prop.two.ni", "power", {"method": "exact", "alpha": "0.025", **inputs}))
    assert row["actual_alpha"] == pytest.approx(size, rel=1e-9, abs=0), row


def test_the_exact_size_where_p1_rounds_to_one_counts_every_table():
    """At OR0 = 1e17, P1 on the boundary is a rounding error from 1, and group 1's binomial window, set from P1 = 1,
    held only N1 − 1 and N1 responses: the tables with fewer, of probability about (N1·(1 − P1))²/2, were left out
    (3.5e-13 of the size at 100,000 per group, where a 60-digit enumeration finds every table with fewer than N1
    responses rejecting)."""
    row = first_row(run("prop.two.ni", "power", {"scale": "or", "method": "exact", "p1": "0.6", "p2": "0.6",
                                                 "n1": "100000", "h1": "lower", "alpha": "0.025", "or0": "1e17"}))
    assert row["actual_alpha"] == pytest.approx(6.666666666664445039e-13, rel=1e-13, abs=0), row


@pytest.mark.parametrize("pid,inputs", [
    ("prop.cc", {"es": "or", "or": "1e20"}),
    ("prop.two", {"es": "or", "or1": "1e14"}),
    ("prop.two.equiv", {"scale": "or", "es": "or", "or1": "1e12", "orl": "0.5", "oru": "1e13"})])
def test_a_given_odds_ratio_is_reported_as_given(pid, inputs):
    """Recomputed from a P1 near 1, the odds ratio lost its digits with 1 − P1 (and was infinite once P1 rounded to
    1): an odds ratio of 10^12 was reported as 9.9999 × 10^11."""
    row = first_row(run(pid, "power", {**base_inputs(pid), **inputs}))
    key = "or1" if pid == "prop.two.equiv" else "or"
    assert row and row[key] == float(inputs["or1" if "or1" in inputs else "or"]), row


def test_a_cluster_ratio_that_leaves_group_2_no_clusters_is_refused():
    """k2 = k1 × the ratio, rounded up, was 0: a power of 100% (from min(1, NaN + NaN)), as the k solve refuses.  A
    product within 1e-9 above a whole number counts as that number, so that one of 1e-9 or less leaves group 2 none,
    which the message of either mode says (the power's spoke of rounding up, which makes such a product 1)."""
    base = base_inputs("prop.two.cluster")
    says = "the cluster ratio is too small: group 2 would have no clusters (k1 × the cluster ratio is 1e-9 or less)"
    for kr in ("1e-12", "5e-11"):
        row = run("prop.two.cluster", "power", {**base, "k1": "10", "kr": kr})["rows"][0]
        assert row.get("_error") == says, row
    row = run("prop.two.cluster", "k", {**base, "kr": "1e-30"})["rows"][0]
    assert row.get("_error") == says, row
    assert first_row(run("prop.two.cluster", "power", {**base, "k1": "11", "kr": "1e-10"}))["k2"] == 1


def test_an_allocation_ratio_that_leaves_group_2_of_an_interval_no_subjects_is_refused():
    """N2 = R × N1, rounded up, was 0: the Wald limits were 0/0, which their clamps made −1 and 1, a width of 2 given
    as an answer (Agresti–Caffo gave a finite interval for the empty group, Newcombe no answer), and a one-sided bound
    pressed against 1 met a target width with 1 and 0 subjects.  As with the cluster ratio, a product of 1e-9 or less
    leaves group 2 no subjects, and both modes say so."""
    base = base_inputs("prop.two.ci")
    says = "the allocation ratio is too small: group 2 would have no subjects (R × N1 is 1e-9 or less)"
    for method in ("wald", "waldcc", "ac", "newcombe"):
        row = run("prop.two.ci", "width", {**base, "method": method, "n1": "10", "r": "5e-11"})["rows"][0]
        assert row.get("_error") == says, (method, row)
    row = run("prop.two.ci", "n", {**base, "r": "1e-30"})["rows"][0]
    assert row.get("_error") == says, row
    row = first_row(run("prop.two.ci", "n", {**base, "method": "wald", "sides": "upper", "p1": "0.3", "p2": "0.2",
                                             "r": "5e-11", "width": "1"}))
    assert row["n1"] == 21 and row["n2"] == 1, row
    assert first_row(run("prop.two.ci", "width", {**base, "n1": "11", "r": "1e-10"}))["n2"] == 1


def test_mcnemars_normal_approximation_where_p10_rounds_to_one():
    """ψ − (P10 − P01)², the variance of a pair's difference, cancelled to 0 once P10 rounded to 1 (ψ = 1 and an odds
    ratio of 1e30): with a critical value of exactly √n the numerator was 0 as well, and min(1, no number) made the
    power 100%, which the n solve took (the one-sided test had no answer).  The power there is 50% (0.5 − 8e-16, the
    formula to 60 digits).  ψ = 1 could also come out a rounding error above 1 from P10 and P01, refused as P10 + P01
    exceeding 1."""
    base = {**base_inputs("prop.paired.mcnemar"), "method": "normal", "es": "or", "psi": "1", "orr": "1e30"}
    for alternative, alpha in (("two", "0.0455002638963584"), ("upper", "0.0227501319481792")):
        row = first_row(run("prop.paired.mcnemar", "power", {**base, "alternative": alternative, "alpha": alpha,
                                                              "n": "4"}))
        assert row["power"] == pytest.approx(0.5, abs=1e-12), (alternative, row)
    for method in ("exact", "normal"):
        row = first_row(run("prop.paired.mcnemar", "power", {**base, "method": method, "orr": "9.08406134939307e-05",
                                                              "alternative": "two", "alpha": "0.05", "n": "50"}))
        assert row and row["psi"] == 1, (method, row)


@pytest.mark.parametrize("odds_ratio", ["1e16", "1e20", "1e30"])
def test_a_matched_case_control_study_with_an_odds_ratio_far_above_one(odds_ratio):
    """1 − P1 and the chance that the case is the exposed member of a set, 1 − e, rounded to 0 once the odds ratio
    passed about 1e16: the variance under H1 was 0 (no answer), and beyond the exposure of a control whose case was
    not exposed was no number, refused as too strong a correlation φ, at φ = 0."""
    base = {**base_inputs("prop.cc.matched"), "phi": "0"}
    near = first_row(run("prop.cc.matched", "power", {**base, "or": "1e10"}))
    row = first_row(run("prop.cc.matched", "power", {**base, "or": odds_ratio}))
    assert row and row["or"] == float(odds_ratio), row
    assert row["power"] == pytest.approx(near["power"], abs=1e-12), (row, near)


def test_text_that_acts_on_a_message_is_written_out():
    """A right-to-left override typed in a field turned the rest of the error line around."""
    res = run("means.two.ni", "n", {**base_inputs("means.two.ni"), "power": "\u202e0.8", "sigma": "1\u2066"})
    assert res["input_errors"]["power"] == "could not read '\\u{202E}0.8' as a number"
    assert "\\u{2066}" in res["input_errors"]["sigma"] and "\u2066" not in res["input_errors"]["sigma"]


def test_an_unknown_choice_is_an_input_error():
    res = run("means.two.t", "power", {**base_inputs("means.two.t"), "alternative": "greater"})
    assert not res["ok"] and "unknown choice 'greater'" in res["input_errors"]["alternative"]
    assert "two, upper, lower" in res["input_errors"]["alternative"]


def test_a_varied_input_without_a_column_gets_one():
    res = run("means.two.cluster", "power", {**base_inputs("means.two.cluster"), "kr": "1 2"})
    assert "kr" in [c["key"] for c in res["columns"]] and [r["kr"] for r in res["rows"]] == [1, 2]
    res = run("means.two.t", "n", {**base_inputs("means.two.t"), "dropout": "0 0.1"})
    col = next(c for c in res["columns"] if c["key"] == "dropout")
    assert col["fmt"] == "pct" and [r["dropout"] for r in res["rows"]] == [0, 0.1]


def test_mixed_model_simulations_do_not_depend_on_the_units():
    """The outcome's location and the time unit change nothing in the analysis, so nothing in the simulation:
    the Satterthwaite degrees of freedom once lost their digits to a large mean and mixed the time unit into
    their finite-difference steps."""
    base = {**base_inputs("sim.lmm.cluster"), "df": "satterthwaite", "k1": "4", "m": "20", "icc": "0.1",
            "reps": "200", "seed": "7"}
    rows = [first_row(run("sim.lmm.cluster", "power", {**base, "mu1": f"{c + 0.8}", "mu2": f"{c}"}, threads=1))
            for c in (0, 310, 10000)]
    assert len({r["power"] for r in rows}) == 1 and max(r["mean_df"] for r in rows) - min(r["mean_df"] for r in rows) < 1e-5
    base = {**base_inputs("sim.lmm.longitudinal"), "df": "satterthwaite", "n1": "12", "re": "slopes", "sd1": "0.5",
            "delta": "0.5", "m": "5", "duration": "4", "reps": "200", "seed": "7"}
    years = first_row(run("sim.lmm.longitudinal", "power", base, threads=1))
    days = first_row(run("sim.lmm.longitudinal", "power", {**base, "delta": repr(0.5 / 365.25), "sd1": repr(0.5 / 365.25),
                                                           "duration": "1461"}, threads=1))
    assert years["power"] == days["power"] and years["mean_df"] == pytest.approx(days["mean_df"], rel=1e-4)


def test_satterthwaites_degrees_of_freedom_keep_their_digits_with_large_clusters():
    """With equal clusters and a positive estimate of the cluster variance the test is the t-test on the cluster
    means, with k1 + k2 − 2 degrees of freedom.  With 10⁷ subjects a cluster the REML criterion is about 10⁸, and its
    rounding swamped the differences its Hessian is taken from (the mean was 6.73 for 6); with an ICC of 10⁻⁶
    besides, the steps took the cluster variance to 0 and back, and 57 of 200 studies of 60 and 30 such clusters
    fell back to the fixed degrees of freedom."""
    base = {**base_inputs("sim.lmm.cluster"), "df": "satterthwaite", "k1": "4", "m": "1e7", "cv": "0", "icc": "0.05",
            "reps": "200", "seed": "7"}
    assert first_row(run("sim.lmm.cluster", "power", base, threads=2))["mean_df"] == pytest.approx(6, rel=1e-3)
    row = first_row(run("sim.lmm.cluster", "power", {**base, "k1": "60", "kr": "0.5", "cv": "1e-3", "icc": "1e-6",
                                                      "mu1": "0", "sigma": "3"}, threads=2))
    assert "df_fallbacks" not in row


def test_satterthwaites_degrees_of_freedom_with_few_clusters_of_millions_are_those_of_a_50_digit_fit():
    """Three clusters a group, of 10⁷ subjects on average with a CV of 5: in 5 of 400 studies the fit ended at a
    cluster variance of 10⁻¹⁰ to 10⁻⁹ σ² where the REML optimum is 0 (the criterion's values cannot tell them
    apart), and Satterthwaite's degrees of freedom fell back to k1 + k2 − 2; elsewhere they were up to 6% off, from
    finite differences of a criterion of 10⁸ and a fit that located θ² only to 6 · 10⁻⁴.  The mean of the 50-digit
    fits' degrees of freedom is 7862226.554949597 (the fallbacks had put it at 6054562)."""
    row = first_row(run("sim.lmm.cluster", "power", {**base_inputs("sim.lmm.cluster"), "k1": "3", "m": "1e7",
                                                     "cv": "5", "icc": "0.01", "mu1": "0.5", "alternative": "upper",
                                                     "reps": "400", "seed": "163117"}, threads=2))
    assert "df_fallbacks" not in row and row["mean_df"] == pytest.approx(7862226.554949597, rel=1e-6), row


def test_clusters_of_one_subject_take_the_degrees_of_freedom_of_the_clusters():
    """With one subject a cluster the cluster and residual variances are not told apart (the criterion is the same for
    every split of their sum), and 159 of 400 studies used k1 + k2 − 2 degrees of freedom with a sentence saying that
    Satterthwaite's could not be computed; the fit now gives the residual variance all of it, where Satterthwaite's
    are k1 + k2 − 2."""
    row = first_row(run("sim.lmm.cluster", "power", {**base_inputs("sim.lmm.cluster"), "k1": "200", "m": "1",
                                                     "cv": "0", "reps": "400"}, threads=2))
    assert "df_fallbacks" not in row and row["mean_df"] == pytest.approx(398, rel=1e-6), row


def simulate_and_approximate(pid: str, inputs: Dict[str, str]):
    sim = first_row(ssc.run(pid, "power", inputs, threads=2, power_curve=False))
    approx = first_row(ssc.run(pid, "power", inputs, power_curve=False, approximate=True))
    return sim, approx


def test_mixed_model_simulations_agree_with_their_approximation_when_one_input_dwarfs_the_others():
    """A follow-up of 10⁶ time units, an intercept SD 10⁵ times the residual SD, a difference of 10¹⁶⁰ or a σ of
    10¹⁶⁰ once made each simulated study's sums of squares lose every digit of its residual variation, or
    overflow: the simulations reported a power near 0 (the approximation near 1, or near α), from fits that
    looked fine or that had all failed.  The approximations of the cluster and stepped-wedge designs squared σ,
    and gave α for a difference of one σ when σ was 10¹⁶⁰.  Values past 10³⁰ are refused now; 10³⁰ (or a residual
    SD of 10⁻³⁰) dwarfs the others as much as the inputs allow."""
    reps = 400
    cases = [("sim.lmm.longitudinal", {"duration": "1000000"}), ("sim.lmm.longitudinal", {"duration": "400000"}),
             ("sim.lmm.longitudinal", {"sd0": "100000"}), ("sim.lmm.longitudinal", {"delta": "1e30"}),
             ("sim.lmm.longitudinal", {"sd1": "1e30"}), ("sim.lmm.longitudinal", {"re": "intercepts", "sd0": "1e30"}),
             ("sim.lmm.longitudinal", {"sde": "1e-30"}),
             ("sim.lmm.longitudinal", {"timing": "custom", "times": "1000000 1000001 1000003"}),
             ("sim.lmm.cluster", {"mu1": "1e30"}), ("sim.lmm.cluster", {"sigma": "1e30"}),
             ("sim.lmm.cluster", {"sigma": "1e30", "mu1": "1e30"}), ("sim.lmm.cluster", {"mu1": "1e12"}),
             ("sim.sw", {"delta": "1e30"}), ("sim.sw", {"sigma": "1e30"}),
             ("sim.sw", {"sigma": "1e30", "delta": "1e30"}),
             ("sim.sw", {"icc": "0.999999", "analysis": "hooper", "cac": "0.5"}), ("sim.mmrm", {"delta": "1e10"}),
             ("sim.mmrm", {"sd": "1e30"})]
    for pid, extra in cases:
        sim, approx = simulate_and_approximate(pid, {**extra, "reps": str(reps)})
        assert sim is not None and approx is not None, (pid, extra)
        assert sim["failures"] == 0, (pid, extra)
        p = approx["power"]
        tolerance = 4 * math.sqrt(max(p * (1 - p), 0.01) / reps) + 0.03  # the approximation's own bias is small
        assert sim["power"] == pytest.approx(p, abs=tolerance), (pid, extra)


def test_mixed_model_simulations_ignore_what_the_analysis_ignores():
    """With complete data the test of the slopes is the two-sample t-test on the subjects' least-squares slopes,
    whatever the intercepts' SD: the simulated power once fell from 98% to 42% at an SD of 10⁵ and came back at
    10⁶.  Two visits 10⁶ time units after the origin are like any two visits one time unit apart: the simulated
    power there was once 1.8% (38% at times 0 and 1).  And beside a σ of 10²⁹ or 10³⁰ (10¹⁶⁰ once, now refused),
    μ1 − μ2 = 0.5 is nothing: the power is that of equal means."""
    base = {"reps": "400"}
    rows = [first_row(run("sim.lmm.longitudinal", "power", {**base, "sd0": sd0}, threads=2))
            for sd0 in ("1", "1e5", "1e6", "1e30")]
    assert len({r["power"] for r in rows}) == 1
    assert max(r["mean_delta"] for r in rows) - min(r["mean_delta"] for r in rows) < 1e-9
    near, far = (first_row(run("sim.lmm.longitudinal", "power", {**base, "timing": "custom", "times": times},
                               threads=2)) for times in ("0 1", "1000000 1000001"))
    assert far["power"] == near["power"] and far["approx_power"] == pytest.approx(near["approx_power"], rel=1e-9)
    null = first_row(run("sim.lmm.cluster", "power", {**base, "mu1": "0"}, threads=2))
    for sigma in ("1e29", "1e30"):
        assert first_row(run("sim.lmm.cluster", "power", {**base, "sigma": sigma}, threads=2))["power"] == null["power"]


def test_a_simulation_whose_studies_fail_by_design_still_has_a_power():
    """Studies whose analysis fails count as not rejecting.  Where the design makes them fail (an MMRM with nearly
    every patient of a group gone, a win ratio with no wins or no losses, a GLMM with a rare outcome in small
    clusters), more than half of them failing once turned the row into an error that said to check for an extreme
    value: the power is given with the count, and refused only when no study could be analyzed."""
    for pid, inputs in (("sim.mmrm", {"drop1": "0.994", "visits": "2", "n1": "10", "reps": "200"}),
                        ("sim.winratio", {"test": "wr", "n1": "2", "reps": "200"}),
                        ("sim.glmm.cluster", {"p1": "0.3", "p2": "0.001", "m": "5", "reps": "200"})):
        res = run(pid, "power", inputs, threads=2)
        row = first_row(res)
        assert row is not None and row["failures"] > 100, (pid, res["rows"][0])
        assert f"{row['failures']:.0f} simulated studies could not be analyzed and count as not rejecting" in \
            row["_statement"], row["_statement"]
    res = run("sim.winratio", "power", {"test": "wr", "n1": "10", "med_d": "600", "med_h": "240", "follow": "1",
                                        "accrual": "1", "reps": "200"}, threads=2)
    assert res["rows"][0].get("_error") == ("the simulated power cannot be computed: none of the 200 simulated "
                                            "studies could be analyzed (each had no wins or no losses, so no finite "
                                            "win ratio)")


def test_a_simulated_study_that_cannot_be_analyzed_says_why():
    """Each failed analysis tells why, and the failures the design brings are counted apart from those of the
    arithmetic.  An MMRM whose visits were correlated within 1e-12 of 1 failed in its arithmetic in most studies,
    and was given a power of 4% to 17% where it is near 80%, or refused because "the model needs patients of both
    groups measured at every visit", though every patient was; an ordinal outcome all in one category was said to
    lie beyond the other group's; and trials in which every pair tied counted without a word."""
    row = first_row(run("sim.mmrm", "power", {"rho": "0.9999999999995", "reps": "400"}, threads=2))
    assert row["failures"] == 0
    assert row["power"] == pytest.approx(row["approx_power"], abs=4 * math.sqrt(0.8 * 0.2 / 400))
    expected = [("sim.mmrm", {"drop1": "0.994", "visits": "2", "n1": "10"}, "a group had no patients left at a visit"),
                ("sim.mmrm", {"drop1": "0.994", "visits": "2", "n1": "10", "analysis": "completers"},
                 "fewer than 2 patients of a group completed the study"),
                ("sim.winratio", {"n1": "10", "med_d": "600", "med_h": "240", "follow": "1", "accrual": "1"},
                 "every pair of patients tied"),
                ("sim.glmm.cluster", {"p1": "0.999", "p2": "0.99", "m": "10", "k1": "5"},
                 "no subject of a group had the outcome, or every one had it")]
    for pid, inputs, why in expected:
        row = first_row(run(pid, "power", {**inputs, "reps": "200"}, threads=2))
        assert row is not None and row["failures"] > 100 and why in row["failures_why"], (pid, inputs)
        assert f"could not be analyzed and count as not rejecting ({row['failures_why']})." in row["_statement"]
    for test in ("wald", "mw"):
        res = run("sim.ordinal", "power", {"p2": "0.999999 0.000001", "n1": "5", "test": test, "reps": "200"},
                  threads=2)
        assert res["rows"][0]["_error"] == ("the simulated power cannot be computed: none of the 200 simulated studies "
                                            "could be analyzed (every outcome fell in one category)")
    res = run("sim.winratio", "power", {"test": "fs", "tier2": "none", "med_d": "1e9", "follow": "1", "n1": "20",
                                        "reps": "200"}, threads=2)
    assert res["rows"][0]["_error"].endswith("could be analyzed (every pair of patients tied)")


def test_the_other_tests_of_a_simulated_study_count_where_its_own_analysis_failed():
    """The tests run alongside on the same simulated studies counted as not rejecting wherever the analysis that
    decides the power failed: with the win ratio deciding, 5 patients a group and hazard ratios of 0.1, the
    Finkelstein–Schoenfeld test had 2.75% on the very trials on which it has 28.75% when it decides."""
    def pair(pid, key, inputs, first, second, power_first, power_second):
        a = first_row(run(pid, "power", {**inputs, key: first}, threads=2))
        b = first_row(run(pid, "power", {**inputs, key: second}, threads=2))
        assert a["failures"] > 0, (pid, inputs)
        assert a[power_second] == b["power"] and b[power_first] == a["power"], (pid, inputs, a, b)
    pair("sim.winratio", "test", {"n1": "5", "hr_d": "0.1", "hr_h": "0.1", "reps": "400"}, "wr", "fs",
         "power_wr", "power_fs")
    pair("sim.mmrm", "analysis", {"drop1": "0.994", "visits": "2", "n1": "10", "reps": "200"}, "mmrm", "completers",
         "power_mmrm", "power_completers")
    pair("sim.ordinal", "test", {"or": "20", "n1": "10", "reps": "400"}, "wald", "lr", "power_wald", "power_lr")


def test_the_statement_says_where_the_other_tests_of_a_simulated_study_failed():
    """A test run alongside that could not be computed in all, or nearly all, of the studies still read as a plain
    0.0%: the Wald test when every study was separated, and the completers' t-test when 198 of 200 trials had fewer
    than 2 completers in a group.  The failures of a first visit were said to leave too few patients to estimate its
    covariances with the earlier visits, which it does not have."""
    row = first_row(run("sim.ordinal", "power", {"or": "1e30", "test": "lr", "reps": "200"}, threads=2))
    assert row["power_wald"] is None and row["power_wald_failures"] == 200 and row["failures"] == 0
    assert ("None of the 200 simulated studies could be analyzed by the Wald test (one group's outcomes all lay at or "
            "above the other's") in row["_statement"]
    assert "the Wald test had a power" not in row["_statement"]
    row = first_row(run("sim.mmrm", "power", {"drop1": "0.994", "visits": "2", "n1": "10", "reps": "200"},
                        threads=2))
    assert row["power_completers_failures"] == 198
    assert ("198 of the 200 simulated trials could not be analyzed by the completers' t-test and count as not "
            "rejecting (fewer than 2 patients of a group completed the study).") in row["_statement"]
    assert "power_mmrm_failures" not in row  # the MMRM decides: its failures are the row's own
    row = first_row(run("sim.mmrm", "power", {"visits": "1", "n1": "3", "drop1": "0.5", "drop2": "0.5",
                                              "reps": "200"}, threads=2))
    assert "only one patient of each group was measured at the first visit" in row["failures_why"]
    assert "earlier visits" not in row["failures_why"]


def test_a_cluster_of_one_subject_is_analyzed_by_logistic_regression():
    """With clusters of one subject the random intercept's SD is not identified; each fit 'succeeded' with the
    standard error of a ridge of maxima, and the simulated power was 1.5% for an approximation of 96%."""
    for inputs, approx in (({"p1": "0.3", "p2": "0.2", "k1": "500", "icc": "0"}, 0.957),
                           ({"p1": "0.5", "p2": "0.2", "k1": "100", "icc": "0.05"}, 0.997),
                           ({"p1": "0.3", "p2": "0.3", "k1": "100", "icc": "0.05"}, 0.05)):
        row = first_row(run("sim.glmm.cluster", "power", {**inputs, "m": "1", "test": "z", "reps": "400"}, threads=2))
        assert row["approx_power"] == pytest.approx(approx, abs=5e-4)
        assert row["power"] == pytest.approx(approx, abs=4 * math.sqrt(approx * (1 - approx) / 400) + 0.005), inputs


def test_the_ordinal_approximation_gives_the_share_of_separated_studies():
    """The approximation gave 100% for the Wald test where nearly every study is separated and it cannot reject;
    it now gives the expected share of separated studies, which a sum over every pair of tables confirms."""
    def separated(p1, n1, p2, n2):
        k, total = len(p1), 0.0
        tables = lambda n: [c for c in itertools.product(range(n + 1), repeat=k) if sum(c) == n]
        prob = lambda c, p: math.factorial(sum(c)) * math.prod(q ** x / math.factorial(x) for q, x in zip(p, c))
        for c1 in tables(n1):
            for c2 in tables(n2):
                used = [j for j in range(k) if c1[j] + c2[j]]
                if len(used) < 2:
                    continue
                g1, g2 = [j for j in range(k) if c1[j]], [j for j in range(k) if c2[j]]
                if g1[-1] <= g2[0] or g2[-1] <= g1[0]:
                    total += prob(c1, p1) * prob(c2, p2)
        return total
    for p2, odds, n1, r in (("0.5 0.5", 1.8, 5, 1), ("0.2 0.5 0.3", 4, 4, 1.5), ("0.1 0.3 0.4 0.2", 0.2, 3, 2)):
        res = ssc.run("sim.ordinal", "power", {"p2": p2, "or": str(odds), "n1": str(n1), "r": str(r)},
                      power_curve=False, approximate=True)
        q2 = [float(x) for x in p2.split()]
        cum, below, q1 = 0.0, 0.0, []
        for j, q in enumerate(q2):
            cum += q
            up = 1.0 if j + 1 == len(q2) else 1 / (1 + math.exp(-(math.log(cum / (1 - cum)) - math.log(odds))))
            q1.append(up - below)
            below = up
        expected = separated(q1, n1, q2, math.ceil(r * n1 - 1e-9))
        assert res["rows"][0]["separated"] == pytest.approx(expected, rel=1e-10), (p2, odds)
    row = ssc.run("sim.ordinal", "power", {"or": "1e5"}, power_curve=False, approximate=True)["rows"][0]
    assert row["power"] == 1 and row["separated"] > 0.99
    assert "In an expected 99.1% of such studies one group's outcomes would all lie at or above" in row["_statement"]
    assert "separated" not in ssc.run("sim.ordinal", "power", {}, power_curve=False, approximate=True)["rows"][0]


def test_an_ordinal_outcome_whose_groups_separate_says_so_both_ways():
    """With odds ratios from 1e5 to 1e30, or 1e-10 to 1e-5, nearly every simulated study has one group at or above
    the other, where the Wald test cannot reject: the power was 0 against an approximation of 1, with no failure
    and no word (only odds ratios of 1e-20 and below, whose group 1 rounds to one category, were refused).  Such
    studies are now failed analyses of the Wald test, and the row says why; the likelihood-ratio test is
    unaffected.  The category probabilities are limited to 0 to 1."""
    for odds in ("1e10", "1e30", "1e-5", "1e-10"):
        res = run("sim.ordinal", "power", {"or": odds, "reps": "200"}, threads=2)
        assert "could be analyzed (one group's outcomes all lay at or above the other's" in res["rows"][0]["_error"]
    row = first_row(run("sim.ordinal", "power", {"or": "1e5", "reps": "200"}, threads=2))
    assert row["failures"] > 190 and row["separated"] == row["failures"] / 200
    assert "one group's outcomes all lay at or above the other's" in row["_statement"]
    row = first_row(run("sim.ordinal", "power", {"or": "1e30", "test": "lr", "reps": "200"}, threads=2))
    assert row["power"] == 1 and row["failures"] == 0
    for text, message in (("-1e-30 0.5 0.5", "values must be at least 0"), ("0.5 1e-31 0.5", "values must be 0 or "
                          "at least 1e-30"), ("1.5 0.5", "values must be at most 1")):
        res = run("sim.ordinal", "power", {"p2": text, "reps": "200"})
        assert res["input_errors"]["p2"] == message
        res = run("sim.ordinal", "power", {"effect": "probs", "p1": text, "p2": "0.5 0.5", "reps": "200"})
        assert res["input_errors"]["p1"] == message


def test_visits_close_together_far_from_time_zero_keep_what_they_tell():
    """Two visits one time unit apart, long before a last visit at time 0: from time 0, their Z'Z looked singular,
    and the analysis dropped the random slope's spread in their difference (and the simulation left it out of the
    data), so that both agreed on 95.5% where the exact GLS power is 70.5%.  The approximation is now that power
    (to the z test) whatever the distance, and the simulation the same as a thousand units before time 0."""
    exact = 0.705418  # 50-digit GLS: each subject's slope shows in its level far from time 0
    base = {"timing": "custom", "sd1": "10", "dropout": "0.2", "df": "z", "reps": "400"}
    near = {delta: first_row(run("sim.lmm.longitudinal", "power", {**base, "times": "-1001 -1000 0", "delta": delta},
                                 threads=2)) for delta in ("5", "0")}
    for distance in (1e6, 1e9, 1e15):
        times = f"{-(distance + 1):.17g} {-distance:.17g} 0"
        for delta in ("5", "0"):
            inputs = {**base, "times": times, "delta": delta}
            approx = first_row(ssc.run("sim.lmm.longitudinal", "power", inputs, power_curve=False, approximate=True))
            sim = first_row(run("sim.lmm.longitudinal", "power", inputs, threads=2))
            assert approx["power"] == pytest.approx(exact if delta == "5" else 0.05, abs=1e-6), (times, delta)
            assert sim["power"] == near[delta]["power"], (times, delta)


def test_robust_and_multiarm_simulations_do_not_depend_on_where_the_means_lie():
    """Each group is simulated around the reference group's mean (each arm around the control's): means of 10¹⁶
    once swamped an SD of 2, which gave the t-tests a type I error of 35% and the multi-arm comparisons none at
    all.  Means 10¹⁶ away give the power of means at 0, study for study."""
    robust = {**base_inputs("sim.means.robust"), "dist": "normal", "sd1": "2", "sd2": "2", "n1": "30", "reps": "400"}
    for diff in (0, 2):
        near = first_row(run("sim.means.robust", "power", {**robust, "mu1": str(diff), "mu2": "0"}, threads=2))
        far = first_row(run("sim.means.robust", "power", {**robust, "mu1": repr(1e16 + diff), "mu2": "1e16"},
                            threads=2))
        assert [far[k] for k in ("power_t", "power_welch", "power_mw")] == \
            [near[k] for k in ("power_t", "power_welch", "power_mw")], diff
    multi = {**base_inputs("sim.multiarm"), "sigma": "2", "reps": "1000"}
    for mu in ((0, 0, 0), (2, 0, 0)):
        near = first_row(run("sim.multiarm", "power", {**multi, "mu0": "0", "mu": " ".join(map(str, mu))}, threads=2))
        far = first_row(run("sim.multiarm", "power", {**multi, "mu0": "1e16",
                                                      "mu": " ".join(repr(1e16 + m) for m in mu)}, threads=2))
        assert far["power"] == near["power"] and far["power_all"] == near["power_all"], mu


def test_satterthwaite_degrees_of_freedom_with_a_residual_sd_far_below_the_random_slopes():
    """With random slopes 10⁸⁰ residual SDs apart over the follow-up, the squares of the finite-difference steps
    passed below the smallest double, and every study fell back to N − 2 degrees of freedom without a word."""
    base = {"n1": "10", "dropout": "0.5", "sd1": "5e29", "duration": "1e30", "delta": "1e30", "reps": "200"}
    rows = [first_row(run("sim.lmm.longitudinal", "power", {**base, "sde": sde}, threads=2))
            for sde in ("1e-20", "1e-30")]
    assert rows[0]["mean_df"] == pytest.approx(rows[1]["mean_df"], rel=1e-5) and rows[0]["mean_df"] < 17, rows


def test_random_slopes_that_dwarf_the_intercepts_away_from_time_zero_give_the_exact_power():
    """With slopes 10⁶ to 10¹⁵ times as variable as the intercepts at time 0 and visits away from it (at 1, 2 and 4,
    in calendar years, 10¹⁰ time units on), a subject's level at the visits all but follows from its slope: measured
    from the first visit, X'V⁻¹X is singular to within 10⁻¹³, which once gave 100% for an exact GLS power of 8.6%, and
    then refused these inputs (the simulated power too) or kept three digits.  Time is now measured from where the
    estimates of the mean level and slope are uncorrelated: the approximation is the exact power (300-digit GLS),
    and the simulation agrees with it within its Monte Carlo error, under the alternative and under the null."""
    cases = [("1 2 4", "1e8", "0.5", "0.3", "10", "2.5e7", 0.0864997134942078),
             ("1 2 4", "1e15", "0.5", "0.3", "10", "2.5e14", 0.0864997134942078),
             ("1 2 4", "1e8", "0", "0", "100", "28284271.2474619", 0.516005273976175),
             ("10000000000 10000000001 10000000003", "1e6", "0.5", "0.3", "10", "250000", 0.0864997134941966)]
    for times, sd1, rho, dropout, n1, delta, exact in cases:
        inputs = {"timing": "custom", "times": times, "sd1": sd1, "rho": rho, "dropout": dropout, "n1": n1,
                  "delta": delta, "df": "z"}
        approx = first_row(ssc.run("sim.lmm.longitudinal", "power", inputs, power_curve=False, approximate=True))
        assert approx["power"] == pytest.approx(exact, abs=1e-9), inputs
    reps = 1000
    base = {"timing": "custom", "times": "1 2 4", "sd1": "1e7", "rho": "0.5", "dropout": "0.3", "n1": "100",
            "df": "z", "reps": str(reps)}
    for delta, exact in (("2828427.125", 0.516005), ("0", 0.05)):
        sim = first_row(run("sim.lmm.longitudinal", "power", {**base, "delta": delta}, threads=2))
        assert sim["failures"] == 0 and sim["approx_power"] == pytest.approx(exact, abs=1e-6)
        assert abs(sim["power"] - exact) < 4 * math.sqrt(exact * (1 - exact) / reps), (delta, sim["power"])


def test_satterthwaite_degrees_of_freedom_are_computed_where_the_random_effects_all_but_coincide():
    """Subjects whose lines all but cross away from the time origin (an intercept SD of 100 and slopes correlated
    0.99 with it, or visits at 10 to 13 with the lines crossing near 0), and slopes 10⁴⁷ times as variable as the
    intercept: the Hessian of the REML criterion came out not positive definite (finite-difference steps near where
    D stops being a covariance matrix, fits on its boundary, a correlation the fit never moved), and every study, or
    a share of them, used N − 2 degrees of freedom while the statement spoke of a Satterthwaite t-test.  Where a
    study still cannot have its own, the row and the statement say how many."""
    base = {"timing": "custom", "reps": "200"}
    for extra, n in (({"times": "0 1 2 3", "sd0": "100", "sd1": "10", "rho": "0.99", "dropout": "0.3"}, 100),
                     ({"times": "10 11 12 13", "sd0": "1", "sd1": "10", "rho": "0.5"}, 20),
                     ({"timing": "equal", "sd1": "5e29", "delta": "6e29", "sde": "1e-30", "dropout": "0.5",
                       "duration": "1e18"}, 6)):
        row = first_row(run("sim.lmm.longitudinal", "power", {**base, **extra, "n1": str(n)}, threads=2))
        assert "df_fallbacks" not in row and row["mean_df"] != 2 * n - 2, (extra, row)
        assert "could not be computed" not in row["_statement"]
    # two visits: the residual and random-effect variances are not told apart, only the covariance of the two
    # measurements, and the Hessian in all four was singular (69 of these studies used N − 2); the fit now holds D
    # at rank 1, which reaches each covariance once
    row = first_row(run("sim.lmm.longitudinal", "power", {**base, "times": "0 1", "n1": "60"}, threads=2))
    assert "df_fallbacks" not in row and row["mean_df"] == pytest.approx(118, abs=1e-3)
    # lines that all but coincide, and almost every subject seen once: in 47 of these studies no subject was measured
    # three times, the ridge of two visits, where 21 used N − 2 and the others took degrees of freedom from anywhere
    # on it; D is held at rank 1 there too
    inputs = {"timing": "equal", "m": "5", "duration": "4", "dropout": "0.999812", "rho": "0.999563", "n1": "50",
              "seed": "10130", "reps": "200"}
    row = first_row(run("sim.lmm.longitudinal", "power", inputs, threads=2))
    assert "df_fallbacks" not in row and "could not be computed" not in row["_statement"], row
    # five subjects a group and 90% dropout: a study with a handful of subjects measured three times can have a
    # Hessian that is not positive definite either way, and the row and the statement say how many used N − 2
    inputs = {"timing": "custom", "times": "0 0.590874 0.815691 1.25445 2.06801", "sd1": "6.372", "sde": "1.175e-06",
              "dropout": "0.899", "n1": "5", "reps": "200"}
    row = first_row(run("sim.lmm.longitudinal", "power", inputs, threads=2))
    k = row.get("df_fallbacks", 0)
    assert 0 < k < 30, row
    assert (f"In {k} of the 200 simulated studies, Satterthwaite's degrees of freedom could not be computed, and the "
            "test used N − 2.") in row["_statement"]


def test_satterthwaites_degrees_of_freedom_are_exact_where_theory_gives_them():
    """Complete, balanced data whose REML fit is interior: the test is the two-sample t-test on the subjects' slopes,
    with N − 2 degrees of freedom.  With visits in calendar years and slopes 10⁶ times the residual SD, Nelder–Mead
    stopped at saddles 0.05 to 0.9 above the optimum, where the degrees of freedom were 9.05 to 10.7 for 10 (or
    N − 2 for want of any); with intercepts and slopes correlated 0.8, the finite differences in D averaged 37.9946
    for 38.  The derivatives are now exact, and a fit that is not at a minimum is taken on from there."""
    row = first_row(run("sim.lmm.longitudinal", "power", {"timing": "custom", "times": "2020 2021 2022 2023 2024",
                                                          "sd1": "1e6", "n1": "6", "reps": "200"}, threads=2))
    assert "df_fallbacks" not in row and row["mean_df"] == pytest.approx(10, rel=1e-5), row
    row = first_row(run("sim.lmm.longitudinal", "power", {"m": "10", "duration": "9", "sd0": "3", "sd1": "1",
                                                          "sde": "0.3", "rho": "0.8", "n1": "20", "reps": "200"},
                        threads=2))
    assert row["mean_df"] == pytest.approx(38, rel=1e-6), row


def test_a_simulated_study_that_cannot_tell_a_groups_slope_apart_counts_as_failed():
    """With two visits and heavy dropout, every subject of a group is sometimes seen only once: nothing then tells
    that group's slope apart, and the study reported an arbitrary estimate (exactly 0 or 1 in the analysis units,
    which moved the mean estimate by 43% of itself over 100 studies) and could even reject.  Such a study now
    could not be analyzed, and counts as not rejecting; the others agree with a 40-digit REML fit."""
    inputs = {"timing": "custom", "times": "-620.357 -426.961", "sd0": "0.873", "sd1": "6.37", "rho": "0.09",
              "sde": "6.6", "dropout": "0.6", "n1": "10", "delta": "0.117", "df": "z", "reps": "100", "r": "1.5"}
    row = first_row(run("sim.lmm.longitudinal", "power", inputs, threads=2))
    assert row["failures"] == 3  # replicates 27, 55 and 93: none of the 10 subjects of group 1 was measured twice
    assert row["mean_delta"] == pytest.approx(-0.0315106, abs=1e-6)


def test_studies_that_leave_the_model_without_an_estimate_count_as_not_rejecting():
    """With heavy dropout and few subjects, a group often has no subject measured twice, or (with two visits) the
    subjects measured twice are too few to determine the covariance of the two measurements: such a study counts as
    not rejecting, as the statement says, and the power is given, where it was withheld with 'check for an extreme
    value' once more than half of the studies were so (146 of 200 here).  When all are, the message says why."""
    row = first_row(run("sim.lmm.longitudinal", "power", {"timing": "custom", "times": "0 1", "dropout": "0.8",
                                                          "n1": "5", "reps": "200"}, threads=2))
    assert row["failures"] == 177 and 0 <= row["power"] < 0.1, row
    assert "177 simulated studies could not be analyzed and count as not rejecting." in row["_statement"]
    res = run("sim.lmm.longitudinal", "power", {"timing": "custom", "times": "0 1", "dropout": "0.99", "n1": "3",
                                                "reps": "100"}, threads=2)
    assert res["rows"][0]["_error"] == ("the simulated power cannot be computed: none of the 100 simulated studies "
                                        "could be analyzed (each had a group with too few subjects measured twice)")


def test_random_slopes_whose_lines_cross_far_from_the_visits_are_analyzed():
    """Intercepts 10¹⁸ times as variable as the slopes and correlated 1 − 10⁻¹³ with them: the subjects' lines cross
    with the least spread 10¹⁸ time units before the visits, and the time origin of the analysis, meant to move there,
    went back to 0 (measured from that far, the visits' times lost their spacing).  X'V⁻¹X then kept too few digits:
    the approximation was refused or up to 1.2 · 10⁻³ off, and 95 of 200 simulated studies failed, which put the power
    at 0.5 for 0.95.  The approximation is now the exact GLS power (400 digits) in any units, with dropout or calendar
    years, and every simulated study is analyzed, with N − 2 degrees of freedom as theory gives for complete data."""
    base = {"timing": "custom", "times": "0 1 2", "sd0": "1e18", "sd1": "1", "sde": "1e-10", "rho": "0.9999999999999",
            "n1": "30", "delta": "0.6454972243679028", "df": "z"}
    for extra in ({}, {"sd0": "1e7", "sd1": "1e-10", "sde": "1e-20", "delta": "6.454972243679028e-11"},
                  {"dropout": "0.3"}, {"times": "2020 2021 2022"}):
        approx = first_row(ssc.run("sim.lmm.longitudinal", "power", {**base, **extra}, power_curve=False,
                                   approximate=True))
        assert approx["power"] == pytest.approx(0.70541800111380, abs=1e-9), extra
    sim = first_row(run("sim.lmm.longitudinal", "power", {**base, "df": "satterthwaite", "reps": "400"}, threads=2))
    assert sim["failures"] == 0 and "df_fallbacks" not in sim and sim["mean_df"] == pytest.approx(58, abs=1e-3), sim
    p = sim["approx_power"]
    assert abs(sim["power"] - p) < 4 * math.sqrt(p * (1 - p) / 400), sim
    # two visits that span time 0, and lines that cross 4 · 10²² time units after them
    span = {"timing": "custom", "times": "-52856.7 0.0263229", "sd0": "1.903e+22", "sd1": "0.5", "sde": "9.6091e-20",
            "rho": "-0.999999999999999", "dropout": "0.5", "n1": "100", "r": "0.5", "delta": "0.21650635094610987",
            "df": "z"}
    approx = first_row(ssc.run("sim.lmm.longitudinal", "power", span, power_curve=False, approximate=True))
    assert approx["power"] == pytest.approx(0.70541800111380, abs=1e-9)


def test_a_fit_that_stops_short_of_a_singular_random_effects_covariance_takes_it():
    """Visits spanning 400 time units with 60% dropout: in one of 100 studies the subjects' lines all but cross at one
    time, and the fit stopped where the criterion is 3 · 10⁻¹⁰ of itself above its optimum on the boundary, with
    Satterthwaite's degrees of freedom 221 for 148; with three visits, in one of 100 studies the fit stopped 0.04 above
    it, at intercepts and slopes correlated 0.9992 (22.78 for 22.86).  The means of 100 studies' degrees of freedom
    from 80- and 60-digit fits are 147.94767275595697 and 22.926241777233322."""
    base = {"timing": "custom", "df": "satterthwaite", "reps": "100", "r": "1.5"}
    for inputs, df in (({"times": "-152.877 -131.082 34.576 35.931 246.46", "sd0": "0.213", "sd1": "0.42",
                         "rho": "-0.13", "sde": "5.97", "dropout": "0.6", "n1": "60", "delta": "0.0236"},
                        147.94767275595697),
                       ({"times": "0.876 1.304 2.696", "sd0": "0.632", "sd1": "8.23", "rho": "0.61", "sde": "0.112",
                         "dropout": "0.6", "n1": "10", "delta": "0.116"}, 22.926241777233322)):
        row = first_row(run("sim.lmm.longitudinal", "power", {**base, **inputs}, threads=2))
        assert "df_fallbacks" not in row and row["mean_df"] == pytest.approx(df, rel=1e-6), row


def test_the_win_ratio_approximation_is_alpha_with_no_effect():
    """Identical groups win as often as they lose; the pair probabilities, estimated from random pairs, once
    carried a fixed Monte Carlo offset in favor of treatment."""
    base = {**base_inputs("sim.winratio"), "hr_d": "1", "hr_h": "1", "tier3": "none", "n1": "300"}
    for side in ("upper", "lower", "two"):
        row = first_row(ssc.run("sim.winratio", "power", {**base, "alternative": side}, approximate=True))
        assert row["power"] == pytest.approx(0.05, abs=1e-12), side
    a = first_row(ssc.run("sim.winratio", "power", {**base, "hr_d": "0.8", "hr_h": "0.75", "alternative": "upper",
                                                    "med_d": "60", "med_h": "24"}, approximate=True))
    b = first_row(ssc.run("sim.winratio", "power", {**base, "hr_d": "1.25", "hr_h": repr(1 / 0.75), "alternative": "lower",
                                                    "med_d": "75", "med_h": "32"}, approximate=True))
    assert a["power"] == pytest.approx(b["power"], abs=1e-12)


@pytest.mark.parametrize("pid,inputs,key,first", [
    ("prop.one.ci", {"method": "wilsoncc", "sides": "lower", "p": "0.05", "width": "0.01"}, "n", 1137),
    ("prop.one.ci", {"method": "wilsoncc", "sides": "lower", "p": "0.01", "width": "0.005"}, "n", 725),
    ("prop.two.ci", {"method": "ac", "sides": "upper", "p1": "0.9", "p2": "0.05", "r": "1", "width": "0.1"}, "n1", 8),
    ("prop.two.ci", {"method": "newcombecc", "sides": "upper", "p1": "0.9", "p2": "0.05", "r": "2", "width": "0.06"},
     "n1", 59),
])
def test_an_interval_stays_within_the_target_width_from_the_reported_size(pid, inputs, key, first):
    """At a few degenerate small sizes a bound falls on the wrong side of the estimate, or the width dips under
    the target before it rises again: the answer is the size from which the width stays within the target."""
    row = first_row(run(pid, "n", {**base_inputs(pid), **inputs}))
    assert row[key] == first and 0 < row["width"] <= float(inputs["width"])
    res = run(pid, "power" if "power" in modes(pid) else "width", {**base_inputs(pid), **inputs,
                                                                    key: f"{first} to {first + 300} by 1"})
    assert all(0 < r["width"] <= float(inputs["width"]) for r in res["rows"])


def test_ancova_searches_from_its_own_smallest_groups():
    row = first_row(run("means.two.ancova", "n", {**base_inputs("means.two.ancova"), "alloc": "n1", "n1": "50",
                                                  "mu1": "2.5", "mu2": "0", "r2": "0.5", "ncov": "6"}))
    assert row["n2"] == 3
    row = first_row(run("means.two.ancova", "n", {**base_inputs("means.two.ancova"), "mu1": "10", "ncov": "0"}))
    assert row["n1"] == 3 and row["power"] <= 1


def test_clearer_messages_and_statements():
    res = run("reg.slope", "b1", {**base_inputs("reg.slope"), "sdtype": "sy", "sigma": "1.5", "sx": "3", "b0": "0.6",
                                  "n": "62", "alternative": "lower"})
    assert "B0·σx must be smaller than σy" in res["rows"][0]["_error"]
    res = run("means.two.t", "n", {**base_inputs("means.two.t"), "alloc": "n1", "n1": "5", "mu1": "0.2"})
    assert "with N1 = 5, however large N2 is" in res["rows"][0]["_error"]
    res = run("means.two.cluster", "power", {**base_inputs("means.two.cluster"), "m": "11.5"})
    assert "whole" in res["input_errors"]["m"]
    assert "achieve " in first_row(run("diag.auc", "n", base_inputs("diag.auc")))["_statement"]
    assert "1 analysis (" in first_row(run("gs.means", "n", {**base_inputs("gs.means"), "looks": "1"}))["_statement"]
    stmt = first_row(run("means.two.ancova", "power", {**base_inputs("means.two.ancova"), "r2": "0"}))["_statement"]
    assert "instead of" not in stmt and "explaining none of the variance" in stmt
    assert "(target 80" in first_row(run("prop.two.cluster", "k", {**base_inputs("prop.two.cluster"), "power": "0.8"}))["_statement"]
    stmt = first_row(run("means.one.ci", "n", {**base_inputs("means.one.ci"), "tolprob": "0.9", "e": "0.5"}))["_statement"]
    assert "at most 0.5 with probability" in stmt and "(target 90%)" in stmt
    stmt = first_row(run("means.two.ancova", "power", {**base_inputs("means.two.ancova"), "alpha": "0.03395",
                                                       "mu1": "0.8758", "mu2": "0.1", "sigma": "0.5587", "r2": "0.363",
                                                       "n1": "50"}))["_statement"]
    assert ">99.9999%" in stmt and "100.0" not in stmt
    row = first_row(run("means.anova2", "power", {**base_inputs("means.anova2"), "means": "0 0.5 0.5 1.5", "n": "10"}))
    assert row["power"] == min(row["pa"], row["pb"], row["pab"])


def test_inputs_the_calculation_does_not_use_are_reported():
    """A program or a language model that sends n2 with equal groups must not believe it was used."""
    res = run("means.two.t", "power", {**base_inputs("means.two.t"), "n2": "99", "zzz": "1", "alloc": "equal"})
    assert "n2" in res["ignored_inputs"] and res["unknown_inputs"] == ["zzz"] and first_row(res)["n2"] == first_row(res)["n1"]
    res = run("means.two.t", "n", {**base_inputs("means.two.t"), "alloc": "sep"})
    assert res["replaced_choices"] == {"alloc": {"given": "sep", "used": "equal"}}
    res = run("means.two.t", "n", {"mu1": "0.5"})
    assert not {"ignored_inputs", "unknown_inputs", "replaced_choices"} & set(res)


@pytest.mark.parametrize("pid,solve,scale,size", [
    ("means.one.t", "mu1", "sigma", {"n": "1e8"}), ("means.two.t", "mu1", "sigma", {"n1": "1e8"}),
    ("reg.slope", "b1", "sigma", {"n": "1e8"}), ("means.two.cluster", "mu1", "sigma", {"k1": "1e7"})])
def test_an_effect_found_with_a_large_sample_scales_with_the_sd(pid, solve, scale, size):
    """The search's tolerance came from the bracket (the SD), far above the effect for a large N: at an SD of 1e-6
    and N = 1e8 the detectable mean was 2.804e-10 at 80.1% power, for 2.802e-10."""
    base = {**base_inputs(pid), **size, scale: "1"}
    unit = first_row(run(pid, solve, base))[solve]
    for sd in (1e-6, 1e-4, 1e3):
        row = first_row(run(pid, solve, {**base, scale: repr(sd)}))
        assert row[solve] / sd == pytest.approx(unit, rel=1e-10) and row["power"] == pytest.approx(0.8, abs=1e-10)


@pytest.mark.parametrize("pid,solve,key,base_key,sds", [
    ("means.one.t", "mu1", "mu1", "mu0", ["sigma"]), ("means.paired.t", "mud", "mud", "mud0", ["sd1", "sd2"]),
    ("means.two.t", "mu1", "mu1", "mu2", ["sigma"]), ("means.two.welch", "mu1", "mu1", "mu2", ["sd1", "sd2"]),
    ("means.two.repeated", "mu1", "mu1", "mu2", ["sigma"]), ("means.two.ancova", "mu1", "mu1", "mu2", ["sigma"]),
    ("means.two.cluster", "mu1", "mu1", "mu2", ["sigma"]), ("reg.slope", "b1", "b1", "b0", ["sigma"])])
def test_an_effect_lost_next_to_a_large_value_is_refused(pid, solve, key, base_key, sds):
    """μ2 = 1e16 keeps nothing of a difference of 0.9: μ1 rounded back to μ2, and the answer claimed 80% power
    (or reported 5% or 100%) for a difference of 0.  The refusal carries no power."""
    base = base_inputs(pid)
    for big in ("1e16", "-1e17", "1e30"):
        row = run(pid, solve, {**base, base_key: big})["rows"][0]
        msg = row.get("_error", "")
        assert "cannot be told apart from" in msg and "power" not in row, (big, row)
        assert "lost to rounding" in msg or "keeps fewer than four significant digits" in msg, (big, msg)
        assert "measure the values from a nearby reference instead" in msg, msg
    msg = run(pid, solve, {**base, base_key: "1e6", **{sd: "1e-30" for sd in sds}})["rows"][0].get("_error", "")
    assert "cannot be told apart from" in msg, msg
    # a large value that still carries the difference to four digits or more gives the answer it gives at 0
    zero = first_row(run(pid, solve, {**base, base_key: "0"}))
    for big in (1e6, 1e10):
        row = first_row(run(pid, solve, {**base, base_key: repr(big)}))
        assert row[key] - big == pytest.approx(zero[key], rel=1e-5), big
        assert row["power"] == pytest.approx(zero["power"], abs=1e-6)


@pytest.mark.parametrize("pid,solve,inputs,key,base_key", [
    ("means.two.t", "mu1", {"mu2": "1e9", "sigma": "1", "n1": "10000"}, "mu1", "mu2"),
    ("means.one.t", "mu1", {"mu0": "120", "sigma": "1e-9", "n": "100"}, "mu1", "mu0"),
    ("means.paired.t", "mud", {"mud0": "5e6", "n": "1e9"}, "mud", "mud0"),
    ("rates.one.poisson", "lam1", {"lam0": "1e20", "n": "50"}, "lam1", "lam0"),
    ("corr.one", "rho1", {"rho0": "0.99999999", "n": "1e9", "alternative": "upper"}, "rho1", "rho0"),
    ("corr.one", "rho1", {"rho0": "-0.99999999999", "alternative": "lower"}, "rho1", "rho0"),
    ("agree.icc", "rho1", {"rho0": "0.99999999999"}, "rho1", "rho0"),
    # the double next to the root keeps the difference, though one a search on the doubles ends at is half a
    # double off it (and the difference it took for the root missed it by as much)
    ("corr.one", "rho1", {"rho0": "-0.9999999999999", "n": "100", "alternative": "upper", "alpha": "1e-20",
                          "power": "0.5", "method": "fisher"}, "rho1", "rho0"),
    ("agree.icc", "rho1", {"k": "2", "rho0": "0.99999999999", "n": "1000", "alternative": "two", "power": "0.5"},
     "rho1", "rho0"),
    ("rates.one.poisson", "lam1", {"lam0": "1.19023e23", "n": "50"}, "lam1", "lam0")])
def test_an_effect_kept_to_its_digits_next_to_a_large_value_is_given(pid, solve, inputs, key, base_key):
    """A difference that the sum keeps to five digits or more was refused as lost to rounding (the doubles next to
    the answer lay more than a millionth of it apart); the power there is the target to the precision shown."""
    target = float({**base_inputs(pid), **inputs}["power"])
    row = first_row(run(pid, solve, {**base_inputs(pid), **inputs}))
    assert row[key] != float(inputs[base_key]) and row["power"] == pytest.approx(target, abs=5e-5), row
    again = first_row(run(pid, "power", {**base_inputs(pid), **inputs, key: repr(row[key])}))
    assert again["power"] == pytest.approx(target, abs=5e-5), again


def test_effects_next_to_a_correlation_near_one():
    """ρ1 next to ρ0 near 1: the search ended at 1 − 1e-9, below ρ0 (an answer below ρ0, or 'the power must exceed
    the significance level'), the critical values stopped at 1 − 1e-12 (4 subjects for 100% power), and the
    statement printed both as 1."""
    base = base_inputs("corr.one")
    for k in (11, 13, 15):
        # as ρ nears 1 the power depends on the distance in Fisher's z alone: that of the doubles given (those next
        # to 1 − 1e-14 are a hundredth of 1e-14 apart), repeated at ρ0 = 1 − 1e-8
        rho0 = 1 - 10.0 ** -k
        rho1 = math.tanh(math.atanh(rho0) + 0.4)
        dz = math.atanh(rho1) - math.atanh(rho0)
        row = first_row(run("corr.one", "power", {**base, "rho0": repr(rho0), "rho1": repr(rho1), "n": "50"}))
        ref = first_row(run("corr.one", "power", {**base, "rho0": repr(1 - 1e-8), "n": "50",
                                                  "rho1": repr(math.tanh(math.atanh(1 - 1e-8) + dz))}))
        assert row["power"] == pytest.approx(ref["power"], abs=1e-7), (k, row, ref)
        assert first_row(run("corr.one", "n", {**base, "rho0": repr(rho0), "rho1": repr(rho1)}))["n"] > 4
    row = first_row(run("corr.one", "rho1", {**base, "rho0": "0.999999", "n": "50"}))
    assert 0.999999 < row["rho1"] < 1 and row["power"] == pytest.approx(0.8, abs=1e-6)
    assert "ρ1 = 0.99999956 when the null correlation is ρ0 = 0.999999" in row["_statement"]
    for rho0, alt in (("0.9999999999999", "two"), ("0.9999999999999999", "upper"), ("-0.9999999999999", "lower")):
        msg = run("corr.one", "rho1", {**base, "rho0": rho0, "alternative": alt})["rows"][0].get("_error", "")
        assert msg.startswith("ρ1 cannot be told apart from ρ0"), (rho0, msg)
    msg = run("corr.one", "rho1", {**base, "rho0": "0.9999999999999999", "alternative": "upper"})["rows"][0]["_error"]
    assert msg.endswith("no number lies between it and 1"), msg
    row = first_row(run("agree.icc", "rho1", {**base_inputs("agree.icc"), "rho0": "0.999999"}))
    assert row["power"] == pytest.approx(0.8, abs=1e-9) and "of 0.9999996 when" in row["_statement"]
    msg = run("agree.icc", "rho1", {**base_inputs("agree.icc"), "rho0": "0.9999999999999"})["rows"][0]["_error"]
    assert msg.startswith("ρ1 cannot be told apart from ρ0 = 0.9999999999999"), msg


def test_a_correlation_whose_target_cannot_be_reached_says_so():
    """Near 1 a target that no ρ1 reaches was blamed on rounding: 'the difference of 9.999e-13 it needs is lost to
    rounding', where the power at the last double below 1 is 9e-27."""
    for method in ("exact", "fisher"):
        msg = run("corr.one", "rho1", {**base_inputs("corr.one"), "rho0": "0.999999999999", "n": "4", "alpha": "1e-30",
                                       "alternative": "upper", "method": method})["rows"][0]["_error"]
        assert msg == "the target power cannot be reached with this N", msg


@pytest.mark.parametrize("inputs", [
    {"k": "2", "rho0": "-0.5", "n": "2", "power": "0.999999", "alternative": "upper"},
    {"k": "2", "rho0": "0.3", "n": "2", "power": "0.999999", "alternative": "upper"},
    {"k": "3", "rho0": "0.3", "n": "2", "alpha": "1e-6", "power": "0.5", "alternative": "lower"},
    {"k": "2", "rho0": "-0.999999", "n": "3", "alpha": "1e-20", "power": "0.5", "alternative": "upper"},
    {"k": "3", "rho0": "-0.499", "n": "2", "alpha": "1e-20", "power": "0.9", "alternative": "upper"},
    {"k": "2", "rho0": "0.999999999", "n": "2", "alternative": "two"},
    {"k": "3", "rho0": "-0.49999999", "n": "30", "alternative": "lower"}])
def test_an_icc_near_its_limits_reaches_the_target_power(inputs):
    """The search ended at ρ0 + (1 − ρ0), which rounds to 1 (an ICC of 1 at 100% power), or at the singular lower
    limit, and its tolerance, relative to ρ1 − ρ0, left answers near a limit at 49% power for 50%."""
    row = first_row(run("agree.icc", "rho1", {**base_inputs("agree.icc"), **inputs}))
    k = float(inputs["k"])
    assert -1 / (k - 1) < row["rho1"] < 1 and 1 + (k - 1) * row["rho1"] > 0, row
    target = float(inputs.get("power", "0.8"))
    assert target - 1e-9 <= row["power"] <= target + min(5e-5, 0.5 * (1 - target)), row
    assert first_row(run("agree.icc", "power", {**base_inputs("agree.icc"), **inputs,
                                                "rho1": repr(row["rho1"])}))["power"] == row["power"]
    msg = run("agree.icc", "rho1", {**base_inputs("agree.icc"), "k": "4", "rho0": "0.3", "n": "2", "alpha": "1e-20",
                                    "power": "0.5", "alternative": "lower"})["rows"][0]["_error"]
    assert msg == "the target power cannot be reached for any admissible value of ρ1", msg


def test_an_icc_next_to_its_lower_limit_has_its_exact_power():
    """1 + (k − 1)ρ lost its digits next to −1/(k − 1) for k − 1 not a power of 2 (the rounding of (k − 1)ρ alone),
    and answers missed the target by more than the precision shown, refusals named wrong powers, and a reachable
    target was called unreachable.  Reference powers: the F distribution in mpmath at 50 digits."""
    base = base_inputs("agree.icc")
    for k, rho0, rho1, n, alt, ref in (("11", "-0.0999999999999", "-0.09999999999999999", "2", "lower", 0.99999906336),
                                       ("4", "-0.333333333330171", "-0.333333333333249", "3", "lower", 0.799913449556),
                                       ("7", "-0.16666666666616547", "-0.16666666666610178", "1000", "upper",
                                        0.79998557)):
        row = first_row(run("agree.icc", "power", {**base, "k": k, "rho0": rho0, "rho1": rho1, "n": n,
                                                   "alternative": alt}))
        assert row["power"] == pytest.approx(ref, abs=1e-9), (k, row["power"])
    # the target of the first is reached at the last ICC above −1/(k − 1); the second lies between two doubles
    # whose powers both miss it; the third is met within the precision shown
    row = first_row(run("agree.icc", "rho1", {**base, "k": "11", "rho0": "-0.0999999999999", "n": "2",
                                              "power": "0.999999", "alternative": "lower"}))
    assert row["rho1"] == -0.09999999999999999 and row["power"] >= 0.999999 - 5e-7, row
    msg = run("agree.icc", "rho1", {**base, "k": "4", "rho0": "-0.333333333330171", "n": "3",
                                    "alternative": "lower"})["rows"][0]["_error"]
    assert msg.startswith("the power cannot be brought to the target") and "0.799913 at -0.333333333333249" in msg
    row = first_row(run("agree.icc", "rho1", {**base, "k": "7", "rho0": "-0.16666666666616547", "n": "1000"}))
    assert row["rho1"] == -0.16666666666610178, row


def test_a_target_near_100_percent_is_met_to_the_finer_precision_shown():
    """Above 1 − 1e-9 the precision shown is finer than the 1e-9 by which a power counted as reaching the target: an
    answer could miss it by more, and a refusal read 'moves from nan at nan to …'.  Where the power itself is not
    computed to that precision (it moves up and down by 3e-12 from one double to the next) the message says so."""
    msg = run("corr.one", "rho1", {**base_inputs("corr.one"), "rho0": "0.99999999999", "n": "1e9",
                                   "power": "0.999999999999", "alternative": "lower"})["rows"][0]["_error"]
    assert "nan" not in msg and msg.startswith("the power cannot be computed to within 5e-13 of the target"), msg
    for target in ("0.9999999999", "0.99999999999"):
        row = first_row(run("agree.icc", "rho1", {**base_inputs("agree.icc"), "rho0": "0.3", "n": "1000",
                                                  "power": target}))
        assert abs(row["power"] - float(target)) <= 0.5 * (1 - float(target)), row


def test_the_exact_correlation_test_keeps_its_critical_value_at_a_tiny_alpha():
    """With ρ0 = 0 the critical value r = t/√(N − 2 + t²) rounded to 1 for a tiny alpha with few subjects, and that
    side was dropped (a power of 0 for every ρ1, and N = 5 where 4 reach the target).  Reference powers: the exact
    density of r integrated with mpmath at 50 digits."""
    base = {**base_inputs("corr.one"), "rho0": "0", "rho1": "0.9999999999999999", "alternative": "upper"}
    for n, alpha, ref in ((4, "1e-20", 0.00027015514), (4, "1e-30", 2.7021598e-14), (5, "1e-30", 3.5596341e-6),
                          (6, "1e-30", 0.99252438)):
        row = first_row(run("corr.one", "power", {**base, "n": str(n), "alpha": alpha}))
        assert row["power"] == pytest.approx(ref, rel=1e-7), (n, alpha, row["power"])
        low = first_row(run("corr.one", "power", {**base, "n": str(n), "alpha": alpha, "rho1": "-0.9999999999999999",
                                                  "alternative": "lower"}))
        assert low["power"] == row["power"]
    size = first_row(run("corr.one", "power", {**base, "rho1": "0", "n": "4", "alpha": "1e-20"}))["power"]
    assert size == pytest.approx(1e-20, rel=1e-9)
    assert first_row(run("corr.one", "n", {**base, "alpha": "1e-20", "power": "0.0002"}))["n"] == 4


@pytest.mark.parametrize("pid,solve,inputs,key", [
    ("means.one.t", "mu1", {"n": "2", "alpha": "1e-12"}, "mu1"),
    ("means.one.t", "mu1", {"n": "2", "alpha": "1e-30"}, "mu1"),
    ("means.two.t", "mu1", {"n1": "2", "alpha": "1e-30"}, "mu1"),
    ("reg.slope", "b1", {"n": "3", "alpha": "1e-30"}, "b1"),
    ("means.crossover", "delta", {"n": "4", "alpha": "1e-30"}, "delta"),
    ("means.two.ancova", "mu1", {"n1": "3", "alpha": "1e-30"}, "mu1"),
    ("means.two.cluster", "mu1", {"k1": "2", "alpha": "1e-30"}, "mu1"),
    ("means.two.slopes", "delta", {"n1": "2", "alpha": "1e-30"}, "delta"),
    ("means.anova", "f", {"n": "2", "alpha": "1e-30", "estype": "f", "k": "3"}, "f"),
    ("reg.multiple", "r2t", {"n": "10", "alpha": "1e-30"}, "r2t"),
    ("rates.one.poisson", "lam1", {"lam0": "1e-14", "n": "50"}, "lam1"),
    ("rates.one.poisson", "lam1", {"lam0": "1e-30", "n": "50", "test": "z"}, "lam1")])
def test_an_effect_search_reaches_what_the_inputs_allow(pid, solve, inputs, key):
    """Searches stopped at 1e9·σ (1e6 for f and the ANCOVA, 1e12·λ0 for a rate, an f² of 1e9) and said that no
    admissible value reaches the target, where the power mode gives 97% at μ1 = 1e12 (n = 2, α = 1e-12)."""
    row = first_row(run(pid, solve, {**base_inputs(pid), **inputs}))
    target = float(base_inputs(pid)["power"])
    assert row["power"] == pytest.approx(target, abs=1e-6), row
    again = first_row(run(pid, "power", {**base_inputs(pid), **inputs, key: repr(row[key])}))
    assert again["power"] == pytest.approx(target, abs=5e-5), again


def test_a_delta_r2_next_to_its_limit_is_found_or_refused_as_between_two_doubles():
    """The search for f² stopped at 1e15, where nine ΔR² below 1 − R²c remain, and called the target unreachable:
    with 3 subjects and α = 1e-8 the power is 0.75 at ΔR² = 0.9999999999999994 and 0.80 at the next double up."""
    base = {**base_inputs("reg.multiple"), "n": "3", "b": "1", "a": "0", "r2c": "0"}
    for alpha, before, after in (("1e-8", "0.75181 at 0.9999999999999994", "0.9999999999999996"),
                                 ("5e-9", "0.638713 at 0.9999999999999998", "0.9999999999999999")):
        msg = run("reg.multiple", "r2t", {**base, "alpha": alpha})["rows"][0]["_error"]
        assert msg.startswith("the power cannot be brought to the target at any ΔR²") and before in msg, msg
        assert msg.endswith(after), msg
    row = first_row(run("reg.multiple", "r2t", {**base, "n": "4", "alpha": "1e-12"}))
    assert row["power"] == pytest.approx(0.8, abs=5e-5), row


def test_a_confidence_level_is_shown_with_its_digits():
    for pid, conf, says in (("means.one.ci", "0.999", "99.9%"), ("means.two.ci", "0.975", "97.5%"),
                            ("corr.ci", "0.9999999999999", "99.99999999999%"), ("prop.one.ci", "0.95", "95%"),
                            ("prop.two.ci", "0.5", "50%")):
        stmt = first_row(run(pid, "n", {**base_inputs(pid), "conf": conf}))["_statement"]
        assert " " + says + " " in stmt, stmt


@pytest.mark.parametrize("pid,side", [("means.one.ci", "one"), ("means.two.ci", "one"), ("prop.one.ci", "lower"),
                                      ("prop.two.ci", "upper")])
def test_a_one_sided_interval_needs_a_confidence_above_one_half(pid, side):
    """At 50% or less a one-sided bound lies on the estimate or on the wrong side: half-widths of 0 or −2.3e14, or
    'the target width cannot be reached'."""
    for solve in ("n", "width"):
        for conf in ("0.5", "0.01", "1e-15"):
            msg = run(pid, solve, {**base_inputs(pid), "sides": side, "conf": conf})["rows"][0].get("_error", "")
            assert msg.startswith("a one-sided interval needs a confidence level above 0.5"), (solve, conf, msg)
    assert first_row(run(pid, "n", {**base_inputs(pid), "sides": side, "conf": "0.6"}))


def test_degenerate_effect_sizes_and_sizes_are_refused():
    """A total of 1 subject for a chi-square test on 2 degrees of freedom; a group the allocation leaves empty."""
    base = base_inputs("prop.chisq")
    for w in ("1e15", "1e30"):
        assert first_row(run("prop.chisq", "n", {**base, "w": w, "df": "2"}))["n"] == 3
    msg = run("prop.chisq", "power", {**base, "n": "2", "df": "2"})["rows"][0]["_error"]
    assert msg == "N must be at least 3: a chi-square test on 2 degrees of freedom compares at least 3 categories"
    # n × a tiny relative size, rounded up, is 1: the slack for rounding, 1e-9 absolute, made it 0 below 1e-9, and
    # the search for n went on until n × the size passed 1e-9 (1,000,001 at 1e-15, with a power of 100%)
    for pid, key, weights, n in (("means.anova", "ratios", "1 1 {}", 20), ("prop.trend", "weights", "1 1 1 {}", 141)):
        base = {**base_inputs(pid), "alloc": "ratio"}
        first = first_row(run(pid, "n", {**base, key: weights.format("1e-6")}))
        for w in ("1e-10", "1e-12", "1e-15", "1e-20", "1e-30"):
            row = first_row(run(pid, "n", {**base, key: weights.format(w)}))
            assert row["n"] == first["n"] and row["sizes"].endswith(", 1"), (w, row)
            row = first_row(run(pid, "power", {**base, key: weights.format(w), "n": str(n)}))
            assert row["sizes"].endswith(", 1"), (w, row)
    # every group of 1 leaves the F-test no degrees of freedom for the error
    base = {**base_inputs("means.anova"), "alloc": "ratio", "estype": "f", "k": "3", "ratios": "1e-30 1e-30 1e-30"}
    for solve, extra in (("power", {"n": "1000"}), ("f", {"n": "1000"}), ("n", {})):
        row = run("means.anova", solve, {**base, **extra})["rows"][0]
        assert row["_error"].startswith("every group has 1 subject (n × its allocation ratio, rounded up)"), row
        assert "no degrees of freedom for the error" in row["_error"] and "n" not in row and "sizes" not in row, row


def test_a_tiny_ratio_gives_its_group_one():
    """n times a ratio of 1e-10 or 1e-12, rounded up, is 1, but the slack for rounding made it 0: messages said
    '0 cluster' and '0 control', and the interval for two means took 1,001 subjects where 6 give its width."""
    msg = run("means.two.cluster", "power", {**base_inputs("means.two.cluster"), "kr": "1e-10"})["rows"][0]["_error"]
    assert msg.startswith("group 2 would have 1 cluster (k1 × the cluster ratio, rounded up)"), msg
    msg = run("diag.auc", "power", {**base_inputs("diag.auc"), "ratio": "1e-12"})["rows"][0]["_error"]
    assert msg.startswith("there would be 1 control (cases × the ratio, rounded up)"), msg
    row = first_row(run("means.two.ci", "n", {**base_inputs("means.two.ci"), "r": "1e-12", "e": "3", "sigma": "1"}))
    assert (row["n1"], row["n2"]) == (6, 1), row


def test_equal_means_far_from_zero_are_no_effect():
    """The mean of equal means of 0.1 with σ = 1e-30 (or of 1e30 with σ = 1) rounded, and left a spread that the
    ANOVA took for an effect of f = 1e13: n = 2 at 100% power.  The two-way ANOVA's additive means 0.1 0.2 0.3 0.4
    leave an interaction of 3e-17 as doubles, 3e13 times σ = 1e-30."""
    for pid, means in (("means.anova", "0.1 0.1 0.1"), ("means.anova", "1e30 1e30 1e30"),
                       ("nonpar.kruskal", "0.7 0.7 0.7"), ("means.rm.anova", "0.7 0.7 0.7")):
        sigma = "1" if means.startswith("1e30") else "1e-30"
        msg = run(pid, "n", {**base_inputs(pid), "means": means, "sigma": sigma})["rows"][0]["_error"]
        assert "means are all equal" in msg, (pid, means, msg)
    base = {**base_inputs("means.anova2"), "a": "2", "b": "3", "sigma": "1e-30"}
    msg = run("means.anova2", "n", {**base, "means": "0.1 0.1 0.1 0.1 0.1 0.1"})["rows"][0]["_error"]
    assert msg == "the cell means are all equal; there is no effect to detect", msg
    base = {**base, "b": "2", "means": "0.1 0.2 0.3 0.4"}
    msg = run("means.anova2", "n", {**base, "goal": "AB"})["rows"][0]["_error"]
    assert msg == "the cell means have no effect of the chosen kind", msg
    row = first_row(run("means.anova2", "power", {**base, "n": "10"}))
    assert row["fab"] == 0 and row["pab"] == pytest.approx(0.05) and row["pa"] > 0.99 and row["pb"] > 0.99
    # scores far from 0 that are equally spaced give the trend of 0, 1, 2, 3
    big = "1e20 1.00000000000000016384e20 1.00000000000000032768e20 1.00000000000000049152e20"
    near = first_row(run("prop.trend", "n", {**base_inputs("prop.trend"), "scores": "custom", "dvals": "0 1 2 3"}))
    far = first_row(run("prop.trend", "n", {**base_inputs("prop.trend"), "scores": "custom", "dvals": big}))
    assert far["n"] == near["n"]
    msg = run("prop.trend", "n", {**base_inputs("prop.trend"), "scores": "custom", "dvals": big,
                                  "props": "0.2 0.2 0.2 0.2"})["rows"][0]["_error"]
    assert msg == "the proportions have no linear trend across the scores", msg


def test_two_way_effects_beyond_the_rounding_of_the_means_are_kept():
    """Effects were set to none up to 8 times what the rounding of the means could make, measured at the largest
    cell: real effects of 7 to 17 doubles' spacing gave 'the cell means are all equal', and a row of 1e6 removed the
    effects of B and A × B that a row of 1 left at f = 235.7."""
    base = {**base_inputs("means.anova2"), "a": "2", "b": "2", "goal": "all"}
    row = first_row(run("means.anova2", "power", {**base, "means": "1000000 1000000 1000000 1000000.000000002",
                                                  "sigma": "1e-9", "n": "10"}))
    # the stored cells differ by 1.979e-9: f = 1.979e-9/4/1e-9 for each effect, and λ = 40 f² on (1, 36) df
    assert row["fa"] == row["fb"] == row["fab"] == pytest.approx(0.4947651177644729, rel=1e-12), row
    assert row["pab"] == pytest.approx(0.86097794658, rel=1e-9), row
    assert first_row(run("means.anova2", "n", {**base, "means": "9e29 9e29 9e29 9.00000000000001e29",
                                                "sigma": "1"}))["n"] == 2
    for shift in ("1", "1000000"):
        row = first_row(run("means.anova2", "power", {**base, "b": "3", "means": f"{shift} {shift} {shift} 0 0 1e-9",
                                                      "sigma": "1e-12", "n": "10"}))
        assert row["fb"] == row["fab"] == pytest.approx(1000 / math.sqrt(18), rel=1e-12), (shift, row)
    # cells one spacing apart differ by no more than the rounding of numbers of their size
    msg = run("means.anova2", "n", {**base, "means": "1 1 1 1.0000000000000002", "sigma": "1e-17"})["rows"][0]
    assert msg["_error"] == ("the cell means differ by no more than the rounding of numbers of their size; there is "
                             "no effect to detect"), msg


@pytest.mark.parametrize("pid,solve,inputs,says", [
    ("agree.icc", "rho1", {"k": "3", "rho0": "-0.49999999", "alternative": "lower"},
     "intraclass correlation of -0.4999999956 when the value under the null hypothesis is -0.49999999,"),
    ("corr.ci", "width", {"rho": "0.9999999", "n": "100"},
     "(0.99999985 to 0.99999993) when the sample correlation is 0.9999999."),
    ("means.one.ci", "n", {"tolprob": "0.9999999", "e": "0.5"}, "(target 99.99999%)"),
    ("means.paired.t", "mud", {"mud0": "1e9"}, "μD = 1000000000.9 (null value 1000000000)"),
    ("means.one.t", "mu1", {"mu0": "1e9"}, "μ0 = 1000000000 and the actual mean μ1 = 1000000000.66,"),
    ("means.two.t", "mu1", {"mu2": "1e9"}, "(μ1 = 1000000001.1, μ2 = 1000000000)"),
    ("reg.slope", "b1", {"b0": "1e9"}, "slope of 1000000000.4 (null value 1000000000)"),
    ("rates.one.poisson", "lam1", {"lam0": "1e9"}, "event rate of 1.000013e+09 when the rate under H0 is 1000000000,"),
    ("corr.two", "n", {"rho1": "0.9999999", "rho2": "0.9999998"}, "correlations of 0.9999999 and 0.9999998 "),
    ("means.two.cluster", "power", {"icc": "0.99999"}, "ICC = 0.99999,"),
    # a value that is not whole never shows as one, or as the limit 1 of a correlation
    ("corr.one", "rho1", {"rho0": "0.9999999", "n": "4", "alternative": "lower", "power": "0.5"},
     "when the null correlation is ρ0 = 0.9999999,"),
    ("corr.one", "power", {"rho0": "0.99990884", "rho1": "0.9999999", "alternative": "upper"},
     "ρ1 = 0.9999999 when the null correlation is ρ0 = 0.99991,"),
    ("rates.one.poisson", "lam1", {"lam0": "0.9999999", "n": "1e9"}, "when the rate under H0 is 0.9999999,"),
    ("means.one.t", "mu1", {"mu0": "0.9999999", "sigma": "0.000133"}, "null mean μ0 = 0.9999999 and")])
def test_values_side_by_side_are_shown_apart(pid, solve, inputs, says):
    """Statements showed two values that differ as one ('-0.5 when the value under the null hypothesis is -0.5',
    'μD = 1e+09 (null value 1000000000)'), or a value as the limit next to it ('(1 to 1)', 'target 100%')."""
    stmt = first_row(run(pid, solve, {**base_inputs(pid), **inputs}))["_statement"]
    assert says in stmt, stmt


@pytest.mark.parametrize("bound", ["obf", "wt_obf", "pocock", "hsd"])
@pytest.mark.parametrize("alt", ["two", "upper"])
def test_group_sequential_boundaries_at_the_smallest_alpha(bound, alt):
    """At an α of 1e-20 to 1e-30 the boundaries lost accuracy: the final nominal p exceeded α and the inflation
    factor dropped below 1."""
    for looks in ("5", "10"):
        for alpha in (1e-20, 1e-30):
            row = first_row(run("gs.means", "n", {**base_inputs("gs.means"), "boundary": bound, "alternative": alt,
                                                 "looks": looks, "alpha": repr(alpha)}))
            table = row["_tables"][0]
            keys = [c["key"] for c in table["columns"]]
            p = table["rows"][-1][keys.index("p")] * (2 if alt == "two" else 1)
            assert p <= alpha and row["inflation"] >= 1, (looks, alpha, p / alpha, row["inflation"])
