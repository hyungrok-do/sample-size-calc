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


@pytest.mark.parametrize("pid,solve,alpha", [("gs.logrank", "power", "1e-10"), ("gs.means", "power", "1e-12"),
                                             ("gs.means", "n", "1e-15"), ("gs.props", "n", "1e-15")])
def test_a_group_sequential_design_at_a_small_alpha(pid, solve, alpha):
    """The other rejection region adds less than a rounding error to the fixed design's power at a small alpha,
    which left its drift unbracketed ("the answer cannot be computed")."""
    row = first_row(run(pid, solve, {**base_inputs(pid), "alpha": alpha}))
    assert row and row["inflation"] == pytest.approx(1, abs=1e-3), row


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
    ("ratio", {"p1": "1e-20", "p2": "1e-20", "r0": "0.8"}, 0.026439090706363813)])
def test_score_tests_at_proportions_near_0_or_1_by_the_normal_approximation(scale, inputs, power):
    """The restricted estimates lost their digits near 0 and 1 (the closed-form root of the difference's cubic is
    off by 1e-16 or more whatever its size, and the smaller root of the ratio's quadratic cancelled to 0): at
    P1 = P2 = 1e-20 the power was 50% where it is about α.  The references evaluate the same formula to 120 digits."""
    row = first_row(run("prop.two.ni", "power", {"scale": scale, "method": "normal", "n1": "100", "h1": "upper",
                                                 "alpha": "0.025", **inputs}))
    assert row["power"] == pytest.approx(power, rel=1e-12), row


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
    ({"scale": "ratio", "p1": "0.3", "p2": "0.2", "n1": "30", "r0": "1e-20", "h1": "upper"}, 6e-20)])
def test_the_exact_size_is_taken_at_the_boundary_itself(inputs, size):
    """P1 was kept within 1e-9 of 0 and 1, so that the size at OR0 = 1e30 was 3.9e-9 where a 150-digit enumeration
    gives 2.6e-30."""
    row = first_row(run("prop.two.ni", "power", {"method": "exact", "alpha": "0.025", **inputs}))
    assert row["actual_alpha"] == pytest.approx(size, rel=1e-9), row


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


def test_random_slopes_that_dwarf_the_intercepts_away_from_time_zero_are_refused_not_misreported():
    """With slopes 10⁸ or 10¹⁵ times as variable as the intercepts at time 0 and visits at 1, 2 and 4, a subject's
    level at the visits all but follows from its slope, and the slopes' variance came out of rounding noise: an
    approximate power of 100% where the exact GLS power is 8.6%, and simulated powers of 63%.  Such inputs now say
    that the answer cannot be computed, and give no simulated power; where enough digits are left (slopes 10⁶ times
    as variable, visits 10¹⁰ time units from 0), the approximation is the exact power."""
    refused = "the answer cannot be computed for these inputs; check for an extreme value"
    for sd1 in (1e8, 1e15):
        inputs = {"timing": "custom", "times": "1 2 4", "sd1": f"{sd1:g}", "rho": "0.5", "dropout": "0.3",
                  "delta": f"{0.25 * sd1:g}", "n1": "10", "df": "z", "reps": "200"}
        approx = ssc.run("sim.lmm.longitudinal", "power", inputs, power_curve=False, approximate=True)
        assert approx["rows"][0].get("_error") == refused, (sd1, approx["rows"][0])
        sim = run("sim.lmm.longitudinal", "power", inputs, threads=2)["rows"][0]
        assert sim.get("_error") == refused and "power" not in sim, (sd1, sim)
    inputs = {"timing": "custom", "times": "10000000000 10000000001 10000000003", "sd1": "1e6", "rho": "0.5",
              "dropout": "0.3", "delta": "250000", "n1": "10", "df": "z"}
    approx = first_row(ssc.run("sim.lmm.longitudinal", "power", inputs, power_curve=False, approximate=True))
    assert approx["power"] == pytest.approx(0.0864997, abs=1e-4)  # 50-digit GLS


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
    (or reported 5% or 100%) for a difference of 0."""
    base = base_inputs(pid)
    for big in ("1e16", "-1e17", "1e30"):
        msg = run(pid, solve, {**base, base_key: big})["rows"][0].get("_error", "")
        assert "cannot be told apart from" in msg and "lost to rounding" in msg, (big, msg)
    msg = run(pid, solve, {**base, base_key: "1e6", **{sd: "1e-30" for sd in sds}})["rows"][0].get("_error", "")
    assert "cannot be told apart from" in msg, msg
    # a large value that still carries the difference gives the answer it gives at 0
    zero = first_row(run(pid, solve, {**base, base_key: "0"}))
    row = first_row(run(pid, solve, {**base, base_key: "1e6"}))
    assert row[key] - 1e6 == pytest.approx(zero[key], rel=1e-6)
    assert row["power"] == pytest.approx(zero["power"], abs=1e-6)


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
    assert "ρ1 = 0.9999996 when the null correlation is ρ0 = 0.999999" in row["_statement"]
    for rho0, alt in (("0.9999999999999", "two"), ("0.9999999999999999", "upper"), ("-0.99999999999", "lower")):
        msg = run("corr.one", "rho1", {**base, "rho0": rho0, "alternative": alt})["rows"][0].get("_error", "")
        assert msg.startswith("ρ1 cannot be told apart from ρ0"), (rho0, msg)
    row = first_row(run("agree.icc", "rho1", {**base_inputs("agree.icc"), "rho0": "0.999999"}))
    assert row["power"] == pytest.approx(0.8, abs=1e-9) and "of 0.9999996 when" in row["_statement"]
    msg = run("agree.icc", "rho1", {**base_inputs("agree.icc"), "rho0": "0.9999999999"})["rows"][0]["_error"]
    assert msg.startswith("ρ1 cannot be told apart from ρ0 = 0.9999999999"), msg


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
    for pid, key, weights in (("means.anova", "ratios", "1 1 1e-20"), ("prop.trend", "weights", "1 1 1 1e-20")):
        base = {**base_inputs(pid), "alloc": "ratio", key: weights}
        msg = run(pid, "n", base)["rows"][0]["_error"]
        assert "is too small: the group would have too few subjects" in msg, msg
        msg = run(pid, "power", base)["rows"][0]["_error"]
        assert "would have no subjects" in msg and "each group needs at least 1" in msg, msg


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
