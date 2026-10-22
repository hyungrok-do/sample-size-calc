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
    """The t, Welch and Mann–Whitney tests are unchanged by a common change of location and scale, so SDs of 1e300
    give the power of SDs of 1 (instead of overflowing to power 0)."""
    def row(scale):
        return first_row(run("sim.means.robust", "power", {**base_inputs("sim.means.robust"), "sd1": f"{scale:g}",
                                                           "sd2": f"{2 * scale:g}", "mu1": f"{scale:g}", "mu2": "0",
                                                           "n1": "40", "reps": "200"}, threads=1))
    one, huge = row(1), row(1e300)
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
