"""Relations that the case-control, prediction-model, paired diagnostic and co-primary methods, and the ratio and
odds-ratio scales of the two-proportion non-inferiority and equivalence tests, must satisfy; and checks that these
methods answer quickly and do not hold on to memory.

The reference-value tests (core/tests) compare selected results with other software.  These tests instead tie each
new method to an older one, to a closed form or to itself under a change of labels:

* an unmatched case-control study is a comparison of two proportions with the cases' exposure from the odds ratio
  and the controls as a ratio allocation; a matched one with one uncorrelated control per case is Schlesselman's
  matched-pairs formula, and has power alpha at an odds ratio of 1;
* comparing two diagnostic tests on the same subjects is McNemar's test among those with the condition (the
  sensitivities) and among those without (the specificities), and both goals have the product of the two powers;
* co-primary endpoints: uncorrelated endpoints multiply their powers, a correlation between equal endpoints helps,
  and an endpoint that is certain to succeed leaves the others' power unchanged;
* non-inferiority and equivalence on the ratio and odds-ratio scales: alpha at the boundary, the same answer with
  the groups swapped and the ratio inverted, and limits symmetric on the log scale behave symmetrically;
* the prediction-model size never falls with more parameters or rises with a better model, and converted inputs
  equal the Cox–Snell input of their (rounded) converted value.
"""
from __future__ import annotations

import json
import math
import os
import random
import subprocess
import sys
import time
from statistics import NormalDist
from typing import Dict, Optional

import pytest

import ssc
import ssc.core as core

Z = NormalDist()


# ----------------------------------------------------------------------------- helpers
def defaults(pid: str) -> Dict[str, str]:
    return {p["id"]: str(p["default"]) for p in ssc.schema(pid)["params"]}


def text(v) -> str:
    return repr(v) if isinstance(v, float) else str(v)


def solve(pid: str, mode: str, **inputs) -> dict:
    """The first result row for the defaults updated with `inputs`; fails with the engine's message if none."""
    res = ssc.run(pid, mode, {**defaults(pid), **{k: text(v) for k, v in inputs.items()}}, threads=1,
                  power_curve=False)
    assert res.get("ok"), (pid, mode, inputs, res.get("error"), res.get("input_errors"))
    row = res["rows"][0]
    assert "_error" not in row, (pid, mode, inputs, row["_error"])
    return row


def try_solve(pid: str, mode: str, **inputs) -> Optional[dict]:
    """Like solve, but None when the engine has no answer (a message instead of a row)."""
    res = ssc.run(pid, mode, {**defaults(pid), **{k: text(v) for k, v in inputs.items()}}, threads=1,
                  power_curve=False)
    if not res.get("ok") or "_error" in res["rows"][0]:
        return None
    return res["rows"][0]


def power(pid: str, **inputs) -> float:
    return solve(pid, "power", **inputs)["power"]


def single_endpoint_power(test: str, n1: int, n2: int, delta: float, alpha: float = 0.025) -> float:
    """One endpoint of a co-primary design on its own: a one-sided two-sample test."""
    return power("means.two.t", alternative="upper", test=test, alpha=alpha, alloc="sep", n1=n1, n2=n2, mu1=delta,
                 mu2=0, sigma=1)


def odds_to_p(odds_ratio: float, p0: float) -> float:
    return odds_ratio * p0 / (1 - p0 + odds_ratio * p0)


# ============================================================================= unmatched case-control
TWO_PROP_TESTS = ["zpool", "zpoolcc", "zunpool", "fisher"]


@pytest.mark.parametrize("test", TWO_PROP_TESTS)
@pytest.mark.parametrize("method", ["exact", "normal"])
def test_case_control_is_two_proportions_with_controls_as_a_ratio_allocation(test, method):
    """Cases are group 1 with the exposure P1 whose odds are OR times those of the controls; M controls per case
    are a ratio allocation N2 = ⌈M·N1⌉ (one control per case: equal groups).  Same test, same power calculation:
    the same power, size and smallest number of cases."""
    for alt, odds_ratio in (("two", 2.0), ("upper", 2.0), ("lower", 0.5)):
        for m in (1, 2, 3, 1.5, 0.7):
            alloc = {"alloc": "equal"} if m == 1 else {"alloc": "ratio", "r": m}
            common = dict(test=test, method=method, alternative=alt, alpha=0.05)
            cc = solve("prop.cc", "power", **common, ncase=60, m=m, es="or", p0=0.3, **{"or": odds_ratio})
            two = solve("prop.two", "power", **common, n1=60, **alloc, es="or", or1=odds_ratio, p2=0.3)
            assert cc["nctrl"] == two["n2"] and cc["n"] == two["n"], (alt, m)
            assert cc["power"] == pytest.approx(two["power"], rel=1e-9), (alt, m)
            assert cc["p1"] == pytest.approx(two["p1"], rel=1e-12)
            if method == "exact":
                assert cc["actual_alpha"] == pytest.approx(two["actual_alpha"], rel=1e-9), (alt, m)
            if m in (1, 2):  # the searches find the same smallest design
                cc = solve("prop.cc", "n", **common, power=0.8, m=m, es="or", p0=0.3, **{"or": odds_ratio})
                two = solve("prop.two", "n", **common, power=0.8, **alloc, es="or", or1=odds_ratio, p2=0.3)
                assert (cc["ncase"], cc["nctrl"]) == (two["n1"], two["n2"]), (alt, m)
                assert cc["power"] == pytest.approx(two["power"], rel=1e-9)


def test_case_control_exposure_input_is_the_odds_ratio_input():
    for p1 in (0.12, 0.46, 0.8):
        by_p1 = solve("prop.cc", "power", es="p1", p1=p1, p0=0.3, ncase=80, m=2)
        by_or = solve("prop.cc", "power", es="or", p0=0.3, ncase=80, m=2, **{"or": p1 / (1 - p1) / (0.3 / 0.7)})
        assert by_p1["or"] == pytest.approx(by_or["or"], rel=1e-12)
        assert by_p1["power"] == pytest.approx(by_or["power"], rel=1e-9)


# ============================================================================= matched case-control
def schlesselman_pairs(odds_ratio: float, p0: float, alpha: float, target: float, sides: int) -> float:
    """Matched pairs needed (Schlesselman 1982): discordant pairs from the binomial test of P = OR/(1 + OR) against
    1/2, divided by the probability that a pair is discordant."""
    p1 = odds_to_p(odds_ratio, p0)
    big_p = odds_ratio / (1 + odds_ratio)
    discordant = p0 * (1 - p1) + p1 * (1 - p0)
    za, zb = Z.inv_cdf(1 - alpha / sides), Z.inv_cdf(target)
    return (za / 2 + zb * math.sqrt(big_p * (1 - big_p))) ** 2 / (big_p - 0.5) ** 2 / discordant


def schlesselman_power(n: int, odds_ratio: float, p0: float, alpha: float, alt: str) -> float:
    p1 = odds_to_p(odds_ratio, p0)
    big_p = odds_ratio / (1 + odds_ratio)
    shift = math.sqrt(n * (p0 * (1 - p1) + p1 * (1 - p0))) * (big_p - 0.5)
    sd = math.sqrt(big_p * (1 - big_p))
    if alt == "two":
        z = Z.inv_cdf(1 - alpha / 2)
        return Z.cdf((shift - z / 2) / sd) + Z.cdf((-shift - z / 2) / sd)
    z = Z.inv_cdf(1 - alpha)
    return Z.cdf(((shift if alt == "upper" else -shift) - z / 2) / sd)


def test_matched_pairs_without_correlation_follow_schlesselman():
    """One control per case and no correlation of exposure within pairs (M = 1, φ = 0): Dupont's power is the
    normal approximation behind Schlesselman's formula, and the number of pairs is the first that reaches it."""
    for odds_ratio in (0.4, 0.5, 1.5, 2, 3, 5):
        for p0 in (0.05, 0.3, 0.7):
            for alpha, target in ((0.01, 0.9), (0.05, 0.8)):
                for alt in ("two", "upper", "lower"):
                    if (alt == "upper") == (odds_ratio < 1) and alt != "two":
                        continue
                    common = dict(mctl=1, phi=0, p0=p0, alpha=alpha, alternative=alt, **{"or": odds_ratio})
                    n = solve("prop.cc.matched", "n", power=target, **common)["n"]
                    where = (odds_ratio, p0, alpha, alt, n)
                    assert power("prop.cc.matched", n=n, **common) == pytest.approx(
                        schlesselman_power(n, odds_ratio, p0, alpha, alt), abs=1e-12), where
                    assert schlesselman_power(n, odds_ratio, p0, alpha, alt) >= target, where
                    assert schlesselman_power(n - 1, odds_ratio, p0, alpha, alt) < target, where
                    if alt != "two":  # one tail only: exactly the formula, rounded up
                        assert n == math.ceil(schlesselman_pairs(odds_ratio, p0, alpha, target, 1)), where


@pytest.mark.parametrize("alt", ["two", "upper", "lower"])
def test_matched_power_at_an_odds_ratio_of_one_is_alpha(alt):
    for mctl in (1, 2, 4):
        for phi in (-0.1, 0, 0.2):
            for alpha in (0.01, 0.05, 0.1):
                got = power("prop.cc.matched", n=80, mctl=mctl, phi=phi, p0=0.3, alpha=alpha, alternative=alt,
                            **{"or": 1})
                assert got == pytest.approx(alpha, abs=1e-12), (mctl, phi, alpha)


def test_matched_power_grows_with_controls_and_falls_with_correlation():
    """More controls per case help (with diminishing returns); exposures that agree within sets (φ > 0) leave
    fewer informative sets."""
    for phi in (0, 0.2):
        by_m = [power("prop.cc.matched", n=80, mctl=m, phi=phi, p0=0.3, **{"or": 2}) for m in range(1, 8)]
        assert all(a < b for a, b in zip(by_m, by_m[1:])), by_m
        gains = [b - a for a, b in zip(by_m, by_m[1:])]
        assert all(a > b for a, b in zip(gains, gains[1:])), gains
    for mctl in (1, 2, 4):
        by_phi = [power("prop.cc.matched", n=80, mctl=mctl, phi=phi, p0=0.3, **{"or": 2})
                  for phi in (-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4)]
        assert all(a > b for a, b in zip(by_phi, by_phi[1:])), by_phi


# ============================================================================= paired diagnostic tests
def paired_cells(a: float, b: float, dep: str, rho: float, disc: float):
    """Test 1 right only and test 2 right only, from the accuracies a and b and their dependence."""
    if dep == "corr":
        both = a * b + rho * math.sqrt(a * (1 - a) * b * (1 - b))
    elif dep == "disc":
        both = (a + b - disc) / 2
    else:
        both = a * b
    return a - both, b - both


def mcnemar(mode: str, alt: str, method: str, alpha: float, cells, **more) -> dict:
    p10, p01 = cells
    return solve("prop.paired.mcnemar", mode, alternative=alt, method=method, alpha=alpha, es="props", p10=p10,
                 p01=p01, dropout=0, **more)


def accuracies(rng: random.Random, alt: str):
    hi = rng.uniform(0.75, 0.95)
    lo = hi - rng.uniform(0.05, 0.2)
    return (hi, lo) if alt != "lower" else (lo, hi)


def test_paired_diagnostic_power_is_mcnemar_in_each_group():
    """The sensitivities are compared by McNemar's test on ⌊N·prevalence⌋ pairs, the specificities on
    ⌊N·(1 − prevalence)⌋; with both goals the power is the product (different subjects)."""
    rng = random.Random(2026)
    for _ in range(60):
        goal, alt = rng.choice(["se", "sp", "both"]), rng.choice(["two", "upper", "lower"])
        method, dep = rng.choice(["exact", "normal"]), rng.choice(["indep", "corr", "disc"])
        prev, n, alpha = rng.choice([0.1, 0.2, 0.37, 0.5]), rng.randint(50, 900), rng.choice([0.01, 0.05, 0.1])
        (se1, se2), (sp1, sp2), rho = accuracies(rng, alt), accuracies(rng, alt), rng.choice([0, 0.2, 0.4])
        dse, dsp = abs(se1 - se2) + 0.1, abs(sp1 - sp2) + 0.1
        row = solve("diag.paired", "power", goal=goal, alternative=alt, method=method, dep=dep, prev=prev, n=n,
                    alpha=alpha, se1=se1, se2=se2, sp1=sp1, sp2=sp2, rse=rho, rsp=rho, dse=dse, dsp=dsp)
        where = (goal, alt, method, dep, prev, n)
        parts = []
        if goal != "sp":
            pairs = math.floor(n * prev + 1e-9)
            want = mcnemar("power", alt, method, alpha, paired_cells(se1, se2, dep, rho, dse), n=pairs)["power"]
            assert row["nd"] == pairs and row["pse"] == pytest.approx(want, abs=1e-12), where
            parts.append(want)
        if goal != "se":
            pairs = math.floor(n * (1 - prev) + 1e-9)
            want = mcnemar("power", alt, method, alpha, paired_cells(sp1, sp2, dep, rho, dsp), n=pairs)["power"]
            assert row["nn"] == pairs and row["psp"] == pytest.approx(want, abs=1e-12), where
            parts.append(want)
        assert row["power"] == pytest.approx(math.prod(parts), abs=1e-12), where


def test_paired_diagnostic_n_is_the_first_total_with_mcnemars_pairs():
    """⌊N·prevalence⌋ takes every whole value as N grows, so the smallest N is the first whose share of subjects
    with (or without) the condition reaches the smallest number of pairs McNemar's test needs."""
    rng = random.Random(7)
    for _ in range(40):
        goal, alt = rng.choice(["se", "sp"]), rng.choice(["two", "upper", "lower"])
        method, dep = rng.choice(["exact", "normal"]), rng.choice(["indep", "corr", "disc"])
        prev, alpha, target = rng.choice([0.1, 0.2, 0.37, 0.5, 0.75]), rng.choice([0.01, 0.05]), rng.choice([0.8, 0.9])
        a1, a2 = accuracies(rng, alt)
        rho, disc = rng.choice([0, 0.2, 0.4]), abs(a1 - a2) + rng.choice([0.02, 0.1])
        row = solve("diag.paired", "n", goal=goal, alternative=alt, method=method, dep=dep, prev=prev, alpha=alpha,
                    power=target, **{goal + "1": a1, goal + "2": a2, "r" + goal: rho, "d" + goal: disc})
        pairs = mcnemar("n", alt, method, alpha, paired_cells(a1, a2, dep, rho, disc), power=target)["n"]
        share = prev if goal == "se" else 1 - prev
        first = math.ceil(pairs / share - 1e-9)
        while math.floor(first * share + 1e-9) < pairs:
            first += 1
        while math.floor((first - 1) * share + 1e-9) >= pairs:
            first -= 1
        assert row["n"] == first, (goal, alt, method, dep, prev, alpha, target, a1, a2, pairs)


# ============================================================================= co-primary endpoints
def coprimary(**inputs) -> dict:
    return solve("means.coprimary", "power", alpha=0.025, **inputs)


@pytest.mark.parametrize("test", ["z", "t"])
def test_uncorrelated_coprimary_endpoints_multiply_their_powers(test):
    for deltas in ((0.3, 0.5), (0.3, 0.5, 0.4), (0.3, 0.5, 0.4, 0.25)):
        for n1, n2 in ((12, 12), (40, 70)):
            row = coprimary(test=test, alloc="sep", n1=n1, n2=n2, deltas=" ".join(map(str, deltas)), rho=0)
            each = [single_endpoint_power(test, n1, n2, d) for d in deltas]
            assert row["power"] == pytest.approx(math.prod(each), abs=1e-12), (deltas, n1, n2)
            assert row["pmin"] == pytest.approx(min(each), abs=1e-12)


@pytest.mark.parametrize("test", ["z", "t"])
def test_coprimary_power_grows_with_the_correlation(test):
    """With equal effects the joint power rises from (nearly) nothing at the most negative common correlation to
    the single endpoint's power as ρ → 1."""
    for k in (2, 3, 4):
        for delta in (0.2, 0.5):
            rhos = [x / 20 for x in range(-19, 20) if x / 20 > -1 / (k - 1) + 1e-9]
            seen = [coprimary(test=test, n1=60, deltas=" ".join([str(delta)] * k), rho=rho)["power"] for rho in rhos]
            assert all(a <= b + 1e-12 for a, b in zip(seen, seen[1:])), (k, delta, seen)
            assert seen[-1] < single_endpoint_power(test, 60, 60, delta) + 1e-12


@pytest.mark.parametrize("test", ["z", "t"])
def test_coprimary_power_lies_between_the_product_and_the_weakest(test):
    """Slepian's inequality: positively correlated endpoints succeed together at least as often as independent
    ones, negatively correlated ones at most as often; never more often than the weakest alone."""
    rng = random.Random(11)
    for _ in range(30):
        k = rng.choice([2, 3, 4])
        deltas = [round(rng.uniform(0.15, 0.6), 3) for _ in range(k)]
        rho = round(rng.uniform(-1 / (k - 1) + 0.02, 0.95), 3)
        n1, n2 = rng.randint(10, 150), rng.randint(10, 150)
        row = coprimary(test=test, alloc="sep", n1=n1, n2=n2, deltas=" ".join(map(str, deltas)), rho=rho)
        each = [single_endpoint_power(test, n1, n2, d) for d in deltas]
        where = (deltas, rho, n1, n2)
        tol = 1e-9 if k == 2 or rho >= 0 else 3e-5  # three or four endpoints, ρ < 0: a numerical MVN probability
        assert row["power"] <= min(each) + tol, where
        if rho >= 0:
            assert row["power"] >= math.prod(each) - tol, where
        else:
            assert row["power"] <= math.prod(each) + tol, where


def test_coprimary_effects_as_differences_and_sds_and_in_any_order():
    std = coprimary(test="t", n1=50, es="std", deltas="0.3 0.4 0.5", rho=0.4)["power"]
    raw = coprimary(test="t", n1=50, es="raw", diffs="3 2 10", sds="10 5 20", rho=0.4)["power"]
    shuffled = coprimary(test="t", n1=50, es="std", deltas="0.5 0.3 0.4", rho=0.4)["power"]
    assert raw == pytest.approx(std, abs=1e-12) and shuffled == pytest.approx(std, abs=1e-9)


@pytest.mark.parametrize("test,k", [
    ("z", 2), ("t", 2), ("z", 3), ("t", 3), ("z", 4), ("t", 4)])
def test_a_certain_coprimary_endpoint_leaves_the_others_power(test, k):
    """An endpoint with an enormous effect always succeeds, so the joint power is that of the other endpoint alone
    (z: exactly the one-sided z-test; t: each endpoint's exact t power, joined by a Gaussian copula)."""
    for rho in (-0.4, 0, 0.5, 0.9):
        if rho <= -1 / (k - 1):
            continue
        for n1, n2 in ((10, 10), (30, 45), (200, 200)):
            row = try_solve("means.coprimary", "power", test=test, alpha=0.025, alloc="sep", n1=n1, n2=n2,
                            deltas=" ".join(["0.4"] + ["50"] * (k - 1)), rho=rho)
            assert row is not None, (rho, n1, n2)
            assert row["power"] == pytest.approx(single_endpoint_power(test, n1, n2, 0.4), abs=1e-9), (rho, n1, n2)


@pytest.mark.parametrize("test", ["z", "t"])
def test_coprimary_n_with_a_certain_endpoint_is_the_single_endpoint_n(test):
    """δ = 0.8 next to δ = 0.2: at the N the weak endpoint needs (about 394 per group) the strong one's power is
    1 to double precision, so N is the weak endpoint's own."""
    single = solve("means.two.t", "n", alternative="upper", test=test, alpha=0.025, power=0.8, mu1=0.2, mu2=0,
                   sigma=1)["n1"]
    for rho in (-0.3, 0, 0.5, 0.9):
        row = try_solve("means.coprimary", "n", test=test, alpha=0.025, power=0.8, deltas="0.2 0.8", rho=rho)
        assert row is not None and row["n1"] == single, (rho, row and row["n1"], single)


# ============================================================================= NI and equivalence: ratio, odds ratio
SCALE_KEYS = {"ratio": ("r0", "ratio", "r1"), "or": ("or0", "or", "or1")}
LIMIT_KEYS = {"diff": ("dl", "du", "diff", "d1"), "ratio": ("rl", "ru", "ratio", "r1"), "or": ("orl", "oru", "or", "or1")}


def boundary_p1(scale: str, bound: float, p2: float) -> float:
    return {"diff": p2 + bound, "ratio": bound * p2, "or": odds_to_p(bound, p2)}[scale]


@pytest.mark.parametrize("scale", ["ratio", "or"])
def test_non_inferiority_power_at_the_boundary_is_alpha(scale):
    """At the boundary the restricted estimates are the true proportions, so the score statistic is standard
    normal: the normal power is alpha whatever the proportions, groups and side."""
    key, es, actual = SCALE_KEYS[scale]
    for bound in (0.7, 0.8, 1.25) if scale == "ratio" else (0.5, 0.67, 1.5, 2):
        for h1 in ("upper", "lower"):
            for p2 in (0.1, 0.4, 0.7):
                for alpha in (0.025, 0.05):
                    for n1, n2 in ((50, 50), (80, 160)):
                        got = power("prop.two.ni", scale=scale, h1=h1, method="normal", alpha=alpha, alloc="sep",
                                    n1=n1, n2=n2, es=es, p2=p2, **{key: bound, actual: bound})
                        assert got == pytest.approx(alpha, abs=1e-12), (bound, h1, p2, alpha, n1, n2)


def size_cases(seed: int, count: int):
    """Random designs (scale, alpha, group sizes, P2), with the generator for the caller's further choices."""
    rng = random.Random(seed)
    for _ in range(count):
        scale, alpha = rng.choice(["diff", "ratio", "or"]), rng.choice([0.025, 0.05])
        n1 = rng.randint(40, 300)
        n2 = rng.choice([n1, rng.randint(40, 300)])
        yield scale, alpha, n1, n2, round(rng.uniform(0.05, 0.95), 3), rng


def test_exact_non_inferiority_size_stays_near_alpha():
    """The score tests are not exact, but with at least 10 expected of each outcome in each group their exact
    rejection rate at the boundary stays within 0.01 of alpha on every scale (the difference scale reaches about
    0.007, the ratio scale 0.008); the reported actual alpha is the exact power at the boundary."""
    checked = 0
    for scale, alpha, n1, n2, p2, rng in size_cases(2026, 150):
        bound = rng.choice({"diff": [-0.15, -0.1, -0.05, 0.05, 0.1], "ratio": [0.7, 0.8, 0.9, 1.1, 1.25],
                            "or": [0.5, 0.67, 1.5, 2]}[scale])
        p1 = boundary_p1(scale, bound, p2)
        if not 0 < p1 < 1 or min(n1 * p1, n1 * (1 - p1), n2 * p2, n2 * (1 - p2)) < 10:
            continue
        key = {"diff": "d0", "ratio": "r0", "or": "or0"}[scale]
        h1 = rng.choice(["upper", "lower"])
        row = solve("prop.two.ni", "power", scale=scale, h1=h1, method="exact", alpha=alpha, alloc="sep", n1=n1,
                    n2=n2, es="props", p1=p1, p2=p2, **{key: bound})
        where = (scale, bound, h1, alpha, n1, n2, p2)
        assert row["power"] == pytest.approx(row["actual_alpha"], abs=1e-12), where
        assert row["actual_alpha"] <= alpha + 0.01, where
        checked += 1
    assert checked >= 60


def test_equivalence_power_at_a_limit_is_at_most_alpha():
    """Just inside a limit the TOST rejects at most as often as the one test for that limit: alpha with the
    normal approximation, and about alpha when the other limit is far away; exactly, within 0.01 of alpha."""
    checked = 0
    for scale, alpha, n1, n2, p2, rng in size_cases(7, 150):
        lk, uk, es, actual = LIMIT_KEYS[scale]
        lo, hi = rng.choice({"diff": [(-0.1, 0.1), (-0.15, 0.05)], "ratio": [(0.8, 1.25), (0.67, 1.5), (0.8, 1.9)],
                             "or": [(0.67, 1.5), (0.5, 2), (0.67, 50)]}[scale])
        at = rng.choice([lo, hi])
        at += 1e-9 * (1 if at == lo else -1) * (1 if scale == "diff" else at)
        p1 = boundary_p1(scale, at, p2)
        if scale == "ratio" and hi * p2 >= 1 or not 0 < p1 < 1:
            continue
        if min(n1 * p1, n1 * (1 - p1), n2 * p2, n2 * (1 - p2)) < 10:
            continue
        common = dict(scale=scale, alpha=alpha, alloc="sep", n1=n1, n2=n2, es=es, p2=p2, **{lk: lo, uk: hi, actual: at})
        normal = power("prop.two.equiv", method="normal", **common)
        exact = power("prop.two.equiv", method="exact", **common)
        where = (scale, lo, hi, at, alpha, n1, n2, p2)
        assert normal <= alpha + 1e-6, where
        assert exact <= alpha + 0.01, where
        checked += 1
    assert checked >= 50
    # with the other limit far away the TOST is the one test: alpha
    for scale, lo, hi in (("ratio", 0.8, 1.9), ("or", 0.67, 50)):
        lk, uk, es, actual = LIMIT_KEYS[scale]
        for at in (lo * (1 + 1e-9), hi * (1 - 1e-9)):
            got = power("prop.two.equiv", scale=scale, method="normal", alpha=0.05, n1=150, es=es, p2=0.4,
                        **{lk: lo, uk: hi, actual: at})
            assert got == pytest.approx(0.05, abs=2e-4), (scale, at)


SWAP_SIDE = {"upper": "lower", "lower": "upper"}


@pytest.mark.parametrize("method", ["normal", "exact"])
def test_swapping_the_groups_inverts_the_ratio_and_the_odds_ratio(method):
    """P1/P2 > R0 is P2/P1 < 1/R0, and the score statistics change only their sign (the restricted estimates are
    the same pair), so the swapped design has the same power, normal or exact."""
    rng = random.Random(3)
    checked = 0
    for _ in range(80):
        scale, alpha = rng.choice(["ratio", "or"]), rng.choice([0.025, 0.05])
        n1, n2 = rng.randint(20, 200), rng.randint(20, 200)
        p1, p2 = round(rng.uniform(0.1, 0.9), 3), round(rng.uniform(0.1, 0.9), 3)
        h1, bound = rng.choice(["upper", "lower"]), rng.choice([0.5, 0.67, 0.8, 0.9, 1.1, 1.25, 1.5, 2])
        key = SCALE_KEYS[scale][0]
        common = dict(scale=scale, method=method, alpha=alpha, alloc="sep", es="props")
        a = try_solve("prop.two.ni", "power", h1=h1, n1=n1, n2=n2, p1=p1, p2=p2, **{key: bound}, **common)
        b = try_solve("prop.two.ni", "power", h1=SWAP_SIDE[h1], n1=n2, n2=n1, p1=p2, p2=p1, **{key: 1 / bound}, **common)
        if a is None or b is None:  # the boundary R0·P2 (or P1/R0) falls outside 0 to 1 in one of them
            continue
        assert a["power"] == pytest.approx(b["power"], abs=1e-12), (scale, h1, bound, n1, n2, p1, p2)
        lk, uk = LIMIT_KEYS[scale][:2]
        lo, hi = rng.choice([(0.8, 1.25), (0.67, 1.5), (0.5, 2), (0.7, 1.2), (0.9, 1.5)])
        q2 = round(rng.uniform(0.1, 0.6 if scale == "ratio" else 0.9), 3)
        q1 = round(boundary_p1(scale, rng.uniform(lo, hi), q2), 4)
        a = try_solve("prop.two.equiv", "power", n1=n1, n2=n2, p1=q1, p2=q2, **{lk: lo, uk: hi}, **common)
        b = try_solve("prop.two.equiv", "power", n1=n2, n2=n1, p1=q2, p2=q1, **{lk: 1 / hi, uk: 1 / lo}, **common)
        if a is None or b is None:  # RU·P2 at or above 1 in one of them
            continue
        assert a["power"] == pytest.approx(b["power"], abs=1e-12), (scale, lo, hi, n1, n2, q1, q2)
        checked += 1
    assert checked >= 20


@pytest.mark.parametrize("method", ["normal", "exact"])
def test_limits_symmetric_on_the_log_scale_behave_symmetrically(method):
    """Limits 1/U and U map onto themselves: on the ratio scale swapping the groups, on the odds-ratio scale also
    counting the other outcome (P → 1 − P inverts the odds ratio), leaves the power unchanged.  A non-inferiority
    margin OR0 on one outcome is 1/OR0 on the other with the other side."""
    rng = random.Random(5)
    for _ in range(25):
        u = rng.choice([1.25, 1.5, 2.0])
        n1, n2, alpha = rng.randint(30, 200), rng.randint(30, 200), rng.choice([0.025, 0.05])
        common = dict(method=method, alpha=alpha, alloc="sep", es="props")
        q2 = round(rng.uniform(0.15, 0.5), 3)
        q1 = round(q2 * rng.uniform(1 / u, u), 4)
        if q1 * u < 1 and q2 * u < 1:
            a = power("prop.two.equiv", scale="ratio", n1=n1, n2=n2, p1=q1, p2=q2, rl=1 / u, ru=u, **common)
            b = power("prop.two.equiv", scale="ratio", n1=n2, n2=n1, p1=q2, p2=q1, rl=1 / u, ru=u, **common)
            assert a == pytest.approx(b, abs=1e-12), ("ratio", u, n1, n2, q1, q2)
        p2 = round(rng.uniform(0.1, 0.9), 3)
        p1 = round(odds_to_p(rng.uniform(1 / u, u), p2), 4)
        lims = dict(scale="or", orl=1 / u, oru=u, **common)
        base = power("prop.two.equiv", n1=n1, n2=n2, p1=p1, p2=p2, **lims)
        flipped = power("prop.two.equiv", n1=n1, n2=n2, p1=round(1 - p1, 4), p2=round(1 - p2, 4), **lims)
        swapped = power("prop.two.equiv", n1=n2, n2=n1, p1=p2, p2=p1, **lims)
        assert flipped == pytest.approx(base, abs=1e-12) and swapped == pytest.approx(base, abs=1e-12), (u, p1, p2)
        h1, or0 = rng.choice(["upper", "lower"]), rng.choice([0.5, 0.67, 1.5, 2])
        ni = dict(scale="or", n1=n1, n2=n2, **common)
        a = power("prop.two.ni", h1=h1, or0=or0, p1=p1, p2=p2, **ni)
        b = power("prop.two.ni", h1=SWAP_SIDE[h1], or0=1 / or0, p1=round(1 - p1, 4), p2=round(1 - p2, 4), **ni)
        assert a == pytest.approx(b, abs=1e-12), (h1, or0, p1, p2)


# ============================================================================= prediction models
PRED_CASES = [
    dict(outcome="binary", r2in="cs", r2cs=0.1, prev=0.174),
    dict(outcome="binary", r2in="nag", r2nag=0.2, prev=0.3),
    dict(outcome="binary", r2in="c", cstat=0.75, prev=0.1),
    dict(outcome="survival", r2in="cs", r2cs=0.05, rate=0.065, fup=2.07),
    dict(outcome="survival", r2in="nag", r2nag=0.3, rate=0.2, fup=3),
    dict(outcome="continuous", r2=0.2),
]


def prediction_n(**inputs) -> dict:
    return solve("reg.prediction", "n", **inputs)


@pytest.mark.parametrize("case", PRED_CASES, ids=lambda c: f"{c['outcome']}-{c.get('r2in', 'r2')}")
def test_prediction_n_never_falls_with_more_parameters(case):
    sizes = [prediction_n(p=p, **case)["n"] for p in (1, 2, 3, 5, 8, 13, 20, 30, 50, 80, 120)]
    assert all(a <= b for a, b in zip(sizes, sizes[1:])), sizes
    assert sizes[-1] > sizes[0]


@pytest.mark.parametrize("outcome,key,values,extra", [
    ("binary", "r2cs", (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.45), dict(prev=0.3)),
    ("binary", "r2nag", (0.01, 0.05, 0.1, 0.2, 0.3, 0.45), dict(r2in="nag", prev=0.2)),
    ("binary", "cstat", (0.55, 0.6, 0.7, 0.8, 0.9, 0.95), dict(r2in="c", prev=0.2)),
    ("survival", "r2cs", (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.45), dict(rate=0.065, fup=2.07)),
    ("continuous", "r2", (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.45), {}),
])
def test_prediction_n_never_rises_with_a_better_model(outcome, key, values, extra):
    for p in (5, 24, 60):
        sizes = [prediction_n(outcome=outcome, p=p, **extra, **{key: v})["n"] for v in values]
        assert all(a >= b for a, b in zip(sizes, sizes[1:])), (p, sizes)


def test_prediction_n_never_falls_with_stricter_targets():
    for outcome, key, values, extra in (("binary", "shrink", (0.8, 0.85, 0.9, 0.95), {}),
                                        ("binary", "delta", (0.1, 0.05, 0.02, 0.01), dict(r2cs=0.3, prev=0.5)),
                                        ("binary", "margin", (0.1, 0.05, 0.03, 0.02), dict(p=2)),
                                        ("continuous", "shrink", (0.8, 0.85, 0.9, 0.95), {}),
                                        ("continuous", "delta", (0.1, 0.05, 0.02, 0.01), {}),
                                        ("continuous", "mmoe", (1.3, 1.2, 1.1, 1.05), dict(p=2))):
        sizes = [prediction_n(outcome=outcome, **extra, **{key: v})["n"] for v in values]
        assert all(a <= b for a, b in zip(sizes, sizes[1:])), (outcome, key, sizes)
        assert sizes[-1] > sizes[0], (outcome, key, sizes)


def largest_cox_snell(outcome: str, prev: float = 0, rate: float = 0, fup: float = 0) -> float:
    """The largest Cox–Snell R² used to convert a Nagelkerke R² (per person for a time-to-event outcome)."""
    if outcome == "binary":
        return -math.expm1(2 * (prev * math.log(prev) + (1 - prev) * math.log1p(-prev)))
    events = rate * fup
    return -math.expm1(2 * (events * math.log(events) - events))


def test_a_nagelkerke_input_is_the_cox_snell_input_of_its_rounded_conversion():
    """As pmsampsize, a Nagelkerke R² becomes the Cox–Snell R² R²Nag × max R²CS rounded to 3 decimals; giving that
    Cox–Snell R² directly gives the same sizes."""
    rng = random.Random(4)
    for _ in range(80):
        outcome, p, nag = rng.choice(["binary", "survival"]), rng.randint(1, 60), round(rng.uniform(0.02, 0.6), 3)
        if outcome == "binary":
            extra = dict(prev=round(rng.uniform(0.02, 0.6), 3))
        else:
            extra = dict(rate=round(rng.uniform(0.01, 0.3), 3), fup=round(rng.uniform(0.5, 5), 2), tpt=2)
        shrink = rng.choice([0.85, 0.9, 0.95])
        a = prediction_n(outcome=outcome, p=p, r2in="nag", r2nag=nag, shrink=shrink, **extra)
        converted = round(nag * largest_cox_snell(outcome, **{k: v for k, v in extra.items() if k != "tpt"}), 3)
        assert a["r2cs"] == pytest.approx(converted, abs=1e-12), (outcome, nag, extra)
        b = prediction_n(outcome=outcome, p=p, r2in="cs", r2cs=a["r2cs"], shrink=shrink, **extra)
        assert (a["n"], a["n_shrink"], a["n_opt"]) == (b["n"], b["n_shrink"], b["n_opt"]), (outcome, p, nag, extra)


def cox_snell_from_c(c: float, prev: float, points: int = 4001) -> float:
    """Riley et al. (2022): a linear predictor normal with variance 1 among those with and without the outcome,
    means √2·Φ⁻¹(C) apart; R²CS = 1 − exp(−2I) with I the outcome's mutual information with it (trapezoid rule)."""
    mu = math.sqrt(2) * Z.inv_cdf(c)
    lp, lq = math.log(prev), math.log1p(-prev)

    def lse(a, b):
        return max(a, b) + math.log1p(math.exp(-abs(a - b)))

    lo, hi = -12.0, 12.0 + mu
    h = (hi - lo) / (points - 1)
    total = 0.0
    for i in range(points):
        x = lo + i * h
        w = h / 2 if i in (0, points - 1) else h
        total += w * (prev * Z.pdf(x - mu) * -lse(lp, lq - mu * x + mu * mu / 2)
                      + (1 - prev) * Z.pdf(x) * -lse(lp + mu * x - mu * mu / 2, lq))
    return -math.expm1(-2 * total)


def test_a_c_statistic_input_is_the_cox_snell_input_of_its_rounded_conversion():
    for c in (0.6, 0.7, 0.8, 0.9):
        for prev in (0.05, 0.174, 0.5):
            a = prediction_n(outcome="binary", p=10, r2in="c", cstat=c, prev=prev)
            assert a["r2cs"] == pytest.approx(cox_snell_from_c(c, prev), abs=5e-5 + 1e-7), (c, prev)
            b = prediction_n(outcome="binary", p=10, r2in="cs", r2cs=a["r2cs"], prev=prev)
            assert (a["n"], a["n_shrink"], a["n_opt"], a["n_risk"]) == (b["n"], b["n_shrink"], b["n_opt"], b["n_risk"])


def test_the_residual_sd_criterion_is_234_plus_p():
    """Riley et al. (2020): a residual SD within a factor of 1.1 needs 234 + p subjects."""
    for p in (1, 2, 5, 10, 24, 50, 100, 500):
        row = prediction_n(outcome="continuous", p=p, mmoe=1.1)
        assert row["n_sd"] == 234 + p, (p, row["n_sd"])
        assert row["n"] >= row["n_sd"]


# ============================================================================= time and memory
def resource_cases():
    """Every solve mode of the new methods at their defaults, their examples, and the ratio and odds-ratio scales
    of non-inferiority and equivalence with both power calculations."""
    cases = []
    for pid in ("prop.cc", "prop.cc.matched", "reg.prediction", "diag.paired", "means.coprimary"):
        schema = ssc.schema(pid)
        for mode in [o["value"] for o in schema["solve_for"]]:
            cases.append((pid, mode, defaults(pid)))
        cases.append((pid, schema["example"]["solve_for"], {**defaults(pid), **schema["example"]["inputs"]}))
    for pid in ("prop.two.ni", "prop.two.equiv"):
        for scale in ("ratio", "or"):
            for method in ("normal", "exact"):
                for mode in [o["value"] for o in ssc.schema(pid)["solve_for"]]:
                    cases.append((pid, mode, {**defaults(pid), "scale": scale, "method": method}))
    return cases


RESOURCE_CASES = resource_cases()


@pytest.mark.parametrize("pid,mode,inputs", RESOURCE_CASES,
                         ids=[f"{c[0]}:{c[1]}:{i}" for i, c in enumerate(RESOURCE_CASES)])
def test_new_methods_answer_quickly(pid, mode, inputs):
    """At most 0.25 s each on a 4-core machine (the example of prop.cc, eight exact searches, is the slowest),
    curves included: 5 s is for a slowdown by an order of magnitude, not for noise."""
    best = math.inf
    for _ in range(2):
        start = time.perf_counter()
        res = ssc.run(pid, mode, inputs, power_curve=True)
        best = min(best, time.perf_counter() - start)
        assert res["ok"] and all("_error" not in r for r in res["rows"]), (pid, mode, res["rows"][0])
        if best < 5:
            break
    assert best < 5, f"{pid} [{mode}] took {best:.1f} s"


MEMORY = """
import gc, json, sys, time
import ssc.core as core

def mem():
    d = {}
    for line in open("/proc/self/status"):
        if line.startswith(("VmRSS", "VmHWM")):
            d[line.split(":")[0]] = int(line.split()[1]) / 1024
    return d

cases = json.loads(sys.stdin.read())
def run(case, curve):
    res = core.run(case[0], case[1], case[2], power_curve=curve)
    assert res["ok"] and all("_error" not in r for r in res["rows"]), case

counts = []
for case in cases:                      # warm up (and size each case's repetitions by its time)
    run(case, True)
    start = time.perf_counter()
    run(case, False)
    counts.append(max(3, min(300, int(0.08 / max(time.perf_counter() - start, 1e-6)))))
gc.collect()
before = mem()
for case, count in zip(cases, counts):
    run(case, True)
    for _ in range(count):
        run(case, False)
gc.collect()
print(json.dumps({"calls": sum(counts) + len(cases), "before": before, "after": mem()}))
"""


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="reads /proc/self/status")
def test_new_methods_do_not_hold_on_to_memory():
    """Thousands of runs of the new methods (cheap ones hundreds of times) in a process of their own.  Measured on
    Linux: 0.11 MB more after about 5,500 runs, while a native leak of 1 KB a call injected into the same loop
    adds 5.9 MB (and 256 bytes a call 1.8 MB), so a bound of 3 MB catches a leak of about 600 bytes a call."""
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([os.path.dirname(os.path.dirname(core.__file__)),
                                                       os.environ.get("PYTHONPATH", "")]))
    out = subprocess.run([sys.executable, "-c", MEMORY], input=json.dumps(RESOURCE_CASES), capture_output=True,
                         text=True, env=env, timeout=600)
    assert out.returncode == 0, out.stderr
    d = json.loads(out.stdout.strip().splitlines()[-1])
    assert d["calls"] >= 3 * len(RESOURCE_CASES), d
    assert d["after"]["VmRSS"] - d["before"]["VmRSS"] < 3, d
    assert d["after"]["VmHWM"] - d["before"]["VmHWM"] < 8, d
