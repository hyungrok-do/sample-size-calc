#!/usr/bin/env python3
"""Independent exact references for the observation-program examples.

Uses combinatorial binomial/multinomial formulas and exhaustive Markov paths,
not the Lean compiler's state-merging interpreter. No floating-point arithmetic.
"""
from fractions import Fraction as Q
from itertools import product
import json
from math import comb, factorial
from pathlib import Path


def binary(n, p):
    return {k: comb(n, k) * p**k * (1-p)**(n-k) for k in range(n+1)}


def paired(n, minus, zero, plus):
    law = {k: Q(0) for k in range(-n, n+1)}
    for a in range(n+1):
        for b in range(n-a+1):
            c = n-a-b
            multiplicity = factorial(n) // (factorial(a)*factorial(b)*factorial(c))
            law[a-b] += multiplicity * plus**a * minus**b * zero**c
    return law


def dependent(n, after_false, after_true):
    law = {k: Q(0) for k in range(n+1)}
    for observations in product([False, True], repeat=n):
        previous = False
        mass = Q(1)
        for observation in observations:
            p = after_true if previous else after_false
            mass *= p if observation else 1-p
            previous = observation
        law[sum(observations)] += mass
    return law


def calibrated(n, null, alternative):
    assert sum(null.values()) == sum(alternative.values()) == 1
    def tail(law, c):
        return sum((p for k, p in law.items() if k >= c), Q(0))
    c = next(c for c in range(n+2) if tail(null, c) <= Q(1, 20))
    return dict(n=n, cutoff=c, typeI=str(tail(null, c)), power=str(tail(alternative, c)))


def main():
    fixture = {
        "provenance": "Exact Fraction arithmetic; combinatorial binomial/multinomial formulas and exhaustive Markov trajectories. Regenerate with scripts/generate_lean_derivation_references.py.",
        "bernoulli-observations": [calibrated(n, binary(n, Q(1, 10)), binary(n, Q(3, 10)))
                                   for n in range(1, 31)],
        "paired-categorical": [calibrated(n, paired(n, Q(3, 20), Q(7, 10), Q(3, 20)),
                                         paired(n, Q(1, 20), Q(3, 5), Q(7, 20)))
                               for n in range(1, 31)],
        "dependent-observations": [calibrated(n, dependent(n, Q(1, 10), Q(7, 10)),
                                             dependent(n, Q(1, 4), Q(4, 5)))
                                   for n in range(1, 13)],
    }
    destination = Path(__file__).resolve().parents[1]/"web/tests/fixtures/lean-derivations.json"
    destination.write_text(json.dumps(fixture, indent=2)+"\n")
    print(f"Wrote 72 independent exact design references to {destination}")


if __name__ == "__main__":
    main()
