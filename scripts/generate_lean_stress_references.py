#!/usr/bin/env python3
"""Exhaustive path references, independent of Lean's state-merging compiler."""
import itertools
import json
from fractions import Fraction
from pathlib import Path
from random import Random

rng = Random(541903)
cases = []
for i in range(32):
    scores = [-2, 0, 1, 3 if i % 4 else 1]
    weights = [rng.randrange(6) for _ in scores]
    other = [rng.randrange(6) for _ in scores]
    if not sum(weights):
        weights[0] = 1
    if not sum(other):
        other[-1] = 1
    n, a, b = 2 + i % 4, rng.choice([-1, 1]), rng.randrange(-2, 3)
    threshold, modulus, residue = rng.randrange(-4, 5), 2 + i % 3, i % 2
    probability, mean, total = Fraction(0), Fraction(0), Fraction(0)
    for path in itertools.product(range(len(scores)), repeat=n):
        state, mass = 0, Fraction(1)
        for choice in path:
            ws = weights if state % 2 == 0 else other
            mass *= Fraction(ws[choice], sum(ws))
            state = a * state + scores[choice] + (b if state % 3 == 0 else -b)
        total += mass
        mean += mass * state
        if state >= threshold and state % modulus != residue:
            probability += mass
    assert total == 1
    fmt = lambda x: f"{x.numerator}/{x.denominator}"
    cases.append(dict(id=i, n=n, a=a, b=b, scores=scores, weights=weights,
                      other=other, threshold=threshold, modulus=modulus,
                      residue=residue, probability=fmt(probability),
                      complement=fmt(1-probability), mean=fmt(mean)))
target = Path(__file__).resolve().parents[1] / "web/tests/fixtures/lean-stress.json"
target.write_text(json.dumps(dict(seed=541903, method="exhaustive paths with fractions.Fraction",
                                cases=cases), indent=2) + "\n")
print(f"Wrote {len(cases)} independent observation-program references to {target}")
