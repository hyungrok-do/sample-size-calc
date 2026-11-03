#!/usr/bin/env python3
"""Independent numerical checks and native compilation of generated C++.

Requires SciPy only for the independent reference. Lean checks are separate.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile
from scipy.optimize import brentq
from scipy.stats import norm


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", type=Path, default=root / "build" / "ssc_cli")
    args = parser.parse_args()

    def run(plan):
        response = subprocess.run([str(args.engine.resolve()), "design", "-"],
                                  input=json.dumps(plan), text=True, capture_output=True, check=True)
        return json.loads(response.stdout)["result"]

    def power(delta, variance, alpha, n):
        z = norm.isf(alpha / 2)
        noncentrality = abs(delta) * math.sqrt(n / variance)
        return norm.cdf(-z - noncentrality) + norm.sf(z - noncentrality)

    checked = 0
    with tempfile.TemporaryDirectory(prefix="ssc-cpp-") as temp:
        normal_source = (root / "examples/design/normal.lean").read_text()
        cases = [
            ("normal", normal_source, 63),
            ("influence", (root / "examples/design/influence.lean").read_text(), 126),
            # Seven-digit counts must survive C++ emission without decimal rounding.
            ("large_count", normal_source.replace("alternative := 0.5", "alternative := 0.004")
             .replace("minN := 1", "minN := 1234567")
             .replace("maxN := 10000", "maxN := 1234579"), 1234567),
        ]
        for kind, source, count in cases:
            for mode in ["sample_size", "power", "effect"]:
                request = {"op": "compile", "source": source, "solveFor": mode}
                if mode != "sample_size":
                    request["n"] = count
                compiled = run(request)
                result = run(compiled["plan"])
                cpp, exe = Path(temp) / "study.cpp", Path(temp) / "study"
                cpp.write_text(compiled["cpp"])
                subprocess.run(["c++", "-std=c++17", "-O2", str(cpp), "-o", str(exe)], check=True)
                standalone = json.loads(subprocess.check_output([str(exe)], text=True))
                key = {"sample_size": "selectedN", "power": "power", "effect": "detectableMagnitude"}[mode]
                assert math.isclose(standalone[key], result[key], abs_tol=2e-12, rel_tol=2e-12), (kind, mode, standalone, result)
                checked += 1

    template = json.loads((root / "examples/design/normal.json").read_text())
    comparisons = 0
    for effect in [-1.5, -.5, .1, .5, 2.0]:
        for sd in [.5, 1.0, 2.0]:
            for alpha in [.001, .01, .05, .2]:
                for target in [.5, .8, .95]:
                    plan = {**template, "parameters": {"alternative": effect, "nullValue": 0, "sd": sd,
                            "alpha": alpha, "targetPower": target}, "search": {"minN": 1, "maxN": 100000}}
                    result = run(plan)
                    root_n = brentq(lambda n: power(effect, 2*sd*sd, alpha, n)-target, 1e-10, 100000)
                    expected = max(1, math.ceil(root_n))
                    assert result["selectedN"] == expected, (plan, result, expected)
                    assert math.isclose(result["power"], power(effect, 2*sd*sd, alpha, expected), abs_tol=2e-13)
                    comparisons += 1
    print(f"{checked} generated C++ programs compiled and agreed; {comparisons} independent SciPy sample-size/power comparisons passed.")


if __name__ == "__main__":
    main()
