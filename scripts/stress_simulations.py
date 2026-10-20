"""Stress test of every simulation procedure: the long version of tests/test_stress.py.

    PYTHONPATH=python python scripts/stress_simulations.py [part ...] [--json report.json]

Parts (all by default):
  fuzz         random inputs from each schema (via tests/test_properties.py): a result or a clear message, quickly
  threads      1, 2, 3, 4, 7 and 16 threads give identical results
  concurrency  every simulation run at once from Python threads gives the same results as one after another
  cancel       a run of a million studies stops promptly when cancelled
  calibration  the type I error at null designs (20,000 studies each), as a z-score against alpha
  bigreps      a million simulated studies: time and memory
  extremes     huge samples, near-degenerate correlations, almost everyone lost: a result or a clear message
"""
import json
import math
import random
import resource
import sys
import threading
import time

ROOT = __import__("pathlib").Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests"))
import ssc  # noqa: E402
import ssc.core as core  # noqa: E402
import test_properties as tp  # noqa: E402

SIMS = [p["id"] for p in ssc.catalog()["procedures"] if ssc.schema(p["id"]).get("simulation")]
report = {}


def run(pid, solve, inputs, threads=0, token=None):
    return core.run(pid, solve, inputs, threads=threads, power_curve=False, token=token)


def strip(res):
    """The rows without timing fields."""
    rows = res.get("rows") or []
    return json.dumps(rows, sort_keys=True)


def fuzz(cases=25):
    out = {}
    for pid in SIMS:
        problems, worst = [], 0.0
        for solve, inputs in tp.fuzz_cases(pid, cases, seed=99):
            inputs = dict(inputs, reps=str(random.Random(str(inputs)).choice([100, 200, 300])))
            t = time.perf_counter()
            res = run(pid, solve, inputs)
            dt = time.perf_counter() - t
            worst = max(worst, dt)
            problems += tp.check_result(pid, solve, inputs, res, dt)
        out[pid] = {"cases": cases, "worst_s": round(worst, 2), "problems": problems[:5]}
        print("fuzz", pid, out[pid], flush=True)
    report["fuzz"] = out


def threads():
    out = {}
    for pid in SIMS:
        bad = []
        for solve, inputs in tp.fuzz_cases(pid, 4, seed=7):
            inputs = dict(inputs, reps="200")
            ref = strip(run(pid, solve, inputs, threads=1))
            for t in (2, 3, 4, 7, 16):
                if strip(run(pid, solve, inputs, threads=t)) != ref:
                    bad.append((solve, inputs, t))
        out[pid] = bad[:3]
        print("threads", pid, "identical" if not bad else bad[:2], flush=True)
    report["threads"] = out


def concurrency():
    jobs = []
    for pid in SIMS:
        ex = ssc.schema(pid)["example"]
        jobs.append((pid, "power", dict(ex["inputs"], reps="300")))
        jobs.append((pid, ex["solve_for"], dict(ex["inputs"], reps="300")))
    sequential = [strip(run(*j)) for j in jobs]
    results = [None] * len(jobs)

    def work(i):
        results[i] = strip(run(*jobs[i]))

    t = time.perf_counter()
    ths = [threading.Thread(target=work, args=(i,)) for i in range(len(jobs))]
    for th in ths:
        th.start()
    for th in ths:
        th.join()
    same = sum(a == b for a, b in zip(sequential, results))
    report["concurrency"] = {"jobs": len(jobs), "identical": same, "seconds": round(time.perf_counter() - t, 1)}
    print("concurrency", report["concurrency"], flush=True)


def cancel(rounds=5):
    out = {}
    for pid in SIMS:
        ex = ssc.schema(pid)["example"]
        lat = []
        ok = True
        for k in range(rounds):
            tok = core.Token()
            box = {}
            th = threading.Thread(target=lambda: box.setdefault("r", run(pid, "power", dict(ex["inputs"], reps="1000000"),
                                                                       token=tok)))
            th.start()
            time.sleep(0.05 + 0.1 * k)
            t = time.perf_counter()
            tok.cancel()
            th.join(10)
            lat.append(time.perf_counter() - t)
            ok = ok and not th.is_alive() and box.get("r", {}).get("cancelled", False)
            tok.close()
        out[pid] = {"cancelled": ok, "max_latency_ms": round(1000 * max(lat), 1)}
        print("cancel", pid, out[pid], flush=True)
    report["cancel"] = out


NULLS = {
    "sim.lmm.cluster": [{"mu1": "0", "k1": "6", "cv": "0.6"}, {"mu1": "0", "k1": "20", "icc": "0.2"}],
    "sim.glmm.cluster": [{"p1": "0.2", "k1": "15"}, {"p1": "0.2", "k1": "30", "m": "40", "cv": "0.5"}],
    "sim.lmm.longitudinal": [{"delta": "0", "n1": "20", "dropout": "0.3"}, {"delta": "0", "re": "intercepts"}],
    "sim.surv.nph": [{"hr1": "1", "hr2": "1"}, {"hr1": "1", "hr2": "1", "test": "maxcombo"},
                     {"hr1": "1", "hr2": "1", "test": "rmst"}],
    "sim.means.robust": [{"mu1": "0", "n1": "15", "dist": "lognormal", "skew": "3"},
                         {"mu1": "0", "dist": "t", "df": "3", "test": "mw"}],
    "sim.multiarm": [{"mu": "0 0 0", "n": "10", "adjust": "holm"}, {"mu": "0 0 0 0 0 0", "adjust": "dunnett"}],
    "sim.ordinal": [{"or": "1", "n1": "25"}, {"or": "1", "test": "lr", "p2": "0.7 0.2 0.1"}],
    "sim.winratio": [{"hr_d": "1", "hr_h": "1", "n1": "60"}, {"hr_d": "1", "hr_h": "1", "test": "wr", "n1": "150"}],
    "sim.mmrm": [{"delta": "0", "n1": "20", "drop1": "0.4"}, {"delta": "0", "visits": "10", "corr": "ar1"}],
    "sim.sw": [{"delta": "0", "steps": "6", "g": "2", "analysis": "hooper", "cac": "0.6"},
               {"delta": "0", "steps": "3", "g": "4", "m": "50"}],
}


def calibration(reps=20000):
    out = {}
    for pid in SIMS:
        rows = []
        for extra in NULLS.get(pid, []):
            inputs = dict(tp.base_inputs(pid), **extra, reps=str(reps))
            r = run(pid, "power", inputs)["rows"][0]
            alpha = float(inputs.get("alpha", 0.05))
            se = math.sqrt(alpha * (1 - alpha) / reps)
            rows.append({"inputs": extra, "type1": r.get("power"), "z": round((r.get("power") - alpha) / se, 2),
                         "error": r.get("_error")})
        out[pid] = rows
        print("calibration", pid, rows, flush=True)
    report["calibration"] = out


def bigreps():
    out = {}
    for pid, inputs in [("sim.means.robust", {"n1": "20"}), ("sim.multiarm", {"n": "20"}), ("sim.ordinal", {"n1": "40"}),
                        ("sim.lmm.cluster", {}), ("sim.surv.nph", {"n1": "50"})]:
        t = time.perf_counter()
        r = run(pid, "power", dict(tp.base_inputs(pid), **inputs, reps="1000000"))["rows"][0]
        out[pid] = {"seconds": round(time.perf_counter() - t, 1), "power": r.get("power"), "se": r.get("power_se"),
                    "maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)}
        print("bigreps", pid, out[pid], flush=True)
    report["bigreps"] = out


EXTREMES = {
    "sim.means.robust": [{"n1": "100000", "reps": "100"}, {"dist": "lognormal", "skew": "20", "n1": "5"},
                         {"dist": "t", "df": "2.01"}, {"dist": "contaminated", "eps": "0.999", "kout": "1000"},
                         {"sd1": "1e-8", "sd2": "1e8"}, {"mu1": "1e12"}],
    "sim.multiarm": [{"mu": "0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5", "n": "2"}, {"r0": "0.1", "n": "200"},
                     {"mu": "100"}, {"n": "100000", "reps": "100"}],
    "sim.ordinal": [{"p2": " ".join(["0.05"] * 20)}, {"p2": "0.999 0.001"}, {"or": "1e6"}, {"or": "1e-6"},
                    {"n1": "200000", "reps": "100"}, {"effect": "probs", "p1": "0 0 0 0 1"}],
    "sim.winratio": [{"n1": "2", "reps": "200"}, {"hr_d": "1e-3"}, {"med_d": "1e-6"}, {"med_d": "1e9"},
                     {"accrual": "0"}, {"eta": "100"}, {"tier3": "score", "margin_s": "1e9"}, {"n1": "2900", "reps": "100"}],
    "sim.mmrm": [{"visits": "12", "drop1": "0.999", "drop2": "0.999"}, {"rho": "0.999"}, {"rho": "-0.2", "visits": "6"},
                 {"rho": "-0.999", "corr": "ar1"}, {"n1": "3", "drop1": "0.9"}, {"n1": "50000", "reps": "100"}],
    "sim.sw": [{"steps": "10", "g": "1", "m": "1"}, {"icc": "0.999"}, {"icc": "0"}, {"cac": "0"},
               {"m": "100000", "reps": "100"}, {"g": "500", "reps": "100"}],
    "sim.lmm.cluster": [{"icc": "0.999"}, {"m": "1", "cv": "5"}, {"k1": "2"}],
    "sim.glmm.cluster": [{"p1": "0.999", "p2": "0.001"}, {"icc": "0.9"}, {"m": "1"}],
    "sim.lmm.longitudinal": [{"dropout": "0.999"}, {"rho": "-0.999"}, {"m": "2"}],
    "sim.surv.nph": [{"hr2": "1e-6"}, {"tchange": "1e9"}, {"eta": "100"}],
}


def extremes():
    out = {}
    for pid in SIMS:
        res_list = []
        for extra in EXTREMES.get(pid, []):
            for solve in tp.modes(pid):
                inputs = dict(tp.base_inputs(pid), **extra)
                inputs.setdefault("reps", "300")
                if "reps" not in extra:
                    inputs["reps"] = "300"
                t = time.perf_counter()
                res = run(pid, solve, inputs)
                dt = time.perf_counter() - t
                probs = tp.check_result(pid, solve, inputs, res, dt)
                row = (res.get("rows") or [{}])[0]
                res_list.append({"extra": extra, "solve": solve, "s": round(dt, 2),
                                 "outcome": row.get("_error") or res.get("error") or res.get("input_errors")
                                 or f"power {row.get('power')}", "problems": probs})
        out[pid] = res_list
        for x in res_list:
            print("extreme", pid, x, flush=True)
    report["extremes"] = out


PARTS = ("fuzz", "threads", "concurrency", "cancel", "calibration", "bigreps", "extremes")

if __name__ == "__main__":
    args = sys.argv[1:]
    out = None
    if "--json" in args:
        i = args.index("--json")
        out = args[i + 1]
        del args[i:i + 2]
    for part in args or PARTS:
        if part not in PARTS:
            sys.exit(f"unknown part {part!r}; choose from {', '.join(PARTS)}")
        t = time.perf_counter()
        globals()[part]()
        print(f"== {part} done in {time.perf_counter() - t:.0f} s", flush=True)
    if out:
        with open(out, "w") as f:
            json.dump(report, f, indent=1, default=str)
