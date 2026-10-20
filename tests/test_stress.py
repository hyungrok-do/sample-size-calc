"""Stress tests of the simulation procedures: any thread count gives the same answer, many simulations at once
give the same answers as one after another, a running simulation stops promptly when cancelled, every test
keeps its type I error at a null design, and a long simulation keeps its memory bounded and gives it back.

The longer sweep (fuzzed and extreme inputs, a million simulated studies, sanitizer builds) is described in
README.md; these are the parts quick enough for every test run.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import threading
import time

import pytest

import ssc
import ssc.core as core

SIMS = [p["id"] for p in ssc.catalog()["procedures"] if ssc.schema(p["id"]).get("simulation")]


def example(pid: str, reps: str = "200") -> dict:
    schema = ssc.schema(pid)
    inputs = {p["id"]: str(p["default"]) for p in schema["params"]}
    inputs.update({k: str(v) for k, v in schema["example"]["inputs"].items()})
    inputs["reps"] = reps
    return inputs


def rows(res: dict) -> str:
    return json.dumps(res.get("rows"), sort_keys=True)


def test_every_simulation_is_covered():
    assert len(SIMS) >= 10


@pytest.mark.parametrize("pid", SIMS)
def test_any_thread_count_gives_the_same_answer(pid):
    solve = ssc.schema(pid)["example"]["solve_for"]
    inputs = example(pid)
    reference = rows(core.run(pid, solve, inputs, threads=1, power_curve=False))
    for threads in (2, 3, 7, 16):
        assert rows(core.run(pid, solve, inputs, threads=threads, power_curve=False)) == reference, threads


def test_simulations_running_at_once_match_one_after_another():
    jobs = [(pid, "power", example(pid)) for pid in SIMS]
    sequential = [rows(core.run(*job, power_curve=False)) for job in jobs]
    results = [None] * len(jobs)

    def work(i):
        results[i] = rows(core.run(*jobs[i], power_curve=False))

    workers = [threading.Thread(target=work, args=(i,)) for i in range(len(jobs))]
    for w in workers:
        w.start()
    for w in workers:
        w.join()
    assert results == sequential


@pytest.mark.parametrize("pid", SIMS)
def test_a_running_simulation_stops_promptly_when_cancelled(pid):
    token = core.Token()
    box = {}
    worker = threading.Thread(target=lambda: box.setdefault(
        "res", core.run(pid, "power", example(pid, reps="1000000"), power_curve=False, token=token)))
    worker.start()
    time.sleep(0.2)
    start = time.perf_counter()
    token.cancel()
    worker.join(10)
    latency = time.perf_counter() - start
    token.close()
    assert not worker.is_alive()
    assert box["res"].get("cancelled"), box["res"]
    assert latency < 2, f"{pid} took {latency:.2f} s to stop"


NULLS = [
    ("sim.lmm.cluster", {"mu1": "0", "k1": "6", "cv": "0.6"}),
    ("sim.glmm.cluster", {"p1": "0.2", "k1": "15", "m": "20"}),
    ("sim.lmm.longitudinal", {"delta": "0", "n1": "20", "dropout": "0.3"}),
    ("sim.surv.nph", {"hr1": "1", "hr2": "1", "test": "maxcombo"}),
    ("sim.means.robust", {"mu1": "0", "n1": "20", "dist": "t", "df": "5", "sd2": "4", "r": "2"}),
    ("sim.multiarm", {"mu": "0 0 0", "n": "10", "adjust": "holm"}),
    ("sim.ordinal", {"or": "1", "n1": "25", "test": "lr"}),
    ("sim.winratio", {"hr_d": "1", "hr_h": "1", "n1": "60"}),
    ("sim.mmrm", {"delta": "0", "n1": "25", "drop1": "0.4"}),
    ("sim.sw", {"delta": "0", "steps": "6", "g": "2", "analysis": "hooper", "cac": "0.6"}),
]


@pytest.mark.parametrize("pid,null", NULLS, ids=[n[0] for n in NULLS])
def test_type_one_error_at_a_null_design(pid, null):
    reps = 6000
    inputs = dict(example(pid), **null, reps=str(reps))
    row = core.run(pid, "power", inputs, power_curve=False)["rows"][0]
    se = math.sqrt(0.05 * 0.95 / reps)
    # small designs may run a little above α (Wald tests with few clusters, Satterthwaite with few subjects),
    # never far from it
    assert -5 * se <= row["power"] - 0.05 <= 5 * se + 0.012, row["power"]


MEMORY = """
import json, sys
import ssc.core as core

def mem():
    d = {}
    for line in open("/proc/self/status"):
        if line.startswith(("VmRSS", "VmHWM")):
            d[line.split(":")[0]] = int(line.split()[1]) / 1024
    return d

core.run("means.two.t", "n", {})
before = mem()
res = core.run("sim.multiarm", "power", {"reps": "400000"}, power_curve=False)
print(json.dumps({"ok": res["ok"], "before": before, "after": mem()}))
"""


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="reads /proc/self/status")
def test_a_long_simulation_keeps_memory_bounded_and_gives_it_back():
    # in a process of its own, so that the peak is this simulation's; 400,000 studies held at once would take
    # over 20 MB, the engine holds a block of 16,384 at a time and returns freed memory to the system
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([os.path.dirname(os.path.dirname(core.__file__)),
                                                       os.environ.get("PYTHONPATH", "")]))
    out = subprocess.run([sys.executable, "-c", MEMORY], capture_output=True, text=True, env=env, timeout=600)
    assert out.returncode == 0, out.stderr
    d = json.loads(out.stdout.strip().splitlines()[-1])
    assert d["ok"]
    assert d["after"]["VmHWM"] - d["before"]["VmHWM"] < 10, d
    assert d["after"]["VmRSS"] - d["before"]["VmRSS"] < 5, d
