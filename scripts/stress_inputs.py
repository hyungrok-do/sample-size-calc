"""What users can type, in every input of every method: the long sweep behind tests/test_extremes.cpp.

    PYTHONPATH=python python scripts/stress_inputs.py [procedure ...] [--limit 5] [--json report.json]

Every number or list input, in every solve mode where it shows, gets a battery of typed texts (blank, letters,
decimal commas, exponents, 1e-300 and 1e300 (refused) and the sizes allowed nearest them, inf and nan, a unicode
minus, ranges good and bad, a list of 3,000 values ...), run the way an interface runs it (the quick approximation
for a simulation, power curve on).  Reported:

  hang      still running after --limit seconds (default 5), and whether Stop (cancelling the token) ends it
  message   a failure without a message, or with an internal one
  blank     an answer or a power that is no number, without a message
  infinite  an infinite value in a row or a chart

Each procedure runs in a process of its own; a run that Stop cannot end ends that process, and the sweep goes
on after it.  The whole sweep takes about an hour on four cores.
"""
import json
import math
import os
import subprocess
import sys
import threading
import time

import ssc
import ssc.core as core

NUMBER = [
    "", " ", "abc", "1,5", "0,8", "1.5.2", "-", ".", "0.", ".5", "5.", "1e3", "1e-3", "1E2", "1e309", "-1e309",
    "inf", "-inf", "Infinity", "nan", "NaN", "-0", "0", "-1", "-0.5", "1e-300", "1e300", "1e-30", "1e30", "-1e30",
    "999999999999", "0.0000001", "1 to 10", "1 to 10 by 0", "10 to 1 by 1", "1 to 1000000 by 1",
    "1 to 2 by 0.0000001",
    "0.5 0.6 0.7", "1;2", "1|2", "−1", "1 000", "10%", "½", "１", "0x10", "1e", "e5", "+5", "++5",
    "5-", "1 to", "to 5", "by 2", "1 by 2", "1 to 5 by", "1 2 to 5", "  7  ", "\t3", "3\n", "1 " * 3000,
    "0.99999999", "1.0000001", "2147483648", "4294967296", "1e10", "-1e10", "0.5 to 0.9 by 0.1",
]
VECTOR = [
    "", " ", "1", "0", "-1", "a b", "1;2", "0.6;0.5 0.4", "1 2 3 4 5 6 7 8 9", "0.5 0.5", "0.2 0.3",
    "-0.1 0.6 0.5", "1 to 5", "0 0 0", "1e309 1", "nan 1", "0.1 " * 50, "1 2 3 " * 400, "0,5 0,5", "1|2|3",
    "0.25 0.25 0.25 0.25", "1 1 1 1 1 1 1 1", "0.001 0.999", "100 200", "−0.5 0.5", "1e-30 1", "0 1e30",
]
INTERNAL = ("missing input", "missing selection", "missing list input", "not bracketed", "internal", "std::",
            "unknown option", "out of range", "exception", "assert")


def visible(spec, inputs, solve):
    """The schema's visibility rule (any clause list whose clauses all hold)."""
    if not spec.get("visible"):
        return True
    return any(all(((solve if c["key"] == "$solve" else inputs.get(c["key"], "")) in c.get("values", []))
                   != bool(c.get("negate")) for c in clauses) for clauses in spec["visible"])


def cases(pid):
    schema = ssc.schema(pid)
    base = {p["id"]: str(p["default"]) for p in schema["params"]}
    for spec in schema["params"]:
        if spec["type"] not in ("number", "vector"):
            continue
        for o in schema["solve_for"]:
            if visible(spec, base, o["value"]):
                for text in NUMBER if spec["type"] == "number" else VECTOR:
                    yield o["value"], o.get("ykey") or o["value"], spec["id"], text, dict(base, **{spec["id"]: text})


def problems(res):
    """What is wrong with a finished result, as (kind, text) pairs."""
    out = []
    if not res.get("ok"):
        text = " ".join([res.get("error") or ""] + [str(v) for v in (res.get("input_errors") or {}).values()])
        if not text.strip():
            out.append(("message", "failed without a message"))
        elif any(s in text.lower() for s in INTERNAL):
            out.append(("message", text.strip()[:160]))
        return out
    for row in res.get("rows", []):
        if row.get("_error"):
            text = row["_error"] + " " + row.get("_detail", "")
            if any(s in text.lower() for s in INTERNAL):
                out.append(("message", text.strip()[:160]))
            continue
        for key, value in row.items():
            if isinstance(value, float) and math.isinf(value):
                out.append(("infinite", f"{key} = {value}"))
    for plot in res.get("plots") or []:
        for series in plot.get("series") or []:
            if any(isinstance(v, float) and math.isinf(v) for key in ("x", "y") for v in series.get(key) or []):
                out.append(("infinite", f"chart {plot.get('kind')}"))
    return out


def child(pid, start, limit):
    """Runs the cases of one procedure from `start` on, one JSON line per finding."""
    sim = bool(ssc.schema(pid).get("simulation"))
    i = -1
    for i, (solve, answer, key, text, inputs) in enumerate(cases(pid)):
        if i < start:
            continue
        where = {"i": i, "procedure": pid, "solve": solve, "input": key, "text": text[:60]}
        token, box = core.Token(), {}
        worker = threading.Thread(target=lambda: box.setdefault("res", core.run(
            pid, solve, inputs, threads=0, power_curve=True, curve_points=41, approximate=sim, token=token)),
            daemon=True)
        worker.start()
        worker.join(limit)
        if worker.is_alive():
            t = time.perf_counter()
            token.cancel()
            worker.join(15)
            stopped = not worker.is_alive()
            print(json.dumps(dict(where, kind="hang", stops=stopped, stop_s=round(time.perf_counter() - t, 2))),
                  flush=True)
            if not stopped:
                os._exit(3)  # the run cannot be ended: the driver resumes after it
            continue
        token.close()
        res = box["res"]
        found = problems(res)
        if res.get("ok"):
            for row in res.get("rows", []):
                for k in ("power", answer):
                    if not row.get("_error") and k in row and row[k] is None:
                        found.append(("blank", f"{k} is no number"))
        for kind, msg in found:
            print(json.dumps(dict(where, kind=kind, message=msg)), flush=True)
    print(json.dumps({"done": True, "runs": i + 1}), flush=True)


def main():
    args = sys.argv[1:]
    limit = float(args[args.index("--limit") + 1]) if "--limit" in args else 5.0
    out = args[args.index("--json") + 1] if "--json" in args else None
    skip = {"--limit", "--json"}
    ids = [a for j, a in enumerate(args) if a not in skip and (j == 0 or args[j - 1] not in skip)]
    ids = ids or [p["id"] for p in ssc.catalog()["procedures"]]
    report, total = [], 0
    for pid in ids:
        t, start, runs = time.perf_counter(), 0, 0
        found = []
        while True:
            proc = subprocess.run([sys.executable, __file__, "--child", pid, str(start), str(limit)],
                                  stdout=subprocess.PIPE, text=True)
            lines = [json.loads(line) for line in proc.stdout.splitlines() if line.startswith("{")]
            found += [x for x in lines if "kind" in x]
            runs = next((x["runs"] for x in lines if x.get("done")), runs)
            if proc.returncode == 3 and found:  # an uncancellable run: go on after it
                start = found[-1]["i"] + 1
                continue
            if proc.returncode not in (0, 3):
                found.append({"procedure": pid, "kind": "crash", "message": f"exit code {proc.returncode}"})
            break
        report += found
        total += runs
        print(f"{pid:24s} {runs:5d} runs  {len(found):3d} findings  ({time.perf_counter() - t:.0f} s)", flush=True)
        for x in found[:10]:
            print(f"    {x['kind']:8s} [{x.get('solve', '')}] {x.get('input', '')}={x.get('text', '')!r}: "
                  f"{x.get('message') or ('stops' if x.get('stops') else 'DOES NOT STOP')}", flush=True)
    print(f"{total} runs, {len(report)} findings", flush=True)
    if out:
        with open(out, "w") as fh:
            json.dump(report, fh, indent=1, ensure_ascii=False)
    sys.exit(1 if any(x["kind"] in ("crash", "message", "blank", "infinite") or
                      (x["kind"] == "hang" and not x.get("stops")) for x in report) else 0)


if __name__ == "__main__":
    if sys.argv[1:2] == ["--child"]:
        child(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]))
    else:
        main()
