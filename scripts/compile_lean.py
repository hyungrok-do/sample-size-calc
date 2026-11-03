#!/usr/bin/env python3
"""Check a study with Lean, compare its exported plan, and emit standalone C++.

Usage: python scripts/compile_lean.py examples/design/normal.lean --cpp study.cpp
Requires `lake build` in lean/ and a native C++ build. Runs the supplied Lean
program on this computer. The browser compiler accepts the documented subset.
"""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--cpp", type=Path, required=True)
    parser.add_argument("--solve", choices=["sample_size", "power", "effect"], default="sample_size")
    parser.add_argument("--n", type=int)
    parser.add_argument("--engine", type=Path, default=root / "build" / "ssc_cli")
    args = parser.parse_args()
    source = args.source.read_text(encoding="utf-8")
    request = {"op": "compile", "source": source, "solveFor": args.solve}
    if args.n is not None:
        request["n"] = args.n
    # Reject unsupported syntax before evaluating the local source in Lean.
    native = subprocess.run([str(args.engine.resolve()), "design", "-"], input=json.dumps(request),
                            capture_output=True, text=True, check=False, timeout=60)
    if native.returncode:
        raise SystemExit(native.stdout or native.stderr)
    compiled = json.loads(native.stdout)["result"]
    with tempfile.TemporaryDirectory(prefix="ssc-study-") as temp:
        path = Path(temp) / "Study.lean"
        path.write_text(source + "\n\ndef main : IO Unit := IO.println (DesignCompiler.compileInput study).compress\n", encoding="utf-8")
        result = subprocess.run(["lake", "env", "lean", "--run", str(path)], cwd=root / "lean",
                                capture_output=True, text=True, check=False, timeout=120)
        if result.returncode:
            raise SystemExit(result.stdout + result.stderr)
        checked = json.loads(result.stdout.strip())
    expected = dict(compiled["plan"])
    expected.pop("solve", None)
    if checked != expected:
        raise SystemExit("Lean and C++ compiler plans differ; no source emitted.")
    args.cpp.write_text(compiled["cpp"], encoding="utf-8")
    print(json.dumps({"leanKernelChecked": True, "compilerPlansAgree": True,
                      "cpp": str(args.cpp), "plan": compiled["plan"]}, indent=2))


if __name__ == "__main__":
    main()
