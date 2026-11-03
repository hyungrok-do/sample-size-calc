#!/usr/bin/env python3
"""Check a study with Lean, compare its exported plan, and emit standalone C++.

Usage: python scripts/compile_lean.py examples/design/normal.lean --cpp study.cpp
Requires `lake build` in lean/ and a native C++ build. Runs the supplied Lean
program on this computer. The browser compiler accepts the documented subset.
"""
import argparse
import json
import re
from pathlib import Path
import subprocess
import tempfile


def native_source(source):
    """Supply the implicit module, preserving the already-validated declaration."""
    position = 0
    while position < len(source):
        if source[position].isspace():
            position += 1
        elif source.startswith('--', position):
            end = source.find('\n', position)
            position = len(source) if end < 0 else end
        elif source.startswith('/-', position):
            depth = 1
            position += 2
            while depth:
                if source.startswith('/-', position):
                    depth += 1
                    position += 2
                elif source.startswith('-/', position):
                    depth -= 1
                    position += 2
                else:
                    position += 1
        else:
            token = re.match(r'[A-Za-z_.]+', source[position:])
            if token is None:
                raise ValueError('Expected a validated study declaration')
            if token.group() == 'def':
                return 'import DesignCompiler.Input\nopen DesignCompiler\n' + source[position:]
            position += len(token.group())
    raise ValueError('Missing study declaration')


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
        path.write_text(native_source(source) + "\n\ndef main : IO Unit := IO.println (DesignCompiler.compileInput study).compress\n", encoding="utf-8")
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
