#!/usr/bin/env python3
"""Derive a Lean-style study, export a Lean proof and/or standalone C++17.

--check-lean invokes Lean's proof kernel for supported exact algebraic derivations.
Example: python scripts/compile_study.py examples/programs/normal-power.study --cpp /tmp/study.cpp --run
"""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--cpp', type=Path)
    parser.add_argument('--lean', type=Path, help='Save the generated exact algebraic derivation as a Lean proof')
    parser.add_argument('--check-lean', action='store_true', help='Check the generated proof with the pinned Lean/Mathlib environment')
    parser.add_argument('--engine', type=Path, default=root / 'build/ssc_cli')
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--cxx', default='c++')
    args = parser.parse_args()
    if not args.cpp and not args.lean:
        parser.error('Supply --lean and/or --cpp')
    if args.run and not args.cpp:
        parser.error('--run requires --cpp')
    if args.check_lean and not args.lean:
        parser.error('--check-lean requires --lean')
    request = {'op': 'compile', 'source': args.source.read_text(encoding='utf-8')}
    response = subprocess.run([str(args.engine.resolve()), 'design', '-'], input=json.dumps(request), text=True, capture_output=True, timeout=60)
    answer = json.loads(response.stdout)
    if response.returncode or not answer.get('ok'):
        raise SystemExit(answer.get('error', response.stderr))
    compiled = answer['result']
    if args.cpp:
        args.cpp.write_text(compiled['cpp'], encoding='utf-8')
    proof_checked = False
    if args.lean:
        proof = compiled.get('derivation', {}).get('leanSource')
        if not proof:
            raise SystemExit('No exact algebraic Lean proof is available for this program; its derivation describes the numerical or simulation method.')
        args.lean.write_text(proof, encoding='utf-8')
        if args.check_lean:
            checked = subprocess.run(['lake', 'env', 'lean', str(args.lean.resolve())], cwd=root / 'lean',
                                     capture_output=True, text=True, timeout=120)
            if checked.returncode or 'sorryAx' in checked.stdout or 'sorry' in checked.stderr:
                raise SystemExit(checked.stdout + checked.stderr)
            proof_checked = True
    if args.run:
        with tempfile.TemporaryDirectory(prefix='ssc-study-build-') as temporary:
            executable = Path(temporary) / 'study'
            subprocess.run([args.cxx, '-std=c++17', '-O2', str(args.cpp.resolve()), '-o', str(executable)], check=True, timeout=120)
            subprocess.run([str(executable)], check=True, timeout=120)
    else:
        print(json.dumps({'cpp': str(args.cpp) if args.cpp else None, 'lean': str(args.lean) if args.lean else None,
                          'compiler': compiled['compiler'], 'leanDerivationKernelChecked': proof_checked,
                          'scope': 'Exact real-number algebra only; not binary64 execution or stochastic-model validation.',
                          'derivation': compiled.get('derivation', {}).get('steps')}, indent=2))


if __name__ == '__main__':
    main()
