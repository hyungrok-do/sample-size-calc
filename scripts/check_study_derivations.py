#!/usr/bin/env python3
"""Kernel-check generated algebraic derivations, including all-root/domain claims."""
import json
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    'quadratic': ('x^2+2*x', '8', '0,10'),
    'linear': ('2*x+1', '8', '-10,10'),
    'rational': ('x/3 + 0.1', '0.7', '0,10'),
    'both_sides': ('x*x+4', '5*x', '0,10'),
    'negative_coefficient': ('-2*x*x+4', '0', '-10,10'),
    'irrational': ('x*x', '2', '0,10'),
    'repeated_root': ('x*x-2*x+1', '0', '0,10'),
    'no_real_roots': ('x*x+2*x+2', '0', '-10,10'),
    'roots_outside_domain': ('x*x', '4', '3,10'),
    'every_value': ('x-x', '0', '-10,10'),
    'constant_false': ('1', '2', '-10,10'),
    'folded_constant': ('a*x + a', '1.4', '-10,10'),
}


def main():
    with tempfile.TemporaryDirectory(prefix='ssc-proof-check-') as directory:
        for name, (left, right, bounds) in CASES.items():
            source = 'def a : Real := 0.1 + 0.2\n' + f'solve for x : Real in [{bounds}] where {left} = {right}'
            result = subprocess.run([str(ROOT/'build/ssc_cli'), 'design', '-'],
                                    input=json.dumps({'op': 'compile', 'source': source}), text=True,
                                    capture_output=True, check=True, timeout=30)
            derivation = json.loads(result.stdout)['result']['derivation']
            assert derivation['kind'] == 'exact_algebraic', (name, derivation)
            proof = Path(directory)/(name+'.lean')
            proof.write_text(derivation['leanSource'])
            checked = subprocess.run(['lake', 'env', 'lean', str(proof)], cwd=ROOT/'lean',
                                     text=True, capture_output=True, timeout=120)
            assert checked.returncode == 0, (name, checked.stdout, checked.stderr)
            assert 'sorryAx' not in checked.stdout + checked.stderr
            audits = re.findall(r"depends on axioms: \[([^\]]*)\]", checked.stdout)
            assert len(audits) >= 2, (name, checked.stdout)
            for audit in audits:
                assert {a.strip() for a in audit.split(',') if a.strip()} <= {'propext', 'Quot.sound', 'Classical.choice'}
            print(f'{name}: Lean kernel checked ({len(audits)} audited theorems)', flush=True)


if __name__ == '__main__':
    main()
