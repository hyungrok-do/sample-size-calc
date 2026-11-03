#!/usr/bin/env python3
"""Compare typed-plan execution, emitted native C++, and independent references."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile

from scipy import integrate, optimize, stats


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine', type=Path, default=root / 'build/ssc_cli')
    parser.add_argument('--cxx', default='c++')
    args = parser.parse_args()

    def engine(request):
        result = subprocess.run([str(args.engine.resolve()), 'design', '-'], input=json.dumps(request),
                                text=True, capture_output=True, timeout=60, check=True)
        value = json.loads(result.stdout)
        assert value['ok'], value
        return value['result']

    def close(a, b):
        assert abs(a-b) <= 2e-10 * max(1, abs(b)), (a, b)

    sources = {p.stem: p.read_text() for p in (root / 'examples/programs').glob('*.study')}
    # Shared stochastic arguments and left-to-right evaluation must survive code generation.
    sources['shared-draw'] = '''
def event (x : Real) : Bool := x-x = 0
simulation { method := simulation, replicates := 1000, seed := 21 }
evaluate probability(event(normal(0,1)))
'''
    sources['array-sharing'] = '''
def event (n : Nat) : Bool :=
  let a := sample(n, fun i => exp(normal(0,1)))
  mean(a) = mean(a)
simulation { replicates := 1000, seed := 21 }
evaluate probability(event(5))
'''
    sources['argument-order'] = '''
def event (x : Real) (y : Real) : Bool := x > y
simulation { method := simulation, replicates := 1000, seed := 21 }
evaluate probability(event(normal(0,1), normal(1,2)))
'''
    sources['real-literals'] = 'evaluate 50000 * 50000 + 2000000000 + 2000000000'
    results = {}
    with tempfile.TemporaryDirectory(prefix='ssc-compiler-check-') as temp:
        for name, source in sources.items():
            compiled = engine({'op': 'compile', 'source': source})
            expected = engine(compiled['plan'])
            cpp, executable = Path(temp) / (name+'.cpp'), Path(temp) / name
            cpp.write_text(compiled['cpp'])
            subprocess.run([args.cxx, '-std=c++17', '-O2', str(cpp), '-o', str(executable)], check=True, timeout=120)
            actual = json.loads(subprocess.check_output([str(executable)], text=True, timeout=120))
            for key, value in actual.items():
                if isinstance(value, (int, float)):
                    close(value, expected[key])
                else:
                    assert value == expected[key], (name, key, value, expected[key])
            results[name] = actual
            print(name, actual['method'], actual['status'], actual.get('solution', actual.get('value')))

    z = stats.norm.isf(.025)
    power = lambda n, d: stats.norm.sf(z-d*(n/2)**.5) + stats.norm.sf(z+d*(n/2)**.5)
    assert results['normal-power']['solution'] == next(n for n in range(2, 10000) if power(n, .5) >= .8)
    close(results['normal-power']['value'], power(63, .5))
    effect = optimize.brentq(lambda d: power(63, d)-.8, .01, 2)
    assert abs(results['detectable-effect']['solution'] - effect) < 2e-8
    assert results['equation']['solution'] == 2
    area_root = optimize.brentq(lambda x: integrate.quad(lambda t: math.exp(-t*t), 0, x)[0]-.7, 0, 3)
    assert abs(results['integration']['solution'] - area_root) < 1e-7
    estimate = results['evaluate-simulation']
    probability = stats.lognorm.sf(1.8, s=1)
    assert estimate['lower'] < probability < estimate['upper']
    assert results['shared-draw']['value'] == results['array-sharing']['value'] == 1
    print('All native generated programs agree with the typed-plan runtime and independent references.')


if __name__ == '__main__':
    main()
