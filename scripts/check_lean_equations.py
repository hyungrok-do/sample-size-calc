#!/usr/bin/env python3
"""Integration checks against the real Lean kernel, emitted C++, and local API."""
import concurrent.futures
import http.client
import json
import math
from pathlib import Path
import subprocess
import tempfile
import threading
import time

from lean_equations import ROOT, compile_source
from serve_lean import CompilerServer

ALLOWED = {'propext', 'Classical.choice', 'Quot.sound'}


def compile_ok(source, status='solved', **kwargs):
    value = compile_source(source, **kwargs)
    assert value['ok'], value.get('error')
    assert value['leanKernelChecked'] is True
    for equation in value['results']:
        assert equation['status'] == status, equation['solution']
        assert set(equation['axioms']) <= ALLOWED
    return value['results'][0]


def cpp_run(equation, *arguments, succeeds=True):
    assert equation['cpp'], equation['lowering']['reason']
    with tempfile.TemporaryDirectory() as directory:
        path, exe = Path(directory)/'equation.cpp', Path(directory)/'equation'
        path.write_text(equation['cpp'])
        subprocess.run(['c++', '-std=c++17', '-O2', str(path), '-o', str(exe)], check=True, timeout=30)
        process = subprocess.run([str(exe), *map(str, arguments)], capture_output=True, text=True, timeout=10)
        if not succeeds:
            assert process.returncode != 0
            return process.stderr
        assert process.returncode == 0, process.stderr
        return json.loads(process.stdout)


def main():
    for unfinished in [
        'theorem invalid : False := by sorry',
        'set_option warningAsError false in\ntheorem invalid : False := by sorry',
        'axiom invalid : False',
        '''open Lean Elab Tactic
elab "bogusProof" : tactic => do
  let goal ← getMainGoal
  goal.assign (mkConst ``True.intro)
  replaceMainGoal []
set_option debug.skipKernelTC true in
theorem invalid : False := by bogusProof''',
    ]:
        result = compile_source(unfinished + '\nsolve for x : ℝ where x = 1')
        assert not result['ok'], result
    examples = {}
    for path in sorted((ROOT/'examples/lean-equations').glob('*.lean')):
        examples[path.stem] = compile_ok(path.read_text(), 'unresolved' if path.stem == 'numerical' else 'solved')
    assert examples['parameters']['solution'] == 'x = (Real.log y - b) / a'
    assert '∀ (a b y : ℝ)' in examples['parameters']['theoremType']
    assert 'a ≠ 0' in examples['parameters']['theoremType'] and '0 < y' in examples['parameters']['theoremType']
    assert any(s['rule'] == 'cubeLink_zero' for s in examples['custom-rule']['steps'])
    for name, expected in [('composition', [(math.log(5)-1)/2]), ('branches', [-1, 5/3]),
                           ('custom-rule', [1.5]), ('quadratic', [-4, 2])]:
        roots = cpp_run(examples[name])['roots']
        assert len(roots) == len(expected)
        assert all(math.isclose(a, b, rel_tol=1e-13, abs_tol=1e-13) for a, b in zip(roots, expected))
    value = cpp_run(examples['numerical'], 0, 1)
    assert abs(math.cos(value['solution']) - value['solution']) <= 1e-10
    assert value['kernelCertified'] is False
    assert 'No sign-changing bracket' in cpp_run(examples['numerical'], 2, 3, succeeds=False)
    assert 'Invalid finite bound' in cpp_run(examples['numerical'], 'nan', 1, succeeds=False)
    assert cpp_run(compile_ok('solve for x : ℝ where Real.exp (x+1) = -2'))['roots'] == []
    assert cpp_run(compile_ok('solve for x : ℝ where 2*x^2 + 3*x = 5'))['roots'] == [-2.5, 1]
    assert compile_ok('solve for x : ℝ where 3*x + 1 = x + 7')['solution'] == 'x = 3'
    conditional = compile_ok('variable (a b y : ℝ)\nsolve for x : ℝ where Real.exp (a*x+b)=y')
    assert '0 < y' in conditional['solution'] and 'a = 0' in conditional['solution']
    assert conditional['lowering']['kind'] == 'analytic_roots'
    assert cpp_run(conditional, 'a=0', 'b=0', 'y=1')['solutionSet'] == 'all_reals'
    assert cpp_run(conditional, 'a=0', 'b=0', 'y=2')['solutionSet'] == 'empty'
    assert cpp_run(conditional, 'a=2', 'b=1', 'y=-1')['solutionSet'] == 'empty'
    # Reusable parameterized code, including every hypothesis of its theorem.
    parameterized = examples['parameters']
    assert [p['name'] for p in parameterized['lowering']['parameters']] == ['a', 'b', 'y']
    assert [h['statement'] for h in parameterized['lowering']['assumptions']] == ['a ≠ 0', '0 < y']
    result = cpp_run(parameterized, 'y=5', 'a=2', 'b=1')
    assert math.isclose(result['roots'][0], (math.log(5)-1)/2, rel_tol=1e-13)
    assert 'Assumption failed' in cpp_run(parameterized, 'a=0', 'b=1', 'y=5', succeeds=False)
    assert 'Assumption failed' in cpp_run(parameterized, 'a=2', 'b=1', 'y=-5', succeeds=False)
    for arguments, message in [
        (['a=2', 'b=1'], 'Missing parameter'),
        (['a=2', 'b=1', 'y=5', 'a=3'], 'Duplicate parameter'),
        (['a=2', 'b=1', 'y=5', 'z=3'], 'Unknown parameter'),
        (['a=nan', 'b=1', 'y=5'], 'Invalid finite'),
        (['a=1e-9999', 'b=1', 'y=5'], 'Invalid finite'),
    ]:
        assert message in cpp_run(parameterized, *arguments, succeeds=False)
    assert examples['precision']['solution'] == 'n = variance / targetVariance'
    assert cpp_run(examples['precision'], 'variance=2', 'targetVariance=0.01')['roots'] == [200]
    assert 'Assumption failed' in cpp_run(examples['precision'], 'variance=2', 'targetVariance=0', succeeds=False)
    # Same compiled conditional expression gives finite, empty, or all-real sets.
    for arguments, kind, roots in [
        (['a=2', 'y=6'], 'finite', [3]),
        (['a=-2', 'y=6'], 'finite', [-3]),
        (['a=0', 'y=6'], 'empty', []),
        (['a=0', 'y=0'], 'all_reals', []),
    ]:
        result = cpp_run(examples['conditional'], *arguments)
        assert result['solutionSet'] == kind and result['roots'] == roots
    square = compile_ok('variable (y : ℝ)\nsolve for x : ℝ where x^2=y')
    assert cpp_run(square, 'y=9')['roots'] == [-3, 3]
    assert cpp_run(square, 'y=0')['roots'] == [0]
    assert cpp_run(square, 'y=-1')['solutionSet'] == 'empty'
    reciprocal = compile_ok('variable (a b : ℝ)\nsolve for x : ℝ where a/x=b')
    assert cpp_run(reciprocal, 'a=6', 'b=2')['roots'] == [3]
    assert cpp_run(reciprocal, 'a=6', 'b=0')['roots'] == [0]  # Lean's totalized division.
    assert cpp_run(reciprocal, 'a=0', 'b=0')['solutionSet'] == 'all_reals'
    assert cpp_run(reciprocal, 'a=0', 'b=2')['solutionSet'] == 'empty'
    intersection = compile_ok('variable (a b : ℝ)\nsolve for x : ℝ where x=a ∧ x=b')
    assert cpp_run(intersection, 'a=2', 'b=2')['roots'] == [2]
    assert cpp_run(intersection, 'a=2', 'b=3')['solutionSet'] == 'empty'
    assert cpp_run(compile_ok('solve for x : ℝ where x=x'))['solutionSet'] == 'all_reals'
    unicode = compile_ok('variable (σ «class» : ℝ)\nsolve for x : ℝ where x=σ + «class»')
    assert cpp_run(unicode, 'σ=2', 'class=3')['roots'] == [5]
    numerical = compile_ok('variable (a : ℝ) (ha : 0 < a)\nsolve for x : ℝ where Real.cos x=a*x', 'unresolved')
    answer = cpp_run(numerical, 0, 1, 'a=2')['solution']
    assert abs(math.cos(answer)-2*answer) <= 1e-10
    assert 'Assumption failed' in cpp_run(numerical, 0, 1, 'a=-1', succeeds=False)
    unsupported_guard = compile_ok('variable (a : ℝ) (h : ∃ y : ℝ, y=a)\nsolve for x : ℝ where x=a')
    assert unsupported_guard['lowering']['kind'] == 'unavailable'
    assert 'No C++ guard' in unsupported_guard['lowering']['reason']
    order_override = compile_ok('local instance : LE ℝ := ⟨fun a b => b<a⟩\nvariable (a : ℝ) (h : 0 ≤ a)\nsolve for x : ℝ where x=a')
    assert order_override['lowering']['kind'] == 'unavailable'
    assert 'standard real arithmetic' in order_override['lowering']['reason']
    overflow = compile_ok('variable (a : ℝ)\nsolve for x : ℝ where x=Real.log (Real.exp a + 1)')
    assert 'Non-finite' in cpp_run(overflow, 'a=1000', succeeds=False)
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory)/'result.json'
        subprocess.run(['python', str(ROOT/'scripts/lean_equations.py'),
                        str(ROOT/'examples/lean-equations/precision.lean'), '--run',
                        '--param', 'variance=2', '--param', 'targetVariance=0.01',
                        '--output', str(output)], check=True, timeout=90)
        assert json.loads(output.read_text())['numericalEvaluation']['roots'] == [200]
    # Actual Lean structures, namespaces, custom notation, and subsequent proof reuse.
    module = '''namespace UserModel
structure Model where
  slope : ℝ
  offset : ℝ
noncomputable def model : Model := ⟨3, 1⟩
local notation "response" => fun x : ℝ => model.slope*x + model.offset
solve answer for x : ℝ where response x = 7
example (x : ℝ) (h : response x = 7) : x = 2 := (answer x).mp h
end UserModel
'''
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory)/'Checked.olean'
        assert compile_ok(module, olean=path)['solution'] == 'x = 2'
        assert path.stat().st_size > 0
    exported = compile_source('solve answer for x : ℝ where 2*x=4')['module']
    assert compile_ok(exported)['solution'] == 'x = 2'
    for invalid in [
        'solve for x : ℝ where missingFunction x = 3',
        '@[equation_rule] theorem bad (x : ℝ) : Real.sin x = 7 ↔ x = 1 := by sorry\nsolve for x : ℝ where Real.sin x = 7',
        'axiom bad : ∀ x : ℝ, Real.sin x = 7 ↔ x = 1\n@[equation_rule] theorem fake (x : ℝ) : Real.sin x = 7 ↔ x = 1 := bad x\nsolve for x : ℝ where Real.sin x = 7',
        'solve for x : ℝ where x + true = 3',
        'def hello := 3',
    ]:
        value = compile_source(invalid)
        assert not value['ok'], invalid
    unsupported = compile_ok('solve for x : ℝ where Real.sin x = x^2', 'unresolved')
    assert unsupported['lowering']['kind'] == 'bracketed_bisection'
    custom_arithmetic = compile_ok('local instance : HAdd ℝ ℝ ℝ := ⟨fun a b => a-b⟩\nsolve for x : ℝ where Real.sin (x+1) = x', 'unresolved')
    assert custom_arithmetic['lowering']['kind'] == 'unavailable'
    assert 'standard real arithmetic' in custom_arithmetic['lowering']['reason']
    # Lean's totalized log/sqrt/division are preserved by the numeric backend.
    total = compile_ok('solve for x : ℝ where x = Real.log (-2) + Real.sqrt (-3) + 1/0')
    assert math.isclose(cpp_run(total)['roots'][0], math.log(2), rel_tol=1e-13)
    slow = '#eval IO.sleep 20000\nsolve for x : ℝ where x=1'
    assert 'exceeded' in compile_source(slow, timeout=0.1)['error']
    server = CompilerServer(0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    def request(method, path, body=None, headers=None):
        connection = http.client.HTTPConnection('127.0.0.1', server.server_port, timeout=70)
        encoded = json.dumps(body) if body is not None else None
        connection.request(method, path, encoded, headers or {'Content-Type': 'application/json'})
        response = connection.getresponse()
        result = response.status, json.loads(response.read())
        connection.close()
        return result

    try:
        assert request('GET', '/api/lean/health')[1]['ready']
        assert request('POST', '/api/lean/compile', {}, {'Origin': 'https://example.org', 'Content-Type': 'application/json'})[0] == 403
        assert request('GET', '/api/lean/health', headers={'Host': 'attacker.example'})[0] == 403
        assert request('DELETE', '/api/lean/compile/early-cancel')[1]['cancelled']
        assert request('POST', '/api/lean/compile', {'source': slow, 'jobId': 'early-cancel'})[1]['error'] == 'Compilation stopped.'
        with concurrent.futures.ThreadPoolExecutor() as pool:
            job = pool.submit(request, 'POST', '/api/lean/compile', {'source': slow, 'jobId': 'cancel-test'})
            deadline = time.monotonic() + 5
            while not server.jobs and time.monotonic() < deadline:
                time.sleep(.01)
            assert request('POST', '/api/lean/compile', {'source': slow, 'jobId': 'busy-test'})[0] == 409
            assert request('DELETE', '/api/lean/compile/cancel-test')[1]['cancelled']
            assert job.result(timeout=5)[1]['error'] == 'Compilation stopped.'
        code, response = request('POST', '/api/lean/compile', {'source': 'solve for x : ℝ where 7*x=21', 'jobId': 'after-cancel'})
        assert code == 200 and response['ok']
        assert response['results'][0]['solution'] == 'x = 3'
    finally:
        server.shutdown()
        server.server_close()
    print('Native Lean: elaboration, extensible proofs, executable parameters and guards, conditional solution sets, reciprocal solving, axiom rejection, C++ evaluation, numerical fallback, CLI, and API cancellation passed.')


if __name__ == '__main__':
    main()
