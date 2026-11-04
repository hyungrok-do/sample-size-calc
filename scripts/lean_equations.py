#!/usr/bin/env python3
"""Compile real Lean modules containing equation-solving commands.

The Lean parser, elaborator, metaprogram and kernel run in the pinned Lake
project. This adapter supplies the default import and transports its results.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import time

from lean_cpp import emit_cpp
from lean_study_cpp import emit_study_cpp

ROOT = Path(__file__).resolve().parents[1]
MARKER = 'SSC_LEAN_REPORT '
MAX_SOURCE_BYTES = 131072
MAX_OUTPUT_BYTES = 4 * 1024 * 1024


def compiler_available():
    return bool(shutil.which('lake') and shutil.which('c++') and
                (ROOT/'lean/.lake/build/lib/lean/DesignCompiler/Planning/Command.olean').exists())


def run_process(command, *, cwd, deadline, cancel=None):
    if cancel is not None and cancel.is_set():
        raise RuntimeError('Compilation stopped.')
    with tempfile.TemporaryFile() as output:
        process = subprocess.Popen(command, cwd=cwd, stdout=output,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        failure = None
        try:
            while process.poll() is None:
                if cancel is not None and cancel.is_set():
                    failure = 'Compilation stopped.'
                elif time.monotonic() > deadline:
                    failure = 'Native compilation exceeded its time limit.'
                elif os.fstat(output.fileno()).st_size > MAX_OUTPUT_BYTES:
                    failure = 'Compiler output exceeded the 4 MB limit.'
                if failure:
                    break
                time.sleep(0.05)
        finally:
            if process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            process.wait()
        if failure:
            raise RuntimeError(failure)
        if os.fstat(output.fileno()).st_size > MAX_OUTPUT_BYTES:
            raise RuntimeError('Compiler output exceeded the 4 MB limit.')
        if cancel is not None and cancel.is_set():
            raise RuntimeError('Compilation stopped.')
        output.seek(0)
        return process.returncode, output.read().decode('utf-8', errors='replace')


def compile_source(source, *, timeout=60, olean=None, cancel=None):
    if not isinstance(source, str) or not 0 < len(source.encode('utf-8')) <= MAX_SOURCE_BYTES:
        return {'ok': False, 'error': 'Supply 1–131072 bytes of Lean source.'}
    if not compiler_available():
        return {'ok': False, 'error': 'Build the native compiler first: cd lean && lake build DesignCompiler.Planning.Command'}
    module = 'import DesignCompiler.Planning.Command\n' + source + '\n'
    with tempfile.TemporaryDirectory(prefix='ssc-lean-') as directory:
        path = Path(directory)/((Path(olean).stem if olean else 'Source') + '.lean')
        path.write_text(module + '\n#export_compiler\n', encoding='utf-8')
        command = ['lake', 'env', 'lean', '-DmaxHeartbeats=2000000', '-DmaxRecDepth=2048', '-R', directory, str(path)]
        if olean:
            command += ['-o', str(Path(olean).resolve())]
        deadline = time.monotonic() + timeout
        try:
            returncode, stdout = run_process(command, cwd=ROOT/'lean', deadline=deadline, cancel=cancel)
        except (RuntimeError, OSError) as error:
            return {'ok': False, 'error': str(error)}
        reports, diagnostics = [], []
        for line in stdout.splitlines():
            if line.startswith(MARKER):
                reports.append(line[len(MARKER):])
            else:
                diagnostics.append(line.replace(str(path), 'source.lean'))
        if returncode or len(reports) != 1 or 'sorryAx' in stdout:
            return {'ok': False, 'error': '\n'.join(diagnostics).strip() or 'Lean did not produce a checked result.'}
        try:
            equations = json.loads(reports[0])
        except json.JSONDecodeError:
            return {'ok': False, 'error': 'Invalid native compiler response.'}
        if not equations:
            return {'ok': False, 'error': 'Add a solve command to the Lean module.'}
        try:
            for index, equation in enumerate(equations):
                if equation.get('kind') == 'statistical_design':
                    equation['cpp'] = emit_study_cpp(equation)
                    cpp, executable = Path(directory)/f'study{index}.cpp', Path(directory)/f'study{index}'
                    cpp.write_text(equation['cpp'], encoding='utf-8')
                    code, output = run_process(['c++', '-std=c++17', '-O2', str(cpp), '-o', str(executable)],
                                               cwd=ROOT, deadline=deadline, cancel=cancel)
                    if code:
                        raise RuntimeError('Statistical C++ compilation failed: ' + output)
                    code, output = run_process([str(executable)], cwd=ROOT, deadline=deadline, cancel=cancel)
                    if code:
                        raise RuntimeError('Statistical evaluation failed: ' + output)
                    equation['evaluation'] = json.loads(output)
                else:
                    equation['cpp'] = emit_cpp(equation['lowering'])
        except (RuntimeError, OSError, ValueError, OverflowError) as error:
            return {'ok': False, 'error': str(error)}
        return {'ok': True, 'compiler': 'lean-equation-compiler-v1', 'language': 'Lean 4',
                'leanKernelChecked': True, 'results': equations,
                'module': module, 'diagnostics': '\n'.join(diagnostics).strip()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--output', type=Path, help='Save the compiled derivations as JSON')
    parser.add_argument('--olean', type=Path, help='Save the compiled Lean module')
    parser.add_argument('--cpp', type=Path, help='Export C++ for a module with exactly one solve command')
    parser.add_argument('--run', action='store_true', help='Compile and execute the exported C++')
    parser.add_argument('--lower', type=float, help='Lower bound for numerical fallback')
    parser.add_argument('--upper', type=float, help='Upper bound for numerical fallback')
    parser.add_argument('--param', action='append', default=[], metavar='NAME=VALUE',
                        help='Real parameter value for --run; repeat for each parameter')
    args = parser.parse_args()
    if args.param and not args.run:
        parser.error('--param requires --run; exported C++ accepts NAME=VALUE at runtime.')
    result = compile_source(args.source.read_text(encoding='utf-8'), olean=args.olean)
    if args.cpp or args.run:
        if not result['ok']:
            parser.exit(1, result['error'] + '\n')
        if len(result['results']) != 1:
            parser.error('C++ export requires exactly one solve command per input file.')
        equation = result['results'][0]
        if equation.get('kind') == 'statistical_design' and (args.param or args.lower is not None or args.upper is not None):
            parser.error('Statistical planning values are specified in the Lean study; edit the source and recompile.')
        if not equation['cpp']:
            parser.error(equation['lowering']['reason'])
        if args.cpp:
            args.cpp.write_text(equation['cpp'], encoding='utf-8')
        if args.run:
            with tempfile.TemporaryDirectory(prefix='ssc-cpp-') as directory:
                cpp, executable = Path(directory)/'equation.cpp', Path(directory)/'equation'
                cpp.write_text(equation['cpp'], encoding='utf-8')
                subprocess.run(['c++', '-std=c++17', '-O2', str(cpp), '-o', str(executable)], check=True, timeout=30)
                command = [str(executable)]
                if equation.get('lowering', {}).get('kind') == 'bracketed_bisection':
                    if args.lower is None or args.upper is None:
                        parser.error('Numerical fallback needs --lower and --upper.')
                    command += [str(args.lower), str(args.upper)]
                command += args.param
                evaluation = subprocess.run(command, capture_output=True, text=True, timeout=10)
                if evaluation.returncode:
                    parser.exit(2, evaluation.stderr.strip() + '\n')
                result['numericalEvaluation'] = json.loads(evaluation.stdout)
    encoded = json.dumps(result, indent=2, ensure_ascii=False)
    if args.output:
        args.output.write_text(encoded + '\n', encoding='utf-8')
    else:
        print(encoded)
    raise SystemExit(0 if result['ok'] else 1)


if __name__ == '__main__':
    main()
