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

ROOT = Path(__file__).resolve().parents[1]
MARKER = 'SSC_LEAN_REPORT '
MAX_SOURCE_BYTES = 131072
MAX_OUTPUT_BYTES = 4 * 1024 * 1024


def compiler_available():
    return bool(shutil.which('lake') and
                (ROOT/'lean/.lake/build/lib/lean/DesignCompiler/Solve.olean').exists())


def compile_source(source, *, timeout=60, olean=None, cancel=None):
    if not isinstance(source, str) or not 0 < len(source.encode('utf-8')) <= MAX_SOURCE_BYTES:
        return {'ok': False, 'error': 'Supply 1–131072 bytes of Lean source.'}
    if not compiler_available():
        return {'ok': False, 'error': 'Build the native compiler first: cd lean && lake build DesignCompiler.Solve'}
    module = 'import DesignCompiler.Solve\n' + source + '\n'
    with tempfile.TemporaryDirectory(prefix='ssc-lean-') as directory:
        path = Path(directory)/((Path(olean).stem if olean else 'Source') + '.lean')
        path.write_text(module + '\n#export_equations\n', encoding='utf-8')
        command = ['lake', 'env', 'lean', '-DmaxHeartbeats=2000000', '-DmaxRecDepth=2048', '-R', directory, str(path)]
        if olean:
            command += ['-o', str(Path(olean).resolve())]
        with tempfile.TemporaryFile() as output:
            process = subprocess.Popen(command, cwd=ROOT/'lean', stdout=output,
                                       stderr=subprocess.STDOUT, start_new_session=True)
            started = time.monotonic()
            failure = None
            try:
                while process.poll() is None:
                    if cancel is not None and cancel.is_set():
                        failure = 'Compilation stopped.'
                    elif time.monotonic() - started > timeout:
                        failure = f'Lean compilation exceeded {timeout} seconds.'
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
                return {'ok': False, 'error': failure}
            if os.fstat(output.fileno()).st_size > MAX_OUTPUT_BYTES:
                return {'ok': False, 'error': 'Compiler output exceeded the 4 MB limit.'}
            output.seek(0)
            stdout = output.read().decode('utf-8', errors='replace')
        reports, diagnostics = [], []
        for line in stdout.splitlines():
            if line.startswith(MARKER):
                reports.append(line[len(MARKER):])
            else:
                diagnostics.append(line.replace(str(path), 'source.lean'))
        if process.returncode or len(reports) != 1 or 'sorryAx' in stdout:
            return {'ok': False, 'error': '\n'.join(diagnostics).strip() or 'Lean did not produce a checked result.'}
        try:
            equations = json.loads(reports[0])
        except json.JSONDecodeError:
            return {'ok': False, 'error': 'Invalid native compiler response.'}
        if not equations:
            return {'ok': False, 'error': 'Add a solve command to the Lean module.'}
        for equation in equations:
            equation['cpp'] = emit_cpp(equation['lowering'])
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
    args = parser.parse_args()
    result = compile_source(args.source.read_text(encoding='utf-8'), olean=args.olean)
    if args.cpp or args.run:
        if not result['ok']:
            parser.exit(1, result['error'] + '\n')
        if len(result['results']) != 1:
            parser.error('C++ export requires exactly one solve command per input file.')
        equation = result['results'][0]
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
                if equation['lowering']['kind'] == 'bracketed_bisection':
                    if args.lower is None or args.upper is None:
                        parser.error('Numerical fallback needs --lower and --upper.')
                    command += [str(args.lower), str(args.upper)]
                result['numericalEvaluation'] = json.loads(subprocess.check_output(command, text=True, timeout=10))
    encoded = json.dumps(result, indent=2, ensure_ascii=False)
    if args.output:
        args.output.write_text(encoded + '\n', encoding='utf-8')
    else:
        print(encoded)
    raise SystemExit(0 if result['ok'] else 1)


if __name__ == '__main__':
    main()
