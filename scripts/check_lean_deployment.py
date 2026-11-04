#!/usr/bin/env python3
"""Gateway contract tests, plus optional native-worker or real Docker execution."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
from urllib.error import HTTPError
from urllib.request import Request, urlopen
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'packaging/lean'))
from app import Application, MAX_REQUEST_BYTES
from containers import DockerCompiler
from lean_equations import compile_source


def check_proof_diagnostics():
    # Lean exits successfully after sorry, so exit status alone is insufficient.
    for warning in ["warning: declaration uses 'sorry'", 'warning: declaration uses `sorry`', 'sorryAx']:
        output = warning + '\nSSC_LEAN_REPORT [{}]\n'
        with patch('lean_equations.compiler_available', return_value=True), \
                patch('lean_equations.run_process', return_value=(0, output)):
            result = compile_source('theorem invalid : False := by sorry')
        assert not result['ok'], result
    print('Legacy adapter rejects unfinished-proof diagnostics even with exit status zero.')


def request(app, path, *, method='GET', value=None, headers=None):
    body = json.dumps(value).encode() if value is not None else b''
    env = {'REQUEST_METHOD': method, 'PATH_INFO': path, 'HTTP_HOST': 'localhost:8080',
           'CONTENT_TYPE': 'application/json', 'CONTENT_LENGTH': str(len(body)),
           'wsgi.input': io.BytesIO(body)}
    env.update(headers or {})
    response = {}
    def start(status, headers):
        response.update(status=int(status.split()[0]), headers=dict(headers))
    output = b''.join(app(env, start))
    response['body'] = json.loads(output) if output and response['headers']['Content-Type'].startswith('application/json') else output
    return response


class ControlledCompiler:
    ready = True
    def __init__(self):
        self.started = threading.Event()
    def available(self):
        return self.ready
    def compile(self, source, cancel):
        self.started.set()
        assert cancel.wait(timeout=10), 'cancellation was not delivered'
        return {'ok': False, 'error': 'Compilation stopped.'}


def check_protocol():
    compiler = ControlledCompiler()
    app = Application(compiler, ROOT/'web', 'http://localhost:8080', max_jobs=1)
    r = request(app, '/api/lean/health')
    assert r['status'] == 200 and r['body']['execution'] == 'hosted'
    for headers in [{'HTTP_HOST': 'attacker.invalid'}, {'HTTP_ORIGIN': 'https://attacker.invalid'},
                    {'HTTP_SEC_FETCH_SITE': 'cross-site'}, {'HTTP_SEC_FETCH_SITE': 'same-site'}]:
        assert request(app, '/api/lean/health', headers=headers)['status'] == 403
    for page in ['/', '/dag.html', '/lean.html']:
        r = request(app, page)
        assert r['status'] == 200 and "frame-ancestors 'none'" in r['headers']['Content-Security-Policy']
        assert r['headers']['Cross-Origin-Embedder-Policy'] == 'require-corp'
        assert request(app, page, method='HEAD')['body'] == b''
    for path in ['/../README.md', '/.env', '/engine/../../README.md', '/api/nope']:
        assert request(app, path)['status'] == 404
    assert request(app, '/api/lean/compile', method='POST', value={},
                   headers={'CONTENT_LENGTH': str(MAX_REQUEST_BYTES+1)})['status'] == 413
    assert request(app, '/api/lean/compile', method='POST', value={},
                   headers={'CONTENT_TYPE': 'text/plain'})['status'] == 415
    for value in [{}, [], {'source': '', 'jobId': 'x'}, {'source': 'x', 'jobId': '../bad'},
                  {'source': 'x' * 131073, 'jobId': 'x'}]:
        assert request(app, '/api/lean/compile', method='POST', value=value)['status'] == 400
    with ThreadPoolExecutor() as pool:
        future = pool.submit(request, app, '/api/lean/compile', method='POST', value={'source': 'x', 'jobId': 'one'})
        assert compiler.started.wait(2)
        for job in ['one', 'two']:
            r = request(app, '/api/lean/compile', method='POST', value={'source': 'x', 'jobId': job})
            assert r['status'] == 429 and r['headers']['Retry-After'] == '2'
        assert request(app, '/api/lean/compile/one', method='DELETE')['status'] == 200
        assert future.result(timeout=3)['body']['error'] == 'Compilation stopped.'
    request(app, '/api/lean/compile/early', method='DELETE')
    assert not request(app, '/api/lean/compile', method='POST', value={'source': 'x', 'jobId': 'early'})['body']['ok']
    assert not app.jobs
    compiler.ready = False
    assert request(app, '/api/lean/health')['status'] == 503
    try:
        Application(compiler, ROOT/'web', 'http://remote.example.org')
        raise AssertionError('remote plaintext configuration was accepted')
    except ValueError:
        pass
    print('Gateway: origin/host, input limits, static assets, CSP, capacity, cancellation and readiness passed.', flush=True)


def check_native():
    env = dict(os.environ, PYTHONPATH=str(ROOT/'scripts'), SSC_JOB_TIMEOUT='60')
    source = (ROOT/'examples/lean-studies/repeated-measures.lean').read_text()
    with tempfile.TemporaryDirectory() as directory:
        frozen = str(Path(directory)/'lean-env.json')
        subprocess.run(['lake', 'env', sys.executable, str(ROOT/'packaging/lean/freeze_environment.py'), frozen],
                       cwd=ROOT/'lean', check=True, env=env)
        env['SSC_LEAN_ENV_FILE'] = frozen
        # A runtime regression that calls either package tool must fail this test.
        for name in ['lake', 'git']:
            forbidden = Path(directory)/name
            forbidden.write_text('#!/bin/sh\necho "Package initialization at runtime is forbidden" >&2\nexit 97\n')
            forbidden.chmod(0o755)
        env['PATH'] = directory + os.pathsep + env['PATH']
        worker = subprocess.run([sys.executable, str(ROOT/'packaging/lean/worker.py')],
                                input=json.dumps({'source': source}), capture_output=True, text=True,
                                env=env, check=True, timeout=90)
    result = json.loads(worker.stdout)
    assert result['ok'], result
    assert result['results'][0]['evaluation']['selectedN'] == 51
    assert result['results'][0]['evaluation']['totalParticipants'] == 102
    print('Worker protocol: actual Lean kernel and generated C++ returned 51 per arm.', flush=True)


def check_docker():
    compiler = DockerCompiler(os.environ.get('SSC_LEAN_WORKER_IMAGE', 'ssc-lean-worker:local'),
                              deployment='ssc-deployment-test')
    compiler.cleanup_previous_jobs()
    assert compiler.available(), 'Build the worker image first.'
    # Probe the actual worker namespace, rather than merely checking flag strings.
    probe = '''import os, socket
assert os.getuid() == 65532
assert not os.path.exists('/var/run/docker.sock')
try:
    open('/opt/compiler/write-test', 'w')
    raise AssertionError('root filesystem is writable')
except OSError:
    pass
s = socket.socket(); s.settimeout(1)
try:
    s.connect(('1.1.1.1', 443))
    raise AssertionError('worker has external networking')
except OSError:
    pass
status = open('/proc/self/status').read()
assert 'NoNewPrivs:\\t1' in status
assert 'CapEff:\\t0000000000000000' in status
print('isolation passed')
'''
    name = 'ssc-lean-isolation-test'
    try:
        create = compiler.create_command(name)
        subprocess.run(create[:-1] + ['--entrypoint=python', create[-1], '-c', probe], check=True, capture_output=True)
        config = json.loads(subprocess.check_output(['docker', 'inspect', name]))[0]['HostConfig']
        assert config['ReadonlyRootfs'] and config['NetworkMode'] == 'none'
        assert not config.get('Binds')
        assert config['Memory'] == 2048*1024*1024 and config['PidsLimit'] == 128
        assert config['NanoCpus'] == 2_000_000_000
        r = subprocess.run(['docker', 'start', '--attach', name], capture_output=True, text=True, timeout=15)
        assert r.returncode == 0 and 'isolation passed' in r.stdout, r.stderr
    finally:
        compiler.remove(name)
    source = (ROOT/'examples/lean-studies/repeated-measures.lean').read_text()
    first = compiler.compile('#eval IO.FS.writeFile "/tmp/job-private" "one"\n' + source, threading.Event())
    assert first['ok'], first
    second = compiler.compile('#eval (show IO Unit from do\n  if (← System.FilePath.pathExists "/tmp/job-private") then throw (IO.userError "leaked job file"))\n' + source, threading.Event())
    assert second['ok'], second
    assert second['results'][0]['evaluation']['selectedN'] == 51
    cancel = threading.Event()
    with ThreadPoolExecutor() as pool:
        pending = pool.submit(compiler.compile, '#eval IO.sleep 60000\n' + source, cancel)
        time.sleep(2)
        cancel.set()
        assert 'stopped' in pending.result(timeout=20)['error']
    short = DockerCompiler(compiler.image, deployment='ssc-deployment-test', timeout=1)
    assert 'time limit' in short.compile('#eval IO.sleep 60000\n' + source, threading.Event())['error']
    remaining = subprocess.check_output(['docker', 'ps', '--all', '--quiet', '--filter',
                                         'label=ssc.lean.deployment=ssc-deployment-test'], text=True)
    assert not remaining.strip(), remaining
    print('Docker: actual isolation, two fresh Lean/C++ jobs, timeout, cancellation and cleanup passed.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--native', action='store_true')
    parser.add_argument('--docker', action='store_true')
    args = parser.parse_args()
    check_proof_diagnostics()
    check_protocol()
    if args.native:
        check_native()
    if args.docker:
        check_docker()
