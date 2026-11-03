#!/usr/bin/env python3
"""Serve the editor and native Lean compiler on loopback for trusted local code.

Lean modules can run IO. This is a local development service, not a sandbox
for untrusted code and not a public compilation service.
"""
import argparse
from collections import deque
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
import re
import threading
import time
from urllib.parse import urlsplit

from lean_equations import ROOT, MAX_SOURCE_BYTES, compile_source, compiler_available


class CompilerServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, port):
        self.jobs = {}
        self.jobs_lock = threading.Condition()
        self.cancelled = deque(maxlen=256)
        super().__init__(('127.0.0.1', port), Handler)


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(ROOT/'web'), **kwargs)

    def json(self, status, value):
        body = json.dumps(value, ensure_ascii=False).encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type', 'application/json; charset=utf-8')
        self.send_header('Cache-Control', 'no-store')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        try:
            self.wfile.write(body)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def allowed(self):
        port = self.server.server_port
        hosts = {f'127.0.0.1:{port}', f'localhost:{port}'}
        host = self.headers.get('Host')
        origin = self.headers.get('Origin')
        return (host in hosts and (not origin or origin == f'http://{host}') and
                self.headers.get('Sec-Fetch-Site') not in {'cross-site', 'same-site'})

    def do_GET(self):
        if not self.allowed():
            return self.json(403, {'ok': False, 'error': 'Use this service from its loopback origin.'})
        if urlsplit(self.path).path == '/api/lean/health':
            return self.json(200, {'ok': True, 'ready': compiler_available(),
                                   'compiler': 'lean-equation-compiler-v1'})
        if urlsplit(self.path).path.startswith('/api/'):
            return self.json(404, {'ok': False, 'error': 'Unknown API route.'})
        return super().do_GET()

    def do_POST(self):
        if not self.allowed():
            return self.json(403, {'ok': False, 'error': 'Cross-origin compilation is disabled.'})
        if self.path != '/api/lean/compile':
            return self.json(404, {'ok': False, 'error': 'Unknown API route.'})
        if self.headers.get('Content-Type', '').split(';')[0] != 'application/json':
            return self.json(415, {'ok': False, 'error': 'Use application/json.'})
        try:
            length = int(self.headers.get('Content-Length', 0))
            if not 0 < length <= MAX_SOURCE_BYTES * 6 + 1024:
                raise ValueError('Request is too large or empty.')
            self.connection.settimeout(10)
            request = json.loads(self.rfile.read(length))
            source, job_id = request['source'], request['jobId']
            if not isinstance(job_id, str) or not re.fullmatch(r'[a-zA-Z0-9-]{1,80}', job_id):
                raise ValueError('Invalid job identifier.')
        except (ValueError, KeyError, TypeError, OSError) as error:
            return self.json(400, {'ok': False, 'error': str(error)})
        cancel = threading.Event()
        with self.server.jobs_lock:
            if job_id in self.server.cancelled:
                return self.json(200, {'ok': False, 'error': 'Compilation stopped.'})
            if self.server.jobs:
                return self.json(409, {'ok': False, 'error': 'Another Lean compilation is running.'})
            self.server.jobs[job_id] = cancel
        try:
            result = compile_source(source, cancel=cancel)
            self.json(200, result)
        finally:
            with self.server.jobs_lock:
                self.server.jobs.pop(job_id, None)
                self.server.jobs_lock.notify_all()

    def do_DELETE(self):
        if not self.allowed():
            return self.json(403, {'ok': False, 'error': 'Cross-origin cancellation is disabled.'})
        match = re.fullmatch(r'/api/lean/compile/([a-zA-Z0-9-]{1,80})', self.path)
        if not match:
            return self.json(404, {'ok': False, 'error': 'Unknown job.'})
        with self.server.jobs_lock:
            self.server.cancelled.append(match[1])
            event = self.server.jobs.get(match[1])
            if event:
                event.set()
                deadline = time.monotonic() + 5
                while match[1] in self.server.jobs and time.monotonic() < deadline:
                    self.server.jobs_lock.wait(timeout=0.1)
        self.json(200, {'ok': True, 'cancelled': True})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port', type=int, default=8765)
    args = parser.parse_args()
    with CompilerServer(args.port) as server:
        print(f'Lean editor: http://127.0.0.1:{server.server_port}/lean.html', flush=True)
        print('Trusted local Lean code only. Stop with Ctrl+C.', flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            for event in server.jobs.values():
                event.set()


if __name__ == '__main__':
    main()
