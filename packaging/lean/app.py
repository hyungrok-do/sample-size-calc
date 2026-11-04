"""Same-origin WSGI gateway. Production entrypoint: gunicorn 'app:create_app()'."""
from collections import deque
from http import HTTPStatus
import json
import logging
import mimetypes
import os
from pathlib import Path
import re
import threading
from urllib.parse import urlsplit

from containers import DockerCompiler

MAX_SOURCE_BYTES = 131072
MAX_REQUEST_BYTES = MAX_SOURCE_BYTES * 6 + 1024
COMPILER = 'lean-equation-compiler-v1'


class Application:
    def __init__(self, compiler, web_root, origin, *, max_jobs=2):
        parsed = urlsplit(origin)
        if (parsed.scheme not in {'https', 'http'} or not parsed.netloc or
                parsed.path or parsed.query or parsed.fragment or parsed.username or parsed.password):
            raise ValueError('SSC_LEAN_PUBLIC_ORIGIN must be a single origin without a trailing slash.')
        if parsed.scheme == 'http' and parsed.hostname not in {'localhost', '127.0.0.1', '::1'}:
            raise ValueError('Remote deployments require an HTTPS origin.')
        if not 1 <= max_jobs <= 8:
            raise ValueError('Use 1–8 concurrent compiler jobs.')
        self.compiler, self.web_root, self.origin = compiler, Path(web_root).resolve(), origin
        self.host, self.max_jobs = parsed.netloc.lower(), max_jobs
        self.lock, self.jobs, self.cancelled = threading.Lock(), {}, deque(maxlen=256)

    def response(self, start, status, value, *, headers=(), head=False, content_type='application/json; charset=utf-8'):
        body = value if isinstance(value, bytes) else json.dumps(value, ensure_ascii=False).encode()
        start(f'{status} {HTTPStatus(status).phrase}', [
            ('Content-Type', content_type), ('Content-Length', str(len(body))),
            ('Cache-Control', 'no-store'), ('X-Content-Type-Options', 'nosniff'),
            ('Cross-Origin-Opener-Policy', 'same-origin'),
            ('Cross-Origin-Embedder-Policy', 'require-corp'),
            ('Cross-Origin-Resource-Policy', 'same-origin'),
            ('Referrer-Policy', 'same-origin'), *headers])
        return [] if head else [body]

    def __call__(self, environ, start):
        method, path = environ['REQUEST_METHOD'], environ.get('PATH_INFO', '/')
        def respond(status, value, **kwargs):
            return self.response(start, status, value, **kwargs)
        if environ.get('HTTP_HOST', '').lower() != self.host:
            return respond(403, {'ok': False, 'error': 'Unexpected request host.'})
        if path.startswith('/api/'):
            if (environ.get('HTTP_ORIGIN') not in {None, self.origin} or
                    environ.get('HTTP_SEC_FETCH_SITE') in {'cross-site', 'same-site'}):
                return respond(403, {'ok': False, 'error': 'Use the compiler from its own site.'})
            if path == '/api/lean/health' and method == 'GET':
                ready = self.compiler.available()
                return respond(200 if ready else 503, {
                    'ok': True, 'ready': ready, 'compiler': COMPILER,
                    'execution': 'hosted', 'maxConcurrentJobs': self.max_jobs})
            if path == '/api/lean/compile' and method == 'POST':
                if environ.get('CONTENT_TYPE', '').split(';')[0] != 'application/json':
                    return respond(415, {'ok': False, 'error': 'Use application/json.'})
                try:
                    length = int(environ.get('CONTENT_LENGTH') or 0)
                    if not 0 < length <= MAX_REQUEST_BYTES:
                        return respond(413, {'ok': False, 'error': 'Request is too large or empty.'})
                    request = json.loads(environ['wsgi.input'].read(length))
                    source, job_id = request['source'], request['jobId']
                    if not isinstance(source, str) or not 0 < len(source.encode()) <= MAX_SOURCE_BYTES:
                        raise ValueError('Supply 1–131072 bytes of Lean source.')
                    if not isinstance(job_id, str) or not re.fullmatch(r'[a-zA-Z0-9-]{1,80}', job_id):
                        raise ValueError('Invalid job identifier.')
                except (ValueError, KeyError, TypeError):
                    return respond(400, {'ok': False, 'error': 'Invalid compiler request.'})
                cancel = threading.Event()
                with self.lock:
                    if job_id in self.cancelled:
                        return respond(200, {'ok': False, 'error': 'Compilation stopped.'})
                    if job_id in self.jobs or len(self.jobs) >= self.max_jobs:
                        return respond(429, {'ok': False, 'error': 'Compiler is busy. Try again shortly.'},
                                       headers=[('Retry-After', '2')])
                    self.jobs[job_id] = cancel
                try:
                    result = self.compiler.compile(source, cancel)
                    return respond(200 if result.get('ok') or cancel.is_set() else 422, result)
                except Exception:
                    logging.exception('Compiler adapter failed')
                    return respond(503, {'ok': False, 'error': 'Compiler temporarily unavailable.'})
                finally:
                    with self.lock:
                        self.jobs.pop(job_id, None)
            match = re.fullmatch(r'/api/lean/compile/([a-zA-Z0-9-]{1,80})', path)
            if match and method == 'DELETE':
                with self.lock:
                    self.cancelled.append(match[1])
                    if event := self.jobs.get(match[1]):
                        event.set()
                return respond(200, {'ok': True, 'cancelled': True})
            return respond(404, {'ok': False, 'error': 'Unknown API route.'})
        if method not in {'GET', 'HEAD'}:
            return respond(405, {'ok': False, 'error': 'Method not allowed.'})
        relative = 'index.html' if path == '/' else path.lstrip('/')
        if any(part.startswith('.') for part in Path(relative).parts):
            return respond(404, {'ok': False, 'error': 'Not found.'})
        file = (self.web_root / relative).resolve()
        if not file.is_relative_to(self.web_root) or not file.is_file():
            return respond(404, {'ok': False, 'error': 'Not found.'})
        body = file.read_bytes()
        headers = []
        if file.suffix == '.html':
            policy = re.search(rb'<meta http-equiv="Content-Security-Policy" content="([^"]+)"', body)
            if policy:
                headers.append(('Content-Security-Policy', policy[1].decode() + "; frame-ancestors 'none'"))
        content_type = {'.mjs': 'text/javascript', '.js': 'text/javascript', '.wasm': 'application/wasm'}.get(
            file.suffix, mimetypes.guess_type(file.name)[0] or 'application/octet-stream')
        return respond(200, body, headers=headers, head=method == 'HEAD', content_type=content_type)


def create_app():
    compiler = DockerCompiler(os.environ['SSC_LEAN_WORKER_IMAGE'],
                              deployment=os.environ.get('SSC_LEAN_DEPLOYMENT', 'sample-size-lean'),
                              timeout=int(os.environ.get('SSC_LEAN_JOB_TIMEOUT', '60')),
                              memory_mb=int(os.environ.get('SSC_LEAN_JOB_MEMORY_MB', '2048')),
                              cpus=int(os.environ.get('SSC_LEAN_JOB_CPUS', '2')),
                              runtime=os.environ.get('SSC_LEAN_RUNTIME') or None)
    compiler.cleanup_previous_jobs()
    return Application(compiler, os.environ.get('SSC_LEAN_WEB_ROOT', '/srv/site'),
                       os.environ['SSC_LEAN_PUBLIC_ORIGIN'],
                       max_jobs=int(os.environ.get('SSC_LEAN_MAX_JOBS', '2')))
