"""Run compiler jobs in fresh, resource-limited Docker containers.

Only this control process has the Docker socket. Workers receive source on stdin,
never a host bind mount, socket, password, or caller-controlled Docker argument.
"""
import json
import os
import re
import subprocess
import tempfile
import threading
import time
import uuid

MAX_RESULT_BYTES = 8 * 1024 * 1024


class DockerCompiler:
    def __init__(self, image, *, deployment='sample-size-lean', timeout=60,
                 memory_mb=2048, cpus=2, runtime=None):
        if not re.fullmatch(r'[a-zA-Z0-9_.-]{1,64}', deployment):
            raise ValueError('Invalid deployment name.')
        if not 1 <= timeout <= 300 or not 512 <= memory_mb <= 16384 or not 1 <= cpus <= 16:
            raise ValueError('Invalid compiler resource limits.')
        self.image, self.deployment = image, deployment
        self.timeout, self.memory_mb, self.cpus = timeout, memory_mb, cpus
        self.runtime = runtime
        self.image_id = None
        self.health_lock = threading.Lock()
        self.last_health = 0
        self.ready = False

    def available(self):
        with self.health_lock:
            if time.monotonic() - self.last_health < 5:
                return self.ready
            try:
                result = subprocess.run(['docker', 'image', 'inspect', self.image],
                                        capture_output=True, text=True, timeout=5, check=True)
                item = json.loads(result.stdout)[0]
                self.image_id = item['Id']
                self.ready = (item['Config'].get('Labels') or {}).get('ssc.lean.worker') == '1'
            except (OSError, subprocess.SubprocessError, ValueError, KeyError, IndexError, TypeError):
                self.ready = False
            self.last_health = time.monotonic()
            return self.ready

    def remove(self, name):
        try:
            subprocess.run(['docker', 'rm', '--force', name], capture_output=True, timeout=10)
        except (OSError, subprocess.SubprocessError):
            pass

    def cleanup_previous_jobs(self):
        # One controller per deployment. Recover workers left by a killed gateway.
        result = subprocess.run(['docker', 'ps', '--all', '--quiet', '--filter',
                                 f'label=ssc.lean.deployment={self.deployment}'],
                                capture_output=True, text=True, timeout=10, check=True)
        for container in result.stdout.split():
            self.remove(container)

    def create_command(self, name):
        command = ['docker', 'create', '--interactive', '--init', '--pull=never',
                   '--name', name, '--label', f'ssc.lean.deployment={self.deployment}',
                   '--network=none', '--read-only', '--user=65532:65532',
                   '--cap-drop=ALL', '--security-opt=no-new-privileges',
                   '--memory', f'{self.memory_mb}m', '--memory-swap', f'{self.memory_mb}m',
                   '--cpus', str(self.cpus), '--pids-limit=128', '--shm-size=16m',
                   '--ulimit=nofile=256:256', '--ulimit=fsize=67108864:67108864',
                   '--log-driver=none', '--tmpfs=/tmp:rw,exec,nosuid,nodev,size=256m,mode=1777',
                   '--env', f'SSC_JOB_TIMEOUT={self.timeout}',
                   '--env', f'SSC_JOB_THREADS={self.cpus}']
        if self.runtime:
            command += ['--runtime', self.runtime]
        return command + [self.image_id or self.image]

    def compile(self, source, cancel):
        if cancel.is_set():
            return {'ok': False, 'error': 'Compilation stopped.'}
        if not self.available():
            return {'ok': False, 'error': 'Compiler worker image is unavailable.'}
        name = f'ssc-lean-{uuid.uuid4().hex}'
        process = None
        deadline = time.monotonic() + self.timeout + 15
        try:
            subprocess.run(self.create_command(name), capture_output=True, timeout=10, check=True)
            if cancel.is_set():
                return {'ok': False, 'error': 'Compilation stopped.'}
            with tempfile.TemporaryFile() as source_file, tempfile.TemporaryFile() as output, tempfile.TemporaryFile() as errors:
                source_file.write(json.dumps({'source': source}, ensure_ascii=False).encode())
                source_file.seek(0)
                process = subprocess.Popen(['docker', 'start', '--attach', '--interactive', name],
                                           stdin=source_file, stdout=output, stderr=errors)
                # A file-backed stdin avoids blocking the supervising thread on a
                # full pipe if Docker or the worker stalls before reading source.
                while process.poll() is None:
                    if cancel.is_set():
                        raise RuntimeError('Compilation stopped.')
                    if time.monotonic() > deadline:
                        raise RuntimeError('Hosted compilation exceeded its time limit.')
                    if os.fstat(output.fileno()).st_size + os.fstat(errors.fileno()).st_size > MAX_RESULT_BYTES:
                        raise RuntimeError('Compiler response exceeded its output limit.')
                    time.sleep(0.05)
                if cancel.is_set():
                    raise RuntimeError('Compilation stopped.')
                if os.fstat(output.fileno()).st_size + os.fstat(errors.fileno()).st_size > MAX_RESULT_BYTES:
                    raise RuntimeError('Compiler response exceeded its output limit.')
                if process.returncode:
                    raise RuntimeError('Compiler worker failed or exceeded its resource limits.')
                output.seek(0)
                result = json.load(output)
                if not isinstance(result, dict) or not isinstance(result.get('ok'), bool):
                    raise ValueError('Invalid compiler response.')
                return result
        except (RuntimeError, OSError, subprocess.SubprocessError, ValueError) as error:
            # Do not return Docker daemon paths or configuration to clients.
            message = str(error) if isinstance(error, RuntimeError) else 'Compiler worker could not complete this request.'
            return {'ok': False, 'error': message}
        finally:
            self.remove(name)
            if process is not None and process.poll() is None:
                process.kill()
                process.wait(timeout=5)
