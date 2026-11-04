"""One request per disposable container. No server credentials enter this process."""
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, '/opt/compiler/scripts')
from lean_equations import MAX_SOURCE_BYTES, compile_source


def main():
    try:
        # Immutable build output: runtime compilation must not invoke Lake/Git,
        # refresh dependencies, or write to the installed Lean/Mathlib bundle.
        environment = json.loads(Path(os.environ.get('SSC_LEAN_ENV_FILE', '/opt/compiler/lean-env.json')).read_text())
        if not environment.get('LEAN_PATH'):
            raise ValueError('The compiler image has no frozen module path.')
        raw = sys.stdin.buffer.read(MAX_SOURCE_BYTES * 6 + 1025)
        if len(raw) > MAX_SOURCE_BYTES * 6 + 1024:
            raise ValueError('Request is too large.')
        request = json.loads(raw)
        result = compile_source(request['source'],
                                timeout=int(os.environ.get('SSC_JOB_TIMEOUT', '60')),
                                threads=int(os.environ.get('SSC_JOB_THREADS', '2')),
                                use_lake=False, compiler_env=environment)
    except (ValueError, KeyError, TypeError, OSError) as error:
        result = {'ok': False, 'error': str(error)}
    print(json.dumps(result, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
