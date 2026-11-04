"""One request per disposable container. No server credentials enter this process."""
import json
import os
import sys

sys.path.insert(0, '/opt/compiler/scripts')
from lean_equations import MAX_SOURCE_BYTES, compile_source


def main():
    try:
        raw = sys.stdin.buffer.read(MAX_SOURCE_BYTES * 6 + 1025)
        if len(raw) > MAX_SOURCE_BYTES * 6 + 1024:
            raise ValueError('Request is too large.')
        request = json.loads(raw)
        result = compile_source(request['source'],
                                timeout=int(os.environ.get('SSC_JOB_TIMEOUT', '60')),
                                threads=int(os.environ.get('SSC_JOB_THREADS', '2')))
    except (ValueError, KeyError, TypeError) as error:
        result = {'ok': False, 'error': str(error)}
    print(json.dumps(result, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
