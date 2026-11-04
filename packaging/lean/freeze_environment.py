"""Run via `lake env python` at build time, never in a compilation request."""
import json
import os
from pathlib import Path
import sys

if 'LEAN_PATH' not in os.environ:
    raise SystemExit('Run this script through lake env after building the library.')
variables = ['LEAN_PATH', 'LEAN_SRC_PATH', 'LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH']
environment = {name: os.pathsep.join(str(Path(part or '.').resolve())
                                   for part in os.environ[name].split(os.pathsep))
               for name in variables if name in os.environ}
Path(sys.argv[1]).write_text(json.dumps(environment), encoding='utf-8')
