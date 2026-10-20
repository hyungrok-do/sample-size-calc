"""Set up, or bring up to date, the repository's virtual environment with the calculator in it.

    python3 scripts/setup_venv.py            # .venv in the repository (run it again after a git pull)
    python3 scripts/setup_venv.py --tests    # also what the tests need (pytest)
    python3 scripts/setup_venv.py --fresh    # delete .venv and start over

The environment is the repository's own (.venv), so the commands (ssc-calc,
ssc-mcp) and the scripts run from one place, whatever else is installed or
active (conda's base, the system's Python).  Installing builds the C++ engine,
so a C++17 compiler and CMake are needed (see the README).  Afterwards:

    source .venv/bin/activate          # Windows: .venv\\Scripts\\activate
    ssc-calc --help
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VENV = ROOT / ".venv"
WINDOWS = sys.platform == "win32"


def python_in(env: Path) -> Path:
    return env / ("Scripts/python.exe" if WINDOWS else "bin/python")


def run(cmd, **kwargs) -> None:
    print("  $ " + " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, **kwargs)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tests", action="store_true", help="also install what the tests need (pytest)")
    ap.add_argument("--fresh", action="store_true", help="delete the environment first and build it again")
    args = ap.parse_args()
    if sys.version_info < (3, 9):
        print(f"Python 3.9 or later is needed; this is {sys.version.split()[0]}.", file=sys.stderr)
        return 1
    if os.environ.get("CONDA_DEFAULT_ENV"):
        print(f"Note: conda's environment '{os.environ['CONDA_DEFAULT_ENV']}' is active. The .venv is made from "
              "its Python, which is fine; to use it, activate only .venv (run `conda deactivate` first, "
              "so the two do not mix).")
    if args.fresh and VENV.exists():
        if Path(sys.prefix).resolve() == VENV.resolve():
            print("Run this with a Python outside .venv to start over (deactivate first).", file=sys.stderr)
            return 1
        print(f"Deleting {VENV}")
        shutil.rmtree(VENV)
    py = python_in(VENV)
    if not py.exists():
        print(f"Creating {VENV} with Python {sys.version.split()[0]}")
        venv.EnvBuilder(with_pip=True).create(VENV)
    target = f"{ROOT}[test]" if args.tests else str(ROOT)
    try:
        run([py, "-m", "pip", "install", "--quiet", "--upgrade", "pip"])
        print("Building the engine and installing the package (a minute or two)…", flush=True)
        run([py, "-m", "pip", "install", target])
        run([py, "-c", "import ssc; from ssc import core; "
                       "print(f'Installed: ssc {ssc.__version__}, engine {core.version()}, "
                       "{len(core.catalog()[\"procedures\"])} methods')"])
    except subprocess.CalledProcessError:
        print("\nThe installation failed. Building the engine needs a C++17 compiler and CMake 3.15 or later "
              "(macOS: xcode-select --install; see the README). `--fresh` starts the environment over.",
              file=sys.stderr)
        return 1
    activate = r".venv\Scripts\activate" if WINDOWS else "source .venv/bin/activate"
    print(f"\nReady. In this folder:\n\n    {activate}\n    ssc-calc --help         # the command line\n\n"
          "After a git pull, run this script again to bring the installation up to date.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
