"""``python -m ssc`` is the ``ssc-calc`` command line."""
from __future__ import annotations

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
