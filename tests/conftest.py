"""Test setup: use the source tree when the package is not installed."""
from __future__ import annotations

import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# An installed package (pip install .) carries the compiled core; in a source
# tree, ssc._native finds the library in the build*/ directories.
if importlib.util.find_spec("ssc") is None:
    sys.path.insert(0, os.path.join(ROOT, "python"))
