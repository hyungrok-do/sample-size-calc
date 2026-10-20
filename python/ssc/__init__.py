"""Sample Size Calculator: power analysis and sample size determination.

The numerical engine is a C++ library (multi-threaded); this package provides
the Python bindings (``ssc.run`` etc.), the ``ssc-calc`` command line (also ``python -m ssc``) and the
``ssc-mcp`` server for language-model clients.
"""
from __future__ import annotations

from .core import Token, catalog, evaluate, hardware_threads, run, schema, search, version

__all__ = ["Token", "catalog", "evaluate", "hardware_threads", "run", "schema", "search", "version"]
__version__ = "0.1.0"
