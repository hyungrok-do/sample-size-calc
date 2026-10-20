"""Locate and load the C++ core (libssc_core) through ctypes.

The core exposes a small C ABI, so one compiled library serves every Python 3
version.  ctypes releases the GIL during foreign calls, so long computations
run in a worker thread without freezing the user interface.

Search order for the library:
  1. the path in the environment variable ``SSC_CORE_LIB``;
  2. ``ssc/_native/`` inside the installed package (wheel / pip install);
  3. ``build*/`` directories of a source checkout (developer builds).
"""
from __future__ import annotations

import ctypes
import os
import sys
from pathlib import Path
from typing import Optional

_LIB_NAMES = {
    "win32": ["libssc_core.dll", "ssc_core.dll"],
    "darwin": ["libssc_core.dylib", "libssc_core.so"],
}
_DEFAULT_NAMES = ["libssc_core.so"]


def _candidate_names() -> list[str]:
    for prefix, names in _LIB_NAMES.items():
        if sys.platform.startswith(prefix):
            return names
    return _DEFAULT_NAMES


def _candidates() -> list[Path]:
    paths: list[Path] = []
    env = os.environ.get("SSC_CORE_LIB")
    if env:
        paths.append(Path(env))
    here = Path(__file__).resolve().parent
    for name in _candidate_names():
        paths.append(here / "_native" / name)
    # developer checkout: <repo>/python/ssc/_native.py -> <repo>/build*/
    repo = here.parent.parent
    if repo.is_dir():
        for build in sorted(repo.glob("build*")):
            for name in _candidate_names():
                paths.append(build / name)
                paths.append(build / "Release" / name)
                paths.append(build / "RelWithDebInfo" / name)
    return paths


class CoreNotFoundError(RuntimeError):
    pass


def load_library() -> ctypes.CDLL:
    tried = []
    for path in _candidates():
        if path.is_file():
            if sys.platform.startswith("win") and hasattr(os, "add_dll_directory"):
                try:
                    os.add_dll_directory(str(path.parent))
                except OSError:
                    pass
            lib = ctypes.CDLL(str(path))
            _declare(lib)
            return lib
        tried.append(str(path))
    raise CoreNotFoundError(
        "Could not find the compiled core library (libssc_core). Build it with\n"
        "    cmake -S . -B build && cmake --build build --config Release\n"
        "or install the package with  pip install .  — or set SSC_CORE_LIB.\n"
        "Searched:\n  " + "\n  ".join(tried)
    )


def _declare(lib: ctypes.CDLL) -> None:
    c_void_p, c_char_p, c_int, c_double = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_double
    lib.ssc_version.restype = c_char_p
    lib.ssc_version.argtypes = []
    lib.ssc_hardware_threads.restype = c_int
    lib.ssc_hardware_threads.argtypes = []
    lib.ssc_catalog_json.restype = c_void_p
    lib.ssc_catalog_json.argtypes = []
    lib.ssc_procedure_json.restype = c_void_p
    lib.ssc_procedure_json.argtypes = [c_char_p]
    lib.ssc_search_json.restype = c_void_p
    lib.ssc_search_json.argtypes = [c_char_p, c_int]
    lib.ssc_run_json.restype = c_void_p
    lib.ssc_run_json.argtypes = [c_char_p, c_void_p]
    lib.ssc_free.restype = None
    lib.ssc_free.argtypes = [c_void_p]
    lib.ssc_token_new.restype = c_void_p
    lib.ssc_token_new.argtypes = []
    lib.ssc_token_cancel.restype = None
    lib.ssc_token_cancel.argtypes = [c_void_p]
    lib.ssc_token_progress.restype = c_double
    lib.ssc_token_progress.argtypes = [c_void_p]
    lib.ssc_token_free.restype = None
    lib.ssc_token_free.argtypes = [c_void_p]
    lib.ssc_eval.restype = c_int
    lib.ssc_eval.argtypes = [c_char_p, ctypes.POINTER(c_double), c_int, ctypes.POINTER(c_double)]


_lib: Optional[ctypes.CDLL] = None


def lib() -> ctypes.CDLL:
    global _lib
    if _lib is None:
        _lib = load_library()
    return _lib


def take_string(ptr: Optional[int]) -> Optional[str]:
    """Copy a char* returned by the core into Python and free it."""
    if not ptr:
        return None
    try:
        return ctypes.string_at(ptr).decode("utf-8", "replace")
    finally:
        lib().ssc_free(ptr)
