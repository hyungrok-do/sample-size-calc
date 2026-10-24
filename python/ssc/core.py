"""Pythonic access to the C++ computational core.

>>> import ssc
>>> res = ssc.run("means.two.t", "n", {"power": "0.8 0.9", "mu1": "1", "mu2": "0", "sigma": "1"})
>>> [row["n1"] for row in res["rows"]]
[17, 23]

All computation happens in C++; this module only converts to and from JSON.
"""
from __future__ import annotations

import ctypes
import json
import threading
from typing import Any, Mapping, Optional, Sequence

from . import _native


class Token:
    """Cancellation / progress handle for a running calculation.

    Thread-safe: it may be cancelled, read or closed from any thread.  Closing waits for a cancel already under
    way (a time limit firing as the calculation ends), and a token closed while a calculation uses it is freed
    when that calculation returns.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._busy = 0  # calculations using it
        self._closed = False
        self._ptr = _native.lib().ssc_token_new()
        if not self._ptr:
            raise MemoryError("could not allocate token")

    @property
    def handle(self) -> int:
        return self._ptr

    def cancel(self) -> None:
        with self._lock:
            if self._ptr:
                _native.lib().ssc_token_cancel(self._ptr)

    def progress(self) -> float:
        with self._lock:
            return float(_native.lib().ssc_token_progress(self._ptr)) if self._ptr else 0.0

    def close(self) -> None:
        with self._lock:
            self._closed = True
            if self._busy or not self._ptr:
                return  # freed when the last calculation using it returns
            ptr, self._ptr = self._ptr, None
        _native.lib().ssc_token_free(ptr)

    def _use(self) -> int:
        """The handle for a calculation, which must call _done() when it returns."""
        with self._lock:
            if not self._ptr:
                raise ValueError("the token is closed")
            self._busy += 1
            return self._ptr

    def _done(self) -> None:
        with self._lock:
            self._busy -= 1
            if self._busy or not self._closed or not self._ptr:
                return
            ptr, self._ptr = self._ptr, None
        _native.lib().ssc_token_free(ptr)

    def __enter__(self) -> "Token":
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    def __del__(self) -> None:  # pragma: no cover - best effort
        try:
            self.close()
        except Exception:
            pass


def version() -> str:
    return _native.lib().ssc_version().decode("ascii")


def hardware_threads() -> int:
    return int(_native.lib().ssc_hardware_threads())


# What the engine answers when a calculation runs out of memory (core/src/capi.cpp); run() gives the same answer when
# there was not even memory for that, and the engine returned NULL.
_NO_MEMORY = ("There is not enough memory for this calculation: use fewer combinations of input values or shorter "
              "lists.")
# What catalog(), search() and schema() raise MemoryError with when the engine returned NULL for want of memory
_NO_MEMORY_FOR_ANSWER = "there is not enough memory for the engine's answer"


def _answer(ptr: Optional[int]) -> dict:
    """The JSON a call returns, which is NULL only when there was not enough memory for it."""
    text = _native.take_string(ptr)
    if text is None:
        raise MemoryError(_NO_MEMORY_FOR_ANSWER)
    return json.loads(text)


def catalog() -> dict:
    """Categories and the list of procedures."""
    return _answer(_native.lib().ssc_catalog_json())


def search(query: str, limit: int = 20) -> dict:
    """Procedures ranked for a free-text query.

    Tolerates typos, plurals, hyphenation, abbreviations and synonyms ("NI",
    "log rank", "equivalance").  Returns ``{"query", "corrected", "partial",
    "hits": [{"id", "name", "category", "score", "terms"}, ...]}``; ``corrected``
    is the query with misspellings fixed ("" if none) and ``partial`` is true
    when no procedure matches every word.  ``limit=0`` returns every match.
    """
    return _answer(_native.lib().ssc_search_json(query.encode("utf-8", "replace"), int(limit)))


def schema(procedure_id: str) -> dict:
    """Full parameter schema of one procedure."""
    raw = str(procedure_id).encode("utf-8", "replace")   # a lone surrogate (from undecodable bytes) becomes "?"
    text = _native.take_string(_native.lib().ssc_procedure_json(raw))
    if text is None:  # an unknown id, or a known one without memory for its schema
        name = raw.decode("utf-8")
        if any(p["id"] == name for p in catalog()["procedures"]):
            raise MemoryError(_NO_MEMORY_FOR_ANSWER)
        raise KeyError(f"unknown procedure: {name}")
    return json.loads(text)


def _as_text(v: Any) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (list, tuple)):
        return " ".join(_as_text(x) for x in v)
    return str(v)


def run(
    procedure_id: str,
    solve_for: str,
    inputs: Mapping[str, Any],
    *,
    threads: int = 0,
    power_curve: bool = True,
    curve_points: int = 41,
    approximate: bool = False,
    token: Optional[Token] = None,
    max_scenarios: Optional[int] = None,
) -> dict:
    """Run a calculation and return the result dictionary.

    ``inputs`` maps parameter ids to text such as ``"0.8 0.9"`` or
    ``"10 to 100 by 10"`` (numbers and lists are converted automatically).
    ``threads=0`` uses every hardware thread.  ``approximate=True`` asks a
    simulation procedure for its quick analytic approximation instead of
    running the simulation; other procedures ignore it.  ``max_scenarios``
    lowers the limit on the number of input combinations (20,000); it cannot
    raise it.
    """
    request = {
        "procedure": procedure_id,
        "solve_for": solve_for,
        "inputs": {k: _as_text(v) for k, v in inputs.items()},
        "options": {"threads": int(threads), "power_curve": bool(power_curve), "curve_points": int(curve_points),
                    "approximate": bool(approximate)},
    }
    if max_scenarios is not None:
        request["options"]["max_scenarios"] = int(max_scenarios)
    payload = json.dumps(request).encode("utf-8")
    handle = token._use() if token is not None else None
    try:
        text = _native.take_string(_native.lib().ssc_run_json(payload, handle))
    finally:
        if token is not None:
            token._done()
    if text is None:  # not even the engine's short answer for running out of memory could be allocated
        return {"ok": False, "error": _NO_MEMORY}
    return json.loads(text)


def evaluate(name: str, *args: float) -> float:
    """Evaluate a distribution function of the core, e.g. ``evaluate("pnt", 2.0, 10, 1.5, 1)``."""
    arr = (ctypes.c_double * max(len(args), 1))(*args)
    out = ctypes.c_double()
    rc = _native.lib().ssc_eval(name.encode("ascii"), arr, len(args), ctypes.byref(out))
    if rc != 0:
        raise ValueError(f"unknown function or wrong number of arguments: {name}{tuple(args)}")
    return out.value


def run_many(requests: Sequence[Mapping[str, Any]], threads: int = 0) -> list[dict]:
    """Convenience: run several requests one after another (each is parallel inside)."""
    return [run(r["procedure"], r["solve_for"], r.get("inputs", {}), threads=threads) for r in requests]
