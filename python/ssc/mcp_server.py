"""``ssc-mcp``: the calculator as a Model Context Protocol server.

Language-model clients (Claude, ChatGPT, Codex, coding agents) get three read-only tools, the same loop a person
follows on the web page or the command line: find the method that fits a study, see its inputs, calculate.

* ``search_methods`` - methods matching a description, best first
* ``describe_method`` - a method's inputs, choices, defaults and when each input applies
* ``calculate`` - sample size, power or detectable effect, with the plain-language statement, notes, and anything
  in the request that was not used (an input that the chosen options leave out, a choice not available here)

The engine does the arithmetic; the model only chooses the method and fills in the inputs.  It runs over stdio
for a local client (``ssc-mcp``) or Streamable HTTP for a hosted one (``ssc-mcp --transport http``), where
``/healthz`` answers for health checks and ``/.well-known/openai-apps-challenge`` serves the domain-verification
token in ``SSC_OPENAI_APPS_CHALLENGE``.  Limits (``SSC_MCP_TIME_LIMIT`` seconds, ``SSC_MCP_MAX_SCENARIOS``,
``SSC_MCP_MAX_REPS``, ``SSC_MCP_THREADS``) keep a hosted server's requests in proportion.

Needs the optional ``mcp`` package: ``pip install "sample-size-calculator[mcp]"``.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
import threading
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Union

from . import core
from .inputs import unused_inputs, when_used

REPO_URL = "https://github.com/hyungrok-do/sample-size-calc"
ICON_URL = "https://raw.githubusercontent.com/hyungrok-do/sample-size-calc/main/packaging/mcp/icon.png"
MAX_ROWS = 50  # scenarios returned in full; a model rarely needs more, and each costs tokens

INSTRUCTIONS = """\
Power analysis and sample size for study designs (t-tests, proportions, survival/logrank, ANOVA, regression,
correlation, nonparametric, equivalence and non-inferiority, cluster-randomized and stepped-wedge, group-sequential,
agreement and diagnostic accuracy, confidence-interval precision, and simulation for mixed models, MMRM, win ratio,
non-proportional hazards). The engine does the arithmetic; never compute a sample size yourself.

1. search_methods with a short description of the design and outcome ("two independent means", "logrank
   non-inferiority", "cluster randomized proportions").
2. describe_method for the chosen id: its solve_for options and inputs, with each input's id, choices, default and
   when it applies.
3. calculate with method_id, solve_for (usually "n", "power", or an effect such as "mu1") and the inputs. Inputs not
   given take their defaults; read them back in the result. Several values ("0.8 0.9") or a range ("10 to 100 by
   10") give several scenarios.

Then quote the result's statement. If the result lists ignored_inputs or replaced_choices, an input you gave was not
used: change the choice that makes it apply (e.g. alloc=sep for separate group sizes) and calculate again. Read the
notes for the method's caveats. Simulation methods return a quick analytic approximation unless simulate=true.
Errors explain what to change (an effect on the wrong side of a one-sided test, a target power at or below alpha)."""


# ----------------------------------------------------------------------------- limits
@dataclass
class Limits:
    time_limit: Optional[float] = None  # seconds per calculation
    max_scenarios: int = 20000
    max_reps: Optional[int] = None      # Monte Carlo replicates per scenario
    threads: int = 0                    # 0 = every hardware thread

    @staticmethod
    def for_transport(transport: str) -> "Limits":
        """Generous on the user's own computer, in proportion on a shared server; the environment overrides."""
        lim = Limits(600, 20000, None, 0) if transport == "stdio" else Limits(60, 200, 5000, 2)

        def env(name: str, cast, current):
            text = os.environ.get(name, "").strip()
            return cast(text) if text else current

        lim.time_limit = env("SSC_MCP_TIME_LIMIT", float, lim.time_limit)
        lim.max_scenarios = env("SSC_MCP_MAX_SCENARIOS", int, lim.max_scenarios)
        lim.max_reps = env("SSC_MCP_MAX_REPS", int, lim.max_reps)
        lim.threads = env("SSC_MCP_THREADS", int, lim.threads)
        return lim


# ----------------------------------------------------------------------------- helpers (no MCP needed)
def _tidy(v: Any) -> Any:
    """Six significant digits for numbers a model reads (0.801281887689584 -> 0.801282)."""
    if isinstance(v, float):
        if not math.isfinite(v):
            return None
        return int(v) if v.is_integer() and abs(v) < 1e15 else float(f"{v:.6g}")
    return v


def _sizes(p: dict) -> str:
    """The sizes the engine allows a number that is not a count: 0 (where the input's limits allow it), or 1e-30 to
    1e30 in size."""
    lo = -math.inf if p.get("min") is None else float(p["min"])
    hi = math.inf if p.get("max") is None else float(p["max"])
    zero = (lo < 0 if p.get("min_open") else lo <= 0) and (hi > 0 if p.get("max_open") else hi >= 0)
    return ("0 or " if zero else "") + "1e-30 to 1e30 in size"


def _as_input(v: Any) -> str:
    """A model's value as engine text: numbers, '0.8 0.9', [0.8, 0.9], or [[0, 0.5, 1], [0, 0, 1]] for lists."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (list, tuple)):
        if v and all(isinstance(x, (list, tuple)) for x in v):
            return "; ".join(" ".join(_as_input(x) for x in inner) for inner in v)
        return " ".join(_as_input(x) for x in v)
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def _when(rule: list) -> str:
    return when_used(rule) if rule else "always"


def _methods() -> Dict[str, dict]:
    return {p["id"]: p for p in core.catalog()["procedures"]}


def _unknown_method(method_id: str) -> dict:
    import difflib

    close = difflib.get_close_matches(method_id, list(_methods()), n=3, cutoff=0.6)
    if not close:
        close = [h["id"] for h in core.search(method_id.replace(".", " "), 3).get("hits", [])]
    return {"ok": False, "error": f"unknown method id '{method_id}'", "did_you_mean": close or None,
            "hint": "use search_methods to find the id"}


def search_methods(description: str = "", limit: int = 5) -> dict:
    """Methods matching a description of a study, best match first (every method when empty)."""
    procs = _methods()
    if not description.strip():
        return {"methods": [{"id": p["id"], "name": p["name"], "category": p["category"]} for p in procs.values()]}
    res = core.search(description, max(1, min(int(limit), 20)))
    methods = []
    for h in res.get("hits", []):
        sch = core.schema(h["id"])
        methods.append({"id": h["id"], "name": h["name"], "category": h["category"], "summary": sch["summary"],
                        "solve_for": [{"value": o["value"], "label": o["label"]} for o in sch["solve_for"]],
                        "simulation": bool(sch.get("simulation"))})
    out: Dict[str, Any] = {"methods": methods}
    if res.get("corrected") and res["corrected"] != description:
        out["searched_for"] = res["corrected"]
    if not methods:
        out["hint"] = "no method matched; try other words, or an empty description to list every method"
    return out


def describe_method(method_id: str) -> dict:
    """A method's solve options and inputs: ids, labels, choices, defaults, ranges and when each applies."""
    if method_id not in _methods():
        return _unknown_method(method_id)
    sch = core.schema(method_id)
    inputs = []
    for p in sch["params"]:
        item: Dict[str, Any] = {"id": p["id"], "label": p["label"], "type": p["type"], "default": p["default"],
                                "used_when": _when(p.get("visible") or [])}
        if p.get("help"):
            item["help"] = p["help"]
        if p["type"] == "select":
            item["choices"] = [{"value": o["value"], "label": o["label"],
                                **({"available_when": _when(o["visible"])} if o.get("visible") else {})}
                               for o in p.get("options", [])]
        elif p["type"] in ("number", "vector"):
            rules = []
            if p.get("min") is not None:
                rules.append(f"{'>' if p.get('min_open') else '≥'} {_tidy(float(p['min']))}")
            if p.get("max") is not None:
                rules.append(f"{'<' if p.get('max_open') else '≤'} {_tidy(float(p['max']))}")
            rules.append("whole numbers" if p.get("integer") else _sizes(p))
            if p["type"] == "vector":
                rules.append("a list of numbers separated by spaces; several lists separated by ';'")
            elif p.get("multi"):
                rules.append("several values allowed")
            item["values"] = ", ".join(rules)
        inputs.append(item)
    return {"id": sch["id"], "name": sch["name"], "category": sch["category"], "summary": sch["summary"],
            "simulation": bool(sch.get("simulation")),
            "solve_for": [{"value": o["value"], "label": o["label"]} for o in sch["solve_for"]],
            "inputs": inputs, "example": sch.get("example"), "notes": sch.get("notes") or []}


def _table(columns: List[dict], rows: List[dict]) -> List[dict]:
    out = []
    for r in rows:
        item: Dict[str, Any] = {}
        for c in columns:
            if c["key"] in r and r[c["key"]] is not None:
                label = c["label"] if c["label"] not in item else f"{c['label']} ({c['key']})"
                item[label] = _tidy(r[c["key"]])
        out.append(item)
    return out


def calculate(method_id: str, solve_for: str, inputs: Optional[Dict[str, Any]] = None, simulate: bool = False,
              limits: Optional[Limits] = None) -> dict:
    """Run a calculation and return what a model needs to report it."""
    lim = limits or Limits()
    if method_id not in _methods():
        return _unknown_method(method_id)
    sch = core.schema(method_id)
    solves = [o["value"] for o in sch["solve_for"]]
    if solve_for not in solves:
        return {"ok": False, "error": f"{method_id} cannot solve for '{solve_for}'; choose one of {', '.join(solves)}"}
    given = {str(k): _as_input(v) for k, v in (inputs or {}).items()}
    params = {p["id"] for p in sch["params"]}
    unknown = sorted(set(given) - params)
    if unknown:
        return {"ok": False, "error": f"unknown input(s) for {method_id}: {', '.join(unknown)}",
                "inputs_of_this_method": sorted(params), "hint": "describe_method lists them"}
    warnings: List[str] = []
    is_sim = bool(sch.get("simulation"))
    if is_sim and simulate and lim.max_reps and "reps" in given:
        # every number in the text, however it is written ("20000," "2e4"): one too many is replaced by the cap
        found = [float(x) for x in re.findall(r"\d+(?:\.\d*)?(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?",
                                               given["reps"])]
        if any(x > lim.max_reps for x in found):
            given["reps"] = str(lim.max_reps)
            warnings.append(f"replicates capped at {lim.max_reps} on this server")
    token = core.Token()
    timer = threading.Timer(lim.time_limit, token.cancel) if lim.time_limit else None
    try:
        if timer:
            timer.start()
        res = core.run(method_id, solve_for, given, threads=lim.threads, power_curve=False,
                       approximate=is_sim and not simulate, token=token, max_scenarios=lim.max_scenarios)
    finally:
        if timer:
            timer.cancel()
            timer.join()  # a cancel already under way finishes before the token goes
        token.close()
    if res.get("cancelled"):
        return {"ok": False, "error": f"the calculation took longer than {lim.time_limit:g} seconds and was stopped",
                "hint": "fewer scenarios, fewer replicates, or the quick approximation (simulate=false)"}
    warnings = unused_inputs(res, sch) + warnings
    out: Dict[str, Any] = {"ok": bool(res.get("ok")), "method": sch["name"], "method_id": method_id,
                           "solve_for": solve_for}
    if not res.get("ok"):
        out["error"] = res.get("error", "the calculation failed")
        if res.get("input_errors"):
            out["input_errors"] = res["input_errors"]
        if warnings:
            out["warnings"] = warnings
        return out
    rows = res.get("rows", [])
    if is_sim:
        out["simulated"] = bool(simulate)
        if not simulate:
            warnings.append("quick analytic approximation; calculate again with simulate=true for the Monte Carlo "
                            "result (slower)")
    scenarios = []
    table = _table(res.get("columns", []), rows[:MAX_ROWS])
    for r, values in zip(rows[:MAX_ROWS], table):
        sc: Dict[str, Any] = {"values": values}
        if r.get("_error"):
            sc["error"] = r["_error"]
        if r.get("_statement"):
            sc["statement"] = r["_statement"]
        if r.get("_notes"):
            sc["notes"] = r["_notes"]
        scenarios.append(sc)
    out["scenarios"] = scenarios
    if len(rows) > MAX_ROWS:
        warnings.append(f"{len(rows)} scenarios; the first {MAX_ROWS} are shown")
    drop = res.get("dropout") or {}
    if drop.get("rows"):
        out["allowing_for_dropout"] = _table(drop["columns"], drop["rows"][:MAX_ROWS])
    if res.get("notes"):
        out["method_notes"] = res["notes"]
    if warnings:
        out["warnings"] = warnings
    return out


def _json(value: dict) -> str:
    """Compact JSON text: what every client reads, at the fewest tokens."""
    import json

    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


# ----------------------------------------------------------------------------- the server
def build_server(limits: Optional[Limits] = None):
    """The MCP server with its three tools (needs the ``mcp`` package)."""
    import anyio
    from mcp.server.mcpserver import MCPServer
    from mcp.types import Icon, ToolAnnotations

    lim = limits or Limits()
    server = MCPServer(name="sample-size-calculator", title="Sample Size Calculator", instructions=INSTRUCTIONS,
                       website_url=REPO_URL, version=core.version(),
                       icons=[Icon(src=ICON_URL, mimeType="image/png", sizes=["512x512"])])

    def hints(title: str) -> ToolAnnotations:
        # computes from the inputs alone: reads and changes nothing, reaches nothing outside, same answer each time
        return ToolAnnotations(title=title, readOnlyHint=True, destructiveHint=False, idempotentHint=True,
                               openWorldHint=False)

    @server.tool(name="search_methods", title="Find a sample-size method", annotations=hints("Find a sample-size method"),
                 structured_output=False)
    async def search_tool(description: str = "", limit: int = 5) -> str:
        """Find power and sample-size methods that fit a study, best match first. Describe the design and outcome
        in a few words, e.g. "two independent means", "paired proportions", "logrank survival non-inferiority",
        "cluster randomized binary outcome", "stepped wedge". An empty description lists every method."""
        return _json(await anyio.to_thread.run_sync(lambda: search_methods(description, limit)))

    @server.tool(name="describe_method", title="Show a method's inputs", annotations=hints("Show a method's inputs"),
                 structured_output=False)
    async def describe_tool(method_id: str) -> str:
        """Show a method's inputs before calculating: what it can solve for, and for each input its id, label,
        choices, default, allowed values and when it applies (some inputs only apply with certain choices, such
        as separate group sizes only with alloc=sep). Also an example and the method's caveats."""
        return _json(await anyio.to_thread.run_sync(lambda: describe_method(method_id)))

    @server.tool(name="calculate", title="Calculate sample size or power",
                 annotations=hints("Calculate sample size or power"), structured_output=False)
    async def calculate_tool(method_id: str, solve_for: str,
                             inputs: Optional[Dict[str, Union[str, float, int, bool, List[Any]]]] = None,
                             simulate: bool = False) -> str:
        """Calculate a sample size, power, or the smallest detectable effect. solve_for is one of the method's
        solve options (often "n" or "power"). inputs maps input ids to values; omitted inputs take their
        defaults, which the result shows. Several values ("0.8 0.9") or a range ("10 to 100 by 10") give one
        scenario each. Report the statement of each scenario. If warnings say an input was ignored or a choice
        replaced, fix the choice and calculate again. For simulation methods, simulate=true runs the Monte Carlo
        simulation instead of the quick approximation."""
        return _json(await anyio.to_thread.run_sync(lambda: calculate(method_id, solve_for, inputs, simulate, lim)))

    from starlette.requests import Request
    from starlette.responses import PlainTextResponse

    @server.custom_route("/healthz", methods=["GET"])
    async def health(_request: Request) -> PlainTextResponse:
        return PlainTextResponse("ok")

    @server.custom_route("/.well-known/openai-apps-challenge", methods=["GET"])
    async def openai_challenge(_request: Request) -> PlainTextResponse:
        token = os.environ.get("SSC_OPENAI_APPS_CHALLENGE", "").strip()
        return PlainTextResponse(token, status_code=200 if token else 404)

    return server


def _fewer_malloc_arenas() -> None:
    """glibc keeps freed memory in an arena per thread (up to 8 a core), and an HTTP server has many threads: under
    eight clients its memory rose from 73 MB to about 400 MB and stayed, though what it used stayed near 20 MB.
    Two arenas kept it near 100 MB.  (Only this server's own process; elsewhere the library leaves malloc alone.)"""
    if not sys.platform.startswith("linux"):
        return
    try:
        import ctypes

        ctypes.CDLL(None).mallopt(-8, 2)  # M_ARENA_MAX; a C library without mallopt (musl) is left as it is
    except (OSError, AttributeError):
        pass


def main(argv: Optional[List[str]] = None) -> int:
    _fewer_malloc_arenas()
    ap = argparse.ArgumentParser(prog="ssc-mcp", description="The sample size calculator as an MCP server.")
    ap.add_argument("--transport", choices=["stdio", "http"], default="stdio",
                    help="stdio for a local client (default); http (Streamable HTTP) to host it")
    ap.add_argument("--host", default=os.environ.get("HOST", "127.0.0.1"), help="with --transport http")
    ap.add_argument("--port", type=int, default=int(os.environ.get("PORT", "8000")), help="with --transport http")
    ap.add_argument("--path", default="/mcp", help="the MCP endpoint's path (default /mcp)")
    ap.add_argument("--allowed-host", action="append", default=[], metavar="HOST",
                    help="a Host header to accept (the public domain); repeat for several")
    args = ap.parse_args(argv)
    try:
        server = build_server(Limits.for_transport(args.transport))
    except ImportError:
        print("ssc-mcp needs the MCP package: pip install \"sample-size-calculator[mcp]\" (or pip install mcp)",
              file=sys.stderr)
        return 1
    if args.transport == "stdio":
        server.run("stdio")
        return 0
    import anyio
    from mcp.server.transport_security import TransportSecuritySettings

    allowed = args.allowed_host + [h.strip() for h in os.environ.get("SSC_MCP_ALLOWED_HOSTS", "").split(",") if h.strip()]
    if allowed:  # behind a public name: accept only it (DNS-rebinding protection)
        security = TransportSecuritySettings(enable_dns_rebinding_protection=True, allowed_hosts=allowed)
    elif args.host in ("127.0.0.1", "localhost", "::1"):
        security = None  # the SDK's protection for a local server
    else:
        security = TransportSecuritySettings(enable_dns_rebinding_protection=False)
    anyio.run(lambda: server.run_streamable_http_async(host=args.host, port=args.port, streamable_http_path=args.path,
                                                       stateless_http=True, json_response=True,
                                                       transport_security=security))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
