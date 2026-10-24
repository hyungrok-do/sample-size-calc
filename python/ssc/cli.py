"""``ssc-calc``: run sample size and power calculations from the command line.

Examples::

    ssc-calc list
    ssc-calc list --search survival
    ssc-calc search equivalance two means
    ssc-calc show means.two.t
    ssc-calc run means.two.t --solve n power="0.8 0.9" mu1=1 mu2=0 sigma=1
    ssc-calc run gs.means --example --json
    ssc-calc eval pnt 2.0 10 1.5 1

Inputs are ``id=value`` pairs using the parameter ids shown by ``show``; a
value may be a list ("0.8 0.9") or a range ("10 to 100 by 10").  All
computation happens in the C++ engine; Ctrl-C cancels a running calculation.
"""
from __future__ import annotations

import argparse
import json
import sys
import textwrap
import threading
import time
from typing import Dict, List, Optional, Sequence

from . import __version__, core
from .formatting import ReportOptions, fmt_cell, result_table_text, sub_table_text, table_csv, text_table
from .inputs import unused_inputs, when_used


def _wrap(text: str, indent: str = "", width: int = 88, hang: str = "") -> str:
    return textwrap.fill(text, width=width, initial_indent=indent, subsequent_indent=indent + hang,
                         break_on_hyphens=False)  # keep "two-sided" together


# ---------------------------------------------------------------------- list
def cmd_list(args: argparse.Namespace) -> int:
    cat = core.catalog()
    order = {c: i for i, c in enumerate(cat["categories"])}
    procs = cat["procedures"]
    if (args.search or "").strip():
        found = {h["id"] for h in core.search(args.search, 0)["hits"]}
        procs = [p for p in procs if p["id"] in found]
    if args.category:
        procs = [p for p in procs if p["category"].lower() == args.category.lower()]
    if args.json:
        print(json.dumps(procs, indent=2, ensure_ascii=False))
        return 0
    if not procs:
        print("No procedures match.", file=sys.stderr)
        return 1
    procs.sort(key=lambda p: (order.get(p["category"], 99), p.get("subcategory", ""), p["name"]))
    width = max(len(p["id"]) for p in procs)
    current = None
    for p in procs:
        if p["category"] != current:
            current = p["category"]
            print(("\n" if current != procs[0]["category"] else "") + current)
        print(f"  {p['id'].ljust(width)}  {p['name']}")
    return 0


# ---------------------------------------------------------------------- show
def cmd_show(args: argparse.Namespace) -> int:
    sch = core.schema(args.procedure)
    if args.json:
        print(json.dumps(sch, indent=2, ensure_ascii=False))
        return 0
    print(sch["name"])
    print(f"{sch['category']} › {sch['subcategory']}   (id: {sch['id']})\n")
    print(_wrap(sch["summary"]))
    print("\nSolve for:")
    for o in sch["solve_for"]:
        print(f"  {o['value']:<10} {o['label']}")
    print("\nParameters:")
    group = None
    for p in sch["params"]:
        if p["group"] != group:
            group = p["group"]
            print(f"  [{group}]")
        default = p["default"]
        line = f"    {p['id']:<12} {p['label']}"
        if p["type"] == "select":
            line += f"  (default {default})"
            print(line)
            for o in p["options"]:
                cond = f"   [when {when_used(o['visible'])}]" if o.get("visible") else ""
                print(f"        {o['value']:<12} {o['label']}{cond}")
        else:
            print(line + (f"  (default {default})" if default != "" else ""))
        if p.get("visible"):
            print(f"        shown when {when_used(p['visible'])}")
    ex = sch.get("example") or {}
    if ex:
        pairs = " ".join(f'{k}="{v}"' if " " in v else f"{k}={v}" for k, v in ex.get("inputs", {}).items())
        print(f"\nExample:\n  ssc-calc run {sch['id']} --solve {ex.get('solve_for', 'n')} {pairs}")
    return 0


# ---------------------------------------------------------------------- run
def _decimals(text: str) -> int:
    """--decimals: a whole number from 0 to 12."""
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{text}' is not a whole number") from None
    if not 0 <= value <= 12:
        raise argparse.ArgumentTypeError("give a number of decimals from 0 to 12")
    return value


def _parse_inputs(pairs: Sequence[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in pairs:
        if "=" not in item:
            raise SystemExit(f"ssc-calc: expected id=value, got {item!r}")
        key, value = item.split("=", 1)
        out[key.strip()] = value.strip()
    return out


def _run_cancellable(procedure: str, solve: str, inputs: Dict[str, str], *, threads: int, curve: bool,
                     show_progress: bool) -> Optional[dict]:
    """Run in a worker thread so that Ctrl-C can cancel the C++ computation."""
    token = core.Token()
    box: Dict[str, object] = {}

    def work() -> None:
        try:
            box["result"] = core.run(procedure, solve, inputs, threads=threads, power_curve=curve, token=token)
        except BaseException as exc:  # re-raised in the main thread
            box["error"] = exc

    worker = threading.Thread(target=work, daemon=True)
    started = time.monotonic()
    worker.start()
    shown = False
    try:
        while worker.is_alive():
            worker.join(0.1)
            if show_progress and time.monotonic() - started > 1.0:
                sys.stderr.write(f"\rcomputing… {100 * token.progress():5.1f}%")
                sys.stderr.flush()
                shown = True
    except KeyboardInterrupt:
        token.cancel()
        worker.join()
        if shown:
            sys.stderr.write("\r" + " " * 24 + "\r")
        print("Cancelled.", file=sys.stderr)
        token.close()
        return None
    if shown:
        sys.stderr.write("\r" + " " * 24 + "\r")
    token.close()
    if "error" in box:
        raise box["error"]  # type: ignore[misc]
    return box["result"]  # type: ignore[return-value]


def _print_result(res: dict, opts: ReportOptions, *, statements: str, details: bool) -> None:
    print(res.get("name", ""))
    print(f"Solve for: {res.get('solve_label', res.get('solve_for'))}   "
          f"({res.get('scenario_count', 0)} scenario(s), {res.get('threads', 1)} thread(s), "
          f"{res.get('elapsed_ms', 0):.0f} ms)\n")
    rows = res.get("rows", [])
    print(result_table_text(res.get("columns", []), rows, opts, indent="  "))
    if statements != "none":
        stm = [(i, r["_statement"]) for i, r in enumerate(rows) if r.get("_statement")]
        if statements == "first":
            stm = stm[:1]
        if stm:
            print()
            for i, text in stm:
                prefix = f"[{i + 1}] " if len(rows) > 1 else ""
                print(_wrap(prefix + text, "  ", hang=" " * len(prefix)))
    if details:
        for i, r in enumerate(rows):
            for t in r.get("_tables", []):
                title = t.get("title", "")
                print(f"\n  {title}" + (f" — row {i + 1}" if len(rows) > 1 else ""))
                print(sub_table_text(t, opts, indent="  "))
    drop = res.get("dropout")
    if drop and drop.get("rows"):
        print("\n  Sample sizes allowing for dropout")
        cols = drop["columns"]
        body = [[fmt_cell(r.get(c["key"]), c.get("fmt", "num"), opts) for c in cols] for r in drop["rows"]]
        print(text_table([c["label"] for c in cols], body, indent="  "))
    notes = list(res.get("notes", []))
    for r in rows:
        for n in r.get("_notes", []):
            if n not in notes:
                notes.append(n)
    if notes:
        print("\nNotes:")
        for n in notes:
            print(_wrap("• " + n, "  ", hang="  "))


def cmd_run(args: argparse.Namespace) -> int:
    sch = core.schema(args.procedure)
    inputs: Dict[str, str] = {}
    solve = args.solve
    if args.example:
        ex = sch.get("example") or {}
        inputs.update(ex.get("inputs", {}))
        solve = solve or ex.get("solve_for")
    inputs.update(_parse_inputs(args.inputs))
    valid = {p["id"] for p in sch["params"]}
    unknown = sorted(set(inputs) - valid)
    if unknown:
        print(f"ssc-calc: unknown parameter(s) for {args.procedure}: {', '.join(unknown)} "
              f"(see 'ssc-calc show {args.procedure}')", file=sys.stderr)
        return 2
    solve = solve or sch["solve_for"][0]["value"]
    if solve not in {o["value"] for o in sch["solve_for"]}:
        print(f"ssc-calc: cannot solve for {solve!r}; choose one of "
              f"{', '.join(o['value'] for o in sch['solve_for'])}", file=sys.stderr)
        return 2
    res = _run_cancellable(args.procedure, solve, inputs, threads=args.threads,
                           curve=args.json and not args.no_curve,
                           show_progress=not args.json and not args.csv and sys.stderr.isatty())
    if res is None:
        return 130
    for line in unused_inputs(res, sch):
        print(f"ssc-calc: {line}", file=sys.stderr)
    if args.docx and res.get("ok"):
        from .docx import write_docx

        try:
            write_docx(args.docx, res, ReportOptions(prob_decimals=args.decimals, statements=args.statements))
        except OSError as e:
            print(f"ssc-calc: cannot write {args.docx}: {e.strerror or e}", file=sys.stderr)
            return 1
        print(f"Wrote {args.docx}", file=sys.stderr)
    # every scenario refused (the target power at or below alpha, say): the run computed nothing
    rows = res.get("rows") or []
    nothing = bool(res.get("ok")) and bool(rows) and all(r.get("_error") for r in rows)
    if args.json:
        print(json.dumps(res, indent=2, ensure_ascii=False))
        return 0 if res.get("ok") and not nothing else 1
    if not res.get("ok"):
        print(f"Error: {res.get('error', 'calculation failed')}", file=sys.stderr)
        for key, msg in (res.get("input_errors") or {}).items():
            print(f"  {key}: {msg}", file=sys.stderr)
        return 1
    if nothing:
        print(f"Error: {rows[0]['_error']}", file=sys.stderr)
    if args.csv:
        sys.stdout.write(table_csv(res))
        return 1 if nothing else 0
    _print_result(res, ReportOptions(prob_decimals=args.decimals), statements=args.statements, details=not args.brief)
    return 1 if nothing else 0


# ---------------------------------------------------------------------- search
def cmd_search(args: argparse.Namespace) -> int:
    query = " ".join(args.words)
    res = core.search(query, args.limit)
    if args.json:
        print(json.dumps(res, indent=2, ensure_ascii=False))
        return 0 if res["hits"] else 1
    if not res["hits"]:
        print(f"No procedure matches \"{query}\".", file=sys.stderr)
        return 1
    if res.get("corrected"):
        print(f"Showing results for \"{res['corrected']}\"")
    elif res.get("partial"):
        print("No procedure matches every word; the closest matches:")
    width = max(len(h["id"]) for h in res["hits"])
    for h in res["hits"]:
        print(f"  {h['id'].ljust(width)}  {h['name']}")
    return 0


# ---------------------------------------------------------------------- eval / version
def cmd_eval(args: argparse.Namespace) -> int:
    try:
        value = core.evaluate(args.function, *[float(a) for a in args.args])
    except ValueError as exc:
        print(f"ssc-calc: {exc}", file=sys.stderr)
        return 2
    print(repr(value))
    return 0


def cmd_version(args: argparse.Namespace) -> int:
    print(f"ssc-calc {__version__} (engine {core.version()}, {core.hardware_threads()} hardware threads)")
    return 0


# ---------------------------------------------------------------------- main
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ssc-calc",
        description="Sample size and power calculations (command-line interface to the C++ engine).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Run 'ssc-calc show <procedure>' to see a procedure's parameters.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("list", help="list procedures")
    p.add_argument("--search", "-s", help="only procedures matching these words (as in 'search')")
    p.add_argument("--category", "-c", help="only this category")
    p.add_argument("--json", action="store_true", help="machine-readable output")
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("search", help="find procedures by keywords, best match first",
                       description="Find procedures by name, topic, keyword or abbreviation (\"NI\", \"ICC\"). "
                                   "Misspelled and partly typed words are matched too.")
    p.add_argument("words", nargs="+", help="e.g. logrank, paired t test, non-inferiority proportions")
    p.add_argument("--limit", "-n", type=int, default=10, help="number of results (default 10; 0 = all)")
    p.add_argument("--json", action="store_true", help="machine-readable output")
    p.set_defaults(func=cmd_search)

    p = sub.add_parser("show", help="describe a procedure's parameters")
    p.add_argument("procedure")
    p.add_argument("--json", action="store_true", help="print the full JSON schema")
    p.set_defaults(func=cmd_show)

    p = sub.add_parser("run", help="run a calculation", formatter_class=argparse.RawDescriptionHelpFormatter,
                       description="Run a calculation; inputs are id=value pairs.\n\nExample:\n"
                                   "  ssc-calc run means.two.t --solve n power=\"0.8 0.9\" mu1=1 mu2=0 sigma=1")
    p.add_argument("procedure")
    p.add_argument("inputs", nargs="*", metavar="id=value")
    p.add_argument("--solve", help="quantity to solve for (default: the procedure's first option)")
    p.add_argument("--example", action="store_true", help="start from the procedure's example inputs")
    p.add_argument("--threads", type=int, default=0, help="worker threads (0 = all cores)")
    p.add_argument("--decimals", type=_decimals, default=4, help="decimals for probabilities, 0 to 12 (default 4)")
    p.add_argument("--statements", choices=["first", "all", "none"], default="all",
                   help="summary sentences to print (default all)")
    p.add_argument("--brief", action="store_true", help="omit per-scenario detail tables")
    fmt = p.add_mutually_exclusive_group()
    fmt.add_argument("--json", action="store_true", help="print the full JSON result")
    fmt.add_argument("--csv", action="store_true", help="print the result table as CSV")
    p.add_argument("--no-curve", action="store_true", help="with --json: skip the power-curve data")
    p.add_argument("--docx", metavar="FILE", help="also write the report as a Word document (.docx)")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("eval", help="evaluate a distribution function, e.g. 'pnt 2 10 1.5 1'")
    p.add_argument("function")
    p.add_argument("args", nargs="*")
    p.set_defaults(func=cmd_eval)

    p = sub.add_parser("version", help="print version information")
    p.set_defaults(func=cmd_version)

    return parser


def _unicode_streams() -> None:
    """Make sure output can carry Greek letters and table rules: a redirected
    stream on Windows defaults to a legacy code page."""
    for stream in (sys.stdout, sys.stderr):
        enc = (getattr(stream, "encoding", None) or "").lower().replace("-", "").replace("_", "")
        if enc not in ("utf8", "utf16", "utf32") and hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8", errors="replace")
            except (ValueError, OSError):
                pass


def main(argv: Optional[Sequence[str]] = None) -> int:
    _unicode_streams()
    parser = build_parser()
    # id=value inputs may follow options (--solve n power=0.8): argparse before
    # Python 3.12 does not hand those to the positional list, so collect them here
    args, extra = parser.parse_known_args(argv)
    if extra:
        if getattr(args, "command", None) == "run" and all("=" in e and not e.startswith("-") for e in extra):
            args.inputs = list(args.inputs) + list(extra)
        else:
            parser.error("unrecognized arguments: " + " ".join(extra))
    try:
        return int(args.func(args))
    except KeyError as exc:
        print(f"ssc-calc: {exc.args[0] if exc.args else exc}", file=sys.stderr)
        return 2
    except MemoryError as exc:  # e.g. no memory even for the engine's answer (core.catalog, search, schema)
        print(f"ssc-calc: {str(exc) or 'there is not enough memory'}", file=sys.stderr)
        return 1
    except BrokenPipeError:  # e.g. piped into `head`
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
