"""Number formatting and plain-text/CSV rendering of results (no Qt needed).

Used by the ``ssc-calc`` command-line tool and its CSV and Word exports.
"""
from __future__ import annotations

import csv
import io
import math
from dataclasses import dataclass
from typing import Any, List, Optional, Sequence


@dataclass
class ReportOptions:
    prob_decimals: int = 4
    statements: str = "first"  # "first" | "all" | "none"
    show_definitions: bool = True
    show_references: bool = True
    show_notes: bool = True
    show_dropout: bool = True
    show_plots: bool = True


def smart(v: float, max_dec: int = 4) -> str:
    """Compact human-readable number: integers without decimals, up to
    ``max_dec`` decimals otherwise, scientific notation for extremes."""
    if not math.isfinite(v):
        return "—"
    if float(v).is_integer() and abs(v) < 1e15:
        return f"{v:,.0f}" if abs(v) >= 10000 else f"{v:.0f}"
    a = abs(v)
    if a >= 1e6 or a < 1e-4:
        return f"{v:.4g}"
    s = f"{v:.{max_dec}f}".rstrip("0").rstrip(".")
    return s if s not in ("-0", "") else "0"


def fmt_cell(v: Any, fmt: str, opts: Optional[ReportOptions] = None) -> str:
    """Display text of a result cell according to its column format."""
    if v is None:
        return "—"
    if isinstance(v, str):
        return v
    if isinstance(v, bool):
        return "Yes" if v else "No"
    try:
        x = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not math.isfinite(x):
        return "—"
    if fmt == "prob":
        return f"{x:.{(opts or ReportOptions()).prob_decimals}f}"
    if fmt == "int":
        return f"{x:,.0f}"
    if fmt == "pct":
        return smart(100 * x, 2) + "%"
    return smart(x)


def raw_cell(v: Any) -> str:
    """Full-precision text for CSV / clipboard export."""
    if v is None:
        return ""
    if isinstance(v, float):
        return repr(v) if math.isfinite(v) else ""
    return str(v)


def table_tsv(result: dict, raw: bool = True, opts: Optional[ReportOptions] = None) -> str:
    """Numeric results as tab-separated text (spreadsheet friendly)."""
    cols = result.get("columns", [])
    lines = ["\t".join(c["label"] for c in cols)]
    for r in result.get("rows", []):
        if r.get("_error"):
            lines.append("\t".join([f"Not computed: {r['_error']}"] + [""] * (len(cols) - 1)))
            continue
        lines.append("\t".join(raw_cell(r.get(c["key"])) if raw else fmt_cell(r.get(c["key"]), c.get("fmt", "num"), opts)
                               for c in cols))
    return "\n".join(lines) + "\n"


def table_csv(result: dict) -> str:
    """Numeric results as CSV with full precision and an error column."""
    cols = result.get("columns", [])
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow([c["label"] for c in cols] + ["error"])
    for r in result.get("rows", []):
        w.writerow([raw_cell(r.get(c["key"])) for c in cols] + [r.get("_error", "")])
    return buf.getvalue()


def text_table(headers: Sequence[str], rows: Sequence[Sequence[str]], *, indent: str = "", left: Sequence[int] = ()) -> str:
    """Monospaced table: numbers right-aligned, columns listed in ``left`` left-aligned."""
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            if i < len(widths):
                widths[i] = max(widths[i], len(cell))
    def line(cells: Sequence[str]) -> str:
        parts: List[str] = []
        for i, cell in enumerate(cells):
            parts.append(cell.ljust(widths[i]) if i in left else cell.rjust(widths[i]))
        return indent + "  ".join(parts).rstrip()
    out = [line(headers), indent + "  ".join("─" * w for w in widths)]
    out.extend(line(list(r) + [""] * (len(headers) - len(r))) for r in rows)
    return "\n".join(out)


def result_table_text(columns: Sequence[dict], rows: Sequence[dict], opts: Optional[ReportOptions] = None,
                      *, indent: str = "") -> str:
    """The main result table of a calculation as monospaced text."""
    headers = [c["label"] for c in columns]
    left = [i for i, c in enumerate(columns) if c.get("fmt") == "text"]
    body = []
    for r in rows:
        if r.get("_error"):
            body.append([f"not computed: {r['_error']}"])
        else:
            body.append([fmt_cell(r.get(c["key"]), c.get("fmt", "num"), opts) for c in columns])
    if any(r.get("_error") for r in rows):
        # error rows span the table: print them left-aligned after the header
        ok = [b for b, r in zip(body, rows) if not r.get("_error")]
        text = text_table(headers, ok, indent=indent, left=left)
        errs = [indent + b[0] for b, r in zip(body, rows) if r.get("_error")]
        return text + ("\n" + "\n".join(errs) if errs else "")
    return text_table(headers, body, indent=indent, left=left)


def sub_table_text(table: dict, opts: Optional[ReportOptions] = None, *, indent: str = "") -> str:
    """A per-row detail table (``_tables`` entry) as monospaced text."""
    cols = table.get("columns", [])
    rows = [[fmt_cell(v, cols[i].get("fmt", "num") if i < len(cols) else "num", opts) for i, v in enumerate(r)]
            for r in table.get("rows", [])]
    return text_table([c["label"] for c in cols], rows, indent=indent)
