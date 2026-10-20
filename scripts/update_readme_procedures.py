"""Regenerate the procedure list in README.md from the engine's catalog.

    PYTHONPATH=python python scripts/update_readme_procedures.py
"""
from __future__ import annotations

import re
from pathlib import Path

import ssc

START, END = "<!-- procedures:start -->", "<!-- procedures:end -->"


def table() -> str:
    cat = ssc.catalog()
    order = {c: i for i, c in enumerate(cat["categories"])}
    procs = sorted(cat["procedures"], key=lambda p: (order[p["category"]], p.get("subcategory", ""), p["name"]))
    lines, current = [], None
    for p in procs:
        if p["category"] != current:
            current = p["category"]
            lines += ["", f"**{current}**", "", "| Procedure | Id | Solve for |", "|---|---|---|"]
        solves = ", ".join(o["label"] for o in ssc.schema(p["id"])["solve_for"])
        lines.append(f"| {p['name']} | `{p['id']}` | {solves} |")
    return f"{len(procs)} procedures in {len(cat['categories'])} categories.\n" + "\n".join(lines)


def main() -> None:
    path = Path(__file__).resolve().parent.parent / "README.md"
    text = path.read_text(encoding="utf-8")
    new = re.sub(re.escape(START) + ".*?" + re.escape(END), START + "\n" + table() + "\n" + END, text, flags=re.S)
    count = len(ssc.catalog()["procedures"])
    new = re.sub(r"\*\*\d+ procedures\*\*", f"**{count} procedures**", new)
    path.write_text(new, encoding="utf-8")
    print(f"README.md: {count} procedures")


if __name__ == "__main__":
    main()
