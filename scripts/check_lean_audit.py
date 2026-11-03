#!/usr/bin/env python3
"""Reject unexpected axioms in the explicit Lean theorem audit."""
import re
import sys
from pathlib import Path

text = Path(sys.argv[1]).read_text()
allowed = {"propext", "Classical.choice", "Quot.sound"}
declarations = re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]", text)
if len(declarations) < 10:
    raise SystemExit("Incomplete Lean axiom audit")
for declaration, raw in declarations:
    axioms = {x.strip() for x in raw.split(",") if x.strip()}
    unexpected = axioms - allowed
    if unexpected:
        raise SystemExit(f"{declaration}: unexpected axioms {sorted(unexpected)}")
if "sorryAx" in text or "error:" in text:
    raise SystemExit("Lean audit contains a placeholder or an error")
print(f"{len(declarations)} declarations audited; only standard classical axioms.")
