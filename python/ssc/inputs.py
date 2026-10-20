"""Inputs in words, from a procedure's description and a result (no Qt, no MCP): when an input applies, and
what a run did not use of the inputs it was given. Shared by the command line and the MCP server."""
from __future__ import annotations

from typing import Dict, List


def applies(spec: dict, values: Dict[str, str], solve: str) -> bool:
    """Whether an input (or a choice's option) applies, given the other choices and what is solved for (the
    schema's visibility rule)."""
    alternatives = spec.get("visible") or []
    if not alternatives:
        return True
    for clauses in alternatives:
        ok = True
        for c in clauses:
            current = solve if c.get("key") == "$solve" else values.get(c.get("key", ""), "")
            hit = current in (c.get("values") or [])
            if hit == bool(c.get("negate")):
                ok = False
                break
        if ok:
            return True
    return False


def when_used(rule: List[List[dict]]) -> str:
    """A visibility rule in words: 'solving for ≠ n and alloc = sep or solving for = n and alloc = n2'."""
    def clause(c: dict) -> str:
        key = "solving for" if c["key"] == "$solve" else c["key"]
        op = "not " if c.get("negate") else ""
        return f"{key} {op}in {{{', '.join(c['values'])}}}" if len(c["values"]) > 1 else \
            f"{key} {'≠' if c.get('negate') else '='} {c['values'][0]}"
    return " or ".join(" and ".join(clause(c) for c in conj) for conj in rule)


def unused_inputs(res: dict, schema: dict) -> List[str]:
    """What the run did not use of the given inputs, in words: an input the chosen options leave out (with
    when it applies), or a choice not available here and the one used instead."""
    params = {p["id"]: p for p in schema.get("params", [])}
    out = []
    for pid in res.get("ignored_inputs") or []:
        rule = params.get(pid, {}).get("visible") or []
        out.append(f"ignored {pid}: not used with these choices" + (f" (used when {when_used(rule)})" if rule else ""))
    for pid, r in (res.get("replaced_choices") or {}).items():
        out.append(f"{pid}={r.get('given')} is not available here; used {pid}={r.get('used')}")
    return out
