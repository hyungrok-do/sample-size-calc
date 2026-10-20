"""Who depends on whom: the engine is used by the command line and the MCP server alike, and nothing goes
through the MCP server, which is an optional outer layer."""
from __future__ import annotations

import ast
import os
import subprocess
import sys
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1] / "python" / "ssc"


def imports_of(path: Path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            yield from (a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            yield ("." * node.level) + (node.module or "")
            yield from (("." * node.level) + (node.module + "." if node.module else "") + a.name
                        for a in node.names)


def test_only_the_mcp_server_uses_the_mcp_sdk_and_nothing_uses_the_server():
    for path in PACKAGE.rglob("*.py"):
        if path.name == "mcp_server.py":
            continue
        for name in imports_of(path):
            assert not (name == "mcp" or name.startswith("mcp.")), (path, name)
            assert "mcp_server" not in name, (path, name)
    server = set(imports_of(PACKAGE / "mcp_server.py"))
    assert not any(n.startswith(".cli") for n in server), server  # the engine, not the command line


def test_the_package_runs_without_the_mcp_package():
    """With the MCP SDK made unimportable, the engine and the command line load and calculate, and ssc-mcp explains
    what to install."""
    code = """
import sys
sys.modules["mcp"] = None
import ssc, ssc.cli, ssc.inputs, ssc.docx
row = ssc.run("means.two.t", "n", {"mu1": 0.5})["rows"][0]
assert row["n1"] == 64, row
assert ssc.cli.main(["run", "means.two.t", "mu1=0.5", "--brief"]) == 0
assert "ssc.mcp_server" not in sys.modules
from ssc import mcp_server
assert mcp_server.main(["--transport", "stdio"]) == 1
print("ok")
"""
    import ssc  # the package these tests use: installed (with its compiled core) or the source tree (conftest.py)
    env = {**os.environ, "PYTHONPATH": str(Path(ssc.__file__).resolve().parents[1]) + os.pathsep +
           os.environ.get("PYTHONPATH", "")}
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env, timeout=120)
    assert out.returncode == 0 and out.stdout.strip().endswith("ok"), out.stdout + out.stderr
    assert "pip install" in out.stderr
