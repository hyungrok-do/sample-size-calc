"""The MCP server (ssc-mcp): its tools, what a model reads back, and the limits of a hosted server."""
from __future__ import annotations

import json

import pytest

pytest.importorskip("mcp")
anyio = pytest.importorskip("anyio")

from mcp import Client  # noqa: E402

from ssc import mcp_server  # noqa: E402
from ssc.mcp_server import Limits, build_server, calculate, describe_method, search_methods  # noqa: E402


def call(tool: str, args: dict, limits: Limits = Limits()) -> dict:
    async def go():
        async with Client(build_server(limits)) as client:
            res = await client.call_tool(tool, args)
            return json.loads(res.content[0].text)

    return anyio.run(go)


def test_three_read_only_tools_with_titles():
    async def go():
        async with Client(build_server()) as client:
            return (await client.list_tools()).tools, client.instructions

    tools, instructions = anyio.run(go)
    assert [t.name for t in tools] == ["search_methods", "describe_method", "calculate"]
    for t in tools:
        hints = t.annotations.model_dump(by_alias=True, exclude_none=True)
        assert t.title and hints["title"] == t.title
        assert hints["readOnlyHint"] is True and hints["destructiveHint"] is False
        assert hints["openWorldHint"] is False and hints["idempotentHint"] is True
        assert len(t.name) <= 64 and t.description
    assert "search_methods" in instructions and "never compute a sample size yourself" in instructions


def test_find_describe_calculate():
    found = call("search_methods", {"description": "two independent means", "limit": 3})
    assert found["methods"][0]["id"] == "means.two.t"
    assert len(search_methods("")["methods"]) == 66
    method = call("describe_method", {"method_id": "means.two.t"})
    n2 = next(i for i in method["inputs"] if i["id"] == "n2")
    assert "alloc = sep" in n2["used_when"] and "whole numbers" in n2["values"]
    alloc = next(i for i in method["inputs"] if i["id"] == "alloc")
    assert {"value": "sep", "label": "Enter N1 and N2 individually", "available_when": "solving for ≠ n"} in alloc["choices"]
    res = call("calculate", {"method_id": "means.two.t", "solve_for": "n",
                             "inputs": {"mu1": 5, "mu2": 0, "sigma": 12, "power": [0.8, 0.9], "dropout": 0.1}})
    assert res["ok"] and [s["values"]["N1"] for s in res["scenarios"]] == [92, 123]
    assert res["scenarios"][1]["statement"].startswith("Group sample sizes of 123 and 123")
    assert res["allowing_for_dropout"][0]["N1'"] == 103


def test_what_was_not_used_is_said():
    res = calculate("means.two.t", "power", {"n1": 40, "n2": 99, "mu1": 0.5})
    assert res["ok"] and res["scenarios"][0]["values"]["N2"] == 40
    assert any(w.startswith("ignored n2") and "alloc = sep" in w for w in res["warnings"])
    res = calculate("means.two.t", "n", {"alloc": "sep", "mu1": 0.5})
    assert any("alloc=sep is not available here; used alloc=equal" in w for w in res["warnings"])
    res = calculate("means.two.t", "power", {"nn": 40})
    assert not res["ok"] and "unknown input(s) for means.two.t: nn" in res["error"]
    assert "n1" in res["inputs_of_this_method"]


def test_mistakes_get_directions():
    assert describe_method("means.twot")["did_you_mean"][0] == "means.two.t"
    assert "choose one of n, power, mu1" in calculate("means.two.t", "width", {})["error"]
    res = calculate("means.two.t", "n", {"mu1": 0.5, "power": 0.04})
    assert "must exceed alpha" in res["scenarios"][0]["error"]
    res = calculate("means.two.t", "power", {"alternative": "greater"})
    assert not res["ok"] and "use one of two, upper, lower" in res["input_errors"]["alternative"]


def test_simulations_approximate_unless_asked():
    res = calculate("sim.mmrm", "n", {"delta": 0.4})
    assert res["ok"] and res["simulated"] is False and any("approximation" in w for w in res["warnings"])
    res = calculate("sim.means.robust", "power", {"n1": 30, "reps": 200}, simulate=True)
    assert res["ok"] and res["simulated"] is True and res["scenarios"][0]["values"]["Studies"] == 200
    capped = calculate("sim.means.robust", "power", {"n1": 30, "reps": 9000}, simulate=True,
                       limits=Limits(max_reps=300))
    assert capped["scenarios"][0]["values"]["Studies"] == 300 and any("capped at 300" in w for w in capped["warnings"])


def test_a_hosted_server_keeps_requests_in_proportion():
    res = calculate("means.two.t", "n", {"mu1": "0.01 to 1 by 0.001"}, limits=Limits(max_scenarios=200))
    assert not res["ok"] and "limit 200" in res["error"]
    slow = calculate("sim.mmrm", "n", {"delta": 0.2, "reps": 5000}, simulate=True, limits=Limits(time_limit=0.05))
    assert not slow["ok"] and "longer than 0.05 seconds" in slow["error"]
    assert Limits.for_transport("http").max_scenarios == 200 and Limits.for_transport("stdio").max_reps is None


def test_health_and_domain_verification_routes(monkeypatch):
    from starlette.testclient import TestClient

    app = build_server().streamable_http_app()
    with TestClient(app) as client:
        assert client.get("/healthz").text == "ok"
        assert client.get("/.well-known/openai-apps-challenge").status_code == 404
        monkeypatch.setenv("SSC_OPENAI_APPS_CHALLENGE", "token-123")
        page = client.get("/.well-known/openai-apps-challenge")
        assert page.status_code == 200 and page.text == "token-123"


def test_inputs_as_a_model_writes_them():
    assert mcp_server._as_input([0.8, 0.9]) == "0.8 0.9"
    assert mcp_server._as_input([[0, 0.5, 1], [0, 0, 1]]) == "0 0.5 1; 0 0 1"
    assert mcp_server._as_input(True) == "true" and mcp_server._as_input(5.0) == "5"
    res = calculate("means.anova", "n", {"means": [[0, 0.5, 1], [0, 0, 1]]})
    assert [s["values"]["n"] for s in res["scenarios"]] == [21, 16]


def test_the_replicate_cap_reads_every_way_of_writing_a_number():
    lim = Limits.for_transport("http")
    for reps in ("20000,", "2e4,", "30000"):
        out = calculate("sim.means.robust", "power", {"reps": reps, "n1": "10"}, simulate=True, limits=lim)
        assert out["scenarios"][0]["values"]["Studies"] == lim.max_reps, reps
        assert out["warnings"] == [f"replicates capped at {lim.max_reps} on this server"]
