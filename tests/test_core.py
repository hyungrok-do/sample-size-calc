"""Python bindings of the C++ engine."""
from __future__ import annotations

import json
import re
import threading
import time
from pathlib import Path

import pytest

import ssc


def test_version_and_threads():
    assert ssc.version()
    assert ssc.hardware_threads() >= 1


def test_every_copy_of_the_version_agrees():
    # a release tag must match it too
    assert ssc.version() == ssc.__version__  # core/src/capi.cpp and python/ssc/__init__.py
    root = Path(__file__).resolve().parents[1]
    sources = [("pyproject.toml", r'^version = "(.+)"$'), ("CMakeLists.txt", r"^project\(\S+ VERSION (\S+)")]
    for name, pattern in sources:
        if (root / name).exists():  # in the source tree
            match = re.search(pattern, (root / name).read_text(encoding="utf-8"), re.M)
            assert match and match.group(1) == ssc.__version__, name


def test_catalog_lists_procedures_by_category():
    cat = ssc.catalog()
    ids = {p["id"] for p in cat["procedures"]}
    assert {"means.two.t", "prop.two", "corr.one", "surv.logrank", "gs.means"} <= ids
    for p in cat["procedures"]:
        assert p["category"] in cat["categories"]
        assert p["name"] and p["summary"]


def test_schema_describes_parameters():
    s = ssc.schema("means.two.t")
    assert s["id"] == "means.two.t"
    assert {"n", "power"} <= {o["value"] for o in s["solve_for"]}
    ids = [p["id"] for p in s["params"]]
    assert "alpha" in ids and "sigma" in ids
    with pytest.raises(KeyError):
        ssc.schema("no.such.procedure")


def test_two_sample_t_test_matches_power_t_test():
    res = ssc.run("means.two.t", "n", {"power": [0.8, 0.9], "mu1": 1, "mu2": 0, "sigma": 1})
    assert res["ok"]
    assert [r["n1"] for r in res["rows"]] == [17, 23]
    # the value of power.t.test in R (strict two-sided power)
    assert res["rows"][0]["power"] == pytest.approx(0.80703671514740, rel=1e-10)
    assert res["rows"][0]["_statement"].startswith("Group sample sizes of 17 and 17")


def test_ranges_expand_into_scenarios_and_curves():
    res = ssc.run("means.two.t", "power", {"n1": "10 to 50 by 10", "mu1": "0.5 1"})
    assert res["ok"] and len(res["rows"]) == 10
    powers = [r["power"] for r in res["rows"] if r["mu1"] == 0.5]
    assert powers == sorted(powers)
    assert res["plots"], "a power curve is computed by default"


def test_input_errors_are_reported():
    res = ssc.run("means.two.t", "n", {"power": 2})
    assert not res["ok"]
    assert "power" in res["input_errors"]


def test_row_errors_do_not_stop_other_scenarios():
    res = ssc.run("means.two.t", "n", {"mu1": "0 1", "mu2": 0})
    assert res["ok"]
    assert "_error" in res["rows"][0] and "_error" not in res["rows"][1]


def test_evaluate_distribution_functions():
    assert ssc.evaluate("pnorm", 1.959963984540054, 1) == pytest.approx(0.975, abs=1e-15)
    assert ssc.evaluate("qt", 0.975, 10, 1) == pytest.approx(2.2281388519862739, rel=1e-14)  # t quantile, 10 df
    with pytest.raises(ValueError):
        ssc.evaluate("no_such_function", 1.0)


def test_thread_count_does_not_change_results():
    inputs = {"n1": "20 to 200 by 20", "p1": "0.5 0.6", "p2": 0.4, "test": "fisher"}
    a = ssc.run("prop.two", "power", inputs, threads=1)
    b = ssc.run("prop.two", "power", inputs, threads=4)
    assert [r["power"] for r in a["rows"]] == [r["power"] for r in b["rows"]]


def test_cancel_from_another_thread():
    heavy = {"n1": "1000 to 7000 by 1", "p1": "0.5", "p2": 0.45, "test": "fisher"}
    out = {}
    with ssc.Token() as token:
        worker = threading.Thread(target=lambda: out.update(res=ssc.run("prop.two", "power", heavy, token=token)))
        started = time.monotonic()
        worker.start()
        time.sleep(0.3)
        token.cancel()
        worker.join(60)
        assert not worker.is_alive()
        assert time.monotonic() - started < 30
    assert out["res"].get("cancelled") is True
    assert not out["res"]["ok"]


def test_closing_a_token_waits_for_a_cancel_under_way(monkeypatch):
    # a time limit firing as the calculation ends: the cancel had read the token and was about to write to it
    # when the other thread freed it (heap corruption on the MCP server)
    from ssc import _native

    lib = _native.lib()
    real_cancel, real_free = lib.ssc_token_cancel, lib.ssc_token_free
    inside, events = threading.Event(), []

    def slow_cancel(ptr):
        inside.set()
        time.sleep(0.2)  # releases the GIL, as the foreign call does
        events.append("cancelled")
        return real_cancel(ptr)

    def free(ptr):
        events.append("freed")
        return real_free(ptr)

    monkeypatch.setattr(lib, "ssc_token_cancel", slow_cancel)
    monkeypatch.setattr(lib, "ssc_token_free", free)
    token = ssc.Token()
    canceller = threading.Thread(target=token.cancel)
    canceller.start()
    inside.wait(5)
    token.close()
    canceller.join(5)
    assert events == ["cancelled", "freed"]


def test_a_token_closed_during_a_calculation_is_freed_when_it_returns():
    out = {}
    token = ssc.Token()
    worker = threading.Thread(target=lambda: out.update(res=ssc.run(
        "prop.two", "power", {"n1": "1000 to 3000 by 1", "test": "fisher"}, token=token)))
    worker.start()
    time.sleep(0.2)
    token.cancel()
    token.close()  # while the engine still reads it
    worker.join(60)
    assert out["res"].get("cancelled") is True
    with pytest.raises(ValueError):
        ssc.run("means.two.t", "n", {}, token=token)


def test_every_procedure_runs_its_example():
    for p in ssc.catalog()["procedures"]:
        ex = ssc.schema(p["id"])["example"]
        res = ssc.run(p["id"], ex["solve_for"], ex["inputs"], power_curve=False)
        assert res["ok"], (p["id"], res.get("error"), res.get("input_errors"))
        assert not [r for r in res["rows"] if "_error" in r], p["id"]


def test_search_ranks_corrects_and_limits():
    res = ssc.search("two sample t test")
    assert res["hits"][0]["id"] == "means.two.t"
    assert {"id", "name", "category", "score", "terms"} <= set(res["hits"][0])
    typo = ssc.search("equivalance two means")
    assert typo["corrected"] == "equivalence two means" and typo["hits"][0]["id"] == "means.two.equiv"
    assert ssc.search("κ")["hits"][0]["id"] == "agree.kappa"
    assert ssc.search("NI proportions")["hits"][0]["id"] == "prop.two.ni"
    assert len(ssc.search("test", 3)["hits"]) == 3
    assert len(ssc.search("test", 0)["hits"]) > 20
    assert ssc.search("qqqzzz")["hits"] == [] and ssc.search("")["hits"] == []
    assert ssc.search("survival phase ii")["partial"] is True



# ----------------------------------------------------------------------------- found by black-box fuzzing
def test_large_stepped_wedges_neither_overflow_nor_run_on():
    row = ssc.run("cluster.sw", "power", {"g": "1000000000"})["rows"][0]
    assert row["clusters"] == 4e9 and row["N"] == 4e11          # was a negative 32-bit wrap-around
    assert ssc.run("cluster.sw", "power", {"steps": "50000"})["input_errors"] == {
        "steps": "values must be at most 1000"}
    assert ssc.run("means.anova", "f", {"estype": "f", "n": "3", "k": "100000000"})["input_errors"] == {
        "k": "values must be at most 1000"}


def test_groups_from_a_ratio_or_a_percentage_meet_the_smallest_group():
    rows = [ssc.run(*args)["rows"][0] for args in (
        ("means.two.ni", "power", {"alloc": "ratio", "n1": "7", "r": "0.1", "test": "t"}),
        ("prop.two", "power", {"alloc": "ratio", "n1": "7", "r": "0.1"}),
        ("means.two.t", "power", {"alloc": "pct", "ntotal": "40", "pct1": "99"}),
        ("means.two.cluster", "power", {"k1": "5", "kr": "0.1"}),
        ("diag.auc", "power", {"n1": "5", "ratio": "0.1"}))]
    assert all(r.get("_error") for r in rows), rows
    assert rows[0]["_error"] == "group 2 would have N2 = 1 (R × N1, rounded up); each group needs at least 2"
    assert ssc.run("means.two.ni", "power", {"alloc": "ratio", "n1": "20", "r": "0.5"})["rows"][0]["n2"] == 10


def test_a_group_from_a_tiny_ratio_has_one_subject():
    # n × the ratio rounded up: a slack of 1e-9 for the rounding of the product left group 2 empty below 1e-9
    row = ssc.run("means.two.t", "power", {"alloc": "ratio", "n1": "10", "r": "5e-11"})["rows"][0]
    assert row["_error"] == "group 2 would have N2 = 1 (R × N1, rounded up); each group needs at least 2"
    tiny = ssc.run("surv.logrank", "power", {"alloc": "ratio", "n1": "10", "r": "1e-12"})["rows"][0]
    typed = ssc.run("surv.logrank", "power", {"alloc": "sep", "n1": "10", "n2": "1"})["rows"][0]
    assert tiny["n2"] == 1 and tiny["power"] == typed["power"]


def test_no_memory_even_for_the_answer_is_said_as_the_engine_says_it(monkeypatch):
    # the engine returns NULL only when not even its short answer for running out of memory can be allocated
    from ssc import _native

    lib = _native.lib()
    monkeypatch.setattr(lib, "ssc_run_json", lambda request, token: None)
    res = ssc.run("means.two.t", "n", {})
    assert res == {"ok": False, "error": "There is not enough memory for this calculation: use fewer combinations of "
                                         "input values or shorter lists."}
    # NULL for a schema: an unknown id, or a known one without memory for it
    monkeypatch.setattr(lib, "ssc_procedure_json", lambda procedure: None)
    with pytest.raises(MemoryError, match="not enough memory"):
        ssc.schema("means.two.t")
    with pytest.raises(KeyError, match="unknown procedure: no.such.procedure"):
        ssc.schema("no.such.procedure")
    monkeypatch.setattr(lib, "ssc_catalog_json", lambda: None)
    with pytest.raises(MemoryError):
        ssc.catalog()
    with pytest.raises(MemoryError):
        ssc.schema("no.such.procedure")  # whether it is known cannot be found out either


def test_text_that_is_not_valid_unicode_is_an_error_not_a_crash():
    res = ssc.run("means.two.t", "n", {"mu1": "\ud800"})
    assert not res["ok"] and "mu1" in res["input_errors"]
    assert ssc.run("\ud800", "n", {})["error"] == "unknown procedure: �"
    with pytest.raises(KeyError):
        ssc.schema("\ud800")
    assert isinstance(ssc.search("\ud800"), dict)
    from ssc import _native
    raw = _native.take_string(_native.lib().ssc_run_json(b'{"procedure":"\xc3","solve_for":"n","inputs":{}}', None))
    assert json.loads(raw)["error"] == "unknown procedure: �"     # valid UTF-8 out of invalid bytes in


def test_declared_limits_are_the_enforced_ones_and_power_stays_a_probability():
    seed = next(p for p in ssc.schema("sim.sw")["params"] if p["id"] == "seed")
    assert seed["max"] == 1e9
    assert ssc.run("sim.sw", "power", {"seed": "1000000000"}, approximate=True)["ok"]
    assert ssc.run("sim.multiarm", "power", {"n": "50000"}, approximate=True)["rows"][0]["power"] == 1
