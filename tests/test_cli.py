"""The ssc-calc command-line tool and the plain-text formatting helpers."""
from __future__ import annotations

import csv
import io
import json

import pytest

from ssc import cli
from ssc.formatting import ReportOptions, fmt_cell, smart, text_table


def test_list_and_search(capsys):
    assert cli.main(["list"]) == 0
    out = capsys.readouterr().out
    assert "means.two.t" in out and "Group-Sequential" in out
    assert cli.main(["list", "--search", "logrnk"]) == 0  # misspelled
    out = capsys.readouterr().out
    assert "surv.logrank" in out and "means.two.t" not in out


def test_search_command(capsys):
    assert cli.main(["search", "kruskall", "walis"]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0] == 'Showing results for "kruskal wallis"' and "nonpar.kruskal" in lines[1]
    assert cli.main(["search", "qqqzzz"]) == 1
    assert cli.main(["search", "NI", "proportions", "--json", "--limit", "2"]) == 0
    res = json.loads(capsys.readouterr().out)
    assert [h["id"] for h in res["hits"]][0] == "prop.two.ni" and len(res["hits"]) == 2


def test_show_describes_a_procedure(capsys):
    assert cli.main(["show", "prop.two"]) == 0
    out = capsys.readouterr().out
    assert "Solve for:" in out and "fisher" in out and "Example:" in out


def test_run_prints_table_and_statement(capsys):
    rc = cli.main(["run", "means.two.t", "--solve", "n", "power=0.8 0.9", "mu1=1", "mu2=0", "sigma=1"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "Target Power" in out
    assert "Group sample sizes of 17 and 17" in out and "of 23 and 23" in out


def test_run_csv_and_json(capsys):
    assert cli.main(["run", "means.two.t", "--csv", "n1=17", "--solve", "power"]) == 0
    rows = list(csv.reader(io.StringIO(capsys.readouterr().out)))
    assert rows[0][0] == "Power" and abs(float(rows[1][0]) - 0.80703671514740) < 1e-10
    assert cli.main(["run", "gs.means", "--example", "--json", "--no-curve"]) == 0
    res = json.loads(capsys.readouterr().out)
    assert res["ok"] and res["rows"][0]["_tables"][0]["title"] == "Boundaries"


def test_run_reports_bad_input(capsys):
    assert cli.main(["run", "means.two.t", "foo=1"]) == 2
    assert "unknown parameter" in capsys.readouterr().err
    assert cli.main(["run", "means.two.t", "power=1.5"]) == 1
    assert "power" in capsys.readouterr().err
    assert cli.main(["run", "means.two.t", "--solve", "nonsense"]) == 2
    assert cli.main(["show", "no.such.procedure"]) == 2


def test_no_memory_for_the_engine_s_answer_is_an_error(capsys, monkeypatch):
    from ssc import _native, core

    lib = _native.lib()
    monkeypatch.setattr(lib, "ssc_run_json", lambda request, token: None)
    assert cli.main(["run", "means.two.t"]) == 1
    assert "Error: There is not enough memory for this calculation" in capsys.readouterr().err
    monkeypatch.setattr(lib, "ssc_procedure_json", lambda procedure: None)
    assert cli.main(["run", "means.two.t"]) == 1  # a known procedure without memory for its schema
    assert capsys.readouterr().err == "ssc-calc: there is not enough memory for the engine's answer\n"
    assert cli.main(["show", "no.such.procedure"]) == 2
    assert capsys.readouterr().err == "ssc-calc: unknown procedure: no.such.procedure\n"
    monkeypatch.setattr(lib, "ssc_catalog_json", lambda: None)
    assert cli.main(["list"]) == 1
    assert capsys.readouterr().err == "ssc-calc: there is not enough memory for the engine's answer\n"

    def no_memory():
        raise MemoryError  # Python's own, without a message

    monkeypatch.setattr(core, "catalog", no_memory)
    assert cli.main(["list"]) == 1
    assert capsys.readouterr().err == "ssc-calc: there is not enough memory\n"


def test_eval(capsys):
    assert cli.main(["eval", "pnorm", "0", "1"]) == 0
    assert float(capsys.readouterr().out) == 0.5


def test_number_formatting():
    o = ReportOptions(prob_decimals=4)
    assert fmt_cell(0.80703671, "prob", o) == "0.8070"
    assert fmt_cell(12345, "int", o) == "12,345"
    assert fmt_cell(0.1, "pct", o) == "10%"
    assert fmt_cell(None, "num", o) == "—"
    assert fmt_cell(float("nan"), "num", o) == "—"
    assert smart(1.5) == "1.5" and smart(2.0) == "2" and smart(1e-7) == "1e-07"


def test_text_table_alignment():
    t = text_table(["Name", "N"], [["a", "5"], ["bb", "123"]], left=[0]).splitlines()
    assert t[0] == "Name    N"
    assert t[2] == "a       5"
    assert t[3] == "bb    123"


def test_a_run_that_computes_nothing_fails(capsys):
    """Every scenario refused (a target power below alpha): the result is printed, with a failing exit code."""
    assert cli.main(["run", "means.two.t", "power=0.04"]) == 1
    assert "must exceed alpha" in capsys.readouterr().err
    assert cli.main(["run", "means.two.t", "--csv", "power=0.04"]) == 1
    assert cli.main(["run", "means.two.t", "power=0.04 0.8"]) == 0  # one scenario computed


def test_an_input_that_is_not_used_is_said(capsys):
    assert cli.main(["run", "means.two.t", "--solve", "power", "n1=40", "n2=99", "mu1=0.5"]) == 0
    err = capsys.readouterr().err
    assert "ignored n2: not used with these choices (used when" in err and "alloc = sep" in err


def test_decimals_must_be_sensible(capsys):
    for bad in ("-1", "100000", "x"):
        with pytest.raises(SystemExit) as exc:
            cli.main(["run", "means.two.t", "--decimals", bad])
        assert exc.value.code == 2
    assert "decimals" in capsys.readouterr().err


def test_text_that_is_not_valid_unicode(capsys):
    assert cli.main(["run", "means.two.t", "mu1=\udcff"]) != 0
    out = capsys.readouterr()
    assert "mu1" in out.out + out.err
    assert cli.main(["show", "\udcff"]) != 0
