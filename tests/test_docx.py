"""The report as a Word document: a well-formed package with the same content as the HTML report."""
from __future__ import annotations

import datetime as dt
import io
import struct
import zipfile
import zlib
import xml.etree.ElementTree as ET

import pytest

from ssc import cli, core
from ssc.docx import build_docx, png_size, write_docx
from ssc.formatting import ReportOptions

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
WHEN = dt.datetime(2026, 1, 2, 3, 4)


def png(width: int, height: int) -> bytes:
    """A blank PNG image of the given size."""
    def chunk(kind: bytes, data: bytes) -> bytes:
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data))
    raw = b"".join(b"\x00" + b"\xff" * (3 * width) for _ in range(height))
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw)) + chunk(b"IEND", b""))


def parts(data: bytes) -> dict:
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        assert z.testzip() is None
        return {name: z.read(name) for name in z.namelist()}


def document(data: bytes) -> ET.Element:
    """word/document.xml, after checking that every XML part parses."""
    p = parts(data)
    for name, body in p.items():
        if name.endswith(".xml") or name.endswith(".rels"):
            ET.fromstring(body)
    return ET.fromstring(p["word/document.xml"])


def text(root: ET.Element) -> str:
    return "\n".join("".join(t.text or "" for t in p.iter(W + "t")) for p in root.iter(W + "p"))


def test_a_report_with_statements_dropout_and_references():
    res = core.run("means.two.t", "n", {"power": "0.8 0.9", "mu1": 1, "mu2": 0, "sigma": 1, "dropout": 0.1},
                   power_curve=False)
    assert res["ok"]
    data = build_docx(res, ReportOptions(statements="all"), now=WHEN)
    p = parts(data)
    assert {"[Content_Types].xml", "_rels/.rels", "word/document.xml", "word/styles.xml",
            "word/_rels/document.xml.rels", "docProps/core.xml", "docProps/app.xml"} <= set(p)
    assert b"<dc:title>Two-Sample T-Tests" in p["docProps/core.xml"]
    root = document(data)
    words = text(root)
    assert words.startswith(res["name"])
    for title in ("Results", "In Plain Words", "Allowing for Dropout", "Definitions", "References"):
        assert title in words
    assert "Group sample sizes of 17 and 17" in words and "of 23 and 23" in words
    assert "Generated 2026-01-02 03:04 by Sample Size Calculator." in words
    tables = list(root.iter(W + "tbl"))
    results = tables[0]
    rows = list(results.iter(W + "tr"))
    assert len(rows) == 1 + len(res["rows"])
    assert rows[0].find(f"{W}trPr/{W}tblHeader") is not None  # the header repeats on each page
    head = ["".join(t.text for t in c.iter(W + "t")) for c in rows[0].iter(W + "tc")]
    assert head == [c["label"] for c in res["columns"]]


def test_wide_tables_turn_the_page_and_extra_tables_follow():
    narrow = core.run("corr.one", "n", {"r1": 0.3}, power_curve=False)
    assert len(narrow["columns"]) <= 8
    assert document(build_docx(narrow, now=WHEN)).find(f"{W}body/{W}sectPr/{W}pgSz").get(W + "orient") is None
    res = core.run("gs.means", core.schema("gs.means")["example"]["solve_for"],
                   core.schema("gs.means")["example"]["inputs"], power_curve=False)
    assert res["ok"] and len(res["columns"]) > 8
    root = document(build_docx(res, now=WHEN))
    assert root.find(f"{W}body/{W}sectPr/{W}pgSz").get(W + "orient") == "landscape"
    words = text(root)
    assert "Boundaries" in words and "Stopping Probabilities" in words


def test_charts_are_embedded_and_fit_two_to_a_page():
    res = core.run("corr.one", "power", {"n": 50, "r1": 0.3}, power_curve=False)
    images = [png(760, 420), png(400, 400)]
    data = build_docx(res, images=images, now=WHEN)
    p = parts(data)
    assert p["word/media/chart1.png"] == images[0] and p["word/media/chart2.png"] == images[1]
    rels = p["word/_rels/document.xml.rels"].decode()
    assert 'Id="rIdImage1"' in rels and 'Target="media/chart2.png"' in rels
    root = document(data)
    assert "Charts" in text(root)
    WP = "{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}"
    extents = [(int(e.get("cx")), int(e.get("cy"))) for e in root.iter(WP + "extent")]
    assert len(extents) == 2
    page = 15840 - 2 * 1440  # portrait Letter's text height, in twips
    for (cx, cy), (w, h) in zip(extents, [(760, 420), (400, 400)]):
        assert abs(cx / cy - w / h) < 0.01
        assert 2 * cy / 635 < page - 1000
    # without charts in the report options, the images are left out
    hidden = parts(build_docx(res, ReportOptions(show_plots=False), images=images, now=WHEN))
    assert not any(n.startswith("word/media/") for n in hidden)


def test_errors_and_awkward_text_still_make_a_document():
    res = core.run("means.two.t", "n", {"power": 1.5})
    assert not res["ok"]
    words = text(document(build_docx(res, now=WHEN)))
    assert "Error:" in words and "power" in words
    odd = {"ok": True, "name": "A <b> & \x01 \"test\"", "solve_label": "N", "columns": [{"key": "n", "label": "N"}],
           "rows": [{"n": 3, "_statement": "Fine & <dandy>"}, {"_error": "bad \x02 row"}]}
    words = text(document(build_docx(odd, now=WHEN)))
    assert 'A <b> &  "test"' in words and "Fine & <dandy>" in words and "bad  row" in words


def test_png_size_and_write(tmp_path):
    assert png_size(png(13, 7)) == (13, 7)
    with pytest.raises(ValueError):
        png_size(b"GIF89a")
    res = core.run("prop.one", "n", {"p0": 0.5, "p1": 0.7}, power_curve=False)
    path = tmp_path / "report.docx"
    write_docx(str(path), res)
    assert zipfile.is_zipfile(path)


def test_the_command_line_writes_word(tmp_path, capsys):
    path = tmp_path / "gs.docx"
    assert cli.main(["run", "gs.means", "--example", "--docx", str(path)]) == 0
    assert f"Wrote {path}" in capsys.readouterr().err
    assert "Boundaries" in text(document(path.read_bytes()))
    assert cli.main(["run", "gs.means", "--example", "--docx", str(tmp_path / "no" / "such" / "dir.docx")]) == 1
    assert "cannot write" in capsys.readouterr().err
