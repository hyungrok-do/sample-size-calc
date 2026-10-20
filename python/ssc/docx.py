"""The report as a Word document (.docx), written with the standard library only.

A .docx file is a zip of XML parts (Office Open XML).  ``build_docx`` writes the report of a result (``ssc-calc run
--docx``): the method, what was solved for, the results table, the plain-language statements, dropout, the extra
tables (such as group-sequential boundaries), the charts (PNG images the caller renders), definitions, notes and
references.  Word, LibreOffice, Google Docs and Pages open it, ready to be
edited into a protocol or a grant.  No library is needed (python-docx would bring lxml, several megabytes).
"""
from __future__ import annotations

import datetime as _dt
import io
import re
import struct
import zipfile
from typing import Any, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape

from .formatting import ReportOptions, fmt_cell

INK = "1F1A33"
ACCENT = "4B3DC4"
MUTED = "6E6888"
ERROR = "D03B3B"
HEAD_BG = "F0EEFB"
ROW_ALT = "F9F8FD"
RULE = "E4E1EE"

EMU_PER_TWIP = 635            # English metric units (images) per twentieth of a point (page layout)
WIDE_TABLE = 8                # tables with more columns than this turn the page sideways
TEXT_COLUMN = 2200            # a text column's width, in twips (about 1.5 inches), as the HTML report's

_W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
_R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"
_A = "http://schemas.openxmlformats.org/drawingml/2006/main"
_PIC = "http://schemas.openxmlformats.org/drawingml/2006/picture"
_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"
_NOT_XML = re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f￾￿]")   # characters XML 1.0 cannot carry


def _text(value: Any) -> str:
    return escape(_NOT_XML.sub("", str(value)))


# ----------------------------------------------------------------------------- paragraphs and runs
def run(text: Any, *, bold: bool = False, italic: bool = False, color: str = "", size: float = 0) -> str:
    props = ("<w:b/>" if bold else "") + ("<w:i/>" if italic else "") + \
        (f'<w:color w:val="{color}"/>' if color else "") + \
        (f'<w:sz w:val="{round(size * 2)}"/><w:szCs w:val="{round(size * 2)}"/>' if size else "")
    return (f"<w:r>{f'<w:rPr>{props}</w:rPr>' if props else ''}"
            f'<w:t xml:space="preserve">{_text(text)}</w:t></w:r>')


def para(*runs: str, style: str = "", align: str = "", after: Optional[int] = None, hanging: int = 0,
         keep_next: bool = False) -> str:
    props = (f'<w:pStyle w:val="{style}"/>' if style else "") + ("<w:keepNext/>" if keep_next else "") + \
        (f'<w:spacing w:before="0" w:after="{after}"/>' if after is not None else "") + \
        (f'<w:ind w:left="{hanging}" w:hanging="{hanging}"/>' if hanging else "") + \
        (f'<w:jc w:val="{align}"/>' if align else "")
    return f"<w:p>{f'<w:pPr>{props}</w:pPr>' if props else ''}{''.join(runs)}</w:p>"


def heading(text: str) -> str:
    return para(run(text), style="Heading1")


# ----------------------------------------------------------------------------- tables
def _cell(content: str, *, fill: str = "", span: int = 1, align: str = "", width: int = 0) -> str:
    props = (f'<w:tcW w:w="{width}" w:type="dxa"/>' if width else "") + \
        (f'<w:gridSpan w:val="{span}"/>' if span > 1 else "") + \
        (f'<w:shd w:val="clear" w:color="auto" w:fill="{fill}"/>' if fill else "")
    return (f"<w:tc>{f'<w:tcPr>{props}</w:tcPr>' if props else ''}"
            f"{para(content, align=align, after=0)}</w:tc>")


def table(columns: Sequence[dict], rows: Sequence[dict], opts: ReportOptions, width: int) -> str:
    """A results table: a shaded header row (repeated on each page), numbers on the right, text on the left,
    rows banded, and a scenario that could not be computed as one cell with its message."""
    n = max(1, len(columns))
    grid = "".join(f'<w:gridCol w:w="{width // n}"/>' for _ in range(n))
    out = [f'<w:tbl><w:tblPr><w:tblStyle w:val="ReportTable"/><w:tblW w:w="0" w:type="auto"/></w:tblPr>'
           f"<w:tblGrid>{grid}</w:tblGrid>", '<w:tr><w:trPr><w:tblHeader/></w:trPr>']
    for c in columns:
        text = c.get("fmt") == "text"
        out.append(_cell(run(c["label"], bold=True), fill=HEAD_BG, align="left" if text else "right",
                         width=TEXT_COLUMN if text else 0))
    out.append("</w:tr>")
    for i, r in enumerate(rows):
        fill = ROW_ALT if i % 2 else ""
        out.append("<w:tr>")
        if r.get("_error"):
            out.append(_cell(run(f"Not computed — {r['_error']}", color=ERROR), fill=fill, span=n, align="left"))
        else:
            for c in columns:
                text = c.get("fmt") == "text"
                out.append(_cell(run(fmt_cell(r.get(c["key"]), c.get("fmt", "num"), opts)), fill=fill,
                                 align="left" if text else "right", width=TEXT_COLUMN if text else 0))
        out.append("</w:tr>")
    out.append("</w:tbl>")
    out.append(para(after=60))  # Word needs a paragraph between a table and what follows
    return "".join(out)


def _definitions(items: Sequence[Tuple[str, str]], width: int) -> str:
    first, second = width // 4, width - width // 4
    out = [f'<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/></w:tblPr>'
           f'<w:tblGrid><w:gridCol w:w="{first}"/><w:gridCol w:w="{second}"/></w:tblGrid>']
    for sym, meaning in items:
        out.append(f"<w:tr>{_cell(run(sym, bold=True))}{_cell(run(meaning))}</w:tr>")
    out.append("</w:tbl>")
    out.append(para(after=60))
    return "".join(out)


# ----------------------------------------------------------------------------- images
def png_size(data: bytes) -> Tuple[int, int]:
    """Width and height of a PNG image, in pixels."""
    if data[:8] != b"\x89PNG\r\n\x1a\n" or data[12:16] != b"IHDR":
        raise ValueError("not a PNG image")
    return struct.unpack(">II", data[16:24])


def _picture(index: int, rel: str, data: bytes, max_w: int, max_h: int) -> str:
    """A chart as an inline picture, as large as fits in max_w × max_h twips."""
    w, h = png_size(data)
    scale = min(max_w / max(1, w), max_h / max(1, h))
    cx, cy = round(w * scale) * EMU_PER_TWIP, round(h * scale) * EMU_PER_TWIP
    return (f'<w:p><w:pPr><w:spacing w:before="60" w:after="120"/></w:pPr><w:r><w:drawing>'
            f'<wp:inline distT="0" distB="0" distL="0" distR="0"><wp:extent cx="{cx}" cy="{cy}"/>'
            f'<wp:docPr id="{index}" name="Chart {index}"/>'
            f'<wp:cNvGraphicFramePr><a:graphicFrameLocks noChangeAspect="1"/></wp:cNvGraphicFramePr>'
            f'<a:graphic><a:graphicData uri="{_PIC}"><pic:pic>'
            f'<pic:nvPicPr><pic:cNvPr id="{index}" name="chart{index}.png"/><pic:cNvPicPr/></pic:nvPicPr>'
            f'<pic:blipFill><a:blip r:embed="{rel}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>'
            f'<pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{cx}" cy="{cy}"/></a:xfrm>'
            f'<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr>'
            f"</pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>")


# ----------------------------------------------------------------------------- the report
def _body(result: dict, opts: ReportOptions, images: Sequence[bytes], width: int, height: int,
          stamp: str) -> Tuple[str, int]:
    """The document's paragraphs and tables, and how many images it shows."""
    parts: List[str] = [para(run(result.get("name", "Results")), style="Title")]
    if not result.get("ok", False):
        parts.append(para(run("Error: ", bold=True, color=ERROR),
                          run(result.get("error", "Calculation failed."), color=ERROR)))
        for key, msg in (result.get("input_errors") or {}).items():
            parts.append(para(run(f"•  {key}: ", color=ERROR), run(msg), hanging=280))
        return "".join(parts), 0
    rows = result.get("rows", [])
    n = len(rows)
    threads = result.get("threads", 1)
    parts.append(para(run("Solved for ", color=MUTED), run(result.get("solve_label", ""), bold=True, color=MUTED),
                      run(f"  ·  {n} scenario{'s' if n != 1 else ''}  ·  {result.get('elapsed_ms', 0.0):.1f} ms "
                          f"on {threads} thread{'s' if threads != 1 else ''}", color=MUTED)))
    parts.append(heading("Results"))
    parts.append(table(result.get("columns", []), rows, opts, width))
    bad = sum(1 for r in rows if r.get("_error"))
    if bad:
        parts.append(para(run(f"{bad} scenario{'s' if bad != 1 else ''} could not be computed (see the messages "
                              "in the table).", color=ERROR)))
    if opts.statements != "none":
        stmts = [(i, r["_statement"]) for i, r in enumerate(rows) if r.get("_statement")]
        if opts.statements == "first":
            stmts = stmts[:1]
        if stmts:
            parts.append(heading("In Plain Words"))
            for i, s in stmts:
                parts.append(para(*([run(f"Scenario {i + 1}. ", color=MUTED)] if n > 1 else []), run(s)))
    for i, notes in [(i, r["_notes"]) for i, r in enumerate(rows) if r.get("_notes")][:20]:
        for note in notes:
            prefix = f"Scenario {i + 1}: " if len(rows) > 1 else ""
            parts.append(para(run(f"{prefix}{note}", color=MUTED)))
    dropout = result.get("dropout")
    if opts.show_dropout and dropout and dropout.get("rows"):
        parts.append(heading("Allowing for Dropout"))
        parts.append(table(dropout["columns"], dropout["rows"], opts, width))
        parts.append(para(run("Inflated sizes N′ = N / (1 − dropout rate), rounded up; D = N′ − N is the expected "
                              "number of dropouts.", color=MUTED)))
    for i, r in enumerate(rows):
        for t in r.get("_tables", []) or []:
            cols = t.get("columns", [])
            trows = [{c["key"]: v for c, v in zip(cols, row)} for row in t.get("rows", [])]
            parts.append(heading(t.get("title", "") + (f" — scenario {i + 1}" if n > 1 else "")))
            parts.append(table(cols, trows, opts, width))
    shown = 0
    if opts.show_plots and images:
        parts.append(heading("Charts"))
        for k, data in enumerate(images, start=1):  # two to a page, below the heading
            parts.append(_picture(k, f"rIdImage{k}", data, min(width, 9000), (height - 1000) // 2 - 200))
            shown += 1
    if opts.show_definitions and result.get("definitions"):
        parts.append(heading("Definitions"))
        parts.append(_definitions(result["definitions"], width))
    if opts.show_notes and result.get("notes"):
        parts.append(heading("Notes"))
        for note in result["notes"]:
            parts.append(para(run("•  "), run(note), hanging=280))
    if opts.show_references and result.get("references"):
        parts.append(heading("References"))
        for k, ref in enumerate(result["references"], start=1):
            parts.append(para(run(f"{k}.  "), run(ref), hanging=360))
    parts.append(para(run(f"Generated {stamp} by Sample Size Calculator.", color=MUTED, size=9)))
    return "".join(parts), shown


def _widest(result: dict) -> int:
    counts = [len(result.get("columns", []))]
    dropout = result.get("dropout") or {}
    counts.append(len(dropout.get("columns", [])))
    for r in result.get("rows", []):
        counts += [len(t.get("columns", [])) for t in r.get("_tables", []) or []]
    return max(counts)


STYLES = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="{_W}">
<w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:eastAsia="Calibri" w:cs="Calibri"/>
<w:color w:val="{INK}"/><w:sz w:val="21"/><w:szCs w:val="21"/><w:lang w:val="en-US"/></w:rPr></w:rPrDefault>
<w:pPrDefault><w:pPr><w:spacing w:after="120" w:line="264" w:lineRule="auto"/></w:pPr></w:pPrDefault></w:docDefaults>
<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:qFormat/></w:style>
<w:style w:type="paragraph" w:styleId="Title"><w:name w:val="Title"/><w:basedOn w:val="Normal"/>
<w:next w:val="Normal"/><w:qFormat/><w:pPr><w:spacing w:after="60"/></w:pPr>
<w:rPr><w:b/><w:sz w:val="34"/><w:szCs w:val="34"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:basedOn w:val="Normal"/>
<w:next w:val="Normal"/><w:qFormat/><w:pPr><w:keepNext/><w:spacing w:before="300" w:after="100"/>
<w:outlineLvl w:val="0"/></w:pPr><w:rPr><w:b/><w:color w:val="{ACCENT}"/><w:sz w:val="25"/><w:szCs w:val="25"/></w:rPr>
</w:style>
<w:style w:type="table" w:default="1" w:styleId="TableNormal"><w:name w:val="Normal Table"/><w:semiHidden/>
<w:tblPr><w:tblInd w:w="0" w:type="dxa"/><w:tblCellMar><w:top w:w="0" w:type="dxa"/><w:left w:w="108" w:type="dxa"/>
<w:bottom w:w="0" w:type="dxa"/><w:right w:w="108" w:type="dxa"/></w:tblCellMar></w:tblPr></w:style>
<w:style w:type="table" w:styleId="ReportTable"><w:name w:val="Report Table"/><w:basedOn w:val="TableNormal"/>
<w:rPr><w:sz w:val="18"/><w:szCs w:val="18"/></w:rPr><w:tblPr><w:tblBorders>
<w:bottom w:val="single" w:sz="4" w:space="0" w:color="{RULE}"/>
<w:insideH w:val="single" w:sz="4" w:space="0" w:color="{RULE}"/></w:tblBorders>
<w:tblCellMar><w:top w:w="60" w:type="dxa"/><w:left w:w="110" w:type="dxa"/><w:bottom w:w="60" w:type="dxa"/>
<w:right w:w="110" w:type="dxa"/></w:tblCellMar></w:tblPr></w:style>
</w:styles>"""


def build_docx(result: dict, opts: Optional[ReportOptions] = None, images: Sequence[bytes] = (),
               now: Optional[_dt.datetime] = None) -> bytes:
    """The report of a calculation result as a .docx file's bytes; ``images`` are the charts, as PNG data."""
    opts = opts or ReportOptions()
    now = now or _dt.datetime.now()
    wide = _widest(result) > WIDE_TABLE
    # US Letter, in twips: portrait with 1-inch margins, or sideways with 0.75-inch ones for wide tables
    page_w, page_h, margin = (15840, 12240, 1080) if wide else (12240, 15840, 1440)
    width = page_w - 2 * margin
    body, shown = _body(result, opts, images, width, page_h - 2 * margin, now.strftime("%Y-%m-%d %H:%M"))
    orient = ' w:orient="landscape"' if wide else ""
    document = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                f'<w:document xmlns:w="{_W}" xmlns:r="{_R}" xmlns:wp="{_WP}" xmlns:a="{_A}" xmlns:pic="{_PIC}">'
                f"<w:body>{body}<w:sectPr><w:pgSz w:w=\"{page_w}\" w:h=\"{page_h}\"{orient}/>"
                f'<w:pgMar w:top="{margin}" w:right="{margin}" w:bottom="{margin}" w:left="{margin}" '
                f'w:header="720" w:footer="720" w:gutter="0"/></w:sectPr></w:body></w:document>')
    rels = [f'<Relationship Id="rIdStyles" Type="{_REL}/styles" Target="styles.xml"/>']
    rels += [f'<Relationship Id="rIdImage{k}" Type="{_REL}/image" Target="media/chart{k}.png"/>'
             for k in range(1, shown + 1)]
    title = _text(result.get("name", "Results"))
    created = now.astimezone(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    parts = {
        "[Content_Types].xml":
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Default Extension="png" ContentType="image/png"/>'
            '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.'
            'wordprocessingml.document.main+xml"/>'
            '<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.'
            'wordprocessingml.styles+xml"/>'
            '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.'
            'core-properties+xml"/>'
            '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.'
            'extended-properties+xml"/></Types>',
        "_rels/.rels":
            f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="{_PKG_REL}">'
            f'<Relationship Id="rId1" Type="{_REL}/officeDocument" Target="word/document.xml"/>'
            f'<Relationship Id="rId2" Type="{_PKG_REL}/metadata/core-properties" Target="docProps/core.xml"/>'
            f'<Relationship Id="rId3" Type="{_REL}/extended-properties" Target="docProps/app.xml"/>'
            "</Relationships>",
        "docProps/core.xml":
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
            'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            f"<dc:title>{title}</dc:title><dc:creator>Sample Size Calculator</dc:creator>"
            f'<dcterms:created xsi:type="dcterms:W3CDTF">{created}</dcterms:created></cp:coreProperties>',
        "docProps/app.xml":
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties">'
            "<Application>Sample Size Calculator</Application></Properties>",
        "word/document.xml": document,
        "word/styles.xml": STYLES,
        "word/_rels/document.xml.rels":
            f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="{_PKG_REL}">'
            + "".join(rels) + "</Relationships>",
    }
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, text in parts.items():
            z.writestr(name, text.encode("utf-8"))
        for k, data in enumerate(images[:shown], start=1):
            z.writestr(f"word/media/chart{k}.png", data, compress_type=zipfile.ZIP_STORED)  # PNG is compressed
    return buf.getvalue()


def write_docx(path: str, result: dict, opts: Optional[ReportOptions] = None, images: Sequence[bytes] = ()) -> None:
    with open(path, "wb") as fh:
        fh.write(build_docx(result, opts, images))
