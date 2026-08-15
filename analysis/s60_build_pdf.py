"""Assemble the report, the figures and the figure legends into one PDF.

Reads results/REPORT.md and results/FIGURE_LEGENDS.md, renders them to A4, and
places each figure on its own page with its legend underneath. Output goes to
results/TNFRSF9_mast_cells_AD_report.pdf.

Text is set in Liberation Sans, which is metric-compatible with the Nimbus Sans
used in the figures. The two superscript signs used in the text are not in that
font's glyph table, so they are emitted as reportlab markup instead.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether,
                                NextPageTemplate, PageBreak, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = RES / "figures"
OUT = RES / "TNFRSF9_mast_cells_AD_report.pdf"

FONT_DIR = Path("/usr/share/fonts/truetype/liberation")
PAGE_W, PAGE_H = A4
MARGIN = 18 * mm
BODY_W = PAGE_W - 2 * MARGIN

FIGURES = [
    ("Figure 1", "fig1_umap_atlas.png"),
    ("Figure 2", "fig2_quantification.png"),
    ("Figure 3", "fig3_spatial.png"),
    ("Figure 4", "fig4_bulk.png"),
    ("Figure 5", "fig5_replication.png"),
    ("Figure 6", "fig6_invitro.png"),
    ("Supplementary Figure 1", "figS1_controls.png"),
]


# ------------------------------------------------------------------ fonts ---
def register_fonts() -> None:
    for name, fn in (("LSans", "LiberationSans-Regular.ttf"),
                     ("LSans-Bold", "LiberationSans-Bold.ttf"),
                     ("LSans-Italic", "LiberationSans-Italic.ttf"),
                     ("LMono", "LiberationMono-Regular.ttf")):
        pdfmetrics.registerFont(TTFont(name, str(FONT_DIR / fn)))
    pdfmetrics.registerFontFamily("LSans", normal="LSans", bold="LSans-Bold",
                                  italic="LSans-Italic", boldItalic="LSans-Bold")


# ----------------------------------------------------------------- styles ---
def styles() -> dict:
    ss = getSampleStyleSheet()
    base = dict(fontName="LSans", leading=11.6, spaceAfter=5, textColor=colors.HexColor("#1A1A1A"))
    return {
        "title": ParagraphStyle("title", ss["Normal"], fontName="LSans-Bold",
                                fontSize=15, leading=18.5, spaceAfter=5,
                                textColor=colors.HexColor("#1A1A1A")),
        "subtitle": ParagraphStyle("subtitle", ss["Normal"], fontName="LSans",
                                   fontSize=10.2, leading=13.5, spaceAfter=12,
                                   textColor=colors.HexColor("#55555F")),
        "h2": ParagraphStyle("h2", ss["Normal"], fontName="LSans-Bold", fontSize=11,
                             leading=13.5, spaceBefore=13, spaceAfter=5,
                             textColor=colors.HexColor("#1A1A1A")),
        "h3": ParagraphStyle("h3", ss["Normal"], fontName="LSans-Bold", fontSize=9.4,
                             leading=12, spaceBefore=9, spaceAfter=3.5,
                             textColor=colors.HexColor("#33333A")),
        "body": ParagraphStyle("body", ss["Normal"], fontSize=8.7, alignment=TA_JUSTIFY, **base),
        "bullet": ParagraphStyle("bullet", ss["Normal"], fontSize=8.7, leftIndent=11,
                                 bulletIndent=3, alignment=TA_LEFT, **base),
        "numbered": ParagraphStyle("numbered", ss["Normal"], fontSize=8.7, leftIndent=14,
                                   bulletIndent=3, alignment=TA_LEFT, **base),
        "cell": ParagraphStyle("cell", ss["Normal"], fontName="LSans", fontSize=7.3,
                               leading=9, textColor=colors.HexColor("#1A1A1A")),
        "cellh": ParagraphStyle("cellh", ss["Normal"], fontName="LSans-Bold", fontSize=7.3,
                                leading=9, textColor=colors.HexColor("#1A1A1A")),
        "legend": ParagraphStyle("legend", ss["Normal"], fontSize=8.0, leading=10.4,
                                 alignment=TA_JUSTIFY, spaceAfter=4,
                                 fontName="LSans", textColor=colors.HexColor("#1A1A1A")),
        "legtitle": ParagraphStyle("legtitle", ss["Normal"], fontName="LSans-Bold",
                                   fontSize=9.2, leading=11.5, spaceBefore=7, spaceAfter=4,
                                   textColor=colors.HexColor("#1A1A1A")),
        "foot": ParagraphStyle("foot", ss["Normal"], fontName="LSans", fontSize=7.2,
                               leading=9, textColor=colors.HexColor("#77777F")),
    }


# ------------------------------------------------------------- md -> flow ---
SUPERS = {"⁺": "<super>+</super>", "⁻": "<super>-</super>"}


def inline(s: str) -> str:
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    s = re.sub(r"`([^`]+)`", r'<font name="LMono" size="7.9">\1</font>', s)
    s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"<i>\1</i>", s)
    for k, v in SUPERS.items():
        s = s.replace(k, v)
    return s


def make_table(rows: list[list[str]], st: dict) -> Table:
    header, body = rows[0], rows[1:]
    ncol = len(header)
    data = [[Paragraph(inline(c), st["cellh"]) for c in header]]
    data += [[Paragraph(inline(c), st["cell"]) for c in r] for r in body]
    # first column carries the labels and needs more room
    first = min(0.34, max(0.18, 1.6 / ncol))
    widths = [BODY_W * first] + [BODY_W * (1 - first) / (ncol - 1)] * (ncol - 1) \
        if ncol > 1 else [BODY_W]
    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.4),
        ("LEFTPADDING", (0, 0), (-1, -1), 3.5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3.5),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.HexColor("#3A3A44")),
        ("LINEBELOW", (0, 1), (-1, -2), 0.25, colors.HexColor("#D8D8DE")),
        ("LINEBELOW", (0, -1), (-1, -1), 0.6, colors.HexColor("#3A3A44")),
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#F2F2F5")),
    ]))
    return t


def md_to_flow(md: str, st: dict, *, skip_h1: bool = False) -> list:
    flow, i = [], 0
    lines = md.split("\n")
    while i < len(lines):
        ln = lines[i]
        s = ln.strip()

        if not s or s == "---":
            i += 1
            continue

        if s.startswith("|") and i + 1 < len(lines) and re.match(r"^\|[\s:|-]+\|$", lines[i + 1].strip()):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not re.match(r"^[\s:|-]+$", "|".join(cells)):
                    rows.append(cells)
                i += 1
            if rows:
                flow += [Spacer(1, 2.5), make_table(rows, st), Spacer(1, 6)]
            continue

        if s.startswith("#"):
            level = len(s) - len(s.lstrip("#"))
            text = s.lstrip("#").strip()
            if level == 1:
                if not skip_h1:
                    flow.append(Paragraph(inline(text), st["title"]))
            elif level == 2:
                flow.append(Paragraph(inline(text), st["h2"]))
            else:
                flow.append(Paragraph(inline(text), st["h3"]))
            i += 1
            continue

        if re.match(r"^([-*]|\d+\.)\s+", s):
            ordered = bool(re.match(r"^\d+\.\s+", s))
            items = []
            while i < len(lines) and re.match(r"^\s*([-*]|\d+\.)\s+", lines[i]):
                marker = re.match(r"^\s*(\S+)\s+", lines[i]).group(1)
                txt = re.sub(r"^\s*([-*]|\d+\.)\s+", "", lines[i])
                i += 1
                # a wrapped line is any indented continuation, whatever the
                # indent width the source happens to use (2 for "*", 3 for "1.")
                while i < len(lines) and lines[i].strip() \
                        and len(lines[i]) - len(lines[i].lstrip()) >= 2 \
                        and not re.match(r"^\s*([-*]|\d+\.)\s+", lines[i]):
                    txt += " " + lines[i].strip()
                    i += 1
                items.append((marker, txt))
            style = st["numbered"] if ordered else st["bullet"]
            for marker, it in items:
                flow.append(Paragraph(inline(it), style,
                                      bulletText=marker if ordered else "•"))
            flow.append(Spacer(1, 3))
            continue

        para = [s]
        i += 1
        while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith(("#", "|", "-", "*", ">")) \
                and lines[i].strip() != "---":
            para.append(lines[i].strip())
            i += 1
        flow.append(Paragraph(inline(" ".join(para)), st["body"]))
    return flow


# --------------------------------------------------------------- legends ---
def parse_legends(md: str) -> tuple[str, dict]:
    preamble, blocks = [], {}
    parts = re.split(r"^## ", md, flags=re.M)
    for chunk in parts[1:]:
        head, _, rest = chunk.partition("\n")
        key = head.split("|")[0].strip()
        blocks[key] = (head.strip(), rest.strip())
    intro = re.split(r"^## ", md, flags=re.M)[0]
    intro = re.sub(r"^# .*$", "", intro, flags=re.M)
    intro = re.sub(r"^\s*-{3,}\s*$", "", intro, flags=re.M).strip()
    return intro, blocks


def main() -> int:
    register_fonts()
    st = styles()

    report = (RES / "REPORT.md").read_text(encoding="utf-8")
    legends_md = (RES / "FIGURE_LEGENDS.md").read_text(encoding="utf-8")
    intro, legends = parse_legends(legends_md)

    # strip the repository navigation lines from the report header
    # drop the repository navigation paragraph and the section-numbering note
    report = re.sub(r"^Figures `results/figures/`.*?(?=\n\s*\n)", "", report,
                    flags=re.M | re.S)
    report = re.sub(r"^\*Sections below are numbered.*?\*\s*$", "", report,
                    flags=re.M | re.S)

    title_line = report.split("\n")[0].lstrip("# ").strip()
    sub_line = ""
    m = re.search(r"^### (.+)$", report, flags=re.M)
    if m:
        sub_line = m.group(1).strip()
        report = report.replace(m.group(0), "", 1)
    report = "\n".join(report.split("\n")[1:])

    flow: list = [
        Paragraph(inline(title_line), st["title"]),
        Paragraph(inline(sub_line), st["subtitle"]) if sub_line else Spacer(1, 2),
        Paragraph("Analysis of six public human transcriptomic cohorts. "
                  "All figures, legends and statistics are reproducible from "
                  '<font name="LMono" size="7.9">analysis/</font> in the project repository.',
                  st["foot"]),
        Spacer(1, 9),
    ]
    flow += md_to_flow(report, st, skip_h1=True)

    # ---- figures, one per page, legend beneath ---------------------------
    # no page break here: the conventions paragraph runs on from the report so
    # the figures are not preceded by a nearly empty separator page
    flow.append(KeepTogether([Paragraph("Figures", st["h2"]),
                              Paragraph(inline(intro), st["legend"])]))

    for name, fname in FIGURES:
        path = FIG / fname
        if not path.exists():
            continue
        head, body = legends.get(name, (name, ""))
        from PIL import Image as PILImage
        with PILImage.open(path) as im:
            w, h = im.size
        draw_w = BODY_W
        draw_h = draw_w * h / w
        max_h = PAGE_H - 2 * MARGIN - 82 * mm      # leave room for the legend
        if draw_h > max_h:
            draw_h = max_h
            draw_w = draw_h * w / h
        flow.append(PageBreak())
        flow.append(Image(str(path), width=draw_w, height=draw_h))
        flow.append(Spacer(1, 7))
        flow.append(Paragraph(inline(head.replace("|", "│")), st["legtitle"]))
        flow += md_to_flow(body, st | {"body": st["legend"]}, skip_h1=True)

    doc = BaseDocTemplate(str(OUT), pagesize=A4,
                          leftMargin=MARGIN, rightMargin=MARGIN,
                          topMargin=MARGIN, bottomMargin=MARGIN,
                          title="TNFRSF9 in cutaneous mast cells in atopic dermatitis",
                          author="MastCells_TNFRS9")

    def footer(canv, _doc):
        canv.saveState()
        canv.setFont("LSans", 7)
        canv.setFillColor(colors.HexColor("#8A8A92"))
        canv.drawString(MARGIN, 11 * mm, "TNFRSF9 in cutaneous mast cells in atopic dermatitis")
        canv.drawRightString(PAGE_W - MARGIN, 11 * mm, str(canv.getPageNumber()))
        canv.restoreState()

    frame = Frame(MARGIN, MARGIN, BODY_W, PAGE_H - 2 * MARGIN, id="body",
                  leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
    doc.addPageTemplates([PageTemplate(id="all", frames=[frame], onPage=footer)])
    doc.build(flow)
    print(f"wrote {OUT.relative_to(ROOT)}  ({OUT.stat().st_size/1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
