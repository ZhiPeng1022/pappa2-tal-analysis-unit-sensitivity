"""Render the v3 manuscript Markdown to a submission-ready PDF."""

import os
import re
from pathlib import Path

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    Paragraph,
    Preformatted,
    SimpleDocTemplate,
    Spacer,
)


ROOT = Path(
    os.environ.get(
        "PAPPA2_PROJECT_ROOT",
        Path(__file__).resolve().parents[2],
    )
)
SOURCE = ROOT / "submission" / "manuscript" / "PAPPA2_manuscript_draft_v3_TAL_atlas.md"
OUTPUT = ROOT / "submission" / "manuscript" / "PAPPA2_manuscript_v3.pdf"
STOP_MARKER = "## Chinese drafting notes"


def escape(text):
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def markup(text):
    text = escape(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"`(.+?)`", r"<font name='Courier'>\1</font>", text)
    return text


def build():
    styles = getSampleStyleSheet()
    body = ParagraphStyle(
        "body", parent=styles["BodyText"], fontSize=10, leading=14, spaceAfter=6
    )
    h1 = ParagraphStyle(
        "h1", parent=styles["Heading1"], fontSize=16, leading=20, spaceAfter=10
    )
    h2 = ParagraphStyle(
        "h2", parent=styles["Heading2"], fontSize=13, leading=16, spaceAfter=8
    )
    h3 = ParagraphStyle(
        "h3", parent=styles["Heading3"], fontSize=11, leading=14, spaceAfter=6
    )
    bullet = ParagraphStyle(
        "bullet", parent=body, leftIndent=14, bulletIndent=4, spaceAfter=3
    )
    code = ParagraphStyle("code", parent=styles["Code"], fontSize=8, leading=10)

    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    story = []
    in_code = False
    code_buffer = []

    for line in lines:
        if line.strip().startswith(STOP_MARKER):
            break
        if line.strip().startswith("```"):
            if in_code:
                story.append(Preformatted("\n".join(code_buffer), code))
                story.append(Spacer(1, 4))
                code_buffer = []
            in_code = not in_code
            continue
        if in_code:
            code_buffer.append(line)
            continue
        if line.startswith("# "):
            story.append(Paragraph(markup(line[2:].strip()), h1))
        elif line.startswith("## "):
            story.append(Paragraph(markup(line[3:].strip()), h2))
        elif line.startswith("### "):
            story.append(Paragraph(markup(line[4:].strip()), h3))
        elif line.strip().startswith("- "):
            story.append(
                Paragraph(markup(line.strip()[2:]), bullet, bulletText="\u2022")
            )
        elif re.match(r"^\d+\.\s", line.strip()):
            number = re.match(r"^(\d+)\.\s", line.strip()).group(1)
            story.append(
                Paragraph(
                    markup(re.sub(r"^\d+\.\s", "", line.strip())),
                    bullet,
                    bulletText=number + ".",
                )
            )
        elif line.strip():
            story.append(Paragraph(markup(line.strip()), body))

    document = SimpleDocTemplate(
        str(OUTPUT),
        pagesize=A4,
        leftMargin=2 * cm,
        rightMargin=2 * cm,
        topMargin=2 * cm,
        bottomMargin=2 * cm,
        title="PAPPA2 TAL neighborhood analysis",
    )
    document.build(story)
    print(OUTPUT)
    print("size:", OUTPUT.stat().st_size)


if __name__ == "__main__":
    build()
