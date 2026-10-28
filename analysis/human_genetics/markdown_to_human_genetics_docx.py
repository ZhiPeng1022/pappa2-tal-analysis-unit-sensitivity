"""Export the Human Genetics Research Article Markdown to DOCX."""

from __future__ import annotations

import re
from pathlib import Path

from docx import Document
from docx.shared import Cm, Pt

import markdown_to_bmc_docx as base


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
SOURCE = (
    ROOT
    / "outputs"
    / "PAPPA2_manuscript_Human_Genetics_Research_Article_v1.md"
)
OUTPUT = (
    ROOT
    / "outputs"
    / "PAPPA2_manuscript_Human_Genetics_Research_Article_v1.docx"
)
FIGURE_DIR = ROOT / "outputs" / "figures"

EXTRA_GENES = [
    "ARHGEF3",
    "C4orf19",
    "KCNIP4",
    "MAST4",
    "NEAT1",
    "PPP2R3A",
    "SOD2",
    "SPP1",
    "TBC1D1",
    "TFPI",
    "TNS3",
    "HIF1A",
    "ADAMTS1",
    "NRP1",
    "UMOD",
    "SLC12A1",
    "KCNJ1",
    "JAG1",
    "STAT6",
]


def main() -> None:
    base.SOURCE = str(SOURCE)
    base.OUTPUT = str(OUTPUT)
    base.FIGURE_DIR = str(FIGURE_DIR)
    base.BODY_SIZE = Pt(10)
    base.GENE_SYMBOLS = sorted(
        set(base.GENE_SYMBOLS + EXTRA_GENES),
        key=len,
        reverse=True,
    )
    base.GENE_RE = re.compile(
        r"\b(" + "|".join(base.GENE_SYMBOLS) + r")\b"
    )

    base.convert(str(SOURCE), str(OUTPUT))

    document = Document(str(OUTPUT))
    for section in document.sections:
        section.page_width = Cm(21.0)
        section.page_height = Cm(29.7)
        section.left_margin = Cm(2.5)
        section.right_margin = Cm(2.5)
        section.top_margin = Cm(2.5)
        section.bottom_margin = Cm(2.5)
    document.save(str(OUTPUT))

    print(f"Written: {OUTPUT}")
    print(f"Bytes: {OUTPUT.stat().st_size}")


if __name__ == "__main__":
    main()
