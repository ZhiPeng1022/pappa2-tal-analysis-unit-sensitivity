"""Verify Human Genetics DOCX layout and journal-level consistency."""

from __future__ import annotations

import re
import zipfile
from pathlib import Path


DOCX = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji\outputs"
    r"\PAPPA2_manuscript_Human_Genetics_Research_Article_v1.docx"
)

EXPECTED_TITLE = (
    "Analysis-unit sensitivity distinguishes donor-level expression "
    "enrichment from within-donor coupling in human kidney thick ascending "
    "limb states"
)


def main() -> None:
    archive = zipfile.ZipFile(DOCX)
    document_xml = archive.read("word/document.xml").decode("utf-8")
    styles_xml = archive.read("word/styles.xml").decode("utf-8")
    text = " ".join(
        re.findall(r"<w:t[^>]*>(.*?)</w:t>", document_xml)
    )
    text = (
        text.replace("&amp;", "&")
        .replace("&lt;", "<")
        .replace("&gt;", ">")
    )

    checks = {
        "title_present": EXPECTED_TITLE in text,
        "line_numbers": "<w:lnNumType" in document_xml,
        "double_spacing": (
            'w:line="480"' in styles_xml
            or 'w:line="480"' in document_xml
        ),
        "no_manual_page_breaks": not re.search(
            r'<w:br[^>]*type="page"',
            document_xml,
        ),
        "a4_page_size": (
            'w:w="11906"' in document_xml
            and 'w:h="16838"' in document_xml
        ),
        "eight_figures": all(
            f"Fig. {number}" in text for number in range(1, 9)
        ),
        "supplementary_tables_s1_s23": all(
            f"Supplementary Table S{number}" in text
            for number in range(1, 24)
        ),
        "no_additional_file_terms": (
            "Additional file" not in text
            and "Additional files" not in text
        ),
        "no_numeric_citations": not re.search(
            r"\[[0-9][0-9,\-]*\]",
            text,
        ),
        "gene_italic_runs": document_xml.count("<w:i/>") >= 20,
    }

    footers = [
        name for name in archive.namelist() if name.startswith("word/footer")
    ]
    checks["page_number_footer"] = any(
        "PAGE" in archive.read(name).decode("utf-8")
        for name in footers
    )

    failed = [name for name, passed in checks.items() if not passed]
    for name, passed in checks.items():
        print(f"{name}: {passed}")
    print(f"paragraphs={len(re.findall(r'<w:p[ >]', document_xml))}")
    print(f"drawings={document_xml.count('<w:drawing>')}")
    print(f"italic_runs={document_xml.count('<w:i/>')}")

    if failed:
        raise SystemExit("Failed checks: " + ", ".join(failed))

    print("All DOCX checks passed")


if __name__ == "__main__":
    main()
