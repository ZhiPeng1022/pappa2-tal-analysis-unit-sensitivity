"""Programmatic layout checks for the Human Genetics PDF preview."""

from __future__ import annotations

import re
from pathlib import Path

import pymupdf


PDF = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji\outputs"
    r"\PAPPA2_manuscript_Human_Genetics_Research_Article_v1.pdf"
)

EXPECTED_TITLE = (
    "Analysis-unit sensitivity distinguishes donor-level expression "
    "enrichment from within-donor coupling in human kidney thick ascending "
    "limb states"
)


def main() -> None:
    document = pymupdf.open(PDF)
    errors: list[str] = []
    blank_pages: list[int] = []
    image_count = 0
    min_width = 595.0
    max_width = 596.0

    for index, page in enumerate(document, start=1):
        text = page.get_text("text").strip()
        if not text:
            blank_pages.append(index)
        if not min_width <= page.rect.width <= max_width:
            errors.append(
                f"Page {index} width is {page.rect.width:.2f} pt, not A4"
            )
        images = page.get_images(full=True)
        image_count += len(images)

        for block in page.get_text("blocks"):
            x0, y0, x1, y1 = block[:4]
            if x0 < 35 or x1 > page.rect.width - 35:
                errors.append(
                    f"Text block outside horizontal margin on page {index}: "
                    f"{x0:.1f}-{x1:.1f}"
                )

    first_text = document[0].get_text("text")
    first_text_without_lines = re.sub(
        r"^\s*\d+\s*$",
        "",
        first_text,
        flags=re.MULTILINE,
    )
    normalized_first = " ".join(first_text_without_lines.split())
    normalized_title = " ".join(EXPECTED_TITLE.split())
    if normalized_title not in normalized_first:
        errors.append("Expected title is not on page 1")
    if blank_pages:
        errors.append(f"Blank pages: {blank_pages}")
    if image_count < 8:
        errors.append(f"Only {image_count} embedded images found")

    print(f"pages={document.page_count}")
    print(f"page_width_pt={document[0].rect.width:.2f}")
    print(f"page_height_pt={document[0].rect.height:.2f}")
    print(f"embedded_images={image_count}")
    print(f"first_page_characters={len(first_text)}")
    print(f"last_page_characters={len(document[-1].get_text('text'))}")

    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        raise SystemExit(1)

    print("All PDF checks passed")


if __name__ == "__main__":
    main()
