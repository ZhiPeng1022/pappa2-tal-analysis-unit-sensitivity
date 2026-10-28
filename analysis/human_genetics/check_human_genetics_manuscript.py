"""Consistency checks for the Human Genetics working manuscript."""

from __future__ import annotations

import re
from pathlib import Path


MANUSCRIPT = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji\outputs"
    r"\PAPPA2_manuscript_Human_Genetics_Research_Article_v1.md"
)

EXPECTED_CITATIONS = [
    "Bikbov et al. 2020",
    "Hoste et al. 2018",
    "Ferenbach and Bonventre 2015",
    "Bonventre and Yang 2011",
    "Lake et al. 2023",
    "Knepper et al. 1999",
    "Brezis et al. 1984",
    "Oxvig 2015",
    "Dauber et al. 2016",
    "Conover et al. 2011",
    "Asghari et al. 2025",
    "Kennedy et al. 2024",
    "Nephroseq v5 2026",
    "Skinnider et al. 2021",
    "Dann et al. 2022",
    "Cippa et al. 2018",
    "Kidney Precision Medicine Project 2026a",
    "Kidney Precision Medicine Project 2026b",
    "Human Kidney Transcriptome Cohort 2026",
    "Abedini et al. 2024",
    "Human Kidney Multi-omic Atlas 2026",
    "Reck et al. 2025",
    "Human Kidney Multiome Resource 2026",
    "Han et al. 2023",
]

EXPECTED_TITLE = (
    "Analysis-unit sensitivity distinguishes donor-level expression "
    "enrichment from within-donor coupling in human kidney thick ascending "
    "limb states"
)


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    errors: list[str] = []
    warnings: list[str] = []

    normalized_text = " ".join(text.split())
    if EXPECTED_TITLE not in normalized_text:
        errors.append("Expected title not found")

    abstract_match = re.search(
        r"## Abstract\s+(.*?)\s+## Introduction",
        text,
        flags=re.DOTALL,
    )
    if not abstract_match:
        errors.append("Abstract section not found")
        abstract = ""
    else:
        abstract = abstract_match.group(1)
        abstract_without_keywords = abstract.split("**Keywords:**", 1)[0]
        cleaned = re.sub(r"\*\*[^*]+\*\*", " ", abstract_without_keywords)
        abstract_words = re.findall(r"\b[\w'-]+\b", cleaned)
        print(f"abstract_words={len(abstract_words)}")
        if not 150 <= len(abstract_words) <= 250:
            errors.append("Abstract is outside 150-250 words")

        keywords_match = re.search(
            r"\*\*Keywords:\*\*\s*(.+?)(?:\n\s*\n|\Z)",
            abstract,
            flags=re.DOTALL,
        )
        if not keywords_match:
            errors.append("Keywords line not found")
        else:
            keywords = [
                keyword.strip()
                for keyword in keywords_match.group(1).split(",")
                if keyword.strip()
            ]
            print(f"keywords={len(keywords)}")
            if not 4 <= len(keywords) <= 6:
                errors.append("Keyword count is outside 4-6")

    body = text.split("## References", 1)[0]
    numeric = re.findall(r"\[[0-9][0-9,\-]*\]", body)
    if numeric:
        errors.append(f"Numeric citations remain: {sorted(set(numeric))}")

    missing_citations = [
        citation
        for citation in EXPECTED_CITATIONS
        if citation not in body
    ]
    if missing_citations:
        errors.append(
            "Reference entries without in-text citation: "
            + ", ".join(missing_citations)
        )

    if "Additional file" in text or "Additional files" in text:
        errors.append("BMC Additional file terminology remains")
    if "10.5281/zenodo.23733571" not in text:
        errors.append("Final exact-version Zenodo DOI is missing")
    if "23698349" in text:
        errors.append("Old working Zenodo DOI remains")

    for number in range(1, 21):
        label = f"Supplementary Table S{number}"
        if text.count(label) < 1:
            errors.append(f"Missing {label}")

    heading_levels = [
        len(match.group(1))
        for match in re.finditer(r"^(#{1,6})\s+", text, flags=re.MULTILINE)
    ]
    if heading_levels and max(heading_levels) > 3:
        errors.append("More than three displayed heading levels are used")

    for figure in range(1, 9):
        if f"Fig. {figure}" not in text and f"Fig. {figure}." not in text:
            warnings.append(f"Fig. {figure} not found in text")

    print(f"characters={len(text)}")
    print(f"lines={text.count(chr(10)) + 1}")
    print(f"max_heading_level={max(heading_levels) if heading_levels else 0}")

    if warnings:
        print("WARNINGS")
        for warning in warnings:
            print(f"- {warning}")
    if errors:
        print("ERRORS")
        for error in errors:
            print(f"- {error}")
        raise SystemExit(1)

    print("All checks passed")


if __name__ == "__main__":
    main()
