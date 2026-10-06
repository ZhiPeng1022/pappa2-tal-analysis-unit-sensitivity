import os
from pathlib import Path
import re

import pdfplumber


ROOT = Path(
    os.environ.get(
        "PAPPA2_PROJECT_ROOT",
        Path(__file__).resolve().parents[2],
    )
)
SUBMISSION = ROOT / "submission"
MANUSCRIPT = (
    SUBMISSION
    / "manuscript"
    / "PAPPA2_manuscript_draft_v9_precise.md"
)
PDF = SUBMISSION / "manuscript" / "PAPPA2_manuscript_v9.pdf"
SUPPLEMENTARY = SUBMISSION / "supplementary_v6"
FIGURE_DIR = SUBMISSION / "figures"
REPORT = ROOT / "outputs" / "PAPPA2_submission_consistency_report_v1.md"


def figure_legends(text):
    return sorted(
        {
            int(number)
            for number in re.findall(
                r"\*\*Figure\s+(\d+)\.",
                text,
            )
        }
    )


def supplementary_tables(text):
    return sorted(
        {
            int(number)
            for number in re.findall(
                r"\|\s*Table\s+S(\d+)\s*\|",
                text,
            )
        }
    )


def supplementary_figures(text):
    return sorted(
        {
            int(number)
            for number in re.findall(
                r"\|\s*Figure\s+S(\d+)\s*\|",
                text,
            )
        }
    )


def check_file_group(prefix, expected_numbers, directory):
    missing = []
    found = []
    for number in expected_numbers:
        matches = sorted(directory.glob(f"{prefix}{number}_*"))
        if len(matches) == 1:
            found.append(matches[0].name)
        else:
            missing.append((number, [path.name for path in matches]))
    return found, missing


def main():
    if not MANUSCRIPT.exists():
        raise FileNotFoundError(MANUSCRIPT)
    text = MANUSCRIPT.read_text(encoding="utf-8")
    legend_numbers = figure_legends(text)
    index_text = (
        SUPPLEMENTARY / "Supplementary_index_v6.md"
    ).read_text(encoding="utf-8")
    table_numbers = supplementary_tables(index_text)
    supplementary_figure_numbers = supplementary_figures(index_text)

    v5_figures = sorted(FIGURE_DIR.glob("*_v5.png"))
    table_files, missing_tables = check_file_group(
        "Table_S",
        table_numbers,
        SUPPLEMENTARY,
    )
    figure_s_files, missing_figure_s = check_file_group(
        "Figure_S",
        supplementary_figure_numbers,
        SUPPLEMENTARY,
    )
    old_reference_hits = []
    for path in SUBMISSION.rglob("*"):
        if path.is_file() and path.suffix.lower() in {
            ".md",
            ".csv",
            ".txt",
        }:
            file_text = path.read_text(
                encoding="utf-8",
                errors="ignore",
            )
            for pattern in [
                "evidence_matrix_v1",
                "figure_assembly_plan_v1",
            ]:
                if pattern in file_text:
                    old_reference_hits.append(
                        f"{path.relative_to(ROOT)} -> {pattern}"
                    )

    page_count = None
    pdf_ok = False
    if PDF.exists():
        with pdfplumber.open(PDF) as pdf:
            page_count = len(pdf.pages)
            pdf_text = "\n".join(
                (page.extract_text() or "")
                for page in pdf.pages
            )
        pdf_ok = (
            "References" in pdf_text
            and "analysis-unit" in pdf_text.lower()
        )

    checks = [
        (
            "Manuscript contains 8 figure legends",
            legend_numbers == list(range(1, 9)),
            f"found {legend_numbers}",
        ),
        (
            "Supplementary index contains Tables S1-S14",
            table_numbers == list(range(1, 15)),
            f"found {table_numbers}",
        ),
        (
            "Supplementary index contains Figure S1",
            supplementary_figure_numbers == [1],
            f"found {supplementary_figure_numbers}",
        ),
        (
            "All v5 figure files exist",
            len(v5_figures) == 8,
            f"found {len(v5_figures)}",
        ),
        (
            "All supplementary table files exist",
            not missing_tables,
            f"missing {missing_tables}",
        ),
        (
            "All supplementary figure files exist",
            not missing_figure_s,
            f"missing {missing_figure_s}",
        ),
        (
            "No v1 matrix/figure-plan references in submission",
            not old_reference_hits,
            f"hits {old_reference_hits}",
        ),
        (
            "v9 PDF exists and contains references",
            pdf_ok,
            f"pages={page_count}",
        ),
        (
            "v9 text uses complete 106,851-cell primary analyses",
            "106,851" in text
            and "analysis-unit sensitivity" in text
            and "logMS normalization" in text,
            "complete-data wording present",
        ),
        (
            "Authors and no-funding statement are present",
            "Peng Zhi" in text
            and "Xiaohe Yan" in text
            and "Chaofan Zhang" in text
            and "received no specific grant" in text,
            "author metadata present",
        ),
        (
            "v9 text does not use the v4 title",
            "Donor-aware analysis-unit comparison reveals"
            not in text,
            "old title absent",
        ),
    ]

    lines = [
        "# PAPPA2 submission consistency report v1",
        "",
        "Date: 2026-09-30",
        "",
        "## Summary",
        "",
        f"- Figure legends: {legend_numbers}",
        f"- Supplementary tables: S{table_numbers}",
        f"- Supplementary figures: S{supplementary_figure_numbers}",
        f"- v9 manuscript with v5 figure files: {len(v5_figures)}",
        f"- v9 PDF pages: {page_count}",
        "",
        "## Checks",
        "",
        "| Check | Status | Detail |",
        "|---|---|---|",
    ]
    for label, passed, detail in checks:
        lines.append(
            f"| {label} | {'PASS' if passed else 'FAIL'} | {detail} |"
        )

    lines.extend(
        [
            "",
            "## v5 figure files",
            "",
            "```text",
            *[path.name for path in v5_figures],
            "```",
            "",
            "## Supplementary files",
            "",
            "```text",
            *table_files,
            *figure_s_files,
            "```",
        ]
    )
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Saved:", REPORT)
    for label, passed, detail in checks:
        print(
            ("PASS" if passed else "FAIL"),
            label,
            detail,
        )


if __name__ == "__main__":
    main()
