"""检查 BMC 版投稿材料与仓库派生数据是否一致。

检查内容：
1. BMC 版主稿存在，且含 Fig. 1-8 图注；
2. 正文图引用按首次出现顺序 1-8；
3. 正文附加文件引用按首次出现顺序 1-13；
4. 仓库 data/derived/supplementary 下为 13 个 Additional file，且无旧 Table_S 残留；
5. 8 张 v5 主图存在；
6. BMC 版 DOCX 含连续行号且无分页符。

用法：
    python manuscript/check_submission_consistency.py
"""

import os
import re
import sys
import zipfile
from pathlib import Path

ROOT = Path(
    os.environ.get(
        "PAPPA2_PROJECT_ROOT",
        Path(__file__).resolve().parents[2],
    )
)
SUBMISSION = ROOT / "submission"
REPO = ROOT / "repo"

MANUSCRIPT_MD = (
    SUBMISSION / "manuscript" / "PAPPA2_manuscript_BMC_Medical_Genomics_v1.md"
)
MANUSCRIPT_DOCX = (
    SUBMISSION / "manuscript" / "PAPPA2_manuscript_BMC_Medical_Genomics_v1.docx"
)
FIGURE_DIR = SUBMISSION / "figures"
REPO_SUPPLEMENTARY = REPO / "data" / "derived" / "supplementary"
REPORT = ROOT / "outputs" / "PAPPA2_bmc_submission_consistency_report.md"

N_ADDITIONAL_FILES = 13


def unique_in_order(pattern, text):
    """按出现顺序返回去重后的编号列表。"""
    seen = []
    for match in re.finditer(pattern, text):
        value = int(match.group(1))
        if value not in seen:
            seen.append(value)
    return seen


def additional_citations_in_order(text):
    """按首次出现顺序返回正文引用的附加文件编号。

    只扫描正文（Introduction 到 Conclusions），并解析
    "Additional files 2 and 3" 这类一次引用多个编号的写法。
    """
    body = text.split("## List of abbreviations")[0]
    seen = []
    for match in re.finditer(r"Additional files? ([^)]+)", body):
        for number in re.findall(r"\d+", match.group(1)):
            value = int(number)
            if value not in seen:
                seen.append(value)
    return seen


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:  # noqa: BLE001
        pass

    if not MANUSCRIPT_MD.exists():
        raise FileNotFoundError(MANUSCRIPT_MD)
    text = MANUSCRIPT_MD.read_text(encoding="utf-8")

    body = text.split("## List of abbreviations")[0]
    legend_numbers = sorted(
        int(number) for number in re.findall(r"\*\*Fig\. (\d)\. ", text)
    )
    figure_citations = unique_in_order(r"\(Fig\. (\d)", body)
    additional_citations = additional_citations_in_order(text)

    repo_files = sorted(
        path.name
        for path in REPO_SUPPLEMENTARY.glob("Additional file *")
    )
    leftover = sorted(
        path.name for path in REPO_SUPPLEMENTARY.glob("Table_S*")
    )
    removed = [
        name
        for name in ("Figure_S1_analysis_units.png",
                     "Table_S9_full_TAL_pseudobulk.csv")
        if (REPO_SUPPLEMENTARY / name).exists()
    ]

    v5_figures = sorted(FIGURE_DIR.glob("*_v5.png"))

    has_lines = False
    no_page_break = False
    if MANUSCRIPT_DOCX.exists():
        document = zipfile.ZipFile(MANUSCRIPT_DOCX).read(
            "word/document.xml"
        ).decode("utf-8")
        has_lines = "lnNumType" in document
        no_page_break = not re.search(r'<w:br[^>]*type="page"', document)

    checks = [
        ("主稿含 Fig. 1-8 图注", legend_numbers == list(range(1, 9)),
         f"found {legend_numbers}"),
        ("正文图引用顺序 1-8", figure_citations == list(range(1, 9)),
         f"found {figure_citations}"),
        (f"正文附加文件引用顺序 1-{N_ADDITIONAL_FILES}",
         additional_citations == list(range(1, N_ADDITIONAL_FILES + 1)),
         f"found {additional_citations}"),
        (f"仓库含 {N_ADDITIONAL_FILES} 个 Additional file",
         len(repo_files) == N_ADDITIONAL_FILES,
         f"found {len(repo_files)}"),
        ("仓库无旧 Table_S 命名残留", not leftover, f"leftover {leftover}"),
        ("已删除的两个文件确已移除", not removed, f"still present {removed}"),
        ("8 张 v5 主图存在", len(v5_figures) == 8,
         f"found {len(v5_figures)}"),
        ("BMC DOCX 有行号且无分页符", has_lines and no_page_break,
         f"lines={has_lines}, no_break={no_page_break}"),
    ]

    lines = [
        "# PAPPA2 BMC submission consistency report",
        "",
        f"Manuscript: `{MANUSCRIPT_MD.name}`",
        "",
        "| Check | Status | Detail |",
        "|---|---|---|",
    ]
    for label, passed, detail in checks:
        lines.append(f"| {label} | {'PASS' if passed else 'FAIL'} | {detail} |")
    lines.extend(["", "## Additional files in repository", "", "```text"])
    lines.extend(repo_files)
    lines.extend(["```", ""])

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    ok = True
    for label, passed, detail in checks:
        ok = ok and passed
        print(f"{'PASS' if passed else 'FAIL'}  {label}  -> {detail}")
    print("报告：", REPORT)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
