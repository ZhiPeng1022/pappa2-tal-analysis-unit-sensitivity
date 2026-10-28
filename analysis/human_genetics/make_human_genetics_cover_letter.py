"""Generate the Human Genetics Research Article cover letter."""

from __future__ import annotations

from pathlib import Path

from docx import Document
from docx.shared import Cm, Pt


OUT_MD = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji\outputs"
    r"\PAPPA2_cover_letter_Human_Genetics_Research_Article_v1.md"
)
OUT_DOCX = OUT_MD.with_suffix(".docx")

DATE = "28 October 2026"
SALUTATION = "Dear Editors,"

SECTIONS = [
    (
        "Article type. ",
        "We submit our manuscript, \"Analysis-unit sensitivity distinguishes "
        "donor-level expression enrichment from within-donor coupling in "
        "human kidney thick ascending limb states\", for consideration as a "
        "Research Article.",
    ),
    (
        "The question. ",
        "Single-cell kidney studies commonly describe a disease-associated "
        "state as a property of a cell type, although donor composition, local "
        "neighborhood structure, and within-donor coupling answer different "
        "questions. We asked whether these analysis units lead to the same "
        "subtype-level conclusion in the human thick ascending limb (TAL).",
    ),
    (
        "What we found. ",
        "Using 106,851 KPMP single-nucleus TAL cells, we compared donor "
        "pseudobulk, Milo neighborhood abundance, neighborhood-expression "
        "association, and within-donor cell-rank association across four TAL "
        "subtypes. The units were not interchangeable. In aTAL2, donor-level "
        "and neighborhood-expression directions were positive, whereas the "
        "within-donor association was negative in combined disease and AKI "
        "(FDR = 0.0086 for each). Complete and subset pseudobulk directions "
        "agreed in 48 of 60 comparisons, but 12 changed direction. The analysis "
        "unit is therefore not a neutral technical choice; it changes which "
        "biological quantity is being estimated. A donor-level decomposition "
        "showed that the overall PAPPA2 difference was mainly a "
        "subtype-expression effect (+0.641) rather than a subtype-composition "
        "effect (+0.010). In a separate prespecified 400-gene panel, "
        "donor-versus-within reversal occurred in 18.6% of target genes versus "
        "12.7% of expression-matched random sets (empirical P = 0.000999).",
    ),
    (
        "Independent support. ",
        "We tested the key directions in two non-KPMP kidney resources. "
        "GSE211785 supported donor-level PAPPA2 disease elevation, and GSE254185 "
        "showed positive donor-level and neighborhood-abundance directions "
        "alongside negative neighborhood-expression and within-donor directions "
        "in injured TAL. A donor-stratified signature transfer achieved 93.9% "
        "accuracy among assigned held-out cells and 93.9% aTAL2 recall. The "
        "external data support a donor-level direction and the general "
        "analysis-unit problem, but do not provide exact aTAL2 replication.",
    ),
    (
        "Kidney regulatory context. ",
        "In NephQTL2, PAPPA2 was a strong tubulointerstitial eGene "
        "(FDR = 1.64e-14) without glomerular eQTL support. We also provide a "
        "prespecified boundary analysis: PAPPA2 was not among FDR-significant "
        "PTWAS genes, and same-cohort snATAC gene activity did not show higher "
        "disease accessibility. We treat this as regulatory context rather than "
        "genetic validation and state that boundary explicitly.",
    ),
    (
        "Why Human Genetics. ",
        "The manuscript connects human kidney single-cell genomics, donor-aware "
        "statistical inference, external cohort transfer, and kidney eQTL "
        "context. Its contribution is not a new algorithm or a claim that PAPPA2 "
        "causes disease. It is a rigorous demonstration that the unit of "
        "inference changes subtype-level conclusions and should be reported "
        "explicitly. We believe this fits the journal's readership in human "
        "genetics, functional genomics, and disease mechanism.",
    ),
    (
        "Declarations. ",
        "The manuscript is original, has not been published, and is not under "
        "consideration elsewhere. All authors have read and approved the "
        "submission. The authors declare no competing interests and no specific "
        "funding. The study used de-identified public or controlled-access data; "
        "no new human samples were collected. No clinical trial was performed. "
        "The manuscript is not intended for a special issue. Because we have no "
        "publication fund, we request consideration under the traditional "
        "subscription route and do not require Open Choice.",
    ),
]

CLOSING = [
    "We have not requested the exclusion of any reviewer.",
    "Thank you for considering our manuscript.",
]

SIGNATURE = [
    "Yours sincerely,",
    "Peng Zhi",
    "Department of General Practice, Henan Provincial Chest Hospital",
    "Zhengzhou, Henan, China",
    "Email: zhiyoyo0409@163.com",
    "ORCID: 0009-0003-2201-5131",
    "On behalf of all authors",
]


def main() -> None:
    markdown = [DATE, SALUTATION]
    markdown.extend(
        f"**{lead.strip()}** {text}" for lead, text in SECTIONS
    )
    markdown.extend(CLOSING)
    markdown.extend(SIGNATURE)
    OUT_MD.write_text("\n\n".join(markdown) + "\n", encoding="utf-8")

    document = Document()
    section = document.sections[0]
    section.page_width = Cm(21.0)
    section.page_height = Cm(29.7)
    for attr in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(section, attr, Cm(2.5))

    style = document.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(11)
    style.paragraph_format.space_after = Pt(8)
    style.paragraph_format.line_spacing = 1.15

    document.add_paragraph(DATE)
    document.add_paragraph(SALUTATION)
    for lead, text in SECTIONS:
        paragraph = document.add_paragraph()
        paragraph.add_run(lead).bold = True
        paragraph.add_run(text)
    for line in CLOSING:
        document.add_paragraph(line)
    for line in SIGNATURE:
        document.add_paragraph(line)

    document.save(OUT_DOCX)
    print(f"Wrote {OUT_MD}")
    print(f"Wrote {OUT_DOCX}")


if __name__ == "__main__":
    main()
