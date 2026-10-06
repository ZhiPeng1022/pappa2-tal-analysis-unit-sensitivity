import os
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(
    os.environ.get(
        "PAPPA2_PROJECT_ROOT",
        Path(__file__).resolve().parents[2],
    )
)
OUTPUT_ROOT = Path(
    os.environ.get(
        "PAPPA2_INTERNAL_ROOT",
        PROJECT_ROOT / "internal",
    )
)
MD_OUTPUT = (
    OUTPUT_ROOT / "PAPPA2_manuscript_evidence_matrix_v2.md"
)
CSV_OUTPUT = (
    OUTPUT_ROOT / "PAPPA2_manuscript_evidence_matrix_v2.csv"
)


ROWS = [
    {
        "id": "E01",
        "position": "Introduction, Methods",
        "claim": (
            "Donor-aware analysis-unit comparison can separate "
            "donor composition, neighborhood association, and "
            "within-donor state association in human TAL cells."
        ),
        "level": "Method framework",
        "status": "Primary framework",
        "evidence": (
            "Complete KPMP TAL source: 106,851 cells from four "
            "target subtypes; frozen neighborhood subset: 4,446 "
            "cells; four inference units compared."
        ),
        "allowed": (
            "We compared four analysis units in a donor-aware TAL "
            "framework."
        ),
        "forbidden": (
            "The framework is not claimed to outperform all existing "
            "neighborhood methods."
        ),
        "next": (
            "Preserve unit definitions and report all four units in "
            "the methods and figures."
        ),
    },
    {
        "id": "E02",
        "position": "Results 3.1",
        "claim": (
            "Disease-related TAL neighborhood states differ by subtype."
        ),
        "level": "B",
        "status": "Primary result 1",
        "evidence": (
            "all_disease: C-TAL-A 62/128, C-TAL-B 32/131, aTAL2 "
            "58/68, frTAL 67/84 significant neighborhoods. CKD aTAL2 "
            "median logFC -1.125 versus frTAL +1.244."
        ),
        "allowed": (
            "TAL subtypes show different disease-related neighborhood "
            "states and directions."
        ),
        "forbidden": (
            "Do not describe a single shared TAL response or uniform "
            "expansion."
        ),
        "next": (
            "Keep the formal Milo result separate from fixed-seed "
            "sensitivity results."
        ),
    },
    {
        "id": "E03",
        "position": "Results 3.3, Discussion",
        "claim": (
            "The four analysis units are not interchangeable."
        ),
        "level": "B/C",
        "status": "Primary method result",
        "evidence": (
            "Complete and subset pseudobulk agreed in 48/60 "
            "directions; 12 conflicts involved CXCL12, HES1, NOTCH3, "
            "and EPHB2. aTAL2 showed opposite PAPPA2 directions "
            "between donor pseudobulk and within-donor cell ranks."
        ),
        "allowed": (
            "Analysis-unit choice changes the question and can change "
            "the conclusion."
        ),
        "forbidden": (
            "Do not combine neighborhood and within-donor estimates "
            "as if they were one statistical test."
        ),
        "next": (
            "Use the concordance table as the central method figure."
        ),
    },
    {
        "id": "E04",
        "position": "Results 3.2",
        "claim": (
            "The 4,446-cell derivative subset is not representative "
            "of the complete TAL source."
        ),
        "level": "B",
        "status": "Coverage result",
        "evidence": (
            "aTAL2 Reference: 6,624 cells from 83 donors in complete "
            "data versus 102 cells from 38 donors in the subset. "
            "frTAL Reference: 6,656 versus 142 cells."
        ),
        "allowed": (
            "The subset preserves broad direction but not all gene- "
            "and subtype-level conclusions."
        ),
        "forbidden": (
            "Do not describe the 4,446-cell matrix as representative "
            "of all TAL cells."
        ),
        "next": (
            "Report complete and subset coverage in the main figure "
            "and supplement."
        ),
    },
    {
        "id": "E05",
        "position": "Results 3.2",
        "claim": (
            "Complete TAL donor pseudobulk provides a donor-level "
            "reference analysis."
        ),
        "level": "B",
        "status": "Primary donor-level result",
        "evidence": (
            "20/60 unadjusted and 25/60 age- and sex-adjusted tests "
            "had FDR < 0.05 across five genes, four subtypes, and "
            "three comparisons."
        ),
        "allowed": (
            "Complete donor pseudobulk is used as the donor-level "
            "reference for direction comparison."
        ),
        "forbidden": (
            "Do not treat donor pseudobulk as equivalent to "
            "neighborhood or within-donor analysis."
        ),
        "next": (
            "Keep unadjusted and adjusted estimates together."
        ),
    },
    {
        "id": "E06",
        "position": "Results 3.2, Discussion",
        "claim": (
            "Subset-to-complete direction changes identify sensitive "
            "gene and subtype conclusions."
        ),
        "level": "C",
        "status": "Method comparison",
        "evidence": (
            "48/60 directions agreed; 12 conflicts were concentrated "
            "in CXCL12, HES1, NOTCH3, and EPHB2. None involved "
            "PAPPA2."
        ),
        "allowed": (
            "A small derivative matrix can change individual test "
            "conclusions while preserving broad trends."
        ),
        "forbidden": (
            "Do not call the subset an external validation cohort."
        ),
        "next": (
            "Use the conflict table for method interpretation and "
            "supplementary reporting."
        ),
    },
    {
        "id": "E07",
        "position": "Results 3.4",
        "claim": (
            "PAPPA2 is a directionally stable case gene in complete "
            "TAL donor pseudobulk."
        ),
        "level": "B",
        "status": "Case result",
        "evidence": (
            "PAPPA2 was positive in 12/12 comparison-subtype groups "
            "with FDR < 0.05. Leave-one-donor analysis preserved "
            "direction in 12/12 groups; age/sex adjustment preserved "
            "the direction."
        ),
        "allowed": (
            "PAPPA2 marks a disease-associated TAL state."
        ),
        "forbidden": (
            "Do not infer that PAPPA2 causes the disease state."
        ),
        "next": (
            "Keep PAPPA2 as a case gene, not the sole study focus."
        ),
    },
    {
        "id": "E08",
        "position": "Results 3.3",
        "claim": (
            "Neighborhood-level PAPPA2 expression aligns with "
            "disease-related Milo logFC in most TAL groups."
        ),
        "level": "B/C",
        "status": "Supporting layer",
        "evidence": (
            "11/12 full-data Spearman correlations were positive; "
            "median 0.711; 10/12 donor-stratified intervals excluded "
            "zero; median purity-adjusted rho 0.723."
        ),
        "allowed": (
            "Neighborhood-level association is directionally stable "
            "after donor resampling and purity adjustment."
        ),
        "forbidden": (
            "Do not treat overlapping neighborhoods as independent "
            "cells or as causal evidence."
        ),
        "next": (
            "Retain resampling intervals and purity adjustment."
        ),
    },
    {
        "id": "E09",
        "position": "Results 3.3",
        "claim": (
            "Within-donor cell-rank analysis shows a distinct "
            "association structure."
        ),
        "level": "B/C",
        "status": "Supporting layer",
        "evidence": (
            "9/12 group medians were positive; median group-level rho "
            "0.205; 7/12 passed FDR. All 9 non-aTAL2 groups remained "
            "positive across cell-count thresholds."
        ),
        "allowed": (
            "Within-donor association reduces donor-level "
            "compositional differences."
        ),
        "forbidden": (
            "Do not interpret the permutation test as independent "
            "cell-level causality."
        ),
        "next": (
            "Keep permutation, threshold sensitivity, and donor "
            "counts together."
        ),
    },
    {
        "id": "E10",
        "position": "Results 3.3, Discussion",
        "claim": (
            "aTAL2 shows cross-unit discordance rather than one "
            "stable reverse mechanism."
        ),
        "level": "C/D",
        "status": "Method case",
        "evidence": (
            "Complete donor pseudobulk was positive in all three "
            "aTAL2 comparisons, whereas within-donor cell-rank "
            "medians were negative in all three; only the combined "
            "comparison passed FDR."
        ),
        "allowed": (
            "aTAL2 is a reproducible example of analysis-unit "
            "discordance."
        ),
        "forbidden": (
            "Do not call aTAL2 a stable reverse mechanism or a "
            "validated biological axis."
        ),
        "next": (
            "Make aTAL2 the central contrast in the method figure."
        ),
    },
    {
        "id": "E11",
        "position": "Results 3.5",
        "claim": (
            "PAPPA2-positive TAL fraction is associated with lower "
            "eGFR in KPMP donors."
        ),
        "level": "B",
        "status": "Internal clinical support",
        "evidence": (
            "Primary OR 0.819 per 10 percentage points (95% CI "
            "0.682-0.985, p=0.0335); disease OR 0.789; CKD OR 0.733; "
            "AKI was not associated."
        ),
        "allowed": (
            "KPMP internal models support an association with lower "
            "eGFR."
        ),
        "forbidden": (
            "Do not call the result independent clinical validation."
        ),
        "next": (
            "Report binned eGFR and the CKD concentration separately."
        ),
    },
    {
        "id": "E12",
        "position": "Results 3.6",
        "claim": (
            "PAPPA2 is directionally distinct among the tested "
            "control genes."
        ),
        "level": "B/C",
        "status": "Specificity support",
        "evidence": (
            "PAPPA2 had 7 FDR-supported within-donor comparisons, "
            "CXCL12 had 2, and NOTCH3, HES1, and EPHB2 had none. "
            "Nephroseq tubulointerstitial negative results: PAPPA2 "
            "7, HES1 6, EPHB2 3, NOTCH3 2, CXCL12 1."
        ),
        "allowed": (
            "PAPPA2 shows the most consistent negative GFR direction "
            "among the tested genes."
        ),
        "forbidden": (
            "Do not claim PAPPA2 is the only or validated GFR gene."
        ),
        "next": (
            "Retain control genes and avoid effect-size superiority "
            "claims."
        ),
    },
    {
        "id": "E13",
        "position": "Results 3.7",
        "claim": (
            "KPMP Visium provides same-consortium spatial "
            "cross-modality support."
        ),
        "level": "D",
        "status": "Supporting, not independent",
        "evidence": (
            "20/23 sections had co-occurrence OR > 1; DKD 11/11; "
            "donor-deduplicated 10/10; DKD versus Reference p=0.078; "
            "Reference 5/6."
        ),
        "allowed": (
            "KPMP spatial data show directional cross-modality "
            "consistency."
        ),
        "forbidden": (
            "Do not call GSE183456 an external cohort or a fully "
            "independent validation."
        ),
        "next": (
            "Keep the KPMP consortium role in the figure legend."
        ),
    },
    {
        "id": "E14",
        "position": "Results 3.8",
        "claim": (
            "GSE137570 provides limited external directional support."
        ),
        "level": "D",
        "status": "Limited external support",
        "evidence": (
            "24 donors: rho=-0.365, p=0.079; adjusted p=0.351; "
            "EPHB2 rho=-0.646, p=0.0006."
        ),
        "allowed": (
            "The external direction is consistent but underpowered."
        ),
        "forbidden": (
            "Do not claim significant or independent validation."
        ),
        "next": (
            "Keep this result clearly labelled as underpowered."
        ),
    },
    {
        "id": "E15",
        "position": "Results 3.8",
        "claim": (
            "Nephroseq provides multi-dataset negative PAPPA2-GFR "
            "association support."
        ),
        "level": "A-",
        "status": "External association support",
        "evidence": (
            "10/50 analyses were significant; 9 were negative, and "
            "7 came from tubulointerstitium."
        ),
        "allowed": (
            "Nephroseq supports a negative PAPPA2-GFR direction "
            "across precomputed kidney datasets."
        ),
        "forbidden": (
            "Do not count analysis rows as independent cohorts or "
            "call the result causal."
        ),
        "next": (
            "Deduplicate by dataset in any future update."
        ),
    },
    {
        "id": "E16",
        "position": "Discussion",
        "claim": (
            "Notch and EPHB2-STAT6 pathways are not consistently "
            "reproduced in TAL."
        ),
        "level": "E",
        "status": "Negative result",
        "evidence": (
            "JAG1, HES4, and STAT6 did not show consistent TAL "
            "changes; NOTCH3 and HES1 were subtype-dependent; EPHB2 "
            "directions were mixed."
        ),
        "allowed": (
            "The data do not support a unified TAL pathway model."
        ),
        "forbidden": (
            "Do not write that JAG1-NOTCH3-HES1 or EPHB2-STAT6 is "
            "established."
        ),
        "next": (
            "Keep as a limitation or negative result."
        ),
    },
    {
        "id": "E17",
        "position": "Discussion, Limitations",
        "claim": (
            "Current data support association and state structure, "
            "not causality."
        ),
        "level": "E",
        "status": "Boundary condition",
        "evidence": (
            "No perturbation, protein, or functional experiment; "
            "neighborhood and expression measurements are not "
            "independent causal instruments."
        ),
        "allowed": (
            "Use associated with and linked to."
        ),
        "forbidden": (
            "Do not use causes, drives, mediates, or regulates."
        ),
        "next": (
            "Apply the wording rule to the manuscript and figures."
        ),
    },
    {
        "id": "E18",
        "position": "Limitations",
        "claim": (
            "Mendelian randomization is not feasible with the "
            "available public protein resources."
        ),
        "level": "E",
        "status": "Data-availability limitation",
        "evidence": (
            "GTEx kidney/blood eQTL: 0; UKB-PPP: PAPPA2 absent; "
            "deCODE SomaScan: PAPPA2 present but only trans-pQTL; "
            "GWAS Catalog: no reported PAPPA2-eGFR/CKD association."
        ),
        "allowed": (
            "MR was not performed because no usable PAPPA2 "
            "cis-pQTL was available."
        ),
        "forbidden": (
            "Do not describe this as a negative MR result."
        ),
        "next": (
            "Revisit only if a public PAPPA2 cis-pQTL becomes "
            "available."
        ),
    },
    {
        "id": "E19",
        "position": "Data and Code Availability",
        "claim": (
            "The complete-TAL extraction and analysis-unit analyses "
            "are reproducible."
        ),
        "level": "A",
        "status": "Reproducibility",
        "evidence": (
            "Extraction, donor pseudobulk, robustness, and "
            "concordance scripts; versioned outputs; supplementary "
            "tables S9-S12."
        ),
        "allowed": (
            "Scripts and result files reproduce the reported "
            "complete-TAL analyses."
        ),
        "forbidden": (
            "Do not omit the difference between complete TAL and the "
            "4,446-cell subset."
        ),
        "next": (
            "Convert absolute local paths to relative paths before "
            "archival."
        ),
    },
    {
        "id": "E20",
        "position": "Abstract, Conclusion",
        "claim": (
            "Donor-aware analysis-unit comparison reveals "
            "subtype-specific TAL neighborhood states and exposes "
            "cross-unit discordance."
        ),
        "level": "Conditional synthesis",
        "status": "Main claim",
        "evidence": (
            "Complete-TAL pseudobulk, neighborhood association, "
            "within-donor cell-rank analysis, robustness checks, "
            "and PAPPA2 case results."
        ),
        "allowed": (
            "TAL states are subtype-specific, analysis units are not "
            "interchangeable, and PAPPA2 marks one associated state."
        ),
        "forbidden": (
            "Do not present PAPPA2 as a validated causal driver or "
            "the sole main finding."
        ),
        "next": (
            "Keep all figures and claims within this boundary."
        ),
    },
]


def render_markdown(frame):
    lines = [
        "# PAPPA2 论文证据矩阵 v2",
        "",
        "日期：2026-09-30",
        "对应主稿："
        "`submission/manuscript/PAPPA2_manuscript_draft_v4_TAL_method.md`",
        "",
        "## 一、证据等级",
        "",
        "- `A`：完整数据、外部关联或可复现性证据，但仍需注意"
        "数据集重叠和适用范围。",
        "- `B`：主队列中供者层级、调整后或完整数据结果。",
        "- `C`：跨分析单元或跨数据切分的方向一致性/冲突证据。",
        "- `D`：探索性、预计算或同联盟空间支持。",
        "- `E`：因果证据不足、阴性结果或限制条件。",
        "",
        "## 二、核心证据",
        "",
        "| ID | 位置 | 主张 | 等级 | 状态 | 关键证据 | 可写 | 不可写 | "
        "下一步 |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in frame.itertuples(index=False):
        lines.append(
            "| {id} | {position} | {claim} | {level} | {status} | "
            "{evidence} | {allowed} | {forbidden} | {next} |".format(
                id=row.id,
                position=row.position,
                claim=row.claim,
                level=row.level,
                status=row.status,
                evidence=row.evidence,
                allowed=row.allowed,
                forbidden=row.forbidden,
                next=row.next,
            )
        )

    lines.extend(
        [
            "",
            "## 三、当前主结论",
            "",
            "```text",
            "供者感知的分析单元比较显示，",
            "TAL 亚类的疾病相关邻域状态具有明显异质性；",
            "供者伪批量、邻域层面和供者内细胞秩相关不能互相替代；",
            "aTAL2 出现可报告的分析单元冲突；",
            "PAPPA2 标记其中一个疾病相关 TAL 状态，",
            "并与较低 eGFR/GFR 方向相关，但不支持因果关系。",
            "```",
            "",
            "## 四、关键词约束",
            "",
            "允许：`associated with`、`linked to`、"
            "`directionally consistent`、`analysis-unit discordance`。",
            "",
            "禁止：`causes`、`drives`、`mediates`、`regulates`、"
            "`validated mechanism`、`independent external cohort`。",
            "",
            "## 五、关键来源",
            "",
            "```text",
            "outputs/KPMP_full_TAL_focus_expression.csv.gz",
            "outputs/KPMP_full_TAL_focus_gene_pseudobulk.csv",
            "outputs/KPMP_full_vs_subset_TAL_focus_pseudobulk_comparison.csv",
            "outputs/KPMP_TAL_analysis_unit_concordance.csv",
            "outputs/KPMP_full_TAL_focus_pseudobulk_lodo.csv",
            "outputs/KPMP_full_TAL_focus_pseudobulk_bootstrap.csv",
            "outputs/KPMP_full_TAL_focus_pseudobulk_region_strata.csv",
            "```",
        ]
    )
    return "\n".join(lines) + "\n"


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(ROWS)
    frame.to_csv(
        CSV_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    MD_OUTPUT.write_text(
        render_markdown(frame),
        encoding="utf-8",
    )
    print("Saved:", MD_OUTPUT)
    print("Saved:", CSV_OUTPUT)
    print("Rows:", len(frame))


if __name__ == "__main__":
    main()
