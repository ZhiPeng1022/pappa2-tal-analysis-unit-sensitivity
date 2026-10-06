import numpy as np
import pandas as pd

from project_paths import OUTPUT_ROOT

COMPARISONS = [
    "all_disease",
    "AKI_vs_Reference",
    "CKD_vs_Reference",
]
FOCUS_GENES = [
    "PAPPA2",
    "CXCL12",
    "NOTCH3",
    "HES1",
    "EPHB2",
]

PSEUDOBULK_PATH = (
    OUTPUT_ROOT / "KPMP_full_TAL_focus_gene_pseudobulk.csv"
)
NHOOD_ASSOCIATION_PATH = (
    OUTPUT_ROOT / "KPMP_full_TAL_nhood_focus_logFC_associations.csv"
)
WITHIN_DONOR_PATH = (
    OUTPUT_ROOT / "KPMP_full_TAL_within_donor_cell_rank.csv"
)
NHOOD_ROBUSTNESS_PATH = (
    OUTPUT_ROOT / "KPMP_full_TAL_nhood_focus_robustness.csv"
)
OLD_SUBSET_PATH = (
    OUTPUT_ROOT / "KPMP_TAL_analysis_unit_concordance.csv"
)

MILO_SUMMARY_PATH = (
    OUTPUT_ROOT / "KPMP_full_TAL_milo_subclass_summary.csv"
)
CONCORDANCE_PATH = (
    OUTPUT_ROOT / "KPMP_full_TAL_analysis_unit_concordance.csv"
)
REPORT_PATH = (
    OUTPUT_ROOT / "PAPPA2_full_data_reanalysis_summary.md"
)


def direction(values):
    values = pd.to_numeric(values, errors="coerce")
    return np.where(
        values > 0,
        1,
        np.where(values < 0, -1, 0),
    )


def markdown_table(frame):
    columns = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in frame.iterrows():
        values = []
        for value in row:
            if isinstance(value, float) and np.isfinite(value):
                text = f"{value:.4g}"
            else:
                text = str(value)
            values.append(text.replace("|", "\\|"))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def build_milo_summary():
    rows = []
    for comparison in COMPARISONS:
        path = OUTPUT_ROOT / (
            "KPMP_milo_full_TAL_final_"
            f"{comparison}_nhood_results.csv"
        )
        results = pd.read_csv(path)
        results["significant"] = (
            results["SpatialFDR"] < 0.05
        )
        for subclass, group in results.groupby(
            "subclass",
            dropna=False,
            observed=True,
        ):
            significant = group["significant"]
            rows.append(
                {
                    "comparison": comparison,
                    "subclass": subclass,
                    "total_annotated_nhoods": len(group),
                    "significant_nhoods": int(
                        significant.sum()
                    ),
                    "positive_nhoods": int(
                        (
                            significant
                            & (group["logFC"] > 0)
                        ).sum()
                    ),
                    "negative_nhoods": int(
                        (
                            significant
                            & (group["logFC"] < 0)
                        ).sum()
                    ),
                    "median_logFC": group["logFC"].median(),
                    "median_subclass_fraction": group[
                        "subclass_fraction"
                    ].median(),
                    "median_nhood_purity": group[
                        "nhood_purity"
                    ].median(),
                }
            )
    summary = pd.DataFrame(rows).sort_values(
        ["comparison", "subclass"]
    )
    summary.to_csv(
        MILO_SUMMARY_PATH,
        index=False,
        encoding="utf-8-sig",
    )
    return summary


def build_concordance():
    pseudobulk = pd.read_csv(PSEUDOBULK_PATH)
    pseudobulk = pseudobulk[
        pseudobulk["comparison"].isin(COMPARISONS)
        & pseudobulk["gene"].isin(FOCUS_GENES)
    ].copy()
    pseudobulk["full_pb_direction"] = direction(
        pseudobulk["mean_case_minus_reference"]
    )
    pseudobulk = pseudobulk.rename(
        columns={
            "mean_case_minus_reference": "full_pb_mean_diff",
            "FDR": "full_pb_FDR",
            "significant_FDR_0_05": "full_pb_significant",
            "adjusted_condition_beta": "full_pb_adjusted_beta",
            "adjusted_FDR": "full_pb_adjusted_FDR",
            "adjusted_significant_FDR_0_05": (
                "full_pb_adjusted_significant"
            ),
        }
    )

    nhood = pd.read_csv(NHOOD_ASSOCIATION_PATH)
    nhood = nhood[
        nhood["comparison"].isin(COMPARISONS)
        & nhood["gene"].isin(FOCUS_GENES)
    ].copy()
    nhood["nhood_direction"] = direction(
        nhood["spearman_rho"]
    )
    nhood = nhood.rename(
        columns={
            "spearman_rho": "nhood_rho",
            "spearman_p_naive": "nhood_p_naive",
            "significant_nhoods": "nhood_significant_count",
            "positive_significant_nhoods": (
                "nhood_positive_significant_count"
            ),
            "negative_significant_nhoods": (
                "nhood_negative_significant_count"
            ),
        }
    )

    within = pd.read_csv(WITHIN_DONOR_PATH)
    within = within[
        within["comparison"].isin(COMPARISONS)
        & within["gene"].isin(FOCUS_GENES)
    ].copy()
    within["within_donor_direction"] = direction(
        within["observed_disease_minus_reference"]
    )
    within = within.rename(
        columns={
            "observed_median_rho_all": "within_donor_rho_all",
            "observed_median_rho_reference": (
                "within_donor_rho_reference"
            ),
            "observed_median_rho_disease": (
                "within_donor_rho_disease"
            ),
            "observed_disease_minus_reference": (
                "within_donor_disease_minus_reference"
            ),
            "FDR_p_disease_minus_reference": (
                "within_donor_difference_FDR"
            ),
        }
    )

    robustness = pd.read_csv(NHOOD_ROBUSTNESS_PATH)
    robustness = robustness[
        robustness["comparison"].isin(COMPARISONS)
        & robustness["gene"].isin(FOCUS_GENES)
    ][
        [
            "comparison",
            "subclass",
            "gene",
            "rho_full",
            "lodo_same_sign_fraction",
            "bootstrap_ci_lower",
            "bootstrap_ci_upper",
            "bootstrap_fraction_positive",
            "partial_rho_controlling_purity",
        ]
    ]

    old = pd.read_csv(OLD_SUBSET_PATH)
    old = old[
        [
            "comparison",
            "subclass",
            "gene",
            "subset_pb_mean_diff",
            "subset_pb_FDR",
            "subset_pb_direction",
        ]
    ]

    keep_pb = [
        "comparison",
        "subclass",
        "gene",
        "reference_donors",
        "case_donors",
        "reference_mean",
        "case_mean",
        "full_pb_mean_diff",
        "full_pb_FDR",
        "full_pb_significant",
        "full_pb_adjusted_beta",
        "full_pb_adjusted_FDR",
        "full_pb_adjusted_significant",
        "full_pb_direction",
    ]
    keep_nhood = [
        "comparison",
        "subclass",
        "gene",
        "n_nhoods",
        "n_unique_cells_in_subclass",
        "median_nhood_cells",
        "nhood_significant_count",
        "nhood_positive_significant_count",
        "nhood_negative_significant_count",
        "median_logFC",
        "median_expression",
        "nhood_rho",
        "nhood_p_naive",
        "nhood_direction",
    ]
    keep_within = [
        "comparison",
        "subclass",
        "gene",
        "eligible_donors",
        "eligible_reference_donors",
        "eligible_disease_donors",
        "within_donor_rho_all",
        "within_donor_rho_reference",
        "within_donor_rho_disease",
        "within_donor_disease_minus_reference",
        "within_donor_difference_FDR",
        "within_donor_direction",
    ]

    concordance = (
        pseudobulk[keep_pb]
        .merge(
            nhood[keep_nhood],
            on=["comparison", "subclass", "gene"],
            how="inner",
            validate="one_to_one",
        )
        .merge(
            within[keep_within],
            on=["comparison", "subclass", "gene"],
            how="inner",
            validate="one_to_one",
        )
        .merge(
            robustness,
            on=["comparison", "subclass", "gene"],
            how="left",
            validate="one_to_one",
        )
        .merge(
            old,
            on=["comparison", "subclass", "gene"],
            how="left",
            validate="one_to_one",
        )
    )

    direction_columns = [
        "full_pb_direction",
        "nhood_direction",
        "within_donor_direction",
    ]
    concordance["direction_agreement_count"] = (
        concordance[direction_columns]
        .eq(
            concordance["full_pb_direction"],
            axis=0,
        )
        .sum(axis=1)
    )
    concordance["direction_conflict"] = (
        concordance[direction_columns].nunique(axis=1)
        > 1
    )
    concordance["full_vs_subset_same_direction"] = (
        concordance["full_pb_direction"]
        == concordance["subset_pb_direction"]
    )
    concordance["full_vs_nhood_same_direction"] = (
        concordance["full_pb_direction"]
        == concordance["nhood_direction"]
    )
    concordance[
        "full_vs_within_donor_same_direction"
    ] = (
        concordance["full_pb_direction"]
        == concordance["within_donor_direction"]
    )
    concordance[
        "nhood_vs_within_donor_same_direction"
    ] = (
        concordance["nhood_direction"]
        == concordance["within_donor_direction"]
    )
    concordance["within_donor_rho"] = concordance[
        "within_donor_disease_minus_reference"
    ]
    concordance = concordance.sort_values(
        ["comparison", "subclass", "gene"]
    )
    concordance.to_csv(
        CONCORDANCE_PATH,
        index=False,
        encoding="utf-8-sig",
    )
    return concordance


def write_report(milo_summary, concordance):
    pappa2 = concordance[
        concordance["gene"] == "PAPPA2"
    ].copy()
    robustness = pd.read_csv(NHOOD_ROBUSTNESS_PATH)
    pappa2_robustness = robustness[
        robustness["gene"] == "PAPPA2"
    ].copy()
    pb_positive = int(
        (pappa2["full_pb_mean_diff"] > 0).sum()
    )
    pb_fdr = int(
        pappa2["full_pb_significant"].sum()
    )
    bootstrap_positive = int(
        (
            pappa2_robustness["bootstrap_ci_lower"] > 0
        ).sum()
    )
    lines = [
        "# PAPPA2 全量数据再分析摘要 v1",
        "",
        "日期：2026-09-30",
        "",
        "## 方法与输入",
        "",
        "- 完整 TAL：106,851 个细胞，223 位供者；",
        "- 使用完整 h5ad 中已有 UMAP；",
        "- Milo 使用 miloR 2.8.1，固定种子 20260927；",
        "- 归一化方法：logMS；",
        "- AKI 比较剔除 1 位没有邻域贡献的供者 163-5；",
        "- 全量 neighborhood membership 已导出；",
        "- 全量供者内秩相关使用 1,000 次置换。",
        "",
        "## 全量 Milo 汇总",
        "",
        markdown_table(milo_summary),
        "",
        "## PAPPA2 四层结果",
        "",
        pappa2[
            [
                "comparison",
                "subclass",
                "full_pb_mean_diff",
                "full_pb_FDR",
                "nhood_rho",
                "nhood_direction",
                "within_donor_rho_reference",
                "within_donor_rho_disease",
                "within_donor_disease_minus_reference",
                "within_donor_difference_FDR",
                "direction_conflict",
            ]
        ]
        .pipe(markdown_table),
        "",
        "## PAPPA2 完整邻域稳健性",
        "",
        markdown_table(
            pappa2_robustness[
                [
                    "comparison",
                    "subclass",
                    "rho_full",
                    "lodo_same_sign_fraction",
                    "bootstrap_ci_lower",
                    "bootstrap_ci_upper",
                    "bootstrap_fraction_positive",
                    "partial_rho_controlling_purity",
                ]
            ]
        ),
        "",
        "## PAPPA2 方向汇总",
        "",
        f"- 供者伪批量正向：{pb_positive}/12；",
        f"- 供者伪批量 FDR < 0.05：{pb_fdr}/12；",
        f"- 邻域 bootstrap 区间完全为正：{bootstrap_positive}/12。",
        "",
        "## 关键解释",
        "",
        "1. 完整数据中 PAPPA2 的供者伪批量和邻域表达关联仍为正向；",
        "2. aTAL2 的供者内秩相关在 all_disease 和 AKI 中相对负向，"
        "并在疾病组与 Reference 组差异检验中通过 FDR；",
        "3. CKD 的 aTAL2 方向差为负，但差异检验未通过 FDR；",
        "4. 因此应写成“跨分析单元冲突在完整数据中仍存在，但显著程度"
        "依疾病比较而异”，不能写成三个比较全部显著；",
        "5. 全量 Milo 与旧子集 Milo 不能直接按显著数量比较，应"
        "分别报告输入细胞数、供者数和效应量。",
        "",
    ]
    REPORT_PATH.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main():
    milo_summary = build_milo_summary()
    concordance = build_concordance()
    write_report(milo_summary, concordance)

    print("Saved:", MILO_SUMMARY_PATH)
    print("Saved:", CONCORDANCE_PATH)
    print("Saved:", REPORT_PATH)
    print("\nFull Milo summary:")
    print(milo_summary.to_string(index=False))
    print("\nPAPPA2 full-data direction conflicts:")
    print(
        concordance.loc[
            (concordance["gene"] == "PAPPA2")
            & concordance["direction_conflict"],
            [
                "comparison",
                "subclass",
                "full_pb_mean_diff",
                "nhood_rho",
                "within_donor_disease_minus_reference",
                "within_donor_difference_FDR",
            ],
        ]
        .round(4)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
