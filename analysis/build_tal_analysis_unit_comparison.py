import os
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(
    os.environ.get(
        "PAPPA2_PROJECT_ROOT",
        Path(__file__).resolve().parents[2],
    )
)
OUTPUT_ROOT = Path(
    os.environ.get(
        "PAPPA2_OUTPUT_ROOT",
        PROJECT_ROOT / "outputs",
    )
)

FULL_PSEUDOBULK = (
    OUTPUT_ROOT / "KPMP_full_TAL_focus_gene_pseudobulk.csv"
)
SUBSET_PSEUDOBULK = (
    OUTPUT_ROOT / "KPMP_TAL_focus_gene_pseudobulk.csv"
)
NEIGHBORHOOD_ROBUSTNESS = (
    OUTPUT_ROOT
    / "KPMP_TAL_seed20260927_nhood_focus_robustness.csv"
)
WITHIN_DONOR = (
    OUTPUT_ROOT
    / "KPMP_TAL_seed20260927_within_donor_cell_rank.csv"
)

OUTPUT_CSV = (
    OUTPUT_ROOT / "KPMP_TAL_analysis_unit_concordance.csv"
)
OUTPUT_MD = (
    OUTPUT_ROOT
    / "KPMP_TAL_analysis_unit_concordance_summary.md"
)

FOCUS_GENES = [
    "PAPPA2",
    "CXCL12",
    "NOTCH3",
    "HES1",
    "EPHB2",
]
SUBCLASS_ORDER = [
    "C-TAL-A",
    "C-TAL-B",
    "aTAL2",
    "frTAL",
]
COMPARISON_ORDER = [
    "all_disease",
    "AKI_vs_Reference",
    "CKD_vs_Reference",
]


def sign_value(values):
    values = np.asarray(values, dtype=float)
    return np.where(
        values > 0,
        1,
        np.where(values < 0, -1, 0),
    )


def load_full_pseudobulk():
    frame = pd.read_csv(FULL_PSEUDOBULK)
    frame = frame[
        frame["gene"].isin(FOCUS_GENES)
    ].copy()
    frame = frame.rename(
        columns={
            "mean_case_minus_reference": (
                "full_pb_mean_diff"
            ),
            "FDR": "full_pb_FDR",
            "significant_FDR_0_05": (
                "full_pb_significant"
            ),
            "adjusted_condition_beta": (
                "full_pb_adjusted_beta"
            ),
            "adjusted_FDR": "full_pb_adjusted_FDR",
            "adjusted_significant_FDR_0_05": (
                "full_pb_adjusted_significant"
            ),
        }
    )
    frame["full_pb_direction"] = sign_value(
        frame["full_pb_mean_diff"]
    )
    return frame[
        [
            "comparison",
            "subclass",
            "gene",
            "reference_donors",
            "case_donors",
            "full_pb_mean_diff",
            "full_pb_FDR",
            "full_pb_significant",
            "full_pb_adjusted_beta",
            "full_pb_adjusted_FDR",
            "full_pb_adjusted_significant",
            "full_pb_direction",
        ]
    ]


def load_subset_pseudobulk():
    frame = pd.read_csv(SUBSET_PSEUDOBULK)
    frame = frame[
        frame["gene"].isin(FOCUS_GENES)
    ].copy()
    frame["comparison"] = frame["comparison"].replace(
        {"disease_vs_Reference": "all_disease"}
    )
    frame = frame.rename(
        columns={
            "mean_case_minus_reference": (
                "subset_pb_mean_diff"
            ),
            "FDR": "subset_pb_FDR",
            "significant_FDR_0_05": (
                "subset_pb_significant"
            ),
        }
    )
    frame["subset_pb_direction"] = sign_value(
        frame["subset_pb_mean_diff"]
    )
    return frame[
        [
            "comparison",
            "subclass",
            "gene",
            "subset_pb_mean_diff",
            "subset_pb_FDR",
            "subset_pb_significant",
            "subset_pb_direction",
        ]
    ]


def load_neighborhood():
    frame = pd.read_csv(NEIGHBORHOOD_ROBUSTNESS)
    frame = frame[
        frame["gene"].isin(FOCUS_GENES)
    ].copy()
    frame["nhood_direction"] = sign_value(
        frame["rho_full"]
    )
    frame["nhood_ci_excludes_zero"] = (
        (frame["bootstrap_ci_lower"] > 0)
        | (frame["bootstrap_ci_upper"] < 0)
    )
    frame["nhood_stable_direction"] = (
        frame["nhood_ci_excludes_zero"]
        & (frame["nhood_direction"] != 0)
    )
    return frame[
        [
            "comparison",
            "subclass",
            "gene",
            "nhoods",
            "unique_cells",
            "donors",
            "rho_full",
            "bootstrap_ci_lower",
            "bootstrap_ci_upper",
            "bootstrap_fraction_positive",
            "lodo_same_sign_fraction",
            "partial_rho_controlling_purity",
            "nhood_direction",
            "nhood_ci_excludes_zero",
            "nhood_stable_direction",
        ]
    ]


def load_within_donor():
    frame = pd.read_csv(WITHIN_DONOR)
    frame = frame[
        frame["gene"].isin(FOCUS_GENES)
    ].copy()
    frame["within_donor_direction"] = sign_value(
        frame["observed_median_rho_all"]
    )
    frame = frame.rename(
        columns={
            "observed_median_rho_all": (
                "within_donor_rho"
            ),
            "FDR_p_all": "within_donor_FDR",
            "eligible_donors": (
                "within_donor_eligible_donors"
            ),
            "eligible_reference_donors": (
                "within_donor_reference_donors"
            ),
            "eligible_disease_donors": (
                "within_donor_disease_donors"
            ),
        }
    )
    frame["within_donor_significant"] = (
        frame["within_donor_FDR"] < 0.05
    )
    return frame[
        [
            "comparison",
            "subclass",
            "gene",
            "within_donor_rho",
            "within_donor_FDR",
            "within_donor_significant",
            "within_donor_eligible_donors",
            "within_donor_reference_donors",
            "within_donor_disease_donors",
            "within_donor_direction",
        ]
    ]


def add_layer_summary(frame):
    direction_columns = [
        "subset_pb_direction",
        "full_pb_direction",
        "nhood_direction",
        "within_donor_direction",
    ]
    frame["directions_available"] = frame[
        direction_columns
    ].notna().sum(axis=1)
    frame["positive_directions"] = (
        frame[direction_columns] == 1
    ).sum(axis=1)
    frame["negative_directions"] = (
        frame[direction_columns] == -1
    ).sum(axis=1)
    frame["direction_consensus"] = np.select(
        [
            (
                frame["positive_directions"]
                > frame["negative_directions"]
            ),
            (
                frame["negative_directions"]
                > frame["positive_directions"]
            ),
        ],
        [1, -1],
        default=0,
    )
    frame["direction_agreement_count"] = np.where(
        frame["direction_consensus"] == 1,
        frame["positive_directions"],
        np.where(
            frame["direction_consensus"] == -1,
            frame["negative_directions"],
            0,
        ),
    )
    frame["direction_conflict"] = (
        (frame["positive_directions"] > 0)
        & (frame["negative_directions"] > 0)
    )

    frame["full_vs_subset_same_direction"] = (
        frame["full_pb_direction"]
        == frame["subset_pb_direction"]
    )
    frame["full_vs_nhood_same_direction"] = (
        frame["full_pb_direction"]
        == frame["nhood_direction"]
    )
    frame["full_vs_within_donor_same_direction"] = (
        frame["full_pb_direction"]
        == frame["within_donor_direction"]
    )
    frame["nhood_vs_within_donor_same_direction"] = (
        frame["nhood_direction"]
        == frame["within_donor_direction"]
    )

    frame["donor_level_support_count"] = (
        frame["full_pb_significant"].astype(int)
        + frame["within_donor_significant"].astype(int)
    )
    frame["evidence_layer_support_count"] = (
        frame["full_pb_significant"].astype(int)
        + frame["nhood_stable_direction"].astype(int)
        + frame["within_donor_significant"].astype(int)
    )
    return frame


def markdown_table(frame):
    headers = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]
    for _, row in frame.iterrows():
        values = []
        for value in row:
            if pd.isna(value):
                values.append("")
            elif isinstance(value, (float, np.floating)):
                values.append(f"{value:.4f}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def build_summary(frame):
    pappa2 = frame[
        frame["gene"] == "PAPPA2"
    ].sort_values(
        [
            "comparison",
            "subclass",
        ]
    )
    gene_summary = (
        frame.groupby("gene", observed=True)
        .agg(
            tests=("gene", "size"),
            full_pb_significant=(
                "full_pb_significant",
                "sum",
            ),
            nhood_stable=(
                "nhood_stable_direction",
                "sum",
            ),
            within_donor_significant=(
                "within_donor_significant",
                "sum",
            ),
            full_vs_subset_direction_agreement=(
                "full_vs_subset_same_direction",
                "sum",
            ),
            full_vs_nhood_direction_agreement=(
                "full_vs_nhood_same_direction",
                "sum",
            ),
            full_vs_within_donor_direction_agreement=(
                "full_vs_within_donor_same_direction",
                "sum",
            ),
        )
        .reset_index()
        .sort_values("gene")
    )
    comparison_summary = (
        frame.groupby(
            ["comparison", "subclass"],
            observed=True,
        )
        .agg(
            tests=("gene", "size"),
            direction_conflicts=(
                "direction_conflict",
                "sum",
            ),
            full_pb_significant=(
                "full_pb_significant",
                "sum",
            ),
            nhood_stable=(
                "nhood_stable_direction",
                "sum",
            ),
            within_donor_significant=(
                "within_donor_significant",
                "sum",
            ),
        )
        .reset_index()
    )
    pappa2_view = pappa2[
        [
            "comparison",
            "subclass",
            "full_pb_mean_diff",
            "full_pb_FDR",
            "subset_pb_mean_diff",
            "subset_pb_FDR",
            "rho_full",
            "bootstrap_ci_lower",
            "bootstrap_ci_upper",
            "within_donor_rho",
            "within_donor_FDR",
            "direction_consensus",
            "direction_conflict",
        ]
    ].copy()

    return f"""# TAL 分析单元一致性汇总

日期：2026-09-30

## 一、分析单元

本表比较四个分析单元：

```text
subset_pb = 旧 4,446 细胞子集的供者伪批量
full_pb = 完整 106,851 细胞 TAL 的供者伪批量
nhood = 固定种子 Milo 邻域平均表达与 logFC 的关联
within_donor = 供者内部细胞表达与邻域暴露的秩相关
```

邻域层面的稳定性使用 500 次供者分层重采样区间；
完整供者伪批量主结果为 Mann-Whitney 检验，年龄和性别调整结果
另行保存在 `KPMP_full_TAL_focus_gene_pseudobulk.csv`。

## 二、PAPPA2 的四层比较

{markdown_table(pappa2_view)}

## 三、不同基因的层级支持数量

`full_vs_*_direction_agreement` 表示与完整供者伪批量方向一致的
检验数，共 12 个比较组合。

{markdown_table(gene_summary)}

## 四、比较和亚类层面的冲突

{markdown_table(comparison_summary)}

## 五、可以直接写入方法论文的结果

1. 完整 TAL 的供者伪批量可用于减少旧子集选择偏倚；
2. 不同分析单元回答不同问题，不能互相替代；
3. 邻域层面、供者伪批量和供者内细胞层面的一致性应作为主要
   方法指标，而不是只报告某一个层面的显著性；
4. 分析单元冲突本身可以作为方法学结果，用来提示组成差异、
   供者异质性或子集选择偏倚；
5. PAPPA2 在完整 TAL 的四个亚类中均与疾病状态升高相关，但
   这仍然不是因果证据。

## 六、必须保留的限制

1. 旧 4,446 细胞子集不是完整 TAL 的随机样本；
2. 邻域层面分析仍基于固定种子 Milo 结果，不应写成正式
   Milo 差异丰度结论；
3. 供者内细胞秩相关保留邻域重叠，不能当作独立样本检验；
4. 年龄和性别调整模型仍是观察性关联；
5. 当前没有外部独立队列或功能实验。
"""


def main():
    full = load_full_pseudobulk()
    subset = load_subset_pseudobulk()
    nhood = load_neighborhood()
    within = load_within_donor()

    merged = (
        full.merge(
            subset,
            on=["comparison", "subclass", "gene"],
            how="outer",
            validate="one_to_one",
        )
        .merge(
            nhood,
            on=["comparison", "subclass", "gene"],
            how="outer",
            validate="one_to_one",
        )
        .merge(
            within,
            on=["comparison", "subclass", "gene"],
            how="outer",
            validate="one_to_one",
        )
    )
    merged = add_layer_summary(merged)
    merged["comparison"] = pd.Categorical(
        merged["comparison"],
        categories=COMPARISON_ORDER,
        ordered=True,
    )
    merged["subclass"] = pd.Categorical(
        merged["subclass"],
        categories=SUBCLASS_ORDER,
        ordered=True,
    )
    merged = merged.sort_values(
        ["comparison", "subclass", "gene"]
    )
    merged.to_csv(
        OUTPUT_CSV,
        index=False,
        encoding="utf-8-sig",
    )
    OUTPUT_MD.write_text(
        build_summary(merged),
        encoding="utf-8",
    )

    print("Saved:", OUTPUT_CSV)
    print("Saved:", OUTPUT_MD)
    print("\nPAPPA2 summary:")
    pappa2 = merged[
        merged["gene"] == "PAPPA2"
    ]
    print(
        pappa2[
            [
                "comparison",
                "subclass",
                "full_pb_mean_diff",
                "full_pb_FDR",
                "subset_pb_mean_diff",
                "subset_pb_FDR",
                "rho_full",
                "bootstrap_ci_lower",
                "bootstrap_ci_upper",
                "within_donor_rho",
                "within_donor_FDR",
                "direction_conflict",
            ]
        ]
        .round(4)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
