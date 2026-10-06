from pathlib import Path

import anndata as ad
import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2")
OUTPUT_DIR = ROOT / "outputs"
INPUT_PATH = Path(
    r"D:\CodexData\RAMP3\augur_tal_full_input"
) / "KPMP_augur_tal_full_input.h5ad"

SUBCLASS_ORDER = [
    "C-TAL-A",
    "C-TAL-B",
    "aTAL2",
    "frTAL",
]
COMPARISON_ORDER = [
    "disease_vs_Reference",
    "AKI_vs_Reference",
    "CKD_vs_Reference",
]
MILO_COMPARISON_ORDER = [
    "all_disease",
    "AKI_vs_Reference",
    "CKD_vs_Reference",
]
FOCUS_GENES = [
    "PAPPA2",
    "JAG1",
    "NOTCH3",
    "HES1",
    "HES4",
    "CXCL12",
    "LHFPL3",
    "NDUFS6",
    "EPHB2",
    "STAT6",
]


def ordered_categorical(frame, column, values):
    result = frame.copy()
    result[column] = pd.Categorical(
        result[column],
        categories=values,
        ordered=True,
    )
    return result


def summarize_milo():
    frames = []
    for path in sorted(
        OUTPUT_DIR.glob(
            "KPMP_milo_TAL_*_nhood_results.csv"
        )
    ):
        frames.append(pd.read_csv(path))
    nhoods = pd.concat(frames, ignore_index=True)
    significant = nhoods[
        nhoods["SpatialFDR"] < 0.05
    ].copy()
    significant["direction"] = significant["logFC"].map(
        lambda value: (
            "positive" if value > 0 else "negative"
        )
    )

    summary = (
        significant.groupby(
            ["comparison", "subclass"],
            observed=True,
        )
        .agg(
            significant_nhoods=("Nhood", "size"),
            positive_nhoods=(
                "direction",
                lambda values: (values == "positive").sum(),
            ),
            negative_nhoods=(
                "direction",
                lambda values: (values == "negative").sum(),
            ),
            median_logFC=("logFC", "median"),
            median_subclass_fraction=(
                "subclass_fraction",
                "median",
            ),
        )
        .reset_index()
    )
    totals = (
        nhoods.groupby(
            ["comparison", "subclass"],
            observed=True,
        )
        .size()
        .rename("total_annotated_nhoods")
        .reset_index()
    )
    summary = totals.merge(
        summary,
        on=["comparison", "subclass"],
        how="left",
    )
    summary = summary.fillna(
        {
            "significant_nhoods": 0,
            "positive_nhoods": 0,
            "negative_nhoods": 0,
        }
    )
    summary = ordered_categorical(
        summary,
        "comparison",
        MILO_COMPARISON_ORDER,
    )
    summary = ordered_categorical(
        summary,
        "subclass",
        SUBCLASS_ORDER,
    )
    return summary.sort_values(
        ["comparison", "subclass"]
    ).reset_index(drop=True)


def summarize_augur():
    auc = pd.read_csv(
        OUTPUT_DIR / "KPMP_Augur_TAL_full_all_AUC.csv"
    )
    adata = ad.read_h5ad(INPUT_PATH, backed="r")
    cell_counts = (
        adata.obs.groupby(
            ["condition", "subclass"],
            observed=True,
        )
        .size()
        .rename("cells")
        .reset_index()
    )

    rows = []
    for _, auc_row in auc.iterrows():
        comparison = auc_row["comparison"]
        subclass = auc_row["cell_type"]
        if comparison == "disease_vs_Reference":
            case_conditions = ["AKI", "CKD"]
        elif comparison == "AKI_vs_Reference":
            case_conditions = ["AKI"]
        else:
            case_conditions = ["CKD"]

        reference_count = cell_counts[
            (cell_counts["condition"] == "Reference")
            & (cell_counts["subclass"] == subclass)
        ]["cells"].sum()
        case_count = cell_counts[
            (cell_counts["condition"].isin(case_conditions))
            & (cell_counts["subclass"] == subclass)
        ]["cells"].sum()

        subsample_path = (
            OUTPUT_DIR
            / (
                "KPMP_Augur_TAL_full_"
                f"{comparison}_subsample_results.csv"
            )
        )
        subsamples = pd.read_csv(subsample_path)
        estimates = subsamples[
            (subsamples["cell_type"] == subclass)
            & (subsamples["metric"] == "roc_auc")
        ]["estimate"]
        rows.append(
            {
                "comparison": comparison,
                "subclass": subclass,
                "reference_cells": int(reference_count),
                "case_cells": int(case_count),
                "auc": auc_row["auc"],
                "auc_sd": estimates.std(),
                "auc_q025": estimates.quantile(0.025),
                "auc_q975": estimates.quantile(0.975),
                "subsample_fold_estimates": len(estimates),
            }
        )
    summary = pd.DataFrame(rows)
    summary = ordered_categorical(
        summary,
        "comparison",
        COMPARISON_ORDER,
    )
    summary = ordered_categorical(
        summary,
        "subclass",
        SUBCLASS_ORDER,
    )
    return summary.sort_values(
        ["comparison", "subclass"]
    ).reset_index(drop=True)


def summarize_limited_vs_full():
    limited = pd.read_csv(
        OUTPUT_DIR / "KPMP_Augur_TAL_all_AUC.csv"
    ).rename(
        columns={
            "cell_type": "subclass",
            "auc": "limited_auc",
            "cells": "limited_cells",
        }
    )
    full = pd.read_csv(
        OUTPUT_DIR / "KPMP_Augur_TAL_full_all_AUC.csv"
    ).rename(
        columns={
            "cell_type": "subclass",
            "auc": "full_auc",
            "cells": "full_cells",
        }
    )
    comparison = limited.merge(
        full,
        on=["comparison", "subclass"],
        how="inner",
    )
    comparison["auc_change_full_minus_limited"] = (
        comparison["full_auc"]
        - comparison["limited_auc"]
    )
    comparison = ordered_categorical(
        comparison,
        "comparison",
        COMPARISON_ORDER,
    )
    comparison = ordered_categorical(
        comparison,
        "subclass",
        SUBCLASS_ORDER,
    )
    return comparison.sort_values(
        ["comparison", "subclass"]
    )[
        [
            "comparison",
            "subclass",
            "limited_cells",
            "full_cells",
            "limited_auc",
            "full_auc",
            "auc_change_full_minus_limited",
        ]
    ].reset_index(drop=True)


def summarize_importance():
    frames = []
    for path in sorted(
        OUTPUT_DIR.glob(
            "KPMP_Augur_TAL_full_"
            "*_feature_importance.csv"
        )
    ):
        frames.append(pd.read_csv(path))
    importance = pd.concat(frames, ignore_index=True)
    aggregated = (
        importance.groupby(
            ["comparison", "cell_type", "gene_symbol"],
            observed=True,
        )["importance"]
        .agg(["mean", "median", "max", "count"])
        .reset_index()
        .rename(
            columns={
                "cell_type": "subclass",
                "mean": "mean_importance",
                "median": "median_importance",
                "max": "max_importance",
                "count": "importance_rows",
            }
        )
    )
    aggregated["rank"] = (
        aggregated.groupby(
            ["comparison", "subclass"],
            observed=True,
        )["mean_importance"]
        .rank(
            method="first",
            ascending=False,
        )
        .astype(int)
    )
    top = (
        aggregated.sort_values(
            ["comparison", "subclass", "rank"]
        )
        .groupby(
            ["comparison", "subclass"],
            observed=True,
        )
        .head(10)
        .reset_index(drop=True)
    )

    focus_rows = []
    for _, auc_row in summarize_augur().iterrows():
        comparison = auc_row["comparison"]
        subclass = auc_row["subclass"]
        subclass_data = aggregated[
            (aggregated["comparison"] == comparison)
            & (aggregated["subclass"] == subclass)
        ]
        for gene in FOCUS_GENES:
            hit = subclass_data[
                subclass_data["gene_symbol"] == gene
            ]
            if len(hit) == 0:
                focus_rows.append(
                    {
                        "comparison": comparison,
                        "subclass": subclass,
                        "gene": gene,
                        "present_in_model": False,
                        "rank": pd.NA,
                        "mean_importance": pd.NA,
                    }
                )
            else:
                row = hit.iloc[0]
                focus_rows.append(
                    {
                        "comparison": comparison,
                        "subclass": subclass,
                        "gene": gene,
                        "present_in_model": True,
                        "rank": row["rank"],
                        "mean_importance": row[
                            "mean_importance"
                        ],
                    }
                )
    focus = pd.DataFrame(focus_rows)
    return top, focus


def markdown_table(frame, columns=None, digits=4):
    selected = frame if columns is None else frame[columns]
    headers = [str(column) for column in selected.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]
    for _, row in selected.iterrows():
        values = []
        for value in row:
            if pd.isna(value):
                values.append("")
            elif isinstance(value, float):
                values.append(f"{value:.{digits}f}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def build_summary_markdown(
    milo,
    augur,
    limited,
    focus_importance,
):
    milo_view = milo.copy()
    milo_view["median_logFC"] = milo_view[
        "median_logFC"
    ].round(4)
    milo_view["median_subclass_fraction"] = milo_view[
        "median_subclass_fraction"
    ].round(4)

    augur_view = augur.copy()
    augur_view["AUC"] = (
        augur_view["auc"].round(4).astype(str)
        + " (SD "
        + augur_view["auc_sd"].round(4).astype(str)
        + ")"
    )
    augur_view = augur_view[
        [
            "comparison",
            "subclass",
            "reference_cells",
            "case_cells",
            "AUC",
        ]
    ]

    limited_view = limited.copy()
    for column in [
        "limited_auc",
        "full_auc",
        "auc_change_full_minus_limited",
    ]:
        limited_view[column] = limited_view[column].round(4)

    pap_present = focus_importance[
        (focus_importance["gene"] == "PAPPA2")
        & (focus_importance["present_in_model"])
    ]
    other_present = focus_importance[
        (focus_importance["gene"] != "PAPPA2")
        & (focus_importance["present_in_model"])
    ]
    text = f"""# PAPPA2 项目 TAL 亚类 Milo 与 Augur 汇总

日期：2026-09-27

## 一、本次完成

- TAL 层面 Milo 已完成 3 个比较：疾病合并、AKI、CKD；
- TAL 层面 Augur 已完成原始交接口径和四亚类全覆盖口径；
- 四亚类全覆盖 Augur 使用 4,446 个 TAL 细胞、2,500 个特征；
- 十个主线相关基因已保留在 TAL 专用表达矩阵中；
- 已完成供者层面的候选基因伪批量比较。

原始交接脚本使用通用 Augur 输入。该输入只包含 583 个 TAL 细胞，
其中 Reference 的 aTAL2 只有 9 个、frTAL 只有 19 个，均低于
min_cells = 20，因此原口径不能得到四个亚类 AUC。

本文档把四亚类全覆盖结果作为 TAL 层主分析，把原有 583 细胞结果
作为敏感性参考。两者不能混用。

## 二、TAL Milo 结果

正式显著性口径：SpatialFDR < 0.05。

{markdown_table(milo_view, ["comparison", "subclass", "total_annotated_nhoods", "significant_nhoods", "positive_nhoods", "negative_nhoods", "median_logFC", "median_subclass_fraction"])}

主要观察：

- 疾病合并时，四个亚类都有显著邻域；
- AKI 中 aTAL2 的显著邻域最多，且全部为正方向；
- CKD 中 frTAL 的显著邻域最多，且全部为正方向；
- C-TAL-A 在 AKI 中同时存在正向和负向邻域，不能写成整体扩张；
- aTAL2 在 CKD 中正负混合，和 AKI 的方向不一致。

## 三、TAL Augur 四亚类结果

子采样设置为 20 次、每次 30 个细胞、3 折交叉验证。
AUC 后的 SD 来自 60 个“子采样 x 折”的 ROC AUC 估计，只用于
描述波动，不是供者层面的置信区间。

{markdown_table(augur_view)}

主要观察：

- AKI 中最容易和 Reference 区分的是 C-TAL-A，AUC = 0.753；
- AKI 的四个亚类都有高于 0.5 的平均 AUC；
- CKD 中平均 AUC 最高的是 frTAL 0.590，其次是 C-TAL-B 0.573；
- CKD 的 C-TAL-A 平均 AUC 为 0.513，接近随机；
- CKD 的 aTAL2 平均 AUC 为 0.479，低于随机；
- 因此 CKD 的 TAL 状态变化不能写成四个亚类一致。

## 四、原始口径与全覆盖口径差异

{markdown_table(limited_view)}

差异最大的例子：

- 疾病合并 C-TAL-B：原口径 0.473，全覆盖 0.606；
- AKI C-TAL-A：原口径 0.614，全覆盖 0.753；
- CKD C-TAL-A：原口径 0.543，全覆盖 0.513。

这说明原 583 细胞输入受抽样子集影响明显，不适合作为最终 TAL
亚类结论。推荐后续使用 `KPMP_Augur_TAL_full_*` 文件。

## 五、Augur 特征重要性

`feature_perc = 0.5` 表示每个子采样只保留一部分特征进入模型。
十个主线候选基因中，只有 PAPPA2 一致进入模型重要性表。

- PAPPA2 记录数：{len(pap_present)}；
- 其他候选基因进入模型的记录数：{len(other_present)}；
- 因此 Augur 不能直接验证 JAG1-NOTCH3-HES1 或 EPHB2-STAT6；
- 将候选基因保留在输入矩阵，不等于它们一定进入最终模型。

## 六、候选基因供者层面表达结果

详细文件：

`KPMP_TAL_focus_gene_pseudobulk.csv`

分析方法：

- 先对每个细胞做 CP10K 归一化和 log1p；
- 再按供者和亚类取平均；
- 疾病组与 Reference 使用双侧 Mann-Whitney 检验；
- FDR 对 120 个“基因 x 亚类 x 比较”检验统一校正。

主要结果：

- PAPPA2 在四个亚类、AKI、CKD 和疾病合并比较中均显著升高；
- CXCL12 在 CKD 的 C-TAL-A、C-TAL-B、frTAL 中显著升高；
- NOTCH3 只在 CKD 的 frTAL 中显著升高；
- HES1 在 CKD 的 aTAL2 中显著降低；
- EPHB2 在 AKI 的 C-TAL-B 和 aTAL2 中显著升高；
- STAT6 没有显著变化；
- JAG1 和 HES4 没有显著变化。

解释边界：

- PAPPA2 的亚类结果最稳定，四个亚类都呈升高方向；
- JAG1-NOTCH3-HES1 没有被 TAL 亚类表达结果整体复现；
- EPHB2 在 AKI 中有信号，但 STAT6 未变化，不能称为 EPHB2-STAT6 轴成立；
- 这些是表达关联，不是通路激活的直接测量；
- 伪批量检验使用供者作为观测单位，比直接检验单个细胞更保守。

## 七、当前结论

可以写的结论：

- TAL 疾病相关邻域变化具有明显亚类差异；
- AKI 的转录可分性以 C-TAL-A 最明显；
- CKD 的信号主要出现在 frTAL 和 C-TAL-B；
- PAPPA2 是四个 TAL 亚类中最一致的疾病升高基因；
- CXCL12 在 CKD 的多个 TAL 亚类中升高。

不能写的结论：

- 不能写四个 TAL 亚类反应一致；
- 不能写 JAG1-NOTCH3-HES1 已在 TAL 亚类中复现；
- 不能写 EPHB2-STAT6 已成立；
- 不能把 AUC 写成因果效应；
- 不能把 Milo 邻域正方向直接解释为细胞类型扩增。

## 八、建议下一步

1. 以 `KPMP_Augur_TAL_full_*` 为 TAL Augur 正式结果；
2. 检查 Milo 邻域纯度较低的少数邻域，必要时做敏感性过滤；
3. 若要连接表达和 Milo，需保存每个邻域的细胞成员，再计算邻域
   平均表达或候选基因评分；
4. 在 TAL 结果稳定前，不重跑 SigXTalk，也不增加新工具。
"""
    return text


def build_handoff_v7(
    milo,
    augur,
    limited,
    focus_importance,
):
    text = f"""# PAPPA2 项目交接文档 v7

日期：2026-09-27

## 一、新会话第一句话

请读取：

```text
C:\\Users\\Elsa\\Documents\\Codex\\2026-09-25\\9-2\\outputs\\PAPPA2_project_handoff_v7.md
```

然后继续 PAPPA2 项目的 TAL 亚类敏感性分析。

## 二、当前状态

- TAL Milo 已完成；
- TAL Augur 原始口径已完成；
- TAL 四亚类全覆盖 Augur 已完成；
- 十个主线候选基因伪批量比较已完成；
- 当前没有后台任务；
- 不需要重新下载 Visium、NicheNet 或 SCENIC 数据。

## 三、TAL Milo 正式结果

- 疾病合并：219 个显著邻域；
- AKI 对 Reference：102 个显著邻域；
- CKD 对 Reference：148 个显著邻域；
- 四个亚类均出现显著邻域；
- AKI 的 aTAL2 和 CKD 的 frTAL 以正方向为主；
- C-TAL-A 和 aTAL2 存在明显的比较间方向差异。

文件：

```text
outputs\\KPMP_milo_TAL_all_disease_nhood_results.csv
outputs\\KPMP_milo_TAL_AKI_vs_Reference_nhood_results.csv
outputs\\KPMP_milo_TAL_CKD_vs_Reference_nhood_results.csv
outputs\\KPMP_milo_TAL_subclass_summary.csv
```

## 四、TAL Augur 正式结果

最终 TAL 层结果应使用：

```text
outputs\\KPMP_Augur_TAL_full_all_AUC.csv
outputs\\KPMP_Augur_TAL_full_feature_importance_mapped.csv
outputs\\KPMP_Augur_TAL_full_AUC_summary.csv
```

输入：

```text
D:\\CodexData\\RAMP3\\augur_tal_full_input\\KPMP_augur_tal_full_input.h5ad
```

规模：

- 4,446 个 TAL 细胞；
- 2,500 个特征；
- 四个目标亚类；
- Reference 的四个亚类均大于 20 个细胞。

主要 AUC：

{markdown_table(augur, ["comparison", "subclass", "auc", "auc_sd"])}

判断：

- AKI 的 C-TAL-A 可分性最强；
- CKD 的 frTAL 和 C-TAL-B 略高于随机；
- CKD 的 aTAL2 平均 AUC 低于 0.5；
- 不能把四个亚类合并成一个方向性结论。

原有 583 细胞口径只作为敏感性参考，不应作为最终结果。

## 五、候选基因

详细结果：

```text
outputs\\KPMP_TAL_focus_gene_pseudobulk.csv
```

当前判断：

- PAPPA2 在四个亚类中稳定升高；
- CXCL12 在 CKD 的多个亚类中升高；
- NOTCH3 只在 CKD 的 frTAL 中升高；
- HES1 在 CKD 的 aTAL2 中降低；
- HES4 和 JAG1 没有显著变化；
- EPHB2 在部分 AKI 亚类升高，但 STAT6 没有显著变化。

不能说 JAG1-NOTCH3-HES1 或 EPHB2-STAT6 已在 TAL 层级复现。

## 六、关键限制

- Milo 使用已有 UMAP，没有重新计算 PCA；
- Milo 是在 TAL 亚类内部重新建图，只能解释亚类邻域状态；
- Augur AUC 不是因果效应；
- Augur 子采样结果共享同一批细胞，SD 不是独立置信区间；
- 伪批量比较是表达关联，不是通路活性；
- 没有实验验证。

## 七、下一步顺序

第一步：

- 检查 Milo 邻域纯度；
- 对纯度较低的邻域做敏感性分析；
- 不重跑全体细胞 Milo。

第二步：

- 保存 TAL Milo 每个邻域的细胞成员；
- 计算候选基因或通路评分在每个邻域中的平均表达；
- 判断 PAPPA2、CXCL12、NOTCH3、HES1 和 EPHB2 与邻域方向的关联。

第三步：

- 只有在 TAL 邻域表达和亚类结果稳定后，才决定是否重跑 SigXTalk；
- 不增加新的重型工具。

## 八、关键文件

```text
outputs\\PAPPA2_TAL_milo_augur_summary.md
outputs\\KPMP_milo_TAL_subclass_summary.csv
outputs\\KPMP_Augur_TAL_full_AUC_summary.csv
outputs\\KPMP_Augur_TAL_full_top_features.csv
outputs\\KPMP_Augur_TAL_full_focus_gene_importance.csv
outputs\\KPMP_TAL_focus_gene_pseudobulk.csv
scripts\\prepare_augur_input.py
scripts\\run_pyaugur_tal_kpmp.py
scripts\\run_milo_tal_kpmp.R
scripts\\analyze_tal_focus_genes.py
scripts\\summarize_tal_milo_augur.py
```

## 九、不要重复

- 不要重新下载数据；
- 不要重新跑全体细胞 Milo 和 Augur；
- 不要把原始 583 细胞 Augur 当作最终 TAL 层级结果；
- 不要把表达关联写成因果；
- 不要用 SigXTalk 覆盖 NicheNet 和 SCENIC。
"""
    return text


def main():
    milo = summarize_milo()
    augur = summarize_augur()
    limited = summarize_limited_vs_full()
    top_features, focus_importance = summarize_importance()

    milo.to_csv(
        OUTPUT_DIR / "KPMP_milo_TAL_subclass_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    augur.to_csv(
        OUTPUT_DIR / "KPMP_Augur_TAL_full_AUC_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    limited.to_csv(
        OUTPUT_DIR / "KPMP_Augur_TAL_limited_vs_full_AUC.csv",
        index=False,
        encoding="utf-8-sig",
    )
    top_features.to_csv(
        OUTPUT_DIR / "KPMP_Augur_TAL_full_top_features.csv",
        index=False,
        encoding="utf-8-sig",
    )
    focus_importance.to_csv(
        OUTPUT_DIR
        / "KPMP_Augur_TAL_full_focus_gene_importance.csv",
        index=False,
        encoding="utf-8-sig",
    )

    summary_path = (
        OUTPUT_DIR / "PAPPA2_TAL_milo_augur_summary.md"
    )
    summary_path.write_text(
        build_summary_markdown(
            milo,
            augur,
            limited,
            focus_importance,
        ),
        encoding="utf-8",
    )
    handoff_path = (
        OUTPUT_DIR / "PAPPA2_project_handoff_v7.md"
    )
    handoff_path.write_text(
        build_handoff_v7(
            milo,
            augur,
            limited,
            focus_importance,
        ),
        encoding="utf-8",
    )
    print("Saved:", summary_path)
    print("Saved:", handoff_path)


if __name__ == "__main__":
    main()
