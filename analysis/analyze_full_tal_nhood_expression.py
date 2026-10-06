import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from project_paths import OUTPUT_ROOT, WORK_ROOT

INPUT_EXPRESSION = (
    OUTPUT_ROOT / "KPMP_full_TAL_focus_expression.csv.gz"
)
MILO_PREFIX = "KPMP_milo_full_TAL_final"

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

EXPRESSION_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_nhood_focus_expression.csv"
)
ASSOCIATION_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_nhood_focus_logFC_associations.csv"
)
VALIDATION_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_nhood_focus_validation.csv"
)


def load_cell_expression():
    if not INPUT_EXPRESSION.exists():
        raise FileNotFoundError(INPUT_EXPRESSION)
    expression = pd.read_csv(
        INPUT_EXPRESSION,
        dtype={"cell_id": str, "donor_id": str},
    )
    if not expression["cell_id"].is_unique:
        raise ValueError("表达表中的 cell_id 不唯一。")
    expression_columns = {
        gene: f"{gene}_lognorm" for gene in FOCUS_GENES
    }
    missing_columns = [
        column
        for column in expression_columns.values()
        if column not in expression.columns
    ]
    if missing_columns:
        raise ValueError(
            f"表达表缺少候选基因列: {missing_columns}"
        )

    result = expression.set_index("cell_id")[
        list(expression_columns.values())
    ].rename(columns={
        column: gene
        for gene, column in expression_columns.items()
    })
    result.index.name = "cell_id"
    return result


def read_membership(comparison):
    path = OUTPUT_ROOT / (
        f"{MILO_PREFIX}_{comparison}_nhood_membership.csv.gz"
    )
    membership = pd.read_csv(
        path,
        dtype={"cell_id": str, "donor_id": str},
    )
    required = {
        "cell_id",
        "subclass",
        "donor_id",
        "condition",
        "Nhood",
        "weight",
    }
    missing = required - set(membership.columns)
    if missing:
        raise ValueError(
            f"{path.name} 缺少列: {sorted(missing)}"
        )
    return membership


def read_nhood_results(comparison):
    path = OUTPUT_ROOT / (
        f"{MILO_PREFIX}_{comparison}_nhood_results.csv"
    )
    results = pd.read_csv(path)
    required = {
        "Nhood",
        "logFC",
        "SpatialFDR",
        "subclass",
        "subclass_fraction",
        "nhood_purity",
    }
    missing = required - set(results.columns)
    if missing:
        raise ValueError(
            f"{path.name} 缺少列: {sorted(missing)}"
        )
    if not results["Nhood"].is_unique:
        raise ValueError(f"{path.name} 中 Nhood 不唯一。")
    return results


def detect_nhood_offset(membership, results):
    """兼容 Milo 成员表和结果表可能不同的 Nhood 起始编号。"""
    membership_ids = set(
        membership["Nhood"].astype(int).tolist()
    )
    result_ids = set(results["Nhood"].astype(int).tolist())
    matches = []
    for offset in range(-2, 3):
        shifted_ids = {
            nhood + offset for nhood in membership_ids
        }
        exact_match = shifted_ids == result_ids
        matches.append((offset, exact_match))
    exact_offsets = [
        offset for offset, exact_match in matches
        if exact_match
    ]
    if len(exact_offsets) != 1:
        raise ValueError(
            "无法唯一确定成员表 Nhood 到结果表 Nhood 的偏移量，"
            f"候选偏移: {exact_offsets}"
        )
    return exact_offsets[0]


def build_nhood_expression(
    comparison,
    membership,
    results,
    cell_expression,
):
    nhood_offset = detect_nhood_offset(membership, results)
    membership = membership.copy()
    membership["result_nhood"] = (
        membership["Nhood"].astype(int) + nhood_offset
    )

    if not np.allclose(
        membership["weight"].to_numpy(dtype=float),
        1.0,
    ):
        raise ValueError(
            "当前脚本只支持 weight 全为 1 的 Milo 成员表。"
        )

    merged = membership.merge(
        cell_expression.reset_index(),
        on="cell_id",
        how="left",
        validate="many_to_one",
        indicator=True,
    )
    missing_cells = merged["_merge"] != "both"
    if missing_cells.any():
        examples = merged.loc[
            missing_cells, "cell_id"
        ].drop_duplicates().head(10).tolist()
        raise ValueError(
            f"{comparison} 有成员细胞不在表达矩阵中，例如: "
            f"{examples}"
        )
    merged = merged.drop(columns="_merge")

    base = (
        merged.groupby("result_nhood", observed=True)
        .agg(
            n_cells=("cell_id", "size"),
            n_donors=("donor_id", "nunique"),
        )
    )
    all_means = (
        merged.groupby("result_nhood", observed=True)[
            FOCUS_GENES
        ]
        .mean()
        .add_suffix("_mean_all")
    )

    merged["condition_group"] = np.where(
        merged["condition"].astype(str) == "Reference",
        "Reference",
        "Disease",
    )
    condition_means = (
        merged.groupby(
            ["result_nhood", "condition_group"],
            observed=True,
        )[FOCUS_GENES]
        .mean()
        .unstack("condition_group")
    )
    condition_means.columns = [
        f"{gene}_mean_{condition.lower()}"
        for gene, condition in condition_means.columns
    ]

    condition_counts = (
        merged.groupby(
            ["result_nhood", "condition_group"],
            observed=True,
        )
        .size()
        .unstack("condition_group")
        .rename(
            columns={
                "Reference": "n_reference_cells",
                "Disease": "n_disease_cells",
            }
        )
    )

    expression = (
        base.join(all_means)
        .join(condition_means)
        .join(condition_counts)
        .reset_index()
    )
    expression = expression.merge(
        results,
        left_on="result_nhood",
        right_on="Nhood",
        how="inner",
        validate="one_to_one",
    )

    classified_membership = membership[
        ["cell_id", "result_nhood"]
    ].merge(
        results[["Nhood", "subclass"]],
        left_on="result_nhood",
        right_on="Nhood",
        how="inner",
        validate="many_to_one",
    )
    unique_cells_by_subclass = (
        classified_membership.groupby(
            "subclass",
            observed=True,
        )["cell_id"]
        .nunique()
        .to_dict()
    )
    expression["unique_cells_in_subclass"] = (
        expression["subclass"].map(
            unique_cells_by_subclass
        )
    )

    if len(expression) != len(results):
        raise ValueError(
            f"{comparison} 合并后邻域数不一致: "
            f"{len(expression)} vs {len(results)}"
        )

    long_rows = []
    for gene in FOCUS_GENES:
        columns = [
            "comparison",
            "result_nhood",
            "subclass",
            "n_cells",
            "n_donors",
            "n_reference_cells",
            "n_disease_cells",
            "unique_cells_in_subclass",
            "logFC",
            "SpatialFDR",
            "subclass_fraction",
            "nhood_purity",
        ]
        gene_table = expression[columns].copy()
        gene_table["gene"] = gene
        gene_table["mean_expression_all"] = expression[
            f"{gene}_mean_all"
        ]
        gene_table["mean_expression_reference"] = expression[
            f"{gene}_mean_reference"
        ]
        gene_table["mean_expression_disease"] = expression[
            f"{gene}_mean_disease"
        ]
        long_rows.append(gene_table)

    return pd.concat(long_rows, ignore_index=True), nhood_offset


def correlate_one_group(group):
    x = group["mean_expression_all"].to_numpy(dtype=float)
    y = group["logFC"].to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]

    if (
        len(x) < 3
        or np.ptp(x) == 0
        or np.ptp(y) == 0
    ):
        pearson_r = np.nan
        pearson_p = np.nan
        spearman_rho = np.nan
        spearman_p = np.nan
    else:
        pearson_r, pearson_p = pearsonr(x, y)
        spearman_rho, spearman_p = spearmanr(x, y)

    significant = group["SpatialFDR"] < 0.05
    return pd.Series(
        {
            "n_nhoods": len(group),
            "n_unique_cells_in_subclass": int(
                group[
                    "unique_cells_in_subclass"
                ].iloc[0]
            ),
            "median_nhood_cells": group[
                "n_cells"
            ].median(),
            "significant_nhoods": int(significant.sum()),
            "positive_significant_nhoods": int(
                (
                    significant
                    & (group["logFC"] > 0)
                ).sum()
            ),
            "negative_significant_nhoods": int(
                (
                    significant
                    & (group["logFC"] < 0)
                ).sum()
            ),
            "median_logFC": group["logFC"].median(),
            "median_expression": group[
                "mean_expression_all"
            ].median(),
            "pearson_r": pearson_r,
            "pearson_p_naive": pearson_p,
            "spearman_rho": spearman_rho,
            "spearman_p_naive": spearman_p,
        }
    )


def build_associations(nhood_expression):
    rows = []
    group_columns = [
        "comparison",
        "subclass",
        "gene",
    ]
    for keys, group in nhood_expression.groupby(
        group_columns,
        observed=True,
    ):
        comparison, subclass, gene = keys
        values = correlate_one_group(group)
        rows.append(
            {
                "comparison": comparison,
                "subclass": subclass,
                "gene": gene,
                **values.to_dict(),
            }
        )
    associations = pd.DataFrame(rows)
    associations["abs_spearman_rho"] = associations[
        "spearman_rho"
    ].abs()
    return associations.sort_values(
        ["comparison", "subclass", "gene"]
    )


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    cell_expression = load_cell_expression()

    expression_tables = []
    validation_rows = []
    for comparison in COMPARISONS:
        membership = read_membership(comparison)
        results = read_nhood_results(comparison)
        nhood_expression, offset = build_nhood_expression(
            comparison,
            membership,
            results,
            cell_expression,
        )
        expression_tables.append(nhood_expression)
        membership_cells = membership[
            "cell_id"
        ].nunique()
        missing_membership_cells = (
            set(membership["cell_id"])
            - set(cell_expression.index)
        )
        validation_rows.append(
            {
                "comparison": comparison,
                "input_cells": len(cell_expression),
                "membership_rows": len(membership),
                "membership_nhoods": membership[
                    "Nhood"
                ].nunique(),
                "membership_cells": membership_cells,
                "memberships_per_cell_mean": (
                    len(membership) / membership_cells
                ),
                "result_nhoods": len(results),
                "mapped_nhoods": nhood_expression[
                    "result_nhood"
                ].nunique(),
                "nhood_offset": offset,
                "missing_membership_cells": len(
                    missing_membership_cells
                ),
                "all_membership_cells_found": (
                    len(missing_membership_cells) == 0
                ),
            }
        )

    nhood_expression = pd.concat(
        expression_tables,
        ignore_index=True,
    )
    associations = build_associations(nhood_expression)
    validation = pd.DataFrame(validation_rows)

    nhood_expression.to_csv(
        EXPRESSION_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    associations.to_csv(
        ASSOCIATION_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    validation.to_csv(
        VALIDATION_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )

    print("Saved:", EXPRESSION_OUTPUT)
    print("Saved:", ASSOCIATION_OUTPUT)
    print("Saved:", VALIDATION_OUTPUT)
    print("\nValidation:")
    print(validation.to_string(index=False))
    print("\nTop absolute Spearman correlations:")
    print(
        associations.sort_values(
            "abs_spearman_rho",
            ascending=False,
        )
        .head(15)[
            [
                "comparison",
                "subclass",
                "gene",
                "n_nhoods",
                "spearman_rho",
                "spearman_p_naive",
            ]
        ]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
