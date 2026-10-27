"""Four-unit independent validation prototype for GSE254185.

The four descriptive units are:

1. donor pseudobulk disease-versus-control expression
2. neighborhood abundance
3. neighborhood-expression association
4. within-donor cell-rank association

The neighborhood implementation uses public UMAP coordinates and a kNN graph.
It is a transparent approximation of Milo, not a claim to reproduce Milo
exactly.
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.neighbors import NearestNeighbors

from validate_gse254185 import (
    BASE,
    GENES,
    STATE_ORDER,
    add_injury_exposure,
    bh_adjust,
    donor_pseudobulk,
    load_expression,
    load_metadata,
)


OUT_DIR = BASE / "gse254185_four_units"


def add_neighborhood_features(expression: pd.DataFrame, k: int = 30) -> pd.DataFrame:
    expression = expression.copy()
    expression["nhood_logfc"] = np.nan
    for gene in GENES:
        if gene in expression.columns:
            expression[f"nhood_mean_{gene}"] = np.nan

    for library, block in expression.groupby("Library", sort=False):
        if len(block) < 10:
            continue
        coords = block[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)
        valid = np.isfinite(coords).all(axis=1)
        if valid.sum() < 10:
            continue

        n_neighbors = min(k + 1, valid.sum())
        model = NearestNeighbors(n_neighbors=n_neighbors, algorithm="auto")
        model.fit(coords[valid])
        _, neighbors = model.kneighbors(coords[valid])
        if n_neighbors > 1:
            neighbors = neighbors[:, 1:]

        valid_index = block.index[valid]
        group = block["group"].to_numpy()
        is_uuo = group == "UUO"
        total_uuo = max(int(is_uuo.sum()), 1)
        total_normal = max(int((~is_uuo).sum()), 1)

        n_uuo = is_uuo[neighbors].sum(axis=1)
        n_normal = (~is_uuo)[neighbors].sum(axis=1)
        logfc = np.log2(
            ((n_uuo + 0.5) / (total_uuo + 1.0))
            / ((n_normal + 0.5) / (total_normal + 1.0))
        )
        expression.loc[valid_index, "nhood_logfc"] = logfc

        for gene in GENES:
            column = f"nhood_mean_{gene}"
            if column not in expression.columns:
                continue
            values = block[gene].to_numpy(dtype=float)
            expression.loc[valid_index, column] = np.nanmean(
                values[neighbors], axis=1
            )
    return expression


def neighborhood_abundance(expression: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for state in STATE_ORDER:
        block = expression[expression["state"].eq(state)]
        values = block["nhood_logfc"].dropna()
        if len(values) < 20:
            continue
        rows.append(
            {
                "state": state,
                "n_centers": len(values),
                "mean_logfc": values.mean(),
                "median_logfc": values.median(),
                "fraction_positive": (values > 0).mean(),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["direction"] = np.where(
            result["mean_logfc"] > 0, "positive", "negative"
        )
    return result


def neighborhood_expression(expression: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for state in STATE_ORDER:
        state_block = expression[expression["state"].eq(state)]
        for gene in GENES:
            column = f"nhood_mean_{gene}"
            if column not in state_block.columns:
                continue
            valid = state_block[[column, "nhood_logfc"]].dropna()
            if (
                len(valid) < 20
                or valid[column].nunique() < 2
                or valid["nhood_logfc"].nunique() < 2
            ):
                continue
            rho, p_value = spearmanr(valid[column], valid["nhood_logfc"])
            rows.append(
                {
                    "state": state,
                    "gene": gene,
                    "n_centers": len(valid),
                    "rho": rho,
                    "p_value": p_value,
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = bh_adjust(result["p_value"])
        result["direction"] = np.where(result["rho"] > 0, "positive", "negative")
    return result


def within_donor_rank(expression: pd.DataFrame) -> pd.DataFrame:
    donor_rows = []
    for state in STATE_ORDER:
        state_block = expression[expression["state"].eq(state)]
        for donor, block in state_block.groupby("donor", sort=False):
            for gene in GENES:
                valid = block[[gene, "nhood_logfc"]].dropna()
                if (
                    len(valid) < 20
                    or valid[gene].nunique() < 2
                    or valid["nhood_logfc"].nunique() < 3
                ):
                    continue
                rho, p_value = spearmanr(valid[gene], valid["nhood_logfc"])
                donor_rows.append(
                    {
                        "state": state,
                        "gene": gene,
                        "donor": donor,
                        "group": block["group"].iloc[0],
                        "n_cells": len(valid),
                        "rho": rho,
                        "p_value": p_value,
                    }
                )

    donor_result = pd.DataFrame(donor_rows)
    if donor_result.empty:
        return donor_result

    summary_rows = []
    for (state, gene), block in donor_result.groupby(["state", "gene"], sort=False):
        uuo = block.loc[block["group"].eq("UUO"), "rho"]
        normal = block.loc[block["group"].eq("Normal"), "rho"]
        if len(uuo) < 2 or len(normal) < 2:
            continue
        p_value = mannwhitneyu(uuo, normal, alternative="two-sided").pvalue
        summary_rows.append(
            {
                "state": state,
                "gene": gene,
                "n_uuo_donors": len(uuo),
                "n_normal_donors": len(normal),
                "uuo_mean_rho": uuo.mean(),
                "normal_mean_rho": normal.mean(),
                "mean_rho_difference": uuo.mean() - normal.mean(),
                "p_value": p_value,
            }
        )
    result = pd.DataFrame(summary_rows)
    if not result.empty:
        result["fdr"] = bh_adjust(result["p_value"])
        result["direction"] = np.where(
            result["mean_rho_difference"] > 0, "positive", "negative"
        )
    return result


def build_comparison(
    pseudobulk: pd.DataFrame,
    abundance: pd.DataFrame,
    neighborhood: pd.DataFrame,
    within: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for state in STATE_ORDER:
        for gene in GENES:
            p = pseudobulk[
                pseudobulk["state"].eq(state) & pseudobulk["gene"].eq(gene)
            ]
            a = abundance[abundance["state"].eq(state)]
            n = neighborhood[
                neighborhood["state"].eq(state) & neighborhood["gene"].eq(gene)
            ]
            w = within[within["state"].eq(state) & within["gene"].eq(gene)]
            if p.empty or a.empty or n.empty or w.empty:
                continue
            rows.append(
                {
                    "state": state,
                    "gene": gene,
                    "pseudobulk_direction": p["direction"].iloc[0],
                    "abundance_direction": a["direction"].iloc[0],
                    "neighborhood_expression_direction": n["direction"].iloc[0],
                    "within_donor_direction": w["direction"].iloc[0],
                    "pseudobulk_difference": p["mean_difference"].iloc[0],
                    "abundance_logfc": a["mean_logfc"].iloc[0],
                    "neighborhood_expression_rho": n["rho"].iloc[0],
                    "within_donor_rho_difference": w["mean_rho_difference"].iloc[0],
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        direction_columns = [
            "pseudobulk_direction",
            "abundance_direction",
            "neighborhood_expression_direction",
            "within_donor_direction",
        ]
        result["n_positive"] = (
            result[direction_columns].eq("positive").sum(axis=1)
        )
        result["n_negative"] = (
            result[direction_columns].eq("negative").sum(axis=1)
        )
        result["discordant_four_unit"] = (
            result["n_positive"].gt(0) & result["n_negative"].gt(0)
        )
    return result


def exact_group_permutation(values: pd.Series, labels: pd.Series) -> tuple[float, float, int]:
    values = values.to_numpy(dtype=float)
    labels = labels.to_numpy()
    n_uuo = int((labels == "UUO").sum())
    if len(values) < 4 or n_uuo < 2 or n_uuo == len(values):
        return np.nan, np.nan, 0
    observed = values[labels == "UUO"].mean() - values[labels == "Normal"].mean()
    null = []
    for chosen in combinations(range(len(values)), n_uuo):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(chosen)] = True
        null.append(values[mask].mean() - values[~mask].mean())
    null = np.asarray(null)
    p_value = (np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1)
    return observed, p_value, len(null)


def donor_level_means(expression: pd.DataFrame, state: str, gene: str) -> pd.DataFrame:
    block = expression[expression["state"].eq(state)]
    means = (
        block.groupby(["donor", "group"], observed=True)[gene]
        .agg(["mean", "size"])
        .reset_index()
        .rename(columns={"mean": "mean_expression", "size": "n_cells"})
    )
    return means[means["n_cells"] >= 5]


def exact_pseudobulk_tests(expression: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for state in STATE_ORDER:
        means = donor_level_means(expression, state, "PAPPA2")
        if means["group"].nunique() < 2:
            continue
        observed, p_value, n_perm = exact_group_permutation(
            means["mean_expression"], means["group"]
        )
        rows.append(
            {
                "state": state,
                "observed_difference": observed,
                "p_exact": p_value,
                "n_permutations": n_perm,
                "n_uuo_donors": int(means["group"].eq("UUO").sum()),
                "n_normal_donors": int(means["group"].eq("Normal").sum()),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr_exact"] = bh_adjust(result["p_exact"])
    return result


def leave_one_donor_out(expression: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for state in STATE_ORDER:
        means = donor_level_means(expression, state, "PAPPA2")
        for donor in means["donor"]:
            kept = means[~means["donor"].eq(donor)]
            if kept["group"].nunique() < 2:
                continue
            observed = (
                kept.loc[kept["group"].eq("UUO"), "mean_expression"].mean()
                - kept.loc[kept["group"].eq("Normal"), "mean_expression"].mean()
            )
            rows.append(
                {
                    "state": state,
                    "removed_donor": donor,
                    "difference": observed,
                }
            )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return (
        result.groupby("state", observed=True)["difference"]
        .agg(["min", "max", "mean"])
        .reset_index()
        .rename(
            columns={
                "min": "leave_one_out_min",
                "max": "leave_one_out_max",
                "mean": "leave_one_out_mean",
            }
        )
    )


def exact_within_tests(expression: pd.DataFrame, exposure: str) -> pd.DataFrame:
    rows = []
    for state in STATE_ORDER:
        block = expression[expression["state"].eq(state)]
        donor_values = []
        for donor, donor_block in block.groupby("donor", sort=False):
            valid = donor_block[["PAPPA2", exposure]].dropna()
            if (
                len(valid) < 20
                or valid[exposure].nunique() < 3
                or valid["PAPPA2"].nunique() < 2
            ):
                continue
            rho, _ = spearmanr(valid["PAPPA2"], valid[exposure])
            donor_values.append(
                {"donor": donor, "group": donor_block["group"].iloc[0], "rho": rho}
            )
        if not donor_values:
            continue
        donor_table = pd.DataFrame(donor_values)
        if donor_table["group"].nunique() < 2:
            continue
        observed, p_value, n_perm = exact_group_permutation(
            donor_table["rho"], donor_table["group"]
        )
        rows.append(
            {
                "exposure": exposure,
                "state": state,
                "observed_rho_difference": observed,
                "p_exact": p_value,
                "n_permutations": n_perm,
                "n_uuo_donors": int(donor_table["group"].eq("UUO").sum()),
                "n_normal_donors": int(donor_table["group"].eq("Normal").sum()),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr_exact"] = bh_adjust(result["p_exact"])
    return result


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta = load_metadata()
    expression = add_injury_exposure(load_expression(meta))
    expression = add_neighborhood_features(expression)

    pseudobulk = donor_pseudobulk(expression)
    abundance = neighborhood_abundance(expression)
    neighborhood = neighborhood_expression(expression)
    within = within_donor_rank(expression)
    comparison = build_comparison(pseudobulk, abundance, neighborhood, within)
    pseudobulk_exact = exact_pseudobulk_tests(expression)
    pseudobulk_leave_one_out = leave_one_donor_out(expression)
    within_injury_exact = exact_within_tests(expression, "injury_exposure")
    within_nhood_exact = exact_within_tests(expression, "nhood_logfc")

    expression.to_pickle(OUT_DIR / "GSE254185_four_unit_cells.pkl")
    pseudobulk.to_csv(OUT_DIR / "GSE254185_pseudobulk.csv", index=False)
    abundance.to_csv(OUT_DIR / "GSE254185_neighborhood_abundance.csv", index=False)
    neighborhood.to_csv(
        OUT_DIR / "GSE254185_neighborhood_expression.csv", index=False
    )
    within.to_csv(OUT_DIR / "GSE254185_within_donor_rank.csv", index=False)
    comparison.to_csv(OUT_DIR / "GSE254185_four_unit_comparison.csv", index=False)
    pseudobulk_exact.to_csv(OUT_DIR / "GSE254185_pseudobulk_exact.csv", index=False)
    pseudobulk_leave_one_out.to_csv(
        OUT_DIR / "GSE254185_pseudobulk_leave_one_out.csv", index=False
    )
    within_injury_exact.to_csv(
        OUT_DIR / "GSE254185_within_injury_exact.csv", index=False
    )
    within_nhood_exact.to_csv(
        OUT_DIR / "GSE254185_within_nhood_exact.csv", index=False
    )

    print("Neighborhood abundance:")
    print(abundance.to_string(index=False))
    print("\nPAPPA2 four-unit comparison:")
    print(
        comparison[comparison["gene"].eq("PAPPA2")].to_string(index=False)
    )
    print("\nPAPPA2 exact donor-pseudobulk tests:")
    print(pseudobulk_exact.to_string(index=False))
    print("\nPAPPA2 exact within-donor tests (injury exposure):")
    print(within_injury_exact.to_string(index=False))
    print("\nPAPPA2 exact within-donor tests (neighborhood logFC):")
    print(within_nhood_exact.to_string(index=False))


if __name__ == "__main__":
    main()
