"""Second independent validation using the public GSE211785 H5AD object.

The non-KPMP Susztak object is read directly from the local D: copy. The
analysis is restricted to RNA TAL cells and uses disease versus control
donors, with exact or Monte Carlo donor-label permutations.
"""

from __future__ import annotations

import os
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors

from validate_gse254185 import (
    GENES,
    add_injury_exposure,
)
from validate_gse254185_four_units import add_neighborhood_features


H5AD = Path(
    os.environ.get(
        "GSE211785_H5AD",
        str(Path(__file__).resolve().parent / "GSE211785_PreSCVI.h5ad"),
    )
)
OUT_DIR = Path(
    os.environ.get(
        "PAPPA2_VALIDATION_OUTPUT",
        str(Path(__file__).resolve().parent / "gse211785_validation"),
    )
)
STATES = ["C_TAL", "M_TAL"]


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return pd.Series(out, index=values.index)


def load_expression() -> pd.DataFrame:
    adata = ad.read_h5ad(H5AD, backed="r")
    obs = adata.obs
    mask = (
        obs["tech"].isin(["SC_RNA", "SN_RNA"])
        & obs["Cluster_Idents"].isin(STATES)
        & obs["sample"].notna()
    )
    subset = adata[mask, GENES].to_memory()
    obs = subset.obs.copy()
    counts = subset.layers["counts"].toarray()
    totals = obs["nCount_RNA"].to_numpy(dtype=float)
    safe_totals = np.where(totals > 0, totals, np.nan)
    normalized = np.log1p(counts / safe_totals[:, None] * 10000.0)

    frame = pd.DataFrame(normalized, columns=GENES, index=obs.index)
    frame["state"] = obs["Cluster_Idents"].astype(str).to_numpy()
    frame["donor"] = obs["sample"].astype(str).to_numpy()
    frame["group_name"] = obs["group"].astype(str).to_numpy()
    frame["group"] = np.where(
        frame["group_name"].eq("Control"), "Normal", "UUO"
    )
    frame["Library"] = obs["tech"].astype(str).to_numpy()
    frame["sex"] = obs["sex"].astype(str).to_numpy()
    frame["age"] = pd.to_numeric(obs["age"], errors="coerce").to_numpy()
    umap = adata.obsm["X_umap"][mask.to_numpy()]
    frame["UMAP_1"] = umap[:, 0]
    frame["UMAP_2"] = umap[:, 1]
    frame["nCount_RNA"] = totals
    return frame


def permutation_test(
    values: np.ndarray,
    labels: np.ndarray,
    n_permutations: int = 20000,
    seed: int = 20261028,
) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    n_case = int((labels == "UUO").sum())
    n_control = int((labels == "Normal").sum())
    if n_case < 2 or n_control < 2:
        return np.nan, np.nan, 0
    observed = values[labels == "UUO"].mean() - values[labels == "Normal"].mean()
    total = len(values)
    combination_count = 1
    # Use exact enumeration only for small donor sets.
    for i in range(1, n_case + 1):
        combination_count = combination_count * (total - i + 1) // i
        if combination_count > n_permutations:
            break
    rng = np.random.default_rng(seed)
    if combination_count <= n_permutations:
        choices = list(combinations(range(total), n_case))
    else:
        choices = [
            tuple(rng.choice(total, n_case, replace=False))
            for _ in range(n_permutations)
        ]
    null = np.empty(len(choices), dtype=float)
    for i, chosen in enumerate(choices):
        mask = np.zeros(total, dtype=bool)
        mask[list(chosen)] = True
        null[i] = values[mask].mean() - values[~mask].mean()
    p_value = (np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1)
    return observed, p_value, len(null)


def donor_pseudobulk(expression: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for state in STATES:
        block = expression[expression["state"].eq(state)]
        donor_means = (
            block.groupby(["donor", "group", "group_name"], observed=True)["PAPPA2"]
            .agg(["mean", "size"])
            .reset_index()
            .rename(columns={"mean": "mean_expression", "size": "n_cells"})
        )
        donor_means = donor_means[donor_means["n_cells"] >= 5]
        for comparison, subset in [
            ("Disease_vs_Control", donor_means),
            (
                "DKD_vs_Control",
                donor_means[donor_means["group_name"].isin(["Control", "DKD"])],
            ),
            (
                "HKD_vs_Control",
                donor_means[donor_means["group_name"].isin(["Control", "HKD"])],
            ),
        ]:
            if subset["group"].nunique() < 2:
                continue
            observed, p_value, n_perm = permutation_test(
                subset["mean_expression"].to_numpy(),
                subset["group"].to_numpy(),
            )
            rows.append(
                {
                    "state": state,
                    "comparison": comparison,
                    "observed_difference": observed,
                    "p_permutation": p_value,
                    "n_permutations": n_perm,
                    "n_case": int(subset["group"].eq("UUO").sum()),
                    "n_control": int(subset["group"].eq("Normal").sum()),
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = bh_adjust(result["p_permutation"])
    return result


def leave_one_out(expression: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for state in STATES:
        block = expression[expression["state"].eq(state)]
        donor_means = (
            block.groupby(["donor", "group"], observed=True)["PAPPA2"]
            .agg(["mean", "size"])
            .reset_index()
            .rename(columns={"mean": "mean_expression", "size": "n_cells"})
        )
        donor_means = donor_means[donor_means["n_cells"] >= 5]
        for donor in donor_means["donor"]:
            kept = donor_means[~donor_means["donor"].eq(donor)]
            if kept["group"].nunique() < 2:
                continue
            difference = (
                kept.loc[kept["group"].eq("UUO"), "mean_expression"].mean()
                - kept.loc[kept["group"].eq("Normal"), "mean_expression"].mean()
            )
            rows.append({"state": state, "removed_donor": donor, "difference": difference})
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    return (
        result.groupby("state", observed=True)["difference"]
        .agg(["min", "max", "mean"])
        .reset_index()
    )


def within_donor(expression: pd.DataFrame, exposure: str) -> pd.DataFrame:
    rows = []
    for state in STATES:
        block = expression[expression["state"].eq(state)]
        donor_rows = []
        for donor, donor_block in block.groupby("donor", observed=True):
            valid = donor_block[["PAPPA2", exposure]].dropna()
            if (
                len(valid) < 20
                or valid[exposure].nunique() < 3
                or valid["PAPPA2"].nunique() < 2
            ):
                continue
            rho, _ = spearmanr(valid["PAPPA2"], valid[exposure])
            donor_rows.append(
                {
                    "donor": donor,
                    "group": donor_block["group"].iloc[0],
                    "rho": rho,
                }
            )
        donor_table = pd.DataFrame(donor_rows)
        if donor_table.empty or donor_table["group"].nunique() < 2:
            continue
        observed, p_value, n_perm = permutation_test(
            donor_table["rho"].to_numpy(),
            donor_table["group"].to_numpy(),
        )
        rows.append(
            {
                "state": state,
                "exposure": exposure,
                "observed_rho_difference": observed,
                "p_permutation": p_value,
                "n_permutations": n_perm,
                "n_case": int(donor_table["group"].eq("UUO").sum()),
                "n_control": int(donor_table["group"].eq("Normal").sum()),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = bh_adjust(result["p_permutation"])
    return result


def neighborhood_summary(expression: pd.DataFrame, exposure: str) -> pd.DataFrame:
    rows = []
    for state in STATES:
        block = expression[expression["state"].eq(state)]
        for column in [f"nhood_mean_{gene}" for gene in GENES]:
            if column not in block.columns:
                continue
            valid = block[[column, exposure]].dropna()
            if (
                len(valid) < 20
                or valid[column].nunique() < 2
                or valid[exposure].nunique() < 2
            ):
                continue
            rho, p_value = spearmanr(valid[column], valid[exposure])
            rows.append(
                {
                    "state": state,
                    "gene": column.removeprefix("nhood_mean_"),
                    "n_centers": len(valid),
                    "rho": rho,
                    "p_value": p_value,
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = bh_adjust(result["p_value"])
    else:
        result = pd.DataFrame(
            columns=["state", "gene", "n_centers", "rho", "p_value", "fdr"]
        )
    return result


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expression = load_expression()
    expression = add_injury_exposure(expression)
    expression = add_neighborhood_features(expression)

    pseudobulk = donor_pseudobulk(expression)
    leave_out = leave_one_out(expression)
    within_injury = within_donor(expression, "injury_exposure")
    within_nhood = within_donor(expression, "nhood_logfc")
    neighborhood = neighborhood_summary(expression, "nhood_logfc")

    expression.to_pickle(OUT_DIR / "GSE211785_TAL_expression.pkl")
    pseudobulk.to_csv(OUT_DIR / "GSE211785_pseudobulk.csv", index=False)
    leave_out.to_csv(OUT_DIR / "GSE211785_leave_one_out.csv", index=False)
    within_injury.to_csv(OUT_DIR / "GSE211785_within_injury.csv", index=False)
    within_nhood.to_csv(OUT_DIR / "GSE211785_within_nhood.csv", index=False)
    neighborhood.to_csv(OUT_DIR / "GSE211785_neighborhood_expression.csv", index=False)

    print("TAL cells:", len(expression))
    print(expression.groupby(["group_name", "state"]).size().to_string())
    print("\nDonor pseudobulk:")
    print(pseudobulk.to_string(index=False))
    print("\nLeave-one-donor-out:")
    print(leave_out.to_string(index=False))
    print("\nWithin donor, injury exposure:")
    print(within_injury.to_string(index=False))
    print("\nWithin donor, neighborhood logFC:")
    print(within_nhood.to_string(index=False))
    print("\nNeighborhood expression, PAPPA2:")
    if neighborhood.empty:
        print("No estimable neighborhood-expression tests")
    else:
        print(
            neighborhood[neighborhood["gene"].eq("PAPPA2")].to_string(index=False)
        )


if __name__ == "__main__":
    main()
