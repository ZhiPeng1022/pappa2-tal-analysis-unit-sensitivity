"""Test PAPPA2 direction conflict in predicted aTAL2-like external cells."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors


BASE = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji\work\signature_transfer")
GSE254_CELLS = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji\work\independent_validation_search"
    r"\gse254185_four_units\GSE254185_four_unit_cells.pkl"
)
GSE211_H5AD = Path(r"D:\Codex_temp\PAPPA2_independent_validation\GSE211785_PreSCVI.h5ad")


def permutation_test(
    values: np.ndarray,
    labels: np.ndarray,
    n_permutations: int = 20000,
    seed: int = 20261028,
) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    case = labels == "UUO"
    n_case = int(case.sum())
    if n_case < 2 or n_case == len(labels):
        return np.nan, np.nan, 0
    observed = values[case].mean() - values[~case].mean()
    total = len(values)
    combination_count = 1
    for i in range(1, n_case + 1):
        combination_count = combination_count * (total - i + 1) // i
        if combination_count > n_permutations:
            break
    rng = np.random.default_rng(seed)
    if combination_count <= n_permutations:
        chosen_list = list(combinations(range(total), n_case))
    else:
        chosen_list = [
            tuple(rng.choice(total, n_case, replace=False))
            for _ in range(n_permutations)
        ]
    null = np.empty(len(chosen_list))
    for i, chosen in enumerate(chosen_list):
        mask = np.zeros(total, dtype=bool)
        mask[list(chosen)] = True
        null[i] = values[mask].mean() - values[~mask].mean()
    p_value = (np.sum(np.abs(null) >= abs(observed)) + 1) / (len(null) + 1)
    return observed, p_value, len(null)


def summarize_frame(
    frame: pd.DataFrame,
    dataset: str,
    margin_threshold: float,
) -> list[dict]:
    subset = frame[
        frame["predicted_subtype"].eq("aTAL2")
        & frame["marker_margin"].ge(margin_threshold)
    ].copy()
    rows = []
    if subset.empty:
        return rows

    donor_means = (
        subset.groupby(["donor", "group"], observed=True)["PAPPA2"]
        .agg(["mean", "size"])
        .reset_index()
    )
    donor_means = donor_means[donor_means["size"] >= 5]
    if donor_means["group"].nunique() >= 2:
        observed, p_value, n_perm = permutation_test(
            donor_means["mean"].to_numpy(),
            donor_means["group"].to_numpy(),
        )
        rows.append(
            {
                "dataset": dataset,
                "margin_threshold": margin_threshold,
                "analysis": "donor pseudobulk PAPPA2",
                "n_case": int(donor_means["group"].eq("UUO").sum()),
                "n_control": int(donor_means["group"].eq("Normal").sum()),
                "effect": observed,
                "p_permutation": p_value,
                "n_permutations": n_perm,
            }
        )

    for exposure in ["injury_exposure", "nhood_logfc"]:
        if exposure not in subset.columns:
            continue
        donor_rows = []
        for donor, block in subset.groupby("donor", observed=True):
            valid = block[["PAPPA2", exposure]].dropna()
            if (
                len(valid) < 10
                or valid[exposure].nunique() < 3
                or valid["PAPPA2"].nunique() < 2
            ):
                continue
            rho, _ = spearmanr(valid["PAPPA2"], valid[exposure])
            donor_rows.append((donor, block["group"].iloc[0], rho))
        if len(donor_rows) < 4:
            continue
        values = np.asarray([row[2] for row in donor_rows])
        labels = np.asarray([row[1] for row in donor_rows])
        if len(set(labels)) < 2:
            continue
        observed, p_value, n_perm = permutation_test(values, labels)
        rows.append(
            {
                "dataset": dataset,
                "margin_threshold": margin_threshold,
                "analysis": f"within-donor PAPPA2 vs {exposure}",
                "n_case": int((labels == "UUO").sum()),
                "n_control": int((labels == "Normal").sum()),
                "effect": observed,
                "p_permutation": p_value,
                "n_permutations": n_perm,
            }
        )
    return rows


def prepare_gse254() -> pd.DataFrame:
    pred = pd.read_csv(BASE / "GSE254185_signature_predictions.csv")
    cells = pd.read_pickle(GSE254_CELLS)
    cells = cells.reset_index().rename(columns={"index": "cell_id"})
    cells["cell_key"] = cells["Library"].astype(str) + "|" + cells["barcode"].astype(str)
    pred["cell_key"] = pred["Library"].astype(str) + "|" + pred["barcode"].astype(str)
    merged = pred.merge(
        cells[
            [
                "cell_key",
                "injury_exposure",
                "nhood_logfc",
            ]
        ],
        on="cell_key",
        how="inner",
    )
    # Keep the prediction-normalized PAPPA2, which is on the same scale as the
    # external marker matrix.
    merged = merged.drop_duplicates("cell_key")
    return merged


def prepare_gse211() -> pd.DataFrame:
    pred = pd.read_csv(BASE / "GSE211785_signature_predictions.csv")
    adata = ad.read_h5ad(GSE211_H5AD, backed="r")
    obs = adata.obs
    mask = (
        obs["tech"].isin(["SC_RNA", "SN_RNA"])
        & obs["Cluster_Idents"].isin(["C_TAL", "M_TAL"])
        & obs["sample"].notna()
    )
    obs_tal = obs.loc[mask].copy()
    if len(pred) != len(obs_tal):
        raise RuntimeError(f"Prediction/obs length mismatch: {len(pred)} vs {len(obs_tal)}")
    pred = pred.copy()
    pred["donor"] = obs_tal["sample"].astype(str).to_numpy()
    pred["group"] = np.where(
        obs_tal["group"].astype(str).eq("Control").to_numpy(),
        "Normal",
        "UUO",
    )
    pred["tech"] = obs_tal["tech"].astype(str).to_numpy()
    umap = adata.obsm["X_umap"][mask.to_numpy()]
    pred["UMAP_1"] = umap[:, 0]
    pred["UMAP_2"] = umap[:, 1]

    # Local disease logFC within each technical library.
    pred["nhood_logfc"] = np.nan
    for tech, block in pred.groupby("tech", observed=True):
        coords = block[["UMAP_1", "UMAP_2"]].to_numpy(float)
        valid = np.isfinite(coords).all(1)
        if valid.sum() < 20:
            continue
        n_neighbors = min(31, valid.sum())
        model = NearestNeighbors(n_neighbors=n_neighbors).fit(coords[valid])
        _, neighbors = model.kneighbors(coords[valid])
        if n_neighbors > 1:
            neighbors = neighbors[:, 1:]
        labels = block["group"].to_numpy()
        is_case = labels == "UUO"
        total_case = max(int(is_case.sum()), 1)
        total_control = max(int((~is_case).sum()), 1)
        n_case = is_case[neighbors].sum(axis=1)
        n_control = (~is_case)[neighbors].sum(axis=1)
        local = np.log2(
            ((n_case + 0.5) / (total_case + 1))
            / ((n_control + 0.5) / (total_control + 1))
        )
        pred.loc[block.index[valid], "nhood_logfc"] = local
    return pred


def main() -> None:
    rows = []
    gse254 = prepare_gse254()
    gse211 = prepare_gse211()
    for threshold in [0.0, 0.1, 0.2]:
        rows.extend(summarize_frame(gse254, "GSE254185", threshold))
        rows.extend(summarize_frame(gse211, "GSE211785", threshold))
    result = pd.DataFrame(rows)
    result.to_csv(BASE / "signature_transfer_validation_summary.csv", index=False)
    print("Predicted aTAL2 counts:")
    for name, frame in [("GSE254185", gse254), ("GSE211785", gse211)]:
        subset = frame[frame["predicted_subtype"].eq("aTAL2")]
        print(name, len(subset), "donors", subset["donor"].nunique())
    print("\nValidation summary:")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
