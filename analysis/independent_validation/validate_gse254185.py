"""Preliminary independent validation of analysis-unit sensitivity in GSE254185.

This script uses the public per-library 10x GEX matrices and the public cell
metadata. It compares two analysis units:

1. donor pseudobulk disease-versus-control expression
2. within-donor cell-rank association with local neighborhood injury exposure

The result is intended as a feasibility test, not a replacement for the full
original Milo workflow.
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.neighbors import NearestNeighbors


BASE = Path(
    os.environ.get(
        "PAPPA2_VALIDATION_DIR",
        str(Path(__file__).resolve().parent),
    )
)
META_PATH = BASE / "GSE254185_cell_meta_data.tsv.gz"
OUT_DIR = BASE / "gse254185_preliminary"

GENES = [
    "PAPPA2",
    "CXCL12",
    "NOTCH3",
    "HES1",
    "EPHB2",
    "HAVCR1",
    "VCAM1",
    "SOX9",
    "CLU",
    "KRT8",
    "KRT18",
]
INJURY_GENES = ["HAVCR1", "VCAM1", "SOX9", "CLU", "KRT8", "KRT18"]
STATE_ORDER = ["cTAL1", "cTAL2", "mTAL", "TAL Injured", "TAL Inflammatory"]


def bh_adjust(p_values: pd.Series) -> pd.Series:
    p = p_values.to_numpy(dtype=float)
    n = len(p)
    adjusted = np.full(n, np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return pd.Series(adjusted, index=p_values.index)
    order = np.argsort(p[valid])
    ranked = p[valid][order]
    q = ranked * valid.sum() / np.arange(1, valid.sum() + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(valid.sum())
    tmp[order] = q
    adjusted[valid] = tmp
    return pd.Series(adjusted, index=p_values.index)


def load_metadata() -> pd.DataFrame:
    meta = pd.read_csv(META_PATH, sep="\t", index_col=0)
    meta["barcode"] = meta.index.str.replace(r"^l\d+_", "", regex=True)
    meta["group"] = np.where(
        meta["Sample"].astype(str).str.startswith("UUO"), "UUO", "Normal"
    )
    meta["donor"] = meta["Sample"].astype(str)
    meta["state"] = meta["Annotation.Lvl2"].astype(str)
    meta = meta[meta["Annotation.Lvl1"].astype(str).eq("TAL")].copy()
    meta = meta[meta["state"].isin(STATE_ORDER)].copy()
    return meta


def load_expression(meta: pd.DataFrame) -> pd.DataFrame:
    records = []
    for library_number in range(1, 7):
        library = f"Library{library_number}"
        path = BASE / f"GSE254185_{library}_gex_matrix.h5"
        with h5py.File(path, "r") as handle:
            group = handle["unknown"]
            genes = np.array([x.decode() for x in group["gene_names"][:]])
            barcodes = np.array([x.decode() for x in group["barcodes"][:]])
            shape = tuple(int(x) for x in group["shape"][:])
            matrix = csc_matrix(
                (group["data"][:], group["indices"][:], group["indptr"][:]),
                shape=shape,
            )

        total_counts = np.asarray(matrix.sum(axis=0)).ravel()
        gene_positions = []
        available_genes = []
        for gene in GENES:
            positions = np.where(genes == gene)[0]
            if len(positions):
                gene_positions.append(int(positions[0]))
                available_genes.append(gene)

        selected_counts = matrix[gene_positions, :].toarray().T
        barcode_index = {barcode: i for i, barcode in enumerate(barcodes)}

        lib_meta = meta[meta["Library"].eq(library)].copy()
        mapped_rows = []
        mapped_positions = []
        for row_id, row in lib_meta.iterrows():
            position = barcode_index.get(row["barcode"])
            if position is not None:
                mapped_rows.append(row_id)
                mapped_positions.append(position)

        if not mapped_rows:
            continue

        positions = np.asarray(mapped_positions, dtype=int)
        counts = selected_counts[positions]
        totals = total_counts[positions]
        safe_totals = np.where(totals > 0, totals, np.nan)
        normalized = np.log1p(counts / safe_totals[:, None] * 10000.0)

        block = lib_meta.loc[mapped_rows].copy()
        for j, gene in enumerate(available_genes):
            block[gene] = normalized[:, j]
        block["total_counts"] = totals
        records.append(block)

    expression = pd.concat(records, ignore_index=False)
    expression = expression.dropna(subset=["PAPPA2"])
    return expression


def add_injury_exposure(expression: pd.DataFrame) -> pd.DataFrame:
    expression = expression.copy()
    expression["injury_score"] = np.nan
    available = [gene for gene in INJURY_GENES if gene in expression.columns]
    for library, block in expression.groupby("Library", sort=False):
        index = block.index
        z = []
        for gene in available:
            values = block[gene].astype(float)
            sd = values.std(ddof=0)
            if not np.isfinite(sd) or sd == 0:
                continue
            z.append((values - values.mean()) / sd)
        if z:
            expression.loc[index, "injury_score"] = np.mean(np.vstack(z), axis=0)

    coords = expression[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)
    valid = np.isfinite(coords).all(axis=1) & np.isfinite(
        expression["injury_score"].to_numpy(dtype=float)
    )
    exposure = np.full(len(expression), np.nan)
    if valid.sum() > 31:
        model = NearestNeighbors(n_neighbors=31, algorithm="auto")
        model.fit(coords[valid])
        _, neighbors = model.kneighbors(coords[valid])
        neighbor_scores = expression.loc[valid, "injury_score"].to_numpy(dtype=float)
        exposure[valid] = neighbor_scores[neighbors[:, 1:]].mean(axis=1)
    expression["injury_exposure"] = exposure
    return expression


def donor_pseudobulk(expression: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for state in STATE_ORDER:
        for gene in GENES:
            if gene not in expression.columns:
                continue
            block = expression[expression["state"].eq(state)]
            if block.empty:
                continue
            donor_means = (
                block.groupby(["donor", "group"], observed=True)[gene]
                .agg(["mean", "size"])
                .reset_index()
                .rename(columns={"mean": "mean_expression", "size": "n_cells"})
            )
            donor_means = donor_means[donor_means["n_cells"] >= 5]
            if donor_means["group"].nunique() < 2:
                continue
            uuo = donor_means.loc[donor_means["group"].eq("UUO"), "mean_expression"]
            normal = donor_means.loc[
                donor_means["group"].eq("Normal"), "mean_expression"
            ]
            if len(uuo) < 2 or len(normal) < 2:
                continue
            p_value = mannwhitneyu(uuo, normal, alternative="two-sided").pvalue
            rows.append(
                {
                    "state": state,
                    "gene": gene,
                    "n_uuo_donors": len(uuo),
                    "n_normal_donors": len(normal),
                    "uuo_mean": uuo.mean(),
                    "normal_mean": normal.mean(),
                    "mean_difference": uuo.mean() - normal.mean(),
                    "p_value": p_value,
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = bh_adjust(result["p_value"])
        result["direction"] = np.where(
            result["mean_difference"] > 0, "positive", "negative"
        )
    return result


def within_donor_rank(expression: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for state in STATE_ORDER:
        state_block = expression[expression["state"].eq(state)].copy()
        if state_block.empty:
            continue
        for donor, block in state_block.groupby("donor", sort=False):
            if len(block) < 20:
                continue
            if block["injury_exposure"].nunique(dropna=True) < 3:
                continue
            for gene in GENES:
                if gene not in block.columns:
                    continue
                valid = block[[gene, "injury_exposure"]].dropna()
                if (
                    len(valid) < 20
                    or valid["injury_exposure"].nunique() < 3
                    or valid[gene].nunique() < 2
                ):
                    continue
                rho, p_value = spearmanr(valid[gene], valid["injury_exposure"])
                rows.append(
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
    donor_result = pd.DataFrame(rows)
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


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta = load_metadata()
    print(f"TAL cells: {len(meta)}")
    print(meta.groupby(["group", "state"]).size().to_string())

    expression = load_expression(meta)
    expression = add_injury_exposure(expression)
    print(f"TAL cells with expression: {len(expression)}")

    pseudobulk = donor_pseudobulk(expression)
    within = within_donor_rank(expression)

    expression.to_pickle(OUT_DIR / "GSE254185_TAL_expression.pkl")
    pseudobulk.to_csv(OUT_DIR / "GSE254185_donor_pseudobulk.csv", index=False)
    within.to_csv(OUT_DIR / "GSE254185_within_donor_rank.csv", index=False)

    print("\nPAPPA2 donor pseudobulk:")
    print(
        pseudobulk[pseudobulk["gene"].eq("PAPPA2")][
            [
                "state",
                "n_uuo_donors",
                "n_normal_donors",
                "mean_difference",
                "p_value",
                "fdr",
                "direction",
            ]
        ].to_string(index=False)
    )
    print("\nPAPPA2 within-donor rank:")
    print(
        within[within["gene"].eq("PAPPA2")][
            [
                "state",
                "n_uuo_donors",
                "n_normal_donors",
                "mean_rho_difference",
                "p_value",
                "fdr",
                "direction",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
