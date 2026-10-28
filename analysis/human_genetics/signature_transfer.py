"""Signature transfer from KPMP TAL subtypes to non-KPMP cohorts.

This is a marker-based transfer analysis, not an exact-label replication.
Markers are learned from a donor-stratified KPMP training set, checked on
held-out KPMP donors, and then projected to GSE254185 and GSE211785.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from sklearn.metrics import accuracy_score, f1_score


KPMP_H5AD = Path(os.environ.get("KPMP_H5AD", r"D:\CodexData\RAMP3\KPMP_new_snRNA.h5ad"))
GSE254185_DIR = Path(
    os.environ.get(
        "GSE254185_DIR",
        r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji\work\independent_validation_search",
    )
)
GSE211785_H5AD = Path(
    os.environ.get(
        "GSE211785_H5AD",
        r"D:\Codex_temp\PAPPA2_independent_validation\GSE211785_PreSCVI.h5ad",
    )
)
OUT_DIR = Path(
    os.environ.get(
        "SIGNATURE_TRANSFER_OUT",
        r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji\work\signature_transfer",
    )
)

LABELS = ["C-TAL-A", "C-TAL-B", "aTAL2", "frTAL"]
SEED = 20261028
N_MARKERS = 35
CHUNK_SIZE = 10000


def label_indices(obs: pd.DataFrame) -> np.ndarray:
    return np.flatnonzero(obs["SubclassLevel3"].isin(LABELS).to_numpy())


def gene_symbols(adata: ad.AnnData) -> np.ndarray:
    return adata.var["feature_name"].astype(str).to_numpy()


def symbol_to_index(symbols: np.ndarray, gene: str) -> int | None:
    positions = np.flatnonzero(symbols == gene)
    return int(positions[0]) if len(positions) else None


def split_donors(obs: pd.DataFrame, labels: np.ndarray) -> tuple[set[str], set[str]]:
    rng = np.random.default_rng(SEED)
    donors = np.asarray(sorted(obs["donor_id"].astype(str).unique()))
    rng.shuffle(donors)
    n_test = max(10, int(round(len(donors) * 0.3)))
    test = set(donors[:n_test])
    train = set(donors[n_test:])
    # Ensure every class appears on both sides where possible.
    donor_labels = pd.DataFrame(
        {"donor_id": obs["donor_id"].astype(str).to_numpy(), "label": labels}
    ).drop_duplicates()
    for label in LABELS:
        present = donor_labels.loc[donor_labels["label"].eq(label), "donor_id"]
        if not test.intersection(set(present)):
            candidate = next((d for d in present if d not in test), None)
            if candidate:
                test.add(candidate)
                train.discard(candidate)
    return train, test


def accumulate_kpmp() -> tuple[
    ad.AnnData,
    np.ndarray,
    pd.DataFrame,
    np.ndarray,
    dict[str, np.ndarray],
]:
    adata = ad.read_h5ad(KPMP_H5AD, backed="r")
    obs = adata.obs
    indices_all = label_indices(obs)
    obs_tal = obs.iloc[indices_all].copy()
    labels = obs_tal["SubclassLevel3"].astype(str).to_numpy()
    train_donors, test_donors = split_donors(obs_tal, labels)
    is_train = obs_tal["donor_id"].astype(str).isin(train_donors).to_numpy()

    n_genes = adata.n_vars
    sum_all = np.zeros((len(LABELS), n_genes), dtype=np.float64)
    sumsq_all = np.zeros((len(LABELS), n_genes), dtype=np.float64)
    detect_all = np.zeros((len(LABELS), n_genes), dtype=np.float64)
    n_all = np.zeros(len(LABELS), dtype=np.float64)
    sum_train = np.zeros((len(LABELS), n_genes), dtype=np.float64)
    detect_train = np.zeros((len(LABELS), n_genes), dtype=np.float64)
    n_train = np.zeros(len(LABELS), dtype=np.float64)

    for start in range(0, len(indices_all), CHUNK_SIZE):
        chunk_indices = indices_all[start : start + CHUNK_SIZE]
        chunk_labels = labels[start : start + CHUNK_SIZE]
        chunk_train = is_train[start : start + CHUNK_SIZE]
        x = adata[chunk_indices, :].to_memory().X.tocsr().astype(np.float64)
        for label_index, label in enumerate(LABELS):
            row_mask = chunk_labels == label
            if not row_mask.any():
                continue
            block = x[row_mask]
            sum_all[label_index] += np.asarray(block.sum(axis=0)).ravel()
            sumsq_all[label_index] += np.asarray(block.multiply(block).sum(axis=0)).ravel()
            detect_all[label_index] += np.asarray(block.getnnz(axis=0)).ravel()
            n_all[label_index] += block.shape[0]

            train_mask = row_mask & chunk_train
            if train_mask.any():
                train_block = x[train_mask]
                sum_train[label_index] += np.asarray(
                    train_block.sum(axis=0)
                ).ravel()
                detect_train[label_index] += np.asarray(
                    train_block.getnnz(axis=0)
                ).ravel()
                n_train[label_index] += train_block.shape[0]
        print(f"processed {min(start + CHUNK_SIZE, len(indices_all))}/{len(indices_all)}")

    mean_all = sum_all / n_all[:, None]
    var_all = np.maximum(sumsq_all / n_all[:, None] - mean_all**2, 1e-8)
    sd_all = np.sqrt(var_all)
    mean_train = sum_train / n_train[:, None]
    detect_fraction_train = detect_train / n_train[:, None]

    rows = []
    for i, label in enumerate(LABELS):
        other = np.max(np.delete(mean_train, i, axis=0), axis=0)
        specificity = mean_train[i] - other
        valid = (
            (mean_train[i] > 0.03)
            & (detect_fraction_train[i] > 0.03)
            & (specificity > 0)
        )
        candidates = np.flatnonzero(valid)
        candidates = candidates[np.argsort(specificity[candidates])[::-1][:N_MARKERS]]
        for position in candidates:
            rows.append(
                {
                    "subtype": label,
                    "gene_index": int(position),
                    "gene": gene_symbols(adata)[position],
                    "train_mean": float(mean_train[i, position]),
                    "train_detection_fraction": float(
                        detect_fraction_train[i, position]
                    ),
                    "specificity": float(specificity[position]),
                }
            )

    markers = pd.DataFrame(rows)
    stats = {
        "all_mean": mean_all,
        "all_sd": sd_all,
        "train_donors": sorted(train_donors),
        "test_donors": sorted(test_donors),
        "indices_all": indices_all,
        "obs_tal": obs_tal,
        "labels": labels,
        "is_train": is_train,
        "markers": markers,
    }
    return adata, indices_all, obs_tal, labels, stats


def marker_index_map(
    symbols: np.ndarray, markers: pd.DataFrame
) -> dict[str, list[int]]:
    result = {}
    for label in LABELS:
        genes = markers.loc[markers["subtype"].eq(label), "gene"].tolist()
        result[label] = [
            index
            for gene in genes
            if (index := symbol_to_index(symbols, gene)) is not None
        ]
    return result


def zscore_block(values: np.ndarray, mean: np.ndarray, sd: np.ndarray) -> np.ndarray:
    return np.clip((values - mean) / sd, -5, 5)


def scores_from_z(
    z: np.ndarray, marker_indices: dict[str, list[int]]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    scores = np.column_stack(
        [
            np.mean(z[:, marker_indices[label]], axis=1)
            if marker_indices[label]
            else np.full(z.shape[0], np.nan)
            for label in LABELS
        ]
    )
    order = np.argsort(scores, axis=1)
    best = order[:, -1]
    second = order[:, -2]
    best_scores = scores[np.arange(len(scores)), best]
    second_scores = scores[np.arange(len(scores)), second]
    return best, best_scores - second_scores, scores


def internal_validation(
    adata: ad.AnnData,
    stats: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    markers = stats["markers"]
    symbols = gene_symbols(adata)
    marker_map = marker_index_map(symbols, markers)
    all_gene_indices = sorted({i for values in marker_map.values() for i in values})
    gene_to_col = {gene_index: position for position, gene_index in enumerate(all_gene_indices)}
    obs_tal = stats["obs_tal"]
    test_indices = stats["indices_all"][~stats["is_train"]]
    test_labels = stats["labels"][~stats["is_train"]]
    mean = stats["all_mean"][:, all_gene_indices].T
    sd = stats["all_sd"][:, all_gene_indices].T
    # Use the global per-gene KPMP mean/sd for held-out KPMP cells.
    gene_mean = np.nanmean(mean, axis=1)
    gene_sd = np.nanmean(sd, axis=1)
    mapped = {
        label: [gene_to_col[i] for i in values]
        for label, values in marker_map.items()
    }
    predicted_parts = []
    margin_parts = []
    for start in range(0, len(test_indices), 2000):
        chunk = adata[test_indices[start : start + 2000], :].to_memory()
        x = chunk.X.tocsr()[:, all_gene_indices].toarray()
        z = zscore_block(x, gene_mean, gene_sd)
        best, margin, _ = scores_from_z(z, mapped)
        predicted_parts.append(np.asarray(LABELS)[best])
        margin_parts.append(margin)
    predicted = np.concatenate(predicted_parts)
    margins = np.concatenate(margin_parts)
    result = pd.DataFrame(
        {
            "true_subtype": test_labels,
            "predicted_subtype": predicted,
            "margin": margins,
        }
    )
    summary_rows = []
    for threshold in [0.0, 0.1, 0.2, 0.3]:
        assigned = result["margin"].ge(threshold)
        true_assigned = test_labels[assigned]
        pred_assigned = predicted[assigned]
        a_tal2_true = true_assigned == "aTAL2"
        a_tal2_pred = pred_assigned == "aTAL2"
        tp = int((a_tal2_true & a_tal2_pred).sum())
        precision = tp / int(a_tal2_pred.sum()) if a_tal2_pred.sum() else np.nan
        recall = tp / int(a_tal2_true.sum()) if a_tal2_true.sum() else np.nan
        summary_rows.append(
            {
                "n_test_cells": len(result),
                "margin_threshold": threshold,
                "n_assigned": int(assigned.sum()),
                "coverage": float(assigned.mean()),
                "accuracy_assigned": accuracy_score(
                    true_assigned, pred_assigned
                ),
                "macro_f1_assigned": f1_score(
                    true_assigned,
                    pred_assigned,
                    labels=LABELS,
                    average="macro",
                ),
                "aTAL2_precision": precision,
                "aTAL2_recall": recall,
                "train_donors": len(stats["train_donors"]),
                "test_donors": len(stats["test_donors"]),
            }
        )
    return pd.DataFrame(summary_rows), result


def load_gse254185(markers: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    meta = pd.read_csv(GSE254185_DIR / "GSE254185_cell_meta_data.tsv.gz", sep="\t", index_col=0)
    meta["barcode"] = meta.index.str.replace(r"^l\d+_", "", regex=True)
    meta = meta[meta["Annotation.Lvl1"].astype(str).eq("TAL")].copy()
    meta["group"] = np.where(
        meta["Sample"].astype(str).str.startswith("UUO"), "UUO", "Normal"
    )
    marker_genes = markers["gene"].drop_duplicates().tolist()
    records = []
    for library_number in range(1, 7):
        library = f"Library{library_number}"
        path = GSE254185_DIR / f"GSE254185_{library}_gex_matrix.h5"
        with h5py.File(path, "r") as handle:
            group = handle["unknown"]
            genes = np.array([x.decode() for x in group["gene_names"][:]])
            barcodes = np.array([x.decode() for x in group["barcodes"][:]])
            shape = tuple(int(x) for x in group["shape"][:])
            matrix = csc_matrix(
                (group["data"][:], group["indices"][:], group["indptr"][:]),
                shape=shape,
            )
        # The matrix is genes x cells.
        gene_indices = [i for i, g in enumerate(genes) if g in set(marker_genes + ["PAPPA2"])]
        selected_genes = [genes[i] for i in gene_indices]
        selected = matrix[gene_indices, :].toarray().T
        total = np.asarray(matrix.sum(axis=0)).ravel()
        safe_total = np.where(total > 0, total, np.nan)
        norm = np.log1p(selected / safe_total[:, None] * 10000)
        barcode_index = {b: i for i, b in enumerate(barcodes)}
        lib_meta = meta[meta["Library"].eq(library)]
        mapped_rows = []
        mapped_positions = []
        for row_id, row in lib_meta.iterrows():
            position = barcode_index.get(row["barcode"])
            if position is not None:
                mapped_rows.append(row_id)
                mapped_positions.append(position)
        if not mapped_rows:
            continue
        block = lib_meta.loc[mapped_rows].copy()
        values = norm[np.asarray(mapped_positions)]
        for j, gene in enumerate(selected_genes):
            block[gene] = values[:, j]
        records.append(block)
    frame = pd.concat(records, ignore_index=False)
    frame["donor"] = frame["Sample"].astype(str)
    return frame, marker_genes


def predict_external(
    frame: pd.DataFrame,
    markers: pd.DataFrame,
    label: str,
) -> pd.DataFrame:
    marker_map = {}
    for subtype in LABELS:
        genes = [
            gene
            for gene in markers.loc[markers["subtype"].eq(subtype), "gene"]
            if gene in frame.columns
        ]
        marker_map[subtype] = genes
    genes = sorted({gene for values in marker_map.values() for gene in values})
    x = frame[genes].to_numpy(dtype=float)
    means = np.nanmean(x, axis=0)
    sds = np.nanstd(x, axis=0)
    sds = np.where(sds > 0, sds, 1)
    z = np.clip((x - means) / sds, -5, 5)
    positions = {gene: i for i, gene in enumerate(genes)}
    mapped = {
        subtype: [positions[g] for g in values if g in positions]
        for subtype, values in marker_map.items()
    }
    best, margin, score_matrix = scores_from_z(z, mapped)
    result = frame.copy()
    result["predicted_subtype"] = np.asarray(LABELS)[best]
    result["marker_margin"] = margin
    result["prediction_dataset"] = label
    return result


def summarize_predictions(predictions: pd.DataFrame, label: str) -> pd.DataFrame:
    return (
        predictions.groupby(["predicted_subtype", "group"], observed=True)
        .agg(
            cells=("predicted_subtype", "size"),
            donors=("donor", "nunique")
            if "donor" in predictions.columns
            else ("sample", "nunique")
            if "sample" in predictions.columns
            else ("predicted_subtype", "size"),
            pappa2_mean=("PAPPA2", "mean"),
        )
        .reset_index()
        .assign(dataset=label)
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adata, indices, obs_tal, labels, stats = accumulate_kpmp()
    markers = stats["markers"]
    markers.to_csv(OUT_DIR / "KPMP_TAL_markers.csv", index=False)
    validation, validation_cells = internal_validation(adata, stats)
    validation.to_csv(OUT_DIR / "KPMP_internal_validation.csv", index=False)
    validation_cells.to_csv(
        OUT_DIR / "KPMP_internal_validation_predictions.csv", index=False
    )
    print("Internal validation:")
    print(validation.to_string(index=False))

    marker_genes = markers["gene"].drop_duplicates().tolist()
    frame254, _ = load_gse254185(markers)
    pred254 = predict_external(frame254, markers, "GSE254185")
    pred254.to_csv(OUT_DIR / "GSE254185_signature_predictions.csv", index=False)
    summary254 = summarize_predictions(pred254, "GSE254185")
    summary254.to_csv(OUT_DIR / "GSE254185_signature_summary.csv", index=False)

    # GSE211785: use the full non-KPMP H5AD and the same marker symbols.
    a211 = ad.read_h5ad(GSE211785_H5AD, backed="r")
    obs211 = a211.obs
    mask211 = (
        obs211["tech"].isin(["SC_RNA", "SN_RNA"])
        & obs211["Cluster_Idents"].isin(["C_TAL", "M_TAL"])
        & obs211["sample"].notna()
    )
    symbols211 = a211.var_names.astype(str).to_numpy()
    marker_indices = [
        i for i in (symbol_to_index(symbols211, gene) for gene in marker_genes)
        if i is not None
    ]
    marker_symbols = [symbols211[i] for i in marker_indices]
    pappa2_index = symbol_to_index(symbols211, "PAPPA2")
    selected_indices = marker_indices + ([pappa2_index] if pappa2_index is not None else [])
    selected_names = marker_symbols + (["PAPPA2"] if pappa2_index is not None else [])
    mask_indices = np.flatnonzero(mask211.to_numpy())
    frames211 = []
    for start in range(0, len(mask_indices), 2000):
        chunk_indices = mask_indices[start : start + 2000]
        chunk = a211[chunk_indices, :].to_memory()
        matrix = chunk.X.tocsr()[:, selected_indices].toarray()
        frame_chunk = pd.DataFrame(
            matrix,
            columns=selected_names,
            index=chunk.obs.index,
        )
        frame_chunk["donor"] = chunk.obs["sample"].astype(str).to_numpy()
        frame_chunk["group"] = chunk.obs["group"].astype(str).to_numpy()
        frame_chunk["true_state"] = chunk.obs["Cluster_Idents"].astype(str).to_numpy()
        frames211.append(frame_chunk)
    frame211 = pd.concat(frames211, ignore_index=False)
    pred211 = predict_external(frame211, markers, "GSE211785")
    pred211.to_csv(OUT_DIR / "GSE211785_signature_predictions.csv", index=False)
    summary211 = summarize_predictions(pred211, "GSE211785")
    summary211.to_csv(OUT_DIR / "GSE211785_signature_summary.csv", index=False)

    print("\nGSE254185 predicted summary:")
    print(summary254.to_string(index=False))
    print("\nGSE211785 predicted summary:")
    print(summary211.to_string(index=False))


if __name__ == "__main__":
    main()
