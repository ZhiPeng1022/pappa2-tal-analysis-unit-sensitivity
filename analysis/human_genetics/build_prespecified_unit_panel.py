"""Build the frozen 400-gene analysis-unit panel and 2,000-gene control pool."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import rankdata


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
KPMP_H5AD = Path(r"D:\CodexData\RAMP3\KPMP_new_snRNA.h5ad")
MARKERS_CSV = ROOT / "work" / "signature_transfer" / "KPMP_TAL_markers.csv"
HALLMARK_GMT = ROOT / "work" / "MSigDB_Hallmark_2020.gmt"
OUT_DIR = ROOT / "outputs" / "prespecified_unit_panel"
SCAN_CACHE = OUT_DIR / "TAL_scan_cache.npz"

LABELS = ["C-TAL-A", "C-TAL-B", "aTAL2", "frTAL"]
SEED = 20261028
CHUNK_SIZE = 10000
PAPPA2 = "PAPPA2"
TARGET_PER_GROUP = 100
CONTROL_POOL_SIZE = 2000

HALLMARK_NAMES = {
    "TNF-alpha Signaling via NF-kB",
    "Inflammatory Response",
    "Hypoxia",
    "Oxidative Phosphorylation",
    "Epithelial Mesenchymal Transition",
}


def symbol_to_unique_index(
    feature_names: np.ndarray,
) -> dict[str, int]:
    frame = pd.DataFrame({"symbol": feature_names})
    counts = frame["symbol"].value_counts()
    unique_symbols = set(counts[counts.eq(1)].index)
    return {
        str(symbol): int(index)
        for index, symbol in enumerate(feature_names)
        if symbol in unique_symbols
        and symbol
        and symbol != "nan"
    }


def split_donors(obs: pd.DataFrame, labels: np.ndarray) -> tuple[set[str], set[str]]:
    rng = np.random.default_rng(SEED)
    donors = np.asarray(sorted(obs["donor_id"].astype(str).unique()))
    rng.shuffle(donors)
    n_test = max(10, int(round(len(donors) * 0.3)))
    test = set(donors[:n_test])
    train = set(donors[n_test:])
    donor_labels = (
        pd.DataFrame(
            {
                "donor_id": obs["donor_id"].astype(str).to_numpy(),
                "label": labels,
            }
        )
        .drop_duplicates()
    )
    for label in LABELS:
        present = donor_labels.loc[
            donor_labels["label"].eq(label),
            "donor_id",
        ]
        if not test.intersection(set(present)):
            candidate = next(
                (donor for donor in present if donor not in test),
                None,
            )
            if candidate:
                test.add(candidate)
                train.discard(candidate)
    return train, test


def read_hallmark_genes() -> set[str]:
    genes: set[str] = set()
    with HALLMARK_GMT.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in HALLMARK_NAMES:
                genes.update(
                    gene.strip().upper()
                    for gene in fields[2:]
                    if gene.strip()
                )
    if not genes:
        raise ValueError("No Hallmark genes parsed")
    return genes


def scan_tal_expression(
    adata: ad.AnnData,
    indices: np.ndarray,
    obs_tal: pd.DataFrame,
) -> dict[str, np.ndarray]:
    donor_ids = np.asarray(
        sorted(obs_tal["donor_id"].astype(str).unique())
    )
    donor_lookup = {donor: index for index, donor in enumerate(donor_ids)}
    donor_codes = np.asarray(
        [
            donor_lookup[donor]
            for donor in obs_tal["donor_id"].astype(str)
        ],
        dtype=np.int32,
    )
    n_donors = len(donor_ids)
    n_genes = adata.n_vars
    sum_lognorm = np.zeros(n_genes, dtype=np.float64)
    detection = np.zeros(n_genes, dtype=np.float64)
    donor_sum_lognorm = np.zeros((n_donors, n_genes), dtype=np.float64)
    donor_cells = np.bincount(donor_codes, minlength=n_donors)
    n_counts = pd.to_numeric(
        obs_tal["nCount_RNA"],
        errors="coerce",
    ).fillna(0).to_numpy(dtype=float)

    for start in range(0, len(indices), CHUNK_SIZE):
        chunk_indices = indices[start : start + CHUNK_SIZE]
        chunk_donor_codes = donor_codes[start : start + CHUNK_SIZE]
        chunk_counts = n_counts[start : start + CHUNK_SIZE]
        x = (
            adata[chunk_indices, :]
            .to_memory()
            .X.tocsr()
            .astype(np.float32)
        )
        scale = 10000.0 / np.maximum(chunk_counts, 1.0)
        x_norm = x.multiply(scale[:, None])
        x_norm.data = np.log1p(x_norm.data)

        sum_lognorm += np.asarray(
            x_norm.sum(axis=0)
        ).ravel()
        detection += np.asarray(
            x.getnnz(axis=0)
        ).ravel()

        n_chunk = x.shape[0]
        donor_indicator = csr_matrix(
            (
                np.ones(n_chunk, dtype=np.float32),
                (
                    chunk_donor_codes,
                    np.arange(n_chunk),
                ),
            ),
            shape=(n_donors, n_chunk),
        )
        donor_sum_lognorm += np.asarray(
            (donor_indicator @ x_norm).toarray(),
            dtype=np.float64,
        )
        print(
            f"processed {min(start + CHUNK_SIZE, len(indices)):,}"
            f"/{len(indices):,} TAL cells",
            flush=True,
        )

    return {
        "donor_ids": donor_ids,
        "donor_codes": donor_codes,
        "donor_sum_lognorm": donor_sum_lognorm,
        "donor_cells": donor_cells,
        "sum_lognorm": sum_lognorm,
        "detection": detection,
        "n_cells": len(indices),
    }


def select_marker_group(
    markers: pd.DataFrame,
    unique_symbols: set[str],
) -> list[str]:
    selected: list[str] = []
    for subtype in LABELS:
        candidates = (
            markers.loc[markers["subtype"].eq(subtype)]
            .sort_values("specificity", ascending=False)["gene"]
            .astype(str)
            .tolist()
        )
        added = 0
        for gene in candidates:
            if gene in unique_symbols and gene not in selected:
                selected.append(gene)
                added += 1
                if added == 25:
                    break
    if len(selected) != 100:
        raise ValueError(
            f"Expected 100 marker genes, found {len(selected)}"
        )
    return selected


def rank_correlation_with_vector(
    matrix: np.ndarray,
    vector: np.ndarray,
) -> np.ndarray:
    ranked_matrix = rankdata(matrix, axis=0).astype(np.float64)
    ranked_vector = rankdata(vector).astype(np.float64)
    x = ranked_matrix - ranked_matrix.mean(axis=0)
    y = ranked_vector - ranked_vector.mean()
    numerator = x.T @ y
    denominator = np.sqrt(
        np.sum(x**2, axis=0) * np.sum(y**2)
    )
    return np.divide(
        numerator,
        denominator,
        out=np.full(matrix.shape[1], np.nan),
        where=denominator > 0,
    )


def select_expression_matched_pool(
    candidate_symbols: list[str],
    gene_means: np.ndarray,
    gene_detection: np.ndarray,
    symbol_to_index: dict[str, int],
    target_means: np.ndarray,
    target_detection: np.ndarray,
) -> list[str]:
    candidates = np.asarray(
        [
            symbol
            for symbol in candidate_symbols
            if symbol in symbol_to_index
        ]
    )
    idx = np.asarray([symbol_to_index[gene] for gene in candidates])
    candidate_mean = gene_means[idx]
    candidate_detection = gene_detection[idx]
    candidate_position = {
        gene: position for position, gene in enumerate(candidates)
    }

    quantiles = np.quantile(target_means, np.linspace(0, 1, 11))
    quantiles[0] = -np.inf
    quantiles[-1] = np.inf
    rng = np.random.default_rng(SEED)
    chosen: list[str] = []
    for bin_index in range(10):
        mask = (
            (candidate_mean >= quantiles[bin_index])
            & (candidate_mean <= quantiles[bin_index + 1])
        )
        block_symbols = candidates[mask]
        block_mean = candidate_mean[mask]
        block_detection = candidate_detection[mask]
        if len(block_symbols) < 210:
            global_distance = (
                np.abs(
                    (candidate_mean - np.median(target_means))
                    / (np.nanstd(target_means) + 1e-8)
                )
                + np.abs(
                    (candidate_detection - np.median(target_detection))
                    / (np.nanstd(target_detection) + 1e-8)
                )
            )
            nearest = np.argsort(global_distance)
            block_symbols = candidates[nearest]
            block_detection = candidate_detection[nearest]
        target_detection_center = np.median(
            target_detection[
                (target_means >= quantiles[bin_index])
                & (target_means <= quantiles[bin_index + 1])
            ]
        )
        if not np.isfinite(target_detection_center):
            target_detection_center = np.median(target_detection)
        order = np.argsort(
            np.abs(block_detection - target_detection_center)
        )
        selected = [
            gene
            for gene in block_symbols[order]
            if gene not in chosen
        ][:210]
        chosen.extend(selected)

    chosen = list(dict.fromkeys(chosen))
    if len(chosen) < 2100:
        remaining = [
            gene for gene in candidates if gene not in set(chosen)
        ]
        remaining_idx = np.asarray(
            [candidate_position[gene] for gene in remaining]
        )
        distance = np.min(
            (
                (
                    candidate_mean[remaining_idx, None]
                    - target_means[None, :]
                )
                / (np.nanstd(target_means) + 1e-8)
            )
            ** 2
            + (
                (
                    candidate_detection[remaining_idx, None]
                    - target_detection[None, :]
                )
                / (np.nanstd(target_detection) + 1e-8)
            )
            ** 2,
            axis=1,
        )
        needed = 2100 - len(chosen)
        chosen.extend(
            np.asarray(remaining)[np.argsort(distance)][:needed].tolist()
        )
    rng.shuffle(chosen)
    return chosen[:2100]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(KPMP_H5AD, backed="r")
    obs = adata.obs
    subclass = obs["SubclassLevel3"].astype(str)
    indices = np.flatnonzero(subclass.isin(LABELS).to_numpy())
    obs_tal = obs.iloc[indices].copy()
    labels = obs_tal["SubclassLevel3"].astype(str).to_numpy()
    feature_names = adata.var["feature_name"].astype(str).to_numpy()
    biotypes = adata.var["feature_biotype"].astype(str).to_numpy()
    symbol_to_index = symbol_to_unique_index(feature_names)
    unique_symbols = set(symbol_to_index)

    if SCAN_CACHE.exists():
        cached = np.load(SCAN_CACHE, allow_pickle=False)
        scan = {
            "donor_ids": cached["donor_ids"].astype(str),
            "donor_codes": cached["donor_codes"].astype(np.int32),
            "donor_sum_lognorm": cached["donor_sum_lognorm"],
            "donor_cells": cached["donor_cells"],
            "sum_lognorm": cached["sum_lognorm"],
            "detection": cached["detection"],
            "n_cells": int(cached["n_cells"]),
        }
        print(f"Loaded scan cache: {SCAN_CACHE}", flush=True)
    else:
        scan = scan_tal_expression(adata, indices, obs_tal)
        np.savez(
            SCAN_CACHE,
            donor_ids=scan["donor_ids"].astype(str),
            donor_codes=scan["donor_codes"],
            donor_sum_lognorm=scan["donor_sum_lognorm"],
            donor_cells=scan["donor_cells"],
            sum_lognorm=scan["sum_lognorm"],
            detection=scan["detection"],
            n_cells=np.asarray(scan["n_cells"]),
        )
    mean_lognorm = scan["sum_lognorm"] / scan["n_cells"]
    detection_fraction = scan["detection"] / scan["n_cells"]
    donor_mean = (
        scan["donor_sum_lognorm"]
        / np.maximum(scan["donor_cells"], 1)[:, None]
    )
    donor_ids = scan["donor_ids"]
    train_donors, test_donors = split_donors(obs_tal, labels)
    train_mask = np.asarray(
        [donor in train_donors for donor in donor_ids]
    )
    if train_mask.sum() != 156 or (~train_mask).sum() != 67:
        raise ValueError(
            f"Unexpected donor split: {train_mask.sum()} train, "
            f"{(~train_mask).sum()} test"
        )

    pappa2_index = symbol_to_index[PAPPA2]
    pappa2_vector = donor_mean[train_mask, pappa2_index]
    eligible = np.flatnonzero(detection_fraction >= 0.10)
    eligible = np.asarray(
        [
            index
            for index in eligible
            if feature_names[index] in unique_symbols
            and feature_names[index] != PAPPA2
        ],
        dtype=np.int64,
    )
    coexpression = rank_correlation_with_vector(
        donor_mean[train_mask][:, eligible],
        pappa2_vector,
    )
    coexpression_frame = pd.DataFrame(
        {
            "gene": feature_names[eligible],
            "pappa2_donor_rho": coexpression,
            "mean_lognorm": mean_lognorm[eligible],
            "detection_fraction": detection_fraction[eligible],
        }
    ).dropna()
    coexpression_frame = coexpression_frame.sort_values(
        "pappa2_donor_rho",
        key=lambda series: series.abs(),
        ascending=False,
    )
    coexpression_frame.to_csv(
        OUT_DIR / "PAPPA2_training_donor_coexpression.csv",
        index=False,
    )

    markers = pd.read_csv(MARKERS_CSV)
    marker_genes = select_marker_group(markers, unique_symbols)
    assigned = set(marker_genes)

    hallmark_symbols = read_hallmark_genes()
    hallmark_candidates = [
        gene
        for gene in hallmark_symbols
        if gene in unique_symbols
        and gene not in assigned
        and detection_fraction[symbol_to_index[gene]] >= 0.10
    ]
    hallmark_frame = pd.DataFrame(
        {
            "gene": hallmark_candidates,
            "mean_lognorm": [
                mean_lognorm[symbol_to_index[gene]]
                for gene in hallmark_candidates
            ],
            "detection_fraction": [
                detection_fraction[symbol_to_index[gene]]
                for gene in hallmark_candidates
            ],
        }
    ).sort_values(
        ["mean_lognorm", "detection_fraction"],
        ascending=False,
    )
    hallmark_genes = hallmark_frame["gene"].head(100).tolist()
    if len(hallmark_genes) != 100:
        raise ValueError("Insufficient Hallmark genes")
    assigned.update(hallmark_genes)
    hallmark_frame.to_csv(
        OUT_DIR / "Hallmark_candidate_ranking.csv",
        index=False,
    )

    coexpression_genes = [
        gene
        for gene in coexpression_frame["gene"]
        if gene not in assigned
    ][:100]
    if len(coexpression_genes) != 100:
        raise ValueError("Insufficient PAPPA2 coexpression genes")
    assigned.update(coexpression_genes)

    target_means = mean_lognorm[
        [symbol_to_index[gene] for gene in assigned]
    ]
    target_detection = detection_fraction[
        [symbol_to_index[gene] for gene in assigned]
    ]
    candidate_controls = [
        symbol
        for symbol in unique_symbols
        if symbol not in assigned
        and detection_fraction[symbol_to_index[symbol]] >= 0.10
    ]
    matched_controls = select_expression_matched_pool(
        candidate_controls,
        mean_lognorm,
        detection_fraction,
        symbol_to_index,
        target_means,
        target_detection,
    )
    internal_controls = matched_controls[:100]
    negative_pool = matched_controls[100:]
    if len(negative_pool) != 2000:
        raise ValueError("Negative-control pool must contain 2,000 genes")

    panel_rows = []
    for gene_set, genes in [
        ("tal_subtype_marker", marker_genes),
        ("injury_inflammation_metabolism", hallmark_genes),
        ("pappa2_coexpression", coexpression_genes),
        ("internal_negative_control", internal_controls),
    ]:
        for rank, gene in enumerate(genes, start=1):
            index = symbol_to_index[gene]
            panel_rows.append(
                {
                    "gene_set": gene_set,
                    "selection_rank": rank,
                    "gene": gene,
                    "feature_index": index,
                    "mean_lognorm": mean_lognorm[index],
                    "detection_fraction": detection_fraction[index],
                }
            )
    panel = pd.DataFrame(panel_rows)
    if len(panel) != 400 or panel["gene"].nunique() != 400:
        raise ValueError("Target panel must contain 400 unique genes")
    panel.to_csv(OUT_DIR / "target_panel_400.csv", index=False)

    control_rows = []
    for rank, gene in enumerate(negative_pool, start=1):
        index = symbol_to_index[gene]
        control_rows.append(
            {
                "selection_rank": rank,
                "gene": gene,
                "feature_index": index,
                "mean_lognorm": mean_lognorm[index],
                "detection_fraction": detection_fraction[index],
            }
        )
    control_pool = pd.DataFrame(control_rows)
    control_pool.to_csv(
        OUT_DIR / "negative_control_pool_2000.csv",
        index=False,
    )

    summary = pd.DataFrame(
        {
            "gene": feature_names,
            "feature_index": np.arange(len(feature_names)),
            "feature_biotype": biotypes,
            "mean_lognorm": mean_lognorm,
            "detection_fraction": detection_fraction,
        }
    )
    summary.to_csv(OUT_DIR / "TAL_gene_expression_summary.csv", index=False)

    panel_hash = hashlib.sha256(
        (OUT_DIR / "target_panel_400.csv").read_bytes()
    ).hexdigest()
    manifest = {
        "date": "2026-10-28",
        "seed": SEED,
        "tal_cells": int(len(indices)),
        "train_donors": int(train_mask.sum()),
        "test_donors": int((~train_mask).sum()),
        "target_panel_size": int(len(panel)),
        "negative_control_pool_size": int(len(control_pool)),
        "target_panel_sha256": panel_hash,
        "group_counts": panel["gene_set"].value_counts().to_dict(),
        "hallmark_sets": sorted(HALLMARK_NAMES),
    }
    (OUT_DIR / "panel_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))
    adata.file.close()


if __name__ == "__main__":
    main()
