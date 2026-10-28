"""Run the frozen 400-gene analysis-unit sensitivity test.

The script computes three gene-specific units:
1. donor pseudobulk;
2. neighborhood-expression association;
3. within-donor cell-rank association.

It then compares the target panel with 1,000 expression-matched random gene
sets drawn from the frozen 2,000-gene control pool.
"""

from __future__ import annotations

import json
import h5py
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix
from scipy.stats import mannwhitneyu, rankdata


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
KPMP_H5AD = Path(r"D:\CodexData\RAMP3\KPMP_new_snRNA.h5ad")
CELL_TABLE = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2\outputs"
    r"\KPMP_full_TAL_focus_expression.csv.gz"
)
MILO_DIR = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2\outputs"
)
MILO_PREFIX = "KPMP_milo_full_TAL_final"
PANEL_DIR = ROOT / "outputs" / "prespecified_unit_panel"
TARGET_PANEL = PANEL_DIR / "target_panel_400.csv"
CONTROL_POOL = PANEL_DIR / "negative_control_pool_2000.csv"
OUT_DIR = ROOT / "outputs" / "prespecified_unit_panel_results"
EXPRESSION_CACHE = OUT_DIR / "panel_expression_lognorm.float32.npy"

LABELS = ["C-TAL-A", "C-TAL-B", "aTAL2", "frTAL"]
COMPARISONS = {
    "all_disease": ["AKI", "CKD"],
    "AKI_vs_Reference": ["AKI"],
    "CKD_vs_Reference": ["CKD"],
}
CONDITION_MAP = {
    "Healthy Reference": "Reference",
    "Reference": "Reference",
    "Acute Kidney Injury": "AKI",
    "AKI": "AKI",
    "Chronic Kidney Disease": "CKD",
    "CKD": "CKD",
}
MIN_CELLS_PER_DONOR = 5
MIN_DONORS_PER_GROUP = 5
MIN_NHOODS = 20
GENE_CHUNK_SIZE = 200
Z_THRESHOLD = 0.5
N_RANDOM_SETS = 1000
SEED = 20261028


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(dtype=float)
    result = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return pd.Series(result, index=values.index)
    order = np.argsort(p[valid])
    ranked = p[valid][order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    adjusted = np.empty_like(q)
    adjusted[order] = np.clip(q, 0, 1)
    result[valid] = adjusted
    return pd.Series(result, index=values.index)


def build_expression_cache(
    h5ad_path: Path,
    gene_indices: np.ndarray,
    target_indices: np.ndarray,
    n_obs: int,
    n_vars: int,
    n_counts: np.ndarray,
) -> Path:
    EXPRESSION_CACHE.parent.mkdir(parents=True, exist_ok=True)
    if (
        EXPRESSION_CACHE.exists()
        and EXPRESSION_CACHE.stat().st_size
        >= len(target_indices) * len(gene_indices) * 4
    ):
        print(f"Using expression cache: {EXPRESSION_CACHE}", flush=True)
        return EXPRESSION_CACHE

    target_position = np.full(n_obs, -1, dtype=np.int64)
    target_position[target_indices] = np.arange(
        len(target_indices),
        dtype=np.int64,
    )
    gene_lookup = np.full(n_vars, -1, dtype=np.int32)
    gene_lookup[gene_indices] = np.arange(
        len(gene_indices),
        dtype=np.int32,
    )
    matrix = np.lib.format.open_memmap(
        EXPRESSION_CACHE,
        mode="w+",
        dtype=np.float32,
        shape=(len(target_indices), len(gene_indices)),
    )

    with h5py.File(h5ad_path, "r") as handle:
        indptr = handle["X/indptr"][:]
        indices_dataset = handle["X/indices"]
        data_dataset = handle["X/data"]
        for start in range(0, n_obs, 20000):
            end = min(start + 20000, n_obs)
            row_starts = indptr[start:end]
            row_ends = indptr[start + 1 : end + 1]
            data_start = int(row_starts[0])
            data_end = int(row_ends[-1])
            nonzero_indices = indices_dataset[data_start:data_end]
            nonzero_data = data_dataset[data_start:data_end]
            row_ids = np.repeat(
                np.arange(end - start, dtype=np.int64),
                row_ends - row_starts,
            )
            out_rows = target_position[start + row_ids]
            out_cols = gene_lookup[nonzero_indices]
            keep = (out_rows >= 0) & (out_cols >= 0)
            if keep.any():
                matrix[
                    out_rows[keep],
                    out_cols[keep],
                ] = nonzero_data[keep]
            if start % 200000 == 0:
                print(
                    f"read all-cell chunk {start:,}/{n_obs:,}",
                    flush=True,
                )

    for start in range(0, len(target_indices), 10000):
        end = min(start + 10000, len(target_indices))
        block = np.asarray(matrix[start:end])
        scale = 10000.0 / np.maximum(n_counts[start:end], 1.0)
        block *= scale[:, None]
        np.log1p(block, out=block)
        matrix[start:end] = block
    matrix.flush()
    print(f"Wrote expression cache: {EXPRESSION_CACHE}", flush=True)
    return EXPRESSION_CACHE


def rank_correlation(
    matrix: np.ndarray,
    exposure: np.ndarray,
    *,
    axis: int = 0,
) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=float)
    exposure = np.asarray(exposure, dtype=float)
    if axis == 1:
        matrix = matrix.T
    x_ranks = rankdata(matrix, axis=0).astype(float)
    y_ranks = rankdata(exposure).astype(float)
    x_centered = x_ranks - x_ranks.mean(axis=0)
    y_centered = y_ranks - y_ranks.mean()
    numerator = x_centered.T @ y_centered
    denominator = np.sqrt(
        np.sum(x_centered**2, axis=0) * np.sum(y_centered**2)
    )
    return np.divide(
        numerator,
        denominator,
        out=np.full(matrix.shape[1], np.nan),
        where=denominator > 0,
    )


def derive_study_design(cell_table: pd.DataFrame) -> dict:
    donor_condition = (
        cell_table.groupby("donor_id", observed=True)["condition"]
        .agg(lambda values: values.mode().iloc[0])
        .to_dict()
    )
    return {
        "donors": np.asarray(sorted(cell_table["donor_id"].unique())),
        "donor_condition": donor_condition,
    }


def build_milo_structures(
    cell_table: pd.DataFrame,
) -> dict[str, dict]:
    cell_to_position = {
        cell_id: position
        for position, cell_id in enumerate(cell_table["cell_id"])
    }
    structures = {}
    for comparison in COMPARISONS:
        membership_path = MILO_DIR / (
            f"{MILO_PREFIX}_{comparison}_nhood_membership.csv.gz"
        )
        result_path = MILO_DIR / (
            f"{MILO_PREFIX}_{comparison}_nhood_results.csv"
        )
        membership = pd.read_csv(
            membership_path,
            dtype={"cell_id": str, "donor_id": str},
        )
        results = pd.read_csv(result_path).rename(
            columns={"subclass": "nhood_subclass"}
        )
        edges = membership.merge(
            results[
                [
                    "Nhood",
                    "nhood_subclass",
                    "logFC",
                    "SpatialFDR",
                    "nhood_purity",
                ]
            ],
            on="Nhood",
            how="inner",
            validate="many_to_one",
        )
        nhood_ids = np.asarray(
            sorted(edges["Nhood"].unique()),
            dtype=int,
        )
        nhood_position = {
            nhood: position for position, nhood in enumerate(nhood_ids)
        }
        cell_positions = edges["cell_id"].map(cell_to_position)
        if cell_positions.isna().any():
            raise ValueError(
                f"{comparison}: membership contains cells outside table"
            )
        row_positions = edges["Nhood"].map(nhood_position).to_numpy()
        col_positions = cell_positions.to_numpy(dtype=int)
        row_counts = np.bincount(
            row_positions,
            minlength=len(nhood_ids),
        ).astype(float)
        raw_weights = coo_matrix(
            (
                np.ones(len(edges), dtype=float),
                (row_positions, col_positions),
            ),
            shape=(len(nhood_ids), len(cell_table)),
        ).tocsr()
        normalized_weights = (
            csr_matrix(
                (
                    1.0 / np.maximum(row_counts, 1.0),
                    (np.arange(len(nhood_ids)), np.arange(len(nhood_ids))),
                ),
                shape=(len(nhood_ids), len(nhood_ids)),
            )
            @ raw_weights
        ).tocsr()
        nhood_metadata = (
            edges[
                [
                    "Nhood",
                    "nhood_subclass",
                    "logFC",
                    "SpatialFDR",
                    "nhood_purity",
                ]
            ]
            .drop_duplicates("Nhood")
            .copy()
        )
        nhood_metadata["nhood_position"] = (
            nhood_metadata["Nhood"].map(nhood_position)
        )
        exposure = (
            edges.groupby("cell_id", observed=True)["logFC"]
            .mean()
            .rename("exposure")
        )
        exposure_by_cell = np.zeros(
            len(cell_table),
            dtype=float,
        )
        exposure_positions = cell_positions.to_numpy(dtype=int)
        exposure_values = exposure.reindex(
            edges["cell_id"].to_numpy()
        ).to_numpy(dtype=float)
        np.add.at(exposure_by_cell, exposure_positions, exposure_values)
        exposure_counts = np.bincount(
            exposure_positions,
            minlength=len(cell_table),
        )
        nonmissing = exposure_counts > 0
        exposure_by_cell[nonmissing] /= exposure_counts[nonmissing]

        subtype_structures = {}
        for subclass in LABELS:
            nhood_subset = nhood_metadata[
                nhood_metadata["nhood_subclass"].eq(subclass)
            ].copy()
            if nhood_subset.empty:
                subtype_structures[subclass] = {
                    "weight_matrix": None,
                    "nhood_logfc": np.asarray([], dtype=float),
                }
                continue
            row_subset = nhood_subset[
                "nhood_position"
            ].to_numpy(dtype=int)
            weight_matrix = normalized_weights[row_subset]
            subtype_structures[subclass] = {
                "weight_matrix": weight_matrix,
                "nhood_logfc": nhood_subset["logFC"].to_numpy(dtype=float),
            }

        structures[comparison] = {
            "edges": edges,
            "exposure_by_cell": exposure_by_cell,
            "subtypes": subtype_structures,
        }
    return structures


def donor_group_indices(
    cell_table: pd.DataFrame,
    donors: np.ndarray,
) -> dict[str, list[np.ndarray]]:
    donor_lookup = {
        donor: position for position, donor in enumerate(donors)
    }
    donor_codes = np.asarray(
        [donor_lookup[donor] for donor in cell_table["donor_id"]],
        dtype=np.int32,
    )
    subclass_values = cell_table["subclass"].astype(str).to_numpy()
    result = {}
    for subclass in LABELS:
        arrays = []
        for donor_index in range(len(donors)):
            indices = np.flatnonzero(
                (donor_codes == donor_index)
                & (subclass_values == subclass)
            )
            arrays.append(indices)
        result[subclass] = arrays
    return result


def within_donor_indices(
    milo_structure: dict,
    cell_table: pd.DataFrame,
    subclass: str,
    donors: np.ndarray,
) -> tuple[list[np.ndarray], np.ndarray]:
    edges = milo_structure["edges"]
    cell_ids = edges.loc[
        edges["nhood_subclass"].eq(subclass),
        "cell_id",
    ].drop_duplicates()
    cell_to_position = {
        cell_id: position
        for position, cell_id in enumerate(cell_table["cell_id"])
    }
    positions = np.asarray(
        [cell_to_position[cell_id] for cell_id in cell_ids],
        dtype=int,
    )
    donor_lookup = {
        donor: position for position, donor in enumerate(donors)
    }
    donor_codes = cell_table.iloc[positions]["donor_id"].map(donor_lookup)
    donor_blocks = []
    for donor_index in range(len(donors)):
        donor_blocks.append(
            positions[
                donor_codes.to_numpy(dtype=int) == donor_index
            ]
        )
    return donor_blocks, milo_structure["exposure_by_cell"]


def compute_effects_for_chunk(
    expression: np.ndarray,
    gene_names: list[str],
    cell_table: pd.DataFrame,
    donors: np.ndarray,
    donor_condition: dict[str, str],
    donor_indices: dict[str, list[np.ndarray]],
    milo_structures: dict[str, dict],
    within_blocks: dict[str, dict[str, tuple[list[np.ndarray], np.ndarray]]],
) -> pd.DataFrame:
    rows = []
    donor_condition_vector = np.asarray(
        [donor_condition[donor] for donor in donors]
    )
    for subclass in LABELS:
        block_arrays = donor_indices[subclass]
        subtype_indices = np.concatenate(
            [indices for indices in block_arrays if len(indices) > 0]
        )
        subtype_detection = (
            expression[subtype_indices] > 0
        ).mean(axis=0)
        donor_means = np.full(
            (len(donors), len(gene_names)),
            np.nan,
            dtype=float,
        )
        donor_detect = np.zeros(
            (len(donors), len(gene_names)),
            dtype=bool,
        )
        for donor_index, indices in enumerate(block_arrays):
            if len(indices) == 0:
                continue
            block = expression[indices]
            donor_means[donor_index] = np.nanmean(block, axis=0)
            donor_detect[donor_index] = np.any(block > 0, axis=0)
        for comparison, case_conditions in COMPARISONS.items():
            disease_mask = np.isin(
                donor_condition_vector,
                case_conditions,
            )
            reference_mask = donor_condition_vector == "Reference"
            nhood = milo_structures[comparison]["subtypes"][subclass]
            nhood_means = (
                np.asarray(
                    nhood["weight_matrix"] @ expression,
                    dtype=float,
                )
                if nhood["weight_matrix"] is not None
                else np.empty((0, len(gene_names)))
            )
            nhood_rho = (
                rank_correlation(
                    nhood_means,
                    nhood["nhood_logfc"],
                )
                if nhood_means.shape[0] >= MIN_NHOODS
                else np.full(len(gene_names), np.nan)
            )

            reference_rhos = []
            disease_rhos = []
            for donor_index, indices in enumerate(
                within_blocks[comparison][subclass][0]
            ):
                if len(indices) < MIN_CELLS_PER_DONOR:
                    continue
                exposure = within_blocks[comparison][subclass][1][indices]
                valid = (
                    np.isfinite(exposure)
                    & np.any(np.isfinite(expression[indices]), axis=1)
                )
                donor_indices_valid = indices[valid]
                exposure_valid = exposure[valid]
                if (
                    len(donor_indices_valid) < MIN_CELLS_PER_DONOR
                    or np.unique(exposure_valid).size < 2
                ):
                    continue
                rho = rank_correlation(
                    expression[donor_indices_valid],
                    exposure_valid,
                )
                if donor_condition[donors[donor_index]] == "Reference":
                    reference_rhos.append(rho)
                else:
                    disease_rhos.append(rho)
            if reference_rhos and disease_rhos:
                reference_rho = np.nanmedian(
                    np.vstack(reference_rhos),
                    axis=0,
                )
                disease_rho = np.nanmedian(
                    np.vstack(disease_rhos),
                    axis=0,
                )
                within_delta = disease_rho - reference_rho
            else:
                within_delta = np.full(len(gene_names), np.nan)
            reference_rho_matrix = (
                np.vstack(reference_rhos)
                if len(reference_rhos) >= MIN_DONORS_PER_GROUP
                else None
            )
            disease_rho_matrix = (
                np.vstack(disease_rhos)
                if len(disease_rhos) >= MIN_DONORS_PER_GROUP
                else None
            )
            donor_variance = np.nanvar(
                donor_means[disease_mask | reference_mask],
                axis=0,
            )
            nhood_variance = (
                np.nanvar(nhood_means, axis=0)
                if nhood_means.shape[0]
                else np.zeros(len(gene_names))
            )
            donor_total_detection = donor_detect[
                disease_mask | reference_mask
            ].sum(axis=0)
            base_valid = (
                (subtype_detection >= 0.10)
                & (donor_total_detection >= 10)
                & (donor_variance > 0)
            )
            within_variance_valid = np.zeros(
                len(gene_names),
                dtype=bool,
            )
            if (
                reference_rho_matrix is not None
                and disease_rho_matrix is not None
            ):
                within_variance_valid = (
                    (np.nanvar(reference_rho_matrix, axis=0) > 0)
                    & (np.nanvar(disease_rho_matrix, axis=0) > 0)
                )

            for gene_index, gene in enumerate(gene_names):
                donor_difference_raw = (
                    np.nanmean(
                        donor_means[disease_mask, gene_index]
                    )
                    - np.nanmean(
                        donor_means[reference_mask, gene_index]
                    )
                )
                valid_donor = bool(
                    base_valid[gene_index]
                    and donor_detect[disease_mask, gene_index].sum() >= 5
                    and donor_detect[reference_mask, gene_index].sum() >= 5
                )
                valid_neighborhood = bool(
                    base_valid[gene_index]
                    and nhood_means.shape[0] >= MIN_NHOODS
                    and nhood_variance[gene_index] > 0
                    and np.nanvar(nhood["nhood_logfc"]) > 0
                )
                valid_within = bool(
                    base_valid[gene_index]
                    and within_variance_valid[gene_index]
                )
                donor_difference = (
                    donor_difference_raw if valid_donor else np.nan
                )
                neighborhood_rho = (
                    nhood_rho[gene_index]
                    if valid_neighborhood
                    else np.nan
                )
                within_effect = (
                    within_delta[gene_index]
                    if valid_within
                    else np.nan
                )
                donor_p = np.nan
                if valid_donor:
                    try:
                        donor_p = mannwhitneyu(
                            donor_means[
                                disease_mask,
                                gene_index,
                            ],
                            donor_means[
                                reference_mask,
                                gene_index,
                            ],
                            alternative="two-sided",
                            nan_policy="omit",
                        ).pvalue
                    except ValueError:
                        donor_p = np.nan

                within_p = np.nan
                if valid_within:
                    try:
                        within_p = mannwhitneyu(
                            np.vstack(disease_rhos)[:, gene_index],
                            np.vstack(reference_rhos)[:, gene_index],
                            alternative="two-sided",
                            nan_policy="omit",
                        ).pvalue
                    except ValueError:
                        within_p = np.nan

                rows.append(
                    {
                        "gene": gene,
                        "subtype": subclass,
                        "comparison": comparison,
                        "donor_effect": donor_difference,
                        "donor_p": donor_p,
                        "neighborhood_rho": neighborhood_rho,
                        "within_delta_rho": within_effect,
                        "within_p": within_p,
                        "subtype_detection_fraction": subtype_detection[gene_index],
                        "donor_valid": valid_donor,
                        "neighborhood_valid": valid_neighborhood,
                        "within_valid": valid_within,
                        "donor_detection_disease": int(
                            donor_detect[disease_mask, gene_index].sum()
                        ),
                        "donor_detection_reference": int(
                            donor_detect[reference_mask, gene_index].sum()
                        ),
                    }
                )
    return pd.DataFrame(rows)


def add_standardized_effects(metrics: pd.DataFrame) -> pd.DataFrame:
    metrics = metrics.copy()
    for effect in ["donor_effect", "neighborhood_rho", "within_delta_rho"]:
        metrics[f"{effect}_z"] = np.nan
        for (_, _), group in metrics.groupby(
            ["comparison", "subtype"],
            observed=True,
        ):
            values = group[effect].to_numpy(dtype=float)
            valid = np.isfinite(values)
            if valid.sum() < 10:
                continue
            center = np.nanmean(values[valid])
            scale = np.nanstd(values[valid])
            if not np.isfinite(scale) or scale == 0:
                continue
            metrics.loc[group.index, f"{effect}_z"] = (
                values - center
            ) / scale
    return metrics


def add_reversal_flags(metrics: pd.DataFrame) -> pd.DataFrame:
    metrics = metrics.copy()
    pairs = {
        "donor_within": ("donor_effect_z", "within_delta_rho_z"),
        "donor_neighborhood": (
            "donor_effect_z",
            "neighborhood_rho_z",
        ),
        "neighborhood_within": (
            "neighborhood_rho_z",
            "within_delta_rho_z",
        ),
    }
    for name, (left, right) in pairs.items():
        left_valid = {
            "donor_effect_z": "donor_valid",
            "neighborhood_rho_z": "neighborhood_valid",
            "within_delta_rho_z": "within_valid",
        }[left]
        right_valid = {
            "donor_effect_z": "donor_valid",
            "neighborhood_rho_z": "neighborhood_valid",
            "within_delta_rho_z": "within_valid",
        }[right]
        valid = (
            metrics[[left, right]].notna().all(axis=1)
            & metrics[left_valid]
            & metrics[right_valid]
        )
        metrics[f"{name}_valid"] = valid
        metrics[f"{name}_reversal"] = (
            valid
            & (
                np.sign(metrics[left])
                * np.sign(metrics[right])
                < 0
            )
            & (metrics[left].abs() >= Z_THRESHOLD)
            & (metrics[right].abs() >= Z_THRESHOLD)
        )
    return metrics


def summarize_genes(
    metrics: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for gene, group in metrics.groupby("gene", observed=True):
        row = {
            "gene": gene,
            "n_groups": len(group),
        }
        for pair in [
            "donor_within",
            "donor_neighborhood",
            "neighborhood_within",
        ]:
            valid = group[f"{pair}_valid"]
            count = int(valid.sum())
            row[f"{pair}_valid_groups"] = count
            row[f"{pair}_reversal_fraction"] = (
                float(
                    group.loc[valid, f"{pair}_reversal"].mean()
                )
                if count
                else np.nan
            )
        rows.append(row)
    summary = pd.DataFrame(rows).merge(
        metadata,
        on="gene",
        how="left",
        validate="one_to_one",
    )
    return summary


def empirical_p(observed: float, random: np.ndarray) -> float:
    random = np.asarray(random, dtype=float)
    valid = random[np.isfinite(random)]
    if not np.isfinite(observed) or len(valid) == 0:
        return np.nan
    return (1.0 + np.sum(valid >= observed)) / (len(valid) + 1.0)


def random_set_test(
    gene_summary: pd.DataFrame,
    panel: pd.DataFrame,
    controls: pd.DataFrame,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, dict]:
    summary_by_gene = gene_summary.set_index("gene")
    target_all = panel["gene"].tolist()
    target_biological = panel.loc[
        ~panel["gene_set"].eq("internal_negative_control"),
        "gene",
    ].tolist()
    control_genes = controls["gene"].tolist()
    target_expression = panel["mean_lognorm"].to_numpy(dtype=float)
    expression_breaks = np.quantile(
        target_expression,
        np.linspace(0, 1, 11),
    )
    target_bins = np.clip(
        np.digitize(target_expression, expression_breaks[1:-1]),
        0,
        9,
    )
    target_bin_counts = np.bincount(
        target_bins,
        minlength=10,
    )
    control_bins = np.clip(
        np.digitize(
            controls["mean_lognorm"].to_numpy(dtype=float),
            expression_breaks[1:-1],
        ),
        0,
        9,
    )
    control_bin_members = {
        bin_index: np.asarray(control_genes)[control_bins == bin_index]
        for bin_index in range(10)
    }
    random_rows = []
    for set_index in range(N_RANDOM_SETS):
        set_genes_400 = rng.choice(
            control_genes,
            size=400,
            replace=False,
        )
        set_genes_300 = rng.choice(
            control_genes,
            size=300,
            replace=False,
        )
        row = {"random_set": set_index + 1}
        stratified_sample = np.concatenate(
            [
                rng.choice(
                    control_bin_members[bin_index],
                    size=int(target_bin_counts[bin_index]),
                    replace=False,
                )
                for bin_index in range(10)
            ]
        )
        for pair in [
            "donor_within",
            "donor_neighborhood",
            "neighborhood_within",
        ]:
            column = f"{pair}_reversal_fraction"
            row[f"{pair}_random_400"] = np.nanmean(
                summary_by_gene.loc[set_genes_400, column]
            )
            row[f"{pair}_random_300"] = np.nanmean(
                summary_by_gene.loc[set_genes_300, column]
            )
            row[f"{pair}_stratified_400"] = np.nanmean(
                summary_by_gene.loc[stratified_sample, column]
            )
        random_rows.append(row)
    random_sets = pd.DataFrame(random_rows)

    result = {"target_panel_size": 400, "random_sets": N_RANDOM_SETS}
    for pair in [
        "donor_within",
        "donor_neighborhood",
        "neighborhood_within",
    ]:
        column = f"{pair}_reversal_fraction"
        observed_400 = float(
            summary_by_gene.loc[target_all, column].mean()
        )
        observed_300 = float(
            summary_by_gene.loc[target_biological, column].mean()
        )
        result[pair] = {
            "target_400_mean": observed_400,
            "random_400_mean": float(
                random_sets[f"{pair}_random_400"].mean()
            ),
            "empirical_p_400": empirical_p(
                observed_400,
                random_sets[f"{pair}_random_400"].to_numpy(),
            ),
            "target_300_biological_mean": observed_300,
            "random_300_mean": float(
                random_sets[f"{pair}_random_300"].mean()
            ),
            "empirical_p_300": empirical_p(
                observed_300,
                random_sets[f"{pair}_random_300"].to_numpy(),
            ),
            "stratified_400_mean": float(
                random_sets[f"{pair}_stratified_400"].mean()
            ),
            "empirical_p_stratified_400": empirical_p(
                observed_400,
                random_sets[f"{pair}_stratified_400"].to_numpy(),
            ),
        }
    return random_sets, result


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    panel = pd.read_csv(TARGET_PANEL)
    controls = pd.read_csv(CONTROL_POOL)
    panel_metadata = panel[
        ["gene", "gene_set"]
    ].rename(columns={"gene_set": "panel_group"})
    combined = pd.concat(
        [
            panel_metadata,
            controls[["gene"]].assign(panel_group="external_control_pool"),
        ],
        ignore_index=True,
    )
    if combined["gene"].duplicated().any():
        raise ValueError("Panel and control genes overlap")
    gene_order = combined["gene"].tolist()

    cell_table = pd.read_csv(
        CELL_TABLE,
        dtype={"cell_id": str, "donor_id": str},
        usecols=[
            "cell_id",
            "donor_id",
            "subclass",
            "condition",
            "n_count_rna",
        ],
    )
    cell_table["condition"] = (
        cell_table["condition"]
        .astype(str)
        .map(CONDITION_MAP)
        .fillna("Other")
    )
    adata = ad.read_h5ad(KPMP_H5AD, backed="r")
    target_indices = np.flatnonzero(
        adata.obs["SubclassLevel3"]
        .astype(str)
        .isin(LABELS)
        .to_numpy()
    )
    if not np.array_equal(
        cell_table["cell_id"].to_numpy(),
        adata.obs.index[target_indices].astype(str).to_numpy(),
    ):
        raise ValueError("Cell table order does not match H5AD TAL rows")
    feature_names = adata.var["feature_name"].astype(str).to_numpy()
    symbol_counts = pd.Series(feature_names).value_counts()
    unique_symbols = set(symbol_counts[symbol_counts.eq(1)].index)
    missing = [gene for gene in gene_order if gene not in unique_symbols]
    if missing:
        raise ValueError(f"Panel genes are not unique: {missing[:10]}")
    feature_lookup = {
        gene: int(np.where(feature_names == gene)[0][0])
        for gene in gene_order
    }
    gene_indices = np.asarray(
        [feature_lookup[gene] for gene in gene_order],
        dtype=int,
    )
    expression_cache = build_expression_cache(
        KPMP_H5AD,
        gene_indices,
        target_indices,
        adata.n_obs,
        adata.n_vars,
        cell_table["n_count_rna"].to_numpy(dtype=float),
    )
    expression_matrix = np.load(
        expression_cache,
        mmap_mode="r",
    )
    adata.file.close()

    design = derive_study_design(cell_table)
    donors = design["donors"]
    milo_structures = build_milo_structures(cell_table)
    donor_indices = donor_group_indices(cell_table, donors)
    within_blocks = {}
    for comparison in COMPARISONS:
        within_blocks[comparison] = {}
        for subclass in LABELS:
            within_blocks[comparison][subclass] = within_donor_indices(
                milo_structures[comparison],
                cell_table,
                subclass,
                donors,
            )

    metric_frames = []
    for start in range(0, len(gene_order), GENE_CHUNK_SIZE):
        chunk_genes = gene_order[start : start + GENE_CHUNK_SIZE]
        expression = np.asarray(
            expression_matrix[
                :,
                start : start + GENE_CHUNK_SIZE,
            ],
            dtype=np.float32,
        )
        metrics = compute_effects_for_chunk(
            expression,
            chunk_genes,
            cell_table,
            donors,
            design["donor_condition"],
            donor_indices,
            milo_structures,
            within_blocks,
        )
        metric_frames.append(metrics)
        print(
            f"processed genes {min(start + GENE_CHUNK_SIZE, len(gene_order)):,}"
            f"/{len(gene_order):,}",
            flush=True,
        )

    metrics = pd.concat(metric_frames, ignore_index=True)
    for column in ["donor_p", "within_p"]:
        metrics[f"{column}_fdr"] = bh_adjust(metrics[column])
    metrics = add_standardized_effects(metrics)
    metrics = add_reversal_flags(metrics)
    metrics.to_csv(
        OUT_DIR / "gene_unit_metrics_long.csv",
        index=False,
    )

    gene_summary = summarize_genes(metrics, combined)
    gene_summary.to_csv(
        OUT_DIR / "gene_unit_summary.csv",
        index=False,
    )

    rng = np.random.default_rng(SEED)
    random_sets, primary = random_set_test(
        gene_summary,
        panel,
        controls,
        rng,
    )
    random_sets.to_csv(
        OUT_DIR / "random_set_results.csv",
        index=False,
    )
    panel_sha = json.loads(
        (PANEL_DIR / "panel_manifest.json").read_text(
            encoding="utf-8"
        )
    )["target_panel_sha256"]
    primary["panel_sha256"] = panel_sha
    primary["primary_endpoint"] = "donor_within"
    primary["z_threshold"] = Z_THRESHOLD
    (OUT_DIR / "primary_result.json").write_text(
        json.dumps(primary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(primary, indent=2))


if __name__ == "__main__":
    main()
