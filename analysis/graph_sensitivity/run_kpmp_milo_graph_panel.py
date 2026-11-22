"""Run the frozen KPMP panel on actual Milo graphs with selected k values."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
WORK = ROOT / "work"
ORIGINAL_OUTPUT = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2\outputs"
)
GRAPH_VARIANT_OUTPUT = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "milo"
PROP_VARIANT_OUTPUT = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "milo_prop"
PCA_VARIANT_OUTPUT = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "milo_pca"
CELL_TABLE_PATH = Path(
    r"D:\CodexData\RAMP3\milo_full_tal_input"
    r"\KPMP_milo_full_tal_cells.csv"
)
FOCUS_TABLE_PATH = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2\outputs"
    r"\KPMP_full_TAL_focus_expression.csv.gz"
)
EXPRESSION_CACHE = (
    ROOT
    / "outputs"
    / "prespecified_unit_panel_results"
    / "panel_expression_lognorm.float32.npy"
)
TARGET_PANEL = (
    ROOT / "outputs" / "prespecified_unit_panel" / "target_panel_400.csv"
)
CONTROL_POOL = (
    ROOT / "outputs" / "prespecified_unit_panel" / "negative_control_pool_2000.csv"
)
OUTPUT_ROOT = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "milo_panel"

sys.path.insert(0, str(WORK))
sys.path.insert(0, str(WORK / "unitaware_framework"))

import run_prespecified_unit_panel as panel  # noqa: E402
from unitaware import paired_reversal_permutation  # noqa: E402


VARIANTS = {
    "k15": {
        "directory": GRAPH_VARIANT_OUTPUT,
        "prefix": "KPMP_milo_graph_k15",
    },
    "k30": {
        "directory": ORIGINAL_OUTPUT,
        "prefix": "KPMP_milo_full_TAL_final",
    },
    "k30_seed": {
        "directory": GRAPH_VARIANT_OUTPUT,
        "prefix": "KPMP_milo_graph_k30_seed20261122",
    },
    "k50": {
        "directory": GRAPH_VARIANT_OUTPUT,
        "prefix": "KPMP_milo_graph_k50",
    },
    "prop005": {
        "directory": PROP_VARIANT_OUTPUT,
        "prefix": "KPMP_milo_graph_k30_prop005",
    },
    "prop020": {
        "directory": PROP_VARIANT_OUTPUT,
        "prefix": "KPMP_milo_graph_k30_prop020",
    },
    "pca30": {
        "directory": PCA_VARIANT_OUTPUT,
        "prefix": "KPMP_milo_graph_pca30_k30",
    },
}


def load_cell_table() -> pd.DataFrame:
    cells = pd.read_csv(
        CELL_TABLE_PATH,
        dtype={"cell_id": str, "donor_id": str},
        usecols=[
            "cell_id",
            "donor_id",
            "condition",
            "subclass",
            "region",
            "UMAP1",
            "UMAP2",
        ],
    )
    focus = pd.read_csv(
        FOCUS_TABLE_PATH,
        dtype={"cell_id": str},
        usecols=["cell_id", "n_count_rna", "PAPPA2_lognorm"],
    )
    if not np.array_equal(
        cells["cell_id"].to_numpy(),
        focus["cell_id"].to_numpy(),
    ):
        raise ValueError("Graph metadata and focus table are not aligned")
    cells["n_count_rna"] = focus["n_count_rna"].to_numpy(dtype=float)
    cells["PAPPA2"] = focus["PAPPA2_lognorm"].to_numpy(dtype=float)
    cells["condition"] = (
        cells["condition"]
        .astype(str)
        .map(panel.CONDITION_MAP)
        .fillna("Other")
    )
    return cells


def attach_case_z(
    case_metrics: pd.DataFrame,
    panel_metrics: pd.DataFrame,
) -> pd.DataFrame:
    scales = (
        panel_metrics.groupby(["comparison", "subtype"], observed=True)[
            ["donor_effect", "within_delta_rho", "neighborhood_rho"]
        ]
        .std(ddof=0)
        .reset_index()
        .rename(
            columns={
                "donor_effect": "scale_donor_effect",
                "within_delta_rho": "scale_within_delta_rho",
                "neighborhood_rho": "scale_neighborhood_rho",
            }
        )
    )
    case_metrics = case_metrics.merge(
        scales,
        on=["comparison", "subtype"],
        how="left",
        validate="many_to_one",
    )
    for effect in [
        "donor_effect",
        "within_delta_rho",
        "neighborhood_rho",
    ]:
        values = case_metrics[effect].to_numpy(dtype=float)
        scale = case_metrics[f"scale_{effect}"].to_numpy(dtype=float)
        case_metrics[f"{effect}_z"] = np.divide(
            values,
            scale,
            out=np.full(len(case_metrics), np.nan),
            where=np.isfinite(values)
            & np.isfinite(scale)
            & (scale > 0),
        )
    case_metrics["donor_within_reversal"] = (
        case_metrics[["donor_effect_z", "within_delta_rho_z"]]
        .notna()
        .all(axis=1)
        & (
            np.sign(case_metrics["donor_effect_z"])
            * np.sign(case_metrics["within_delta_rho_z"])
            < 0
        )
        & (case_metrics["donor_effect_z"].abs() >= panel.Z_THRESHOLD)
        & (case_metrics["within_delta_rho_z"].abs() >= panel.Z_THRESHOLD)
    )
    return case_metrics


def run_variant(
    variant: str,
    cells: pd.DataFrame,
    panel_table: pd.DataFrame,
    controls: pd.DataFrame,
    expression: np.ndarray,
) -> dict:
    settings = VARIANTS[variant]
    panel.MILO_DIR = settings["directory"]
    panel.MILO_PREFIX = settings["prefix"]
    if not (
        settings["directory"]
        / f"{settings['prefix']}_all_disease_nhood_results.csv"
    ).exists():
        raise FileNotFoundError(
            f"Missing Milo results for {variant}: {settings['directory']}"
        )

    print(
        f"Running actual Milo panel variant {variant} "
        f"from {settings['prefix']}",
        flush=True,
    )
    structures = panel.build_milo_structures(cells)
    design = panel.derive_study_design(cells)
    donors = design["donors"]
    donor_indices = panel.donor_group_indices(cells, donors)
    within_blocks = {
        comparison: {
            subclass: panel.within_donor_indices(
                structures[comparison],
                cells,
                subclass,
                donors,
            )
            for subclass in panel.LABELS
        }
        for comparison in panel.COMPARISONS
    }

    gene_order = (
        panel_table["gene"].tolist()
        + controls["gene"].tolist()
    )
    metric_parts = []
    for start in range(0, len(gene_order), panel.GENE_CHUNK_SIZE):
        stop = min(start + panel.GENE_CHUNK_SIZE, len(gene_order))
        expression_chunk = np.asarray(
            expression[:, start:stop],
            dtype=np.float32,
        )
        metric_parts.append(
            panel.compute_effects_for_chunk(
                expression_chunk,
                gene_order[start:stop],
                cells,
                donors,
                design["donor_condition"],
                donor_indices,
                structures,
                within_blocks,
            )
        )
        print(
            f"{variant}: processed genes {stop}/{len(gene_order)}",
            flush=True,
        )
    metrics = pd.concat(metric_parts, ignore_index=True)
    for column in ["donor_p", "within_p"]:
        metrics[f"{column}_fdr"] = panel.bh_adjust(metrics[column])
    metrics = panel.add_standardized_effects(metrics)
    metrics = panel.add_reversal_flags(metrics)

    output_dir = OUTPUT_ROOT / variant
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(
        output_dir / "gene_unit_metrics_long.csv",
        index=False,
    )

    combined = pd.concat(
        [
            panel_table[["gene", "gene_set"]].rename(
                columns={"gene_set": "panel_group"}
            ),
            controls[["gene"]].assign(
                panel_group="external_control_pool"
            ),
        ],
        ignore_index=True,
    )
    summary = panel.summarize_genes(metrics, combined)
    summary.to_csv(
        output_dir / "gene_unit_summary.csv",
        index=False,
    )
    random_sets, random_summary = panel.random_set_test(
        summary,
        panel_table,
        controls,
        np.random.default_rng(panel.SEED + len(variant)),
    )
    random_sets.to_csv(
        output_dir / "random_set_results.csv",
        index=False,
    )
    target_genes = set(panel_table["gene"])
    null_table, pairing_summary = paired_reversal_permutation(
        metrics,
        left_col="donor_effect_z",
        right_col="within_delta_rho_z",
        target_genes=target_genes,
        gene_col="gene",
        strata_cols=("subtype", "comparison"),
        z_threshold=panel.Z_THRESHOLD,
        n_permutations=panel.N_RANDOM_SETS * 10,
        seed=panel.SEED + 1000 + len(variant),
        aggregate="gene",
    )
    null_table.to_csv(
        output_dir / "effect_aware_null.csv",
        index=False,
    )

    pappa2_raw = panel.compute_effects_for_chunk(
        cells["PAPPA2"].to_numpy(dtype=np.float32)[:, None],
        ["PAPPA2"],
        cells,
        donors,
        design["donor_condition"],
        donor_indices,
        structures,
        within_blocks,
    )
    pappa2_raw = attach_case_z(pappa2_raw, metrics)
    pappa2_raw.insert(0, "variant", variant)
    pappa2_raw.to_csv(
        output_dir / "PAPPA2_case_metrics.csv",
        index=False,
    )

    result = {
        "variant": variant,
        "milo_prefix": settings["prefix"],
        "target_panel_size": 400,
        "cells": int(len(cells)),
        "donor_within": random_summary["donor_within"],
        "effect_aware_null_mean": float(pairing_summary["null_mean"]),
        "effect_aware_null_q025": float(pairing_summary["null_q025"]),
        "effect_aware_null_q975": float(pairing_summary["null_q975"]),
        "effect_aware_empirical_p": float(pairing_summary["empirical_p"]),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2), flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--variants",
        default="k15,k30,k50",
        help="Comma-separated variants from k15, k30, k50.",
    )
    args = parser.parse_args()
    variants = [
        value.strip()
        for value in args.variants.split(",")
        if value.strip()
    ]
    unknown = sorted(set(variants) - set(VARIANTS))
    if unknown:
        raise ValueError(f"Unknown variants: {unknown}")

    cells = load_cell_table()
    panel_table = pd.read_csv(TARGET_PANEL)
    controls = pd.read_csv(CONTROL_POOL)
    expression = np.load(EXPRESSION_CACHE, mmap_mode="r")
    expected_genes = len(panel_table) + len(controls)
    if expression.shape != (len(cells), expected_genes):
        raise ValueError("Expression cache shape is not aligned")

    rows = []
    for variant in variants:
        rows.append(
            run_variant(
                variant,
                cells,
                panel_table,
                controls,
                expression,
            )
        )
    pd.DataFrame(rows).to_csv(
        OUTPUT_ROOT / "milo_panel_summary.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
