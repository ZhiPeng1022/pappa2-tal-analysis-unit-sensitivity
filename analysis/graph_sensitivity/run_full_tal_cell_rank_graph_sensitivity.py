"""Compare full-TAL within-donor cell-rank effects across Milo graphs."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
ORIGINAL_REPO = Path(r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2")
ORIGINAL_OUTPUT = ORIGINAL_REPO / "outputs"
GRAPH_OUTPUT = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "milo"
PROP_OUTPUT = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "milo_prop"
PCA_OUTPUT = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "milo_pca"
FOCUS_EXPRESSION = (
    ORIGINAL_OUTPUT / "KPMP_full_TAL_focus_expression.csv.gz"
)
OUTPUT_ROOT = (
    ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "cell_rank"
)

sys.path.insert(0, str(ORIGINAL_REPO / "repo" / "analysis"))
sys.path.insert(0, str(ORIGINAL_REPO / "repo"))

import analyze_full_tal_within_donor_cell_rank as rank  # noqa: E402
import analyze_full_tal_within_donor_nhood_effect as nhood  # noqa: E402


VARIANTS = {
    "k15": {
        "directory": GRAPH_OUTPUT,
        "prefix": "KPMP_milo_graph_k15",
    },
    "k30": {
        "directory": ORIGINAL_OUTPUT,
        "prefix": "KPMP_milo_full_TAL_final",
    },
    "k30_seed": {
        "directory": GRAPH_OUTPUT,
        "prefix": "KPMP_milo_graph_k30_seed20261122",
    },
    "k50": {
        "directory": GRAPH_OUTPUT,
        "prefix": "KPMP_milo_graph_k50",
    },
    "prop005": {
        "directory": PROP_OUTPUT,
        "prefix": "KPMP_milo_graph_k30_prop005",
    },
    "prop020": {
        "directory": PROP_OUTPUT,
        "prefix": "KPMP_milo_graph_k30_prop020",
    },
    "pca30": {
        "directory": PCA_OUTPUT,
        "prefix": "KPMP_milo_graph_pca30_k30",
    },
}
GENE = "PAPPA2"
PERMUTATIONS = 1000
PERMUTATION_SEED = 20261122


def label_permutation_p(
    blocks: list[dict],
    *,
    gene_index: int,
    seed: int,
) -> float:
    observed = np.asarray(
        [
            block["observed_rho"][gene_index]
            for block in blocks
        ],
        dtype=float,
    )
    reference = np.asarray(
        [
            block["condition"] == "Reference"
            for block in blocks
        ],
        dtype=bool,
    )
    reference_count = int(reference.sum())
    disease_count = len(reference) - reference_count
    if reference_count == 0 or disease_count == 0:
        return np.nan
    observed_difference = (
        np.nanmedian(observed[~reference])
        - np.nanmedian(observed[reference])
    )
    rng = np.random.default_rng(seed)
    null = np.empty(PERMUTATIONS, dtype=float)
    for index in range(PERMUTATIONS):
        shuffled = rng.permutation(reference)
        null[index] = (
            np.nanmedian(observed[~shuffled])
            - np.nanmedian(observed[shuffled])
        )
    return float(
        (1.0 + np.sum(np.abs(null) >= abs(observed_difference)))
        / (PERMUTATIONS + 1.0)
    )


def load_expression() -> pd.DataFrame:
    frame = pd.read_csv(
        FOCUS_EXPRESSION,
        dtype={"cell_id": str, "donor_id": str},
    )
    if not frame["cell_id"].is_unique:
        raise ValueError("Full TAL focus expression has duplicate cell_id")
    columns = {gene: f"{gene}_lognorm" for gene in nhood.FOCUS_GENES}
    missing = [
        column
        for column in columns.values()
        if column not in frame.columns
    ]
    if missing:
        raise ValueError(f"Missing focus columns: {missing}")
    expression = frame.set_index("cell_id")[
        list(columns.values())
    ].rename(
        columns={column: gene for gene, column in columns.items()}
    )
    expression.index.name = "cell_id"
    return expression


def run_variant(
    variant: str,
    expression: pd.DataFrame,
) -> pd.DataFrame:
    settings = VARIANTS[variant]
    nhood.OUTPUT_ROOT = settings["directory"]
    nhood.MILO_PREFIX = settings["prefix"]
    print(
        f"Reading full-TAL Milo edges for {variant}",
        flush=True,
    )
    edges_by_comparison = []
    for comparison in nhood.COMPARISONS:
        edges_by_comparison.append(
            nhood.read_and_merge_edges(
                comparison,
                expression,
            )
        )
    edges = pd.concat(edges_by_comparison, ignore_index=True)

    rows = []
    for comparison in nhood.COMPARISONS:
        comparison_edges = edges[
            edges["comparison"].eq(comparison)
        ]
        for subclass in sorted(
            comparison_edges["nhood_subclass"].dropna().unique()
        ):
            group_edges = comparison_edges[
                comparison_edges["nhood_subclass"].eq(subclass)
            ]
            cell_data = rank.build_cell_table(group_edges)
            blocks = rank.prepare_donor_blocks(cell_data)
            reference = np.asarray(
                [
                    block["observed_rho"][
                        nhood.FOCUS_GENES.index(GENE)
                    ]
                    for block in blocks
                    if block["condition"] == "Reference"
                ],
                dtype=float,
            )
            disease = np.asarray(
                [
                    block["observed_rho"][
                        nhood.FOCUS_GENES.index(GENE)
                    ]
                    for block in blocks
                    if block["condition"] != "Reference"
                ],
                dtype=float,
            )
            reference = reference[np.isfinite(reference)]
            disease = disease[np.isfinite(disease)]
            if len(reference) and len(disease):
                difference = (
                    float(np.median(disease))
                    - float(np.median(reference))
                )
                p_value = mannwhitneyu(
                    disease,
                    reference,
                    alternative="two-sided",
                ).pvalue
            else:
                difference = np.nan
                p_value = np.nan
            rows.append(
                {
                    "variant": variant,
                    "gene": GENE,
                    "comparison": comparison,
                    "subclass": subclass,
                    "cells_total": int(len(cell_data)),
                    "donors_total": int(
                        cell_data["donor_id"].nunique()
                    ),
                    "eligible_reference_donors": int(len(reference)),
                    "eligible_disease_donors": int(len(disease)),
                    "median_rho_reference": (
                        float(np.median(reference))
                        if len(reference)
                        else np.nan
                    ),
                    "median_rho_disease": (
                        float(np.median(disease))
                        if len(disease)
                        else np.nan
                    ),
                    "disease_minus_reference": difference,
                    "donor_mannwhitney_p": float(p_value),
                    "label_permutation_p": label_permutation_p(
                        blocks,
                        gene_index=nhood.FOCUS_GENES.index(GENE),
                        seed=(
                            PERMUTATION_SEED
                            + len(variant)
                            + len(comparison)
                            + len(str(subclass))
                        ),
                    ),
                }
            )
            print(
                f"{variant} {comparison} {subclass}: "
                f"delta={difference:.4f}, "
                f"n={len(reference)}/{len(disease)}",
                flush=True,
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    expression = load_expression()
    rows = []
    for variant in VARIANTS:
        result = run_variant(variant, expression)
        result.to_csv(
            OUTPUT_ROOT / f"full_cell_rank_{variant}.csv",
            index=False,
        )
        rows.append(result)
    combined = pd.concat(rows, ignore_index=True)
    combined.to_csv(
        OUTPUT_ROOT / "full_cell_rank_graph_sensitivity.csv",
        index=False,
    )
    summary = (
        combined[
            combined["subclass"].eq("aTAL2")
            & combined["comparison"].eq("all_disease")
        ][
            [
                "variant",
                "median_rho_reference",
                "median_rho_disease",
                "disease_minus_reference",
                "donor_mannwhitney_p",
                "label_permutation_p",
            ]
        ]
    )
    summary.to_csv(
        OUTPUT_ROOT / "aTAL2_PAPPA2_graph_sensitivity.csv",
        index=False,
    )
    (OUTPUT_ROOT / "manifest.json").write_text(
        json.dumps(
            {
                "gene": GENE,
                "variants": list(VARIANTS),
                "comparisons": nhood.COMPARISONS,
                "statistic": (
                    "median donor within-donor cell-rank Spearman rho "
                    "in disease minus reference"
                ),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
