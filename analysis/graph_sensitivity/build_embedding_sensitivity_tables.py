"""Build UMAP-versus-PCA embedding sensitivity tables."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
GRAPH = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1"
PANEL_ROOT = GRAPH / "milo_panel"
CELL_RANK = GRAPH / "cell_rank" / "full_cell_rank_graph_sensitivity.csv"
OUTPUT_PANEL = GRAPH / "panel_calibration_by_embedding.csv"
OUTPUT_CELL = GRAPH / "aTAL2_PAPPA2_by_embedding.csv"
OUTPUT_RUNS = GRAPH / "embedding_sensitivity_runs.csv"
OUTPUT_STABILITY = GRAPH / "embedding_sensitivity_stability.csv"
OUTPUT_REPORT = GRAPH / "embedding_sensitivity_report.md"

sys.path.insert(0, str(ROOT / "work" / "unitaware_framework"))
from unitaware import summarize_calibration_stability  # noqa: E402


VARIANTS = {
    "umap_k30_seed": {
        "panel_variant": "k30_seed",
        "cell_variant": "k30_seed",
        "embedding": "UMAP2D",
    },
    "pca30": {
        "panel_variant": "pca30",
        "cell_variant": "pca30",
        "embedding": "panel_PCA30",
    },
}


def main() -> None:
    panel_rows = []
    for embedding_variant, settings in VARIANTS.items():
        path = (
            PANEL_ROOT
            / settings["panel_variant"]
            / "summary.json"
        )
        payload = json.loads(path.read_text(encoding="utf-8"))
        endpoint = payload["donor_within"]
        panel_rows.append(
            {
                "embedding_variant": embedding_variant,
                "embedding": settings["embedding"],
                "target_reversal_fraction": endpoint["target_400_mean"],
                "expression_matched_mean": endpoint["random_400_mean"],
                "expression_matched_p": endpoint["empirical_p_400"],
                "effect_aware_null_mean": payload[
                    "effect_aware_null_mean"
                ],
                "effect_aware_p": payload[
                    "effect_aware_empirical_p"
                ],
                "target_minus_effect_aware": (
                    endpoint["target_400_mean"]
                    - payload["effect_aware_null_mean"]
                ),
            }
        )
    panel = pd.DataFrame(panel_rows)
    panel.to_csv(OUTPUT_PANEL, index=False)

    cell_rank = pd.read_csv(CELL_RANK)
    cell_rank = cell_rank[
        cell_rank["variant"].isin(
            [settings["cell_variant"] for settings in VARIANTS.values()]
        )
        & cell_rank["comparison"].eq("all_disease")
        & cell_rank["subclass"].eq("aTAL2")
    ].copy()
    cell_rank["embedding_variant"] = cell_rank["variant"].map(
        {
            settings["cell_variant"]: embedding_variant
            for embedding_variant, settings in VARIANTS.items()
        }
    )
    cell_rank["embedding"] = cell_rank["embedding_variant"].map(
        {
            embedding_variant: settings["embedding"]
            for embedding_variant, settings in VARIANTS.items()
        }
    )
    cell_rank.to_csv(OUTPUT_CELL, index=False)

    runs = []
    for row in panel.itertuples(index=False):
        runs.append(
            {
                "endpoint": "expression_matched_bias",
                "embedding_variant": row.embedding_variant,
                "effect": (
                    row.target_reversal_fraction
                    - row.expression_matched_mean
                ),
                "p_value": row.expression_matched_p,
            }
        )
        runs.append(
            {
                "endpoint": "donor_conditioned_panel_residual",
                "embedding_variant": row.embedding_variant,
                "effect": row.target_minus_effect_aware,
                "p_value": row.effect_aware_p,
            }
        )
    for row in cell_rank.itertuples(index=False):
        runs.append(
            {
                "endpoint": "aTAL2_PAPPA2_cell_rank_conflict",
                "embedding_variant": row.embedding_variant,
                "effect": row.disease_minus_reference,
                "p_value": row.label_permutation_p,
            }
        )
    runs = pd.DataFrame(runs)
    runs.to_csv(OUTPUT_RUNS, index=False)

    stability = summarize_calibration_stability(
        runs,
        group_cols=("endpoint",),
        effect_col="effect",
        p_col="p_value",
    )
    stability.to_csv(OUTPUT_STABILITY, index=False)

    report = [
        "# Embedding-space sensitivity",
        "",
        "Both runs fix k = 30, prop = 0.10, and seed 20261122. Only the "
        "neighborhood representation changes: atlas UMAP versus PCA30 built "
        "from the frozen 2,400-gene panel.",
        "",
        "## Panel calibration",
        "",
        panel.to_string(index=False),
        "",
        "## Stability summary",
        "",
        stability.to_string(index=False),
        "",
        "## aTAL2/PAPPA2 cell-rank",
        "",
        cell_rank[
            [
                "embedding_variant",
                "embedding",
                "disease_minus_reference",
                "donor_mannwhitney_p",
                "label_permutation_p",
            ]
        ].to_string(index=False),
        "",
    ]
    OUTPUT_REPORT.write_text("\n".join(report), encoding="utf-8")
    print(stability.to_string(index=False))


if __name__ == "__main__":
    main()
