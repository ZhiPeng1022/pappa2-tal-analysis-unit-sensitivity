"""Build provenance summaries for the actual-Milo graph sensitivity runs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
ORIGINAL_OUTPUT = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2\outputs"
)
GRAPH_OUTPUT = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "milo"
PROP_OUTPUT = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "milo_prop"
PCA_OUTPUT = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "milo_pca"
OUTPUT = (
    ROOT
    / "outputs"
    / "kpmp_graph_sensitivity"
    / "v1"
    / "milo_provenance.csv"
)

VARIANTS = {
    "k15": {
        "directory": GRAPH_OUTPUT,
        "prefix": "KPMP_milo_graph_k15",
        "k": 15,
        "prop": 0.1,
        "seed": 20261122,
        "embedding": "UMAP2D",
    },
    "k30_original": {
        "directory": ORIGINAL_OUTPUT,
        "prefix": "KPMP_milo_full_TAL_final",
        "k": 30,
        "prop": 0.1,
        "seed": 20260927,
        "embedding": "UMAP2D",
    },
    "k30_seed": {
        "directory": GRAPH_OUTPUT,
        "prefix": "KPMP_milo_graph_k30_seed20261122",
        "k": 30,
        "prop": 0.1,
        "seed": 20261122,
        "embedding": "UMAP2D",
    },
    "k50": {
        "directory": GRAPH_OUTPUT,
        "prefix": "KPMP_milo_graph_k50",
        "k": 50,
        "prop": 0.1,
        "seed": 20261122,
        "embedding": "UMAP2D",
    },
    "prop005": {
        "directory": PROP_OUTPUT,
        "prefix": "KPMP_milo_graph_k30_prop005",
        "k": 30,
        "prop": 0.05,
        "seed": 20261122,
        "embedding": "UMAP2D",
    },
    "prop020": {
        "directory": PROP_OUTPUT,
        "prefix": "KPMP_milo_graph_k30_prop020",
        "k": 30,
        "prop": 0.20,
        "seed": 20261122,
        "embedding": "UMAP2D",
    },
    "pca30": {
        "directory": PCA_OUTPUT,
        "prefix": "KPMP_milo_graph_pca30_k30",
        "k": 30,
        "prop": 0.10,
        "seed": 20261122,
        "embedding": "panel_PCA30",
    },
}
COMPARISONS = [
    "all_disease",
    "AKI_vs_Reference",
    "CKD_vs_Reference",
]


def main() -> None:
    rows = []
    for variant, settings in VARIANTS.items():
        directory = settings["directory"]
        prefix = settings["prefix"]
        for comparison in COMPARISONS:
            result_path = directory / (
                f"{prefix}_{comparison}_nhood_results.csv"
            )
            membership_path = directory / (
                f"{prefix}_{comparison}_nhood_membership.csv.gz"
            )
            if not result_path.exists() or not membership_path.exists():
                raise FileNotFoundError(
                    f"Missing Milo inputs for {variant}: {comparison}"
                )
            results = pd.read_csv(result_path)
            membership = pd.read_csv(
                membership_path,
                usecols=["cell_id", "donor_id"],
                dtype={"cell_id": str, "donor_id": str},
            )
            significant = results[
                results["SpatialFDR"].astype(float) < 0.05
            ]
            rows.append(
                {
                    "variant": variant,
                    "embedding": settings["embedding"],
                    "k": settings["k"],
                    "prop": settings["prop"],
                    "seed": settings["seed"],
                    "comparison": comparison,
                    "cells": int(membership["cell_id"].nunique()),
                    "donors": int(membership["donor_id"].nunique()),
                    "neighborhoods": int(results["Nhood"].nunique()),
                    "significant_neighborhoods": int(len(significant)),
                    "positive_significant_neighborhoods": int(
                        (significant["logFC"].astype(float) > 0).sum()
                    ),
                    "negative_significant_neighborhoods": int(
                        (significant["logFC"].astype(float) < 0).sum()
                    ),
                }
            )
    provenance = pd.DataFrame(rows)
    provenance.to_csv(OUTPUT, index=False)

    report = [
        "# Milo run provenance",
        "",
        "All graph-stability runs used actual miloR neighborhoods and logMS "
        "normalization. The k/seed sensitivity uses UMAP2D with `prop = 0.1`; "
        "the fixed `k = 30` proportion sensitivity adds `prop = 0.05` and "
        "`0.20`; the embedding sensitivity uses panel PCA30. The k = 30 "
        "original UMAP run is the prespecified primary analysis; the other "
        "runs are post hoc sensitivity settings.",
        "",
        provenance.to_string(index=False),
        "",
    ]
    OUTPUT.with_suffix(".md").write_text(
        "\n".join(report),
        encoding="utf-8",
    )
    print(provenance.to_string(index=False))


if __name__ == "__main__":
    main()
