"""Build a PCA30 representation from the frozen 2,400-gene TAL panel."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import IncrementalPCA


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
CELL_TABLE = Path(
    r"D:\CodexData\RAMP3\milo_full_tal_input"
    r"\KPMP_milo_full_tal_cells.csv"
)
EXPRESSION = (
    ROOT
    / "outputs"
    / "prespecified_unit_panel_results"
    / "panel_expression_lognorm.float32.npy"
)
OUTPUT_DIR = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1" / "pca_input"
OUTPUT = OUTPUT_DIR / "KPMP_milo_panel_pca30.csv"
SUMMARY = OUTPUT_DIR / "KPMP_milo_panel_pca30_summary.json"
CHUNK_SIZE = 10000
N_COMPONENTS = 30
SEED = 20261122


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    cells = pd.read_csv(
        CELL_TABLE,
        dtype={"cell_id": str, "donor_id": str},
    )
    expression = np.load(EXPRESSION, mmap_mode="r")
    if expression.shape[0] != len(cells):
        raise ValueError("Cell table and expression cache are not aligned")

    n_cells, n_genes = expression.shape
    total = np.zeros(n_genes, dtype=np.float64)
    total_sq = np.zeros(n_genes, dtype=np.float64)
    for start in range(0, n_cells, CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, n_cells)
        block = np.asarray(expression[start:stop], dtype=np.float64)
        total += block.sum(axis=0)
        total_sq += np.square(block).sum(axis=0)
    mean = total / n_cells
    variance = np.maximum(total_sq / n_cells - np.square(mean), 0.0)
    sd = np.sqrt(variance)
    usable = sd > 1e-8
    if usable.sum() < N_COMPONENTS:
        raise ValueError("Too few variable genes for PCA")

    pca = IncrementalPCA(
        n_components=N_COMPONENTS,
        batch_size=CHUNK_SIZE,
    )
    for start in range(0, n_cells, CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, n_cells)
        block = np.asarray(expression[start:stop], dtype=np.float32)
        block = (block - mean.astype(np.float32)) / sd.astype(np.float32)
        pca.partial_fit(block[:, usable])
        if start % 20000 == 0:
            print(f"PCA fit {stop:,}/{n_cells:,}", flush=True)

    scores = np.empty((n_cells, N_COMPONENTS), dtype=np.float32)
    for start in range(0, n_cells, CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, n_cells)
        block = np.asarray(expression[start:stop], dtype=np.float32)
        block = (block - mean.astype(np.float32)) / sd.astype(np.float32)
        scores[start:stop] = pca.transform(
            block[:, usable]
        ).astype(np.float32)

    pca_columns = [f"PCA{index}" for index in range(1, N_COMPONENTS + 1)]
    result = pd.concat(
        [
            cells.reset_index(drop=True),
            pd.DataFrame(scores, columns=pca_columns),
        ],
        axis=1,
    )
    result.to_csv(OUTPUT, index=False)
    summary = {
        "n_cells": int(n_cells),
        "n_genes": int(n_genes),
        "usable_genes": int(usable.sum()),
        "n_components": int(N_COMPONENTS),
        "seed_not_used_for_deterministic_pca": SEED,
        "explained_variance_ratio": [
            float(value)
            for value in pca.explained_variance_ratio_
        ],
        "output": str(OUTPUT),
    }
    SUMMARY.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(OUTPUT)


if __name__ == "__main__":
    main()
