import anndata as ad
import h5py
import numpy as np
import pandas as pd

from project_paths import DATA_ROOT, OUTPUT_ROOT

INPUT_H5AD = DATA_ROOT / "KPMP_new_snRNA.h5ad"

TARGET_SUBCLASSES = [
    "C-TAL-A",
    "C-TAL-B",
    "aTAL2",
    "frTAL",
]
FOCUS_GENES = [
    "PAPPA2",
    "CXCL12",
    "NOTCH3",
    "HES1",
    "EPHB2",
]

CONDITION_MAP = {
    "Healthy Reference": "Reference",
    "Reference": "Reference",
    "Acute Kidney Injury": "AKI",
    "AKI": "AKI",
    "Chronic Kidney Disease": "CKD",
    "CKD": "CKD",
    "Diabetes Mellitus - Resilient": "Resilient",
}

EXPRESSION_OUTPUT = (
    OUTPUT_ROOT / "KPMP_full_TAL_focus_expression.csv.gz"
)
MANIFEST_OUTPUT = (
    OUTPUT_ROOT / "KPMP_full_TAL_focus_expression_manifest.csv"
)


def resolve_gene_indices(adata):
    feature_names = (
        adata.var["feature_name"].astype(str).to_numpy()
    )
    gene_indices = {}
    for gene in FOCUS_GENES:
        matches = np.flatnonzero(feature_names == gene)
        if len(matches) != 1:
            raise ValueError(
                f"{gene} matched {len(matches)} feature rows."
            )
        gene_indices[gene] = int(matches[0])
    return gene_indices


def read_focus_counts(input_path, row_indices, gene_indices):
    count_matrix = np.zeros(
        (len(row_indices), len(FOCUS_GENES)),
        dtype=np.float32,
    )
    focus_columns = [
        gene_indices[gene] for gene in FOCUS_GENES
    ]
    focus_lookup = {
        gene_index: column
        for column, gene_index in enumerate(focus_columns)
    }

    with h5py.File(input_path, "r") as handle:
        indptr = handle["X/indptr"][:]
        index_dataset = handle["X/indices"]
        data_dataset = handle["X/data"]

        for output_row, input_row in enumerate(row_indices):
            start = int(indptr[input_row])
            end = int(indptr[input_row + 1])
            row_indices_data = index_dataset[start:end]
            row_data = data_dataset[start:end]

            for position, gene_index in enumerate(
                row_indices_data
            ):
                column = focus_lookup.get(int(gene_index))
                if column is not None:
                    count_matrix[output_row, column] = (
                        row_data[position]
                    )

            if (output_row + 1) % 10000 == 0:
                print(
                    f"Read {output_row + 1:,} / "
                    f"{len(row_indices):,} cells",
                    flush=True,
                )

    return count_matrix


def build_metadata(adata, target_rows):
    obs = adata.obs.iloc[target_rows].copy()
    condition = (
        obs["ConditionCategory"]
        .astype(str)
        .map(CONDITION_MAP)
        .fillna("Other")
    )
    metadata = pd.DataFrame(
        {
            "cell_id": obs.index.astype(str),
            "donor_id": obs["donor_id"].astype(str).to_numpy(),
            "subclass": (
                obs["SubclassLevel3"].astype(str).to_numpy()
            ),
            "condition": condition.to_numpy(),
            "cell_state_level1": (
                obs["CellStateLevel1"].astype(str).to_numpy()
            ),
            "cell_state_level2": (
                obs["CellStateLevel2"].astype(str).to_numpy()
            ),
            "region": obs["region"].astype(str).to_numpy(),
            "age": obs["Age"].astype(str).to_numpy(),
            "sex": obs["sex"].astype(str).to_numpy(),
            "disease": obs["disease"].astype(str).to_numpy(),
            "egfr_bin": (
                obs[
                    "Baseline eGFR (ml/min/1.73m2) (Binned)"
                ]
                .astype(str)
                .to_numpy()
            ),
            "n_count_rna": pd.to_numeric(
                obs["nCount_RNA"],
                errors="coerce",
            ).to_numpy(dtype=float),
            "n_feature_rna": pd.to_numeric(
                obs["nFeature_RNA"],
                errors="coerce",
            ).to_numpy(dtype=float),
        }
    )
    return metadata


def normalize_focus_counts(metadata, count_matrix):
    library_size = metadata[
        "n_count_rna"
    ].to_numpy(dtype=float)
    scale = 10000.0 / np.maximum(library_size, 1.0)
    return np.log1p(count_matrix * scale[:, None])


def build_manifest(expression):
    rows = []
    for (subclass, condition), group in expression.groupby(
        ["subclass", "condition"],
        observed=True,
    ):
        rows.append(
            {
                "subclass": subclass,
                "condition": condition,
                "cells": len(group),
                "donors": group["donor_id"].nunique(),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["subclass", "condition"]
    )


def main():
    if not INPUT_H5AD.exists():
        raise FileNotFoundError(INPUT_H5AD)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(INPUT_H5AD, backed="r")
    subclass_values = (
        adata.obs["SubclassLevel3"].astype(str).to_numpy()
    )
    target_rows = np.flatnonzero(
        np.isin(subclass_values, TARGET_SUBCLASSES)
    )
    gene_indices = resolve_gene_indices(adata)

    print(
        f"Extracting {len(target_rows):,} TAL cells "
        f"for {len(FOCUS_GENES)} focus genes",
        flush=True,
    )
    print("Gene indices:", gene_indices, flush=True)

    count_matrix = read_focus_counts(
        INPUT_H5AD,
        target_rows,
        gene_indices,
    )
    metadata = build_metadata(adata, target_rows)
    normalized = normalize_focus_counts(
        metadata,
        count_matrix,
    )

    expression = metadata.copy()
    for gene_index, gene in enumerate(FOCUS_GENES):
        expression[f"{gene}_count"] = (
            count_matrix[:, gene_index]
        )
        expression[f"{gene}_lognorm"] = (
            normalized[:, gene_index]
        )

    expression.to_csv(
        EXPRESSION_OUTPUT,
        index=False,
        compression="gzip",
    )
    manifest = build_manifest(expression)
    manifest.to_csv(
        MANIFEST_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    adata.file.close()

    print("\nSaved:", EXPRESSION_OUTPUT)
    print("Saved:", MANIFEST_OUTPUT)
    print("\nManifest:")
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
