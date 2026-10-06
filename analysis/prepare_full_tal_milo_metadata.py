import anndata as ad
import numpy as np
import pandas as pd

from project_paths import DATA_ROOT, OUTPUT_ROOT, PROJECT_ROOT

INPUT_H5AD = DATA_ROOT / "KPMP_new_snRNA.h5ad"
OUTPUT_DIR = DATA_ROOT / "milo_full_tal_input"
OUTPUT_PATH = OUTPUT_DIR / "KPMP_milo_full_tal_cells.csv"
MANIFEST_PATH = (
    OUTPUT_ROOT / "KPMP_full_TAL_milo_input_manifest.csv"
)

TAL_SUBCLASSES = [
    "C-TAL-A",
    "C-TAL-B",
    "aTAL2",
    "frTAL",
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


def main():
    if not INPUT_H5AD.exists():
        raise FileNotFoundError(INPUT_H5AD)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PROJECT_ROOT.joinpath("outputs").mkdir(
        parents=True,
        exist_ok=True,
    )

    adata = ad.read_h5ad(INPUT_H5AD, backed="r")
    subclass = adata.obs["SubclassLevel3"].astype(str)
    target_rows = np.flatnonzero(
        subclass.isin(TAL_SUBCLASSES).to_numpy()
    )
    if len(target_rows) == 0:
        raise ValueError("No target TAL cells found.")

    obs = adata.obs.iloc[target_rows].copy()
    umap = np.asarray(
        adata.obsm["X_umap"][target_rows, :]
    )
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
            "condition": condition.to_numpy(),
            "subclass": (
                obs["SubclassLevel3"].astype(str).to_numpy()
            ),
            "region": obs["region"].astype(str).to_numpy(),
            "UMAP1": umap[:, 0],
            "UMAP2": umap[:, 1],
        }
    )

    manifest = (
        metadata.groupby(
            ["subclass", "condition"],
            observed=True,
        )
        .agg(
            cells=("cell_id", "size"),
            donors=("donor_id", "nunique"),
        )
        .reset_index()
        .sort_values(["subclass", "condition"])
    )

    metadata.to_csv(
        OUTPUT_PATH,
        index=False,
        encoding="utf-8-sig",
    )
    manifest.to_csv(
        MANIFEST_PATH,
        index=False,
        encoding="utf-8-sig",
    )
    adata.file.close()

    print("Saved:", OUTPUT_PATH)
    print("Saved:", MANIFEST_PATH)
    print("Cells:", len(metadata))
    print("Donors:", metadata["donor_id"].nunique())
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
