"""Kidney regulatory support for PAPPA2 and fixed TAL marker genes.

Inputs
------
1. NephQTL2 Supplementary Data 4:
   gene-level eQTL results for glomerular (GLOM) and tubulointerstitial
   (TUBE) kidney compartments.
2. GSE211785 multiome H5AD:
   same-cohort snATAC gene-activity values, used as a chromatin-accessibility
   proxy for C_TAL and M_TAL cells.
3. KPMP fixed marker table:
   used only to select an aTAL2 marker set that is not re-optimized here.

Outputs
-------
outputs/regulatory_support/
  NephQTL2_PAPPA2_gene_level.csv
  NephQTL2_fixed_marker_overlap.csv
  GSE211785_snATAC_donor_pseudobulk.csv
  GSE211785_snATAC_gene_summary.csv
  regulatory_support_summary.csv
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from openpyxl import load_workbook


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
NEPH_XLSX = ROOT / "work" / "nephqtl2_supplementary_data.xlsx"
MARKERS_CSV = ROOT / "work" / "signature_transfer" / "KPMP_TAL_markers.csv"
H5AD = Path(
    r"D:\Codex_temp\PAPPA2_independent_validation"
    r"\GSE211785_PreSCVI.h5ad"
)
OUT_DIR = ROOT / "outputs" / "regulatory_support"

PAPPA2 = "PAPPA2"
CONTROL_GENES = ["CXCL12", "NOTCH3", "HES1", "EPHB2"]
FIXED_ATAC_GENES = [
    PAPPA2,
    "SPP1",
    "KCNIP4",
    "HIF1A",
    "ADAMTS1",
    "NRP1",
    "TFPI",
    "SOD2",
    "NEAT1",
    *CONTROL_GENES,
]
STATES = ["C_TAL", "M_TAL"]
MIN_ATAC_CELLS = 5
N_PERMUTATIONS = 20_000
SEED = 20261028


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(dtype=float)
    result = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return pd.Series(result, index=values.index)
    ranked = np.argsort(p[valid])
    ordered = p[valid][ranked]
    q = ordered * len(ordered) / np.arange(1, len(ordered) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    adjusted = np.empty_like(q)
    adjusted[ranked] = np.clip(q, 0, 1)
    result[valid] = adjusted
    return pd.Series(result, index=values.index)


def permutation_difference(
    values: np.ndarray,
    labels: np.ndarray,
    n_permutations: int = N_PERMUTATIONS,
    seed: int = SEED,
) -> tuple[float, float, int]:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    n_case = int(np.sum(labels == "Disease"))
    n_control = int(np.sum(labels == "Control"))
    if n_case < 2 or n_control < 2:
        return np.nan, np.nan, 0

    observed = (
        values[labels == "Disease"].mean()
        - values[labels == "Control"].mean()
    )
    total = len(values)
    n_choose = min(n_case, n_control)
    combination_count = 1
    for i in range(1, n_choose + 1):
        combination_count = (
            combination_count * (total - i + 1) // i
        )
        if combination_count > n_permutations:
            break

    if combination_count <= n_permutations:
        selections = list(combinations(range(total), n_choose))
    else:
        rng = np.random.default_rng(seed)
        selections = [
            tuple(rng.choice(total, n_choose, replace=False))
            for _ in range(n_permutations)
        ]

    null = np.empty(len(selections), dtype=float)
    for i, selected in enumerate(selections):
        mask = np.zeros(total, dtype=bool)
        mask[list(selected)] = True
        if n_choose != n_case:
            mask = ~mask
        null[i] = values[mask].mean() - values[~mask].mean()

    p_value = (
        np.sum(np.abs(null) >= abs(observed)) + 1
    ) / (len(null) + 1)
    return observed, p_value, len(null)


def leave_one_donor_range(
    donor_table: pd.DataFrame,
    value_column: str,
) -> tuple[float, float]:
    estimates = []
    for donor in donor_table["donor"]:
        kept = donor_table[~donor_table["donor"].eq(donor)]
        if kept["group"].nunique() < 2:
            continue
        estimate = (
            kept.loc[kept["group"].eq("Disease"), value_column].mean()
            - kept.loc[kept["group"].eq("Control"), value_column].mean()
        )
        estimates.append(estimate)
    if not estimates:
        return np.nan, np.nan
    return float(np.min(estimates)), float(np.max(estimates))


def read_nephqtl2_gene_level() -> pd.DataFrame:
    workbook = load_workbook(NEPH_XLSX, read_only=True, data_only=True)
    sheet = workbook["Sup. Data 4"]
    rows = []
    for row in sheet.iter_rows(min_row=7, values_only=True):
        if row[0] is None or len(row) < 8:
            continue
        rows.append(
            {
                "ensembl_id": row[0],
                "gene": row[1],
                "glom_fdr": row[2],
                "glom_egene": row[3],
                "glom_model_size": row[4],
                "tube_fdr": row[5],
                "tube_egene": row[6],
                "tube_model_size": row[7],
            }
        )
    frame = pd.DataFrame(rows)
    for column in [
        "glom_fdr",
        "glom_egene",
        "glom_model_size",
        "tube_fdr",
        "tube_egene",
        "tube_model_size",
    ]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["gene"] = frame["gene"].astype(str).str.strip()
    return frame


def read_nephqtl2_ptwas() -> pd.DataFrame:
    workbook = load_workbook(NEPH_XLSX, read_only=True, data_only=True)
    sheet = workbook["Sup. Data 9"]
    rows = []
    for row in sheet.iter_rows(min_row=7, values_only=True):
        if row[0] is None:
            continue
        rows.append(
            {
                "ensembl_id": row[0],
                "gene": row[1],
                "tissue": row[2],
                "trait": row[3],
                "p_value": row[4],
                "q_value": row[5],
            }
        )
    frame = pd.DataFrame(rows)
    frame["gene"] = frame["gene"].astype(str).str.strip()
    frame["p_value"] = pd.to_numeric(frame["p_value"], errors="coerce")
    frame["q_value"] = pd.to_numeric(frame["q_value"], errors="coerce")
    return frame


def fixed_marker_overlap(
    gene_level: pd.DataFrame,
    ptwas: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    markers = pd.read_csv(MARKERS_CSV)
    atal2_markers = (
        markers.loc[markers["subtype"].eq("aTAL2"), "gene"]
        .astype(str)
        .str.strip()
        .drop_duplicates()
        .tolist()
    )
    fixed_rows = []
    for gene in sorted(set([PAPPA2, *CONTROL_GENES, *atal2_markers])):
        row = {
            "gene": gene,
            "fixed_set": (
                "PAPPA2"
                if gene == PAPPA2
                else "control"
                if gene in CONTROL_GENES
                else "aTAL2_top35"
            ),
        }
        eqtl = gene_level[gene_level["gene"].eq(gene)]
        if eqtl.empty:
            row.update(
                {
                    "tube_fdr": np.nan,
                    "tube_egene": np.nan,
                    "tube_model_size": np.nan,
                    "glom_fdr": np.nan,
                }
            )
        else:
            row.update(
                {
                    "tube_fdr": eqtl["tube_fdr"].iloc[0],
                    "tube_egene": eqtl["tube_egene"].iloc[0],
                    "tube_model_size": eqtl["tube_model_size"].iloc[0],
                    "glom_fdr": eqtl["glom_fdr"].iloc[0],
                }
            )
        hits = ptwas[
            ptwas["gene"].eq(gene)
            & ptwas["tissue"].isin(["TUBE", "GLOM"])
            & ptwas["trait"].isin(["eGFR", "UACR"])
        ]
        for _, hit in hits.iterrows():
            key = f"ptwas_{hit['tissue'].lower()}_{str(hit['trait']).lower()}_q"
            row[key] = hit["q_value"]
        fixed_rows.append(row)

    fixed = pd.DataFrame(fixed_rows)
    group_summary = (
        fixed.groupby("fixed_set", dropna=False)
        .agg(
            n_genes=("gene", "size"),
            n_tube_egene=("tube_egene", lambda x: int(np.nansum(x))),
            median_tube_fdr=("tube_fdr", "median"),
        )
        .reset_index()
    )
    return fixed, group_summary


def load_atac_gene_activity() -> pd.DataFrame:
    adata = ad.read_h5ad(H5AD, backed="r")
    obs = adata.obs
    mask = (
        obs["tech"].astype(str).eq("SN_ATAC")
        & obs["Cluster_Idents"].astype(str).isin(STATES)
        & obs["sample"].notna()
        & obs["group"].astype(str).isin(["Control", "DKD", "HKD"])
    )
    genes = [gene for gene in FIXED_ATAC_GENES if gene in adata.var_names]
    subset = adata[mask, genes].to_memory()
    matrix = subset.X.toarray() if hasattr(subset.X, "toarray") else np.asarray(subset.X)
    frame = pd.DataFrame(matrix, columns=genes, index=subset.obs.index)
    frame["donor"] = subset.obs["sample"].astype(str).to_numpy()
    frame["state"] = subset.obs["Cluster_Idents"].astype(str).to_numpy()
    frame["group_name"] = subset.obs["group"].astype(str).to_numpy()
    frame["group"] = np.where(
        frame["group_name"].eq("Control"),
        "Control",
        "Disease",
    )
    return frame


def atac_donor_analysis(
    activity: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    donor_tables = []
    summary_rows = []
    genes = [
        column
        for column in activity.columns
        if column not in {"donor", "state", "group_name", "group"}
    ]
    for state in STATES:
        state_block = activity[activity["state"].eq(state)]
        for gene in genes:
            donor_table = (
                state_block.groupby(
                    ["donor", "group", "group_name"],
                    observed=True,
                )[gene]
                .agg(["mean", "size"])
                .reset_index()
                .rename(
                    columns={
                        "mean": "mean_gene_activity",
                        "size": "n_cells",
                    }
                )
            )
            donor_table = donor_table[
                donor_table["n_cells"] >= MIN_ATAC_CELLS
            ].copy()
            donor_table["state"] = state
            donor_table["gene"] = gene
            donor_tables.append(donor_table)
            if donor_table["group"].nunique() < 2:
                continue
            difference, p_value, n_perm = permutation_difference(
                donor_table["mean_gene_activity"].to_numpy(),
                donor_table["group"].to_numpy(),
            )
            loo_min, loo_max = leave_one_donor_range(
                donor_table,
                "mean_gene_activity",
            )
            summary_rows.append(
                {
                    "state": state,
                    "gene": gene,
                    "difference_disease_minus_control": difference,
                    "p_exact_or_monte_carlo": p_value,
                    "n_permutations": n_perm,
                    "n_disease_donors": int(
                        donor_table["group"].eq("Disease").sum()
                    ),
                    "n_control_donors": int(
                        donor_table["group"].eq("Control").sum()
                    ),
                    "loo_min": loo_min,
                    "loo_max": loo_max,
                }
            )
    donors = pd.concat(donor_tables, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    summary["fdr"] = bh_adjust(summary["p_exact_or_monte_carlo"])
    return donors, summary


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    gene_level = read_nephqtl2_gene_level()
    ptwas = read_nephqtl2_ptwas()
    fixed_markers, marker_summary = fixed_marker_overlap(
        gene_level,
        ptwas,
    )

    pappa2_rows = gene_level[gene_level["gene"].eq(PAPPA2)].copy()
    pappa2_rows.to_csv(
        OUT_DIR / "NephQTL2_PAPPA2_gene_level.csv",
        index=False,
    )
    fixed_markers.to_csv(
        OUT_DIR / "NephQTL2_fixed_marker_overlap.csv",
        index=False,
    )
    marker_summary.to_csv(
        OUT_DIR / "NephQTL2_fixed_marker_group_summary.csv",
        index=False,
    )

    activity = load_atac_gene_activity()
    donor_table, atac_summary = atac_donor_analysis(activity)
    donor_table.to_csv(
        OUT_DIR / "GSE211785_snATAC_donor_pseudobulk.csv",
        index=False,
    )
    atac_summary.to_csv(
        OUT_DIR / "GSE211785_snATAC_gene_summary.csv",
        index=False,
    )

    pappa2_atac = atac_summary[atac_summary["gene"].eq(PAPPA2)].copy()
    summary_rows = []
    if not pappa2_rows.empty:
        row = pappa2_rows.iloc[0]
        summary_rows.append(
            {
                "module": "NephQTL2_eQTL",
                "analysis_unit": "tubulointerstitial_eGene",
                "result": str(row["tube_egene"]),
                "effect_or_fdr": row["tube_fdr"],
                "interpretation": (
                    "PAPPA2 is a significant TUBE eGene; GLOM is negative"
                ),
            }
        )
    for _, row in pappa2_atac.iterrows():
        summary_rows.append(
            {
                "module": "GSE211785_snATAC",
                "analysis_unit": f"{row['state']}_donor_pseudobulk",
                "result": row["difference_disease_minus_control"],
                "effect_or_fdr": row["fdr"],
                "interpretation": (
                    "Chromatin-accessibility proxy; same-cohort orthogonal "
                    "modality, not an independent cohort"
                ),
            }
        )
    pd.DataFrame(summary_rows).to_csv(
        OUT_DIR / "regulatory_support_summary.csv",
        index=False,
    )

    print("NephQTL2 PAPPA2 rows:")
    print(pappa2_rows.to_string(index=False))
    print("\nNephQTL2 fixed marker groups:")
    print(marker_summary.to_string(index=False))
    print("\nGSE211785 snATAC PAPPA2:")
    print(pappa2_atac.to_string(index=False))
    print("\nFixed marker genes tested in ATAC:")
    print(np.sort(atac_summary["gene"].unique()))


if __name__ == "__main__":
    main()
