"""Build the Human Genetics Supplementary Tables S1-S19."""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
BMC_DIR = (
    ROOT
    / "outputs"
    / "PAPPA2_BMC_Medical_Genomics_v2_package"
    / "04_Additional_Files"
)
SIGNATURE_DIR = ROOT / "work" / "signature_transfer"
REGULATORY_DIR = ROOT / "outputs" / "regulatory_support"
OUT_DIR = ROOT / "outputs" / "Human_Genetics_Supplementary_Tables"

TITLES = {
    1: "Permutation null distributions for the within-donor cell-rank analysis",
    2: "Complete-data Milo neighborhood summary",
    3: "Subset-level Augur separability summary",
    4: "Complete versus subset pseudobulk direction comparison",
    5: "Neighborhood-expression association robustness",
    6: "Donor-stratified bootstrap replicates for neighborhood association",
    7: "Joint analysis-unit concordance table",
    8: "Within-donor cell-rank association results",
    9: "Complete donor pseudobulk results",
    10: "Leave-one-donor-out replicates for neighborhood association",
    11: "KPMP ordinal eGFR models",
    12: "GSE137570 directional eGFR association",
    13: "Nephroseq precomputed GFR associations",
    14: "Independent non-KPMP validation of PAPPA2 and analysis-unit direction",
    15: "Donor-stratified KPMP TAL marker panel",
    16: "Held-out KPMP signature-transfer validation",
    17: "External aTAL2-like signature transfer",
    18: "NephQTL2 gene-level eQTL results for PAPPA2 and fixed markers",
    19: "GSE211785 snATAC gene-activity analysis",
    20: "Composition and expression decomposition",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def copy_table(source: Path, number: int) -> Path:
    destination = OUT_DIR / f"Supplementary_Table_S{number}.csv"
    shutil.copy2(source, destination)
    return destination


def build_s18() -> Path:
    pappa2 = pd.read_csv(
        REGULATORY_DIR / "NephQTL2_PAPPA2_gene_level.csv"
    )
    markers = pd.read_csv(
        REGULATORY_DIR / "NephQTL2_fixed_marker_overlap.csv"
    )
    pappa2.insert(0, "row_type", "PAPPA2_gene_level")
    markers.insert(0, "row_type", "fixed_marker")
    combined = pd.concat([pappa2, markers], ignore_index=True, sort=False)
    destination = OUT_DIR / "Supplementary_Table_S18.csv"
    combined.to_csv(destination, index=False)
    return destination


def build_s20() -> Path:
    summary = pd.read_csv(
        ROOT
        / "outputs"
        / "composition_expression_decomposition"
        / "composition_expression_summary.csv"
    )
    states = pd.read_csv(
        ROOT
        / "outputs"
        / "composition_expression_decomposition"
        / "composition_expression_state_contributions.csv"
    )
    summary.insert(0, "record_type", "decomposition_summary")
    states.insert(0, "record_type", "state_contribution")
    combined = pd.concat([summary, states], ignore_index=True, sort=False)
    destination = OUT_DIR / "Supplementary_Table_S20.csv"
    combined.to_csv(destination, index=False)
    return destination


def build_s17() -> Path:
    summary = pd.read_csv(
        SIGNATURE_DIR / "signature_transfer_validation_summary.csv"
    )
    summary.insert(0, "record_type", "donor_and_within_donor_summary")

    gse254 = pd.read_csv(
        SIGNATURE_DIR / "GSE254185_signature_predictions.csv",
        usecols=[
            "Annotation.Lvl2",
            "group",
            "predicted_subtype",
            "marker_margin",
        ],
    )
    gse254 = gse254[
        gse254["predicted_subtype"].eq("aTAL2")
        & gse254["marker_margin"].ge(0.2)
    ]
    gse254_composition = (
        gse254.groupby(["Annotation.Lvl2", "group"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .rename(columns={"Annotation.Lvl2": "backbone_state"})
    )
    gse254_composition.insert(0, "dataset", "GSE254185")
    gse254_composition.insert(1, "record_type", "backbone_composition")
    gse254_composition["margin_threshold"] = 0.2

    gse211 = pd.read_csv(
        SIGNATURE_DIR / "GSE211785_signature_predictions.csv",
        usecols=[
            "true_state",
            "group",
            "predicted_subtype",
            "marker_margin",
        ],
    )
    gse211 = gse211[
        gse211["predicted_subtype"].eq("aTAL2")
        & gse211["marker_margin"].ge(0.2)
    ]
    gse211_composition = (
        gse211.groupby(["true_state", "group"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .rename(columns={"true_state": "backbone_state"})
    )
    gse211_composition.insert(0, "dataset", "GSE211785")
    gse211_composition.insert(1, "record_type", "backbone_composition")
    gse211_composition["margin_threshold"] = 0.2

    composition = pd.concat(
        [gse254_composition, gse211_composition],
        ignore_index=True,
        sort=False,
    )
    combined = pd.concat(
        [summary, composition],
        ignore_index=True,
        sort=False,
    )
    destination = OUT_DIR / "Supplementary_Table_S17.csv"
    combined.to_csv(destination, index=False)
    return destination


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    files: dict[int, Path] = {}
    for number in range(1, 15):
        files[number] = copy_table(
            BMC_DIR / f"Additional file {number}.csv",
            number,
        )
    files[15] = copy_table(
        SIGNATURE_DIR / "KPMP_TAL_markers.csv",
        15,
    )
    files[16] = copy_table(
        SIGNATURE_DIR / "KPMP_internal_validation.csv",
        16,
    )
    files[17] = build_s17()
    files[18] = build_s18()
    files[19] = copy_table(
        REGULATORY_DIR / "GSE211785_snATAC_gene_summary.csv",
        19,
    )
    files[20] = build_s20()

    manifest_rows = []
    index_lines = [
        "# Human Genetics Supplementary Tables",
        "",
        "All files are comma-separated values (.csv).",
        "",
    ]
    for number in range(1, 21):
        path = files[number]
        digest = sha256(path)
        manifest_rows.append(
            {
                "table": f"S{number}",
                "title": TITLES[number],
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": digest,
            }
        )
        index_lines.extend(
            [
                f"## Supplementary Table S{number}",
                "",
                TITLES[number],
                "",
                f"File: `{path.name}`",
                "",
                f"SHA256: `{digest}`",
                "",
            ]
        )

    pd.DataFrame(manifest_rows).to_csv(
        OUT_DIR / "Supplementary_manifest.csv",
        index=False,
    )
    (OUT_DIR / "Supplementary_index.md").write_text(
        "\n".join(index_lines),
        encoding="utf-8",
    )
    print(f"Built {len(files)} tables in {OUT_DIR}")


if __name__ == "__main__":
    main()
