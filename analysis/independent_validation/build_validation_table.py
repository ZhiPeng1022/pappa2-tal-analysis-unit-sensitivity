"""Build Additional file 14 from the independent-validation result tables."""

import os
from pathlib import Path

import pandas as pd


BASE = Path(
    os.environ.get(
        "PAPPA2_VALIDATION_DATA",
        str(Path(__file__).resolve().parents[2] / "data" / "derived" / "independent_validation"),
    )
)
OUT = BASE / "Additional file 14.csv"
GSE211785 = BASE
GSE254185 = BASE


def add_row(rows, **values):
    rows.append(
        {
            "dataset": values.get("dataset", ""),
            "comparison": values.get("comparison", ""),
            "tal_state": values.get("tal_state", ""),
            "analysis_unit": values.get("analysis_unit", ""),
            "effect": values.get("effect", ""),
            "nominal_p": values.get("nominal_p", ""),
            "fdr": values.get("fdr", ""),
            "n_case": values.get("n_case", ""),
            "n_control": values.get("n_control", ""),
            "n_centers": values.get("n_centers", ""),
            "notes": values.get("notes", ""),
        }
    )


def main() -> None:
    rows = []

    pseudobulk = pd.read_csv(GSE211785 / "GSE211785_pseudobulk.csv")
    for _, row in pseudobulk.iterrows():
        add_row(
            rows,
            dataset="GSE211785",
            comparison=row["comparison"],
            tal_state=row["state"],
            analysis_unit="donor pseudobulk",
            effect=round(float(row["observed_difference"]), 6),
            nominal_p=row["p_permutation"],
            fdr=row["fdr"],
            n_case=row["n_case"],
            n_control=row["n_control"],
            notes="two-sided donor-label permutation",
        )

    for filename, unit in [
        ("GSE211785_within_injury.csv", "within-donor rank, injury score"),
        ("GSE211785_within_nhood.csv", "within-donor rank, neighborhood logFC"),
    ]:
        table = pd.read_csv(GSE211785 / filename)
        for _, row in table.iterrows():
            add_row(
                rows,
                dataset="GSE211785",
                comparison="Disease vs Control",
                tal_state=row["state"],
                analysis_unit=unit,
                effect=round(float(row["observed_rho_difference"]), 6),
                nominal_p=row["p_permutation"],
                fdr=row["fdr"],
                n_case=row["n_case"],
                n_control=row["n_control"],
                notes="donor-level Spearman difference; two-sided permutation",
            )

    neighborhood = pd.read_csv(
        GSE211785 / "GSE211785_neighborhood_expression.csv"
    )
    neighborhood = neighborhood[neighborhood["gene"].eq("PAPPA2")]
    for _, row in neighborhood.iterrows():
        add_row(
            rows,
            dataset="GSE211785",
            comparison="Disease vs Control",
            tal_state=row["state"],
            analysis_unit="neighborhood-expression association",
            effect=round(float(row["rho"]), 6),
            nominal_p=row["p_value"],
            fdr=row["fdr"],
            n_centers=row["n_centers"],
            notes="descriptive; neighborhood centers overlap",
        )

    pseudobulk_254 = pd.read_csv(GSE254185 / "GSE254185_pseudobulk_exact.csv")
    pseudobulk_254 = pseudobulk_254[
        pseudobulk_254["state"].isin(["mTAL", "TAL Injured"])
    ]
    for _, row in pseudobulk_254.iterrows():
        add_row(
            rows,
            dataset="GSE254185",
            comparison="obstructive nephropathy vs normal",
            tal_state=row["state"],
            analysis_unit="donor pseudobulk",
            effect=round(float(row["observed_difference"]), 6),
            nominal_p=row["p_exact"],
            fdr=row["fdr_exact"],
            n_case=row["n_uuo_donors"],
            n_control=row["n_normal_donors"],
            notes="two-sided donor-label permutation",
        )

    for filename, unit in [
        ("GSE254185_within_injury_exact.csv", "within-donor rank, injury score"),
        ("GSE254185_within_nhood_exact.csv", "within-donor rank, neighborhood logFC"),
    ]:
        table = pd.read_csv(GSE254185 / filename)
        table = table[table["state"].eq("TAL Injured")]
        for _, row in table.iterrows():
            add_row(
                rows,
                dataset="GSE254185",
                comparison="obstructive nephropathy vs normal",
                tal_state=row["state"],
                analysis_unit=unit,
                effect=round(float(row["observed_rho_difference"]), 6),
                nominal_p=row["p_exact"],
                fdr=row["fdr_exact"],
                n_case=row["n_uuo_donors"],
                n_control=row["n_normal_donors"],
                notes="donor-level Spearman difference; two-sided permutation",
            )

    neighborhood_254 = pd.read_csv(
        GSE254185 / "GSE254185_neighborhood_expression.csv"
    )
    neighborhood_254 = neighborhood_254[
        neighborhood_254["gene"].eq("PAPPA2")
        & neighborhood_254["state"].eq("TAL Injured")
    ]
    for _, row in neighborhood_254.iterrows():
        add_row(
            rows,
            dataset="GSE254185",
            comparison="obstructive nephropathy vs normal",
            tal_state=row["state"],
            analysis_unit="neighborhood-expression association",
            effect=round(float(row["rho"]), 6),
            nominal_p=row["p_value"],
            fdr=row["fdr"],
            n_centers=row["n_centers"],
            notes="descriptive; neighborhood centers overlap",
        )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"Wrote {OUT} with {len(rows)} rows")


if __name__ == "__main__":
    main()
