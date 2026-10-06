from pathlib import Path
import gzip
import re

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr


WORK_ROOT = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2\work"
)
SOFT_PATH = WORK_ROOT / "GSE137570_family.soft.gz"
COUNTS_PATH = (
    WORK_ROOT / "GSE137570_Normalized_read_counts.xlsx"
)
OUTPUT_PATH = (
    WORK_ROOT / "GSE137570_PAPPA2_eGFR_validation.csv"
)
DATA_OUTPUT = (
    WORK_ROOT / "GSE137570_PAPPA2_eGFR_validation_data.csv"
)

GENES = [
    "PAPPA2",
    "CXCL12",
    "NOTCH3",
    "HES1",
    "EPHB2",
]


def parse_soft_metadata():
    text = gzip.open(
        SOFT_PATH,
        "rt",
        encoding="utf-8",
        errors="replace",
    ).read()
    blocks = text.split("^SAMPLE = ")[1:]
    rows = []
    for block in blocks:
        lines = block.splitlines()
        row = {
            "gsm": lines[0].strip(),
            "title": next(
                (
                    line.split("=", 1)[1].strip()
                    for line in lines
                    if line.startswith("!Sample_title")
                ),
                "",
            ),
            "age": next(
                (
                    line.split(":", 1)[1].strip()
                    for line in lines
                    if line.startswith(
                        "!Sample_characteristics_ch1 = age:"
                    )
                ),
                "",
            ),
            "sex": next(
                (
                    line.split(":", 1)[1].strip()
                    for line in lines
                    if line.lower().startswith(
                        "!sample_characteristics_ch1 = sex:"
                    )
                ),
                "",
            ),
            "egfr": next(
                (
                    line.split(":", 1)[1].strip()
                    for line in lines
                    if line.lower().startswith(
                        "!sample_characteristics_ch1 = egfr:"
                    )
                ),
                "",
            ),
            "fibrosis": next(
                (
                    line.split(":", 1)[1].strip()
                    for line in lines
                    if line.lower().startswith(
                        "!sample_characteristics_ch1 = "
                        "% tubulointerstitial fibrosis"
                    )
                ),
                "",
            ),
        }
        rows.append(row)
    metadata = pd.DataFrame(rows)
    metadata["cohort"] = np.where(
        metadata["title"].str.startswith("S"),
        "Cohort 1",
        "Cohort 2",
    )
    return metadata


def load_expression():
    counts = pd.read_excel(
        COUNTS_PATH,
        sheet_name="COHORT 1_NORMALIZED READ COUNTS",
    )
    counts = counts.rename(
        columns={counts.columns[0]: "gene"}
    )
    counts["gene"] = counts["gene"].astype(str)
    selected = counts[
        counts["gene"].isin(GENES)
    ].copy()
    missing = sorted(
        set(GENES) - set(selected["gene"])
    )
    if missing:
        raise ValueError(
            f"Missing genes in GSE137570: {missing}"
        )
    return selected.set_index("gene").T.reset_index().rename(
        columns={"index": "title"}
    )


def main():
    metadata = parse_soft_metadata()
    metadata = metadata[
        metadata["cohort"] == "Cohort 1"
    ].copy()
    metadata["age"] = pd.to_numeric(
        metadata["age"],
        errors="coerce",
    )
    metadata["egfr"] = pd.to_numeric(
        metadata["egfr"],
        errors="coerce",
    )
    metadata["fibrosis"] = pd.to_numeric(
        metadata["fibrosis"],
        errors="coerce",
    )
    expression = load_expression()
    data = metadata.merge(
        expression,
        on="title",
        how="inner",
        validate="one_to_one",
    )
    data = data.dropna(
        subset=[
            "age",
            "egfr",
            "fibrosis",
            "sex",
            *GENES,
        ]
    ).copy()

    rows = []
    for gene in GENES:
        rho, rho_p = spearmanr(
            data[gene],
            data["egfr"],
        )
        formula = (
            f"egfr ~ {gene}"
            " + age + C(sex) + fibrosis"
        )
        result = smf.ols(
            formula,
            data=data,
        ).fit(cov_type="HC3")
        coefficient = float(
            result.params[gene]
        )
        confidence = result.conf_int().loc[gene]
        rows.append(
            {
                "gene": gene,
                "donors": len(data),
                "spearman_rho": float(rho),
                "spearman_p": float(rho_p),
                "adjusted_coefficient": coefficient,
                "adjusted_ci_lower": float(
                    confidence.iloc[0]
                ),
                "adjusted_ci_upper": float(
                    confidence.iloc[1]
                ),
                "adjusted_p": float(
                    result.pvalues[gene]
                ),
                "direction": (
                    "higher_gene_lower_eGFR"
                    if coefficient < 0
                    else "higher_gene_higher_eGFR"
                ),
            }
        )

    results = pd.DataFrame(rows)
    results["FDR"] = np.nan
    order = np.argsort(
        results["adjusted_p"].to_numpy()
    )
    ranked = results["adjusted_p"].to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(
        1,
        len(ranked) + 1,
    )
    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    results.loc[
        results.index[order],
        "FDR",
    ] = adjusted

    data.to_csv(
        DATA_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    results.to_csv(
        OUTPUT_PATH,
        index=False,
        encoding="utf-8-sig",
    )
    print("Saved:", DATA_OUTPUT)
    print("Saved:", OUTPUT_PATH)
    print("\nGSE137570 Cohort 1 validation:")
    print(results.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
