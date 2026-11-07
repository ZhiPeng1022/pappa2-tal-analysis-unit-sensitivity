"""Run donor-level exploratory eGFR association in GSE217427."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
AUDIT = ROOT / "work" / "bulk_clinical_cohort_audit"
OUT = ROOT / "outputs" / "gse217427_donor_sensitivity"
MATRIX = AUDIT / "GSE217427_expression_matrix.TPM.gct.gz"
METADATA = AUDIT / "geo_all_kidney_characteristics_long.csv"

SIGNATURE = [
    "PAPPA2",
    "TMPRSS4",
    "LIX1",
    "WNT10A",
    "KRT7",
    "CRIP1",
    "S100A14",
    "LAMC2",
    "LAMB3",
    "SOD3",
    "B3GALT5",
]


def zscore(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        raise ValueError(f"Zero-variance vector: {values.name}")
    return (values - values.mean()) / sd


def fit_ols(
    data: pd.DataFrame,
    formula: str,
    model_name: str,
    role: str,
) -> dict:
    model = smf.ols(formula, data=data).fit()
    robust = model.get_robustcov_results(cov_type="HC3")
    index = list(model.params.index).index("score_sd")
    beta = float(model.params["score_sd"])
    se = float(robust.bse[index])
    return {
        "model": model_name,
        "role": role,
        "n_donors": int(model.nobs),
        "coefficient": beta,
        "robust_se_hc3": se,
        "ci_lower": beta - 1.96 * se,
        "ci_upper": beta + 1.96 * se,
        "p_value": float(robust.pvalues[index]),
        "r_squared": float(model.rsquared),
    }


def read_expression() -> pd.DataFrame:
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as f:
        f.readline()
        f.readline()
        expression = pd.read_csv(f, sep="\t")
    expression["Description"] = expression["Description"].astype(str).str.upper()
    selected = expression[
        expression["Description"].isin(SIGNATURE)
    ].copy()
    selected = selected.drop_duplicates("Description").set_index("Description")
    return selected.drop(columns=["Name"])


def read_metadata() -> pd.DataFrame:
    data = pd.read_csv(METADATA, dtype=str)
    data = data[data["accession"].eq("GSE217427")]
    wide = data.pivot_table(
        index=["gsm", "sample_title"],
        columns="characteristic",
        values="value",
        aggfunc="first",
    ).reset_index()
    wide["egfr"] = pd.to_numeric(wide["egfr"], errors="coerce")
    wide["age"] = pd.to_numeric(wide["age"], errors="coerce")
    wide["bmi"] = pd.to_numeric(wide["bmi"], errors="coerce")
    return wide


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    expression = read_expression()
    metadata = read_metadata()
    expression = expression.rename(columns=lambda column: column.strip())
    metadata["sample_title"] = metadata["sample_title"].str.strip()
    if set(metadata["sample_title"]) != set(expression.columns):
        missing = set(metadata["sample_title"]) - set(expression.columns)
        extra = set(expression.columns) - set(metadata["sample_title"])
        raise ValueError(f"Sample mismatch. Missing={missing}, extra={extra}")

    long = expression.T.reset_index().rename(
        columns={"index": "sample_title"}
    )
    long = long.merge(
        metadata,
        on="sample_title",
        how="inner",
        validate="one_to_one",
    )
    expression_long = long[
        ["gsm", "sample_title", "donor", "tissue", "disease state",
         "age", "Sex", "bmi", "egfr"]
    ].copy()
    expression_long[SIGNATURE] = long[SIGNATURE].apply(
        pd.to_numeric,
        errors="coerce",
    )

    donor_meta = expression_long.groupby("donor").agg(
        n_tissues=("sample_title", "size"),
        age=("age", "first"),
        sex=("Sex", "first"),
        bmi=("bmi", "first"),
        disease_state=("disease state", "first"),
        egfr=("egfr", "first"),
        egfr_min=("egfr", "min"),
        egfr_max=("egfr", "max"),
    )
    if not (donor_meta["egfr_min"] == donor_meta["egfr_max"]).all():
        raise ValueError("eGFR differs between tissues within a donor")

    donor_expression = expression_long.groupby("donor")[SIGNATURE].mean()
    log_expression = np.log2(donor_expression + 1)
    gene_z = log_expression.apply(zscore, axis=0)
    donor_score = gene_z.mean(axis=1)
    donor_data = donor_meta.reset_index()
    donor_data["score_sd"] = zscore(donor_score).values
    donor_data["PAPPA2_sd"] = gene_z["PAPPA2"].values
    donor_data = donor_data.rename(
        columns={"disease_state": "disease"}
    )
    donor_data.to_csv(OUT / "GSE217427_donor_scores.csv", index=False)

    results = [
        fit_ols(
            donor_data,
            "egfr ~ score_sd",
            "donor_average_unadjusted",
            "exploratory primary",
        ),
        fit_ols(
            donor_data,
            "egfr ~ score_sd + age + C(sex) + C(disease)",
            "donor_average_adjusted",
            "exploratory sensitivity",
        ),
    ]
    pd.DataFrame(results).to_csv(
        OUT / "GSE217427_donor_association_models.csv",
        index=False,
    )
    print(pd.DataFrame(results).to_string(index=False))


if __name__ == "__main__":
    main()
