"""Run the locked 11-gene clinical association in GSE175759."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import pearsonr, spearmanr


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
AUDIT = ROOT / "work" / "bulk_clinical_cohort_audit"
OUT = ROOT / "outputs" / "gse175759_locked_validation"
COUNTS = OUT / "GSE175759_11gene_raw_counts.csv"
QC = OUT / "GSE175759_sample_qc.csv"
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
    term: str,
    model_name: str,
    role: str,
) -> dict:
    model = smf.ols(formula, data=data).fit()
    robust = model.get_robustcov_results(cov_type="HC3")
    index = list(model.params.index).index(term)
    beta = float(model.params[term])
    se = float(robust.bse[index])
    ci_lower = beta - 1.96 * se
    ci_upper = beta + 1.96 * se
    return {
        "model": model_name,
        "role": role,
        "n": int(model.nobs),
        "coefficient": beta,
        "robust_se_hc3": se,
        "ci_lower": ci_lower,
        "ci_upper": ci_upper,
        "p_value": float(robust.pvalues[index]),
        "r_squared": float(model.rsquared),
    }


def load_metadata() -> pd.DataFrame:
    characteristics = pd.read_csv(METADATA, dtype=str)
    characteristics = characteristics[
        characteristics["accession"].eq("GSE175759")
    ]
    wide = characteristics.pivot_table(
        index="gsm",
        columns="characteristic",
        values="value",
        aggfunc="first",
    ).reset_index()
    wide = wide.rename(
        columns={
            "estimated gfr (ckd-epi)": "eGFR",
            "technical outlier": "technical_outlier",
        }
    )
    wide["eGFR"] = pd.to_numeric(wide["eGFR"], errors="coerce")
    return wide


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(COUNTS)
    metadata = load_metadata()
    data = (
        counts.merge(
            metadata,
            on="gsm",
            how="inner",
            validate="one_to_one",
        )
    )
    if len(data) != 90:
        raise ValueError(f"Expected 90 samples, found {len(data)}")

    data = data[
        ~data["technical_outlier"].str.contains(
            "Outlier",
            case=False,
            na=False,
        )
    ].copy()
    data = data[data["eGFR"].notna()].copy()
    if len(data) < 80:
        raise ValueError(f"Too few analysable samples: {len(data)}")

    gene_counts = data[SIGNATURE].apply(pd.to_numeric, errors="coerce")
    library_size = pd.to_numeric(
        data["library_size"],
        errors="coerce",
    )
    log_cpm = np.log2(gene_counts.div(library_size, axis=0) * 1e6 + 1)
    gene_z = log_cpm.apply(zscore, axis=0)
    score = gene_z.mean(axis=1)
    data["score_sd"] = zscore(score)
    data["PAPPA2_sd"] = gene_z["PAPPA2"]
    data["log_library_size"] = np.log(library_size)
    data.to_csv(OUT / "GSE175759_sample_scores.csv", index=False)

    primary = fit_ols(
        data,
        "eGFR ~ score_sd",
        "score_sd",
        "primary_unadjusted",
        "prespecified primary",
    )
    diagnosis_adjusted = fit_ols(
        data,
        "eGFR ~ score_sd + C(diagnosis)",
        "score_sd",
        "diagnosis_adjusted",
        "prespecified secondary",
    )
    full_adjusted = fit_ols(
        data,
        "eGFR ~ score_sd + C(diagnosis) + log_library_size",
        "score_sd",
        "diagnosis_and_library_adjusted",
        "sensitivity",
    )
    pd.DataFrame(
        [primary, diagnosis_adjusted, full_adjusted]
    ).to_csv(
        OUT / "primary_and_sensitivity_models.csv",
        index=False,
    )

    pappa2_model = fit_ols(
        data,
        "eGFR ~ PAPPA2_sd",
        "PAPPA2_sd",
        "PAPPA2_only_unadjusted",
        "prespecified secondary",
    )
    pd.DataFrame([pappa2_model]).to_csv(
        OUT / "pappa2_only_model.csv",
        index=False,
    )

    leave_one_out = []
    for omitted in SIGNATURE:
        genes = [gene for gene in SIGNATURE if gene != omitted]
        frame = data.copy()
        frame["score_sd"] = zscore(gene_z[genes].mean(axis=1))
        result = fit_ols(
            frame,
            "eGFR ~ score_sd",
            "score_sd",
            f"leave_out_{omitted}",
            "sensitivity",
        )
        result["omitted_gene"] = omitted
        leave_one_out.append(result)
    pd.DataFrame(leave_one_out).to_csv(
        OUT / "leave_one_gene_out.csv",
        index=False,
    )

    pearson_r, pearson_p = pearsonr(data["score_sd"], data["eGFR"])
    spearman_r, spearman_p = spearmanr(data["score_sd"], data["eGFR"])
    pd.DataFrame(
        [
            {
                "method": "Pearson",
                "estimate": float(pearson_r),
                "p_value": float(pearson_p),
            },
            {
                "method": "Spearman",
                "estimate": float(spearman_r),
                "p_value": float(spearman_p),
            },
        ]
    ).to_csv(
        OUT / "score_egfr_correlations.csv",
        index=False,
    )

    diagnosis_counts = (
        data.groupby("diagnosis")
        .agg(
            n=("gsm", "size"),
            mean_eGFR=("eGFR", "mean"),
            mean_score_sd=("score_sd", "mean"),
        )
        .reset_index()
    )
    diagnosis_counts.to_csv(
        OUT / "diagnosis_summary.csv",
        index=False,
    )
    summary = {
        "cohort": "GSE175759",
        "n_after_exclusions": len(data),
        "excluded_technical_outliers": 3,
        "score_direction_expected": "negative",
        "genes": SIGNATURE,
        "primary_model": primary,
        "diagnosis_adjusted_model": diagnosis_adjusted,
        "diagnosis_and_library_adjusted_model": full_adjusted,
        "pappa2_only_model": pappa2_model,
        "correlations": {
            "pearson": [float(pearson_r), float(pearson_p)],
            "spearman": [float(spearman_r), float(spearman_p)],
        },
        "limitations": [
            "no age",
            "no sex",
            "mixed primary kidney diseases",
            "observational association only",
        ],
    }
    (OUT / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(pd.DataFrame([primary, diagnosis_adjusted, full_adjusted]).to_string(index=False))
    print("\nPAPPA2 only:")
    print(pd.DataFrame([pappa2_model]).to_string(index=False))
    print("\nCorrelations:")
    print(f"Pearson r={pearson_r:.3f}, p={pearson_p:.4g}")
    print(f"Spearman rho={spearman_r:.3f}, p={spearman_p:.4g}")


if __name__ == "__main__":
    main()
