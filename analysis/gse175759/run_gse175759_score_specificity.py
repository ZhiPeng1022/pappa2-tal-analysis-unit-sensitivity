"""Run the post hoc GSE175759 score-specificity analyses."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import pearsonr


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
AUDIT = ROOT / "work" / "bulk_clinical_cohort_audit"
OUT = ROOT / "outputs" / "gse175759_locked_validation"
MAPPING = AUDIT / "GSE217427_expression_matrix.TPM.gct.gz"
SAMPLE_METADATA = AUDIT / "geo_all_kidney_characteristics_long.csv"
BASE_SCORES = OUT / "GSE175759_sample_scores.csv"

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
TAL_CONTENT = ["UMOD", "SLC12A1", "KCNJ1", "CASR", "CLDN16"]
N_NULL = 10_000
RANDOM_SEED = 20261108


def zscore(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        raise ValueError(f"Zero-variance vector: {values.name}")
    return (values - values.mean()) / sd


def load_symbol_map() -> pd.Series:
    with gzip.open(MAPPING, "rt", encoding="utf-8", errors="replace") as f:
        f.readline()
        f.readline()
        annotation = pd.read_csv(
            f,
            sep="\t",
            usecols=["Name", "Description"],
        )
    annotation = annotation.dropna().copy()
    annotation["Name"] = annotation["Name"].astype(str).str.split(".").str[0]
    annotation["Description"] = annotation["Description"].astype(str).str.upper()
    annotation = annotation[
        annotation["Description"].ne("")
        & annotation["Description"].ne("NAN")
    ]
    return annotation.drop_duplicates("Name").set_index("Name")["Description"]


def load_file_list() -> pd.DataFrame:
    rows = []
    for line in (AUDIT / "GSE175759_filelist.txt").read_text(
        encoding="utf-8"
    ).splitlines():
        if not line.startswith("File\t"):
            continue
        fields = line.split("\t")
        file_name = fields[1]
        gsm = file_name.split("_", 1)[0]
        rows.append({"gsm": gsm, "file_name": file_name})
    return pd.DataFrame(rows)


def build_full_matrix() -> pd.DataFrame:
    symbol_map = load_symbol_map()
    files = load_file_list()
    columns = {}
    for row in files.itertuples(index=False):
        counts = pd.read_csv(
            AUDIT / "GSE175759_samples" / row.file_name,
            sep="\t",
            skiprows=1,
            header=None,
            names=["id", "count"],
        )
        counts["base_id"] = counts["id"].astype(str).str.split(".").str[0]
        counts["symbol"] = counts["base_id"].map(symbol_map)
        counts = counts.dropna(subset=["symbol"])
        counts["count"] = pd.to_numeric(counts["count"], errors="coerce").fillna(0)
        columns[row.gsm] = counts.groupby("symbol")["count"].sum()
    matrix = pd.DataFrame(columns)
    return matrix.fillna(0)


def load_metadata() -> pd.DataFrame:
    data = pd.read_csv(SAMPLE_METADATA, dtype=str)
    data = data[data["accession"].eq("GSE175759")]
    wide = data.pivot_table(
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


def fit_terms(
    data: pd.DataFrame,
    formula: str,
    terms: list[str],
    model_name: str,
) -> pd.DataFrame:
    model = smf.ols(formula, data=data).fit()
    robust = model.get_robustcov_results(cov_type="HC3")
    rows = []
    for term in terms:
        index = list(model.params.index).index(term)
        beta = float(model.params[term])
        se = float(robust.bse[index])
        rows.append(
            {
                "model": model_name,
                "term": term,
                "n": int(model.nobs),
                "coefficient": beta,
                "robust_se_hc3": se,
                "ci_lower": beta - 1.96 * se,
                "ci_upper": beta + 1.96 * se,
                "p_value": float(robust.pvalues[index]),
                "r_squared": float(model.rsquared),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw_counts = build_full_matrix()
    raw_counts.to_csv(OUT / "GSE175759_full_gene_raw_counts.csv.gz")

    metadata = load_metadata()
    keep = metadata[
        ~metadata["technical_outlier"].str.contains(
            "Outlier",
            case=False,
            na=False,
        )
    ][["gsm", "eGFR", "diagnosis"]].copy()
    samples = keep["gsm"].tolist()
    counts = raw_counts[samples]
    library_size = counts.sum(axis=0)
    log_cpm = np.log2(counts.div(library_size, axis=1) * 1e6 + 1)
    log_cpm.to_csv(OUT / "GSE175759_full_gene_logcpm.csv.gz")

    expression = log_cpm.T
    expression.index.name = "gsm"
    data = expression.reset_index().merge(
        keep,
        on="gsm",
        how="inner",
        validate="one_to_one",
    )
    if len(data) != 87:
        raise ValueError(f"Expected 87 samples, found {len(data)}")

    missing_base = [
        gene
        for gene in SIGNATURE + TAL_CONTENT
        if gene not in data.columns
    ]
    if missing_base:
        raise ValueError(f"Missing required genes: {missing_base}")

    base_scores = pd.read_csv(BASE_SCORES)
    data = data.merge(
        base_scores[["gsm", "score_sd"]],
        on="gsm",
        how="inner",
        validate="one_to_one",
    )

    # Analysis 1: TAL-content adjustment
    tal_z = data[TAL_CONTENT].apply(zscore, axis=0)
    data["TAL_content_sd"] = zscore(tal_z.mean(axis=1))
    tal_models = pd.concat(
        [
            fit_terms(
                data,
                "eGFR ~ score_sd",
                ["score_sd"],
                "injury_only",
            ),
            fit_terms(
                data,
                "eGFR ~ score_sd + TAL_content_sd",
                ["score_sd", "TAL_content_sd"],
                "TAL_content_adjusted",
            ),
            fit_terms(
                data,
                "eGFR ~ score_sd + TAL_content_sd + C(diagnosis)",
                ["score_sd", "TAL_content_sd"],
                "TAL_content_and_diagnosis_adjusted",
            ),
        ],
        ignore_index=True,
    )
    tal_models.to_csv(OUT / "tal_content_models.csv", index=False)
    score_tal_r, score_tal_p = pearsonr(
        data["score_sd"],
        data["TAL_content_sd"],
    )

    # Analysis 2: expression-stratified matched 11-gene sets
    detection = log_cpm.gt(0).mean(axis=1)
    mean_expression = log_cpm.mean(axis=1)
    eligible = log_cpm[
        detection.ge(0.8)
        & mean_expression.ge(1.0)
        & ~log_cpm.index.isin(SIGNATURE + TAL_CONTENT)
    ].copy()
    eligible_means = eligible.mean(axis=1)
    bin_edges = np.quantile(
        eligible_means.to_numpy(dtype=float),
        np.linspace(0, 1, 21),
    )
    control_bins = np.searchsorted(
        bin_edges[1:-1],
        eligible_means.to_numpy(dtype=float),
        side="right",
    )
    eligible_by_bin = {
        int(bin_id): eligible_means.iloc[
            np.where(control_bins == bin_id)[0]
        ].index.tolist()
        for bin_id in sorted(set(control_bins))
    }
    signature_bins = {
        gene: int(
            np.searchsorted(
                bin_edges[1:-1],
                float(mean_expression.loc[gene]),
                side="right",
            )
        )
        for gene in SIGNATURE
    }

    rng = np.random.default_rng(RANDOM_SEED)
    patient_egfr = data.set_index("gsm").loc[samples, "eGFR"]
    egfr_values = patient_egfr.to_numpy(dtype=float)
    egfr_centered = egfr_values - egfr_values.mean()
    gene_z = eligible.apply(zscore, axis=1)
    null_betas = np.empty(N_NULL, dtype=float)
    for index in range(N_NULL):
        selected = []
        used = set()
        for signature_gene in SIGNATURE:
            candidates = eligible_by_bin[signature_bins[signature_gene]]
            available = [gene for gene in candidates if gene not in used]
            if not available:
                raise RuntimeError(f"Empty matched bin for {signature_gene}")
            gene = available[int(rng.integers(0, len(available)))]
            selected.append(gene)
            used.add(gene)
        score = gene_z.loc[selected].mean(axis=0)
        score = (score - score.mean()) / score.std(ddof=0)
        null_betas[index] = float(
            np.dot(egfr_centered, score.to_numpy(dtype=float))
            / len(samples)
        )

    observed_model = smf.ols(
        "eGFR ~ score_sd",
        data=data,
    ).fit()
    observed_beta = float(observed_model.params["score_sd"])
    empirical_p = (
        1 + int(np.sum(null_betas <= observed_beta))
    ) / (N_NULL + 1)
    percentile = float(
        100 * np.mean(null_betas <= observed_beta)
    )
    pd.DataFrame(
        {
            "matched_set": np.arange(1, N_NULL + 1),
            "coefficient": null_betas,
        }
    ).to_csv(
        OUT / "matched_set_null_coefficients.csv",
        index=False,
    )
    pd.DataFrame(
        [
            {
                "observed_coefficient": observed_beta,
                "null_mean": float(null_betas.mean()),
                "null_sd": float(null_betas.std(ddof=1)),
                "null_2.5_percentile": float(np.percentile(null_betas, 2.5)),
                "null_97.5_percentile": float(np.percentile(null_betas, 97.5)),
                "observed_percentile": percentile,
                "empirical_negative_tail_p": empirical_p,
                "n_null": N_NULL,
                "eligible_control_genes": len(eligible),
            }
        ]
    ).to_csv(
        OUT / "matched_set_null_summary.csv",
        index=False,
    )

    # Analysis 3: PAPPA2 versus the remaining 10 genes
    data["PAPPA2_sd"] = zscore(data["PAPPA2"])
    remaining = [gene for gene in SIGNATURE if gene != "PAPPA2"]
    remaining_z = data[remaining].apply(zscore, axis=0)
    data["score_without_PAPPA2_sd"] = zscore(remaining_z.mean(axis=1))
    incremental_models = pd.concat(
        [
            fit_terms(
                data,
                "eGFR ~ PAPPA2_sd",
                ["PAPPA2_sd"],
                "PAPPA2_only",
            ),
            fit_terms(
                data,
                "eGFR ~ PAPPA2_sd + score_without_PAPPA2_sd",
                ["PAPPA2_sd", "score_without_PAPPA2_sd"],
                "PAPPA2_plus_remaining_10",
            ),
        ],
        ignore_index=True,
    )
    incremental_models.to_csv(
        OUT / "pappa2_incremental_models.csv",
        index=False,
    )
    reduced = smf.ols("eGFR ~ PAPPA2_sd", data=data).fit()
    full = smf.ols(
        "eGFR ~ PAPPA2_sd + score_without_PAPPA2_sd",
        data=data,
    ).fit()
    f_value, f_p, f_df_diff = full.compare_f_test(reduced)
    pd.DataFrame(
        [
            {
                "n": len(data),
                "reduced_r_squared": float(reduced.rsquared),
                "full_r_squared": float(full.rsquared),
                "delta_r_squared": float(full.rsquared - reduced.rsquared),
                "nested_f": float(f_value),
                "nested_f_p": float(f_p),
                "df_difference": float(f_df_diff),
            }
        ]
    ).to_csv(
        OUT / "pappa2_incremental_comparison.csv",
        index=False,
    )

    summary = {
        "n_samples": len(data),
        "tal_content_genes": TAL_CONTENT,
        "score_tal_correlation": {
            "pearson_r": float(score_tal_r),
            "p_value": float(score_tal_p),
        },
        "matched_set_null": {
            "observed_beta": observed_beta,
            "empirical_negative_tail_p": float(empirical_p),
            "observed_percentile": percentile,
            "random_seed": RANDOM_SEED,
        },
        "pappa2_incremental": {
            "delta_r_squared": float(full.rsquared - reduced.rsquared),
            "nested_f": float(f_value),
            "nested_f_p": float(f_p),
        },
    }
    (OUT / "score_specificity_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    print(tal_models.to_string(index=False))
    print(
        f"\nscore-TAL correlation r={score_tal_r:.3f}, "
        f"P={score_tal_p:.4g}"
    )
    print(
        f"matched-set observed beta={observed_beta:.3f}, "
        f"empirical P={empirical_p:.4g}, percentile={percentile:.2f}"
    )
    print("\nPAPPA2 incremental models:")
    print(incremental_models.to_string(index=False))
    print(
        f"delta R2={full.rsquared - reduced.rsquared:.4f}, "
        f"nested F P={f_p:.4g}"
    )


if __name__ == "__main__":
    main()
