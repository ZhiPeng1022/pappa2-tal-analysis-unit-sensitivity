"""Post hoc disease-only sensitivity analyses for GSE175759."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from run_gse175759_locked_clinical_validation import fit_ols


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
OUT = ROOT / "outputs" / "gse175759_locked_validation"
DATA = OUT / "GSE175759_sample_scores.csv"


def main() -> None:
    data = pd.read_csv(DATA)
    disease = data[data["diagnosis"].ne("Control")].copy()
    igan = data[data["diagnosis"].eq("IgAN")].copy()
    rows = [
        fit_ols(
            disease,
            "eGFR ~ score_sd",
            "score_sd",
            "disease_only_unadjusted",
            "post hoc disease-only sensitivity",
        ),
        fit_ols(
            disease,
            "eGFR ~ score_sd + C(diagnosis)",
            "score_sd",
            "disease_only_diagnosis_adjusted",
            "post hoc disease-only sensitivity",
        ),
        fit_ols(
            disease,
            "eGFR ~ score_sd + C(diagnosis) + log_library_size",
            "score_sd",
            "disease_only_diagnosis_and_library_adjusted",
            "post hoc disease-only sensitivity",
        ),
        fit_ols(
            igan,
            "eGFR ~ score_sd",
            "score_sd",
            "IgAN_only_unadjusted",
            "post hoc single-disease sensitivity",
        ),
    ]
    result = pd.DataFrame(rows)
    result.to_csv(
        OUT / "posthoc_disease_subgroup_models.csv",
        index=False,
    )
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
