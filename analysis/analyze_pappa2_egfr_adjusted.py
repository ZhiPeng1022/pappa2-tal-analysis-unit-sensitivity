import os
from pathlib import Path
import re

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.miscmodels.ordinal_model import OrderedModel


PROJECT_ROOT = Path(
    os.environ.get(
        "PAPPA2_PROJECT_ROOT",
        Path(__file__).resolve().parents[2],
    )
)
DATA_ROOT = Path(
    os.environ.get(
        "PAPPA2_DATA_ROOT",
        PROJECT_ROOT / "data" / "RAMP3",
    )
)
OUTPUT_ROOT = Path(
    os.environ.get(
        "PAPPA2_OUTPUT_ROOT",
        PROJECT_ROOT / "outputs",
    )
)
WORK_ROOT = Path(
    os.environ.get(
        "PAPPA2_WORK_ROOT",
        PROJECT_ROOT / "work",
    )
)

STATE_GROUPS = (
    DATA_ROOT / "state_pseudobulk" / "KPMP_PAPPA2_state_groups.csv"
)
DONOR_EGFR = (
    OUTPUT_ROOT / "KPMP_donor_eGFR_midpoint.csv"
)

MODEL_DATA_OUTPUT = (
    WORK_ROOT / "KPMP_PAPPA2_TAL_egfr_model_data.csv"
)
MODEL_OUTPUT = (
    WORK_ROOT / "KPMP_PAPPA2_TAL_egfr_models.csv"
)
LINEAR_OUTPUT = (
    WORK_ROOT / "KPMP_PAPPA2_TAL_egfr_linear_sensitivity.csv"
)

TAL_SUBCLASSES = [
    "C-TAL-A",
    "C-TAL-B",
    "C/M-TAL-A",
    "C/M-TAL-B",
    "M-TAL",
    "aTAL1",
    "aTAL2",
    "dC-TAL",
    "dM-TAL",
    "frTAL",
]

EXPOSURE = "pappa2_tal_fraction_per_10pct"


def age_midpoint(value):
    if value is None or pd.isna(value):
        return np.nan
    match = re.search(r"(\d+)-(\d+)", str(value))
    if not match:
        return np.nan
    return (
        float(match.group(1))
        + float(match.group(2))
    ) / 2.0


def egfr_category(value):
    if value is None or pd.isna(value):
        return np.nan
    text = str(value)
    if text == "Not available":
        return np.nan
    if text.startswith(">"):
        # ">60" is not a precise bin. Assign to the
        # lowest category above 60 for the primary analysis.
        return 2
    match = re.search(r"(\d+)-(\d+)", text)
    if not match:
        return np.nan
    lower = int(match.group(1))
    if lower < 30:
        return 0
    if lower < 60:
        return 1
    if lower < 90:
        return 2
    return 3


def disease_group(value):
    text = str(value)
    if text == "normal":
        return "Reference"
    if text == "acute kidney injury":
        return "AKI"
    if text == "chronic kidney disease":
        return "CKD"
    return text


def build_donor_table():
    groups = pd.read_csv(
        STATE_GROUPS,
        dtype={"donor_id": str},
    )
    if not set(TAL_SUBCLASSES).issubset(
        set(groups["subclass"])
    ):
        missing = sorted(
            set(TAL_SUBCLASSES)
            - set(groups["subclass"])
        )
        raise ValueError(
            f"Missing TAL subclasses: {missing}"
        )

    groups["is_positive"] = (
        groups["target_status"] == "positive"
    )
    donor_rows = []
    for donor_id, subset in groups.groupby(
        "donor_id",
        observed=True,
    ):
        total_cells = int(
            subset["cells"].sum()
        )
        positive_cells = int(
            subset.loc[
                subset["is_positive"], "cells"
            ].sum()
        )
        subclass_counts = (
            subset.groupby(
                "subclass",
                observed=True,
            )["cells"]
            .sum()
        )
        row = {
            "donor_id": donor_id,
            "disease": disease_group(
                subset["disease"].iloc[0]
            ),
            "sex": str(subset["sex"].iloc[0]),
            "age": str(subset["age"].iloc[0]),
            "assay": str(subset["assay"].iloc[0]),
            "total_tal_cells": total_cells,
            "positive_tal_cells": positive_cells,
            "pappa2_tal_positive_fraction": (
                positive_cells / total_cells
                if total_cells > 0
                else np.nan
            ),
        }
        weights = subset["cells"].to_numpy(
            dtype=float
        )
        for column in [
            "mean_nCount_RNA",
            "mean_nFeature_RNA",
            "mean_percent_mt",
        ]:
            values = subset[column].to_numpy(
                dtype=float
            )
            row[f"weighted_{column}"] = np.average(
                values,
                weights=weights,
            )
        for subclass in TAL_SUBCLASSES:
            column_name = (
                subclass.lower()
                .replace("-", "_")
                .replace("/", "_")
            )
            row[f"{column_name}_fraction"] = (
                float(
                    subclass_counts.get(
                        subclass,
                        0,
                    )
                )
                / total_cells
            )
        donor_rows.append(row)

    donor_table = pd.DataFrame(donor_rows)
    donor_table["age_midpoint"] = donor_table[
        "age"
    ].map(age_midpoint)
    donor_table[EXPOSURE] = (
        donor_table["pappa2_tal_positive_fraction"]
        / 0.10
    )

    egfr = pd.read_csv(
        DONOR_EGFR,
        dtype={"donor_id": str},
    )
    egfr["egfr_category"] = egfr[
        "egfr_bin"
    ].map(egfr_category)
    egfr["egfr_gt60_ambiguous"] = (
        egfr["egfr_bin"]
        .astype(str)
        .str.startswith(">")
    )
    donor_table = donor_table.merge(
        egfr,
        on="donor_id",
        how="inner",
        validate="one_to_one",
    )
    donor_table["disease"] = pd.Categorical(
        donor_table["disease"],
        categories=["Reference", "AKI", "CKD"],
    )
    donor_table["sex"] = pd.Categorical(
        donor_table["sex"],
        categories=["female", "male"],
    )
    return donor_table


def formula_for_model(
    model_number,
    scope,
):
    formula = f"egfr_category ~ {EXPOSURE}"
    if model_number >= 2:
        formula += " + age_midpoint + C(sex)"
        if scope in ["all", "disease"]:
            formula += " + C(disease)"
    if model_number >= 3:
        formula += (
            " + c_tal_a_fraction"
            " + c_tal_b_fraction"
            " + atal2_fraction"
            " + frtal_fraction"
        )
    if model_number >= 4:
        formula += (
            " + weighted_mean_nFeature_RNA"
            " + C(assay)"
        )
    return formula


def extract_exposure_result(result):
    parameters = result.params
    if EXPOSURE not in parameters.index:
        raise ValueError(
            f"Exposure not found in model: {list(parameters.index)}"
        )
    coefficient = float(parameters[EXPOSURE])
    confidence = result.conf_int().loc[EXPOSURE]
    p_value = float(result.pvalues[EXPOSURE])
    return {
        "coefficient": coefficient,
        "odds_ratio_per_10pct": float(
            np.exp(coefficient)
        ),
        "confidence_lower": float(
            np.exp(confidence.iloc[0])
        ),
        "confidence_upper": float(
            np.exp(confidence.iloc[1])
        ),
        "p_value": p_value,
    }


def fit_ordinal_models(model_data):
    rows = []
    for scope, subset in [
        ("all", model_data),
        (
            "disease",
            model_data[
                model_data["disease"].isin(
                    ["AKI", "CKD"]
                )
            ],
        ),
        (
            "Reference",
            model_data[
                model_data["disease"] == "Reference"
            ],
        ),
        (
            "AKI",
            model_data[
                model_data["disease"] == "AKI"
            ],
        ),
        (
            "CKD",
            model_data[
                model_data["disease"] == "CKD"
            ],
        ),
    ]:
        subset = subset.dropna(
            subset=[
                "egfr_category",
                EXPOSURE,
                "age_midpoint",
                "sex",
                "disease",
                "c_tal_a_fraction",
                "c_tal_b_fraction",
                "atal2_fraction",
                "frtal_fraction",
                "weighted_mean_nFeature_RNA",
                "assay",
            ]
        ).copy()
        subset["disease"] = subset[
            "disease"
        ].astype(str)
        for model_number in [1, 2, 3, 4]:
            formula = formula_for_model(
                model_number,
                scope,
            )
            try:
                model = OrderedModel.from_formula(
                    formula,
                    data=subset,
                    distr="logit",
                )
                result = model.fit(
                    method="bfgs",
                    maxiter=200,
                    disp=False,
                )
                extracted = extract_exposure_result(
                    result
                )
                rows.append(
                    {
                        "scope": scope,
                        "model": f"model_{model_number}",
                        "formula": formula,
                        "donors": len(subset),
                        "converged": bool(
                            result.mle_retvals[
                                "converged"
                            ]
                        ),
                        **extracted,
                        "direction": (
                            "higher_PAPPA2_lower_eGFR"
                            if extracted["coefficient"] < 0
                            else "higher_PAPPA2_higher_eGFR"
                        ),
                    }
                )
            except Exception as error:
                rows.append(
                    {
                        "scope": scope,
                        "model": f"model_{model_number}",
                        "formula": formula,
                        "donors": len(subset),
                        "converged": False,
                        "coefficient": np.nan,
                        "odds_ratio_per_10pct": np.nan,
                        "confidence_lower": np.nan,
                        "confidence_upper": np.nan,
                        "p_value": np.nan,
                        "direction": "not_estimated",
                        "error": str(error),
                    }
                )
    return pd.DataFrame(rows)


def linear_formula(model_number, scope):
    return formula_for_model(
        model_number,
        scope,
    ).replace(
        "egfr_category ~",
        "egfr_midpoint ~",
    )


def fit_linear_sensitivity(model_data):
    rows = []
    for scope, subset in [
        ("all", model_data),
        (
            "disease",
            model_data[
                model_data["disease"].isin(
                    ["AKI", "CKD"]
                )
            ],
        ),
        (
            "Reference",
            model_data[
                model_data["disease"] == "Reference"
            ],
        ),
        (
            "AKI",
            model_data[
                model_data["disease"] == "AKI"
            ],
        ),
        (
            "CKD",
            model_data[
                model_data["disease"] == "CKD"
            ],
        ),
    ]:
        subset = subset.dropna(
            subset=[
                "egfr_midpoint",
                EXPOSURE,
                "age_midpoint",
                "sex",
                "disease",
                "c_tal_a_fraction",
                "c_tal_b_fraction",
                "atal2_fraction",
                "frtal_fraction",
                "weighted_mean_nFeature_RNA",
                "assay",
            ]
        ).copy()
        subset["disease"] = subset[
            "disease"
        ].astype(str)
        for model_number in [1, 2, 3, 4]:
            formula = linear_formula(
                model_number,
                scope,
            )
            try:
                result = smf.ols(
                    formula,
                    data=subset,
                ).fit(cov_type="HC3")
                coefficient = float(
                    result.params[EXPOSURE]
                )
                confidence = result.conf_int().loc[
                    EXPOSURE
                ]
                rows.append(
                    {
                        "scope": scope,
                        "model": f"model_{model_number}",
                        "formula": formula,
                        "donors": len(subset),
                        "coefficient": coefficient,
                        "confidence_lower": float(
                            confidence.iloc[0]
                        ),
                        "confidence_upper": float(
                            confidence.iloc[1]
                        ),
                        "p_value": float(
                            result.pvalues[EXPOSURE]
                        ),
                        "direction": (
                            "higher_PAPPA2_lower_eGFR"
                            if coefficient < 0
                            else "higher_PAPPA2_higher_eGFR"
                        ),
                    }
                )
            except Exception as error:
                rows.append(
                    {
                        "scope": scope,
                        "model": f"model_{model_number}",
                        "formula": formula,
                        "donors": len(subset),
                        "coefficient": np.nan,
                        "confidence_lower": np.nan,
                        "confidence_upper": np.nan,
                        "p_value": np.nan,
                        "direction": "not_estimated",
                        "error": str(error),
                    }
                )
    return pd.DataFrame(rows)


def main():
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    model_data = build_donor_table()
    ordinal = fit_ordinal_models(model_data)
    linear = fit_linear_sensitivity(model_data)
    model_data.to_csv(
        MODEL_DATA_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    ordinal.to_csv(
        MODEL_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    linear.to_csv(
        LINEAR_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    print("Saved:", MODEL_DATA_OUTPUT)
    print("Saved:", MODEL_OUTPUT)
    print("Saved:", LINEAR_OUTPUT)
    print("\nOrdinal PAPPA2 results:")
    print(
        ordinal[
            [
                "scope",
                "model",
                "donors",
                "odds_ratio_per_10pct",
                "confidence_lower",
                "confidence_upper",
                "p_value",
                "direction",
            ]
        ]
        .round(4)
        .to_string(index=False)
    )
    print("\nLinear midpoint sensitivity:")
    print(
        linear[
            [
                "scope",
                "model",
                "donors",
                "coefficient",
                "confidence_lower",
                "confidence_upper",
                "p_value",
                "direction",
            ]
        ]
        .round(4)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
