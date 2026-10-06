import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu

from project_paths import OUTPUT_ROOT

EXPRESSION_INPUT = (
    OUTPUT_ROOT / "KPMP_full_TAL_focus_expression.csv.gz"
)
SUBSET_INPUT = (
    OUTPUT_ROOT / "KPMP_TAL_focus_gene_pseudobulk.csv"
)

FOCUS_GENES = [
    "PAPPA2",
    "CXCL12",
    "NOTCH3",
    "HES1",
    "EPHB2",
]
SUBCLASS_ORDER = [
    "C-TAL-A",
    "C-TAL-B",
    "aTAL2",
    "frTAL",
]
COMPARISONS = {
    "all_disease": ["AKI", "CKD"],
    "AKI_vs_Reference": ["AKI"],
    "CKD_vs_Reference": ["CKD"],
}

FULL_OUTPUT = (
    OUTPUT_ROOT / "KPMP_full_TAL_focus_gene_pseudobulk.csv"
)
COMPARISON_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_vs_subset_TAL_focus_pseudobulk_comparison.csv"
)


def benjamini_hochberg(p_values):
    values = np.asarray(p_values, dtype=float)
    result = np.full_like(values, np.nan)
    valid = np.isfinite(values)
    valid_values = values[valid]
    if len(valid_values) == 0:
        return result
    order = np.argsort(valid_values)
    ranked = valid_values[order]
    n_values = len(valid_values)
    adjusted = ranked * n_values / np.arange(
        1,
        n_values + 1,
    )
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    valid_result = np.empty_like(adjusted)
    valid_result[order] = adjusted
    result[valid] = valid_result
    return result


def age_to_midpoint(value):
    text = str(value)
    if "-" not in text:
        return np.nan
    lower, upper = text.split("-", maxsplit=1)
    try:
        return (float(lower) + float(upper)) / 2.0
    except ValueError:
        return np.nan


def build_donor_pseudobulk(expression):
    expression = expression[
        expression["condition"].isin(
            ["Reference", "AKI", "CKD"]
        )
    ].copy()
    expression_cols = [
        f"{gene}_lognorm" for gene in FOCUS_GENES
    ]
    donor_data = (
        expression.groupby(
            ["donor_id", "condition", "subclass"],
            observed=True,
        )
        .agg(
            cells=("cell_id", "size"),
            age=("age", "first"),
            sex=("sex", "first"),
            **{
                column: (column, "mean")
                for column in expression_cols
            },
        )
        .reset_index()
    )
    donor_data["age_midpoint"] = donor_data["age"].map(
        age_to_midpoint
    )
    donor_data["sex_binary"] = (
        donor_data["sex"]
        .astype(str)
        .map({"male": 1.0, "female": 0.0})
    )
    donor_data = donor_data.rename(
        columns={
            column: column.removesuffix("_lognorm")
            for column in expression_cols
        }
    )
    return donor_data


def adjusted_condition_model(
    reference_values,
    case_values,
):
    values = np.concatenate(
        [reference_values, case_values]
    )
    # The adjusted result is filled in later by caller.
    return values


def compare_expression(pseudobulk):
    rows = []
    for subclass in SUBCLASS_ORDER:
        subclass_data = pseudobulk[
            pseudobulk["subclass"] == subclass
        ]
        reference = subclass_data[
            subclass_data["condition"] == "Reference"
        ]
        for comparison, case_conditions in COMPARISONS.items():
            disease = subclass_data[
                subclass_data["condition"].isin(
                    case_conditions
                )
            ]
            for gene in FOCUS_GENES:
                reference_values = reference[gene].to_numpy(
                    dtype=float
                )
                disease_values = disease[gene].to_numpy(
                    dtype=float
                )
                if (
                    len(reference_values) == 0
                    or len(disease_values) == 0
                ):
                    p_value = np.nan
                elif np.ptp(
                    np.concatenate(
                        [
                            reference_values,
                            disease_values,
                        ]
                    )
                ) == 0:
                    p_value = 1.0
                else:
                    p_value = mannwhitneyu(
                        disease_values,
                        reference_values,
                        alternative="two-sided",
                    ).pvalue

                rows.append(
                    {
                        "comparison": comparison,
                        "subclass": subclass,
                        "gene": gene,
                        "reference_donors": len(
                            reference_values
                        ),
                        "case_donors": len(
                            disease_values
                        ),
                        "reference_positive_fraction": (
                            np.mean(reference_values > 0)
                            if len(reference_values)
                            else np.nan
                        ),
                        "case_positive_fraction": (
                            np.mean(disease_values > 0)
                            if len(disease_values)
                            else np.nan
                        ),
                        "reference_median": np.median(
                            reference_values
                        ),
                        "case_median": np.median(
                            disease_values
                        ),
                        "reference_mean": np.mean(
                            reference_values
                        ),
                        "case_mean": np.mean(
                            disease_values
                        ),
                        "mean_case_minus_reference": (
                            np.mean(disease_values)
                            - np.mean(reference_values)
                        ),
                        "p_value": p_value,
                    }
                )
    result = pd.DataFrame(rows)
    result["FDR"] = benjamini_hochberg(
        result["p_value"]
    )
    result["significant_FDR_0_05"] = result["FDR"] < 0.05
    result["direction"] = np.sign(
        result["mean_case_minus_reference"]
    )
    return result


def add_adjusted_models(result, pseudobulk):
    adjusted_rows = []
    for _, row in result.iterrows():
        subclass = row["subclass"]
        comparison = row["comparison"]
        gene = row["gene"]
        case_conditions = COMPARISONS[comparison]
        data = pseudobulk[
            (pseudobulk["subclass"] == subclass)
            & (
                pseudobulk["condition"].isin(
                    ["Reference", *case_conditions]
                )
            )
        ][
            [
                "condition",
                "age_midpoint",
                "sex_binary",
                gene,
            ]
        ].dropna()
        data = data.rename(columns={gene: "expression"})
        data["case"] = (
            data["condition"] != "Reference"
        ).astype(float)

        beta = np.nan
        p_value = np.nan
        model_status = "not_fit"
        if (
            data["case"].nunique() == 2
            and len(data) >= 8
            and data["age_midpoint"].notna().sum() >= 8
            and data["sex_binary"].notna().sum() >= 8
        ):
            try:
                model = smf.ols(
                    "expression ~ case + age_midpoint + sex_binary",
                    data=data,
                ).fit()
                beta = model.params["case"]
                p_value = model.pvalues["case"]
                model_status = "fit"
            except Exception:
                model_status = "failed"

        adjusted_rows.append(
            {
                "comparison": comparison,
                "subclass": subclass,
                "gene": gene,
                "adjusted_n_donors": len(data),
                "adjusted_condition_beta": beta,
                "adjusted_p_value": p_value,
                "adjusted_model_status": model_status,
            }
        )

    adjusted = pd.DataFrame(adjusted_rows)
    adjusted["adjusted_FDR"] = benjamini_hochberg(
        adjusted["adjusted_p_value"]
    )
    adjusted["adjusted_significant_FDR_0_05"] = (
        adjusted["adjusted_FDR"] < 0.05
    )
    return result.merge(
        adjusted,
        on=["comparison", "subclass", "gene"],
        how="left",
        validate="one_to_one",
    )


def compare_with_subset(full_result):
    subset = pd.read_csv(SUBSET_INPUT)
    subset = subset[
        subset["gene"].isin(FOCUS_GENES)
    ].copy()
    subset["comparison"] = subset["comparison"].replace(
        {"disease_vs_Reference": "all_disease"}
    )
    subset = subset.rename(
        columns={
            "mean_case_minus_reference": (
                "subset_mean_case_minus_reference"
            ),
            "p_value": "subset_p_value",
            "FDR": "subset_FDR",
            "significant_FDR_0_05": (
                "subset_significant_FDR_0_05"
            ),
        }
    )
    subset["subset_direction"] = np.sign(
        subset["subset_mean_case_minus_reference"]
    )
    merged = full_result.merge(
        subset[
            [
                "comparison",
                "subclass",
                "gene",
                "subset_mean_case_minus_reference",
                "subset_p_value",
                "subset_FDR",
                "subset_significant_FDR_0_05",
                "subset_direction",
            ]
        ],
        on=["comparison", "subclass", "gene"],
        how="inner",
        validate="one_to_one",
    )
    merged["direction_agreement"] = (
        merged["direction"] == merged["subset_direction"]
    )
    merged["full_only_significant"] = (
        merged["significant_FDR_0_05"]
        & ~merged["subset_significant_FDR_0_05"]
    )
    merged["subset_only_significant"] = (
        ~merged["significant_FDR_0_05"]
        & merged["subset_significant_FDR_0_05"]
    )
    return merged


def main():
    if not EXPRESSION_INPUT.exists():
        raise FileNotFoundError(EXPRESSION_INPUT)
    expression = pd.read_csv(EXPRESSION_INPUT)
    pseudobulk = build_donor_pseudobulk(expression)
    full_result = compare_expression(pseudobulk)
    full_result = add_adjusted_models(
        full_result,
        pseudobulk,
    )
    comparison = compare_with_subset(full_result)

    full_result.to_csv(
        FULL_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    comparison.to_csv(
        COMPARISON_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )

    print("Saved:", FULL_OUTPUT)
    print("Saved:", COMPARISON_OUTPUT)
    print(
        "Full significant:",
        int(full_result["significant_FDR_0_05"].sum()),
        "/",
        len(full_result),
    )
    print(
        "Adjusted full significant:",
        int(
            full_result[
                "adjusted_significant_FDR_0_05"
            ].sum()
        ),
        "/",
        len(full_result),
    )
    print(
        "Direction agreement with subset:",
        int(comparison["direction_agreement"].sum()),
        "/",
        len(comparison),
    )
    print("\nPAPPA2 full-data results:")
    pappa2 = full_result[
        full_result["gene"] == "PAPPA2"
    ].sort_values(["comparison", "subclass"])
    print(
        pappa2[
            [
                "comparison",
                "subclass",
                "reference_donors",
                "case_donors",
                "mean_case_minus_reference",
                "p_value",
                "FDR",
                "adjusted_condition_beta",
                "adjusted_p_value",
                "adjusted_FDR",
            ]
        ]
        .round(4)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
