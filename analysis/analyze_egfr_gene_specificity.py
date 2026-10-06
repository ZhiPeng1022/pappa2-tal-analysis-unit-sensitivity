from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.miscmodels.ordinal_model import OrderedModel

from analyze_pappa2_egfr_adjusted import (
    DONOR_EGFR,
    WORK_ROOT,
    egfr_category,
    disease_group,
)


STATE_DIR = Path(
    r"D:\CodexData\RAMP3\state_pseudobulk"
)
GROUPS_PATH = (
    STATE_DIR / "KPMP_PAPPA2_state_groups.csv"
)
GENES_PATH = (
    STATE_DIR / "KPMP_PAPPA2_state_genes.csv"
)
RAW_PATH = (
    STATE_DIR / "KPMP_PAPPA2_state_raw.npz"
)

OUTPUT_PATH = (
    WORK_ROOT
    / "KPMP_PAPPA2_TAL_egfr_gene_specificity.csv"
)
MODEL_DATA_OUTPUT = (
    WORK_ROOT
    / "KPMP_PAPPA2_TAL_egfr_gene_specificity_data.csv"
)

GENES = [
    "PAPPA2",
    "CXCL12",
    "NOTCH3",
    "HES1",
    "EPHB2",
]


def age_midpoint(value):
    text = str(value)
    first, second = text.split("-")
    return (float(first) + float(second)) / 2.0


def benjamini_hochberg(p_values):
    values = np.asarray(p_values, dtype=float)
    result = np.full_like(values, np.nan)
    valid = np.isfinite(values)
    valid_values = values[valid]
    if len(valid_values) == 0:
        return result
    order = np.argsort(valid_values)
    ranked = valid_values[order]
    n_values = len(ranked)
    adjusted = ranked * n_values / np.arange(
        1,
        n_values + 1,
    )
    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    valid_result = np.empty_like(adjusted)
    valid_result[order] = adjusted
    result[valid] = valid_result
    return result


def build_gene_data():
    groups = pd.read_csv(
        GROUPS_PATH,
        dtype={"donor_id": str},
    )
    genes = pd.read_csv(GENES_PATH)
    gene_lookup = dict(
        zip(
            genes["gene"],
            genes["gene_index"].astype(int),
        )
    )
    missing = [
        gene for gene in GENES
        if gene not in gene_lookup
    ]
    if missing:
        raise ValueError(
            f"Missing genes in state matrix: {missing}"
        )

    with np.load(
        RAW_PATH,
        mmap_mode="r",
    ) as state_matrix:
        positive_count = state_matrix[
            "positive_count"
        ]
        donor_tables = []
        for gene in GENES:
            gene_index = gene_lookup[gene]
            gene_positive = np.asarray(
                positive_count[:, gene_index],
                dtype=float,
            )
            table = groups[
                ["donor_id", "cells", "disease"]
            ].copy()
            table["positive_cells"] = gene_positive
            donor = (
                table.groupby(
                    "donor_id",
                    observed=True,
                )
                .agg(
                    cells=("cells", "sum"),
                    positive_cells=(
                        "positive_cells",
                        "sum",
                    ),
                    disease=(
                        "disease",
                        "first",
                    ),
                )
                .reset_index()
            )
            donor["positive_fraction"] = (
                donor["positive_cells"]
                / donor["cells"]
            )
            donor["gene"] = gene
            donor_tables.append(donor)

    fractions = pd.concat(
        donor_tables,
        ignore_index=True,
    )
    fraction_pivot = fractions.pivot(
        index="donor_id",
        columns="gene",
        values="positive_fraction",
    ).reset_index()
    fraction_pivot.columns.name = None

    covariate_rows = []
    for donor_id, subset in groups.groupby(
        "donor_id",
        observed=True,
    ):
        subclass_counts = (
            subset.groupby(
                "subclass",
                observed=True,
            )["cells"]
            .sum()
        )
        total = int(subset["cells"].sum())
        covariate_rows.append(
            {
                "donor_id": donor_id,
                "disease": disease_group(
                    subset["disease"].iloc[0]
                ),
                "sex": str(subset["sex"].iloc[0]),
                "age": str(subset["age"].iloc[0]),
                "total_tal_cells": total,
                "c_tal_a_fraction": (
                    float(
                        subclass_counts.get(
                            "C-TAL-A",
                            0,
                        )
                    )
                    / total
                ),
                "c_tal_b_fraction": (
                    float(
                        subclass_counts.get(
                            "C-TAL-B",
                            0,
                        )
                    )
                    / total
                ),
                "atal2_fraction": (
                    float(
                        subclass_counts.get(
                            "aTAL2",
                            0,
                        )
                    )
                    / total
                ),
                "frtal_fraction": (
                    float(
                        subclass_counts.get(
                            "frTAL",
                            0,
                        )
                    )
                    / total
                ),
            }
        )
    covariates = pd.DataFrame(covariate_rows)
    covariates["age_midpoint"] = covariates[
        "age"
    ].map(age_midpoint)

    egfr = pd.read_csv(
        DONOR_EGFR,
        dtype={"donor_id": str},
    )
    egfr["egfr_category"] = egfr[
        "egfr_bin"
    ].map(egfr_category)

    data = (
        fraction_pivot.merge(
            covariates,
            on="donor_id",
            how="inner",
            validate="one_to_one",
        )
        .merge(
            egfr,
            on="donor_id",
            how="inner",
            validate="one_to_one",
        )
    )
    for gene in GENES:
        data[f"{gene}_per_10pct"] = (
            data[gene] / 0.10
        )
    return data


def formula_for_gene(gene, scope, model_number):
    exposure = f"{gene}_per_10pct"
    formula = f"egfr_category ~ {exposure}"
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
    return formula


def fit_models(data):
    rows = []
    for gene in GENES:
        for scope in ["all", "disease", "CKD", "AKI"]:
            if scope == "all":
                subset = data
            elif scope == "disease":
                subset = data[
                    data["disease"].isin(
                        ["AKI", "CKD"]
                    )
                ]
            else:
                subset = data[
                    data["disease"] == scope
                ]
            subset = subset.copy()
            subset["disease"] = subset[
                "disease"
            ].astype(str)
            subset = subset.dropna(
                subset=[
                    "egfr_category",
                    f"{gene}_per_10pct",
                    "age_midpoint",
                    "sex",
                    "disease",
                    "c_tal_a_fraction",
                    "c_tal_b_fraction",
                    "atal2_fraction",
                    "frtal_fraction",
                ]
            )
            for model_number in [2, 3]:
                formula = formula_for_gene(
                    gene,
                    scope,
                    model_number,
                )
                exposure = f"{gene}_per_10pct"
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
                    coefficient = float(
                        result.params[exposure]
                    )
                    confidence = (
                        result.conf_int().loc[
                            exposure
                        ]
                    )
                    rows.append(
                        {
                            "gene": gene,
                            "scope": scope,
                            "model": (
                                f"model_{model_number}"
                            ),
                            "donors": len(subset),
                            "coefficient": coefficient,
                            "odds_ratio_per_10pct": (
                                float(
                                    np.exp(coefficient)
                                )
                            ),
                            "confidence_lower": float(
                                np.exp(
                                    confidence.iloc[0]
                                )
                            ),
                            "confidence_upper": float(
                                np.exp(
                                    confidence.iloc[1]
                                )
                            ),
                            "p_value": float(
                                result.pvalues[exposure]
                            ),
                            "direction": (
                                "lower_eGFR"
                                if coefficient < 0
                                else "higher_eGFR"
                            ),
                        }
                    )
                except Exception as error:
                    rows.append(
                        {
                            "gene": gene,
                            "scope": scope,
                            "model": (
                                f"model_{model_number}"
                            ),
                            "donors": len(subset),
                            "coefficient": np.nan,
                            "odds_ratio_per_10pct": (
                                np.nan
                            ),
                            "confidence_lower": np.nan,
                            "confidence_upper": np.nan,
                            "p_value": np.nan,
                            "direction": "not_estimated",
                            "error": str(error),
                        }
                    )
    results = pd.DataFrame(rows)
    results["FDR"] = np.nan
    group_columns = ["scope", "model"]
    for keys, group in results.groupby(
        group_columns,
        observed=True,
    ):
        mask = (
            (results["scope"] == keys[0])
            & (results["model"] == keys[1])
        )
        results.loc[mask, "FDR"] = (
            benjamini_hochberg(
                results.loc[mask, "p_value"]
            )
        )
    return results


def main():
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    data = build_gene_data()
    results = fit_models(data)
    data.to_csv(
        MODEL_DATA_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    results.to_csv(
        OUTPUT_PATH,
        index=False,
        encoding="utf-8-sig",
    )
    print("Saved:", MODEL_DATA_OUTPUT)
    print("Saved:", OUTPUT_PATH)
    print("\neGFR gene specificity:")
    print(
        results[
            [
                "gene",
                "scope",
                "model",
                "donors",
                "odds_ratio_per_10pct",
                "confidence_lower",
                "confidence_upper",
                "p_value",
                "FDR",
                "direction",
            ]
        ]
        .round(4)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
