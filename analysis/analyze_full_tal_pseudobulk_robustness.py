import numpy as np
import pandas as pd

from project_paths import OUTPUT_ROOT

EXPRESSION_INPUT = (
    OUTPUT_ROOT / "KPMP_full_TAL_focus_expression.csv.gz"
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
BOOTSTRAP_REPLICATES = 1000
BOOTSTRAP_SEED = 20260930

LODO_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_focus_pseudobulk_lodo.csv"
)
BOOTSTRAP_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_focus_pseudobulk_bootstrap.csv"
)
REGION_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_focus_pseudobulk_region_strata.csv"
)
SUMMARY_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_focus_pseudobulk_robustness_summary.md"
)


def build_donor_pseudobulk(expression):
    expression = expression[
        expression["condition"].isin(
            ["Reference", "AKI", "CKD"]
        )
    ].copy()
    gene_columns = [
        f"{gene}_lognorm" for gene in FOCUS_GENES
    ]
    donor = (
        expression.groupby(
            ["donor_id", "condition", "subclass"],
            observed=True,
        )[gene_columns]
        .mean()
        .reset_index()
    )
    donor = donor.rename(
        columns={
            column: column.removesuffix("_lognorm")
            for column in gene_columns
        }
    )
    return donor


def mean_difference(reference, case, gene):
    reference_values = reference[gene].to_numpy(dtype=float)
    case_values = case[gene].to_numpy(dtype=float)
    return (
        np.mean(case_values)
        - np.mean(reference_values)
    )


def leave_one_donor_out(pseudobulk):
    rows = []
    for subclass in SUBCLASS_ORDER:
        subclass_data = pseudobulk[
            pseudobulk["subclass"] == subclass
        ]
        donors = subclass_data["donor_id"].unique()
        for comparison, case_conditions in COMPARISONS.items():
            for gene in FOCUS_GENES:
                full_reference = subclass_data[
                    subclass_data["condition"]
                    == "Reference"
                ]
                full_case = subclass_data[
                    subclass_data["condition"].isin(
                        case_conditions
                    )
                ]
                full_effect = mean_difference(
                    full_reference,
                    full_case,
                    gene,
                )
                lodo_effects = []
                for donor in donors:
                    reduced = subclass_data[
                        subclass_data["donor_id"] != donor
                    ]
                    reference = reduced[
                        reduced["condition"] == "Reference"
                    ]
                    case = reduced[
                        reduced["condition"].isin(
                            case_conditions
                        )
                    ]
                    if len(reference) == 0 or len(case) == 0:
                        continue
                    lodo_effects.append(
                        mean_difference(
                            reference,
                            case,
                            gene,
                        )
                    )
                lodo_effects = np.asarray(
                    lodo_effects,
                    dtype=float,
                )
                rows.append(
                    {
                        "comparison": comparison,
                        "subclass": subclass,
                        "gene": gene,
                        "full_effect": full_effect,
                        "donors_testable": len(lodo_effects),
                        "lodo_median": np.median(lodo_effects),
                        "lodo_min": np.min(lodo_effects),
                        "lodo_max": np.max(lodo_effects),
                        "lodo_same_sign_fraction": np.mean(
                            np.sign(lodo_effects)
                            == np.sign(full_effect)
                        ),
                        "lodo_max_abs_change": np.max(
                            np.abs(lodo_effects - full_effect)
                        ),
                    }
                )
    return pd.DataFrame(rows)


def bootstrap_donor_effects(pseudobulk):
    generator = np.random.default_rng(BOOTSTRAP_SEED)
    rows = []
    for subclass in SUBCLASS_ORDER:
        subclass_data = pseudobulk[
            pseudobulk["subclass"] == subclass
        ]
        for comparison, case_conditions in COMPARISONS.items():
            case_data = subclass_data[
                subclass_data["condition"].isin(
                    case_conditions
                )
            ]
            reference_data = subclass_data[
                subclass_data["condition"] == "Reference"
            ]
            case_strata = [
                case_data[
                    case_data["condition"] == condition
                ]
                for condition in case_conditions
            ]
            case_strata = [
                stratum
                for stratum in case_strata
                if len(stratum) > 0
            ]
            reference_values = reference_data
            if (
                len(reference_values) == 0
                or len(case_strata) == 0
            ):
                continue

            for gene in FOCUS_GENES:
                effects = []
                for _ in range(BOOTSTRAP_REPLICATES):
                    reference_sample = reference_values.sample(
                        n=len(reference_values),
                        replace=True,
                        random_state=int(
                            generator.integers(0, 2**31 - 1)
                        ),
                    )
                    case_samples = [
                        stratum.sample(
                            n=len(stratum),
                            replace=True,
                            random_state=int(
                                generator.integers(
                                    0,
                                    2**31 - 1,
                                )
                            ),
                        )
                        for stratum in case_strata
                    ]
                    case_sample = pd.concat(
                        case_samples,
                        ignore_index=True,
                    )
                    effects.append(
                        np.mean(
                            case_sample[gene].to_numpy(
                                dtype=float
                            )
                        )
                        - np.mean(
                            reference_sample[gene].to_numpy(
                                dtype=float
                            )
                        )
                    )
                effects = np.asarray(effects, dtype=float)
                rows.append(
                    {
                        "comparison": comparison,
                        "subclass": subclass,
                        "gene": gene,
                        "bootstrap_replicates": len(effects),
                        "bootstrap_median": np.median(effects),
                        "bootstrap_ci_lower": np.percentile(
                            effects,
                            2.5,
                        ),
                        "bootstrap_ci_upper": np.percentile(
                            effects,
                            97.5,
                        ),
                        "bootstrap_fraction_positive": np.mean(
                            effects > 0
                        ),
                    }
                )
    return pd.DataFrame(rows)


def region_stratified_effects(expression):
    gene_columns = [
        f"{gene}_lognorm" for gene in FOCUS_GENES
    ]
    region_data = expression[
        expression["condition"].isin(
            ["Reference", "AKI", "CKD"]
        )
        & expression["region"].isin(
            ["Cortex", "Medulla"]
        )
    ].copy()
    donor_region = (
        region_data.groupby(
            [
                "region",
                "donor_id",
                "condition",
                "subclass",
            ],
            observed=True,
        )[gene_columns]
        .mean()
        .reset_index()
    )
    donor_region = donor_region.rename(
        columns={
            column: column.removesuffix("_lognorm")
            for column in gene_columns
        }
    )
    rows = []
    for region in ["Cortex", "Medulla"]:
        region_subset = donor_region[
            donor_region["region"] == region
        ]
        for subclass in SUBCLASS_ORDER:
            subclass_data = region_subset[
                region_subset["subclass"] == subclass
            ]
            reference = subclass_data[
                subclass_data["condition"] == "Reference"
            ]
            for comparison, case_conditions in COMPARISONS.items():
                case = subclass_data[
                    subclass_data["condition"].isin(
                        case_conditions
                    )
                ]
                for gene in FOCUS_GENES:
                    if len(reference) == 0 or len(case) == 0:
                        effect = np.nan
                    else:
                        effect = mean_difference(
                            reference,
                            case,
                            gene,
                        )
                    rows.append(
                        {
                            "region": region,
                            "comparison": comparison,
                            "subclass": subclass,
                            "gene": gene,
                            "reference_donors": len(reference),
                            "case_donors": len(case),
                            "mean_case_minus_reference": effect,
                        }
                    )
    return pd.DataFrame(rows)


def build_summary(lodo, bootstrap, region):
    pappa2_lodo = lodo[
        lodo["gene"] == "PAPPA2"
    ]
    pappa2_bootstrap = bootstrap[
        bootstrap["gene"] == "PAPPA2"
    ]
    region_pivot = region.pivot_table(
        index=[
            "comparison",
            "subclass",
            "gene",
        ],
        columns="region",
        values="mean_case_minus_reference",
        observed=True,
    ).reset_index()
    region_pivot["region_direction_agreement"] = (
        np.sign(region_pivot["Cortex"])
        == np.sign(region_pivot["Medulla"])
    )
    pappa2_region = region_pivot[
        region_pivot["gene"] == "PAPPA2"
    ]

    def render(frame):
        lines = [
            "| " + " | ".join(frame.columns) + " |",
            "|"
            + "|".join(["---"] * len(frame.columns))
            + "|",
        ]
        for _, row in frame.iterrows():
            values = []
            for value in row:
                if pd.isna(value):
                    values.append("")
                elif isinstance(value, float):
                    values.append(f"{value:.4f}")
                else:
                    values.append(str(value))
            lines.append("| " + " | ".join(values) + " |")
        return "\n".join(lines)

    return f"""# 完整 TAL 供者伪批量稳健性汇总

日期：2026-09-30

## 一、留一供者

{render(pappa2_lodo)}

## 二、供者分层 bootstrap

{render(pappa2_bootstrap)}

## 三、皮质和髓质分层

{render(pappa2_region)}

## 四、当前判断

- 完整 TAL 供者伪批量的方向需要在留一供者后保持；
- bootstrap 区间用于描述供者层面稳定性，不是独立细胞检验；
- 皮质和髓质分层可以检查区域组成是否解释效应；
- 如果区域方向不一致，应写成区域相关状态，而不是全器官统一变化。
"""


def main():
    expression = pd.read_csv(EXPRESSION_INPUT)
    pseudobulk = build_donor_pseudobulk(expression)
    lodo = leave_one_donor_out(pseudobulk)
    bootstrap = bootstrap_donor_effects(pseudobulk)
    region = region_stratified_effects(expression)

    lodo.to_csv(
        LODO_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    bootstrap.to_csv(
        BOOTSTRAP_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    region.to_csv(
        REGION_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    SUMMARY_OUTPUT.write_text(
        build_summary(lodo, bootstrap, region),
        encoding="utf-8",
    )

    print("Saved:", LODO_OUTPUT)
    print("Saved:", BOOTSTRAP_OUTPUT)
    print("Saved:", REGION_OUTPUT)
    print("Saved:", SUMMARY_OUTPUT)
    print("\nPAPPA2 leave-one-donor-out:")
    print(
        lodo[lodo["gene"] == "PAPPA2"][
            [
                "comparison",
                "subclass",
                "full_effect",
                "lodo_median",
                "lodo_min",
                "lodo_max",
                "lodo_same_sign_fraction",
            ]
        ]
        .round(4)
        .to_string(index=False)
    )
    print("\nPAPPA2 region effects:")
    region_pivot = region.pivot_table(
        index=[
            "comparison",
            "subclass",
            "gene",
        ],
        columns="region",
        values="mean_case_minus_reference",
        observed=True,
    ).reset_index()
    print(
        region_pivot[
            region_pivot["gene"] == "PAPPA2"
        ]
        .round(4)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
