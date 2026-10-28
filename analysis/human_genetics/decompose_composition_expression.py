"""Decompose overall TAL expression differences into composition and expression.

The donor is the unit of analysis. For each donor, the overall TAL expression
is the cell-weighted mean across the four TAL subtypes:

    overall_d = sum_s proportion_ds * mean_expression_ds

The disease-reference difference is separated into:

1. composition: changes in subtype proportions at reference state means;
2. expression: changes in subtype-specific means at reference proportions;
3. remainder: interaction plus differences in the within-group covariance
   between subtype proportion and subtype mean.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


INPUT = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2\outputs"
    r"\KPMP_full_TAL_focus_expression.csv.gz"
)
OUT_DIR = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji\outputs"
    r"\composition_expression_decomposition"
)
GENES = ["PAPPA2", "CXCL12", "NOTCH3", "HES1", "EPHB2"]
STATES = ["C-TAL-A", "C-TAL-B", "aTAL2", "frTAL"]
COMPARISONS = {
    "all_disease": ["AKI", "CKD"],
    "AKI_vs_Reference": ["AKI"],
    "CKD_vs_Reference": ["CKD"],
}
N_BOOTSTRAP = 2000
SEED = 20261028


def build_donor_state_table(expression: pd.DataFrame) -> pd.DataFrame:
    expression = expression[
        expression["subclass"].isin(STATES)
        & expression["condition"].isin(["Reference", "AKI", "CKD"])
    ].copy()
    value_columns = [f"{gene}_lognorm" for gene in GENES]
    donor_state = (
        expression.groupby(
            ["donor_id", "condition", "subclass"],
            observed=True,
        )
        .agg(
            n_cells=("cell_id", "size"),
            **{column: (column, "mean") for column in value_columns},
        )
        .reset_index()
    )
    donor_totals = (
        donor_state.groupby(["donor_id", "condition"], observed=True)[
            "n_cells"
        ]
        .sum()
        .rename("n_tal_cells")
        .reset_index()
    )
    donor_state = donor_state.merge(
        donor_totals,
        on=["donor_id", "condition"],
        how="left",
        validate="many_to_one",
    )
    donor_state["proportion"] = (
        donor_state["n_cells"] / donor_state["n_tal_cells"]
    )
    for gene in GENES:
        donor_state[f"{gene}_contribution"] = (
            donor_state["proportion"]
            * donor_state[f"{gene}_lognorm"]
        )
    return donor_state


def group_moments(
    donor_state: pd.DataFrame,
    donors: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    block = donor_state[donor_state["donor_id"].isin(donors)]
    proportions = {}
    means = {}
    for state in STATES:
        state_block = block[block["subclass"].eq(state)]
        present = state_block.drop_duplicates("donor_id")
        total_donors = block["donor_id"].nunique()
        denominator = total_donors if total_donors else np.nan
        proportions[state] = (
            present["proportion"].sum() / denominator
            if denominator
            else np.nan
        )
        means[state] = (
            present[
                [f"{gene}_lognorm" for gene in GENES]
            ].mean()
            if not present.empty
            else pd.Series(np.nan, index=[f"{gene}_lognorm" for gene in GENES])
        )
    proportion_array = {
        state: float(proportions[state]) for state in STATES
    }
    mean_array = {
        gene: np.array(
            [
                means[state][f"{gene}_lognorm"]
                for state in STATES
            ],
            dtype=float,
        )
        for gene in GENES
    }
    return proportion_array, mean_array


def calculate_components(
    reference_state: pd.DataFrame,
    disease_state: pd.DataFrame,
    reference_donors: np.ndarray,
    disease_donors: np.ndarray,
) -> tuple[dict[str, dict[str, float]], pd.DataFrame]:
    reference_props, reference_means = group_moments(
        reference_state,
        reference_donors,
    )
    disease_props, disease_means = group_moments(
        disease_state,
        disease_donors,
    )
    reference_matrix = (
        reference_state[
            reference_state["donor_id"].isin(reference_donors)
        ]
        .groupby("donor_id", observed=True)[
            [f"{gene}_contribution" for gene in GENES]
        ]
        .sum()
    )
    disease_matrix = (
        disease_state[
            disease_state["donor_id"].isin(disease_donors)
        ]
        .groupby("donor_id", observed=True)[
            [f"{gene}_contribution" for gene in GENES]
        ]
        .sum()
    )

    components: dict[str, dict[str, float]] = {}
    state_rows = []
    for gene in GENES:
        reference_prop_vector = np.array(
            [reference_props[state] for state in STATES],
            dtype=float,
        )
        disease_prop_vector = np.array(
            [disease_props[state] for state in STATES],
            dtype=float,
        )
        reference_mean_vector = reference_means[gene]
        disease_mean_vector = disease_means[gene]

        composition = np.nansum(
            (disease_prop_vector - reference_prop_vector)
            * reference_mean_vector
        )
        expression = np.nansum(
            reference_prop_vector
            * (disease_mean_vector - reference_mean_vector)
        )
        observed = (
            disease_matrix[f"{gene}_contribution"].mean()
            - reference_matrix[f"{gene}_contribution"].mean()
        )
        remainder = observed - composition - expression
        components[gene] = {
            "observed_difference": observed,
            "composition_effect": composition,
            "expression_effect": expression,
            "interaction_and_covariance": remainder,
            "reference_proportion_total": reference_prop_vector.sum(),
            "disease_proportion_total": disease_prop_vector.sum(),
        }

        for index, state in enumerate(STATES):
            state_rows.append(
                {
                    "gene": gene,
                    "state": state,
                    "reference_proportion": reference_prop_vector[index],
                    "disease_proportion": disease_prop_vector[index],
                    "proportion_difference": (
                        disease_prop_vector[index]
                        - reference_prop_vector[index]
                    ),
                    "reference_mean": reference_mean_vector[index],
                    "disease_mean": disease_mean_vector[index],
                    "mean_difference": (
                        disease_mean_vector[index]
                        - reference_mean_vector[index]
                    ),
                    "composition_contribution": (
                        disease_prop_vector[index]
                        - reference_prop_vector[index]
                    )
                    * reference_mean_vector[index],
                    "expression_contribution": (
                        reference_prop_vector[index]
                        * (
                            disease_mean_vector[index]
                            - reference_mean_vector[index]
                        )
                    ),
                }
            )
    return components, pd.DataFrame(state_rows)


def bootstrap_components(
    reference_state: pd.DataFrame,
    disease_state: pd.DataFrame,
    reference_donors: np.ndarray,
    disease_donors: np.ndarray,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows = []
    for _ in range(N_BOOTSTRAP):
        sampled_reference = rng.choice(
            reference_donors,
            size=len(reference_donors),
            replace=True,
        )
        sampled_disease = rng.choice(
            disease_donors,
            size=len(disease_donors),
            replace=True,
        )
        components, _ = calculate_components(
            reference_state,
            disease_state,
            sampled_reference,
            sampled_disease,
        )
        for gene, values in components.items():
            rows.append(
                {
                    "gene": gene,
                    **values,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expression = pd.read_csv(
        INPUT,
        usecols=[
            "cell_id",
            "donor_id",
            "subclass",
            "condition",
            *[f"{gene}_lognorm" for gene in GENES],
        ],
    )
    donor_state = build_donor_state_table(expression)
    donor_state.to_csv(
        OUT_DIR / "donor_state_proportions_and_means.csv",
        index=False,
    )

    summary_rows = []
    state_contribution_frames = []
    rng = np.random.default_rng(SEED)
    for comparison, case_conditions in COMPARISONS.items():
        reference_state = donor_state[
            donor_state["condition"].eq("Reference")
        ]
        disease_state = donor_state[
            donor_state["condition"].isin(case_conditions)
        ]
        reference_donors = reference_state["donor_id"].unique()
        disease_donors = disease_state["donor_id"].unique()

        components, state_contributions = calculate_components(
            reference_state,
            disease_state,
            reference_donors,
            disease_donors,
        )
        bootstrap = bootstrap_components(
            reference_state,
            disease_state,
            reference_donors,
            disease_donors,
            rng,
        )
        state_contributions.insert(0, "comparison", comparison)
        state_contribution_frames.append(state_contributions)

        for gene, values in components.items():
            gene_bootstrap = bootstrap[bootstrap["gene"].eq(gene)]
            row = {
                "comparison": comparison,
                "gene": gene,
                **values,
                "n_reference_donors": len(reference_donors),
                "n_disease_donors": len(disease_donors),
            }
            for component in [
                "observed_difference",
                "composition_effect",
                "expression_effect",
                "interaction_and_covariance",
            ]:
                row[f"{component}_ci_lower"] = gene_bootstrap[
                    component
                ].quantile(0.025)
                row[f"{component}_ci_upper"] = gene_bootstrap[
                    component
                ].quantile(0.975)
            summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    summary["absolute_component_sum"] = (
        summary["composition_effect"].abs()
        + summary["expression_effect"].abs()
        + summary["interaction_and_covariance"].abs()
    )
    summary["composition_share_of_absolute_components"] = (
        summary["composition_effect"].abs()
        / summary["absolute_component_sum"].replace(0, np.nan)
    )
    summary["expression_share_of_absolute_components"] = (
        summary["expression_effect"].abs()
        / summary["absolute_component_sum"].replace(0, np.nan)
    )
    summary.to_csv(
        OUT_DIR / "composition_expression_summary.csv",
        index=False,
    )
    pd.concat(state_contribution_frames, ignore_index=True).to_csv(
        OUT_DIR / "composition_expression_state_contributions.csv",
        index=False,
    )

    print(summary[
        [
            "comparison",
            "gene",
            "observed_difference",
            "composition_effect",
            "expression_effect",
            "interaction_and_covariance",
            "composition_share_of_absolute_components",
            "expression_share_of_absolute_components",
        ]
    ].to_string(index=False))


if __name__ == "__main__":
    main()
