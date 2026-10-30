"""Core analysis-unit functions."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


def _bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(dtype=float)
    result = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return pd.Series(result, index=values.index)
    order = np.argsort(p[valid])
    ranked = p[valid][order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    adjusted = np.empty_like(q)
    adjusted[order] = np.clip(q, 0, 1)
    result[valid] = adjusted
    return pd.Series(result, index=values.index)


def _safe_mannwhitney(
    case_values: np.ndarray,
    reference_values: np.ndarray,
) -> float:
    case_values = np.asarray(case_values, dtype=float)
    reference_values = np.asarray(reference_values, dtype=float)
    case_values = case_values[np.isfinite(case_values)]
    reference_values = reference_values[np.isfinite(reference_values)]
    if (
        len(case_values) < 2
        or len(reference_values) < 2
        or np.ptp(case_values) == 0
        or np.ptp(reference_values) == 0
    ):
        return np.nan
    return float(
        mannwhitneyu(
            case_values,
            reference_values,
            alternative="two-sided",
        ).pvalue
    )


@dataclass
class UnitResults:
    donor: pd.DataFrame
    neighborhood: pd.DataFrame
    within_donor: pd.DataFrame
    reversals: pd.DataFrame


def _aligned(
    expression: pd.DataFrame,
    cells: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not expression.index.is_unique:
        raise ValueError("expression index must be unique")
    cells = cells.copy()
    if "cell_id" in cells.columns:
        cells = cells.set_index("cell_id")
    common = expression.index.intersection(cells.index)
    if len(common) != len(cells):
        raise ValueError("cell metadata and expression index do not match")
    return expression.loc[common], cells.loc[common]


def donor_pseudobulk(
    expression: pd.DataFrame,
    cells: pd.DataFrame,
    *,
    state_col: str = "state",
    donor_col: str = "donor_id",
    condition_col: str = "condition",
    case_levels: list[str],
    reference_level: str,
    comparisons: dict[str, list[str]] | None = None,
    min_donors: int = 3,
) -> pd.DataFrame:
    expression, cells = _aligned(expression, cells)
    if comparisons is None:
        comparisons = {"disease_vs_reference": list(case_levels)}
    rows = []
    for state, state_cells in cells.groupby(state_col, observed=True):
        donor_means = (
            expression.loc[state_cells.index].groupby(
                state_cells[donor_col],
                observed=True,
            )
            .mean()
        )
        donor_condition = (
            state_cells.groupby(donor_col, observed=True)[condition_col]
            .agg(lambda values: values.mode().iloc[0])
        )
        for comparison, levels in comparisons.items():
            case_donors = donor_condition[
                donor_condition.isin(levels)
            ].index
            reference_donors = donor_condition[
                donor_condition.eq(reference_level)
            ].index
            if (
                len(case_donors) < min_donors
                or len(reference_donors) < min_donors
            ):
                continue
            for gene in expression.columns:
                case_values = donor_means.loc[case_donors, gene]
                reference_values = donor_means.loc[reference_donors, gene]
                rows.append(
                    {
                        "gene": gene,
                        "state": state,
                        "comparison": comparison,
                        "donor_effect": (
                            case_values.mean()
                            - reference_values.mean()
                        ),
                        "donor_p": _safe_mannwhitney(
                            case_values,
                            reference_values,
                        ),
                        "n_case_donors": len(case_values),
                        "n_reference_donors": len(reference_values),
                    }
                )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["donor_fdr"] = _bh_adjust(result["donor_p"])
    return result


def neighborhood_expression(
    expression: pd.DataFrame,
    cells: pd.DataFrame,
    neighborhood_logfc: pd.DataFrame,
    *,
    state_col: str = "state",
    neighborhood_col: str = "neighborhood_id",
    min_neighborhoods: int = 10,
) -> pd.DataFrame:
    expression, cells = _aligned(expression, cells)
    rows = []
    for state in cells[state_col].drop_duplicates():
        state_cells = cells[cells[state_col].eq(state)]
        merged = expression.loc[state_cells.index].join(
            state_cells[[neighborhood_col]]
        )
        means = merged.groupby(neighborhood_col, observed=True).mean()
        metadata = neighborhood_logfc[
            neighborhood_logfc[state_col].eq(state)
        ].set_index(neighborhood_col)
        for comparison, group in metadata.groupby(
            "comparison",
            observed=True,
        ):
            common = means.index.intersection(group.index)
            if len(common) < min_neighborhoods:
                continue
            logfc = group.loc[common, "logfc"]
            for gene in expression.columns:
                rho = spearmanr(means.loc[common, gene], logfc).statistic
                rows.append(
                    {
                        "gene": gene,
                        "state": state,
                        "comparison": comparison,
                        "neighborhood_effect": rho,
                        "n_neighborhoods": len(common),
                    }
                )
    return pd.DataFrame(rows)


def within_donor(
    expression: pd.DataFrame,
    cells: pd.DataFrame,
    *,
    exposure_col: str = "exposure",
    state_col: str = "state",
    donor_col: str = "donor_id",
    condition_col: str = "condition",
    case_levels: list[str],
    reference_level: str,
    comparisons: dict[str, list[str]] | None = None,
    min_cells: int = 10,
) -> pd.DataFrame:
    expression, cells = _aligned(expression, cells)
    if comparisons is None:
        comparisons = {"disease_vs_reference": list(case_levels)}
    rows = []
    for state, state_cells in cells.groupby(state_col, observed=True):
        donor_rhos = []
        for donor, donor_cells in state_cells.groupby(
            donor_col,
            observed=True,
        ):
            valid = (
                donor_cells[exposure_col].notna()
                & np.isfinite(donor_cells[exposure_col])
            )
            donor_cells = donor_cells[valid]
            if (
                len(donor_cells) < min_cells
                or donor_cells[exposure_col].nunique() < 2
            ):
                continue
            block = expression.loc[donor_cells.index]
            exposure = donor_cells[exposure_col].to_numpy(dtype=float)
            rho = {
                gene: spearmanr(
                    block[gene].to_numpy(dtype=float),
                    exposure,
                ).statistic
                for gene in expression.columns
            }
            donor_rhos.append(
                {
                    "donor_id": donor,
                    "condition": donor_cells[condition_col].iloc[0],
                    **rho,
                }
            )
        if not donor_rhos:
            continue
        rho_table = pd.DataFrame(donor_rhos)
        for comparison, levels in comparisons.items():
            case_table = rho_table[rho_table["condition"].isin(levels)]
            reference_table = rho_table[
                rho_table["condition"].eq(reference_level)
            ]
            if (
                len(case_table) < 3
                or len(reference_table) < 3
            ):
                continue
            for gene in expression.columns:
                case_values = case_table[gene].to_numpy(dtype=float)
                reference_values = reference_table[gene].to_numpy(dtype=float)
                rows.append(
                    {
                        "gene": gene,
                        "state": state,
                        "comparison": comparison,
                        "within_donor_effect": (
                            np.nanmedian(case_values)
                            - np.nanmedian(reference_values)
                        ),
                        "within_donor_p": _safe_mannwhitney(
                            case_values,
                            reference_values,
                        ),
                        "n_case_donors": len(case_values),
                        "n_reference_donors": len(reference_values),
                    }
                )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["within_donor_fdr"] = _bh_adjust(
            result["within_donor_p"]
        )
    return result


def standardize_and_score(
    donor: pd.DataFrame,
    neighborhood: pd.DataFrame,
    within: pd.DataFrame,
    *,
    z_threshold: float = 0.5,
) -> pd.DataFrame:
    merged = donor.merge(
        neighborhood,
        on=["gene", "state", "comparison"],
        how="outer",
        validate="one_to_one",
    ).merge(
        within,
        on=["gene", "state", "comparison"],
        how="outer",
        validate="one_to_one",
    )
    effect_columns = [
        "donor_effect",
        "neighborhood_effect",
        "within_donor_effect",
    ]
    for effect in effect_columns:
        merged[f"{effect}_z"] = np.nan
        for _, group in merged.groupby(
            ["state", "comparison"],
            observed=True,
        ):
            values = group[effect].to_numpy(dtype=float)
            valid = np.isfinite(values)
            if valid.sum() < 3:
                continue
            scale = np.nanstd(values[valid])
            if not np.isfinite(scale) or scale == 0:
                continue
            merged.loc[group.index, f"{effect}_z"] = (
                values / scale
            )
    pairs = {
        "donor_within": ("donor_effect_z", "within_donor_effect_z"),
        "donor_neighborhood": (
            "donor_effect_z",
            "neighborhood_effect_z",
        ),
        "neighborhood_within": (
            "neighborhood_effect_z",
            "within_donor_effect_z",
        ),
    }
    for name, (left, right) in pairs.items():
        valid = merged[[left, right]].notna().all(axis=1)
        merged[f"{name}_valid"] = valid
        merged[f"{name}_reversal"] = (
            valid
            & (
                np.sign(merged[left]) * np.sign(merged[right])
                < 0
            )
            & (merged[left].abs() >= z_threshold)
            & (merged[right].abs() >= z_threshold)
        )
    return merged


def reversal_fraction(
    effects: pd.DataFrame,
    *,
    left_col: str,
    right_col: str,
    z_threshold: float = 0.5,
) -> float:
    """Return the fraction of finite pairs with opposite, large effects."""
    if z_threshold < 0:
        raise ValueError("z_threshold must be non-negative")
    missing = [
        column
        for column in (left_col, right_col)
        if column not in effects.columns
    ]
    if missing:
        raise ValueError(f"Missing effect columns: {missing}")
    valid = effects[[left_col, right_col]].notna().all(axis=1)
    if not valid.any():
        return np.nan
    left = effects.loc[valid, left_col].to_numpy(dtype=float)
    right = effects.loc[valid, right_col].to_numpy(dtype=float)
    reversal = (
        (np.sign(left) * np.sign(right) < 0)
        & (np.abs(left) >= z_threshold)
        & (np.abs(right) >= z_threshold)
    )
    return float(reversal.mean())


def paired_reversal_permutation(
    effects: pd.DataFrame,
    *,
    left_col: str,
    right_col: str,
    target_genes: list[str] | set[str],
    gene_col: str = "gene",
    strata_cols: tuple[str, ...] = ("state", "comparison"),
    z_threshold: float = 0.5,
    n_permutations: int = 10000,
    seed: int = 20261028,
    aggregate: str = "gene",
) -> tuple[pd.DataFrame, dict]:
    """Test whether target gene pairs reverse more often than chance.

    The permutation preserves each target gene's left-unit effect and the
    marginal distribution of the right-unit effect within each stratum. It
    breaks only the gene-wise pairing between the two units. This conditions
    the null on the observed effect-size distribution and is therefore less
    sensitive to large-effect selection than an expression-only random set.
    """
    required = {gene_col, left_col, right_col, *strata_cols}
    missing = sorted(required.difference(effects.columns))
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    if n_permutations < 1:
        raise ValueError("n_permutations must be positive")
    if z_threshold < 0:
        raise ValueError("z_threshold must be non-negative")
    if aggregate not in {"gene", "row"}:
        raise ValueError("aggregate must be 'gene' or 'row'")

    target = set(target_genes)
    if not target:
        raise ValueError("target_genes must not be empty")

    frame = effects.copy()
    valid = frame[[left_col, right_col]].notna().all(axis=1)
    frame = frame.loc[valid]
    target_frame = frame[frame[gene_col].isin(target)]
    if target_frame.empty:
        raise ValueError("No finite target-gene effect pairs")

    rng = np.random.default_rng(seed)
    n_target_rows = 0
    if aggregate == "row":
        null_counts = np.zeros(n_permutations, dtype=np.int64)
        denominator = 0
        observed_counts = 0
    else:
        target_gene_names = sorted(target_frame[gene_col].unique())
        gene_position = {
            gene: index
            for index, gene in enumerate(target_gene_names)
        }
        gene_denominator = np.zeros(
            len(target_gene_names),
            dtype=np.int64,
        )
        gene_observed = np.zeros(
            len(target_gene_names),
            dtype=np.int64,
        )
        null_counts = np.zeros(
            (n_permutations, len(target_gene_names)),
            dtype=np.uint16,
        )

    for _, stratum in frame.groupby(
        list(strata_cols),
        observed=True,
        sort=False,
    ):
        target_rows = stratum[stratum[gene_col].isin(target)]
        if target_rows.empty:
            continue
        left = target_rows[left_col].to_numpy(dtype=float)
        right = target_rows[right_col].to_numpy(dtype=float)
        right_pool = stratum[right_col].to_numpy(dtype=float)
        n_target = len(target_rows)
        n_target_rows += n_target

        observed_flags = (
            (np.sign(left) * np.sign(right) < 0)
            & (np.abs(left) >= z_threshold)
            & (np.abs(right) >= z_threshold)
        )
        if aggregate == "row":
            denominator += n_target
            observed_counts += int(np.sum(observed_flags))
        else:
            positions = np.asarray(
                [
                    gene_position[gene]
                    for gene in target_rows[gene_col]
                ],
                dtype=np.int64,
            )
            gene_denominator[positions] += 1
            gene_observed[positions] += observed_flags.astype(np.int64)

        left_reversal_eligible = np.abs(left) >= z_threshold
        for permutation in range(n_permutations):
            sampled_right = rng.choice(
                right_pool,
                size=n_target,
                replace=False,
            )
            permutation_flags = (
                left_reversal_eligible
                & (np.sign(left) * np.sign(sampled_right) < 0)
                & (np.abs(sampled_right) >= z_threshold)
            )
            if aggregate == "row":
                null_counts[permutation] += int(
                    np.sum(permutation_flags)
                )
            else:
                np.add.at(
                    null_counts[permutation],
                    positions,
                    permutation_flags.astype(np.uint16),
                )

    if aggregate == "row":
        if denominator == 0:
            raise ValueError("No finite target-gene effect pairs")
        observed = observed_counts / denominator
        null_fractions = null_counts / denominator
    else:
        if not np.any(gene_denominator):
            raise ValueError("No finite target-gene effect pairs")
        observed = float(
            np.mean(gene_observed / gene_denominator)
        )
        null_fractions = np.mean(
            null_counts / gene_denominator[None, :],
            axis=1,
        )
    empirical_p = (
        1.0 + float(np.sum(null_fractions >= observed))
    ) / (n_permutations + 1.0)
    null_table = pd.DataFrame(
        {
            "permutation": np.arange(1, n_permutations + 1),
            "reversal_fraction": null_fractions,
        }
    )
    summary = {
        "left_effect": left_col,
        "right_effect": right_col,
        "z_threshold": float(z_threshold),
        "n_target_rows": int(n_target_rows),
        "n_permutations": int(n_permutations),
        "seed": int(seed),
        "aggregate": aggregate,
        "observed_reversal_fraction": observed,
        "null_mean": float(null_fractions.mean()),
        "null_sd": float(null_fractions.std(ddof=1))
        if n_permutations > 1
        else 0.0,
        "null_q025": float(np.quantile(null_fractions, 0.025)),
        "null_median": float(np.quantile(null_fractions, 0.5)),
        "null_q975": float(np.quantile(null_fractions, 0.975)),
        "empirical_p": empirical_p,
    }
    return null_table, summary


def analyze_units(
    expression: pd.DataFrame,
    cells: pd.DataFrame,
    *,
    case_levels: list[str],
    reference_level: str,
    neighborhood_logfc: pd.DataFrame,
    exposure_col: str = "exposure",
    state_col: str = "state",
    donor_col: str = "donor_id",
    condition_col: str = "condition",
    neighborhood_col: str = "neighborhood_id",
    comparisons: dict[str, list[str]] | None = None,
    z_threshold: float = 0.5,
) -> UnitResults:
    donor = donor_pseudobulk(
        expression,
        cells,
        state_col=state_col,
        donor_col=donor_col,
        condition_col=condition_col,
        case_levels=case_levels,
        reference_level=reference_level,
        comparisons=comparisons,
    )
    neighborhood = neighborhood_expression(
        expression,
        cells,
        neighborhood_logfc,
        state_col=state_col,
        neighborhood_col=neighborhood_col,
    )
    within = within_donor(
        expression,
        cells,
        exposure_col=exposure_col,
        state_col=state_col,
        donor_col=donor_col,
        condition_col=condition_col,
        case_levels=case_levels,
        reference_level=reference_level,
        comparisons=comparisons,
    )
    reversals = standardize_and_score(
        donor,
        neighborhood,
        within,
        z_threshold=z_threshold,
    )
    return UnitResults(
        donor=donor,
        neighborhood=neighborhood,
        within_donor=within,
        reversals=reversals,
    )


def decompose_composition(
    expression: pd.DataFrame,
    cells: pd.DataFrame,
    *,
    gene: str,
    state_col: str = "state",
    donor_col: str = "donor_id",
    condition_col: str = "condition",
    case_levels: list[str],
    reference_level: str,
) -> pd.DataFrame:
    expression, cells = _aligned(expression, cells)
    if gene not in expression.columns:
        raise ValueError(f"Missing gene: {gene}")
    frame = cells.copy()
    frame["expression"] = expression[gene]
    donor_state = (
        frame.groupby(
            [donor_col, condition_col, state_col],
            observed=True,
        )
        .agg(
            n_cells=("expression", "size"),
            mean_expression=("expression", "mean"),
        )
        .reset_index()
    )
    totals = (
        donor_state.groupby(
            [donor_col, condition_col],
            observed=True,
        )["n_cells"]
        .sum()
        .rename("n_total")
        .reset_index()
    )
    donor_state = donor_state.merge(
        totals,
        on=[donor_col, condition_col],
    )
    donor_state["proportion"] = (
        donor_state["n_cells"] / donor_state["n_total"]
    )
    donor_state["contribution"] = (
        donor_state["proportion"] * donor_state["mean_expression"]
    )
    reference = donor_state[
        donor_state[condition_col].eq(reference_level)
    ]
    disease = donor_state[
        donor_state[condition_col].isin(case_levels)
    ]
    states = sorted(donor_state[state_col].unique())
    ref_prop = (
        reference.groupby(state_col)["proportion"].mean()
        .reindex(states, fill_value=0)
    )
    dis_prop = (
        disease.groupby(state_col)["proportion"].mean()
        .reindex(states, fill_value=0)
    )
    ref_mean = (
        reference.groupby(state_col)["mean_expression"].mean()
        .reindex(states)
    )
    dis_mean = (
        disease.groupby(state_col)["mean_expression"].mean()
        .reindex(states)
    )
    observed = (
        disease.groupby(donor_col)["contribution"].sum().mean()
        - reference.groupby(donor_col)["contribution"].sum().mean()
    )
    composition = float(
        ((dis_prop - ref_prop) * ref_mean).sum()
    )
    expression_effect = float(
        (ref_prop * (dis_mean - ref_mean)).sum()
    )
    return pd.DataFrame(
        [
            {
                "gene": gene,
                "observed_difference": observed,
                "composition_effect": composition,
                "expression_effect": expression_effect,
                "interaction_and_covariance": (
                    observed - composition - expression_effect
                ),
            }
        ]
    )
