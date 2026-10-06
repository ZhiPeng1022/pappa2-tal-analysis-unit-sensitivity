import numpy as np
import pandas as pd
from scipy.stats import pearsonr, rankdata, spearmanr

from project_paths import OUTPUT_ROOT, WORK_ROOT

INPUT_EXPRESSION = (
    OUTPUT_ROOT / "KPMP_full_TAL_focus_expression.csv.gz"
)
MILO_PREFIX = "KPMP_milo_full_TAL_final"

COMPARISONS = [
    "all_disease",
    "AKI_vs_Reference",
    "CKD_vs_Reference",
]
FOCUS_GENES = [
    "PAPPA2",
    "CXCL12",
    "NOTCH3",
    "HES1",
    "EPHB2",
]

BOOTSTRAP_REPLICATES = 500
BOOTSTRAP_SEED = 20260927

SUMMARY_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_nhood_focus_robustness.csv"
)
BOOTSTRAP_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_nhood_focus_bootstrap_rhos.csv.gz"
)
LODO_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_nhood_focus_lodo_rhos.csv.gz"
)


def load_cell_expression():
    if not INPUT_EXPRESSION.exists():
        raise FileNotFoundError(INPUT_EXPRESSION)
    expression = pd.read_csv(
        INPUT_EXPRESSION,
        dtype={"cell_id": str, "donor_id": str},
    )
    if not expression["cell_id"].is_unique:
        raise ValueError("cell_id is not unique.")
    expression_columns = {
        gene: f"{gene}_lognorm" for gene in FOCUS_GENES
    }
    missing_columns = [
        column
        for column in expression_columns.values()
        if column not in expression.columns
    ]
    if missing_columns:
        raise ValueError(
            f"Missing focus-gene columns: {missing_columns}"
        )
    result = expression.set_index("cell_id")[
        list(expression_columns.values())
    ].rename(columns={
        column: gene
        for gene, column in expression_columns.items()
    })
    result.index.name = "cell_id"
    return result


def read_comparison_inputs(comparison):
    membership_path = OUTPUT_ROOT / (
        f"{MILO_PREFIX}_{comparison}_nhood_membership.csv.gz"
    )
    result_path = OUTPUT_ROOT / (
        f"{MILO_PREFIX}_{comparison}_nhood_results.csv"
    )
    membership = pd.read_csv(
        membership_path,
        dtype={"cell_id": str, "donor_id": str},
    )
    results = pd.read_csv(result_path)
    return membership, results


def build_comparison_edges(
    comparison,
    membership,
    results,
    cell_expression,
):
    results = results.rename(
        columns={
            "Nhood": "result_nhood",
            "subclass": "nhood_subclass",
        }
    )
    membership = membership.copy()
    membership["result_nhood"] = (
        membership["Nhood"].astype(int) + 1
    )

    if not np.allclose(
        membership["weight"].to_numpy(dtype=float),
        1.0,
    ):
        raise ValueError(
            "The sensitivity script requires unit nhood weights."
        )

    edges = membership.merge(
        cell_expression.reset_index(),
        on="cell_id",
        how="left",
        validate="many_to_one",
        indicator=True,
    )
    if (edges["_merge"] != "both").any():
        raise ValueError(
            f"{comparison}: some membership cells are missing "
            "from the expression matrix."
        )
    edges = edges.drop(columns="_merge")

    metadata_columns = [
        "result_nhood",
        "nhood_subclass",
        "logFC",
        "SpatialFDR",
        "subclass_fraction",
        "nhood_purity",
    ]
    edges = edges.merge(
        results[metadata_columns],
        on="result_nhood",
        how="inner",
        validate="many_to_one",
    )
    edges.insert(0, "comparison", comparison)
    return edges


def compute_nhood_rhos(means, log_fc):
    rhos = []
    for gene_index in range(means.shape[1]):
        x = means[:, gene_index]
        y = log_fc
        valid = np.isfinite(x) & np.isfinite(y)
        x_valid = x[valid]
        y_valid = y[valid]
        if (
            len(x_valid) < 3
            or np.ptp(x_valid) == 0
            or np.ptp(y_valid) == 0
        ):
            rhos.append(np.nan)
        else:
            rho, _ = spearmanr(x_valid, y_valid)
            rhos.append(rho)
    return np.asarray(rhos, dtype=float)


def partial_spearman(x, y, control):
    valid = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(control)
    )
    x = x[valid]
    y = y[valid]
    control = control[valid]
    if (
        len(x) < 4
        or np.ptp(x) == 0
        or np.ptp(y) == 0
        or np.ptp(control) == 0
    ):
        return np.nan

    x_rank = rankdata(x)
    y_rank = rankdata(y)
    control_rank = rankdata(control)

    x_coefficients = np.polyfit(
        control_rank,
        x_rank,
        1,
    )
    y_coefficients = np.polyfit(
        control_rank,
        y_rank,
        1,
    )
    x_residual = x_rank - np.polyval(
        x_coefficients,
        control_rank,
    )
    y_residual = y_rank - np.polyval(
        y_coefficients,
        control_rank,
    )
    if np.ptp(x_residual) == 0 or np.ptp(y_residual) == 0:
        return np.nan
    rho, _ = pearsonr(x_residual, y_residual)
    return rho


def build_group_matrices(edges):
    nhoods = np.asarray(
        sorted(
            edges["result_nhood"]
            .astype(int)
            .unique()
        ),
        dtype=int,
    )
    donors = np.asarray(
        sorted(edges["donor_id"].astype(str).unique()),
        dtype=object,
    )
    nhood_lookup = {
        nhood: index
        for index, nhood in enumerate(nhoods)
    }
    donor_lookup = {
        donor: index
        for index, donor in enumerate(donors)
    }

    edge_nhood = edges["result_nhood"].astype(int).map(
        nhood_lookup
    ).to_numpy(dtype=int)
    edge_donor = edges["donor_id"].astype(str).map(
        donor_lookup
    ).to_numpy(dtype=int)

    count_matrix = np.zeros(
        (len(nhoods), len(donors)),
        dtype=float,
    )
    np.add.at(count_matrix, (edge_nhood, edge_donor), 1.0)

    sum_matrices = []
    for gene in FOCUS_GENES:
        sum_matrix = np.zeros_like(count_matrix)
        np.add.at(
            sum_matrix,
            (edge_nhood, edge_donor),
            edges[gene].to_numpy(dtype=float),
        )
        sum_matrices.append(sum_matrix)

    donor_conditions = (
        edges[["donor_id", "condition"]]
        .astype({"donor_id": str, "condition": str})
        .drop_duplicates()
    )
    if donor_conditions["donor_id"].duplicated().any():
        raise ValueError(
            "A donor is assigned to more than one condition."
        )
    donor_condition_lookup = dict(
        zip(
            donor_conditions["donor_id"],
            donor_conditions["condition"],
        )
    )
    donor_condition = np.asarray(
        [
            donor_condition_lookup[str(donor)]
            for donor in donors
        ],
        dtype=object,
    )
    condition_indices = [
        np.flatnonzero(donor_condition == condition)
        for condition in sorted(
            set(donor_condition.tolist())
        )
    ]
    condition_counts = {
        condition: int(
            (donor_condition == condition).sum()
        )
        for condition in sorted(
            set(donor_condition.tolist())
        )
    }
    return (
        nhoods,
        donors,
        count_matrix,
        sum_matrices,
        condition_indices,
        condition_counts,
    )


def weighted_nhood_means(
    count_matrix,
    sum_matrices,
    weights,
):
    denominator = count_matrix @ weights
    means = np.full(
        (count_matrix.shape[0], len(sum_matrices)),
        np.nan,
        dtype=float,
    )
    valid = denominator > 0
    for gene_index, sum_matrix in enumerate(sum_matrices):
        numerator = sum_matrix @ weights
        means[valid, gene_index] = (
            numerator[valid] / denominator[valid]
        )
    return means, valid


def donor_weights_from_indices(
    n_donors,
    condition_indices,
    random_generator,
):
    weights = np.zeros(n_donors, dtype=float)
    for indices in condition_indices:
        probabilities = np.full(
            len(indices),
            1.0 / len(indices),
        )
        sampled = random_generator.multinomial(
            len(indices),
            probabilities,
        )
        weights[indices] = sampled
    return weights


def condition_stratified_rho(
    edges,
    nhoods,
    log_fc,
    group_column,
):
    condition_group = np.where(
        edges["condition"].astype(str) == "Reference",
        "Reference",
        "Disease",
    )
    condition_edges = edges.assign(
        condition_group=condition_group
    )
    counts = (
        condition_edges.groupby(
            ["result_nhood", "condition_group"],
            observed=True,
        )
        .size()
        .unstack("condition_group")
        .reindex(nhoods)
    )
    means = (
        condition_edges.groupby(
            ["result_nhood", "condition_group"],
            observed=True,
        )[FOCUS_GENES]
        .mean()
        .unstack("condition_group")
        .reindex(nhoods)
    )

    result = {}
    for gene in FOCUS_GENES:
        reference_values = means[
            (gene, "Reference")
        ].to_numpy(dtype=float)
        disease_values = means[
            (gene, "Disease")
        ].to_numpy(dtype=float)
        if group_column == "Reference":
            values = reference_values
        elif group_column == "Disease":
            values = disease_values
        else:
            values = disease_values - reference_values
        valid = np.isfinite(values) & np.isfinite(log_fc)
        if (
            valid.sum() < 3
            or np.ptp(values[valid]) == 0
            or np.ptp(log_fc[valid]) == 0
        ):
            result[gene] = np.nan
        else:
            rho, _ = spearmanr(
                values[valid],
                log_fc[valid],
            )
            result[gene] = rho

    return result, counts


def analyze_group(
    comparison,
    subclass,
    group_edges,
    bootstrap_replicates,
    random_generator,
):
    nhood_meta = (
        group_edges[
            [
                "result_nhood",
                "logFC",
                "SpatialFDR",
                "nhood_purity",
                "subclass_fraction",
            ]
        ]
        .drop_duplicates("result_nhood")
        .sort_values("result_nhood")
        .reset_index(drop=True)
    )
    (
        nhoods,
        donors,
        count_matrix,
        sum_matrices,
        condition_indices,
        condition_counts,
    ) = build_group_matrices(group_edges)

    if not np.array_equal(
        nhoods,
        nhood_meta["result_nhood"].to_numpy(dtype=int),
    ):
        raise ValueError(
            f"{comparison}/{subclass}: nhood order mismatch."
        )

    log_fc = nhood_meta["logFC"].to_numpy(dtype=float)
    full_weights = np.ones(len(donors), dtype=float)
    full_means, full_valid = weighted_nhood_means(
        count_matrix,
        sum_matrices,
        full_weights,
    )
    full_rhos = compute_nhood_rhos(full_means, log_fc)

    lodo_rows = []
    lodo_rhos = []
    for donor_index, donor in enumerate(donors):
        weights = np.ones(len(donors), dtype=float)
        weights[donor_index] = 0.0
        means, valid = weighted_nhood_means(
            count_matrix,
            sum_matrices,
            weights,
        )
        rhos = compute_nhood_rhos(means, log_fc)
        lodo_rhos.append(rhos)
        for gene_index, gene in enumerate(FOCUS_GENES):
            lodo_rows.append(
                {
                    "comparison": comparison,
                    "subclass": subclass,
                    "gene": gene,
                    "left_out_donor": donor,
                    "n_valid_nhoods": int(valid.sum()),
                    "rho": rhos[gene_index],
                }
            )
    lodo_rhos = np.asarray(lodo_rhos, dtype=float)

    bootstrap_rows = []
    bootstrap_rhos = []
    for replicate in range(1, bootstrap_replicates + 1):
        weights = donor_weights_from_indices(
            len(donors),
            condition_indices,
            random_generator,
        )
        means, valid = weighted_nhood_means(
            count_matrix,
            sum_matrices,
            weights,
        )
        minimum_valid_nhoods = max(
            3,
            int(np.ceil(0.5 * len(nhoods))),
        )
        if int(valid.sum()) < minimum_valid_nhoods:
            continue
        rhos = compute_nhood_rhos(means, log_fc)
        bootstrap_rhos.append(rhos)
        for gene_index, gene in enumerate(FOCUS_GENES):
            bootstrap_rows.append(
                {
                    "comparison": comparison,
                    "subclass": subclass,
                    "gene": gene,
                    "replicate": replicate,
                    "n_valid_nhoods": int(valid.sum()),
                    "rho": rhos[gene_index],
                }
            )
    bootstrap_rhos = np.asarray(
        bootstrap_rhos,
        dtype=float,
    )

    disease_fraction_counts = (
        condition_stratified_rho(
            group_edges,
            nhoods,
            log_fc,
            "difference",
        )[1]
    )
    counts_reference = disease_fraction_counts[
        "Reference"
    ].to_numpy(dtype=float)
    counts_disease = disease_fraction_counts[
        "Disease"
    ].to_numpy(dtype=float)
    both_conditions = (
        (counts_reference > 0)
        & (counts_disease > 0)
    )
    disease_fraction = np.full(
        len(nhoods),
        np.nan,
        dtype=float,
    )
    disease_fraction[both_conditions] = (
        counts_disease[both_conditions]
        / (
            counts_reference[both_conditions]
            + counts_disease[both_conditions]
        )
    )

    reference_rhos, _ = condition_stratified_rho(
        group_edges,
        nhoods,
        log_fc,
        "Reference",
    )
    disease_rhos, _ = condition_stratified_rho(
        group_edges,
        nhoods,
        log_fc,
        "Disease",
    )
    difference_rhos, _ = condition_stratified_rho(
        group_edges,
        nhoods,
        log_fc,
        "difference",
    )

    purity = nhood_meta[
        "nhood_purity"
    ].to_numpy(dtype=float)
    rows = []
    for gene_index, gene in enumerate(FOCUS_GENES):
        full_rho = full_rhos[gene_index]
        lodo_gene = lodo_rhos[:, gene_index]
        lodo_valid = np.isfinite(lodo_gene)
        if (
            np.isfinite(full_rho)
            and lodo_valid.sum() > 0
        ):
            lodo_same_sign_fraction = np.mean(
                np.sign(lodo_gene[lodo_valid])
                == np.sign(full_rho)
            )
            lodo_max_abs_change = np.max(
                np.abs(lodo_gene[lodo_valid] - full_rho)
            )
        else:
            lodo_same_sign_fraction = np.nan
            lodo_max_abs_change = np.nan

        bootstrap_gene = bootstrap_rhos[
            :, gene_index
        ]
        bootstrap_valid = bootstrap_gene[
            np.isfinite(bootstrap_gene)
        ]
        if len(bootstrap_valid) > 0:
            bootstrap_median = np.median(
                bootstrap_valid
            )
            bootstrap_ci_lower = np.percentile(
                bootstrap_valid,
                2.5,
            )
            bootstrap_ci_upper = np.percentile(
                bootstrap_valid,
                97.5,
            )
            bootstrap_positive_fraction = np.mean(
                bootstrap_valid > 0
            )
            bootstrap_gt_half_fraction = np.mean(
                bootstrap_valid > 0.5
            )
            bootstrap_lt_minus_half_fraction = np.mean(
                bootstrap_valid < -0.5
            )
        else:
            bootstrap_median = np.nan
            bootstrap_ci_lower = np.nan
            bootstrap_ci_upper = np.nan
            bootstrap_positive_fraction = np.nan
            bootstrap_gt_half_fraction = np.nan
            bootstrap_lt_minus_half_fraction = np.nan

        expression = full_means[:, gene_index]
        purity_rho = np.nan
        disease_fraction_rho = np.nan
        partial_purity_rho = partial_spearman(
            expression,
            log_fc,
            purity,
        )
        valid_purity = (
            np.isfinite(expression)
            & np.isfinite(purity)
        )
        if (
            valid_purity.sum() >= 3
            and np.ptp(expression[valid_purity]) > 0
            and np.ptp(purity[valid_purity]) > 0
        ):
            purity_rho, _ = spearmanr(
                expression[valid_purity],
                purity[valid_purity],
            )
        valid_disease_fraction = (
            np.isfinite(expression)
            & np.isfinite(disease_fraction)
        )
        if (
            valid_disease_fraction.sum() >= 3
            and np.ptp(
                expression[
                    valid_disease_fraction
                ]
            )
            > 0
            and np.ptp(
                disease_fraction[
                    valid_disease_fraction
                ]
            )
            > 0
        ):
            disease_fraction_rho, _ = spearmanr(
                expression[
                    valid_disease_fraction
                ],
                disease_fraction[
                    valid_disease_fraction
                ],
            )

        rows.append(
            {
                "comparison": comparison,
                "subclass": subclass,
                "gene": gene,
                "nhoods": len(nhoods),
                "unique_cells": int(
                    group_edges["cell_id"].nunique()
                ),
                "donors": len(donors),
                "donors_reference": condition_counts.get(
                    "Reference",
                    0,
                ),
                "donors_aki": condition_counts.get(
                    "AKI",
                    0,
                ),
                "donors_ckd": condition_counts.get(
                    "CKD",
                    0,
                ),
                "rho_full": full_rho,
                "lodo_valid_donors": int(
                    lodo_valid.sum()
                ),
                "lodo_median": (
                    np.nanmedian(lodo_gene)
                    if lodo_valid.sum() > 0
                    else np.nan
                ),
                "lodo_min": (
                    np.nanmin(lodo_gene)
                    if lodo_valid.sum() > 0
                    else np.nan
                ),
                "lodo_max": (
                    np.nanmax(lodo_gene)
                    if lodo_valid.sum() > 0
                    else np.nan
                ),
                "lodo_same_sign_fraction": (
                    lodo_same_sign_fraction
                ),
                "lodo_max_abs_change": (
                    lodo_max_abs_change
                ),
                "bootstrap_requested": bootstrap_replicates,
                "bootstrap_valid": len(
                    bootstrap_valid
                ),
                "bootstrap_median": bootstrap_median,
                "bootstrap_ci_lower": (
                    bootstrap_ci_lower
                ),
                "bootstrap_ci_upper": (
                    bootstrap_ci_upper
                ),
                "bootstrap_fraction_positive": (
                    bootstrap_positive_fraction
                ),
                "bootstrap_fraction_gt_0_5": (
                    bootstrap_gt_half_fraction
                ),
                "bootstrap_fraction_lt_minus_0_5": (
                    bootstrap_lt_minus_half_fraction
                ),
                "rho_reference_mean_vs_logFC": (
                    reference_rhos[gene]
                ),
                "rho_disease_mean_vs_logFC": (
                    disease_rhos[gene]
                ),
                "rho_condition_difference_vs_logFC": (
                    difference_rhos[gene]
                ),
                "rho_expression_vs_purity": purity_rho,
                "partial_rho_controlling_purity": (
                    partial_purity_rho
                ),
                "rho_expression_vs_disease_fraction": (
                    disease_fraction_rho
                ),
                "nhoods_with_both_conditions": int(
                    both_conditions.sum()
                ),
                "full_valid_nhoods": int(
                    full_valid.sum()
                ),
            }
        )

    summary = pd.DataFrame(rows)
    bootstrap = pd.DataFrame(bootstrap_rows)
    lodo = pd.DataFrame(lodo_rows)
    return summary, bootstrap, lodo


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    cell_expression = load_cell_expression()

    summary_tables = []
    bootstrap_tables = []
    lodo_tables = []
    for comparison_index, comparison in enumerate(
        COMPARISONS
    ):
        membership, results = read_comparison_inputs(
            comparison
        )
        edges = build_comparison_edges(
            comparison,
            membership,
            results,
            cell_expression,
        )
        subclasses = sorted(
            edges["nhood_subclass"].unique()
        )
        for subclass_index, subclass in enumerate(
            subclasses
        ):
            group_edges = edges[
                edges["nhood_subclass"] == subclass
            ].copy()
            group_seed = (
                BOOTSTRAP_SEED
                + comparison_index * 100
                + subclass_index
            )
            random_generator = (
                np.random.default_rng(group_seed)
            )
            print(
                f"Analyzing {comparison} / {subclass}: "
                f"{group_edges['result_nhood'].nunique()} "
                "nhoods",
                flush=True,
            )
            summary, bootstrap, lodo = analyze_group(
                comparison,
                subclass,
                group_edges,
                BOOTSTRAP_REPLICATES,
                random_generator,
            )
            summary_tables.append(summary)
            bootstrap_tables.append(bootstrap)
            lodo_tables.append(lodo)

    summary = pd.concat(
        summary_tables,
        ignore_index=True,
    )
    bootstrap = pd.concat(
        bootstrap_tables,
        ignore_index=True,
    )
    lodo = pd.concat(
        lodo_tables,
        ignore_index=True,
    )
    summary.to_csv(
        SUMMARY_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    bootstrap.to_csv(
        BOOTSTRAP_OUTPUT,
        index=False,
        compression="gzip",
        encoding="utf-8-sig",
    )
    lodo.to_csv(
        LODO_OUTPUT,
        index=False,
        compression="gzip",
        encoding="utf-8-sig",
    )

    print("\nSaved:", SUMMARY_OUTPUT)
    print("Saved:", BOOTSTRAP_OUTPUT)
    print("Saved:", LODO_OUTPUT)
    print(
        "\nPAPPA2 stability summary:"
    )
    pappa2 = summary[
        summary["gene"] == "PAPPA2"
    ].sort_values(["comparison", "subclass"])
    print(
        pappa2[
            [
                "comparison",
                "subclass",
                "nhoods",
                "donors",
                "rho_full",
                "lodo_same_sign_fraction",
                "bootstrap_ci_lower",
                "bootstrap_ci_upper",
                "bootstrap_fraction_positive",
            ]
        ]
        .round(3)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
