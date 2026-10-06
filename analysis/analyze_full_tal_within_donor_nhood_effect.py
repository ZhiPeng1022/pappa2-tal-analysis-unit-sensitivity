import numpy as np
import pandas as pd

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

MIN_CELLS_PER_DONOR_NHOOD = 3
MIN_NHOODS_PER_DONOR = 2
MIN_DONORS_PER_CONDITION = 5
PERMUTATIONS = 500
PERMUTATION_SEED = 20260927

SUMMARY_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_within_donor_nhood_effect.csv"
)
NULL_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_within_donor_nhood_null.csv.gz"
)
DIAGNOSTIC_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_within_donor_nhood_diagnostics.csv"
)

METRICS = [
    "reference_slope",
    "disease_slope",
    "disease_minus_reference",
    "reference_slope_purity_adjusted",
    "disease_slope_purity_adjusted",
    "disease_minus_reference_purity_adjusted",
]


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


def read_and_merge_edges(
    comparison,
    cell_expression,
):
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
    results = pd.read_csv(result_path).rename(
        columns={
            "Nhood": "result_nhood",
            "subclass": "nhood_subclass",
        }
    )
    membership["result_nhood"] = (
        membership["Nhood"].astype(int) + 1
    )
    if not np.allclose(
        membership["weight"].to_numpy(dtype=float),
        1.0,
    ):
        raise ValueError(
            "The analysis requires unit neighborhood weights."
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
            f"{comparison}: missing cells in expression matrix."
        )
    edges = edges.drop(columns="_merge")
    metadata = results[
        [
            "result_nhood",
            "nhood_subclass",
            "logFC",
            "SpatialFDR",
            "subclass_fraction",
            "nhood_purity",
        ]
    ]
    edges = edges.merge(
        metadata,
        on="result_nhood",
        how="inner",
        validate="many_to_one",
    )
    edges.insert(0, "comparison", comparison)
    return edges


def build_donor_nhood_table(edges):
    group_columns = [
        "comparison",
        "nhood_subclass",
        "donor_id",
        "condition",
        "result_nhood",
        "logFC",
        "nhood_purity",
        "subclass_fraction",
    ]
    aggregation = {
        "n_cells": ("cell_id", "size"),
    }
    aggregation.update(
        {
            gene: (gene, "mean")
            for gene in FOCUS_GENES
        }
    )
    table = (
        edges.groupby(
            group_columns,
            observed=True,
        )
        .agg(**aggregation)
        .reset_index()
    )
    return table


def weighted_center(values, weights, donor_codes, n_donors):
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        values = values[:, None]
    totals = np.bincount(
        donor_codes,
        weights=weights,
        minlength=n_donors,
    )
    centered_columns = []
    for column_index in range(values.shape[1]):
        sums = np.bincount(
            donor_codes,
            weights=(
                weights * values[:, column_index]
            ),
            minlength=n_donors,
        )
        means = np.divide(
            sums,
            totals,
            out=np.zeros_like(sums),
            where=totals > 0,
        )
        centered_columns.append(
            values[:, column_index]
            - means[donor_codes]
        )
    return np.column_stack(centered_columns)


def fit_within_donor_models(
    expression_values,
    neighborhood_log_fc,
    disease_indicator,
    purity,
    weights,
    donor_codes,
    n_donors,
):
    n_genes = expression_values.shape[1]
    y_centered = weighted_center(
        expression_values,
        weights,
        donor_codes,
        n_donors,
    )
    x_centered = weighted_center(
        neighborhood_log_fc,
        weights,
        donor_codes,
        n_donors,
    )[:, 0]
    purity_centered = weighted_center(
        purity,
        weights,
        donor_codes,
        n_donors,
    )[:, 0]

    base_x = np.column_stack(
        [
            x_centered,
            disease_indicator * x_centered,
        ]
    )
    adjusted_x = np.column_stack(
        [
            x_centered,
            disease_indicator * x_centered,
            purity_centered,
            disease_indicator * purity_centered,
        ]
    )

    base_beta, base_rank = weighted_least_squares(
        base_x,
        y_centered,
        weights,
    )
    adjusted_beta, adjusted_rank = (
        weighted_least_squares(
            adjusted_x,
            y_centered,
            weights,
        )
    )

    observed = np.full((6, n_genes), np.nan)
    if base_rank == base_x.shape[1]:
        observed[0] = base_beta[0]
        observed[1] = (
            base_beta[0] + base_beta[1]
        )
        observed[2] = base_beta[1]
    if adjusted_rank == adjusted_x.shape[1]:
        observed[3] = adjusted_beta[0]
        observed[4] = (
            adjusted_beta[0] + adjusted_beta[1]
        )
        observed[5] = adjusted_beta[1]
    return observed


def weighted_least_squares(
    design_matrix,
    response_matrix,
    weights,
):
    weighted_design = (
        design_matrix * weights[:, None]
    )
    xtwx = design_matrix.T @ weighted_design
    xtwy = design_matrix.T @ (
        response_matrix * weights[:, None]
    )
    beta, _, rank, _ = np.linalg.lstsq(
        xtwx,
        xtwy,
        rcond=None,
    )
    return beta, rank


def empirical_two_sided_p(observed, null_values):
    null_values = np.asarray(
        null_values,
        dtype=float,
    )
    valid = np.isfinite(null_values)
    if not np.isfinite(observed) or valid.sum() == 0:
        return np.nan
    extreme = np.abs(
        null_values[valid]
    ) >= abs(observed)
    return (1.0 + extreme.sum()) / (
        valid.sum() + 1.0
    )


def analyze_group(
    comparison,
    subclass,
    group_table,
    permutations,
    random_generator,
):
    total_nhoods = group_table[
        "result_nhood"
    ].nunique()
    total_donors = group_table[
        "donor_id"
    ].nunique()
    rows_before_filter = len(group_table)

    group_table = group_table[
        group_table["n_cells"]
        >= MIN_CELLS_PER_DONOR_NHOOD
    ].copy()
    donor_row_counts = group_table.groupby(
        "donor_id",
        observed=True,
    ).size()
    valid_donors = donor_row_counts[
        donor_row_counts >= MIN_NHOODS_PER_DONOR
    ].index
    group_table = group_table[
        group_table["donor_id"].isin(valid_donors)
    ].copy()

    diagnosis = {
        "comparison": comparison,
        "subclass": subclass,
        "total_nhoods": total_nhoods,
        "total_donors": total_donors,
        "rows_before_filter": rows_before_filter,
        "rows_after_filter": len(group_table),
        "nhoods_used": group_table[
            "result_nhood"
        ].nunique(),
        "donors_used": group_table[
            "donor_id"
        ].nunique(),
        "reference_donors_used": group_table.loc[
            group_table["condition"] == "Reference",
            "donor_id",
        ].nunique(),
        "disease_donors_used": group_table.loc[
            group_table["condition"] != "Reference",
            "donor_id",
        ].nunique(),
        "min_cells_per_donor_nhood": (
            MIN_CELLS_PER_DONOR_NHOOD
        ),
        "min_nhoods_per_donor": (
            MIN_NHOODS_PER_DONOR
        ),
    }

    if (
        diagnosis["reference_donors_used"]
        < MIN_DONORS_PER_CONDITION
        or diagnosis["disease_donors_used"]
        < MIN_DONORS_PER_CONDITION
        or diagnosis["nhoods_used"] < 5
    ):
        diagnosis["model_fit"] = False
        observed = np.full(
            (len(METRICS), len(FOCUS_GENES)),
            np.nan,
        )
        null_values = np.full(
            (permutations, len(METRICS), len(FOCUS_GENES)),
            np.nan,
        )
        return observed, null_values, diagnosis

    donor_codes, unique_donors = pd.factorize(
        group_table["donor_id"],
        sort=True,
    )
    disease_indicator = (
        group_table["condition"].astype(str)
        != "Reference"
    ).to_numpy(dtype=float)
    weights = group_table[
        "n_cells"
    ].to_numpy(dtype=float)
    neighborhood_log_fc = group_table[
        "logFC"
    ].to_numpy(dtype=float)
    purity = group_table[
        "nhood_purity"
    ].to_numpy(dtype=float)
    expression_values = group_table[
        FOCUS_GENES
    ].to_numpy(dtype=float)
    n_donors = len(unique_donors)

    observed = fit_within_donor_models(
        expression_values,
        neighborhood_log_fc,
        disease_indicator,
        purity,
        weights,
        donor_codes,
        n_donors,
    )

    unique_nhoods = np.asarray(
        sorted(
            group_table["result_nhood"]
            .astype(int)
            .unique()
        ),
        dtype=int,
    )
    nhood_log_fc = (
        group_table[
            ["result_nhood", "logFC"]
        ]
        .drop_duplicates("result_nhood")
        .set_index("result_nhood")
        .reindex(unique_nhoods)["logFC"]
        .to_numpy(dtype=float)
    )
    row_nhood = group_table[
        "result_nhood"
    ].astype(int).to_numpy()

    null_rows = []
    for permutation in range(permutations):
        shuffled_log_fc = random_generator.permutation(
            nhood_log_fc
        )
        shuffle_lookup = dict(
            zip(unique_nhoods, shuffled_log_fc)
        )
        permuted_log_fc = np.asarray(
            [
                shuffle_lookup[nhood]
                for nhood in row_nhood
            ],
            dtype=float,
        )
        permutation_metrics = fit_within_donor_models(
            expression_values,
            permuted_log_fc,
            disease_indicator,
            purity,
            weights,
            donor_codes,
            n_donors,
        )
        null_rows.append(permutation_metrics)

    null_values = np.asarray(
        null_rows,
        dtype=float,
    )
    diagnosis["model_fit"] = True
    return observed, null_values, diagnosis


def build_null_records(
    comparison,
    subclass,
    null_values,
):
    records = []
    for permutation in range(null_values.shape[0]):
        for metric_index, metric in enumerate(
            METRICS
        ):
            for gene_index, gene in enumerate(
                FOCUS_GENES
            ):
                records.append(
                    {
                        "comparison": comparison,
                        "subclass": subclass,
                        "gene": gene,
                        "metric": metric,
                        "permutation": permutation + 1,
                        "value": null_values[
                            permutation,
                            metric_index,
                            gene_index,
                        ],
                    }
                )
    return records


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    cell_expression = load_cell_expression()

    all_edges = []
    for comparison in COMPARISONS:
        all_edges.append(
            read_and_merge_edges(
                comparison,
                cell_expression,
            )
        )
    edges = pd.concat(all_edges, ignore_index=True)
    donor_nhood = build_donor_nhood_table(edges)

    summary_rows = []
    diagnostic_rows = []
    null_records = []
    for comparison_index, comparison in enumerate(
        COMPARISONS
    ):
        comparison_table = donor_nhood[
            donor_nhood["comparison"] == comparison
        ]
        subclasses = sorted(
            comparison_table[
                "nhood_subclass"
            ].unique()
        )
        for subclass_index, subclass in enumerate(
            subclasses
        ):
            group_table = comparison_table[
                comparison_table["nhood_subclass"]
                == subclass
            ].copy()
            analysis_seed = (
                PERMUTATION_SEED
                + comparison_index * 100
                + subclass_index
            )
            random_generator = (
                np.random.default_rng(analysis_seed)
            )
            print(
                f"Analyzing {comparison} / {subclass}: "
                f"{group_table['result_nhood'].nunique()} "
                "nhoods",
                flush=True,
            )
            observed, null_values, diagnosis = (
                analyze_group(
                    comparison,
                    subclass,
                    group_table,
                    PERMUTATIONS,
                    random_generator,
                )
            )
            diagnostic_rows.append(diagnosis)
            null_records.extend(
                build_null_records(
                    comparison,
                    subclass,
                    null_values,
                )
            )
            for gene_index, gene in enumerate(
                FOCUS_GENES
            ):
                summary_row = {
                    **diagnosis,
                    "gene": gene,
                    "permutations": PERMUTATIONS,
                }
                for metric_index, metric in enumerate(
                    METRICS
                ):
                    metric_observed = observed[
                        metric_index,
                        gene_index,
                    ]
                    metric_null = null_values[
                        :,
                        metric_index,
                        gene_index,
                    ]
                    summary_row[
                        f"observed_{metric}"
                    ] = metric_observed
                    summary_row[
                        f"p_{metric}"
                    ] = empirical_two_sided_p(
                        metric_observed,
                        metric_null,
                    )
                    summary_row[
                        f"null_median_{metric}"
                    ] = np.nanmedian(metric_null)
                    summary_row[
                        f"null_lower95_{metric}"
                    ] = np.nanpercentile(
                        metric_null,
                        2.5,
                    )
                    summary_row[
                        f"null_upper95_{metric}"
                    ] = np.nanpercentile(
                        metric_null,
                        97.5,
                    )
                summary_rows.append(summary_row)

    summary = pd.DataFrame(summary_rows)
    diagnostics = pd.DataFrame(diagnostic_rows)
    null_table = pd.DataFrame(null_records)

    summary.to_csv(
        SUMMARY_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    diagnostics.to_csv(
        DIAGNOSTIC_OUTPUT,
        index=False,
        encoding="utf-8-sig",
    )
    null_table.to_csv(
        NULL_OUTPUT,
        index=False,
        compression="gzip",
        encoding="utf-8-sig",
    )

    print("\nSaved:", SUMMARY_OUTPUT)
    print("Saved:", DIAGNOSTIC_OUTPUT)
    print("Saved:", NULL_OUTPUT)
    print("\nPAPPA2 within-donor results:")
    pappa2 = summary[
        summary["gene"] == "PAPPA2"
    ].sort_values(["comparison", "subclass"])
    print(
        pappa2[
            [
                "comparison",
                "subclass",
                "observed_reference_slope_purity_adjusted",
                "p_reference_slope_purity_adjusted",
                "observed_disease_slope_purity_adjusted",
                "p_disease_slope_purity_adjusted",
                "observed_disease_minus_reference_purity_adjusted",
                "p_disease_minus_reference_purity_adjusted",
            ]
        ]
        .round(4)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
