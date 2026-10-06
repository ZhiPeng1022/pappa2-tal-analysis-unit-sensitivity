from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from analyze_full_tal_within_donor_nhood_effect import (
    COMPARISONS,
    FOCUS_GENES,
    OUTPUT_ROOT,
    load_cell_expression,
    read_and_merge_edges,
)


MIN_CELLS_PER_DONOR = 5
PERMUTATIONS = 1000
PERMUTATION_SEED = 20260927

SUMMARY_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_within_donor_cell_rank.csv"
)
NULL_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_within_donor_cell_rank_null.csv.gz"
)
DIAGNOSTIC_OUTPUT = (
    OUTPUT_ROOT
    / "KPMP_full_TAL_within_donor_cell_rank_diagnostics.csv"
)


def build_cell_table(edges):
    cell_data = (
        edges.groupby("cell_id", observed=True)
        .agg(
            donor_id=("donor_id", "first"),
            condition=("condition", "first"),
            exposure=("logFC", "mean"),
            mean_purity=("nhood_purity", "mean"),
            nhood_memberships=("result_nhood", "nunique"),
            **{
                gene: (gene, "first")
                for gene in FOCUS_GENES
            },
        )
        .reset_index()
    )
    return cell_data


def rank_correlation_with_exposure(
    exposure_rank,
    expression_ranks,
):
    x_centered = (
        exposure_rank - exposure_rank.mean()
    )
    y_centered = (
        expression_ranks
        - expression_ranks.mean(axis=0)
    )
    numerator = x_centered @ y_centered
    denominator = np.sqrt(
        np.sum(x_centered**2)
        * np.sum(y_centered**2, axis=0)
    )
    return np.divide(
        numerator,
        denominator,
        out=np.full(
            expression_ranks.shape[1],
            np.nan,
        ),
        where=denominator > 0,
    )


def prepare_donor_blocks(cell_data):
    blocks = []
    for donor_id, donor_data in cell_data.groupby(
        "donor_id",
        observed=True,
    ):
        donor_data = donor_data.dropna(
            subset=["exposure", *FOCUS_GENES]
        )
        if len(donor_data) < MIN_CELLS_PER_DONOR:
            continue
        if donor_data["exposure"].nunique() < 2:
            continue

        exposure_rank = rankdata(
            donor_data["exposure"].to_numpy(dtype=float)
        )
        expression_ranks = np.column_stack(
            [
                rankdata(
                    donor_data[gene].to_numpy(
                        dtype=float
                    )
                )
                for gene in FOCUS_GENES
            ]
        ).astype(float)
        condition = str(
            donor_data["condition"].iloc[0]
        )
        if (
            donor_data["condition"].astype(str).nunique()
            > 1
        ):
            raise ValueError(
                f"{donor_id} has multiple conditions."
            )

        blocks.append(
            {
                "donor_id": str(donor_id),
                "condition": condition,
                "n_cells": len(donor_data),
                "exposure_rank": exposure_rank,
                "expression_ranks": expression_ranks,
                "observed_rho": (
                    rank_correlation_with_exposure(
                        exposure_rank,
                        expression_ranks,
                    )
                ),
            }
        )
    return blocks


def median_or_nan(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan
    return np.median(values)


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
    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    valid_result = np.empty_like(adjusted)
    valid_result[order] = adjusted
    result[valid] = valid_result
    return result


def analyze_group(
    comparison,
    subclass,
    cell_data,
    permutations,
    random_generator,
):
    blocks = prepare_donor_blocks(cell_data)
    reference_blocks = [
        block for block in blocks
        if block["condition"] == "Reference"
    ]
    disease_blocks = [
        block for block in blocks
        if block["condition"] != "Reference"
    ]

    diagnosis = {
        "comparison": comparison,
        "subclass": subclass,
        "cells_total": len(cell_data),
        "donors_total": cell_data[
            "donor_id"
        ].nunique(),
        "eligible_donors": len(blocks),
        "eligible_reference_donors": len(
            reference_blocks
        ),
        "eligible_disease_donors": len(
            disease_blocks
        ),
        "median_cells_per_eligible_donor": (
            np.median(
                [
                    block["n_cells"]
                    for block in blocks
                ]
            )
            if blocks
            else np.nan
        ),
        "min_cells_per_donor": MIN_CELLS_PER_DONOR,
    }

    if len(blocks) == 0:
        observed = np.full(
            (4, len(FOCUS_GENES)),
            np.nan,
        )
        null_values = np.full(
            (permutations, 3, len(FOCUS_GENES)),
            np.nan,
        )
        difference_null = np.full(
            (permutations, len(FOCUS_GENES)),
            np.nan,
        )
        diagnosis["model_fit"] = False
        return (
            observed,
            null_values,
            difference_null,
            np.full(len(FOCUS_GENES), np.nan),
            diagnosis,
        )

    reference_rhos = np.asarray(
        [
            block["observed_rho"]
            for block in reference_blocks
        ],
        dtype=float,
    )
    disease_rhos = np.asarray(
        [
            block["observed_rho"]
            for block in disease_blocks
        ],
        dtype=float,
    )
    all_rhos = np.asarray(
        [
            block["observed_rho"]
            for block in blocks
        ],
        dtype=float,
    )
    observed = np.vstack(
        [
            np.nanmedian(all_rhos, axis=0),
            (
                np.nanmedian(reference_rhos, axis=0)
                if len(reference_rhos) > 0
                else np.full(
                    len(FOCUS_GENES),
                    np.nan,
                )
            ),
            (
                np.nanmedian(disease_rhos, axis=0)
                if len(disease_rhos) > 0
                else np.full(
                    len(FOCUS_GENES),
                    np.nan,
                )
            ),
            (
                np.nanmedian(disease_rhos, axis=0)
                - np.nanmedian(reference_rhos, axis=0)
                if (
                    len(disease_rhos) > 0
                    and len(reference_rhos) > 0
                )
                else np.full(
                    len(FOCUS_GENES),
                    np.nan,
                )
            ),
        ]
    )

    null_values = np.full(
        (permutations, 3, len(FOCUS_GENES)),
        np.nan,
    )
    for permutation in range(permutations):
        permuted_rhos = []
        for block in blocks:
            permuted_indices = (
                random_generator.permutation(
                    block["n_cells"]
                )
            )
            permuted_rhos.append(
                rank_correlation_with_exposure(
                    block["exposure_rank"],
                    block[
                        "expression_ranks"
                    ][permuted_indices],
                )
            )
        permuted_rhos = np.asarray(
            permuted_rhos,
            dtype=float,
        )
        reference_count = len(reference_blocks)
        null_values[
            permutation,
            0,
        ] = np.nanmedian(permuted_rhos, axis=0)
        if reference_count > 0:
            null_values[
                permutation,
                1,
            ] = np.nanmedian(
                permuted_rhos[:reference_count, :],
                axis=0,
            )
        if len(disease_blocks) > 0:
            null_values[
                permutation,
                2,
            ] = np.nanmedian(
                permuted_rhos[reference_count:, :],
                axis=0,
            )

    difference_null = np.full(
        (permutations, len(FOCUS_GENES)),
        np.nan,
    )
    all_conditions = np.asarray(
        [
            block["condition"] == "Reference"
            for block in blocks
        ],
        dtype=bool,
    )
    if (
        len(reference_blocks) > 1
        and len(disease_blocks) > 1
    ):
        for permutation in range(permutations):
            shuffled_labels = (
                random_generator.permutation(
                    all_conditions
                )
            )
            shuffled_reference = all_rhos[
                shuffled_labels
            ]
            shuffled_disease = all_rhos[
                ~shuffled_labels
            ]
            difference_null[
                permutation
            ] = (
                np.nanmedian(shuffled_disease, axis=0)
                - np.nanmedian(
                    shuffled_reference,
                    axis=0,
                )
            )

    p_values = {
        "p_all": [
            empirical_two_sided_p(
                observed[0, gene_index],
                null_values[:, 0, gene_index],
            )
            for gene_index in range(len(FOCUS_GENES))
        ],
        "p_reference": [
            empirical_two_sided_p(
                observed[1, gene_index],
                null_values[:, 1, gene_index],
            )
            for gene_index in range(len(FOCUS_GENES))
        ],
        "p_disease": [
            empirical_two_sided_p(
                observed[2, gene_index],
                null_values[:, 2, gene_index],
            )
            for gene_index in range(len(FOCUS_GENES))
        ],
        "p_disease_minus_reference": [
            empirical_two_sided_p(
                observed[3, gene_index],
                difference_null[:, gene_index],
            )
            for gene_index in range(len(FOCUS_GENES))
        ],
    }
    diagnosis["model_fit"] = bool(
        len(reference_blocks) >= 3
        and len(disease_blocks) >= 3
    )
    return (
        observed,
        null_values,
        difference_null,
        p_values,
        diagnosis,
    )


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    cell_expression = load_cell_expression()
    edges = pd.concat(
        [
            read_and_merge_edges(
                comparison,
                cell_expression,
            )
            for comparison in COMPARISONS
        ],
        ignore_index=True,
    )

    summary_rows = []
    diagnostic_rows = []
    null_records = []
    for comparison_index, comparison in enumerate(
        COMPARISONS
    ):
        comparison_edges = edges[
            edges["comparison"] == comparison
        ]
        subclasses = sorted(
            comparison_edges[
                "nhood_subclass"
            ].unique()
        )
        for subclass_index, subclass in enumerate(
            subclasses
        ):
            group_edges = comparison_edges[
                comparison_edges["nhood_subclass"]
                == subclass
            ]
            cell_data = build_cell_table(group_edges)
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
                f"{len(cell_data)} cells",
                flush=True,
            )
            (
                observed,
                null_values,
                difference_null,
                p_values,
                diagnosis,
            ) = analyze_group(
                comparison,
                subclass,
                cell_data,
                PERMUTATIONS,
                random_generator,
            )
            diagnostic_rows.append(diagnosis)

            blocks = prepare_donor_blocks(cell_data)
            reference_rhos = np.asarray(
                [
                    block["observed_rho"]
                    for block in blocks
                    if block["condition"] == "Reference"
                ],
                dtype=float,
            )
            disease_rhos = np.asarray(
                [
                    block["observed_rho"]
                    for block in blocks
                    if block["condition"] != "Reference"
                ],
                dtype=float,
            )
            for gene_index, gene in enumerate(
                FOCUS_GENES
            ):
                reference_gene = (
                    reference_rhos[:, gene_index]
                    if len(reference_rhos) > 0
                    else np.asarray([])
                )
                disease_gene = (
                    disease_rhos[:, gene_index]
                    if len(disease_rhos) > 0
                    else np.asarray([])
                )
                all_gene = np.concatenate(
                    [
                        reference_gene[
                            np.isfinite(
                                reference_gene
                            )
                        ],
                        disease_gene[
                            np.isfinite(
                                disease_gene
                            )
                        ],
                    ]
                )
                summary_rows.append(
                    {
                        **diagnosis,
                        "gene": gene,
                        "permutations": PERMUTATIONS,
                        "observed_median_rho_all": (
                            observed[0, gene_index]
                        ),
                        "p_all": p_values["p_all"][
                            gene_index
                        ],
                        "observed_median_rho_reference": (
                            observed[1, gene_index]
                        ),
                        "p_reference": (
                            p_values["p_reference"][
                                gene_index
                            ]
                        ),
                        "observed_median_rho_disease": (
                            observed[2, gene_index]
                        ),
                        "p_disease": (
                            p_values["p_disease"][
                                gene_index
                            ]
                        ),
                        "observed_disease_minus_reference": (
                            observed[3, gene_index]
                        ),
                        "p_disease_minus_reference": (
                            p_values[
                                "p_disease_minus_reference"
                            ][gene_index]
                        ),
                        "positive_all": int(
                            (all_gene > 0).sum()
                        ),
                        "negative_all": int(
                            (all_gene < 0).sum()
                        ),
                        "positive_reference": int(
                            (
                                reference_gene > 0
                            ).sum()
                        ),
                        "negative_reference": int(
                            (
                                reference_gene < 0
                            ).sum()
                        ),
                        "positive_disease": int(
                            (
                                disease_gene > 0
                            ).sum()
                        ),
                        "negative_disease": int(
                            (
                                disease_gene < 0
                            ).sum()
                        ),
                    }
                )

            for permutation in range(
                PERMUTATIONS
            ):
                for metric_index, metric in enumerate(
                    [
                        "all",
                        "reference",
                        "disease",
                    ]
                ):
                    for gene_index, gene in enumerate(
                        FOCUS_GENES
                    ):
                        null_records.append(
                            {
                                "comparison": comparison,
                                "subclass": subclass,
                                "gene": gene,
                                "metric": metric,
                                "permutation": (
                                    permutation + 1
                                ),
                                "value": null_values[
                                    permutation,
                                    metric_index,
                                    gene_index,
                                ],
                            }
                        )
                for gene_index, gene in enumerate(
                    FOCUS_GENES
                ):
                    null_records.append(
                        {
                            "comparison": comparison,
                            "subclass": subclass,
                            "gene": gene,
                            "metric": (
                                "disease_minus_reference"
                            ),
                            "permutation": (
                                permutation + 1
                            ),
                            "value": difference_null[
                                permutation,
                                gene_index,
                            ],
                        }
                    )

    summary = pd.DataFrame(summary_rows)
    for p_column in [
        "p_all",
        "p_reference",
        "p_disease",
        "p_disease_minus_reference",
    ]:
        summary[f"FDR_{p_column}"] = (
            benjamini_hochberg(
                summary[p_column].to_numpy()
            )
        )
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
    print("\nPAPPA2 within-donor cell-rank results:")
    pappa2 = summary[
        summary["gene"] == "PAPPA2"
    ].sort_values(["comparison", "subclass"])
    print(
        pappa2[
            [
                "comparison",
                "subclass",
                "model_fit",
                "eligible_reference_donors",
                "eligible_disease_donors",
                "observed_median_rho_all",
                "p_all",
                "FDR_p_all",
                "observed_median_rho_reference",
                "p_reference",
                "FDR_p_reference",
                "observed_median_rho_disease",
                "p_disease",
                "FDR_p_disease",
                "observed_disease_minus_reference",
                "p_disease_minus_reference",
                "FDR_p_disease_minus_reference",
            ]
        ]
        .round(4)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
