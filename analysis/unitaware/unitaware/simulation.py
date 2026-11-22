"""Simulation helpers for unit-aware sensitivity analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd


def simulate_reversal(
    *,
    n_reference_donors: int = 40,
    n_case_donors: int = 40,
    cells_per_donor: tuple[int, int] = (40, 100),
    donor_expression_effect: float = 0.0,
    donor_effect_sd: float = 0.0,
    within_donor_coupling: float = 0.0,
    within_donor_case_addon: float = 0.0,
    within_coupling_sd: float = 0.0,
    exposure_shift: float = 0.0,
    noise_sd: float = 0.8,
    n_genes: int = 1,
    seed: int = 20261028,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    cell_rows = []
    expression_rows = []
    gene_names = [f"GENE_{index:03d}" for index in range(n_genes)]
    gene_donor_effect = rng.normal(
        donor_expression_effect,
        donor_effect_sd,
        size=n_genes,
    )
    gene_reference_coupling = rng.normal(
        within_donor_coupling,
        within_coupling_sd,
        size=n_genes,
    )
    gene_case_addon = rng.normal(
        within_donor_case_addon,
        within_coupling_sd,
        size=n_genes,
    )
    neighborhoods = []
    donor_index = 0
    for condition, count in [
        ("Reference", n_reference_donors),
        ("Disease", n_case_donors),
    ]:
        case = int(condition == "Disease")
        for _ in range(count):
            donor_id = f"{condition}_{donor_index}"
            donor_index += 1
            n_cells = int(
                rng.integers(
                    cells_per_donor[0],
                    cells_per_donor[1] + 1,
                )
            )
            donor_intercept = rng.normal(0, 0.5, size=n_genes)
            exposure = rng.normal(
                exposure_shift * case,
                1.0,
                size=n_cells,
            )
            expression = (
                donor_intercept[None, :]
                + gene_donor_effect[None, :] * case
                + (
                    gene_reference_coupling[None, :]
                    + gene_case_addon[None, :] * case
                )
                * exposure[:, None]
                + rng.normal(
                    0,
                    noise_sd,
                    size=(n_cells, n_genes),
                )
            )
            for cell_index in range(n_cells):
                cell_id = f"{donor_id}_{cell_index}"
                cell_rows.append(
                    {
                        "cell_id": cell_id,
                        "donor_id": donor_id,
                        "condition": condition,
                        "state": "target",
                        "exposure": exposure[cell_index],
                    }
                )
                expression_rows.append(
                    {
                        "cell_id": cell_id,
                        **{
                            gene: expression[cell_index, gene_index]
                            for gene_index, gene in enumerate(gene_names)
                        },
                    }
                )

    cells = pd.DataFrame(cell_rows).set_index("cell_id")
    expression = pd.DataFrame(expression_rows).set_index("cell_id")

    n_neighborhoods = 50
    exposure_quantiles = pd.qcut(
        cells["exposure"],
        q=n_neighborhoods,
        labels=False,
        duplicates="drop",
    )
    cells["neighborhood_id"] = exposure_quantiles.astype(int)
    disease_fraction = (
        cells["condition"].eq("Disease")
        .groupby(cells["neighborhood_id"], observed=True)
        .mean()
    )
    neighborhood_table = pd.DataFrame(
        {
            "neighborhood_id": sorted(cells["neighborhood_id"].unique()),
            "state": "target",
            "comparison": "disease_vs_reference",
            "logfc": disease_fraction.reindex(
                sorted(cells["neighborhood_id"].unique())
            ).to_numpy(),
        }
    )
    return expression, cells, neighborhood_table
