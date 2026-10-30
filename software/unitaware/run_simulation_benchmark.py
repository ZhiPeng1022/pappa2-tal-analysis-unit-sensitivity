"""Run the unit-aware simulation benchmark used by the Genome Medicine draft."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from unitaware import analyze_units, simulate_reversal


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
OUT_DIR = ROOT / "outputs" / "unitaware_simulation"
SEED = 20261031
N_REPLICATES = 100
N_GENES = 8

SCENARIOS = {
    "null": {
        "donor_expression_effect": 0.0,
        "donor_effect_sd": 0.15,
        "within_donor_coupling": 0.0,
        "within_donor_case_addon": 0.0,
        "within_coupling_sd": 0.10,
        "exposure_shift": 0.0,
    },
    "donor_positive_only": {
        "donor_expression_effect": 1.2,
        "donor_effect_sd": 0.8,
        "within_donor_coupling": 0.0,
        "within_donor_case_addon": 0.0,
        "within_coupling_sd": 0.4,
        "exposure_shift": 0.0,
    },
    "within_positive_only": {
        "donor_expression_effect": 0.0,
        "donor_effect_sd": 0.4,
        "within_donor_coupling": -0.9,
        "within_donor_case_addon": 1.8,
        "within_coupling_sd": 0.5,
        "exposure_shift": 0.0,
    },
    "concordant_positive": {
        "donor_expression_effect": 1.2,
        "donor_effect_sd": 0.8,
        "within_donor_coupling": -0.9,
        "within_donor_case_addon": 1.8,
        "within_coupling_sd": 0.5,
        "exposure_shift": 0.0,
    },
    "opposite_donor_positive_within_negative": {
        "donor_expression_effect": 2.0,
        "donor_effect_sd": 0.8,
        "within_donor_coupling": 0.9,
        "within_donor_case_addon": -1.8,
        "within_coupling_sd": 0.5,
        "exposure_shift": 0.0,
    },
    "opposite_donor_negative_within_positive": {
        "donor_expression_effect": -2.0,
        "donor_effect_sd": 0.8,
        "within_donor_coupling": -0.9,
        "within_donor_case_addon": 1.8,
        "within_coupling_sd": 0.5,
        "exposure_shift": 0.0,
    },
}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    seed_rng = np.random.default_rng(SEED)
    rows = []
    for scenario, parameters in SCENARIOS.items():
        for replicate in range(N_REPLICATES):
            seed = int(seed_rng.integers(1, 2**31 - 1))
            expression, cells, neighborhood = simulate_reversal(
                n_reference_donors=25,
                n_case_donors=25,
                cells_per_donor=(30, 80),
                noise_sd=0.8,
                n_genes=N_GENES,
                seed=seed,
                **parameters,
            )
            result = analyze_units(
                expression,
                cells,
                case_levels=["Disease"],
                reference_level="Reference",
                neighborhood_logfc=neighborhood,
            )
            table = result.reversals
            rows.append(
                {
                    "scenario": scenario,
                    "replicate": replicate + 1,
                    "seed": seed,
                    "mean_donor_effect": table["donor_effect"].mean(),
                    "mean_within_effect": table[
                        "within_donor_effect"
                    ].mean(),
                    "mean_neighborhood_effect": table[
                        "neighborhood_effect"
                    ].mean(),
                    "donor_within_reversal_fraction": table[
                        "donor_within_reversal"
                    ].mean(),
                    "donor_neighborhood_reversal_fraction": table[
                        "donor_neighborhood_reversal"
                    ].mean(),
                    "neighborhood_within_reversal_fraction": table[
                        "neighborhood_within_reversal"
                    ].mean(),
                }
            )
        print(f"completed {scenario}", flush=True)

    runs = pd.DataFrame(rows)
    runs.to_csv(OUT_DIR / "simulation_runs.csv", index=False)
    summary = (
        runs.groupby("scenario", observed=True)
        .agg(
            mean_donor_effect=("mean_donor_effect", "mean"),
            mean_within_effect=("mean_within_effect", "mean"),
            mean_neighborhood_effect=(
                "mean_neighborhood_effect",
                "mean",
            ),
            donor_within_reversal=(
                "donor_within_reversal_fraction",
                "mean",
            ),
            donor_neighborhood_reversal=(
                "donor_neighborhood_reversal_fraction",
                "mean",
            ),
            neighborhood_within_reversal=(
                "neighborhood_within_reversal_fraction",
                "mean",
            ),
        )
        .reindex(SCENARIOS)
        .reset_index()
    )
    summary.to_csv(OUT_DIR / "simulation_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
