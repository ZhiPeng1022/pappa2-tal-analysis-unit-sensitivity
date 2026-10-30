import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from unitaware import (
    analyze_units,
    decompose_composition,
    paired_reversal_permutation,
    reversal_fraction,
    simulate_reversal,
)
from unitaware.cli import main


class UnitAwareCoreTests(unittest.TestCase):
    def test_positive_donor_negative_within_is_detected(self):
        expression, cells, neighborhood = simulate_reversal(
            donor_expression_effect=1.2,
            within_donor_coupling=0.9,
            within_donor_case_addon=-1.8,
            noise_sd=0.35,
            n_genes=20,
            n_reference_donors=30,
            n_case_donors=30,
            seed=7,
        )
        result = analyze_units(
            expression,
            cells,
            case_levels=["Disease"],
            reference_level="Reference",
            neighborhood_logfc=neighborhood,
        )
        self.assertTrue((result.reversals["donor_effect"] > 0).all())
        self.assertLess(
            result.reversals["within_donor_effect"].mean(),
            -0.3,
        )
        self.assertGreaterEqual(
            result.reversals["donor_within_reversal"].mean(),
            0.1,
        )

    def test_concordant_positive_effects_are_not_reversals(self):
        effects = pd.DataFrame(
            {
                "left": [1.0, 2.0, -0.2, np.nan],
                "right": [0.8, 1.2, 0.9, 2.0],
            }
        )
        self.assertEqual(
            reversal_fraction(
                effects,
                left_col="left",
                right_col="right",
                z_threshold=0.5,
            ),
            0.0,
        )
        self.assertTrue(
            np.isnan(
                reversal_fraction(
                    pd.DataFrame({"left": [np.nan], "right": [1.0]}),
                    left_col="left",
                    right_col="right",
                )
            )
        )

    def test_paired_permutation_detects_coupled_opposite_effects(self):
        effects = pd.DataFrame(
            {
                "gene": [f"G{index}" for index in range(8)],
                "state": ["TAL"] * 8,
                "comparison": ["Disease_vs_Reference"] * 8,
                "left": np.linspace(-2.0, 2.0, 8),
                "right": -np.linspace(-2.0, 2.0, 8),
            }
        )
        null_runs, summary = paired_reversal_permutation(
            effects,
            left_col="left",
            right_col="right",
            target_genes=set(effects["gene"]),
            strata_cols=("state", "comparison"),
            n_permutations=199,
            seed=11,
        )
        self.assertEqual(len(null_runs), 199)
        self.assertEqual(
            summary["observed_reversal_fraction"],
            0.75,
        )
        self.assertLess(summary["empirical_p"], 0.05)
        self.assertLess(summary["null_mean"], 1.0)

    def test_composition_expression_decomposition(self):
        cells = pd.DataFrame(
            {
                "donor_id": ["R", "R", "R", "R", "D", "D", "D", "D"],
                "condition": ["Reference"] * 4 + ["Disease"] * 4,
                "state": ["A", "A", "A", "B", "A", "B", "B", "B"],
            },
            index=[f"C{index}" for index in range(8)],
        )
        expression = pd.DataFrame(
            {"GENE": [2.0, 2.0, 2.0, 0.0, 4.0, 0.0, 0.0, 0.0]},
            index=cells.index,
        )
        result = decompose_composition(
            expression,
            cells,
            gene="GENE",
            case_levels=["Disease"],
            reference_level="Reference",
        ).iloc[0]
        self.assertAlmostEqual(result["observed_difference"], -0.5)
        self.assertAlmostEqual(result["composition_effect"], -1.0)
        self.assertAlmostEqual(result["expression_effect"], 1.5)
        self.assertAlmostEqual(
            result["interaction_and_covariance"],
            -1.0,
        )

    def test_simulation_without_effect_is_near_zero_average(self):
        expression, cells, neighborhood = simulate_reversal(
            donor_expression_effect=0.0,
            donor_effect_sd=0.0,
            within_donor_coupling=0.0,
            within_donor_case_addon=0.0,
            within_coupling_sd=0.0,
            noise_sd=0.1,
            n_genes=3,
            n_reference_donors=12,
            n_case_donors=12,
            seed=19,
        )
        result = analyze_units(
            expression,
            cells,
            case_levels=["Disease"],
            reference_level="Reference",
            neighborhood_logfc=neighborhood,
        )
        self.assertLess(
            abs(result.reversals["donor_effect"].mean()),
            0.2,
        )
        self.assertLess(
            abs(result.reversals["within_donor_effect"].mean()),
            0.2,
        )

    def test_cli_simulate_writes_expected_tables(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "simulated"
            argv = [
                "unitaware",
                "simulate",
                "--seed",
                "5",
                "--out",
                str(output),
            ]
            with patch.object(sys, "argv", argv):
                self.assertEqual(main(), 0)
            self.assertTrue((output / "expression.csv").is_file())
            self.assertTrue((output / "cells.csv").is_file())
            self.assertTrue(
                (output / "neighborhood_logfc.csv").is_file()
            )
            with (output / "cells.csv").open(
                newline="",
                encoding="utf-8",
            ) as handle:
                reader = csv.DictReader(handle)
                self.assertIn("donor_id", reader.fieldnames)


if __name__ == "__main__":
    unittest.main()
