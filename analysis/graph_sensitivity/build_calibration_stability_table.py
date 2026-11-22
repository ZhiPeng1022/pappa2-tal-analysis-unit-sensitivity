"""Build the formal calibration-stability audit table."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
GRAPH = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1"
PANEL = GRAPH / "panel_calibration_by_graph.csv"
CELL_RANK = GRAPH / "aTAL2_PAPPA2_by_graph.csv"
OUTPUT = GRAPH / "calibration_stability_table.csv"
OUTPUT_LONG = GRAPH / "calibration_stability_runs.csv"
OUTPUT_MD = GRAPH / "calibration_stability_report.md"

sys.path.insert(0, str(ROOT / "work" / "unitaware_framework"))
from unitaware import summarize_calibration_stability  # noqa: E402


def main() -> None:
    panel = pd.read_csv(PANEL)
    cell_rank = pd.read_csv(CELL_RANK)

    panel_runs = []
    for row in panel.itertuples(index=False):
        panel_runs.append(
            {
                "endpoint": "expression_matched_bias",
                "variant": row.variant,
                "effect": (
                    row.target_reversal_fraction
                    - row.expression_matched_mean
                ),
                "p_value": row.expression_matched_p,
            }
        )
        panel_runs.append(
            {
                "endpoint": "donor_conditioned_panel_residual",
                "variant": row.variant,
                "effect": row.target_minus_effect_aware,
                "p_value": row.effect_aware_p,
            }
        )
    cell_runs = [
        {
            "endpoint": "aTAL2_PAPPA2_cell_rank_conflict",
            "variant": row.variant,
            "effect": row.disease_minus_reference,
            "p_value": row.label_permutation_p,
        }
        for row in cell_rank.itertuples(index=False)
    ]
    runs = pd.DataFrame(panel_runs + cell_runs)
    runs.to_csv(OUTPUT_LONG, index=False)

    stability = summarize_calibration_stability(
        runs,
        group_cols=("endpoint",),
        effect_col="effect",
        p_col="p_value",
        alpha=0.05,
        null_tolerance=0.0,
    )
    stability.to_csv(OUTPUT, index=False)

    lines = [
        "# Calibration stability audit",
        "",
        "Stability is descriptive: it reports direction and decision "
        "consistency across the four actual-Milo settings. It does not create "
        "a new combined p value.",
        "",
        stability.to_string(index=False),
        "",
        "## Interpretation",
        "",
        "- Expression matching has a stable positive excess, showing that the "
        "need for effect-aware calibration is graph- and seed-robust.",
        "- The donor-conditioned panel residual changes sign and significance "
        "across settings, so a universal calibrated-null conclusion is not "
        "supported.",
        "- The aTAL2/PAPPA2 cell-rank conflict remains negative and significant "
        "across settings, providing a stable case-level benchmark.",
        "",
    ]
    OUTPUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(stability.to_string(index=False))


if __name__ == "__main__":
    main()
