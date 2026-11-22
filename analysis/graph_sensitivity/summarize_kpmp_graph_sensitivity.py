"""Summarize actual Milo graph sensitivity for the KPMP benchmark."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
SOURCE = ROOT / "outputs" / "kpmp_graph_sensitivity" / "v1"
PANEL_DIR = SOURCE / "milo_panel"
CELL_RANK_DIR = SOURCE / "cell_rank"
COMPOSITION = (
    ROOT
    / "outputs"
    / "composition_expression_decomposition"
    / "composition_expression_summary.csv"
)


def main() -> None:
    panel_rows = []
    for variant in ["k15", "k30", "k30_seed", "k50"]:
        path = PANEL_DIR / variant / "summary.json"
        if not path.exists():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        endpoint = payload["donor_within"]
        panel_rows.append(
            {
                "variant": variant,
                "target_reversal_fraction": endpoint["target_400_mean"],
                "expression_matched_mean": endpoint["random_400_mean"],
                "expression_matched_p": endpoint["empirical_p_400"],
                "effect_aware_null_mean": payload[
                    "effect_aware_null_mean"
                ],
                "effect_aware_p": payload[
                    "effect_aware_empirical_p"
                ],
                "target_minus_effect_aware": (
                    endpoint["target_400_mean"]
                    - payload["effect_aware_null_mean"]
                ),
            }
        )
    panel = pd.DataFrame(panel_rows)
    panel.to_csv(
        SOURCE / "panel_calibration_by_graph.csv",
        index=False,
    )

    cell_rank = pd.read_csv(
        CELL_RANK_DIR / "aTAL2_PAPPA2_graph_sensitivity.csv"
    )
    cell_rank.to_csv(
        SOURCE / "aTAL2_PAPPA2_by_graph.csv",
        index=False,
    )

    composition = pd.read_csv(COMPOSITION)
    composition = composition[
        composition["gene"].eq("PAPPA2")
        & composition["comparison"].eq("all_disease")
    ].copy()
    composition.to_csv(
        SOURCE / "PAPPA2_composition_boundary.csv",
        index=False,
    )

    report = [
        "# KPMP graph and seed sensitivity",
        "",
        "## Important distinction",
        "",
        "The formal sensitivity analysis uses actual Milo runs. The files "
        "`kpmp_graph_k15_*` and `kpmp_graph_k50_*` in the parent directory "
        "are exploratory fixed-kNN approximations and are not used as formal "
        "Milo results.",
        "",
        "## Frozen 400-gene donor-within calibration",
        "",
        panel.to_string(index=False),
        "",
        "Expression matching consistently underestimated the reversal "
        "fraction. The residual after donor-conditioned effect-aware "
        "calibration, however, was not stable across graph or seed choices: "
        "the original k=30 Milo graph gave no excess, whereas k=15, k=50, "
        "and a new k=30 seed retained positive residuals. The graph-robust "
        "claim is therefore that effect matching is an inadequate null; the "
        "claim that the calibrated kidney excess is universally absent is "
        "not supported.",
        "",
        "## aTAL2/PAPPA2 donor-within cell-rank conflict",
        "",
        cell_rank.to_string(index=False),
        "",
        "The disease-minus-reference within-donor effect remained negative "
        "in every graph variant. The donor-level PAPPA2 difference is graph-"
        "invariant because it is a pseudobulk quantity; the composition-"
        "expression decomposition is likewise unchanged by the neighbor "
        "graph because it is computed from donor subtype proportions and "
        "subtype means.",
        "",
        "## Composition boundary",
        "",
        composition[
            [
                "comparison",
                "gene",
                "observed_difference",
                "composition_effect",
                "expression_effect",
                "interaction_and_covariance",
            ]
        ].to_string(index=False),
        "",
        "## Manuscript implication",
        "",
        "- Retain the aTAL2/PAPPA2 cell-rank conflict as a graph- and seed-"
        "robust observation.",
        "- Retain the expression-matched-versus-effect-aware contrast as the "
        "main methodological finding.",
        "- Do not state that effect-aware calibration abolished the kidney "
        "panel excess as a universal result. Report the graph/seed boundary "
        "explicitly and describe the conclusion as calibration-dependent.",
        "- Keep the original k=30 fixed-seed analysis as the prespecified "
        "primary result, with the graph/seed sensitivity as a robustness "
        "boundary.",
        "",
    ]
    (SOURCE / "milo_graph_sensitivity_report.md").write_text(
        "\n".join(report),
        encoding="utf-8",
    )
    print(panel.to_string(index=False))
    print(cell_rank.to_string(index=False))


if __name__ == "__main__":
    main()
