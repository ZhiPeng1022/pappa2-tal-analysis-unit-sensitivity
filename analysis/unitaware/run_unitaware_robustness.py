"""Run effect-aware robustness analyses for the Genome Medicine revision."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
FRAMEWORK = ROOT / "work" / "unitaware_framework"
sys.path.insert(0, str(FRAMEWORK))

from unitaware import paired_reversal_permutation  # noqa: E402


PANEL_DIR = ROOT / "outputs" / "prespecified_unit_panel"
RESULTS_DIR = ROOT / "outputs" / "prespecified_unit_panel_results"
OUT_DIR = ROOT / "outputs" / "unitaware_robustness"
METRICS = RESULTS_DIR / "gene_unit_metrics_long.csv"
PANEL = PANEL_DIR / "target_panel_400.csv"
CONTROLS = PANEL_DIR / "negative_control_pool_2000.csv"

THRESHOLDS = (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
N_RANDOM_SETS = 1000
N_PERMUTATIONS = 10000
SEED = 20261028

PAIRS = {
    "donor_within": ("donor_effect_z", "within_delta_rho_z"),
    "donor_neighborhood": (
        "donor_effect_z",
        "neighborhood_rho_z",
    ),
    "neighborhood_within": (
        "neighborhood_rho_z",
        "within_delta_rho_z",
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def empirical_p(observed: float, null: np.ndarray) -> float:
    values = np.asarray(null, dtype=float)
    values = values[np.isfinite(values)]
    if not np.isfinite(observed) or len(values) == 0:
        return np.nan
    return float(
        (1.0 + np.sum(values >= observed)) / (len(values) + 1.0)
    )


def gene_reversal_fractions(
    metrics: pd.DataFrame,
    *,
    left_col: str,
    right_col: str,
    threshold: float,
) -> pd.Series:
    valid = metrics[[left_col, right_col]].notna().all(axis=1)
    frame = metrics.loc[valid, ["gene", left_col, right_col]].copy()
    frame["reversal"] = (
        (
            np.sign(frame[left_col])
            * np.sign(frame[right_col])
            < 0
        )
        & (frame[left_col].abs() >= threshold)
        & (frame[right_col].abs() >= threshold)
    )
    return frame.groupby("gene", observed=True)["reversal"].mean()


def draw_random_sets(
    panel: pd.DataFrame,
    controls: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(SEED)
    control_genes = controls["gene"].to_numpy(dtype=object)
    target_expression = panel["mean_lognorm"].to_numpy(dtype=float)
    breaks = np.quantile(
        target_expression,
        np.linspace(0, 1, 11),
    )
    target_bins = np.clip(
        np.digitize(target_expression, breaks[1:-1]),
        0,
        9,
    )
    target_bin_counts = np.bincount(target_bins, minlength=10)
    control_bins = np.clip(
        np.digitize(
            controls["mean_lognorm"].to_numpy(dtype=float),
            breaks[1:-1],
        ),
        0,
        9,
    )
    members = {
        index: control_genes[control_bins == index]
        for index in range(10)
    }

    sample_400 = []
    sample_300 = []
    sample_stratified = []
    for _ in range(N_RANDOM_SETS):
        sample_400.append(
            rng.choice(control_genes, size=400, replace=False)
        )
        sample_300.append(
            rng.choice(control_genes, size=300, replace=False)
        )
        sample_stratified.append(
            np.concatenate(
                [
                    rng.choice(
                        members[index],
                        size=int(target_bin_counts[index]),
                        replace=False,
                    )
                    for index in range(10)
                ]
            )
        )
    return (
        np.vstack(sample_400),
        np.vstack(sample_300),
        np.vstack(sample_stratified),
    )


def threshold_sensitivity(
    metrics: pd.DataFrame,
    panel: pd.DataFrame,
    controls: pd.DataFrame,
) -> pd.DataFrame:
    target_400 = panel["gene"].tolist()
    target_300 = panel.loc[
        ~panel["gene_set"].eq("internal_negative_control"),
        "gene",
    ].tolist()
    sample_400, sample_300, sample_stratified = draw_random_sets(
        panel,
        controls,
    )
    rows = []
    for threshold in THRESHOLDS:
        for pair, (left_col, right_col) in PAIRS.items():
            fractions = gene_reversal_fractions(
                metrics,
                left_col=left_col,
                right_col=right_col,
                threshold=threshold,
            )
            observed_400 = float(
                fractions.reindex(target_400).mean()
            )
            observed_300 = float(
                fractions.reindex(target_300).mean()
            )
            random_400 = np.nanmean(
                fractions.reindex(
                    sample_400.ravel()
                ).to_numpy().reshape(sample_400.shape),
                axis=1,
            )
            random_300 = np.nanmean(
                fractions.reindex(
                    sample_300.ravel()
                ).to_numpy().reshape(sample_300.shape),
                axis=1,
            )
            random_stratified = np.nanmean(
                fractions.reindex(
                    sample_stratified.ravel()
                ).to_numpy().reshape(sample_stratified.shape),
                axis=1,
            )
            rows.append(
                {
                    "threshold": threshold,
                    "pair": pair,
                    "target_400": observed_400,
                    "random_400": float(random_400.mean()),
                    "empirical_p_400": empirical_p(
                        observed_400,
                        random_400,
                    ),
                    "target_300": observed_300,
                    "random_300": float(random_300.mean()),
                    "empirical_p_300": empirical_p(
                        observed_300,
                        random_300,
                    ),
                    "random_stratified_400": float(
                        random_stratified.mean()
                    ),
                    "empirical_p_stratified_400": empirical_p(
                        observed_400,
                        random_stratified,
                    ),
                }
            )
    return pd.DataFrame(rows)


def pairing_permutation(
    metrics: pd.DataFrame,
    panel: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    target_genes = set(panel["gene"])
    summary_rows = []
    null_tables = []
    for pair, (left_col, right_col) in PAIRS.items():
        null_table, summary = paired_reversal_permutation(
            metrics,
            left_col=left_col,
            right_col=right_col,
            target_genes=target_genes,
            gene_col="gene",
            strata_cols=("subtype", "comparison"),
            z_threshold=0.5,
            n_permutations=N_PERMUTATIONS,
            seed=SEED,
        )
        null_table.insert(0, "pair", pair)
        null_tables.append(null_table)
        summary_rows.append({"pair": pair, **summary})
    return (
        pd.DataFrame(summary_rows),
        pd.concat(null_tables, ignore_index=True),
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    metrics = pd.read_csv(METRICS)
    panel = pd.read_csv(PANEL)
    controls = pd.read_csv(CONTROLS)

    thresholds = threshold_sensitivity(
        metrics,
        panel,
        controls,
    )
    thresholds.to_csv(
        OUT_DIR / "threshold_sensitivity.csv",
        index=False,
    )

    permutation_summary, permutation_null = pairing_permutation(
        metrics,
        panel,
    )
    permutation_summary.to_csv(
        OUT_DIR / "pairing_permutation_summary.csv",
        index=False,
    )
    permutation_null.to_csv(
        OUT_DIR / "pairing_permutation_null.csv",
        index=False,
    )

    manifest = {
        "analysis": "unitaware effect-aware robustness",
        "seed": SEED,
        "thresholds": list(THRESHOLDS),
        "n_random_sets": N_RANDOM_SETS,
        "n_pairing_permutations": N_PERMUTATIONS,
        "input_sha256": {
            "gene_unit_metrics_long.csv": sha256(METRICS),
            "target_panel_400.csv": sha256(PANEL),
            "negative_control_pool_2000.csv": sha256(CONTROLS),
        },
        "permutation_definition": (
            "Within subtype and comparison, right-unit effects were "
            "reassigned across genes without replacement. Left-unit "
            "effects and all effect-size marginals were preserved."
        ),
    }
    (OUT_DIR / "robustness_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(thresholds.to_string(index=False))
    print()
    print(permutation_summary.to_string(index=False))


if __name__ == "__main__":
    main()
