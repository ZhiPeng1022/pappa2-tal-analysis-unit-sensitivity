"""Create Supplementary Figure S1 for analysis-unit robustness."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
RESULT_DIR = ROOT / "outputs" / "prespecified_unit_panel_results"
ROBUSTNESS_DIR = ROOT / "outputs" / "unitaware_robustness"
OUTPUT = ROOT / "outputs" / "figures" / "Figure_S1_panel_validation.tif"

BLUE = "#5B8DB8"
RED = "#C96A6A"
ORANGE = "#D99A6C"
TEAL = "#4F9A94"
NEUTRAL = "#7B8794"
LIGHT = "#DCE6EE"


def style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 5.8,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def main() -> None:
    style()
    random = pd.read_csv(RESULT_DIR / "random_set_results.csv")
    primary = pd.read_json(RESULT_DIR / "primary_result.json", typ="series")
    permutation = pd.read_csv(
        ROBUSTNESS_DIR / "pairing_permutation_summary.csv"
    )
    thresholds = pd.read_csv(
        ROBUSTNESS_DIR / "threshold_sensitivity.csv"
    )

    figure = plt.figure(figsize=(170 / 25.4, 64 / 25.4))
    grid = figure.add_gridspec(
        1,
        3,
        width_ratios=[1.0, 1.2, 1.05],
        wspace=0.52,
        left=0.065,
        right=0.985,
        bottom=0.22,
        top=0.90,
    )

    # a, Expression matching alone creates an apparent target-panel excess.
    ax = figure.add_subplot(grid[0, 0])
    values = random["donor_within_random_400"].to_numpy(dtype=float)
    ax.hist(
        values,
        bins=24,
        color=LIGHT,
        edgecolor="white",
        linewidth=0.5,
    )
    target = float(primary["donor_within"]["target_400_mean"])
    biological = float(primary["donor_within"]["target_300_biological_mean"])
    ax.axvline(target, color=RED, lw=1.2)
    ax.axvline(biological, color=BLUE, lw=1.0, ls="--")
    ax.text(
        target,
        ax.get_ylim()[1] * 0.95,
        "400",
        color=RED,
        ha="left",
        va="top",
        rotation=90,
    )
    ax.text(
        biological,
        ax.get_ylim()[1] * 0.95,
        "300",
        color=BLUE,
        ha="right",
        va="top",
        rotation=90,
    )
    ax.set_xlabel("Donor-within reversal fraction")
    ax.set_ylabel("Random sets")
    ax.set_title(
        "a  Expression-matched null",
        loc="left",
        fontweight="bold",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # b, Re-pairing the same effects removes the apparent excess.
    ax = figure.add_subplot(grid[0, 1])
    pair_order = [
        "donor_within",
        "donor_neighborhood",
        "neighborhood_within",
    ]
    pair_labels = [
        "Donor-within",
        "Donor-neighborhood",
        "Neighborhood-within",
    ]
    pair_colors = [BLUE, ORANGE, TEAL]
    permutation = permutation.set_index("pair").reindex(pair_order)
    positions = np.arange(len(pair_order))
    for index, row in enumerate(permutation.itertuples()):
        color = pair_colors[index]
        ax.plot(
            [row.null_q025 * 100, row.null_q975 * 100],
            [index, index],
            color=color,
            lw=3.0,
            alpha=0.30,
            solid_capstyle="round",
            zorder=1,
        )
        ax.scatter(
            row.null_mean * 100,
            index,
            marker="|",
            s=80,
            color=color,
            linewidths=1.4,
            zorder=2,
        )
        ax.scatter(
            row.observed_reversal_fraction * 100,
            index,
            s=28,
            facecolor="white",
            edgecolor=color,
            linewidth=1.2,
            zorder=3,
        )
        label = (
            "P < 0.001"
            if row.empirical_p < 0.001
            else f"P = {row.empirical_p:.3f}"
        )
        ax.text(
            row.observed_reversal_fraction * 100,
            index - 0.25,
            label,
            ha="center",
            va="bottom",
            fontsize=5.8,
            color=NEUTRAL,
        )
    ax.set_yticks(positions)
    ax.set_yticklabels(pair_labels)
    ax.invert_yaxis()
    ax.set_xlabel("Reversal fraction (%)")
    ax.set_title(
        "b  Effect-aware conditional null",
        loc="left",
        fontweight="bold",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(
        0.02,
        0.03,
        "open circle: observed\nbar: null mean and 95% interval",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=5.6,
        color=NEUTRAL,
    )

    # c, Expression-matched separation persists across effect thresholds.
    ax = figure.add_subplot(grid[0, 2])
    for pair, color, label in [
        ("donor_within", BLUE, "Donor-within"),
        ("donor_neighborhood", ORANGE, "Donor-neighborhood"),
    ]:
        block = thresholds[thresholds["pair"].eq(pair)]
        ax.plot(
            block["threshold"],
            block["target_400"] * 100,
            color=color,
            lw=1.2,
            marker="o",
            ms=2.5,
            label=f"{label}, target",
        )
        ax.plot(
            block["threshold"],
            block["random_400"] * 100,
            color=color,
            lw=0.9,
            ls="--",
            marker="o",
            ms=2.0,
            mfc="white",
            label=f"{label}, matched",
        )
    ax.set_xlabel("Standardized-effect threshold")
    ax.set_ylabel("Reversal fraction (%)")
    ax.set_title(
        "c  Threshold sensitivity",
        loc="left",
        fontweight="bold",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper right", handlelength=1.8, borderaxespad=0.2)

    figure.savefig(
        OUTPUT,
        dpi=600,
        format="tiff",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(figure)
    print(OUTPUT)


if __name__ == "__main__":
    main()
