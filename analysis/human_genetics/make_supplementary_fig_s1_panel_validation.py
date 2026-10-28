"""Create Figure 9 for the prespecified analysis-unit panel."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
RESULT_DIR = ROOT / "outputs" / "prespecified_unit_panel_results"
PANEL_DIR = ROOT / "outputs" / "prespecified_unit_panel"
OUTPUT = ROOT / "outputs" / "figures" / "Figure_S1_panel_validation.tif"

BLUE = "#2A6F97"
RED = "#C44E52"
GREEN = "#3D9970"
ORANGE = "#D98C74"
NEUTRAL = "#6C757D"


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
            "legend.fontsize": 6.5,
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
    summary = pd.read_csv(RESULT_DIR / "gene_unit_summary.csv")
    panel = pd.read_csv(PANEL_DIR / "target_panel_400.csv")

    group_order = [
        "pappa2_coexpression",
        "tal_subtype_marker",
        "injury_inflammation_metabolism",
        "internal_negative_control",
        "external_control_pool",
    ]
    group_labels = [
        "PAPPA2 coexpression",
        "TAL subtype markers",
        "Injury / inflammation",
        "Internal controls",
        "External controls",
    ]
    group_values = (
        summary.groupby("panel_group")["donor_within_reversal_fraction"]
        .mean()
        .reindex(group_order)
        .to_numpy()
    )

    figure = plt.figure(figsize=(170 / 25.4, 64 / 25.4))
    grid = figure.add_gridspec(
        1,
        3,
        width_ratios=[1.0, 1.15, 0.95],
        wspace=0.48,
        left=0.065,
        right=0.985,
        bottom=0.20,
        top=0.90,
    )

    ax = figure.add_subplot(grid[0, 0])
    values = random["donor_within_random_400"].to_numpy(dtype=float)
    ax.hist(
        values,
        bins=24,
        color="#C3DCEC",
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
    ax.set_title("a  Random-set null", loc="left", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = figure.add_subplot(grid[0, 1])
    positions = np.arange(len(group_order))
    colors = [BLUE, "#5C8FB3", ORANGE, "#B7C4CE", NEUTRAL]
    ax.barh(
        positions,
        group_values * 100,
        color=colors,
        height=0.66,
    )
    random_mean = float(primary["donor_within"]["random_400_mean"]) * 100
    ax.axvline(random_mean, color=RED, lw=0.9, ls="--")
    ax.set_yticks(positions)
    ax.set_yticklabels(group_labels)
    ax.invert_yaxis()
    ax.set_xlabel("Donor-within reversal (%)")
    ax.set_title("b  Panel groups", loc="left", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for y, value in zip(positions, group_values):
        ax.text(
            value * 100 + 0.25,
            y,
            f"{value * 100:.1f}",
            va="center",
            fontsize=6.2,
        )

    ax = figure.add_subplot(grid[0, 2])
    labels = ["Donor–within", "Donor–neighborhood", "Neighborhood–within"]
    keys = ["donor_within", "donor_neighborhood", "neighborhood_within"]
    target_values = [
        float(primary[key]["target_400_mean"]) * 100
        for key in keys
    ]
    random_values = [
        float(primary[key]["stratified_400_mean"]) * 100
        for key in keys
    ]
    p_values = [
        float(primary[key]["empirical_p_stratified_400"])
        for key in keys
    ]
    y = np.arange(len(keys))
    for index in y:
        ax.plot(
            [random_values[index], target_values[index]],
            [index, index],
            color="#CED4DA",
            lw=2.0,
            zorder=1,
        )
    ax.scatter(random_values, y, color=NEUTRAL, s=22, zorder=2)
    ax.scatter(target_values, y, color=RED, s=26, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Reversal fraction (%)")
    ax.set_title("c  Unit-pair comparison", loc="left", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    x_limit = max(target_values + random_values) + 4
    ax.set_xlim(0, x_limit)
    for index, p_value in enumerate(p_values):
        text = "P = 0.001" if p_value <= 0.001 else f"P = {p_value:.3f}"
        ax.text(
            x_limit - 0.2,
            index - 0.22,
            text,
            ha="right",
            va="bottom",
            fontsize=6.0,
            color=NEUTRAL,
        )

    figure.savefig(
        OUTPUT,
        dpi=300,
        format="tiff",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(figure)
    print(OUTPUT)


if __name__ == "__main__":
    main()
