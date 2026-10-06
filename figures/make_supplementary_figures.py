"""Generate the supplementary schematic explaining the two analysis units."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch


OUT = Path(r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2\outputs")


def box(ax, x, y, width, height, text, face="#EAF2F8", edge="#2A6F97", fontsize=7.5):
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.01,rounding_size=0.02",
        linewidth=0.9,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=fontsize)


def arrow(ax, x1, y1, x2, y2, color="#6C757D"):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=9, linewidth=0.9, color=color
        )
    )


def figure_s1():
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for axis in axes:
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.axis("off")

    # Panel A: neighborhood-level analysis
    axis = axes[0]
    axis.text(0.5, 0.95, "A. Neighborhood-level analysis", ha="center", fontsize=9.5)
    for x, y in [(0.18, 0.62), (0.30, 0.62), (0.24, 0.72), (0.36, 0.72)]:
        axis.add_patch(
            Circle((x, y), 0.10, facecolor="#D6EAF8", edgecolor="#2A6F97", alpha=0.5, linewidth=0.8)
        )
    axis.text(0.27, 0.44, "Overlapping Milo neighborhoods", ha="center", fontsize=7.5)
    box(
        axis,
        0.55,
        0.58,
        0.40,
        0.22,
        "Unit = neighborhood\nmean PAPPA2 expression\nvs Milo logFC",
        face="#E8F4EC",
        edge="#3D9970",
    )
    arrow(axis, 0.46, 0.68, 0.55, 0.68)
    box(
        axis,
        0.10,
        0.14,
        0.80,
        0.22,
        "Statistics: Spearman correlation across neighborhoods\n"
        "500 donor-stratified resamples + purity-adjusted estimate",
        face="#F1F3F5",
        edge="#6C757D",
        fontsize=7,
    )

    # Panel B: within-donor cell-level analysis
    axis = axes[1]
    axis.text(0.5, 0.95, "B. Within-donor cell-level analysis", ha="center", fontsize=9.5)
    box(axis, 0.05, 0.55, 0.38, 0.30, "Donor 1\ncells", face="#D6EAF8", edge="#2A6F97")
    box(axis, 0.05, 0.16, 0.38, 0.30, "Donor 2\ncells", face="#D6EAF8", edge="#2A6F97")
    box(
        axis,
        0.55,
        0.36,
        0.40,
        0.30,
        "Unit = donor\none rho per donor\n(cell PAPPA2 vs\nneighborhood exposure)",
        face="#E8F4EC",
        edge="#3D9970",
    )
    arrow(axis, 0.43, 0.70, 0.55, 0.62)
    arrow(axis, 0.43, 0.31, 0.55, 0.42)
    box(
        axis,
        0.10,
        0.02,
        0.80,
        0.10,
        "Statistics: 1,000 within-donor permutations; 3, 5, 10 cell thresholds",
        face="#F1F3F5",
        edge="#6C757D",
        fontsize=7,
    )

    figure.text(
        0.5,
        0.02,
        "The two analyses use different units and must not be combined.",
        ha="center",
        fontsize=8,
        color="#C44E52",
    )
    figure.tight_layout(rect=[0, 0.05, 1, 1])
    output = OUT / "FigS1_analysis_units.png"
    figure.savefig(output, dpi=300)
    plt.close(figure)
    return output


if __name__ == "__main__":
    print(figure_s1())
