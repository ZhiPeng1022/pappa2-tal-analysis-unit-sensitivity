"""Generate v2 manuscript figures 2 and 3.

Figure 2: TAL subtype disease-associated neighborhood remodeling.
Figure 3: PAPPA2 donor-level pseudobulk elevation.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2")
OUT = ROOT / "outputs"

SUBCLASSES = ["C-TAL-A", "C-TAL-B", "aTAL2", "frTAL"]
COMPARISONS = ["disease_vs_Reference", "AKI_vs_Reference", "CKD_vs_Reference"]
COMPARISON_LABELS = {
    "disease_vs_Reference": "Disease",
    "AKI_vs_Reference": "AKI",
    "CKD_vs_Reference": "CKD",
}

POSITIVE = "#2A6F97"
NEGATIVE = "#C44E52"
NEUTRAL = "#6C757D"
COMPARISON_COLORS = {
    "disease_vs_Reference": "#2A6F97",
    "AKI_vs_Reference": "#E07A5F",
    "CKD_vs_Reference": "#3D9970",
}


def apply_style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
        }
    )


def load_tables():
    milo = pd.read_csv(OUT / "KPMP_milo_TAL_subclass_summary.csv")
    augur = pd.read_csv(OUT / "KPMP_Augur_TAL_full_AUC_summary.csv")
    pseudo = pd.read_csv(OUT / "KPMP_TAL_focus_gene_pseudobulk.csv")
    egfr_models = pd.read_csv(OUT / "KPMP_PAPPA2_TAL_egfr_models.csv")
    egfr_specificity = pd.read_csv(OUT / "KPMP_PAPPA2_TAL_egfr_gene_specificity.csv")
    egfr_linear = pd.read_csv(OUT / "KPMP_PAPPA2_TAL_egfr_linear_sensitivity.csv")
    nephro_genes = pd.read_csv(OUT / "nephroseq_GFR_gene_summary.csv")
    nephro_sig = pd.read_csv(OUT / "nephroseq_GFR_significant_analyses.csv")
    gse_summary = pd.read_csv(OUT / "GSE137570_PAPPA2_eGFR_validation.csv")
    gse_data = pd.read_csv(ROOT / "work" / "GSE137570_PAPPA2_eGFR_validation_data.csv")
    spatial = pd.read_csv(
        OUT / "KPMP_Visium_PAPPA2_branch_signature_sample_associations.csv"
    )
    robustness = pd.read_csv(OUT / "KPMP_TAL_seed20260927_nhood_focus_robustness.csv")
    within = pd.read_csv(OUT / "KPMP_TAL_seed20260927_within_donor_cell_rank.csv")
    threshold = pd.read_csv(
        OUT / "KPMP_TAL_seed20260927_within_donor_cell_rank_threshold_sensitivity.csv"
    )
    return (
        milo,
        augur,
        pseudo,
        egfr_models,
        egfr_specificity,
        egfr_linear,
        nephro_genes,
        nephro_sig,
        gse_summary,
        gse_data,
        spatial,
        robustness,
        within,
        threshold,
    )


def figure2(milo, augur):
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(10.6, 3.5),
        gridspec_kw={"width_ratios": [1.05, 1.0, 1.35]},
    )

    # Panel a: combined-disease significant neighborhoods
    ax = axes[0]
    combined = milo[milo["comparison"] == "all_disease"].set_index("subclass")
    x = np.arange(len(SUBCLASSES))
    positive = [combined.loc[s, "positive_nhoods"] for s in SUBCLASSES]
    negative = [combined.loc[s, "negative_nhoods"] for s in SUBCLASSES]
    total = [combined.loc[s, "total_annotated_nhoods"] for s in SUBCLASSES]
    significant = [combined.loc[s, "significant_nhoods"] for s in SUBCLASSES]
    ax.bar(x, positive, color=POSITIVE, label="Positive")
    ax.bar(x, negative, bottom=positive, color=NEGATIVE, label="Negative")
    for index, (sig, tot) in enumerate(zip(significant, total)):
        ax.text(index, sig + 1.5, f"{sig}/{tot}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(SUBCLASSES, rotation=20, ha="right")
    ax.set_ylabel("Significant neighborhoods")
    ax.set_title("Combined disease")
    ax.set_ylim(0, max(total) + 12)
    ax.legend(frameon=False, loc="upper right")

    # Panel b: median logFC heatmap
    ax = axes[1]
    matrix = np.array(
        [
            [
                milo[
                    (milo["comparison"] == comparison)
                    & (milo["subclass"] == subclass)
                ]["median_logFC"].iloc[0]
                for subclass in SUBCLASSES
            ]
            for comparison in ["AKI_vs_Reference", "CKD_vs_Reference"]
        ]
    )
    limit = 2.4
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(np.arange(len(SUBCLASSES)))
    ax.set_xticklabels(SUBCLASSES, rotation=20, ha="right")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["AKI", "CKD"])
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix[row, column]
            ax.text(
                column,
                row,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=7,
                color="black",
            )
    ax.set_title("Median neighborhood logFC")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.ax.tick_params(labelsize=6.5)

    # Panel c: Augur AUC
    ax = axes[2]
    offsets = np.linspace(-0.24, 0.24, len(COMPARISONS))
    for offset, comparison in zip(offsets, COMPARISONS):
        subset = augur[augur["comparison"] == comparison].set_index("subclass")
        positions = np.arange(len(SUBCLASSES)) + offset
        values = [subset.loc[s, "auc"] for s in SUBCLASSES]
        lower = [subset.loc[s, "auc"] - subset.loc[s, "auc_q025"] for s in SUBCLASSES]
        upper = [subset.loc[s, "auc_q975"] - subset.loc[s, "auc"] for s in SUBCLASSES]
        ax.errorbar(
            positions,
            values,
            yerr=[lower, upper],
            fmt="o",
            markersize=4,
            capsize=2,
            linewidth=1,
            color=COMPARISON_COLORS[comparison],
            label=COMPARISON_LABELS[comparison],
        )
    ax.axhline(0.5, color=NEUTRAL, linestyle="--", linewidth=0.8)
    ax.set_xticks(np.arange(len(SUBCLASSES)))
    ax.set_xticklabels(SUBCLASSES, rotation=20, ha="right")
    ax.set_ylabel("Augur AUC")
    ax.set_title("Disease vs Reference separability")
    ax.set_ylim(0.15, 1.0)
    ax.legend(frameon=False, loc="upper left", ncol=1)

    fig.tight_layout()
    output = OUT / "Fig2_TAL_neighborhood_remodeling.png"
    fig.savefig(output)
    plt.close(fig)
    return output


def figure3(pseudo):
    pappa2 = pseudo[pseudo["gene"] == "PAPPA2"].copy()
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(10.6, 3.5),
        gridspec_kw={"width_ratios": [1.25, 1.0, 1.1]},
    )

    # Panel a: effect sizes
    ax = axes[0]
    x = np.arange(len(SUBCLASSES))
    width = 0.25
    for index, comparison in enumerate(COMPARISONS):
        subset = pappa2[pappa2["comparison"] == comparison].set_index("subclass")
        values = [subset.loc[s, "mean_case_minus_reference"] for s in SUBCLASSES]
        ax.bar(
            x + (index - 1) * width,
            values,
            width,
            color=COMPARISON_COLORS[comparison],
            label=COMPARISON_LABELS[comparison],
        )
    ax.set_xticks(x)
    ax.set_xticklabels(SUBCLASSES, rotation=20, ha="right")
    ax.set_ylabel("Mean case - reference\n(log1p CP10K)")
    ax.set_title("PAPPA2 donor pseudobulk")
    ax.legend(frameon=False, loc="upper left")
    ax.axhline(0, color=NEUTRAL, linewidth=0.7)

    # Panel b: FDR heatmap
    ax = axes[1]
    matrix = np.array(
        [
            [
                pappa2[
                    (pappa2["comparison"] == comparison)
                    & (pappa2["subclass"] == subclass)
                ]["FDR"].iloc[0]
                for subclass in SUBCLASSES
            ]
            for comparison in COMPARISONS
        ]
    )
    image = ax.imshow(-np.log10(matrix), cmap="Blues", aspect="auto")
    ax.set_xticks(np.arange(len(SUBCLASSES)))
    ax.set_xticklabels(SUBCLASSES, rotation=20, ha="right")
    ax.set_yticks(np.arange(len(COMPARISONS)))
    ax.set_yticklabels([COMPARISON_LABELS[c] for c in COMPARISONS])
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix[row, column]
            label = f"{value:.1e}" if value < 0.001 else f"{value:.3f}"
            ax.text(column, row, label, ha="center", va="center", fontsize=6.5)
    ax.set_title("FDR")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.ax.tick_params(labelsize=6.5)

    # Panel c: positive fraction, disease comparison
    ax = axes[2]
    disease = pappa2[pappa2["comparison"] == "disease_vs_Reference"].set_index("subclass")
    reference = [disease.loc[s, "reference_positive_fraction"] * 100 for s in SUBCLASSES]
    case = [disease.loc[s, "case_positive_fraction"] * 100 for s in SUBCLASSES]
    ax.bar(x - 0.18, reference, 0.36, color="#B0B7BD", label="Reference")
    ax.bar(x + 0.18, case, 0.36, color=POSITIVE, label="Disease")
    ax.set_xticks(x)
    ax.set_xticklabels(SUBCLASSES, rotation=20, ha="right")
    ax.set_ylabel("PAPPA2-positive donors (%)")
    ax.set_ylim(0, 105)
    ax.set_title("Combined disease positive fraction")
    ax.legend(frameon=False, loc="lower right")

    fig.tight_layout()
    output = OUT / "Fig3_PAPPA2_pseudobulk.png"
    fig.savefig(output)
    plt.close(fig)
    return output


def figure4(robustness, within, threshold):
    comparison_order = ["all_disease", "AKI_vs_Reference", "CKD_vs_Reference"]
    comparison_labels = {
        "all_disease": "All",
        "AKI_vs_Reference": "AKI",
        "CKD_vs_Reference": "CKD",
    }
    group_labels = [
        f"{comparison_labels[c]} {s}" for c in comparison_order for s in SUBCLASSES
    ]

    robust = robustness[robustness["gene"] == "PAPPA2"].copy()
    robust["comparison"] = pd.Categorical(
        robust["comparison"], categories=comparison_order, ordered=True
    )
    robust["subclass"] = pd.Categorical(
        robust["subclass"], categories=SUBCLASSES, ordered=True
    )
    robust = robust.sort_values(["comparison", "subclass"]).reset_index(drop=True)

    donor = within[within["gene"] == "PAPPA2"].copy()
    donor["comparison"] = pd.Categorical(
        donor["comparison"], categories=comparison_order, ordered=True
    )
    donor["subclass"] = pd.Categorical(
        donor["subclass"], categories=SUBCLASSES, ordered=True
    )
    donor = donor.sort_values(["comparison", "subclass"]).reset_index(drop=True)

    fig, axes = plt.subplots(2, 3, figsize=(11.5, 7.0))
    x = np.arange(len(group_labels))

    # Panel a: neighborhood-level rho with donor-stratified interval
    ax = axes[0, 0]
    for xi, rho, low, high in zip(
        x,
        robust["rho_full"],
        robust["bootstrap_ci_lower"],
        robust["bootstrap_ci_upper"],
    ):
        color = POSITIVE if low > 0 else (NEGATIVE if high < 0 else NEUTRAL)
        ax.plot([xi, xi], [low, high], color=color, linewidth=1.2)
        ax.plot([xi], [rho], marker="o", color=color, markersize=4)
    ax.axhline(0, color=NEUTRAL, linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(group_labels, rotation=90, fontsize=5.5)
    ax.set_ylabel("Spearman rho")
    ax.set_title("Neighborhood-level association\n(PAPPA2 mean expression vs Milo logFC)")

    # Panel b: purity-adjusted correlation
    ax = axes[0, 1]
    ax.scatter(
        robust["rho_full"],
        robust["partial_rho_controlling_purity"],
        s=22,
        color=POSITIVE,
    )
    limits = [-0.3, 0.9]
    ax.plot(limits, limits, color=NEUTRAL, linestyle="--", linewidth=0.8)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_xlabel("Unadjusted rho")
    ax.set_ylabel("Partial rho (purity adjusted)")
    ax.set_title("Neighborhood-level purity adjustment")

    # Panel c: within-donor cell-level association
    ax = axes[0, 2]
    for xi, rho, fdr in zip(
        x, donor["observed_median_rho_all"], donor["FDR_p_all"]
    ):
        color = POSITIVE if rho > 0 else NEGATIVE
        marker = "o" if fdr < 0.05 else "o"
        face = color if fdr < 0.05 else "white"
        ax.plot(
            [xi],
            [rho],
            marker=marker,
            markersize=5,
            markerfacecolor=face,
            markeredgecolor=color,
        )
    ax.axhline(0, color=NEUTRAL, linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(group_labels, rotation=90, fontsize=5.5)
    ax.set_ylabel("Group-level median rho")
    ax.set_title("Within-donor cell-level association\n(filled = FDR < 0.05)")

    # Panel d: cell-count threshold sensitivity
    ax = axes[1, 0]
    thresholds = [3, 5, 10]
    pappa2_threshold = threshold[
        (threshold["gene"] == "PAPPA2")
        & (threshold["min_cells_per_donor"].isin(thresholds))
    ]
    non_atal = []
    atal = []
    for value in thresholds:
        subset = pappa2_threshold[pappa2_threshold["min_cells_per_donor"] == value]
        non_atal.append(
            int(
                (
                    (subset["subclass"] != "aTAL2")
                    & (subset["observed_median_rho_all"] > 0)
                ).sum()
            )
        )
        atal.append(
            int(
                (
                    (subset["subclass"] == "aTAL2")
                    & (subset["observed_median_rho_all"] < 0)
                ).sum()
            )
        )
    positions = np.arange(len(thresholds))
    ax.bar(positions - 0.18, non_atal, 0.36, color=POSITIVE, label="Non-aTAL2 positive")
    ax.bar(positions + 0.18, atal, 0.36, color=NEGATIVE, label="aTAL2 negative")
    ax.set_xticks(positions)
    ax.set_xticklabels([f">={t} cells" for t in thresholds])
    ax.set_ylabel("Number of groups")
    ax.set_ylim(0, 10.5)
    ax.set_title("Within-donor threshold sensitivity")
    ax.legend(frameon=False, loc="upper left")

    # Panel e: aTAL2 direction contrast
    ax = axes[1, 1]
    robust_index = robust.set_index(["comparison", "subclass"])
    donor_index = donor.set_index(["comparison", "subclass"])
    positions = np.arange(len(comparison_order))
    neighborhood_values = [
        robust_index.loc[(comparison, "aTAL2"), "rho_full"]
        for comparison in comparison_order
    ]
    donor_values = [
        donor_index.loc[(comparison, "aTAL2"), "observed_median_rho_all"]
        for comparison in comparison_order
    ]
    ax.bar(
        positions - 0.18,
        neighborhood_values,
        0.36,
        color=POSITIVE,
        label="Neighborhood level",
    )
    ax.bar(
        positions + 0.18,
        donor_values,
        0.36,
        color=NEGATIVE,
        label="Within-donor cells",
    )
    ax.axhline(0, color=NEUTRAL, linewidth=0.8)
    ax.set_xticks(positions)
    ax.set_xticklabels([comparison_labels[c] for c in comparison_order])
    ax.set_ylabel("rho")
    ax.set_title("aTAL2 depends on analysis unit")
    ax.legend(frameon=False, loc="lower left")

    # Panel f: eligible donors
    ax = axes[1, 2]
    reference = donor["eligible_reference_donors"].values
    disease = donor["eligible_disease_donors"].values
    ax.bar(x, reference, color="#B0B7BD", label="Reference")
    ax.bar(x, disease, bottom=reference, color=COMPARISON_COLORS["disease_vs_Reference"], label="Disease")
    ax.set_xticks(x)
    ax.set_xticklabels(group_labels, rotation=90, fontsize=5.5)
    ax.set_ylabel("Eligible donors")
    ax.set_title("Within-donor eligible sample size")
    ax.legend(frameon=False, loc="upper right")

    fig.tight_layout()
    output = OUT / "Fig4_two_level_association.png"
    fig.savefig(output)
    plt.close(fig)
    return output


def _forest(ax, labels, estimates, lower, upper, colors, xlabel, logx=False, reference=1.0):
    y = np.arange(len(labels))[::-1]
    for yi, estimate, low, high, color in zip(y, estimates, lower, upper, colors):
        ax.plot([low, high], [yi, yi], color=color, linewidth=1.4)
        ax.plot([estimate], [yi], marker="o", color=color, markersize=4.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    if logx:
        ax.set_xscale("log")
    ax.axvline(reference, color=NEUTRAL, linestyle="--", linewidth=0.8)
    ax.set_xlabel(xlabel)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def figure5(egfr_models, egfr_specificity, egfr_linear):
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(11.2, 3.9),
        gridspec_kw={"width_ratios": [1.25, 1.0, 1.05]},
    )

    # Panel a: PAPPA2 ordinal eGFR models
    rows = [
        ("All, model 2", "all", "model_2"),
        ("All, model 3", "all", "model_3"),
        ("All, model 4", "all", "model_4"),
        ("Disease, model 2", "disease", "model_2"),
        ("CKD, model 2", "CKD", "model_2"),
        ("AKI, model 2", "AKI", "model_2"),
    ]
    scope_colors = {
        "all": "#2A6F97",
        "disease": "#E07A5F",
        "CKD": "#3D9970",
        "AKI": "#8E7DBE",
    }
    labels, estimates, lower, upper, colors = [], [], [], [], []
    for label, scope, model in rows:
        row = egfr_models[
            (egfr_models["scope"] == scope) & (egfr_models["model"] == model)
        ].iloc[0]
        labels.append(label)
        estimates.append(row["odds_ratio_per_10pct"])
        lower.append(row["confidence_lower"])
        upper.append(row["confidence_upper"])
        colors.append(scope_colors[scope])
    ax = axes[0]
    _forest(
        ax,
        labels,
        estimates,
        lower,
        upper,
        colors,
        "Odds ratio per 10% PAPPA2-positive fraction",
        logx=True,
        reference=1.0,
    )
    ax.set_title("KPMP ordinal eGFR models")

    # Panel b: candidate-gene specificity at all, model 2
    subset = egfr_specificity[
        (egfr_specificity["scope"] == "all") & (egfr_specificity["model"] == "model_2")
    ].sort_values("odds_ratio_per_10pct")
    gene_labels = subset["gene"].tolist()
    gene_colors = [POSITIVE if gene == "PAPPA2" else NEUTRAL for gene in gene_labels]
    ax = axes[1]
    _forest(
        ax,
        gene_labels,
        subset["odds_ratio_per_10pct"].tolist(),
        subset["confidence_lower"].tolist(),
        subset["confidence_upper"].tolist(),
        gene_colors,
        "Odds ratio per 10% positive fraction",
        logx=True,
        reference=1.0,
    )
    ax.set_title("Candidate-gene specificity (all donors)")

    # Panel c: linear sensitivity on eGFR midpoint
    scopes = ["all", "disease", "CKD", "AKI"]
    linear = egfr_linear[egfr_linear["model"] == "model_2"].set_index("scope")
    ax = axes[2]
    _forest(
        ax,
        scopes,
        [linear.loc[s, "coefficient"] for s in scopes],
        [linear.loc[s, "confidence_lower"] for s in scopes],
        [linear.loc[s, "confidence_upper"] for s in scopes],
        [scope_colors[s] for s in scopes],
        "eGFR midpoint change per 10% (mL/min/1.73 m2)",
        logx=False,
        reference=0.0,
    )
    ax.set_title("Linear sensitivity analysis")

    fig.tight_layout()
    output = OUT / "Fig5_KPMP_eGFR_specificity.png"
    fig.savefig(output)
    plt.close(fig)
    return output


def figure6(spatial):
    group_order = ["Reference", "AKI", "DKD"]
    group_colors = {"Reference": "#6C757D", "AKI": "#E07A5F", "DKD": "#2A6F97"}
    spatial = spatial.copy()
    spatial["disease_group"] = pd.Categorical(
        spatial["disease_group"], categories=group_order, ordered=True
    )
    spatial = spatial.sort_values(["disease_group", "sample"]).reset_index(drop=True)

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(11.0, 3.8),
        gridspec_kw={"width_ratios": [1.55, 1.0, 1.05]},
    )

    # Panel a: per-section co-occurrence
    ax = axes[0]
    x = np.arange(len(spatial))
    colors = [group_colors[group] for group in spatial["disease_group"]]
    ax.bar(x, spatial["high_state_cooccurrence_odds_ratio"], color=colors)
    ax.axhline(1.0, color=NEUTRAL, linestyle="--", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(spatial["sample"], rotation=90, fontsize=5)
    ax.set_ylabel("Co-occurrence odds ratio")
    above = int((spatial["high_state_cooccurrence_odds_ratio"] > 1).sum())
    ax.set_title(f"Per-section co-occurrence ({above}/{len(spatial)} > 1)")
    handles = [Patch(facecolor=group_colors[g], label=g) for g in group_order]
    ax.legend(handles=handles, frameon=False, loc="upper right", fontsize=6)
    ax.text(
        0.02,
        0.02,
        "Donor-deduplicated DKD: 10/10 > 1",
        transform=ax.transAxes,
        fontsize=6,
        va="bottom",
    )

    # Panel b: group distribution
    ax = axes[1]
    values = [
        spatial.loc[
            spatial["disease_group"] == group, "high_state_cooccurrence_odds_ratio"
        ].values
        for group in group_order
    ]
    box = ax.boxplot(
        values,
        tick_labels=group_order,
        widths=0.5,
        patch_artist=True,
        showfliers=False,
    )
    for patch, group in zip(box["boxes"], group_order):
        patch.set_facecolor(group_colors[group])
        patch.set_alpha(0.35)
    rng = np.random.default_rng(1)
    for index, (group, group_values) in enumerate(zip(group_order, values), start=1):
        jitter = rng.normal(0, 0.04, len(group_values))
        ax.scatter(
            np.full(len(group_values), index) + jitter,
            group_values,
            s=16,
            color=group_colors[group],
            zorder=3,
        )
        ax.text(
            index,
            max(group_values) + 0.06,
            f"{int((group_values > 1).sum())}/{len(group_values)}",
            ha="center",
            fontsize=6.5,
        )
    ax.axhline(1.0, color=NEUTRAL, linestyle="--", linewidth=0.8)
    ax.set_ylabel("Co-occurrence odds ratio")
    ax.set_title("By disease group")

    # Panel c: signature correlation vs co-occurrence
    ax = axes[2]
    for group in group_order:
        subset = spatial[spatial["disease_group"] == group]
        ax.scatter(
            subset["spearman_main_repair"],
            subset["high_state_cooccurrence_odds_ratio"],
            s=18,
            color=group_colors[group],
            label=group,
        )
    ax.axhline(1.0, color=NEUTRAL, linestyle="--", linewidth=0.8)
    ax.set_xlabel("Main-repair signature Spearman rho")
    ax.set_ylabel("Co-occurrence odds ratio")
    ax.set_title("Correlation vs co-occurrence")
    ax.legend(frameon=False, loc="upper left", fontsize=6)

    fig.tight_layout()
    output = OUT / "Fig6_KPMP_spatial_crossmodality.png"
    fig.savefig(output)
    plt.close(fig)
    return output


def figure7(nephro_genes, nephro_sig, gse_summary, gse_data):
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.0))

    # Panel a: Nephroseq PAPPA2 significant GFR analyses
    ax = axes[0, 0]
    pappa2 = nephro_sig[nephro_sig["gene"] == "PAPPA2"].copy()
    pappa2["tissue_types"] = pappa2["tissue_types"].astype(str).str.strip()
    pappa2 = pappa2.sort_values("r_value")
    labels = [
        f"{row.dataset_name} ({row.tissue_types})" for row in pappa2.itertuples()
    ]
    tissue_colors = {
        "Tubulointerstitium": NEGATIVE,
        "Kidney": POSITIVE,
        "Glomerulus": NEUTRAL,
    }
    colors = [tissue_colors.get(tissue, NEUTRAL) for tissue in pappa2["tissue_types"]]
    _forest(
        ax,
        labels,
        pappa2["r_value"].tolist(),
        pappa2["r_value"].tolist(),
        pappa2["r_value"].tolist(),
        colors,
        "Spearman r (GFR)",
        logx=False,
        reference=0.0,
    )
    ax.set_title("Nephroseq PAPPA2 GFR analyses")
    ax.tick_params(axis="y", labelsize=6.5)

    # Panel b: negative tubulointerstitial analyses per gene
    ax = axes[0, 1]
    genes = ["PAPPA2", "CXCL12", "NOTCH3", "HES1", "EPHB2"]
    subset = nephro_sig.copy()
    subset["tissue_types"] = subset["tissue_types"].astype(str).str.strip()
    negative = (
        subset[
            (subset["tissue_types"] == "Tubulointerstitium") & (subset["r_value"] < 0)
        ]
        .groupby("gene")
        .size()
        .reindex(genes)
        .fillna(0)
    )
    bar_colors = [POSITIVE if gene == "PAPPA2" else NEUTRAL for gene in genes]
    ax.bar(genes, negative.values, color=bar_colors)
    for index, value in enumerate(negative.values):
        ax.text(index, value + 0.1, str(int(value)), ha="center", va="bottom", fontsize=7)
    ax.set_ylabel("Negative tubulointerstitial analyses")
    ax.set_title("Nephroseq negative tubulointerstitial results")
    ax.set_ylim(0, max(negative.values) + 1.5)

    # Panels c and d: GSE137570 public eGFR scatter
    for ax, gene in zip(axes[1, :], ["PAPPA2", "EPHB2"]):
        values = gse_data[gene].astype(float)
        egfr = gse_data["egfr"].astype(float)
        ax.scatter(values, egfr, s=18, color=POSITIVE if gene == "PAPPA2" else NEUTRAL)
        coefficients = np.polyfit(values, egfr, 1)
        x_grid = np.linspace(values.min(), values.max(), 50)
        ax.plot(x_grid, np.polyval(coefficients, x_grid), color=NEGATIVE, linewidth=1)
        summary_row = gse_summary[gse_summary["gene"] == gene].iloc[0]
        ax.set_xlabel(f"{gene} expression")
        ax.set_ylabel("eGFR (mL/min/1.73 m2)")
        ax.set_title(
            f"GSE137570 {gene}: rho={summary_row['spearman_rho']:.3f}, "
            f"p={summary_row['spearman_p']:.3f}"
        )

    fig.tight_layout()
    output = OUT / "Fig7_external_GFR_associations.png"
    fig.savefig(output)
    plt.close(fig)
    return output


def _schematic_box(ax, x, y, width, height, text, face="#EAF2F8", edge="#2A6F97", fontsize=7.2):
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.008,rounding_size=0.015",
        linewidth=0.9,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(box)
    ax.text(
        x + width / 2,
        y + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
    )


def _schematic_arrow(ax, x1, y1, x2, y2, color="#6C757D"):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=0.9,
            color=color,
        )
    )


def figure1():
    fig, axes = plt.subplots(
        1, 3, figsize=(11.5, 4.2), gridspec_kw={"width_ratios": [1.45, 1.25, 1.0]}
    )
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

    # Panel a: data roles
    ax = axes[0]
    ax.set_title("Data roles", fontsize=9)
    boxes = [
        (0.02, 0.48, 0.44, 0.30, "KPMP discovery\nsnRNA-seq\nGSE183277", "#E3EEF7", "#2A6F97"),
        (0.54, 0.48, 0.44, 0.30, "KPMP spatial\nVisium\nGSE183456", "#E3EEF7", "#2A6F97"),
        (0.02, 0.08, 0.44, 0.28, "External public\nGSE137570\n24 donors", "#F7EDE8", "#E07A5F"),
        (0.54, 0.08, 0.44, 0.28, "Nephroseq v5\nprecomputed\nGFR analyses", "#F7EDE8", "#E07A5F"),
    ]
    for x, y, width, height, text, face, edge in boxes:
        _schematic_box(ax, x, y, width, height, text, face=face, edge=edge)
    _schematic_arrow(ax, 0.46, 0.63, 0.54, 0.63)
    _schematic_arrow(ax, 0.46, 0.22, 0.54, 0.22)
    ax.text(0.50, 0.90, "Same KPMP consortium", ha="center", fontsize=7.5, color="#2A6F97")
    ax.text(0.50, 0.42, "External data", ha="center", fontsize=7.5, color="#E07A5F")

    # Panel b: analysis flow
    ax = axes[1]
    ax.set_title("Analysis flow", fontsize=9)
    steps = [
        (0.02, 0.70, 0.44, 0.22, "Milo neighborhood\nabundance"),
        (0.54, 0.70, 0.44, 0.22, "Augur disease\nseparability"),
        (0.02, 0.40, 0.44, 0.22, "Donor pseudobulk\nPAPPA2 expression"),
        (0.54, 0.40, 0.44, 0.22, "Neighborhood-level\nassociation"),
        (0.02, 0.10, 0.44, 0.22, "Within-donor\ncell-level rank"),
        (0.54, 0.10, 0.44, 0.22, "KPMP eGFR\nordinal model"),
    ]
    for x, y, width, height, text in steps:
        _schematic_box(ax, x, y, width, height, text)
    _schematic_arrow(ax, 0.46, 0.81, 0.54, 0.81)
    _schematic_arrow(ax, 0.76, 0.70, 0.76, 0.62)
    _schematic_arrow(ax, 0.54, 0.51, 0.46, 0.51)
    _schematic_arrow(ax, 0.24, 0.40, 0.24, 0.32)
    _schematic_arrow(ax, 0.46, 0.21, 0.54, 0.21)

    # Panel c: gene panel
    ax = axes[2]
    ax.set_title("Candidate and control genes", fontsize=9)
    _schematic_box(
        ax,
        0.10,
        0.68,
        0.80,
        0.20,
        "PAPPA2\ncandidate",
        face="#D6EAF8",
        edge="#2A6F97",
        fontsize=8,
    )
    for index, gene in enumerate(["CXCL12", "NOTCH3", "HES1", "EPHB2"]):
        _schematic_box(
            ax,
            0.10,
            0.52 - index * 0.13,
            0.80,
            0.10,
            gene,
            face="#F1F3F5",
            edge="#6C757D",
            fontsize=7.5,
        )

    fig.tight_layout()
    output = OUT / "Fig1_study_design.png"
    fig.savefig(output)
    plt.close(fig)
    return output


def figure8():
    fig, ax = plt.subplots(figsize=(11.5, 5.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(
        0.5,
        0.96,
        "Subtype-specific TAL neighborhood state model",
        ha="center",
        fontsize=10.5,
        color="#1F2933",
    )

    # Left: TAL subtype states
    _schematic_box(
        ax,
        0.03,
        0.56,
        0.27,
        0.28,
        "TAL subtypes\nC-TAL-A\nC-TAL-B\nfrTAL",
        face="#D6EAF8",
        edge="#2A6F97",
        fontsize=8,
    )
    ax.text(
        0.165,
        0.515,
        "positive neighborhood association",
        ha="center",
        fontsize=7,
        color="#2A6F97",
    )
    _schematic_box(
        ax,
        0.03,
        0.20,
        0.27,
        0.24,
        "aTAL2\n(exploratory\ndecoupling)",
        face="#F7EDE8",
        edge="#E07A5F",
        fontsize=8,
    )

    # Middle: framework
    _schematic_box(
        ax,
        0.36,
        0.36,
        0.27,
        0.48,
        "Donor-aware\nneighborhood\nframework\n\nMilo + Augur\npseudobulk\nwithin-donor rank",
        face="#E8F4EC",
        edge="#3D9970",
        fontsize=8,
    )

    # Right: PAPPA2 case
    _schematic_box(
        ax,
        0.70,
        0.56,
        0.27,
        0.28,
        "PAPPA2-associated\nTAL state",
        face="#D6EAF8",
        edge="#2A6F97",
        fontsize=8,
    )
    _schematic_box(
        ax,
        0.70,
        0.20,
        0.27,
        0.24,
        "Lower kidney\nfunction\n(eGFR)",
        face="#F7EDE8",
        edge="#E07A5F",
        fontsize=8,
    )
    _schematic_arrow(ax, 0.30, 0.70, 0.36, 0.70)
    _schematic_arrow(ax, 0.63, 0.70, 0.70, 0.70)
    _schematic_arrow(ax, 0.835, 0.56, 0.835, 0.44)
    ax.text(0.665, 0.735, "case gene", ha="center", fontsize=7, color="#6C757D")
    ax.text(0.86, 0.50, "associated", ha="center", fontsize=7, color="#6C757D")

    # Bottom: supported and not supported
    _schematic_box(
        ax,
        0.03,
        0.03,
        0.45,
        0.13,
        "Supported: subtype-specific state structure,\n"
        "exploratory aTAL2 decoupling, PAPPA2 association",
        face="#E8F4EC",
        edge="#3D9970",
        fontsize=7,
    )
    _schematic_box(
        ax,
        0.52,
        0.03,
        0.45,
        0.13,
        "Not supported: causality, direct communication,\n"
        "unified JAG1-NOTCH3-HES1 or EPHB2-STAT6 mechanism",
        face="#FDECEC",
        edge="#C44E52",
        fontsize=7,
    )
    fig.tight_layout()
    output = OUT / "Fig8_TAL_state_model.png"
    fig.savefig(output)
    plt.close(fig)
    return output


def main():
    apply_style()
    (
        milo,
        augur,
        pseudo,
        egfr_models,
        egfr_specificity,
        egfr_linear,
        nephro_genes,
        nephro_sig,
        gse_summary,
        gse_data,
        spatial,
        robustness,
        within,
        threshold,
    ) = load_tables()
    outputs = [
        figure1(),
        figure2(milo, augur),
        figure3(pseudo),
        figure4(robustness, within, threshold),
        figure5(egfr_models, egfr_specificity, egfr_linear),
        figure6(spatial),
        figure7(nephro_genes, nephro_sig, gse_summary, gse_data),
        figure8(),
    ]
    for path in outputs:
        print(f"{path} ({path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
