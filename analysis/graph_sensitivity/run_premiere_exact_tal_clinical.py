"""Explore exact TAL states and eGFR categories in the PREMIERE cohort.

The integrated Lake et al. atlas contains KPMP, BUKMAP, and PREMIERE cells.
The primary analysis uses PREMIERE donors whose IDs do not occur in KPMP.
eGFR is an ordered category, not a continuous measurement.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.formula.api import ols


ROOT = Path(r"C:\Users\Elsa\Documents\Codex\2026-10-15\ji")
H5AD = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-09-25\9\work"
    r"\KPMP_integrated.h5ad"
)
OUT = ROOT / "outputs" / "premiere_exact_tal_clinical" / "v1"
TAL_STATES = [
    "C-TAL",
    "M-TAL",
    "aTAL1",
    "aTAL2",
    "dC-TAL",
    "dM-TAL",
    "MD",
]
MIN_TAL_CELLS = 5
MIN_STATE_CELLS = 5
MIN_DONORS = 8
EGFR_RANK = {
    "20-29": 1,
    "30-39": 2,
    "40-49": 3,
    "50-59": 4,
    ">60": 5,
}


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(dtype=float)
    result = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return pd.Series(result, index=values.index)
    order = np.argsort(p[valid])
    ranked = p[valid][order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    adjusted = np.empty_like(q)
    adjusted[order] = np.clip(q, 0, 1)
    result[valid] = adjusted
    return pd.Series(result, index=values.index)


def dataset_label(experiment: str) -> str:
    for prefix in ["BUKMAP", "KPMP", "PREMIERE"]:
        if str(experiment).startswith(prefix):
            return prefix
    return "other"


def age_years(development_stage: str) -> float:
    text = str(development_stage)
    if "-year-old" in text:
        token = text.split("-year-old")[0].split()[-1]
        try:
            return float(token)
        except ValueError:
            return np.nan
    ordinal = {
        "first": 10,
        "second": 20,
        "third": 30,
        "fourth": 40,
        "fifth": 50,
        "sixth": 60,
        "seventh": 70,
        "eighth": 80,
        "ninth": 90,
    }
    for word, value in ordinal.items():
        if word in text:
            return float(value)
    return np.nan


def load_metadata() -> tuple[ad.AnnData, pd.DataFrame, int]:
    adata = ad.read_h5ad(H5AD, backed="r")
    obs = adata.obs.copy()
    obs["dataset"] = obs["experiment"].astype(str).map(dataset_label)
    obs["eGFR_rank"] = obs["eGFR"].astype(str).map(EGFR_RANK)
    obs["age_years"] = obs["development_stage"].map(age_years)
    pap_idx = int(
        np.flatnonzero(
            adata.var["feature_name"].astype(str).to_numpy()
            == "PAPPA2"
        )[0]
    )
    return adata, obs, pap_idx


def read_pappa2(
    adata: ad.AnnData,
    cell_indices: np.ndarray,
    pap_idx: int,
    chunk_size: int = 5000,
) -> np.ndarray:
    result = np.full(len(cell_indices), np.nan, dtype=float)
    for start in range(0, len(cell_indices), chunk_size):
        stop = min(start + chunk_size, len(cell_indices))
        block = adata.X[cell_indices[start:stop], [pap_idx]]
        result[start:stop] = np.asarray(
            block.toarray() if hasattr(block, "toarray") else block
        ).ravel()
    return result


def build_donor_table(
    obs: pd.DataFrame,
    pappa2: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tal = obs.copy()
    tal["PAPPA2"] = pappa2
    tal = tal[
        tal["subclass.l1"].astype(str).eq("TAL")
        & tal["subclass.l3"].astype(str).isin(TAL_STATES)
        & tal["eGFR_rank"].notna()
    ].copy()

    donor_tal_count = (
        tal.groupby("donor_id", observed=True)
        .size()
        .rename("n_tal_cells")
        .reset_index()
    )
    donor_tal_count = donor_tal_count[
        donor_tal_count["n_tal_cells"] >= MIN_TAL_CELLS
    ]
    tal = tal[tal["donor_id"].isin(donor_tal_count["donor_id"])]

    donor_metadata = (
        tal[
            [
                "donor_id",
                "dataset",
                "condition.l1",
                "condition.l2",
                "eGFR",
                "eGFR_rank",
                "sex",
                "age_years",
            ]
        ]
        .drop_duplicates("donor_id")
        .merge(donor_tal_count, on="donor_id", how="left")
    )
    state_counts = (
        tal.groupby(["donor_id", "subclass.l3"], observed=True)
        .agg(
            n_cells=("PAPPA2", "size"),
            PAPPA2_mean=("PAPPA2", "mean"),
        )
        .reset_index()
    )
    state_counts["TAL_fraction"] = (
        state_counts["n_cells"]
        / state_counts["donor_id"].map(
            donor_tal_count.set_index("donor_id")["n_tal_cells"]
        )
    )
    state_counts.loc[
        state_counts["n_cells"] < MIN_STATE_CELLS,
        "PAPPA2_mean",
    ] = np.nan

    total_pappa2 = (
        tal.groupby("donor_id", observed=True)["PAPPA2"]
        .mean()
        .rename("PAPPA2_total_mean")
        .reset_index()
    )
    donor_metadata = donor_metadata.merge(
        total_pappa2,
        on="donor_id",
        how="left",
        validate="one_to_one",
    )
    return donor_metadata, state_counts


def test_association(
    donor_table: pd.DataFrame,
    value_col: str,
    label: str,
) -> dict[str, object]:
    frame = donor_table[["eGFR_rank", value_col]].dropna()
    if len(frame) < MIN_DONORS:
        return {
            "endpoint": label,
            "n_donors": int(len(frame)),
            "rho": np.nan,
            "p": np.nan,
        }
    result = spearmanr(
        frame["eGFR_rank"],
        frame[value_col],
        alternative="two-sided",
    )
    return {
        "endpoint": label,
        "n_donors": int(len(frame)),
        "rho": float(result.statistic),
        "p": float(result.pvalue),
    }


def run_cohort(
    donor_metadata: pd.DataFrame,
    state_counts: pd.DataFrame,
    *,
    cohort_name: str,
    donor_ids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = donor_metadata[
        donor_metadata["donor_id"].isin(donor_ids)
    ].copy()
    states = state_counts[
        state_counts["donor_id"].isin(donor_ids)
    ].copy()
    rows = []
    for state in TAL_STATES:
        state_frame = states[states["subclass.l3"].eq(state)].copy()
        fraction = metadata[["donor_id", "eGFR_rank"]].merge(
            state_frame[["donor_id", "TAL_fraction"]],
            on="donor_id",
            how="left",
        )
        fraction["TAL_fraction"] = fraction["TAL_fraction"].fillna(0.0)
        rows.append(
            {
                "cohort": cohort_name,
                "state": state,
                "endpoint": "TAL_fraction",
                **test_association(
                    fraction,
                    "TAL_fraction",
                    f"{state} TAL fraction",
                ),
            }
        )
        mean = metadata[["donor_id", "eGFR_rank"]].merge(
            state_frame[["donor_id", "PAPPA2_mean"]],
            on="donor_id",
            how="left",
        )
        rows.append(
            {
                "cohort": cohort_name,
                "state": state,
                "endpoint": "PAPPA2_mean",
                **test_association(
                    mean,
                    "PAPPA2_mean",
                    f"PAPPA2 in {state}",
                ),
            }
        )
    rows.append(
        {
            "cohort": cohort_name,
            "state": "all_TAL",
            "endpoint": "PAPPA2_mean",
            **test_association(
                metadata,
                "PAPPA2_total_mean",
                "PAPPA2 in all TAL",
            ),
        }
    )
    result = pd.DataFrame(rows)
    result["p_fdr"] = bh_adjust(result["p"])
    return result, metadata


def adjusted_models(
    donor_metadata: pd.DataFrame,
    state_counts: pd.DataFrame,
    *,
    donor_ids: set[str],
) -> pd.DataFrame:
    """Estimate condition-adjusted ordinal-eGFR associations.

    The coefficient is expressed per one-SD increase in the donor-level
    endpoint. These models are exploratory because eGFR category and
    diagnosis are strongly collinear.
    """
    metadata = donor_metadata[
        donor_metadata["donor_id"].isin(donor_ids)
    ].copy()
    states = state_counts[
        state_counts["donor_id"].isin(donor_ids)
    ].copy()
    rows = []
    for state in TAL_STATES:
        state_frame = states[states["subclass.l3"].eq(state)]
        fraction = metadata.merge(
            state_frame[["donor_id", "TAL_fraction"]],
            on="donor_id",
            how="left",
        )
        fraction["TAL_fraction"] = fraction["TAL_fraction"].fillna(0.0)
        candidates = [
            (
                "TAL_fraction",
                f"{state} TAL fraction",
                fraction,
            )
        ]
        mean = metadata.merge(
            state_frame[["donor_id", "PAPPA2_mean"]],
            on="donor_id",
            how="left",
        )
        candidates.append(
            (
                "PAPPA2_mean",
                f"PAPPA2 in {state}",
                mean,
            )
        )
        for value_col, label, frame in candidates:
            model_frame = frame[
                [
                    "eGFR_rank",
                    value_col,
                    "condition.l1",
                ]
            ].dropna().rename(columns={"condition.l1": "condition"})
            if (
                len(model_frame) < 12
                or model_frame["condition"].nunique() < 2
            ):
                rows.append(
                    {
                        "state": state,
                        "endpoint": label,
                        "n_donors": int(len(model_frame)),
                        "beta_per_sd": np.nan,
                        "p_condition_adjusted": np.nan,
                    }
                )
                continue
            scale = float(model_frame[value_col].std(ddof=0))
            if not np.isfinite(scale) or scale == 0:
                continue
            model_frame = model_frame.assign(
                value_z=(
                    model_frame[value_col]
                    - model_frame[value_col].mean()
                )
                / scale
            )
            fit = ols(
                "eGFR_rank ~ value_z + C(condition)",
                data=model_frame,
            ).fit(cov_type="HC3")
            rows.append(
                {
                    "state": state,
                    "endpoint": label,
                    "n_donors": int(len(model_frame)),
                    "beta_per_sd": float(fit.params["value_z"]),
                    "p_condition_adjusted": float(fit.pvalues["value_z"]),
                    "condition_levels": int(
                        model_frame["condition"].nunique()
                    ),
                }
            )

    all_tal = metadata[
        ["eGFR_rank", "PAPPA2_total_mean", "condition.l1"]
    ].dropna().rename(columns={"condition.l1": "condition"})
    if len(all_tal) >= 12 and all_tal["condition"].nunique() >= 2:
        scale = float(all_tal["PAPPA2_total_mean"].std(ddof=0))
        if np.isfinite(scale) and scale > 0:
            all_tal = all_tal.assign(
                value_z=(
                    all_tal["PAPPA2_total_mean"]
                    - all_tal["PAPPA2_total_mean"].mean()
                )
                / scale
            )
            fit = ols(
                "eGFR_rank ~ value_z + C(condition)",
                data=all_tal,
            ).fit(cov_type="HC3")
            rows.append(
                {
                    "state": "all_TAL",
                    "endpoint": "PAPPA2 in all TAL",
                    "n_donors": int(len(all_tal)),
                    "beta_per_sd": float(fit.params["value_z"]),
                    "p_condition_adjusted": float(fit.pvalues["value_z"]),
                    "condition_levels": int(
                        all_tal["condition"].nunique()
                    ),
                }
            )

    result = pd.DataFrame(rows)
    result["p_condition_adjusted_fdr"] = bh_adjust(
        result["p_condition_adjusted"]
    )
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    adata, obs, pap_idx = load_metadata()
    obs["dataset"] = obs["experiment"].astype(str).map(dataset_label)

    kpmp_donors = set(
        obs.loc[obs["dataset"].eq("KPMP"), "donor_id"].astype(str)
    )
    premiere = obs[
        obs["dataset"].eq("PREMIERE")
        & obs["subclass.l1"].astype(str).eq("TAL")
        & obs["eGFR_rank"].notna()
    ].copy()
    premiere_ids = set(premiere["donor_id"].astype(str))
    exclusive_ids = premiere_ids.difference(kpmp_donors)
    if not exclusive_ids:
        raise ValueError("No PREMIERE donors remain after excluding KPMP")

    cell_indices = np.flatnonzero(obs.index.isin(premiere.index))
    pappa2 = read_pappa2(adata, cell_indices, pap_idx)
    donor_table, state_counts = build_donor_table(
        obs.iloc[cell_indices].copy(),
        pappa2,
    )
    donor_table["premiere_exclusive"] = donor_table["donor_id"].isin(
        exclusive_ids
    )
    donor_table.to_csv(
        OUT / "PREMIERE_TAL_donor_state_table.csv",
        index=False,
    )
    state_counts.to_csv(
        OUT / "PREMIERE_TAL_state_counts.csv",
        index=False,
    )

    primary, primary_metadata = run_cohort(
        donor_table,
        state_counts,
        cohort_name="PREMIERE_exclusive",
        donor_ids=exclusive_ids,
    )
    sensitivity, sensitivity_metadata = run_cohort(
        donor_table,
        state_counts,
        cohort_name="PREMIERE_all",
        donor_ids=premiere_ids,
    )
    results = pd.concat([primary, sensitivity], ignore_index=True)
    results.to_csv(
        OUT / "exact_TAL_state_eGFR_associations.csv",
        index=False,
    )
    adjusted = adjusted_models(
        donor_table,
        state_counts,
        donor_ids=exclusive_ids,
    )
    adjusted.to_csv(
        OUT / "condition_adjusted_models.csv",
        index=False,
    )
    adata.file.close()

    report = [
        "# PREMIERE exact TAL-state and eGFR-category audit",
        "",
        "## Design",
        "",
        f"- PREMIERE TAL donors: {len(premiere_ids)}.",
        f"- PREMIERE donors not present in the KPMP label set: "
        f"{len(exclusive_ids)}.",
        "- eGFR is an ordered category: 20-29, 30-39, 40-49, 50-59, >60.",
        "- Exact TAL labels are atlas labels: "
        + ", ".join(TAL_STATES)
        + ".",
        "- The primary cohort is PREMIERE-exclusive; all PREMIERE is a "
        "sensitivity cohort.",
        "",
        "## Important boundary",
        "",
        "The low-eGFR samples are concentrated in CKD. The association is "
        "therefore a disease-linked ordinal association, not an independent "
        "continuous-eGFR validation and not evidence of causality.",
        "",
        "## Primary cohort",
        "",
        primary.sort_values(["p", "state"]).to_string(index=False),
        "",
        "## Condition-adjusted models",
        "",
        adjusted.sort_values(
            ["p_condition_adjusted", "state"]
        ).to_string(index=False),
        "",
        "Condition adjustment uses diagnosis category and HC3 robust "
        "standard errors. Because diagnosis and eGFR category are strongly "
        "collinear, the adjusted estimate is an exploratory within-diagnosis "
        "residual test, not a fully adjusted clinical model.",
        "",
        "## Sensitivity cohort",
        "",
        sensitivity.sort_values(["p", "state"]).to_string(index=False),
        "",
        "## Donor overlap",
        "",
        f"- PREMIERE/KPMP overlapping donor labels: "
        f"{len(premiere_ids.intersection(kpmp_donors))}.",
        f"- PREMIERE-exclusive donor labels: {len(exclusive_ids)}.",
        "",
    ]
    (OUT / "PREMIERE_exact_TAL_report.md").write_text(
        "\n".join(report),
        encoding="utf-8",
    )
    manifest = {
        "input": str(H5AD),
        "premiere_donors": len(premiere_ids),
        "premiere_exclusive_donors": len(exclusive_ids),
        "primary_donors": int(primary_metadata["donor_id"].nunique()),
        "sensitivity_donors": int(
            sensitivity_metadata["donor_id"].nunique()
        ),
        "state_labels": TAL_STATES,
        "egfr_mapping": EGFR_RANK,
        "primary_endpoint": "Spearman eGFR rank vs donor-level state fraction/PAPPA2",
    }
    (OUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(primary.sort_values(["p", "state"]).to_string(index=False))


if __name__ == "__main__":
    main()
