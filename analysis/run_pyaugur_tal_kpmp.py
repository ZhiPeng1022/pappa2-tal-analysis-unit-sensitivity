import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
from pyaugur import calculate_auc


DEFAULT_INPUT_PATH = (
    Path(r"D:\CodexData\RAMP3\augur_input")
    / "KPMP_augur_input.h5ad"
)
DEFAULT_OUTPUT_DIR = Path(
    r"C:\Users\Elsa\Documents\Codex\2026-09-25\9-2\outputs"
)
DEFAULT_OUTPUT_PREFIX = "KPMP_Augur_TAL"
TAL_SUBCLASSES = [
    "C-TAL-A",
    "C-TAL-B",
    "frTAL",
    "aTAL2",
]
COMPARISONS = {
    "disease_vs_Reference": {
        "Reference": "Reference",
        "AKI": "Disease",
        "CKD": "Disease",
    },
    "AKI_vs_Reference": {
        "Reference": "Reference",
        "AKI": "AKI",
    },
    "CKD_vs_Reference": {
        "Reference": "Reference",
        "CKD": "CKD",
    },
}


def parse_list(value):
    if value is None:
        return []
    return [
        item.strip()
        for item in value.split(",")
        if item.strip()
    ]


def map_feature_names(frame, var_names):
    result = frame.copy()
    parsed = result["gene"].str.extract(r"^gene_(\d+)$")
    result["feature_idx"] = pd.to_numeric(
        parsed[0],
        errors="coerce",
    )
    result = result.dropna(subset=["feature_idx"])
    result["feature_idx"] = result["feature_idx"].astype(int)
    result["gene_symbol"] = result["feature_idx"].map(
        pd.Series(
            var_names.astype(str),
            index=range(len(var_names)),
        )
    )
    return result


def run_comparison(adata, label, mapping, args):
    subset = adata[
        adata.obs["condition"].isin(mapping.keys())
    ].copy()
    subset.obs["label"] = subset.obs["condition"].map(mapping)
    print(
        f"\n{label}: {subset.n_obs} cells",
        flush=True,
    )
    print(
        subset.obs.groupby(
            [args.cell_type_column, "label"],
            observed=True,
        )
        .size()
        .unstack(fill_value=0)
        .to_string(),
        flush=True,
    )
    result = calculate_auc(
        subset,
        label_col="label",
        cell_type_col=args.cell_type_column,
        n_subsamples=args.n_subsamples,
        subsample_size=args.subsample_size,
        folds=args.folds,
        min_cells=args.min_cells,
        var_quantile=args.var_quantile,
        feature_perc=args.feature_perc,
        n_threads=args.n_threads,
        classifier=args.classifier,
        seed=args.seed,
    )
    auc = result["AUC"].copy()
    auc["comparison"] = label
    auc["cells"] = subset.n_obs
    auc.to_csv(
        args.output_dir
        / f"{args.output_prefix}_{label}_AUC.csv",
        index=False,
        encoding="utf-8-sig",
    )
    result["results"].to_csv(
        args.output_dir
        / f"{args.output_prefix}_{label}_subsample_results.csv",
        index=False,
        encoding="utf-8-sig",
    )
    importance = map_feature_names(
        result["feature_importance"],
        subset.var_names,
    )
    importance["comparison"] = label
    importance.to_csv(
        args.output_dir
        / f"{args.output_prefix}_{label}_feature_importance.csv",
        index=False,
        encoding="utf-8-sig",
    )
    print(auc.to_string(index=False), flush=True)
    return auc, importance


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT_PATH,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
    )
    parser.add_argument(
        "--output-prefix",
        default=DEFAULT_OUTPUT_PREFIX,
    )
    parser.add_argument(
        "--cell-type-column",
        default="subclass",
    )
    parser.add_argument(
        "--cell-types",
        type=parse_list,
        default=",".join(TAL_SUBCLASSES),
    )
    parser.add_argument("--n-subsamples", type=int, default=20)
    parser.add_argument("--subsample-size", type=int, default=30)
    parser.add_argument("--folds", type=int, default=3)
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--var-quantile", type=float, default=0.5)
    parser.add_argument("--feature-perc", type=float, default=0.5)
    parser.add_argument("--n-threads", type=int, default=4)
    parser.add_argument("--classifier", default="rf")
    parser.add_argument("--seed", type=int, default=42)
    return parser


def main():
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.input)
    if args.cell_types:
        adata = adata[
            adata.obs[args.cell_type_column].isin(args.cell_types)
        ].copy()
    print(
        "TAL cells:",
        adata.n_obs,
        "cell types:",
        adata.obs[args.cell_type_column].nunique(),
        "genes:",
        adata.n_vars,
        flush=True,
    )

    auc_frames = []
    importance_frames = []
    for label, mapping in COMPARISONS.items():
        auc, importance = run_comparison(
            adata,
            label,
            mapping,
            args,
        )
        auc_frames.append(auc)
        importance_frames.append(importance)

    pd.concat(auc_frames, ignore_index=True).to_csv(
        args.output_dir / f"{args.output_prefix}_all_AUC.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pd.concat(importance_frames, ignore_index=True).to_csv(
        args.output_dir
        / f"{args.output_prefix}_feature_importance_mapped.csv",
        index=False,
        encoding="utf-8-sig",
    )


if __name__ == "__main__":
    main()
