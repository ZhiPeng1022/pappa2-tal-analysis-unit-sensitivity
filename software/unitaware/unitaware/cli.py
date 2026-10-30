"""Command-line entry point for unitaware."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .core import analyze_units
from .simulation import simulate_reversal


def read_table(path: str) -> pd.DataFrame:
    source = Path(path)
    separator = "\t" if source.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(source, sep=separator)


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="unitaware",
        description="Unit-aware single-cell sensitivity analysis",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    simulate = subparsers.add_parser(
        "simulate",
        help="Generate a synthetic unit-reversal dataset",
    )
    simulate.add_argument("--seed", type=int, default=20261028)
    simulate.add_argument("--out", type=str, default="")
    analyze = subparsers.add_parser(
        "analyze",
        help="Run donor, neighborhood, and within-donor unit analysis",
    )
    analyze.add_argument("--expression", required=True)
    analyze.add_argument("--cells", required=True)
    analyze.add_argument("--neighborhood-logfc", required=True)
    analyze.add_argument("--out", required=True)
    analyze.add_argument("--case-level", action="append", required=True)
    analyze.add_argument("--reference-level", default="Reference")
    analyze.add_argument("--state-col", default="state")
    analyze.add_argument("--donor-col", default="donor_id")
    analyze.add_argument("--condition-col", default="condition")
    analyze.add_argument("--neighborhood-col", default="neighborhood_id")
    analyze.add_argument("--exposure-col", default="exposure")
    args = parser.parse_args()

    if args.command == "simulate":
        expression, cells, neighborhood = simulate_reversal(
            seed=args.seed
        )
        if args.out:
            output_dir = Path(args.out)
            output_dir.mkdir(parents=True, exist_ok=True)
            expression.to_csv(output_dir / "expression.csv")
            cells.to_csv(output_dir / "cells.csv")
            neighborhood.to_csv(
                output_dir / "neighborhood_logfc.csv",
                index=False,
            )
        print(f"cells={len(cells)}")
        print(f"genes={expression.shape[1]}")
        print(f"neighborhoods={len(neighborhood)}")
    elif args.command == "analyze":
        expression = read_table(args.expression).set_index("cell_id")
        cells = read_table(args.cells)
        neighborhood = read_table(args.neighborhood_logfc)
        result = analyze_units(
            expression,
            cells,
            case_levels=args.case_level,
            reference_level=args.reference_level,
            neighborhood_logfc=neighborhood,
            exposure_col=args.exposure_col,
            state_col=args.state_col,
            donor_col=args.donor_col,
            condition_col=args.condition_col,
            neighborhood_col=args.neighborhood_col,
        )
        output_dir = Path(args.out)
        output_dir.mkdir(parents=True, exist_ok=True)
        result.donor.to_csv(
            output_dir / "donor_pseudobulk.csv",
            index=False,
        )
        result.neighborhood.to_csv(
            output_dir / "neighborhood_expression.csv",
            index=False,
        )
        result.within_donor.to_csv(
            output_dir / "within_donor.csv",
            index=False,
        )
        result.reversals.to_csv(
            output_dir / "unit_reversals.csv",
            index=False,
        )
        print(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
