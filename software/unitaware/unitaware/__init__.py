"""Unit-aware analysis of single-cell disease associations."""

from .core import (
    analyze_units,
    decompose_composition,
    donor_pseudobulk,
    neighborhood_expression,
    paired_reversal_permutation,
    reversal_fraction,
    standardize_and_score,
    within_donor,
)
from .simulation import simulate_reversal

__all__ = [
    "analyze_units",
    "decompose_composition",
    "donor_pseudobulk",
    "neighborhood_expression",
    "paired_reversal_permutation",
    "reversal_fraction",
    "simulate_reversal",
    "standardize_and_score",
    "within_donor",
]
