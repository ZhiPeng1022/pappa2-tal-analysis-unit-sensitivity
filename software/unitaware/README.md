# unitaware v0.2.0

`unitaware` is a small reference implementation for testing whether
disease-associated cell-state signals depend on the analysis unit.

It distinguishes:

1. donor pseudobulk expression;
2. neighborhood-expression association;
3. within-donor cell-rank association.

Neighborhood abundance is retained as a state-level unit and is not treated as
a gene-level effect.

## Input model

The API expects:

- an expression table indexed by `cell_id`, with genes as columns;
- a cell table with `donor_id`, `condition`, and `state`;
- an optional cell-level neighborhood exposure column;
- an optional neighborhood table with `neighborhood_id`, `state`,
  `comparison`, and `logfc`.

## Example

```python
import pandas as pd
from unitaware import analyze_units

results = analyze_units(
    expression=expression,
    cells=cells,
    case_levels=["Disease"],
    reference_level="Reference",
    exposure_col="exposure",
    neighborhood_logfc=neighborhood_logfc,
)
```

The returned object contains donor, neighborhood, and within-donor long tables
plus a merged reversal table.

## Effect-aware calibration

An expression-matched random set does not necessarily match the two effect-size
distributions being compared. The package therefore also provides a conditional
pairing permutation:

```python
from unitaware import paired_reversal_permutation

null_runs, summary = paired_reversal_permutation(
    effects,
    left_col="donor_effect_z",
    right_col="within_donor_effect_z",
    target_genes=target_genes,
    strata_cols=("state", "comparison"),
    z_threshold=0.5,
    n_permutations=10000,
    seed=20261028,
    aggregate="gene",
)
```

This test preserves each target gene's left-unit effect and the marginal
distribution of the right-unit effect within each stratum, while breaking the
gene-wise pairing between units. It should be reported alongside, not instead
of, conventional matched-set comparisons. `aggregate="gene"` first averages
valid strata within each gene and then averages genes; `aggregate="row"` gives
each finite gene-stratum pair equal weight.

## Command line

Generate a synthetic dataset:

```bash
unitaware simulate --seed 7 --out examples/simulated
```

Analyze an existing dataset:

```bash
unitaware analyze \
  --expression expression.csv \
  --cells cells.csv \
  --neighborhood-logfc neighborhood_logfc.csv \
  --case-level Disease \
  --reference-level Reference \
  --out results
```

The expression table must contain a `cell_id` column plus gene columns. The cell
table must contain `cell_id`, `donor_id`, `condition`, `state`, and an exposure
column. The neighborhood table must contain `neighborhood_id`, `state`,
`comparison`, and `logfc`.

## Scope

This is a reference framework, not a substitute for study-specific
normalization, covariate adjustment, or causal inference.

## License

The software is released under the MIT License. See `LICENSE`.
