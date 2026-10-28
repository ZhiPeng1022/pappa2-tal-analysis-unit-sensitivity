# v1.0.7 - prespecified analysis-unit panel

Date: 2026-10-28

This release adds a prospectively frozen 400-gene analysis-unit panel and a
2,000-gene expression- and detection-matched control pool.

## Added

- Frozen 400-gene target panel, 100 genes per biological or control group.
- 2,000-gene expression-matched external control pool.
- 1,000 random 400-gene and 300-gene null sets.
- Expression-stratified random-set sensitivity analysis.
- Gene-level metrics across donor, neighborhood, and within-donor units.
- Supplementary Tables S21-S23.
- Supplementary Fig. S1.

## Primary result

Donor-versus-within reversal:

```text
target panel:               18.6%
expression-matched random:  12.7%
empirical P:                0.000999
```

Biological 300-gene subset:

```text
target:                     20.7%
random:                     12.7%
empirical P:                0.000999
```

Internal negative controls matched the random control distribution
(12.3% versus 12.6%).

## Scope

The result is a panel-level enrichment. It does not imply that every selected
gene is stable or causal.
