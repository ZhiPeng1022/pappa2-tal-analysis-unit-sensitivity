# PAPPA2 GSE181757 expression-PC sensitivity analysis v1

Date: 2026-11-08

## Status

This analysis tests whether the positive association in the preselected
GSE181757 Validation stratum is explained by global expression structure.
It is a sensitivity analysis, not a new primary endpoint and not evidence of
causality.

## Analysis code

```text
work\run_gse181757_locked_validation.py
```

## Frozen score

The 11-gene score is unchanged:

```text
PAPPA2, TMPRSS4, LIX1, WNT10A, KRT7, CRIP1,
S100A14, LAMC2, LAMB3, SOD3, B3GALT5
```

The score is the equal-weight mean of within-Validation gene-wise z-scores.
No gene selection, removal, or reweighting used the outcome.

## Expression-PC construction

- Samples: 270 Validation samples.
- Probes: top 5,000 variable probes within Validation.
- Excluded probes: all 17 annotated probes for the 11 score genes.
- Transformation: probe-wise standardization across Validation samples.
- Components: first five principal components.
- Outcome use: none in probe filtering, standardization, or PCA.

The exclusion of score-gene probes prevents the adjustment from directly
removing the measured score signal. The resulting components represent
global expression structure outside the score genes.

## Primary result and sensitivity models

| Model | n | Progressors | OR per 1-SD score | 95% CI | P | AUC |
|---|---:|---:|---:|---:|---:|---:|
| Unadjusted primary | 270 | 43 | 1.63 | 1.20-2.22 | 0.0018 | 0.646 |
| Adjusted for PC1-PC3 | 270 | 43 | 1.76 | 1.14-2.71 | 0.0112 | 0.662 |
| Adjusted for PC1-PC5 | 270 | 43 | 1.67 | 1.07-2.63 | 0.0245 | 0.663 |

The first three components explained 48.1% of the variance in the selected
probes; the first five explained 54.3%.

## Score-to-PC correlations

| Component | Pearson r | P |
|---|---:|---:|
| PC1 | -0.635 | 7.13e-32 |
| PC2 | -0.049 | 0.423 |
| PC3 | -0.227 | 1.65e-4 |
| PC4 | -0.182 | 0.0027 |
| PC5 | -0.098 | 0.108 |

The score is strongly related to PC1. This indicates that the score captures
part of a broad transcriptomic axis rather than being statistically
independent of global expression structure.

## Interpretation

The Validation association was not abolished after adjustment for the first
three or five expression principal components. Therefore, the positive
Validation result is not explained by these leading axes of global expression
structure alone.

The analysis does not establish batch robustness because no public batch
labels, age, sex, fibrosis, or transplant-specific clinical covariates were
available. PCA adjustment can also remove biological variation related to
injury, so the adjusted estimates should not be interpreted as a purified
causal effect.

The Discovery stratum remained null in the primary analysis, and the
score-by-stratum interaction was P = 0.046. The evidence should therefore be
described as:

> a positive Validation association that persisted after leading
> expression-PC adjustment, but with unexplained cross-stratum heterogeneity.

It must not be described as stable external clinical validation.

## Output files

```text
outputs\gse181757_locked_validation\validation_expression_pc_sensitivity.csv
outputs\gse181757_locked_validation\validation_expression_pc_variance.csv
outputs\gse181757_locked_validation\validation_expression_pcs.csv
outputs\gse181757_locked_validation\validation_pca_top_probes.csv
outputs\gse181757_locked_validation\validation_score_pc_correlations.csv
```
