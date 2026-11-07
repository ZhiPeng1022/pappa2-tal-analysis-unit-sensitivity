# PAPPA2 GSE175759 score-specificity results v1

Date: 2026-11-08

## Status

Post hoc sensitivity analyses for the already completed GSE175759 clinical
association. These analyses were locked before execution and are not new
primary endpoints.

Protocol:

```text
outputs\PAPPA2_GSE175759_score_specificity_protocol_v1.md
```

Analysis code:

```text
work\run_gse175759_score_specificity.py
```

## Analysis 1: TAL-content adjustment

TAL-content markers:

```text
UMOD, SLC12A1, KCNJ1, CASR, CLDN16
```

The injury score and TAL-content score had Pearson r = 0.287
(P = 0.0071).

| Model | Injury-score effect | 95% CI | P |
|---|---:|---:|---:|
| Injury only | -7.98 | -12.29 to -3.66 | 0.00050 |
| TAL-content adjusted | -8.81 | -13.76 to -3.85 | 0.00079 |
| TAL-content and diagnosis adjusted | -10.13 | -15.70 to -4.56 | 0.00062 |

The TAL-content term itself was not significant in the adjusted model
(2.90; 95% CI -1.44 to 7.23; P = 0.194). The injury-score association is
therefore not explained by the fixed TAL-content score.

## Analysis 2: expression-stratified matched gene sets

The null used 10,000 matched 11-gene sets drawn from 14,876 eligible control
genes after 20-bin expression stratification.

Observed injury-score coefficient:

```text
-7.98 eGFR points per SD
```

Null coefficient distribution:

- mean: -3.25;
- SD: 1.28;
- 2.5th percentile: -5.58;
- 97.5th percentile: -0.54.

The observed score was at the 0.01 percentile of the null distribution
(one-sided empirical P = 0.0002). The observed association is therefore
markedly stronger than expression-matched random 11-gene sets.

## Analysis 3: PAPPA2 versus the remaining score

| Model | Term | Effect | 95% CI | P |
|---|---|---:|---:|---:|
| PAPPA2 only | PAPPA2 | -9.04 | -14.08 to -4.00 | 0.00071 |
| PAPPA2 + remaining 10 genes | PAPPA2 | -7.26 | -12.31 to -2.21 | 0.0060 |
| PAPPA2 + remaining 10 genes | Remaining 10-gene score | -4.10 | -8.11 to -0.10 | 0.0478 |

The remaining 10-gene component increased R-squared from 0.1689 to 0.1972
(delta R-squared = 0.0283). The nested F test was not significant
(P = 0.089).

Interpretation:

```text
The residual 10-gene signal is possible but not robustly established.
The data do not justify a claim that the composite score clearly outperforms
PAPPA2 alone.
```

## Overall interpretation

The primary clinical association survives TAL-content adjustment and performs
better than expression-matched random 11-gene sets. This materially strengthens
the specificity of the clinical association. The incremental value of the
non-PAPPA2 component remains borderline and should be described cautiously.

## Output files

```text
outputs\gse175759_locked_validation\GSE175759_full_gene_logcpm.csv.gz
outputs\gse175759_locked_validation\tal_content_models.csv
outputs\gse175759_locked_validation\matched_set_null_summary.csv
outputs\gse175759_locked_validation\matched_set_null_coefficients.csv
outputs\gse175759_locked_validation\pappa2_incremental_models.csv
outputs\gse175759_locked_validation\pappa2_incremental_comparison.csv
outputs\gse175759_locked_validation\score_specificity_summary.json
```
