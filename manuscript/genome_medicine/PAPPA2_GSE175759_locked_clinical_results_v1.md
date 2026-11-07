# PAPPA2 GSE175759 locked clinical-association results v1

Date: 2026-11-08

## Status

Independent external clinical-association analysis in a primary
kidney-disease cohort. This is not a causal or clinical-utility validation.

Protocol:

```text
outputs\PAPPA2_GSE175759_locked_clinical_protocol_v1.md
```

Analysis code:

```text
work\run_gse175759_locked_clinical_validation.py
work\run_gse175759_posthoc_disease_sensitivity.py
```

## Cohort and quality control

GSE175759 contains 90 microdissected tubulointerstitial kidney samples.
Three technical outliers were excluded before scoring, leaving 87 samples.

- all 11 score genes detected in all 90 original samples;
- eGFR available for all 87 analysed samples;
- library size range: 22.3-58.4 million counts;
- no age or sex variables are publicly available.

Sample composition after exclusion:

| Diagnosis | n | Mean eGFR | Mean score |
|---|---:|---:|---:|
| Control | 22 | 95.8 | -0.61 |
| IgA nephropathy | 43 | 91.5 | 0.11 |
| Minimal change disease | 9 | 86.8 | 0.11 |
| Membranous nephropathy | 4 | 98.6 | 0.51 |
| Focal segmental glomerulosclerosis | 3 | 80.7 | 0.18 |
| Diabetic nephropathy | 3 | 80.4 | 1.58 |
| Lupus nephritis | 3 | 130.7 | 0.14 |

## Primary and prespecified models

All estimates use HC3 robust standard errors.

| Model | n | eGFR change per 1-SD score | 95% CI | P |
|---|---:|---:|---:|---:|
| Primary, unadjusted | 87 | -7.98 | -12.29 to -3.66 | 0.00050 |
| Diagnosis-adjusted | 87 | -8.71 | -13.79 to -3.63 | 0.0012 |
| Diagnosis and library-size adjusted | 87 | -8.56 | -13.43 to -3.70 | 0.00090 |

The direction is consistent with the prespecified expectation: higher
TAL-injury score, lower eGFR.

PAPPA2-only model:

```text
beta = -9.06 eGFR points per 1 SD
95% CI = -14.08 to -4.03
P = 0.00067
```

Correlations:

- Pearson r = -0.363, P = 0.00056;
- Spearman rho = -0.289, P = 0.0065.

## Post hoc disease-only sensitivity

The control group had a lower mean score than the disease groups, so
disease-only models were added after the primary analysis.

| Model | n | eGFR change per 1 SD | 95% CI | P |
|---|---:|---:|---:|---:|
| Disease samples only | 65 | -9.27 | -15.10 to -3.45 | 0.0027 |
| Disease, diagnosis-adjusted | 65 | -9.67 | -15.82 to -3.53 | 0.0031 |
| Disease, diagnosis and library adjusted | 65 | -9.13 | -15.12 to -3.14 | 0.0041 |
| IgA nephropathy only | 43 | -6.07 | -13.36 to 1.21 | 0.110 |

The association is therefore not generated solely by the control samples,
but it is not independently established within IgA nephropathy.

## Leave-one-gene-out

All 11 leave-one-gene-out scores retained a negative coefficient with a
95% confidence interval excluding zero. Estimates ranged from -7.28 to
-8.40 eGFR points per SD.

## Interpretation

GSE175759 provides independent, direction-consistent clinical-association
support in primary kidney disease. It is the strongest external clinical
evidence currently available for the frozen 11-gene score.

Required limitations:

- no age or sex adjustment;
- mixed disease diagnoses;
- observational cross-sectional eGFR;
- no causal or treatment-response inference;
- IgA-only analysis is underpowered;
- a second independent primary kidney-disease cohort with complete score-gene
  coverage was not identified.

Recommended manuscript wording:

> In an independent tubulointerstitial cohort of 87 primary kidney-disease
> samples, a higher 11-gene TAL-injury score was associated with lower eGFR,
> including after diagnosis and library-size adjustment. The absence of age
> and sex and the mixed diagnostic composition limit clinical interpretation.

Do not describe this as stable multi-cohort validation.

## Output files

```text
outputs\gse175759_locked_validation\GSE175759_gene_coverage.csv
outputs\gse175759_locked_validation\GSE175759_sample_scores.csv
outputs\gse175759_locked_validation\primary_and_sensitivity_models.csv
outputs\gse175759_locked_validation\posthoc_disease_subgroup_models.csv
outputs\gse175759_locked_validation\leave_one_gene_out.csv
outputs\gse175759_locked_validation\analysis_summary.json
```
