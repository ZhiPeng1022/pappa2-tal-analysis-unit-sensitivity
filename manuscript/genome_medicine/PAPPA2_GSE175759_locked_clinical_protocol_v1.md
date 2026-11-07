# PAPPA2 GSE175759 locked clinical-association protocol v1

Date locked: 2026-11-08

## Purpose

Test the frozen 11-gene TAL-injury score against continuous eGFR in an
independent primary kidney-disease cohort. This is an observational
clinical-association analysis, not a causal or clinical-utility validation.

## Data

GEO series:

```text
GSE175759
```

Tissue:

```text
manually microdissected tubulointerstitial kidney tissue
```

The series contains 90 samples:

- 46 IgA nephropathy;
- 22 living-donor controls;
- 9 minimal change disease;
- 4 membranous nephropathy;
- 3 focal segmental glomerulosclerosis;
- 3 diabetic nephropathy;
- 3 lupus nephritis.

Three samples are marked by the submitters as:

```text
Outlier_excluded from main analysis
```

They will be excluded before scoring and modeling.

## Frozen score

The score and weights are unchanged:

```text
PAPPA2, TMPRSS4, LIX1, WNT10A, KRT7, CRIP1,
S100A14, LAMC2, LAMB3, SOD3, B3GALT5
```

The Ensembl identifiers are fixed from the GSE217427 annotation:

```text
PAPPA2   ENSG00000116183
TMPRSS4  ENSG00000137648
LIX1     ENSG00000145721
WNT10A   ENSG00000135925
KRT7     ENSG00000135480
CRIP1    ENSG00000213145
S100A14  ENSG00000189334
LAMC2    ENSG00000058085
LAMB3    ENSG00000196878
SOD3     ENSG00000109610
B3GALT5  ENSG00000183778
```

All 11 genes were detected in all 90 samples before exclusions.

## Processing

- Counts are converted to counts per million using the total counts across all
  genes in each sample.
- Values are transformed as `log2(CPM + 1)`.
- Each gene is z-scored within the analysis cohort after outlier exclusion.
- The score is the equal-weight mean of the 11 gene-wise z-scores.
- The composite score is standardized once more so the coefficient is per
  one-standard-deviation increase.

No outcome information is used for gene selection, gene removal, probe choice,
or weighting.

## Primary endpoint

```text
eGFR ~ score_sd
```

Estimated by ordinary least squares with HC3 heteroskedasticity-robust
standard errors.

Expected direction:

```text
higher injury score -> lower eGFR
```

Report:

- eGFR change per one-SD score;
- 95% confidence interval;
- two-sided p value;
- n and model R-squared.

## Pre-specified secondary and sensitivity analyses

- `eGFR ~ score_sd + C(diagnosis)`;
- `eGFR ~ score_sd + C(diagnosis) + log_library_size`;
- PAPPA2-only model using the same primary structure;
- leave-one-gene-out score reconstruction;
- Pearson and Spearman correlations.

These analyses are descriptive and will not be used to select a preferred gene
set or model after seeing results.

## Multiplicity

There is one primary endpoint and one primary model. Secondary and sensitivity
models will be labeled as such and will not receive confirmatory claims.

## Interpretation rules

- A negative association with a confidence interval excluding zero is
  independent external clinical-association support.
- A null or positive association must be reported as a failure to transfer the
  expected clinical direction.
- The absence of age and sex is a material limitation.
- The cohort combines several kidney diseases; diagnosis adjustment is
  mandatory for interpretation.
- Transplant cohorts remain separate and cannot be pooled with this primary
  kidney-disease cohort.
- No result from this cohort establishes causality or clinical utility.
