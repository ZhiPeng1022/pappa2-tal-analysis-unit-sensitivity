# Release v1.0.11

Date: 2026-11-08

## GSE175759 clinical-association analysis

This release adds the independent GSE175759 primary kidney-disease
clinical-association analysis to the Genome Medicine manuscript.

The locked analysis used 87 samples after excluding three
submitter-marked technical outliers. The frozen 11-gene TAL-injury score was
associated with lower eGFR:

```text
-7.98 eGFR points per 1-SD score
95% CI -12.29 to -3.66
P = 0.00050
```

Diagnosis-adjusted:

```text
-8.71 eGFR points per SD
95% CI -13.79 to -3.63
P = 0.0012
```

The disease-only estimate was -9.27 (P = 0.0027). All leave-one-gene-out
scores retained a negative coefficient with confidence intervals excluding
zero.

## Score-specificity analyses

TAL-content-adjusted estimate:

```text
-8.81 eGFR points per SD
95% CI -13.76 to -3.85
P = 0.00079
```

Against 10,000 expression-stratified matched 11-gene sets, the observed score
was at percentile 0.01 (one-sided empirical P = 0.0002; null mean -3.25).

PAPPA2 alone remained strongly associated with eGFR. The remaining 10-gene
component had a borderline incremental coefficient (-4.10; robust P = 0.0478),
but the nested F test was not significant (P = 0.089). The manuscript does not
claim that the composite score clearly outperforms PAPPA2 alone.

## Added files

- Genome Medicine manuscript v7;
- Additional file 31;
- GSE175759 full gene-level logCPM and derived tables;
- clinical-cohort audit;
- GSE175759 analysis scripts;
- cover letter and pre-submission inquiry v4;
- APC waiver request v3.

## Boundary

- One independent primary kidney-disease clinical-association cohort;
- no stable multicohort clinical validation;
- no clinical utility or causal claim;
- age and sex unavailable in GSE175759;
- APC matters are separate from editorial correspondence.

The exact Zenodo DOI for this release is generated asynchronously after the
GitHub release is published.
