# PAPPA2 GSE175759 score-specificity protocol v1

Date locked: 2026-11-08

## Status

These are post hoc sensitivity analyses for the already completed GSE175759
clinical association. They are designed to test alternative explanations for
the observed score-eGFR relationship. They are not new primary endpoints.

## Cohort and base score

Use the locked GSE175759 analysis cohort:

- 90 samples before exclusions;
- three submitter-marked technical outliers excluded;
- 87 samples for analysis;
- frozen 11-gene TAL-injury score;
- continuous eGFR.

The full gene-level matrix will be reconstructed from the 90 public per-sample
count files using the GSE217427 Ensembl-to-symbol annotation. Multiple Ensembl
identifiers mapping to the same symbol will be summed. Expression will be
converted to log2(CPM + 1).

## Analysis 1: TAL-content adjustment

A fixed TAL identity score will be built from:

```text
UMOD, SLC12A1, KCNJ1, CASR, CLDN16
```

Each marker will be z-scored across the 87 samples on log2(CPM + 1), and the
TAL-content score will be the equal-weight mean. The injury and TAL-content
scores will then be standardized.

Models:

```text
eGFR ~ injury_score
eGFR ~ injury_score + TAL_content_score
eGFR ~ injury_score + TAL_content_score + C(diagnosis)
```

All models will use HC3 robust standard errors. The Pearson correlation between
the two scores will be reported.

## Analysis 2: expression-stratified matched gene sets

Eligible control genes will:

- map to a unique gene symbol;
- be detected above zero in at least 80% of the 87 samples;
- have mean log2(CPM + 1) of at least 1;
- exclude the 11 injury-score genes and the five TAL-content markers.

Eligible genes will be divided into 20 equal-frequency bins by mean expression
across the 87 samples. For each of the 11 score genes, one control gene will be
drawn from the same expression bin without replacement. Each matched 11-gene
set will be processed exactly like the observed score: gene-wise z score,
equal-weight mean, standardization, and OLS of eGFR on score.

Number of matched sets:

```text
10,000
```

The observed coefficient will be compared with the one-sided negative tail of
the null distribution. The empirical P value will be:

```text
(1 + number of null coefficients <= observed coefficient) / (10,000 + 1)
```

The observed percentile within the null distribution will also be reported.

## Analysis 3: PAPPA2 versus the remainder of the score

Define:

```text
PAPPA2_sd = standardized PAPPA2 expression
score_without_PAPPA2 = standardized mean of the other 10 gene-wise z scores
```

Models:

```text
eGFR ~ PAPPA2_sd
eGFR ~ PAPPA2_sd + score_without_PAPPA2
```

Both models will use HC3 robust standard errors. A nested F test will compare
the two models. The analysis asks whether the remaining 10-gene component adds
information beyond PAPPA2 alone.

## Interpretation rules

- Failure to survive TAL-content adjustment weakens a TAL-injury-specific
  interpretation.
- A matched-set percentile near the middle of the null distribution means the
  observed score is not distinguishable from expression-matched gene sets.
- A non-significant incremental component means the composite score should not
  be presented as adding information beyond PAPPA2 alone.
- All results, including negative or conflicting findings, will be reported.
