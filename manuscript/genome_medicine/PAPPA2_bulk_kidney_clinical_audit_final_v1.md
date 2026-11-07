# PAPPA2 bulk kidney clinical-cohort audit final v1

Date: 2026-11-08

## Objective

Perform a final broad search for an independent bulk human kidney cohort with
sample-level clinical outcome data suitable for testing the frozen 11-gene
TAL-injury score.

## Search scope

GEO:

- 965 unique GSE records from kidney/renal clinical queries;
- 432 human kidney-related bulk series with at least 10 samples were scanned
  through the GEO sample-metadata text interface;
- all previously audited series were reused from local cache.

ArrayExpress/BioStudies:

- targeted searches for kidney/renal biopsy with eGFR, serum creatinine,
  fibrosis, progression, and human samples;
- GEO mirrors were not counted as independent cohorts;
- non-GEO SDRF files were checked when a relevant record was found.

## A-layer definition

For this audit, an A-layer series required:

1. primary human kidney-disease tissue, not blood, cancer, or a secondary
   liver-disease context;
2. at least approximately 20 samples;
3. sample-level quantitative eGFR or an equivalent prespecified clinical
   endpoint;
4. a public expression matrix;
5. complete or near-complete detection of the frozen 11-gene score.

## Result: one A-layer cohort

### GSE175759

- 90 tubulointerstitial kidney samples;
- 87 after excluding three submitter-marked technical outliers;
- 11/11 score genes detected;
- continuous eGFR available for all analysed samples;
- mixed primary kidney diseases, including IgA nephropathy and controls;
- no public age or sex variables.

Locked association result:

| Model | n | eGFR change per score SD | 95% CI | P |
|---|---:|---:|---:|---:|
| Unadjusted primary | 87 | -7.98 | -12.29 to -3.66 | 0.00050 |
| Diagnosis-adjusted | 87 | -8.71 | -13.79 to -3.63 | 0.0012 |
| Disease samples only | 65 | -9.27 | -15.10 to -3.45 | 0.0027 |

PAPPA2-only coefficient:

```text
-9.06 eGFR points per SD; 95% CI -14.08 to -4.03; P = 0.00067
```

All leave-one-gene-out scores retained a significant negative coefficient.

Detailed report:

```text
outputs\PAPPA2_GSE175759_locked_clinical_results_v1.md
```

## Best rejected candidates

### GSE166239

- 18 FFPE primary kidney biopsies: 6 controls, 6 hypertensive nephropathy,
  6 type 2 diabetic nephropathy;
- public eGFR, age, sex, diagnosis, and proteinuria;
- public logCPM matrix lacks `CRIP1` and `B3GALT5`;
- only 9/11 score genes are represented.

Decision:

```text
Rejected because the frozen 11-gene score cannot be reconstructed.
No eGFR association was run after observing incomplete coverage.
```

### GSE50892

- 17 kidney-biopsy samples from patients with liver disease;
- only 8 samples have GFR values;
- secondary liver-disease context rather than a primary CKD cohort.

Decision:

```text
Rejected: insufficient clinical data and unsuitable disease context.
```

### GSE287697

- 58 samples, but the public matrix is peripheral blood mononuclear cells;
- eGFR and transplant variables are available, but the tissue does not test
  kidney TAL expression.

Decision:

```text
Rejected for this score-validation question.
```

### GSE217427

- 44 tissue samples collapse to 16 donors;
- AKI versus no-AKI context;
- eGFR available;
- 11/11 genes detectable;
- donor-level unadjusted association: beta 0.84, P = 0.897;
- adjusted model: beta -6.61, P = 0.121.

Decision:

```text
Exploratory B-layer support only; null and underpowered.
```

### Other B-layer contexts

- GSE181757: 440 transplant biopsies with split-sample positive but
  inconsistent outcome association;
- GSE137570: small public eGFR/progression resource;
- GSE21374: 282 transplant biopsies with 51 failures and follow-up, but no
  age/sex and not yet used as the primary validation;
- GSE275824: 63 transplant samples from repeated time points and 22 patients;
- GSE290167: donor histology and donor age but no quantitative eGFR;
- GSE25902: repeated transplant histology and IFTA scores, not eGFR.

## ArrayExpress

No additional qualifying primary kidney-disease cohort was identified.

`E-GEOD-93798` states in the paper description that clinical creatinine and
eGFR data were used, but the public SDRF does not expose sample-level eGFR.
Its sample identifiers encode age and sex but not an unambiguous quantitative
eGFR. It was therefore not treated as an A-layer cohort.

`E-GEOD-69074` is a mouse kidney experiment and is not eligible.

## Final decision

The final audit found:

- one strong independent primary kidney-disease clinical-association cohort
  with complete score coverage: GSE175759;
- no second independent primary kidney-disease cohort with complete 11-gene
  coverage and usable sample-level eGFR.

Therefore:

> The study now has one independent external primary-kidney-disease
> clinical-association cohort, but it still does not have stable multicohort
> external clinical validation.

The appropriate manuscript claim is "independent external clinical-association
support", not "validated clinical biomarker" and not "stable multicohort
validation".

The final broad search is complete. Further cohort expansion is not justified
unless a new resource with unpublished clinical annotations becomes
available.
