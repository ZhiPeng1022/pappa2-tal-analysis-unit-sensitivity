# Release v1.0.4

Date: 2026-10-28

## Changes

- Added independent non-KPMP analyses of GSE211785 and GSE254185.
- Added donor-pseudobulk, neighborhood-feature, neighborhood-expression, and
  within-donor rank scripts for the two validation cohorts.
- Added `Additional file 14.csv` with validation effects, permutation p values,
  FDRs, and donor counts.
- Added curated validation result tables under
  `data/derived/independent_validation/`.
- Updated the README, additional-file index, `CITATION.cff`, and
  `.zenodo.json` to version 1.0.4.

## Main validation findings

- GSE211785 reproduced donor-level PAPPA2 disease elevation in C_TAL
  (+0.283; permutation p = 0.0315; six-test FDR = 0.063) and M_TAL
  (+0.200; p = 0.105).
- GSE254185 showed positive donor-level and neighborhood-abundance
  directions but negative neighborhood-expression and within-donor
  directions in TAL Injured.
- The exact aTAL2/frTAL discovery labels were not available in either
  independent cohort. These analyses are independent-cohort validation,
  not experimental or causal validation.

## Archive note

Publishing this release through Zenodo should create a new exact-version DOI.
The manuscript and submission package must be updated to that DOI before
submission.
