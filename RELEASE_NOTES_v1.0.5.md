# v1.0.5 - Human Genetics Research Article candidate

Date: 2026-10-28

This release adds the Human Genetics Research Article working package and the
new analyses required to support it.

## Added

- Human Genetics Research Article manuscript in Markdown, DOCX, and PDF.
- Human Genetics cover letter in Markdown, DOCX, and PDF.
- Eight submission TIFF figures.
- Supplementary Tables S1-S20.
- Donor-stratified signature transfer across KPMP, GSE211785, and GSE254185.
- NephQTL2 tubulointerstitial and glomerular eQTL results for PAPPA2.
- GSE211785 snATAC gene-activity boundary analysis.
- Donor-level composition and expression decomposition.
- Manuscript consistency, DOCX layout, and PDF layout checks.

## Key results

- Signature transfer achieved 93.9% accuracy among assigned held-out cells and
  93.9% aTAL2 recall.
- External cohorts supported donor-level PAPPA2 elevation, but within-donor
  negative estimates were underpowered.
- PAPPA2 was a strong tubulointerstitial eGene (FDR = 1.64e-14) without
  glomerular eQTL support.
- PAPPA2 was not among FDR-significant PTWAS genes, and same-cohort snATAC did
  not show higher disease accessibility.
- In the combined disease comparison, the observed four-subtype TAL PAPPA2
  difference was +0.563, with an expression component of +0.641, a composition
  component of +0.010, and an interaction/covariance term of -0.088.

## Scientific boundary

The release does not claim exact external aTAL2 replication, genetic
colocalization, causality, or a validated biomarker.

## Archive status

This is the source release from which a new exact-version Zenodo DOI will be
generated. The final manuscript must insert that DOI after archiving and be
re-exported before journal submission.
