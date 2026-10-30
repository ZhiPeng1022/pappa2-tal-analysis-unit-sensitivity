# PAPPA2 Genome Medicine upgrade plan v2

Date: 2026-10-29

## Current position

The manuscript remains a high-risk Genome Medicine Research Article, but the
statistical argument is now more defensible. The paper is no longer presented
as evidence of broad panel-level negative coupling. It is presented as:

> a reusable framework showing that analysis units can change single-cell
> disease conclusions, with explicit expression-matched and effect-aware
> calibration.

## Completed in this round

1. Added effect-aware conditional pairing permutation to unitaware v0.2.0.
2. Expanded unit tests from one to six and added an MIT license.
3. Built and verified the Python wheel.
4. Replaced the unsupported panel-enrichment claim in manuscript v4.
5. Updated cover letter, pre-submission inquiry, and APC request to v2.
6. Added Additional files 28 and 29.
7. Regenerated Supplementary Figure S1 at 600 dpi.
8. Updated the manuscript to a 36-page checked DOCX/PDF.
9. Assembled draft package v4 with 69 files and verified ZIP integrity.
10. Searched for new non-KPMP cohorts and identified GSE284715 as the best
    near-term candidate.

## Scientific judgment

Strengths:

- complete KPMP TAL resource;
- multiple analysis units in the same cells;
- case gene plus panel-level framework;
- two existing partial external cohorts;
- reproducible software and effect-aware calibration;
- explicit negative findings and boundaries.

Remaining high-risk gaps:

- no exact aTAL2 external replication;
- no experimental or functional validation;
- no strong external clinical outcome cohort;
- GSE284715 matrix was not downloaded;
- the effect-aware null removes the broad donor-within coupling claim.

## Decision gate

The next decision should be based on editorial feedback, not another blind
journal downgrade.

If Genome Medicine says the work is outside scope or the APC waiver is not
viable:

- use the same corrected package for Human Genetics;
- retain the effect-aware null and negative findings;
- do not return to the old BMC Medical Genomics claim set.

If Genome Medicine is interested:

- finish GSE284715 through the resumable parallel downloader;
- the corrected archive has been released;
- the new DOI has been inserted
  (`https://doi.org/10.5281/zenodo.23811766`);
- submit.

## APC

Current list APC: GBP 3,690 / USD 5,490 / EUR 4,490.

Automatic country-based waiver does not apply to China. A discretionary
waiver or discount must be requested at submission. The package includes a
v2 request letter.
