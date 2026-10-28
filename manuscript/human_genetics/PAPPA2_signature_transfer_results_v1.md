# PAPPA2 TAL signature-transfer results v1

Date: 2026-10-28

## Purpose

Test whether the four KPMP TAL subtype signatures can identify an
`aTAL2`-like state in two non-KPMP kidney resources, and whether PAPPA2
shows the same donor-level positive / within-donor negative direction
contrast.

This is signature transfer, not exact-label replication.

This analysis was completed after the v1.0.4 archive. If it is used in the
Human Genetics Brief Report, it must be included in a subsequent release
(v1.0.5) with a new exact-version DOI; do not describe v1.0.4 as containing
the signature-transfer analysis.

## Discovery and held-out validation

Data:

- KPMP full snRNA-seq: `D:\CodexData\RAMP3\KPMP_new_snRNA.h5ad`
- Four labels: `C-TAL-A`, `C-TAL-B`, `aTAL2`, `frTAL`
- 106,851 TAL cells
- Donor-stratified split:
  - 156 training donors
  - 67 held-out test donors
- Top 35 markers per subtype selected from training donors
- Marker score: mean within-subtype marker z score
- High-confidence threshold: `marker_margin >= 0.2`

Held-out prediction:

| Metric | Result |
|---|---:|
| Test cells | 23,152 |
| Coverage at margin >= 0.2 | 86.66% |
| Accuracy among assigned cells | 93.90% |
| Macro F1 among assigned cells | 93.17% |
| aTAL2 precision | 83.28% |
| aTAL2 recall | 93.93% |

The marker panel transfers well enough within KPMP to justify a cautious
external projection.

Selected aTAL2 markers in the training split include:

```text
SPP1, KCNIP4, HIF1A, ADAMTS1, NRP1, TFPI, KAZN, NLGN1,
SOD2, NEAT1, PRKCE, C4orf19, ITGB6, SASH1
```

These markers should be treated as a signature of an inflammatory/injured
TAL-like state, not as a complete aTAL2 identity.

## External projection

### GSE254185

High-confidence aTAL2-like cells:

- 1,176 cells total at `marker_margin >= 0.2`
- 377 at margin 0; 1,116 at margin 0.1; 1,176 at margin 0.2
- Mostly derived from `TAL Inflammatory` and `TAL Injured`

Backbone-state composition at margin >= 0.2:

| Backbone state | Normal | UUO |
|---|---:|---:|
| TAL Inflammatory | 15 | 550 |
| TAL Injured | 4 | 333 |
| mTAL | 62 | 166 |
| cTAL2 | 2 | 22 |
| cTAL1 | 6 | 13 |

PAPPA2:

| Analysis | Effect | p | Interpretable donors |
|---|---:|---:|---|
| donor pseudobulk, UUO vs Normal | +0.554 | 0.0043 | 5 vs 6 |
| within-donor rank vs injury score | -0.581 | 0.125 | 4 vs 2 |
| within-donor rank vs neighborhood logFC | -0.365 | 0.125 | 4 vs 2 |

The donor-level increase is significant in the permutation test. The
within-donor negative direction remains directionally consistent but is
underpowered after requiring high-confidence aTAL2-like cells.

### GSE211785

High-confidence aTAL2-like cells:

- 3,751 cells at `marker_margin >= 0.2`
- Mostly mapped from `M_TAL`, with a smaller `C_TAL` component

Backbone-state composition at margin >= 0.2:

| Backbone state | Control | DKD | HKD |
|---|---:|---:|---:|
| M_TAL | 18 | 16 | 2,620 |
| C_TAL | 100 | 111 | 886 |

PAPPA2:

| Analysis | Effect | p | Interpretable donors |
|---|---:|---:|---|
| donor pseudobulk, disease vs control | +0.570 | 0.0235 | 10 vs 6 |
| within-donor rank vs neighborhood logFC | -0.095 | 0.405 | 8 vs 4 |

The donor-level PAPPA2 elevation is significant. The within-donor
direction is negative but not significant.

## Interpretation

The signature transfer supports three claims:

1. The four KPMP TAL signatures can be projected to independent kidney
   datasets with high held-out accuracy within KPMP.
2. Predicted aTAL2-like cells show higher PAPPA2 at the donor level in both
   external cohorts.
3. The within-donor direction is negative in both cohorts, producing the
   same qualitative direction conflict as the discovery analysis.

It does not support these claims:

1. The external predictions are not exact `aTAL2` labels.
2. The high-confidence GSE254185 aTAL2-like cells are predominantly
   inflammatory/injured TAL cells, not a clean cortical aTAL2 population.
3. The within-donor estimates are underpowered and not statistically
   significant at the high-confidence threshold.

## Decision

For a **Human Genetics Brief Report**, this is probably enough to prepare a
focused submission:

- one main finding: donor-level PAPPA2 elevation in an aTAL2-like signature
  state;
- one secondary finding: a qualitative donor/within-donor direction conflict;
- held-out KPMP validation and two external cohorts in one figure/table.

For a full Human Genetics Research Article, it is not yet enough. The next
limited addition would be kidney eQTL/ATAC/GWAS enrichment, but only with a
pre-specified analysis and a clear result. Do not add a token genetic
analysis solely to make the paper look more genetic.

## Files

```text
work\signature_transfer\
  KPMP_TAL_markers.csv
  KPMP_internal_validation.csv
  KPMP_internal_validation_predictions.csv
  GSE254185_signature_predictions.csv
  GSE254185_signature_summary.csv
  GSE211785_signature_predictions.csv
  GSE211785_signature_summary.csv
  signature_transfer_validation_summary.csv

work\signature_transfer.py
work\analyze_signature_transfer_validation.py
```
