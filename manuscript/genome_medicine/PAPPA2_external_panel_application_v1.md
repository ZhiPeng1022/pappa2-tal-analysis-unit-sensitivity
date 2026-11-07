# PAPPA2 external panel application v1

Date: 2026-11-01

## Purpose

Test whether the frozen 400-gene analysis-unit panel shows the same
expression-matched versus effect-aware calibration pattern outside KPMP.

## Data

Local processed data from two non-KPMP human kidney cohorts were used. No
additional download was required.

- GSE211785: 29,847 C_TAL and M_TAL cells from 35 donors.
- GSE254185: 9,912 TAL cells from 12 kidneys.

The 400-gene target panel and 2,000-gene external control pool were mapped to
the input gene identifiers:

- 385 target genes were available in both cohorts;
- 1,915 external control genes were available in both cohorts;
- 2,300 genes were analyzed in total.

## Analysis

For each cohort, state, and disease comparison, the analysis estimated:

1. donor pseudobulk difference between disease and Reference donors;
2. within-donor association between gene expression and local neighborhood
   disease logFC.

Both effects were standardized across the 2,300 available genes. A reversal
required opposite signs and absolute standardized effects of at least 0.5 in
both units. Donor pseudobulk required at least five cells per donor, and
within-donor analysis required at least 20 cells per donor. A state-comparison
group required at least three disease and three Reference donors.

Two controls were used:

- expression-selected random sets of 385 genes drawn from the 1,915 external
  control genes, with 1,000 draws;
- a post hoc effect-aware pairing null with 10,000 permutations that preserved
  each gene's left-unit effect and the marginal distribution of the right-unit
  effect within each state and comparison.

The external neighborhood exposure is based on a k-nearest-neighbor graph built
from public UMAP coordinates. It is a transparent approximation of the KPMP
Milo exposure, not an exact reproduction of Milo.

## Results

### GSE211785

- cells: 29,847;
- estimable state-comparison groups: 6;
- target panel: 19.49%;
- external control genes: 14.28%;
- expression-matched random sets: 14.27%, empirical P = 0.000999;
- effect-aware conditional null: 19.20%, empirical P = 0.347.

Panel-group mean reversal fractions:

| Group | Reversal fraction |
|---|---:|
| TAL subtype markers | 24.55% |
| Injury, inflammation, and metabolism | 19.70% |
| PAPPA2 coexpression | 17.02% |
| Internal negative controls | 16.84% |
| External control pool | 14.28% |

### GSE254185

- cells: 9,912;
- estimable state-comparison groups: 2 (`cTAL2` and `TAL Injured`);
- target panel: 18.75%;
- external control genes: 16.06%;
- expression-matched random sets: 16.07%, empirical P = 0.0120;
- effect-aware conditional null: 18.65%, empirical P = 0.484.

Panel-group mean reversal fractions:

| Group | Reversal fraction |
|---|---:|
| PAPPA2 coexpression | 30.32% |
| TAL subtype markers | 21.51% |
| Internal negative controls | 16.49% |
| External control pool | 16.06% |
| Injury, inflammation, and metabolism | 7.50% |

## Interpretation

Both external cohorts reproduced the calibration pattern observed in KPMP:

- expression-matched gene sets suggested that the target panel was enriched
  for cross-unit reversal;
- conditional nulls that preserved the observed effect-size marginals removed
  that apparent enrichment.

The result supports the general framework claim that expression matching alone
does not distinguish biological coupling from effect-size selection. It does
not show that the 400-gene panel is universally enriched for analysis-unit
reversal.

The two external cohorts differ substantially in disease context, cell number,
and state coverage. GSE254185 contributed only two estimable
state-comparison groups, so its result is supportive rather than definitive.

## Output files

```text
outputs\external_panel_application\gse211785\external_panel_result.json
outputs\external_panel_application\gse211785\gene_unit_metrics_long.csv
outputs\external_panel_application\gse211785\gene_unit_summary.csv
outputs\external_panel_application\gse211785\random_control_sets.csv
outputs\external_panel_application\gse211785\conditional_null.csv

outputs\external_panel_application\gse254185\external_panel_result.json
outputs\external_panel_application\gse254185\gene_unit_metrics_long.csv
outputs\external_panel_application\gse254185\gene_unit_summary.csv
outputs\external_panel_application\gse254185\random_control_sets.csv
outputs\external_panel_application\gse254185\conditional_null.csv
```

## Boundary

- The analysis uses donor pseudobulk and within-donor cell rank only.
- Neighborhood exposure is a kNN approximation, not an exact external Milo
  implementation.
- There is no external clinical outcome endpoint.
- The external results do not validate the aTAL2/frTAL subtype labels.
- The result should be presented as generalizability of the calibration
  pattern, not as external biological replication of the target panel.
