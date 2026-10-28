# PAPPA2 prespecified analysis-unit panel results v1

Date: 2026-10-28

## Status

Completed once using the frozen 400-gene panel defined in:

```text
C:\Users\Elsa\Documents\Codex\2026-10-15\ji\outputs\PAPPA2_prespecified_unit_panel_v1.md
```

Frozen panel SHA256:

```text
d4942c72ed3596ddfa00f312c403e3e56d412d2aed9ad6e6501e03967ad0a4bf
```

No genes, thresholds, subtypes, comparisons, or outcomes were changed after
the primary result was observed.

## Design

Target panel:

- 100 TAL subtype markers;
- 100 Hallmark injury, inflammatory, hypoxic, and metabolic genes;
- 100 PAPPA2 donor-pseudobulk coexpression genes;
- 100 internal expression-matched negative controls.

Negative-control distribution:

- 2,000 expression- and detection-matched genes;
- 1,000 random sets of 400 genes;
- an additional expression-stratified sensitivity analysis.

Gene-specific units:

- donor pseudobulk;
- neighborhood-expression association;
- within-donor cell-rank association.

Milo neighborhood abundance remains a state-level unit and is not treated as a
gene-level effect.

## Primary result

The donor-versus-within reversal fraction was higher in the 400-gene target
panel than in expression-matched random gene sets.

| Analysis | Target | Random | Empirical P |
|---|---:|---:|---:|
| Target 400 | 18.57% | 12.68% | 0.000999 |
| Biological 300 | 20.66% | 12.70% | 0.000999 |
| Expression-stratified sensitivity | 18.57% | 12.77% | 0.000999 |

The minimum attainable empirical P with 1,000 random sets is 0.000999.

## Secondary results

### Donor versus neighborhood

| Analysis | Target | Random | Empirical P |
|---|---:|---:|---:|
| Target 400 | 22.06% | 14.31% | 0.000999 |
| Biological 300 | 24.38% | 14.34% | 0.000999 |
| Expression-stratified sensitivity | 22.06% | 14.40% | 0.000999 |

### Neighborhood versus within-donor

| Analysis | Target | Random | Empirical P |
|---|---:|---:|---:|
| Target 400 | 20.11% | 19.45% | 0.0949 |
| Biological 300 | 19.77% | 19.45% | 0.2987 |
| Expression-stratified sensitivity | 20.11% | 19.37% | 0.0709 |

This pair did not pass the corrected secondary threshold.

## Panel group results

Donor-versus-within reversal fractions:

| Panel group | Reversal fraction |
|---|---:|
| PAPPA2 coexpression | 25.08% |
| TAL subtype markers | 21.14% |
| Injury/inflammation/metabolism | 15.75% |
| Internal negative controls | 12.31% |
| External matched controls | 12.69% |

The internal negative-control group was statistically indistinguishable from
random 100-gene sets from the control pool:

```text
internal controls: 12.31%
random 100-gene sets: 12.60%
two-sided empirical P = 0.825
```

This argues against a general expression-level artefact.

## Disease-specific pattern

Donor-versus-within reversal remained higher than matched controls in AKI and
CKD:

| Comparison | Biological 300 | External controls |
|---|---:|---:|
| AKI | 18.83% | 13.63% |
| CKD | 19.50% | 11.88% |
| Combined disease | 18.00% | 12.23% |

## Interpretation

The prespecified analysis supports three claims:

1. Analysis-unit reversal is more frequent in biologically selected TAL genes
   than in expression-matched controls.
2. The result is present for both donor-versus-within and
   donor-versus-neighborhood comparisons.
3. The neighborhood-versus-within contrast does not show a panel-level
   excess.

It does not establish that every selected gene is robust or causal. The result
is a panel-level enrichment test, and individual top genes require separate
validation.

## Output files

```text
C:\Users\Elsa\Documents\Codex\2026-10-15\ji\outputs\prespecified_unit_panel\
  target_panel_400.csv
  negative_control_pool_2000.csv
  panel_manifest.json
  PAPPA2_training_donor_coexpression.csv
  Hallmark_candidate_ranking.csv
  TAL_gene_expression_summary.csv

C:\Users\Elsa\Documents\Codex\2026-10-15\ji\outputs\prespecified_unit_panel_results\
  gene_unit_metrics_long.csv
  gene_unit_summary.csv
  random_set_results.csv
  primary_result.json
```
