# PAPPA2 prespecified analysis-unit panel v1

Date: 2026-10-28

Status: frozen before outcome analysis

## Purpose

Test whether analysis-unit direction reversal is concentrated in biologically
defined TAL disease genes rather than in expression-matched control genes.

The four-unit hierarchy is retained, but only three units are gene-specific:

1. donor pseudobulk expression;
2. neighborhood-expression association;
3. within-donor cell-rank association.

Milo neighborhood abundance is state-level and does not produce a gene-specific
direction. It is therefore not included as a fourth gene-level test.

## Target panel

The target panel has exactly 400 unique genes, 100 per group:

1. TAL subtype markers: top 25 genes for each of C-TAL-A, C-TAL-B, aTAL2, and
   frTAL from the existing donor-stratified KPMP marker table.
2. Injury, inflammatory, hypoxic, and metabolic genes: union of five
   MSigDB Hallmark 2020 sets, filtered to genes detected in KPMP TAL cells,
   then ranked by reference-donor mean expression.
3. PAPPA2 coexpression genes: top genes by absolute Spearman correlation with
   PAPPA2 in donor pseudobulk from the 156 KPMP training donors, excluding
   genes already assigned to groups 1 and 2.
4. Expression-matched negative controls: 100 genes selected from the external
   negative-control pool without using disease or analysis-unit results.

Group assignment is disjoint and priority-ordered. A gene already assigned to
an earlier group cannot enter a later group.

## Negative-control pool

A 2,000-gene pool is selected from genes not in the 400-gene target panel and
matched to the target panel by mean TAL expression and detection decile.
One thousand random sets of 400 genes are sampled from this pool without
replacement within each set.

## Expression filters

A gene is eligible in a subtype and comparison when:

- at least 10% of TAL cells in that subtype have non-zero expression;
- at least 10 donors have non-zero expression;
- the gene has non-zero variance in both disease and Reference donor groups.

Genes failing these filters remain in the frozen list and are reported as
non-evaluable in that group.

## Unit-specific effects

- Donor pseudobulk effect: disease donor mean minus Reference donor mean.
- Neighborhood-expression effect: Spearman correlation across Milo
  neighborhoods between neighborhood mean gene expression and Milo logFC.
- Within-donor effect: median disease donor Spearman rho minus median
  Reference donor Spearman rho for expression versus neighborhood exposure.

Effect sizes are standardized within each subtype and disease comparison using
the pooled target plus negative-control gene distribution.

## Primary endpoint

A donor-versus-within reversal is counted when:

- donor pseudobulk and within-donor effects have opposite signs;
- both absolute standardized effects are at least 0.5.

For each gene, the reversal fraction is computed across evaluable
subtype-comparison groups. The panel statistic is the mean reversal fraction
across the 400 target genes. The empirical P value compares this statistic with
1,000 expression-matched random 400-gene sets.

One-sided empirical P:

```text
(1 + number of random sets >= observed) / 1001
```

The primary test is donor-versus-within reversal.

## Secondary endpoints

- donor-versus-neighborhood reversal;
- neighborhood-versus-within reversal;
- reversal fraction in the 300 biological genes versus the 100 internal
  negative-control genes.

Secondary results are reported with Bonferroni correction across three tests.

## Decision rule

- A significant primary result supports a general analysis-unit reversal
  program in biologically selected TAL genes.
- A null result is retained as a null result.
- No gene set, threshold, subtype, or comparison will be changed after the
  primary result is observed.
