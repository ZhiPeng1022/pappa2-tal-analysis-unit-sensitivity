# PAPPA2 Human Genetics Research Article route v1

Date: 2026-10-28

## Decision

The manuscript will be developed as a full **Research Article**, not a Brief
Report.

Primary target:

```text
Human Genetics (Springer)
```

This is the hybrid Springer journal with ISSN 0340-6717 / 1432-1203, not the
fully open-access BMC/Springer Nature journal Human Genomics. The traditional
subscription route does not charge an APC. Open Choice remains optional.

Official journal information checked on 2026-10-28:

- Research Articles are defined as definitive and comprehensive descriptions
  of major research findings of broad significance.
- Brief Reports have a separate 750-1,000 word limit, no formal section
  headings, and no more than two figures/tables.
- Research Article abstracts should be 150-250 words.
- Provide 4-6 keywords.
- References use author-year style.
- Use no more than three displayed heading levels.

Therefore, the 750-1,000 word limit and two-display-item limit do not apply to
the Research Article version.

## Current evidence inventory

The study already contains full-research-scale discovery data and two
non-KPMP support cohorts:

- KPMP discovery data: 106,851 TAL cells;
- four prespecified analysis units;
- complete-versus-subset comparison: direction changed in 12/60 comparisons;
- aTAL2 cross-unit discordance with FDR-supported internal results;
- donor-stratified KPMP signature validation: 93.90% accuracy among assigned
  held-out cells;
- GSE254185: non-KPMP TAL multiome cohort;
- GSE211785: non-KPMP kidney RNA/ATAC cohort;
- donor-level PAPPA2 disease elevation supported in both external cohorts;
- external within-donor directions are negative but underpowered;
- public eGFR, Nephroseq, and KPMP spatial results provide supporting but not
  definitive evidence;
- complete code and derived-data archive already exist.

What is still missing for a competitive Research Article is not word count. It
is a second form of independent evidence that links the TAL disease state to
human kidney regulatory or genetic variation and makes the result consequential
beyond one case gene.

## Central argument

The article should be framed as a general measurement problem in single-cell
kidney genomics:

> The unit of inference is an effect modifier. Donor composition, spatial
> neighborhood abundance, neighborhood-expression coupling, and within-donor
> cell ranking answer different questions. In human TAL data, changing the
> analysis unit can reverse a subtype-level disease conclusion. This behavior
> is reproducible in non-KPMP cohorts and must therefore be reported explicitly.

PAPPA2 remains the prespecified case gene. It is not the central mechanism and
the paper must not claim causality, a validated biomarker, or exact aTAL2
replication in external cohorts.

## Research Article structure

### Abstract

Use a conventional Background-Methods-Results-Conclusions abstract of
180-230 words. State the analysis-unit problem in the first two sentences,
report the two external cohorts, and end with a reporting recommendation rather
than a causal claim.

### Introduction

End with three prespecified aims:

1. determine whether four established analysis units produce concordant
   subtype-level disease conclusions in KPMP TAL cells;
2. test whether the direction pattern transfers to non-KPMP kidney cohorts;
3. determine whether the disease-state signal has human kidney
   regulatory/genetic support.

### Results modules

1. **Analysis-unit discordance is common in TAL states.**
   Present the full KPMP map, subtype coverage, neighborhood heterogeneity,
   complete-versus-subset direction changes, and all four analysis units.

2. **aTAL2 is the strongest boundary case, not the only finding.**
   Report aTAL2 and frTAL without turning PAPPA2 into a gene-discovery paper.
   Keep the prespecified control genes and all tested states visible.

3. **Donor-stratified transfer recovers an aTAL2-like state outside KPMP.**
   Report internal held-out validation, threshold sensitivity, external
   projection, and backbone-state composition.

4. **External cohorts preserve donor-level association but not the complete
   four-unit pattern.**
   GSE211785 supports donor-level disease elevation. GSE254185 preserves a
   qualitative analysis-unit conflict. State clearly which claims fail.

5. **Human kidney regulatory and genetic evidence triangulates the signal.**
   This module must be added before calling the paper submission-ready as a
   competitive Research Article.

6. **A reporting framework follows from the empirical results.**
   Recommend that single-cell disease studies name the analysis unit, show at
   least one donor-aware and one within-donor result, and avoid interpreting
   all three as the same disease association.

## Required additional module: kidney regulatory/genetic triangulation

The analysis must be prespecified and reported whether positive or null.

### Primary question

Are the PAPPA2-associated TAL disease signal and the aTAL2-like signature
enriched for human kidney regulatory activity or kidney-function genetic risk?

### Inputs

- GTEx kidney cortex eQTL;
- a disease-context human kidney eQTL resource if a suitable public,
  non-KPMP dataset can be verified;
- GSE211785 snATAC-seq to test TAL accessibility at PAPPA2 and prespecified
  signature genes;
- CKDGen or another current, citable GWAS summary statistic for eGFR, CKD, or
  kidney function;
- the prespecified PAPPA2 case gene and the fixed aTAL2 marker panel.

### Analyses

1. Test whether PAPPA2 has a kidney-cortex eQTL and report all credible
   variant-gene associations without selecting only significant variants.
2. Test whether PAPPA2 and the fixed marker panel overlap TAL-accessible
   chromatin in GSE211785.
3. Use a prespecified gene-set test for kidney-function GWAS enrichment.
4. Report effect size, permutation or empirical P value, multiple-testing
   correction, and sensitivity to marker count and gene-window definition.
5. Include negative-control gene sets matched for expression and gene length.

### Interpretation rules

- A positive regulatory result can support a functional context, not
  causality.
- A null GWAS result must be reported and must not be reframed as evidence
  against a disease association.
- Do not add the analysis if the gene set, phenotype, or statistical test is
  chosen only after inspecting the result.
- Do not use a token genetic analysis as decoration.

## Main figure plan

Reduce the BMC eight-figure layout to six main figures. Robustness panels move
to supplementary information.

1. Study design, nephron/TAL location, and four analysis units.
2. TAL state map, donor composition, and neighborhood heterogeneity.
3. Four-unit direction comparison and complete-versus-subset instability.
4. PAPPA2 case analysis and kidney-function association.
5. Donor-stratified transfer and two non-KPMP external cohorts.
6. Kidney regulatory/genetic triangulation and proposed reporting framework.

## Submission sequence

1. Freeze the analysis-unit definitions and external-cohort claims.
2. Complete the kidney regulatory/genetic triangulation with a prespecified
   analysis plan.
3. Rebuild the manuscript in Human Genetics Research Article format.
4. Convert references from BMC numeric style to author-year style.
5. Reduce to six main figures and move lower-value panels to supplementary
   information.
6. Update the public archive to a new exact-version DOI after the new results
   are final.
7. Prepare a journal-specific cover letter that foregrounds the general
   analysis-unit problem and the independent support, not PAPPA2 alone.

## Go/no-go rule

- If the kidney regulatory/genetic module is positive and survives the
  prespecified correction, submit to Human Genetics as a Research Article.
- If the module is null but the two external cohorts remain directionally
  supportive, the paper is still a full Research Article by scope. In that
  case, use Human Genetics as an ambitious first attempt and have a
  lower-risk hybrid or subscription journal ready as the next target.
- Do not compress the paper into a Brief Report merely to fit the current
  evidence. If the evidence is insufficient for Research Article standards,
  change the journal rather than misrepresent the study's scope.

## Regulatory-support progress

Completed on 2026-10-28:

- PAPPA2 is a strong TUBE eGene in NephQTL2 (FDR = 1.64e-14);
- PAPPA2 has no glomerular eQTL support in that resource (FDR = 0.526);
- the fixed aTAL2 top-35 marker set is not broadly enriched for TUBE eGenes;
- PAPPA2 is not among NephQTL2 FDR-significant PTWAS genes for eGFR or UACR;
- GSE211785 snATAC gene activity does not show higher disease accessibility
  for PAPPA2 in C_TAL or M_TAL.

The result provides a kidney tubulointerstitial regulatory context but does not
provide genetic causality or disease-associated chromatin validation. The
manuscript should present this as a boundary analysis and retain the main
emphasis on analysis-unit sensitivity and external cohort transfer.

Detailed report:

```text
C:\Users\Elsa\Documents\Codex\2026-10-15\ji\outputs\PAPPA2_kidney_regulatory_support_v1.md
```
