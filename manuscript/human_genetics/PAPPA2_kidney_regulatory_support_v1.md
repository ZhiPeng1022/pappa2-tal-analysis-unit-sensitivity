# PAPPA2 kidney regulatory support v1

Date: 2026-10-28

## Purpose

Test whether PAPPA2 has a human kidney regulatory context that can add
independent biological support to the analysis-unit sensitivity manuscript.

This module was prespecified to distinguish:

1. tubulointerstitial versus glomerular eQTL support;
2. disease-associated chromatin accessibility in an orthogonal modality;
3. whether the fixed aTAL2 marker panel is broadly enriched for kidney
   tubulointerstitial eGenes.

## Data sources

### NephQTL2 kidney eQTL

Han et al. Nat Commun 2023.

```text
DOI: 10.1038/s41467-023-37691-7
PMID: 37076491
```

The resource contains gene-level eQTL results from 240 glomerular (GLOM) and
311 tubulointerstitial (TUBE) human kidney samples. Supplementary Data 4 was
used for gene-level eQTL results; Supplementary Data 9 was used only to check
whether PAPPA2 appeared among FDR-significant PTWAS associations for eGFR or
UACR.

### GSE211785 snATAC

Abedini et al. Nat Genet 2024.

```text
DOI: 10.1038/s41588-024-01802-x
GEO: GSE211785
```

The non-KPMP multiome object was used as an orthogonal chromatin-accessibility
resource. The available H5AD contains a gene-activity representation for
55,955 snATAC nuclei. C_TAL and M_TAL cells were analyzed at the donor level.
This is same-cohort orthogonal evidence, not an independent cohort.

## Results

### PAPPA2 is a strong kidney tubulointerstitial eGene

| Compartment | Gene-level FDR | eGene | Model size |
|---|---:|---:|---:|
| TUBE | 1.64e-14 | yes | 1.035 |
| GLOM | 0.526 | no | 0.138 |

PAPPA2 therefore has strong kidney tubulointerstitial genetic regulation and no
glomerular eQTL support in NephQTL2.

### The aTAL2 marker panel is not broadly enriched for tubulointerstitial eGenes

The fixed aTAL2 top-35 marker set had 11 TUBE eGenes among 35 genes. Against all
18,027 gene symbols represented in the NephQTL2 gene-level table, this was not
an enrichment:

```text
one-sided Fisher exact P = 0.95
odds ratio = 0.58
```

This is exploratory because eGene detection depends on expression and other
gene properties. It does not support the claim that the aTAL2 marker program as
a whole is enriched for kidney eQTLs.

Several individual markers do have strong TUBE eQTLs, including TBC1D1, SPP1,
C4orf19, KCNIP4, and SOD2. These genes can be reported descriptively, but must
not be presented as independent validation selected after looking at their
eQTL results.

### PAPPA2 is not among FDR-significant NephQTL2 PTWAS genes

PAPPA2 was not present in NephQTL2 Supplementary Data 9, which lists only
TUBE/GLOM eGFR and UACR PTWAS associations with FDR < 0.05. This means that
PAPPA2 did not pass that resource's PTWAS significance filter. It is not a
negative test with a reported effect estimate, because non-significance may
also reflect instrument availability and power.

### Disease-associated snATAC gene activity does not support PAPPA2 elevation

GSE211785 snATAC analysis used four disease donors and seven control donors for
both TAL states. Donor-level exact permutations tested disease versus control.

| State | Difference disease minus control | Exact P | FDR |
|---|---:|---:|---:|
| C_TAL | -0.0293 | 0.843 | 1.0 |
| M_TAL | -0.0163 | 0.937 | 1.0 |

Neither direction was positive, and leave-one-donor ranges crossed zero.

The strongest exploratory snATAC signal was HIF1A in M_TAL
(difference = +0.120, exact P = 0.057, FDR = 0.995). It did not survive
multiple-testing correction and should not be promoted as a finding.

## Interpretation

The module supports a **kidney tubulointerstitial regulatory context** for
PAPPA2. It does not support:

- kidney-function causality;
- disease-associated PAPPA2 chromatin accessibility;
- a broad aTAL2 marker enrichment in kidney eQTL data;
- genetic colocalization with eGFR or UACR.

The correct manuscript claim is:

> PAPPA2 is a strong kidney tubulointerstitial eGene in an independent human
> kidney eQTL resource. This regulatory context is consistent with a
> tubule-compartment signal, but it does not establish disease causality, and
> same-cohort snATAC gene activity did not show higher disease accessibility.

## Decision

Keep this module as a concise boundary-of-evidence section and, if useful, one
small supplementary panel. Do not use it to claim genetic validation or
mechanistic support.

The main Research Article evidence remains:

1. systematic analysis-unit discordance in KPMP TAL cells;
2. donor-stratified signature transfer;
3. two non-KPMP external cohorts;
4. an explicit boundary analysis showing where regulatory and chromatin
   support stops.

This result changes the journal positioning slightly. The paper is stronger as
a rigorous methods-and-biology study than as a PAPPA2 genetics paper.

## Reproduction

Script:

```text
C:\Users\Elsa\Documents\Codex\2026-10-15\ji\work\analyze_kidney_regulatory_support.py
```

Outputs:

```text
C:\Users\Elsa\Documents\Codex\2026-10-15\ji\outputs\regulatory_support\
  NephQTL2_PAPPA2_gene_level.csv
  NephQTL2_fixed_marker_overlap.csv
  NephQTL2_fixed_marker_group_summary.csv
  GSE211785_snATAC_donor_pseudobulk.csv
  GSE211785_snATAC_gene_summary.csv
  regulatory_support_summary.csv
```
