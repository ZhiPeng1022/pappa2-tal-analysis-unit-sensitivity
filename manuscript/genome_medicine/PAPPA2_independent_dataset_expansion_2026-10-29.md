# PAPPA2 independent dataset expansion

Date: 2026-10-29

## Question

Can the analysis-unit result be tested in a non-KPMP human kidney
single-cell or spatial cohort without adding experimental work?

## Search performed

GEO was searched through NCBI E-utilities for human kidney single-cell,
single-nucleus, and spatial studies in disease settings. The broad query
returned 860 records; the first 200 were screened by title, summary, sample
count, and tissue type. Targeted searches were also run for thick ascending
limb, kidney biopsy, AKI, CKD, diabetic kidney disease, lupus nephritis, IgA
nephropathy, transplant rejection, and spatial transcriptomics.

## Highest-priority new candidate

### GSE284715 / GSE284742

Human kidney transplant biopsies:

- 3 stable allografts;
- 3 T-cell-mediated rejection biopsies;
- 22,183 annotated human nuclei in the public metadata table;
- 3,607 TAL-lineage nuclei across cortical TAL, medullary TAL, and TAL injury
  states;
- donor/sample identifiers and disease groups are present;
- the publication explicitly reports TAL injury states;
- `PAPPA2` is independently listed as a `TAL_Injury_h1` marker in the article's
  supplementary marker table.

This is the strongest new candidate for an exact or near-exact TAL analysis
outside KPMP. It tests transplant rejection rather than AKI/CKD, so it would
support analysis-unit portability rather than the original disease contrast.

Data access:

- processed human count matrices are available as per-sample Matrix Market
  files in GEO;
- six human matrix files total approximately 732 MB before barcodes and
  features;
- the combined RAW archive is approximately 2.16 GB.

Current blocker:

- GEO transfer was approximately 6-15 KB/s during this work session;
- a small 0.5 MB metadata file completed, but a 20 MB barcode file did not
  complete within a practical interval.

No quantitative result from GSE284715 is claimed in the current manuscript.

## Other candidates

| Accession | Content | Strength | Main limitation |
|---|---|---|---|
| GSE195460 | Human diabetic kidney disease snRNA/snATAC | Disease versus control, multiome | One control donor; partial overlap with GSE131882 |
| GSE131882 | Early diabetic nephropathy, 3 disease and 3 control libraries | Small, direct disease comparison | Broad tubule labels, older processing |
| GSE286911 | IgA nephropathy snRNA-seq | Biopsy disease progression | Only six newly generated libraries |
| GSE317012 | BK polyomavirus-associated nephropathy biopsies | Serial samples and tubular injury | 1.9 GB processed archive; infection-specific context |
| GSE232222 | CKD snRNA/snATAC plus spatial data | Large CKD resource | Not an untouched external cohort for the current paper without full reannotation |
| GSE271188 | Human glomerular disease snRNA-seq | Disease diversity and Xenium companion | Glomerular sampling may underrepresent TAL |
| GSE315877 | Paired pre/post metabolic surgery kidney biopsies | Five paired donors with TAL cells | No untreated disease-reference contrast |

## Recommendation

1. Prioritize GSE284715 because it has explicit TAL states and a usable
   disease-reference design.
2. Retrieve the six human matrices through a faster institutional mirror or an
   overnight transfer.
3. Recompute donor pseudobulk, TAL-injury marker score, and within-donor
   PAPPA2-marker coupling.
4. Report this as near-exact TAL portability, not as an exact aTAL2
   replication.
5. Keep the current manuscript's statement that exact external aTAL2
   replication remains unavailable until that analysis is complete.

## Boundary

Dataset availability is not validation. The current work did not obtain or
analyze the GSE284715 expression matrix, so no new quantitative claim has been
added to the manuscript or submission package.
