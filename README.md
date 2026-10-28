# PAPPA2 TAL analysis-unit sensitivity and Human Genetics Research Article archive

Public code and derived-data repository:

```text
https://github.com/ZhiPeng1022/pappa2-tal-analysis-unit-sensitivity
```

Zenodo concept DOI:

```text
https://doi.org/10.5281/zenodo.23196128
```

Zenodo exact DOI for the analyses and Human Genetics package:

```text
https://doi.org/10.5281/zenodo.23733571
```

Analysis code for the manuscript:

```text
Analysis-unit sensitivity distinguishes donor-level expression enrichment
from within-donor coupling in human kidney thick ascending limb states
```

## Overview

This repository contains the complete-data analysis code, curated derived
tables, and generated figures for four human thick ascending limb (TAL)
subtypes: C-TAL-A, C-TAL-B, aTAL2, and frTAL.

The analysis compares:

1. Milo neighborhood abundance testing;
2. Augur disease-versus-reference separability;
3. donor-level PAPPA2 pseudobulk testing;
4. complete-data neighborhood-level expression association;
5. within-donor cell-level rank analysis;
6. KPMP eGFR association models;
7. KPMP GSE183456 spatial cross-modality analysis;
8. external public eGFR and Nephroseq association checks;
9. independent non-KPMP validation in GSE211785 and GSE254185;
10. donor-stratified TAL signature transfer;
11. NephQTL2 kidney tubulointerstitial eQTL context;
12. GSE211785 snATAC gene-activity boundary analysis;
13. donor-level composition and expression decomposition.

PAPPA2 is analyzed as a case gene. The study does not claim a new
algorithm, software method, or causal mechanism.

## Data sources

| Data | Accession or source |
|---|---|
| KPMP snRNA-seq | GEO superseries `GSE183279`; snRNA component `GSE183277` |
| KPMP Visium | GEO `GSE183456` |
| Independent multi-omic validation | GEO `GSE211785`; non-KPMP object |
| Independent injury validation | GEO `GSE254185` |
| External public eGFR | GEO `GSE137570` |
| External GFR associations | Nephroseq v5 |
| Candidate external resource, not used | GEO `GSE115098` |

`GSE183456` is part of the KPMP spatial resource and is not an external
cohort. `GSE115098` was not used for external validation because
sample-level GFR annotations are not publicly available.

## Independent validation

Release `v1.0.4` adds two non-KPMP human kidney validation resources.

`GSE211785` contains C_TAL and M_TAL RNA annotations. After requiring at least
five TAL cells per donor, 29,847 TAL cells from 35 donors were available
(16 control, 7 DKD, and 12 HKD). PAPPA2 was higher in disease than control
donors in C_TAL (+0.283; permutation p = 0.0315; six-test FDR = 0.063) and was
directionally positive in M_TAL (+0.200; p = 0.105). Leave-one-donor-out
estimates remained positive.

`GSE254185` contains 9,912 annotated TAL cells from 7 normal and 5
obstructive-nephropathy kidneys. PAPPA2 donor pseudobulk was positive in mTAL
(+0.341; exact p = 0.0121; five-state FDR = 0.060) and in TAL Injured. In
TAL Injured, neighborhood-expression and within-donor directions were negative,
showing a qualitative cross-unit signal. The exact aTAL2 and frTAL labels were
not available in either independent cohort.

These analyses are independent-cohort validation, not experimental or causal
validation. The detailed outputs are in `data/derived/independent_validation/`
and the scripts are in `analysis/independent_validation/`.

## Release v1.0.5: Human Genetics Research Article candidate

Release `v1.0.5` adds the Research Article working manuscript, Human Genetics
cover letter, 20 supplementary tables, eight submission TIFF figures, and the
new donor-stratified signature-transfer and kidney-regulatory analyses.

The signature-transfer analysis trains on 156 KPMP donors, tests on 67
held-out donors, and recovers an aTAL2-like state in GSE211785 and GSE254185.
The external data support a donor-level direction but do not provide exact
aTAL2 label replication or a significant within-donor effect.

NephQTL2 identifies PAPPA2 as a strong tubulointerstitial eGene
(FDR = 1.64e-14) without glomerular eQTL support. PAPPA2 is not among the
FDR-significant PTWAS genes, and GSE211785 snATAC does not show higher disease
accessibility. These results are retained as a regulatory-context and boundary
analysis.

The donor-level decomposition shows that the overall four-subtype TAL PAPPA2
difference is dominated by subtype expression rather than subtype composition.
In the combined disease comparison, the observed difference is +0.563, the
expression component is +0.641, the composition component is +0.010, and the
interaction/covariance term is -0.088.

## Directory structure

```text
repo/
├─ data/
│  └─ derived/
│     ├─ supplementary/       # BMC additional files 1-14 (submission naming)
│     ├─ independent_validation/  # GSE211785 and GSE254185 result tables
│     ├─ signature_transfer/  # Donor-stratified marker transfer summaries
│     ├─ regulatory_support/  # NephQTL2 eQTL and GSE211785 snATAC results
│     ├─ composition_expression/  # Donor-level decomposition
│     ├─ human_genetics_supplementary/  # Supplementary Tables S1-S20
│     └─ reanalysis/          # Full-data reanalysis and manifest summaries
├─ analysis/
│  ├─ run_milo_tal_kpmp.R
│  ├─ prepare_full_tal_milo_metadata.py
│  ├─ analyze_full_tal_donor_pseudobulk.py
│  ├─ analyze_full_tal_pseudobulk_robustness.py
│  ├─ analyze_full_tal_nhood_expression.py
│  ├─ sensitivity_full_tal_nhood_expression.py
│  ├─ analyze_full_tal_within_donor_cell_rank.py
│  ├─ summarize_full_tal_analysis_units.py
│  ├─ run_pyaugur_tal_kpmp.py
│  ├─ summarize_tal_milo_augur.py
│  ├─ analyze_tal_nhood_focus_expression.py
│  ├─ sensitivity_tal_nhood_expression.py
│  ├─ analyze_tal_within_donor_cell_rank.py
│  ├─ analyze_pappa2_egfr_adjusted.py
│  ├─ analyze_egfr_gene_specificity.py
│  ├─ analyze_gse137570_pappa2_egfr.py
│  ├─ query_nephroseq_gfr_genes.py
│  ├─ independent_validation/ # Non-KPMP validation scripts
│  └─ human_genetics/         # Signature transfer, eQTL, decomposition scripts
├─ figures/
│  ├─ make_figures_v2.py
│  ├─ make_supplementary_figures.py
│  ├─ check_figures.py
│  └─ generated/              # Final figure PNG files
├─ manuscript/
│  ├─ human_genetics/         # Research Article manuscript and cover letter
│  ├─ build_manuscript_evidence_matrix.py
│  ├─ search_refs.py
│  ├─ verify_refs.py
│  ├─ markdown_to_docx.py
│  └─ markdown_to_pdf.py
├─ requirements.txt
├─ CITATION.cff
├─ .zenodo.json
└─ README.md
```

Archive metadata are provided in `CITATION.cff` and `.zenodo.json`.
The GitHub repository is the working public archive. Zenodo creates a
versioned DOI archive from tagged GitHub releases; the concept DOI above
resolves to the latest archived version.

## Submission files (BMC Medical Genomics)

The supplementary material follows the BMC requirement that additional files
are named `Additional file N` and cited in the manuscript in numeric order.
The mapping is listed in
`data/derived/supplementary/Additional_files_index.md`.

Two files that were part of the v1.0.2 archive are no longer used:

| Removed file | Reason |
|---|---|
| `Table_S9_full_TAL_pseudobulk.csv` | Identical to `Table_S3_full_pseudobulk.csv` (`Additional file 9.csv`) |
| `Figure_S1_analysis_units.png` | Duplicated the analysis-unit design already shown in main Fig. 3 |

Main Fig. 3 is retained. The statistics shown in the removed supplementary
schematic (500 donor-stratified resamples, purity-adjusted partial
correlation, 1,000 within-donor permutations, 3/5/10 cell thresholds) are
already reported in the Methods section of the manuscript.

## Submission files (Human Genetics)

The Human Genetics Research Article working files are in
`manuscript/human_genetics/`. The main DOCX is the editable manuscript, the
PDF is the review copy, and the TXT source is the Markdown file. The submission
package contains eight TIFF figures and Supplementary Tables S1-S20.

The final manuscript cites `https://doi.org/10.5281/zenodo.23733571`.

## Release v1.0.6: final DOI binding

Release `v1.0.6` updates the manuscript, repository metadata, and citation
files to cite the exact `v1.0.5` archive:

```text
https://doi.org/10.5281/zenodo.23733571
```

No analysis results changed in this release.

## Software

Python packages and versions are listed in `requirements.txt`. The Milo
analysis used R 4.6.1 and `miloR` 2.8.1, with
`SingleCellExperiment` 1.34.0, `edgeR` 4.10.5, and `limma` 3.68.5.

## Random seeds and key parameters

| Parameter | Value |
|---|---|
| Milo seed | `20260927` |
| Milo normalization | `logMS` |
| Complete TAL cells | `106,851` |
| Old comparator subset | `4,446` |
| Donor-stratified resamples | 500 |
| Within-donor permutations | 1,000 |
| Augur subsamples | 20 |
| Augur cells per subsample | 30 |
| Augur cross-validation folds | 3 |
| Cell thresholds | 3, 5, 10 |
| Primary cell threshold | 5 |
| Milo threshold | `SpatialFDR < 0.05` |
| Signature high-confidence margin | 0.2 |
| Signature bootstrap / permutation cap | 20,000 |
| Decomposition bootstrap resamples | 2,000 |

## Path configuration

The complete-data analysis chain uses:

```text
PAPPA2_PROJECT_ROOT
PAPPA2_DATA_ROOT
PAPPA2_OUTPUT_ROOT
PAPPA2_WORK_ROOT
PAPPA2_VALIDATION_DIR
PAPPA2_VALIDATION_OUTPUT
PAPPA2_VALIDATION_DATA
GSE211785_H5AD
```

`analysis/project_paths.py` resolves defaults relative to the repository.
Set `PAPPA2_DATA_ROOT` to the directory containing `KPMP_new_snRNA.h5ad`:

```powershell
$env:PAPPA2_DATA_ROOT = "D:\path\to\RAMP3"
```

The Milo script also accepts:

```text
MILO_METADATA_PATH
MILO_OUTPUT_DIR
MILO_OUTPUT_PREFIX
MILO_SEED
MILO_EXPORT_MEMBERSHIP
MILO_NORM_METHOD
```

All scripts resolve project, output, work, and data paths from the
`PAPPA2_*` environment variables. The optional variables are
`PAPPA2_PROJECT_ROOT`, `PAPPA2_DATA_ROOT`, `PAPPA2_OUTPUT_ROOT`, and
`PAPPA2_WORK_ROOT`.

## Reproducibility notes

1. Neighborhoods overlap, so conventional correlation p values are
   descriptive rather than independent-sample inference.
2. The neighborhood-level and within-donor cell-level analyses use
   different units and should not be combined.
3. Nephroseq analysis rows are not independent cohorts.
4. No causal or direct cell-cell communication claims are supported by
   these analyses.

## Citation

The manuscript citation will be added after publication. For the code and
derived-data archive, cite the Zenodo concept DOI
`https://doi.org/10.5281/zenodo.23196128`.

## Contact

```text
Peng Zhi, zhiyoyo0409@163.com
```
