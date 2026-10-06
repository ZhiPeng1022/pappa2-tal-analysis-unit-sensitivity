# PAPPA2 TAL analysis-unit sensitivity analysis

Public code and derived-data repository:

```text
https://github.com/ZhiPeng1022/pappa2-tal-analysis-unit-sensitivity
```

Analysis code for the manuscript:

```text
Subtype-specific thick ascending limb states and analysis-unit
sensitivity in human kidney disease
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
8. external public eGFR and Nephroseq association checks.

PAPPA2 is analyzed as a case gene. The study does not claim a new
algorithm, software method, or causal mechanism.

## Data sources

| Data | Accession or source |
|---|---|
| KPMP snRNA-seq | GEO superseries `GSE183279`; snRNA component `GSE183277` |
| KPMP Visium | GEO `GSE183456` |
| External public eGFR | GEO `GSE137570` |
| External GFR associations | Nephroseq v5 |
| Candidate external resource, not used | GEO `GSE115098` |

`GSE183456` is part of the KPMP spatial resource and is not an external
cohort. `GSE115098` was not used for external validation because
sample-level GFR annotations are not publicly available.

## Directory structure

```text
repo/
├─ data/
│  └─ derived/
│     ├─ supplementary/       # Supplementary Tables S1-S14 and Figure S1
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
│  └─ query_nephroseq_gfr_genes.py
├─ figures/
│  ├─ make_figures_v2.py
│  ├─ make_supplementary_figures.py
│  ├─ check_figures.py
│  └─ generated/              # Final figure PNG files
├─ manuscript/
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
versioned DOI archive from tagged GitHub releases.

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

## Path configuration

The complete-data analysis chain uses:

```text
PAPPA2_PROJECT_ROOT
PAPPA2_DATA_ROOT
PAPPA2_OUTPUT_ROOT
PAPPA2_WORK_ROOT
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

```text
[Manuscript citation to be added]
```

## Contact

```text
Peng Zhi, zhiyoyo0409@163.com
```
