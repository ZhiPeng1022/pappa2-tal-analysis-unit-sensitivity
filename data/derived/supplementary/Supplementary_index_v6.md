# Supplementary Material Index v6

Manuscript: `Subtype-specific thick ascending limb states and analysis-unit
sensitivity in human kidney disease`

Date: 2026-10-07

Repository: `https://github.com/ZhiPeng1022/pappa2-tal-analysis-unit-sensitivity`

Zenodo concept DOI: `https://doi.org/10.5281/zenodo.23196128`

## Supplementary tables

| Table | File | Content |
|---|---|---|
| Table S1 | `Table_S1_milo_full_summary.csv` | Complete-data Milo neighborhood counts, directions and median logFC |
| Table S2 | `Table_S2_augur_subset_summary.csv` | Secondary subset-level Augur separability summary |
| Table S3 | `Table_S3_full_pseudobulk.csv` | Complete donor pseudobulk results for five genes |
| Table S4 | `Table_S4_full_neighborhood_robustness.csv` | Complete-data neighborhood association with leave-one-donor and bootstrap stability |
| Table S5 | `Table_S5_full_within_donor_cell_rank.csv` | Complete-data within-donor cell-rank results |
| Table S6 | `Table_S6_kpmp_egfr_models.csv` | KPMP ordinal eGFR models |
| Table S7 | `Table_S7_gse137570_validation.csv` | GSE137570 directional eGFR association |
| Table S8 | `Table_S8_nephroseq_significant_analyses.csv` | Nephroseq precomputed GFR associations |
| Table S9 | `Table_S9_full_TAL_pseudobulk.csv` | Duplicate machine-readable source for Table S3 if required |
| Table S10 | `Table_S10_full_vs_subset_direction.csv` | Complete versus 4,446-cell subset direction comparison |
| Table S11 | `Table_S11_full_analysis_unit_concordance.csv` | Joint donor pseudobulk, neighborhood, and within-donor comparison |
| Table S12 | `Table_S12_full_nhood_bootstrap_rhos.csv.gz` | Donor-stratified bootstrap correlation replicates |
| Table S13 | `Table_S13_full_nhood_lodo_rhos.csv.gz` | Leave-one-donor neighborhood-association replicates |
| Table S14 | `Table_S14_full_within_donor_null.csv.gz` | Permutation null distributions for the complete-data cell-rank analysis |

## Supplementary figures and diagnostics

| Item | File | Content |
|---|---|---|
| Figure S1 | `Figure_S1_analysis_units.png` | Analysis-unit schematic |
| Diagnostic S1 | `KPMP_full_TAL_within_donor_cell_rank_diagnostics.csv` | Eligible donors and model-fit diagnostics |
| Summary | `PAPPA2_full_data_reanalysis_summary.md` | Full-data reanalysis summary for internal and response-letter use |
| Data and code manifest | `PAPPA2_code_data_manifest_v3.md` | Full-data analyses, software versions, seeds and archive checklist |
| Reference verification | `PAPPA2_reference_verification_v1.md` | DOI and PMID verification record |

## Important interpretation notes

1. Complete-data Milo used the atlas-provided UMAP and logMS normalization.
2. One AKI donor without neighborhood representation was excluded from the
   AKI differential-abundance test.
3. PAPPA2 was directionally positive in all 12 donor-pseudobulk
   subtype-comparison groups; 11 of 12 passed FDR.
4. aTAL2 cross-unit discordance was significant in combined disease and AKI
   but not CKD.
5. frTAL had a smaller related cross-unit pattern, significant only in the
   combined disease comparison.
6. Nephroseq rows are not independent cohorts.
7. GSE183456 is a KPMP spatial resource, not an external cohort.
