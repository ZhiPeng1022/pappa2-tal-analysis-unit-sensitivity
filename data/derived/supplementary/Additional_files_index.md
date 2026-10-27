# Additional files index

Manuscript: `Subtype-specific thick ascending limb states and analysis-unit
sensitivity in human kidney disease`

Target journal: `BMC Medical Genomics`

Date: 2026-10-28

Repository: `https://github.com/ZhiPeng1022/pappa2-tal-analysis-unit-sensitivity`

Zenodo concept DOI: `https://doi.org/10.5281/zenodo.23196128`

Additional files are numbered in the order in which they are first cited in
the main manuscript, as required by BMC.

## Additional files

| Upload name | Format | Title | Content |
|---|---|---|---|
| `Additional file 1.csv.gz` | .csv.gz | Permutation null distributions for the within-donor cell-rank analysis | Null distributions from 1,000 within-donor expression permutations |
| `Additional file 2.csv` | .csv | Complete-data Milo neighborhood summary | Neighborhood counts, directions and median logFC |
| `Additional file 3.csv` | .csv | Subset-level Augur separability summary | Mean AUC and subsample intervals in the 4,446-cell subset |
| `Additional file 4.csv` | .csv | Complete versus subset pseudobulk direction comparison | Direction of effect in the complete and subset analyses |
| `Additional file 5.csv` | .csv | Neighborhood-expression association robustness | Spearman correlations with leave-one-donor and bootstrap stability |
| `Additional file 6.csv.gz` | .csv.gz | Donor-stratified bootstrap replicates for neighborhood association | Per-replicate correlation values |
| `Additional file 7.csv` | .csv | Joint analysis-unit concordance table | Donor pseudobulk, neighborhood and within-donor directions joined |
| `Additional file 8.csv` | .csv | Within-donor cell-rank association results | Donor-level correlations between expression and neighborhood exposure |
| `Additional file 9.csv` | .csv | Complete donor pseudobulk results | Effect sizes, adjusted estimates and false discovery rates |
| `Additional file 10.csv.gz` | .csv.gz | Leave-one-donor-out replicates for neighborhood association | Per-replicate correlations with one donor removed |
| `Additional file 11.csv` | .csv | KPMP ordinal eGFR models | Ordinal logistic estimates for PAPPA2-positive TAL fraction |
| `Additional file 12.csv` | .csv | GSE137570 directional eGFR association | Unadjusted and adjusted association estimates |
| `Additional file 13.csv` | .csv | Nephroseq precomputed GFR associations | Precomputed GFR results for PAPPA2 and four control genes |
| `Additional file 14.csv` | .csv | Independent non-KPMP validation of PAPPA2 and analysis-unit direction | GSE211785 and GSE254185 effects, permutation p values, FDRs and donor counts |

## Diagnostic and reference files not submitted as additional files

| File | Content |
|---|---|
| `KPMP_full_TAL_within_donor_cell_rank_diagnostics.csv` | Eligible donors and model-fit diagnostics |
| `PAPPA2_full_data_reanalysis_summary.md` | Full-data reanalysis summary |
| `PAPPA2_code_data_manifest_v3.md` | Analyses, software versions, seeds and archive checklist |
| `PAPPA2_reference_verification_v1.md` | DOI and PMID verification record |

## Changes relative to v1.0.2

Two files were removed when the manuscript was adapted to BMC:

| Removed file | Reason |
|---|---|
| `Table_S9_full_TAL_pseudobulk.csv` | Identical to `Table_S3_full_pseudobulk.csv` (now `Additional file 9.csv`) |
| `Figure_S1_analysis_units.png` | Duplicated the analysis-unit design already shown in main Fig. 3 |

Main Fig. 3 is retained. All statistical parameters shown in the removed
supplementary schematic (500 donor-stratified resamples, purity-adjusted
partial correlation, 1,000 within-donor permutations, 3/5/10 cell
thresholds) are already stated in the Methods section of the manuscript.

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
8. GSE211785 and GSE254185 are non-KPMP independent cohorts; they do not
   contain the exact aTAL2/frTAL discovery labels.
