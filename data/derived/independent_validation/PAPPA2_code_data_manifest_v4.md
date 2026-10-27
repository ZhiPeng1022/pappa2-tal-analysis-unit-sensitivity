# PAPPA2 independent-validation manifest v4

Date: 2026-10-28

Release: 1.0.4

## Data sources

| Dataset | Role | Access |
|---|---|---|
| GSE211785 | Non-KPMP donor-level validation | GEO; non-KPMP PreSCVI object |
| GSE254185 | Non-KPMP analysis-unit validation | GEO; per-library GEX matrices and metadata |
| GSE183279 / GSE183456 | KPMP discovery and same-consortium spatial support | GEO |
| GSE137570 | Public directional eGFR support | GEO |
| Nephroseq v5 | Precomputed GFR association support | Public web interface |

## Analysis scripts

| File | Purpose |
|---|---|
| `analysis/independent_validation/validate_gse254185.py` | Load TAL cells, donor pseudobulk, injury exposure, and preliminary within-donor rank |
| `analysis/independent_validation/validate_gse254185_four_units.py` | kNN neighborhood features, four-unit comparison, exact permutation, leave-one-donor checks |
| `analysis/independent_validation/validate_gse211785.py` | C_TAL/M_TAL donor pseudobulk, permutation tests, neighborhood features, within-donor rank |
| `analysis/independent_validation/build_validation_table.py` | Build Additional file 14 from validation result tables |

## Key parameters

| Parameter | Value |
|---|---|
| GSE254185 TAL cells | 9,912 |
| GSE254185 donors | 12 (7 normal, 5 obstructive nephropathy) |
| GSE211785 TAL cells | 29,847 |
| GSE211785 TAL-covered donors | 35 (16 control, 7 DKD, 12 HKD) |
| Neighborhood k | 30 |
| Donor-label permutations | 20,000 unless exact enumeration |
| Expression normalization | CP10K + log1p |
| Multiple testing | Benjamini-Hochberg |

## Result files

- `Additional file 14.csv`
- `GSE211785_pseudobulk.csv`
- `GSE211785_leave_one_out.csv`
- `GSE211785_within_injury.csv`
- `GSE211785_within_nhood.csv`
- `GSE211785_neighborhood_expression.csv`
- `GSE254185_pseudobulk_exact.csv`
- `GSE254185_pseudobulk_leave_one_out.csv`
- `GSE254185_within_injury_exact.csv`
- `GSE254185_within_nhood_exact.csv`
- `GSE254185_four_unit_comparison.csv`
- `GSE254185_neighborhood_abundance.csv`
- `GSE254185_neighborhood_expression.csv`
- `GSE254185_within_donor_rank.csv`

## Limitations

- GSE211785 contains C_TAL and M_TAL rather than the exact aTAL2/frTAL discovery labels.
- GSE254185 uses obstructive nephropathy rather than the KPMP disease mix.
- Neighborhood p values are descriptive because neighborhoods overlap.
- These are independent-cohort analyses, not experimental or causal validation.
