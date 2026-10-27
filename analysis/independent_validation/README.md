# Independent non-KPMP validation

These scripts reproduce the independent analyses added in release `v1.0.4`.
They use public data and do not require access to the controlled KPMP source
files.

## Data

### GSE254185

Download the following public files from the GEO supplementary directory:

```text
GSE254185_cell_meta_data.tsv.gz
GSE254185_Library1_gex_matrix.h5
...
GSE254185_Library6_gex_matrix.h5
```

Put them in a directory and set:

```text
PAPPA2_VALIDATION_DIR=<directory containing the six GEX matrices and metadata>
```

Run:

```text
python validate_gse254185.py
python validate_gse254185_four_units.py
```

### GSE211785

Download the non-KPMP `GSE211785_Susztak_SC_SN_ATAC_merged_PreSCVI_final.h5ad.gz`
file from GEO and decompress it. Do not use the file whose name contains
`KPMP`.

Set:

```text
GSE211785_H5AD=<path to the decompressed PreSCVI h5ad>
```

Run:

```text
python validate_gse211785.py
```

The combined Additional file 14 can then be rebuilt with:

```text
python build_validation_table.py
```

## Interpretation boundaries

- GSE211785 supports donor-level PAPPA2 disease elevation in C_TAL and M_TAL.
- GSE211785 does not use the aTAL2/frTAL labels from the KPMP discovery set.
- GSE254185 provides a qualitative cross-unit signal in injured TAL, but its
  disease context and state labels differ from KPMP.
- Neighborhood-expression p values are descriptive because kNN neighborhoods
  overlap.
- These analyses are independent-cohort validation, not experimental or causal
  validation.
