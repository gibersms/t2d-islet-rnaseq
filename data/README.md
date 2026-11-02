# Data

The data are not stored in this repository. Run `scripts/01_download_data.R` to download them from NCBI GEO into `data/raw/`. The later scripts write intermediate files to `data/processed/`.

| File | Source |
|---|---|
| `GSE164416_raw_counts_GRCh38.p13_NCBI.tsv.gz` | NCBI-generated raw counts for GEO series GSE164416 |
| `Human.GRCh38.p13.annot.tsv.gz` | NCBI human gene annotation table used with those counts |
| `sample_info_all_columns.csv` | Sample characteristics from the GEO series matrix |
