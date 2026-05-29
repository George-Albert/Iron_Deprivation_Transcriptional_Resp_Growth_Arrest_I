# Data sources

## Publication context

Data correspond to the analyses reported in:

Zavala-Urtecho et al. (2022), Frontiers in Microbiology, DOI: 10.3389/fmicb.2022.956602.

## Current repository data locations

### Legacy data layout (kept for backward compatibility)

- `Analyses/inputs/Raw_data.csv`: transcriptomic input table.
- `Analyses/inputs/DESeq2.xlsx`: differential expression summary.
- `Analyses/inputs/Lipid_analysis.xlsx`: lipidomic summary data.
- `GEO_submission/GEO_metadata.xlsx`: GEO metadata file.
- `GEO_submission/Read_counts.xlsx`: read count matrix for submission.

### Standardized data layout (new)

- `data/raw/`: immutable source data.
- `data/processed/`: transformed data and analysis-ready tables.

## Data handling principles

- Do not overwrite raw files.
- Track preprocessing scripts in `analysis/scripts/`.
- Keep provenance notes for derived files in commit messages and docs.
