# Iron Deprivation Transcriptional Responses and Growth Arrest in *Mycobacterium tuberculosis*

[![DOI](https://img.shields.io/badge/DOI-10.3389/fmicb.2022.956602-blue)](https://doi.org/10.3389/fmicb.2022.956602)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains analysis code, input data references, and manuscript assets for the publication:

**Zavala-Urtecho et al. (2022)**  
*Iron deprivation enhances transcriptional responses to in vitro growth arrest of Mycobacterium tuberculosis.*  
Frontiers in Microbiology, 13:956602  
https://doi.org/10.3389/fmicb.2022.956602

## Project goals

- Keep the repository clean and organized.
- Improve reproducibility and traceability of analyses.
- Document data sources, methods, and execution workflow in English.

## Repository structure

- `analysis/` — analysis scripts, reusable functions, and notebooks.
- `data/` — standardized data layout (`raw/` and `processed/`).
- `docs/` — reproducibility, methods, and data documentation.
- `results/` — generated figures, tables, and logs.
- `manuscript/` — manuscript-related assets (main/supplementary/proofs).
- `geo-submission/` — GEO submission metadata files.

Legacy project assets are still available in:
- `Analyses/`
- `Manuscript/`
- `GEO_submission/`

## Quick start

1. Create the software environment:
   ```bash
   conda env create -f environment.yml
   conda activate iron-deprivation-r
   ```
2. Read the reproducibility guide:
   - [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md)
3. Run the workflow helper:
   ```bash
   make help
   ```

## Documentation

- General docs index: [`docs/README.md`](docs/README.md)
- Reproducibility steps: [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md)
- Data inventory and provenance: [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md)
- Methodological details: [`docs/METHODS.md`](docs/METHODS.md)

## Citation

Please cite this work using [`CITATION.cff`](CITATION.cff).

## License

MIT License. See [`LICENSE`](LICENSE).
