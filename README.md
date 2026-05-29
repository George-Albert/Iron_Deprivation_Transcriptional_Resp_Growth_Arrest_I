# Iron Deprivation Enhances Transcriptional Responses to Growth Arrest

[![DOI](https://img.shields.io/badge/DOI-10.3389%2Ffmicb.2022.956602-blue)](https://doi.org/10.3389/fmicb.2022.956602)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.3.0-blue)](https://www.r-project.org/)

## Overview

This repository contains reproducible analysis code, data, and manuscript materials for the published paper:

**"Iron deprivation enhances transcriptional responses to in vitro growth arrest of *Mycobacterium tuberculosis*"**  
Alebouyeh S, Cárdenas-Pestana JA, Vazquez L, Prados-Rosales R, Del Portillo P, Sanz J, Menéndez MC, García MJ  
*Frontiers in Microbiology*, **2022**, 13:956602  
https://doi.org/10.3389/fmicb.2022.956602

## Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/George-Albert/Iron_Deprivation_Transcriptional_Resp_Growth_Arrest_I.git
cd Iron_Deprivation_Transcriptional_Resp_Growth_Arrest_I
```

### 2. Set Up Environment
```bash
# Install Conda environment
make install-env

# Activate environment
conda activate mtb-iron-deprivation
```

### 3. Validate and Run
```bash
# Validate repository structure
make validate

# Run all analyses
make run-legacy-scripts
```

## Repository Structure

```
├── docs/                          # Complete documentation
│   ├── README.md                  # Documentation index
│   ├── REPRODUCIBILITY.md         # Step-by-step reproduction guide
│   ├── DATA_SOURCES.md            # Data descriptions and sources
│   └── METHODS.md                 # Technical methodology details
├── data/
│   ├── raw/                       # Original unmodified data
│   └── processed/                 # Cleaned and processed datasets
├── analysis/
│   ├── scripts/                   # Main R analysis scripts
│   ├── functions/                 # Reusable R functions
│   └── notebooks/                 # RMarkdown exploratory notebooks
├── results/
│   ├── figures/                   # Generated plots
│   ├── tables/                    # Output data tables
│   └── logs/                      # Execution logs
├── manuscript/
│   ├── main/                      # Main figures
│   ├── supplementary/             # Supplementary materials
│   └── proofs/                    # Manuscript proofs
├── geo-submission/                # GEO database submission files
├── Analyses/                      # Legacy analysis directory (preserved)
│   ├── codes/                     # Original R scripts
│   ├── inputs/                    # Original input data
│   └── outputs/                   # Original outputs
├── .gitignore                     # Git exclusion rules
├── CITATION.cff                   # Citation metadata
├── environment.yml                # Conda environment
├── Makefile                       # Automation targets
├── LICENSE                        # MIT License
└── README.md                      # This file
```

## Documentation

For detailed information, visit:
- **[REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md)** - Complete reproduction guide
- **[DATA_SOURCES.md](docs/DATA_SOURCES.md)** - Data descriptions
- **[METHODS.md](docs/METHODS.md)** - Technical methodology
- **[docs/README.md](docs/README.md)** - Documentation index

## Key Findings

- ~714 genes show iron-dependent expression changes during growth arrest
- Stress and metal homeostasis genes more strongly upregulated without iron
- Energy metabolism genes more downregulated without iron
- Lipid remodeling (PDIM) appears independent of iron status

## Citation

```bibtex
@article{Alebouyeh2022,
  title={Iron deprivation enhances transcriptional responses to in vitro growth arrest of {Mycobacterium tuberculosis}},
  author={Alebouyeh, S and Cárdenas-Pestana, JA and Vazquez, L and Prados-Rosales, R and Del Portillo, P and Sanz, J and Menéndez, MC and García, MJ},
  journal={Frontiers in Microbiology},
  volume={13},
  pages={956602},
  year={2022},
  doi={10.3389/fmicb.2022.956602}
}
```

## Make Targets

```bash
make help              # Show all targets
make install-env      # Create conda environment
make validate         # Validate directory structure
make run-legacy-scripts # Run all R scripts
make clean            # Remove artifacts
```

## License

MIT License - See [LICENSE](LICENSE) file

## Contact

For questions or issues, see [docs/README.md](docs/README.md) or open an issue.

---

**Paper DOI**: [10.3389/fmicb.2022.956602](https://doi.org/10.3389/fmicb.2022.956602)
