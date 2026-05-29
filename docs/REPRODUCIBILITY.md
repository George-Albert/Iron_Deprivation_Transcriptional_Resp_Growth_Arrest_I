# Reproducibility guide

## Scope

This document describes how to reproduce the analysis workflow associated with:

Zavala-Urtecho et al. (2022), Frontiers in Microbiology, DOI: 10.3389/fmicb.2022.956602.

## 1. Environment setup

1. Install Conda (Miniconda or Anaconda).
2. Create and activate the environment:

```bash
conda env create -f environment.yml
conda activate iron-deprivation-r
```

## 2. Inspect data inputs

- Standardized data layout:
  - `data/raw/`
  - `data/processed/`
- Legacy data used in the publication remain in:
  - `Analyses/inputs/`
  - `GEO_submission/`

See [`DATA_SOURCES.md`](DATA_SOURCES.md) for details.

## 3. Run analyses

Use `make help` to inspect available targets.

```bash
make validate
make run-legacy-scripts
```

The `run-legacy-scripts` target provides a draft execution flow based on figure-specific scripts under `Analyses/codes/`.

## 4. Collect outputs

Store generated artifacts in:

- `results/figures/`
- `results/tables/`
- `results/logs/`

## 5. Traceability checklist

- Record software versions (`R --version` and package versions).
- Keep raw data immutable.
- Keep transformed/intermediate files in `data/processed/`.
- Keep logs from each run in `results/logs/`.
