# Reproducibility Guide

This project uses Python/PyMC for the non-spatial Bayesian models and R-INLA for
the spatial extensions. The commands below are meant to make a clean checkout
easy to verify before running the longer model fits.

## 1. Python Environment

Create and activate a Python environment, then install:

```bash
pip install -r requirements.txt
```

The Python scripts expect `data/complete_combined_datasets.csv`. Scripts that
save summaries or figures write them to `outputs/`.

## 2. R Environment

Install the R dependencies with:

```bash
Rscript models/r_inla/spatial/install_packages.R
```

The spatial scripts use INLA, Arrow, data.table, Matrix, sf, and ggplot2.

## 3. Sanity Checks

Run syntax checks without fitting the models:

```bash
make check
```

## 4. Main Workflow

Common targets:

```bash
make m5
make s2
make s6
make diagnostics
make maps
```

These Make targets are convenience wrappers around the scripts listed in the
README. The full model suite can also be run script by script. Runtime can be
substantial because the models are fit directly.

## 5. Data

Data provenance is documented in `DATA_AVAILABILITY.md`. Source licenses and
redistribution permissions should be checked before releasing `data/` publicly.
