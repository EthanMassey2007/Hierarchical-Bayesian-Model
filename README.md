# Hierarchical Bayesian Dengue Model

This repository contains a research workflow for modeling weekly dengue case counts across municipalities in Rio de Janeiro, Brazil. The core model is a hierarchical Bayesian negative-binomial regression built with PyMC. It combines epidemiological, weather, socioeconomic, and transportation-connectivity covariates to estimate dengue case dynamics and compare different lag assumptions for prior case counts.

The project is intended for exploratory modeling, lag sensitivity analysis, and train/test evaluation of municipal dengue risk signals.

## Project Overview

The modeling pipeline:

1. Loads weekly municipality-level dengue, weather, and socioeconomic data.
2. Normalizes municipality names and joins them to IBGE municipality codes.
3. Adds transportation/connectivity features from air, road, and fluvial datasets.
4. Creates lagged case-count features by municipality.
5. Standardizes covariates.
6. Fits a hierarchical Bayesian count model with municipality, week, and year effects.
7. Reports posterior diagnostics and predictive metrics.
8. Optionally plots lag comparisons and actual-vs-predicted results.

The main target is weekly `cases` for municipalities in Rio de Janeiro (`TARGET_STATE = "RJ"`).

## Repository Structure

```text
.
|-- README.md
|-- lag_sweep_model_fixed.py
|-- year_split_model.py
|-- parse_data.py
|-- check_zeros.py
|-- lag_sweep_results.csv
`-- data/
    |-- cases.csv
    |-- temperature.csv
    |-- humidity.csv
    |-- rainfall.csv
    |-- idhm.csv
    |-- municipios.csv
    |-- aero_anac_2017_2023.parquet
    |-- fluvi_road_ibge.parquet
    |-- adjacency_matrix_correct.parquet
    |-- complete_combined_datasets.csv
    `-- hub_pop_density.csv
```

## Main Scripts

| Script | Purpose |
| --- | --- |
| `lag_sweep_model_fixed.py` | Fits the hierarchical Bayesian model across one or more lag values on the full configured time period. Saves lag metrics to `lag_sweep_results.csv`. |
| `year_split_model.py` | Fits on a training period and evaluates on a held-out test period. Saves metrics to `lag_sweep_train_test_results.csv` when enabled. |
| `parse_data.py` | Earlier end-to-end model/preprocessing script for building the merged modeling dataset and fitting a single model configuration. |
| `check_zeros.py` | Small utility for checking zero-case rows for Rio de Janeiro in a local case dataset. |

## Data Inputs

The modeling scripts expect the following files under `data/`:

| File | Expected role |
| --- | --- |
| `cases.csv` | Weekly dengue case counts with `municipio`, `year`, `week`, and `cases`. |
| `temperature.csv` | Weekly temperature by municipality. |
| `humidity.csv` | Weekly humidity by municipality. |
| `rainfall.csv` | Weekly rainfall by municipality. |
| `idhm.csv` | Weekly or repeated socioeconomic IDHM values by municipality. |
| `municipios.csv` | Municipality metadata used to map names to IBGE codes and filter Rio de Janeiro. |
| `aero_anac_2017_2023.parquet` | Air passenger/connectivity data. Used by `lag_sweep_model_fixed.py` and optionally by `year_split_model.py`. |
| `fluvi_road_ibge.parquet` | Road and fluvial connectivity data keyed by IBGE municipality codes. |

CSV column names are normalized to lowercase in the scripts. Municipality names are also lowercased, accent-stripped, and cleaned before joins.

## Model Summary

The main PyMC model estimates dengue counts with a negative-binomial likelihood:

```text
cases ~ NegativeBinomial(mu, alpha)
log(mu) = intercept
          + municipality random intercept
          + week-of-year effect
          + year effect
          + standardized covariates
```

Typical covariates include:

- Rainfall
- Humidity
- Temperature
- IDHM
- Lagged log case counts
- Air passenger/connectivity features
- Road connectivity features
- Fluvial connectivity features
- Optional time index in the train/test model

The scripts report MAE, RMSE, WAPE, WAPE-based accuracy, R-squared, Rhat, bulk ESS, and tail ESS.

## Requirements

Use Python 3.10 or newer. The project depends on:

- `numpy`
- `pandas`
- `pyarrow`
- `pymc`
- `arviz`
- `matplotlib`
- `scikit-learn`

Install the dependencies with:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy pandas pyarrow pymc arviz matplotlib scikit-learn
```

For reproducible project use, consider saving the environment after installation:

```bash
python -m pip freeze > requirements.txt
```

## Quick Start

From the repository root:

```bash
source .venv/bin/activate
python lag_sweep_model_fixed.py
```

This will:

- Build the merged modeling dataframe from `data/`.
- Fit the configured lag values.
- Print model diagnostics and predictive metrics.
- Save `lag_sweep_results.csv` when `SAVE_RESULTS_CSV = True`.
- Display plots when `MAKE_PLOTS = True`.

To run the train/test year split workflow:

```bash
python year_split_model.py
```

By default, `year_split_model.py` trains on 2017-2022 and tests on 2023-2025.

## Configuration

The main settings live near the top of each model script.

Common options:

| Setting | Description |
| --- | --- |
| `TARGET_STATE` | State filter for municipalities. Currently `RJ`. |
| `START_YEAR`, `END_YEAR` | Year range for the full-data lag sweep. |
| `TRAIN_START_YEAR`, `TRAIN_END_YEAR` | Training years in `year_split_model.py`. |
| `TEST_START_YEAR`, `TEST_END_YEAR` | Held-out test years in `year_split_model.py`. |
| `START_LAG`, `END_LAG` | Inclusive lag range to test. |
| `MAKE_PLOTS` | Whether to show Matplotlib figures. |
| `USE_FULL_COVARIATE_SET` | Whether to include the expanded transport/connectivity covariates. |
| `APPLY_LOG1P_TO_SKEWED_FEATURES` | Whether to log-transform skewed connectivity features. |
| `DRAWS`, `TUNE`, `CHAINS`, `CORES` | PyMC sampling controls. |
| `TARGET_ACCEPT` | NUTS target acceptance rate. |
| `RANDOM_SEED` | Sampling seed for repeatability. |

For faster smoke tests, reduce the sampling settings, for example:

```python
DRAWS = 200
TUNE = 500
CHAINS = 2
CORES = 2
MAKE_PLOTS = False
```

## Outputs

`lag_sweep_model_fixed.py` writes:

- `lag_sweep_results.csv`

The current tracked result file contains one configured lag:

| lag | accuracy_pct | mae | rmse | wape | r2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2 | 55.90 | 3.11 | 10.95 | 0.441 | 0.933 |

`year_split_model.py` writes:

- `lag_sweep_train_test_results.csv`

when `SAVE_RESULTS_CSV = True`.

Both main workflows print posterior summaries and diagnostics to the console.

## Metric Notes

The reported `accuracy_pct` is based on WAPE:

```text
accuracy_pct = 100 * (1 - sum(abs(actual - predicted)) / sum(actual))
```

This is not classification accuracy. It is a normalized aggregate error score where higher is better and 100% would be a perfect prediction under this metric.

## Reproducibility Notes

- PyMC sampling can take a long time with the default settings: `DRAWS = 1000`, `TUNE = 2000`, `CHAINS = 4`.
- The models use `RANDOM_SEED = 42`, but exact results can still vary across operating systems, dependency versions, BLAS backends, and sampler behavior.
- Some plots are displayed interactively with `plt.show()` and are not automatically saved as image files.
- Large data files are stored directly in `data/`. If publishing this repository publicly, confirm that the datasets can be redistributed.

## Suggested GitHub Additions

For a more production-ready public repository, consider adding:

- `requirements.txt` or `environment.yml`
- `.gitignore` for virtual environments, Python caches, and generated outputs
- A data license or citation section
- A project license such as MIT, Apache-2.0, or GPL
- Saved example plots under a `figures/` directory

## License

No license file is currently included. Add a license before distributing or reusing the project publicly.
