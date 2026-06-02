# Hierarchical Bayesian Dengue Models

This repository contains a model-comparison workflow for weekly dengue case counts across municipalities in Rio de Janeiro, Brazil. The models are hierarchical Bayesian negative-binomial regressions built with PyMC and evaluated with a held-out year split.

The current workflow uses a combined dataset:

```text
data/complete_combined_datasets.csv
```

The model family is organized from a true baseline model through covariate, lagged-weather, interpolation, and lagged-case extensions.

## Research Goal

The main goal is to evaluate whether environmental and socioeconomic covariates, especially rainfall, improve dengue case modeling after accounting for municipality-level differences, seasonality, annual effects, and recent outbreak momentum.

The models are designed to answer questions such as:

- How much predictive signal comes from municipality, week, and year effects alone?
- Do same-week weather and IDHM covariates improve the model?
- Do lagged weather covariates improve the model?
- Does interpolation of missing humidity and temperature values change results?
- How much predictive power is added by lagged dengue cases?
- Does lagged rainfall remain important after lagged cases are included?

## Model Lineup

| Model | Script | Description | Covariates |
| --- | --- | --- | --- |
| M0 | `base_model.py` | Base/null hierarchical model | None |
| M1 | `base_model_covariates.py` | Base + same-week covariates | `rainfall`, `humidity`, `temperature`, `idhm` |
| M2 | `base_model_lag_weather.py` | Base + lagged covariates | `rainfall_lag`, `humidity_lag`, `temperature_lag`, `idhm` |
| M3 | `base_model_interpolation.py` | Base + covariates + interpolation | `rainfall`, interpolated `humidity`, interpolated `temperature`, `idhm` |
| M4 | `base_model_lag_cases.py` | Base + covariates + lagged log cases | `rainfall`, `humidity`, `temperature`, `idhm`, `log_cases_lag` |
| M5 | `base_model_lag_cases_weather.py` | Base + lagged covariates + lagged log cases | `rainfall_lag`, `humidity_lag`, `temperature_lag`, `idhm`, `log_cases_lag` |
| M6 | `base_model_lag_cases_weather_interpolation.py` | Base + lagged covariates + lagged log cases + interpolation | `rainfall_lag`, interpolated `humidity_lag`, interpolated `temperature_lag`, `idhm`, `log_cases_lag` |

## Model Structure

All M0-M6 models use a negative-binomial likelihood for weekly dengue counts:

```text
cases_it ~ NegativeBinomial(mu_it, alpha)
```

with a log link:

```text
log(mu_it) = intercept
             + municipality effect
             + week-of-year effect
             + year effect
             + covariate effects
```

M0 excludes the covariate effects and functions as the true null baseline. M1-M6 add covariates or feature-engineering steps to test whether those additions improve prediction or change posterior covariate effects.

The hierarchical terms are:

- Municipality effects: persistent differences between municipalities.
- Week effects: recurring seasonal patterns.
- Year effects: annual outbreak intensity or reporting differences.

Covariates are standardized before model fitting, so beta coefficients are interpreted per one-standard-deviation increase in the covariate.

## Data Inputs

Primary input:

| File | Role |
| --- | --- |
| `data/complete_combined_datasets.csv` | Main weekly municipality-level modeling dataset. |

Expected core columns:

| Column | Role |
| --- | --- |
| `municipio` | Municipality name. |
| `year` | ISO year. |
| `week` | ISO week number. |
| `cases` | Weekly dengue case count. |
| `rainfall` | Weekly rainfall. |
| `humidity` | Weekly humidity. |
| `temperature` | Weekly temperature. |
| `idhm` | Municipal Human Development Index covariate. |

Additional data files used by support scripts:

| File | Role |
| --- | --- |
| `data/hub_pop_density.csv` | Municipality lookup and population-density support data. |
| `data/aero_anac_2017_2023.parquet` | Air transportation data for expanded/support workflows. |
| `data/fluvi_road_ibge.parquet` | Road and fluvial connectivity data for expanded/support workflows. |
| `data/adjacency_matrix_correct.parquet` | Spatial adjacency data for future spatial modeling. |
| `data/municipios.csv` | Municipality metadata. |

## Missing-Data Handling

All models drop rows missing required core fields:

```text
municipio, year, week, cases
```

Models without interpolation drop rows missing their required covariates.

Interpolation models fill missing `humidity` and `temperature` values using municipality-level linear temporal interpolation:

```text
x(t) = x(t0) + ((t - t0) / (t1 - t0)) * (x(t1) - x(t0))
```

where `t0` and `t1` are observed dates surrounding the missing value. Interpolation is limited to internal gaps of at most 8 weeks:

```text
INTERPOLATION_LIMIT_WEEKS = 8
```

Rainfall is not interpolated because it is assumed to be complete in the current dataset.

For train/test evaluation, interpolation is applied separately to the training and testing periods to avoid using held-out data during training preprocessing.

## Lag Definitions

The lagged-case models use:

```text
CASE_LAG_WEEKS = 4
log_cases_lag = log1p(cases from exactly 4 calendar weeks earlier)
```

The lagged-weather models use:

```text
WEATHER_LAG_WEEKS = 6
rainfall_lag     = rainfall from exactly 6 calendar weeks earlier
humidity_lag     = humidity from exactly 6 calendar weeks earlier
temperature_lag  = temperature from exactly 6 calendar weeks earlier
```

Lag construction is municipality-specific and calendar-date based. The scripts validate that lagged values come from strictly earlier dates.

For held-out test evaluation, rows are removed if their lagged case value would come from inside the held-out test period. This prevents leakage from future test outcomes into test predictions.

## Train/Test Evaluation

The main model scripts use:

```text
TRAIN_START_YEAR = 2017
TRAIN_END_YEAR   = 2022
TEST_START_YEAR  = 2023
TEST_END_YEAR    = 2023
```

The training scaler and categorical encoders are fit only on the training data. The held-out test year is not assigned a learned year effect; unseen test years use a zero year-effect contribution during prediction.

Metrics reported by the scripts include:

- MAE
- RMSE
- WAPE
- WAPE-based accuracy
- R-squared
- R-hat
- Bulk ESS
- Tail ESS

The WAPE-based accuracy is:

```text
accuracy_pct = 100 * (1 - sum(abs(actual - predicted)) / sum(actual))
```

This is not classification accuracy. It is a normalized aggregate count-error metric where higher is better.

## Repository Structure

```text
.
|-- README.md
|-- base_model.py
|-- base_model_covariates.py
|-- base_model_lag_weather.py
|-- base_model_interpolation.py
|-- base_model_lag_cases.py
|-- base_model_lag_cases_weather.py
|-- base_model_lag_cases_weather_interpolation.py
|-- lag_sweep_model_fixed.py
|-- year_split_model.py
|-- parse_data.py
|-- check_zeros.py
|-- lag_sweep_results.csv
`-- data/
    |-- complete_combined_datasets.csv
    |-- hub_pop_density.csv
    |-- aero_anac_2017_2023.parquet
    |-- fluvi_road_ibge.parquet
    |-- adjacency_matrix_correct.parquet
    `-- municipios.csv
```

## Script Guide

| Script | Purpose |
| --- | --- |
| `base_model.py` | M0 true baseline with no covariates. |
| `base_model_covariates.py` | M1 same-week weather and IDHM covariate model. |
| `base_model_lag_weather.py` | M2 lagged-weather covariate model. |
| `base_model_interpolation.py` | M3 same-week covariate model with humidity/temperature interpolation. |
| `base_model_lag_cases.py` | M4 same-week covariates plus lagged log cases. |
| `base_model_lag_cases_weather.py` | M5 lagged weather plus lagged log cases. |
| `base_model_lag_cases_weather_interpolation.py` | M6 lagged weather plus lagged log cases with humidity/temperature interpolation. |
| `lag_sweep_model_fixed.py` | Lag-sweep support script using the combined dataset. |
| `year_split_model.py` | Support script for train/test year-split modeling. |
| `parse_data.py` | Support preprocessing/modeling script using the combined dataset. |
| `check_zeros.py` | Utility for inspecting zero-case rows. |

## Quick Start

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy pandas pyarrow pymc arviz matplotlib scikit-learn
```

Run the true baseline:

```bash
python base_model.py
```

Run the full lagged-cases, lagged-weather, interpolation model:

```bash
python base_model_lag_cases_weather_interpolation.py
```

Run the lag-sweep support workflow:

```bash
python lag_sweep_model_fixed.py
```

## Configuration

The main settings are defined near the top of each model script.

Common settings:

| Setting | Description |
| --- | --- |
| `DATA_START_YEAR`, `DATA_END_YEAR` | Full modeling year range. |
| `TRAIN_START_YEAR`, `TRAIN_END_YEAR` | Training years for held-out evaluation. |
| `TEST_START_YEAR`, `TEST_END_YEAR` | Test years for held-out evaluation. |
| `CASE_LAG_WEEKS` | Number of weeks used for lagged cases. |
| `WEATHER_LAG_WEEKS` | Number of weeks used for lagged weather. |
| `INTERPOLATION_LIMIT_WEEKS` | Maximum internal gap length for interpolation models. |
| `DRAWS`, `TUNE`, `CHAINS`, `CORES` | PyMC sampling settings. |
| `TARGET_ACCEPT` | NUTS target acceptance rate. |
| `MAKE_PLOTS` | Whether to show actual-vs-predicted plots. |
| `SAVE_OUTPUTS` | Whether to save metrics and predictions as CSV files. |
| `RUN_TRAIN_TEST_EVALUATION` | Whether to run the held-out train/test evaluation after the full-data fit. |

PyMC sampling can take a long time. For a fast smoke test, reduce:

```python
DRAWS = 100
TUNE = 200
CHAINS = 2
CORES = 2
MAKE_PLOTS = False
```

For final research results, use larger sampling settings and report convergence diagnostics.

## Interpreting Betas

Because covariates are standardized and the model uses a log link, beta coefficients are interpreted on the log expected-count scale.

For example:

```text
beta[rainfall_lag] = 0.204
exp(0.204) = 1.226
```

This corresponds to about a 22.6% increase in expected weekly cases for a one-standard-deviation increase in lagged rainfall, holding the other model terms constant.

Posterior uncertainty should be assessed with the HDI, R-hat, and ESS. A covariate mean alone is not enough to determine whether an effect is reliable.

## Reproducibility Notes

- The scripts use `RANDOM_SEED = 42`.
- Exact Bayesian sampling results may vary across machines and package versions.
- Rows included in each model may differ if interpolation fills values that another model drops.
- When comparing predictive metrics, confirm whether models are being evaluated on the same held-out rows.
- Generated CSV outputs are disabled by default in most model scripts with `SAVE_OUTPUTS = False`.

## Future Work

Possible extensions include:

- Spatial random effects using municipality adjacency.
- BYM2/CAR spatial models, potentially fit with R-INLA.
- Rainfall nonlinearities, such as quadratic rainfall terms.
- Rainfall interactions with IDHM or municipality vulnerability groups.
- Posterior predictive interval coverage comparisons across M0-M6.

## License

No license file is currently included. Add a license before distributing or reusing this project publicly.
