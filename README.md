# Hierarchical Bayesian Dengue Models

This repository contains a model-comparison workflow for weekly dengue case counts across municipalities in Rio de Janeiro, Brazil. The project starts with non-spatial hierarchical Bayesian negative-binomial models in Python/PyMC and extends the strongest non-spatial model into spatial R-INLA models.

The main research goal is to evaluate whether rainfall and other environmental, socioeconomic, temporal, and spatial features improve dengue case modeling while preserving interpretability.

## Current Modeling Strategy

The model family is organized in two stages:

1. M0-M6: non-spatial hierarchical Bayesian models fit in Python/PyMC.
2. S1-S5: spatial extensions fit in R-INLA.

The current recommended main spatial model is:

```text
S2 = M5 + BYM2 spatial random effect + adjacency-based lagged neighboring cases
```

S3-S5 are best treated as sensitivity or extension models. They test alternative spatial spread and mobility assumptions, but they have not clearly replaced S2 as the clean main model.

## Data Inputs

Primary modeling dataset:

```text
data/complete_combined_datasets.csv
```

Spatial and mobility support files:

| File | Role |
| --- | --- |
| `data/adjacency_matrix_correct.parquet` | Municipality adjacency graph for BYM2 and adjacency-based spatial lag. |
| `data/fluvi_road_ibge.parquet` | Road/fluvial connectivity and municipality coordinate support. |
| `data/aero_anac_2017_2023.parquet` | Air passenger mobility support. |
| `data/hub_pop_density.csv` | Municipality lookup and population-density support data. |
| `data/municipios.csv` | Municipality metadata. |

Because the model scripts are now organized into `base_model/` and `spatial_R/`, the shared data directory remains one level above those script folders:

```text
project root/
|-- data/
|-- base_model/
`-- spatial_R/
```

If a script reports that it cannot find `complete_combined_datasets.csv`, check that its `DATA_DIR` points to the project-root `data/` folder rather than a nested `base_model/data/` or `spatial_R/data/` folder.

Expected core columns in `complete_combined_datasets.csv`:

| Column | Role |
| --- | --- |
| `municipio` | Municipality name. |
| `year` | ISO year. |
| `week` | ISO week number. |
| `cases` | Weekly dengue case count. |
| `rainfall` | Weekly rainfall. |
| `humidity` | Weekly humidity. |
| `temperature` | Weekly temperature. |
| `idhm` | Municipal Human Development Index. |

## Non-Spatial Models

| Model | Script | Description | Main covariates |
| --- | --- | --- | --- |
| M0 | `base_model/base_model.py` | True base/null model | None |
| M1 | `base_model/base_model_covariates.py` | Base + same-week covariates | `rainfall`, `humidity`, `temperature`, `idhm` |
| M2 | `base_model/base_model_lag_weather.py` | Base + lagged covariates | `rainfall_lag`, `humidity_lag`, `temperature_lag`, `idhm` |
| M3 | `base_model/base_model_interpolation.py` | Base + covariates + interpolation | `rainfall`, interpolated `humidity`, interpolated `temperature`, `idhm` |
| M4 | `base_model/base_model_lag_cases.py` | Base + covariates + lagged log cases | `rainfall`, `humidity`, `temperature`, `idhm`, `log_cases_lag` |
| M5 | `base_model/base_model_lag_cases_weather.py` | Base + lagged covariates + lagged log cases | `rainfall_lag`, `humidity_lag`, `temperature_lag`, `idhm`, `log_cases_lag` |
| M6 | `base_model/base_model_lag_cases_weather_interpolation.py` | Base + lagged covariates + lagged log cases + interpolation | `rainfall_lag`, interpolated `humidity_lag`, interpolated `temperature_lag`, `idhm`, `log_cases_lag` |

M5 is the strongest current non-spatial model because lagged own cases substantially improve predictive performance while still allowing interpretation of lagged weather effects.

## Spatial INLA Models

| Model | Script | Description | Added spatial/mobility term |
| --- | --- | --- | --- |
| S1 | `spatial_R/spatial_inla_model_s1.R` | M5 + BYM2 spatial random effect | BYM2 structured/unstructured municipality effect |
| S2 | `spatial_R/spatial_inla_model_s2.R` | S1 + adjacency-based neighboring case lag | `neighbor_log_cases_lag` |
| S3 | `spatial_R/spatial_inla_model_s3.R` | S1 + distance-weighted neighboring case lag | `distance_log_cases_lag` |
| S4 | `spatial_R/spatial_inla_model_s4_road.R` | S2 + road connectivity | `road_connectivity_log` |
| S5 | `spatial_R/spatial_inla_model_s5_air.R` | S2 + air passenger mobility | `air_passenger_week_lag_log` |

### Spatial Model Definitions

S1 adds BYM2 spatial smoothing so adjacent municipalities can share unexplained baseline risk.

S2 adds adjacency-based neighboring dengue pressure:

```text
neighbor_cases_lag[i,t] = mean cases among adjacent municipalities at t - CASE_LAG_WEEKS
neighbor_log_cases_lag[i,t] = log1p(neighbor_cases_lag[i,t])
```

S3 replaces the adjacency-only neighboring case lag with a distance-weighted lag:

```text
distance_cases_lag[i,t] = sum_j w_ij * cases[j,t-4] / sum_j w_ij
w_ij = exp(-distance_ij / 50)
```

Only municipalities within the distance cutoff are included, and each lag value must come from a date strictly before the prediction date.

S4 adds road connectivity:

```text
road_connectivity_log[i] = log1p(sum incident road connectivity for municipality i)
```

S5 adds previous-month air mobility:

```text
air_passenger_week_lag_raw[i,t]
= inbound + outbound aero_pass_week for municipality i during the previous calendar month

air_passenger_week_lag_log[i,t] = log1p(air_passenger_week_lag_raw[i,t])
```

S5 intentionally uses previous-month air mobility instead of same-year annual totals to avoid future-covariate leakage in the 2023 test period.

## Model Form

The Python models use a negative-binomial likelihood:

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

The spatial INLA models use the same count-model idea but add BYM2 spatial structure and, where relevant, spatial lag or mobility covariates:

```text
log(mu_it) = intercept
             + fixed covariate effects
             + BYM2 municipality spatial effect
```

All continuous covariates are standardized before fitting. Coefficients are interpreted as the change in log expected cases for a one-standard-deviation increase in the covariate, holding other terms constant.

For example:

```text
beta[rainfall_lag] = 0.204
exp(0.204) = 1.226
```

This corresponds to about a 22.6% increase in expected weekly cases for a one-standard-deviation increase in lagged rainfall, conditional on the rest of the model.

## Lags And Leakage Prevention

The core lag settings are defined near the top of the scripts:

```text
CASE_LAG_WEEKS = 4
WEATHER_LAG_WEEKS = 6
```

Lagged cases:

```text
log_cases_lag = log1p(cases from exactly CASE_LAG_WEEKS earlier)
```

Lagged weather:

```text
rainfall_lag    = rainfall from exactly WEATHER_LAG_WEEKS earlier
humidity_lag    = humidity from exactly WEATHER_LAG_WEEKS earlier
temperature_lag = temperature from exactly WEATHER_LAG_WEEKS earlier
```

Leakage checks used in the model scripts:

- Own-case lags must come from dates strictly before the current row.
- Neighbor-case lags must come from dates strictly before the current row.
- Test rows are dropped if their own-case or neighbor-case lag would come from inside the held-out test period.
- Scalers are fit on training data only.
- S5 air mobility uses previous-month data only.
- Interpolation models avoid using training-period held-out values during training preprocessing.

## Missing Data And Interpolation

All models drop rows missing required core fields:

```text
municipio, year, week, cases
```

Models without interpolation drop rows missing their required covariates.

Interpolation models fill missing humidity and temperature values using municipality-level linear temporal interpolation:

```text
x(t) = x(t0) + ((t - t0) / (t1 - t0)) * (x(t1) - x(t0))
```

where `t0` and `t1` are observed dates surrounding the missing value. Interpolation is limited to internal gaps of at most 8 weeks:

```text
INTERPOLATION_LIMIT_WEEKS = 8
```

Rainfall is not interpolated because it is assumed to be complete in the current dataset.

## Train/Test Evaluation

The main evaluation split is:

```text
TRAIN_START_YEAR = 2017
TRAIN_END_YEAR   = 2022
TEST_START_YEAR  = 2023
TEST_END_YEAR    = 2023
```

Metrics reported by the scripts include:

- MAE
- RMSE
- WAPE
- WAPE-based accuracy
- R-squared
- DIC and WAIC for INLA full-data fits
- Posterior means and 95% credible intervals

The WAPE-based accuracy is:

```text
accuracy_pct = 100 * (1 - sum(abs(actual - predicted)) / sum(actual))
```

This is not classification accuracy. It is a normalized aggregate count-error metric where higher is better.

## Current Results Snapshot

These are the current observed results from the most recent model runs. Re-run the scripts before reporting final paper numbers.

| Model | Test accuracy | Test WAPE | Test R2 | Notes |
| --- | ---: | ---: | ---: | --- |
| S1 | about 51.6% | not recorded here | not recorded here | BYM2 spatial effect only. |
| S2 | 53.257% | 0.4674 | 0.8684 | Recommended main spatial model so far. |
| S3 | 51.430% | 0.4857 | 0.8492 | Distance-weighted lag did not beat S2. |
| S4-road | 53.435% | 0.4657 | 0.8694 | Small road-mobility gain; road interval crossed zero. |
| S5-air | 53.431% | 0.4657 | 0.8706 | Small air-mobility gain; train-fit interval nearly crossed zero. |

Current full-data S2 fixed effects:

| Term | Posterior mean |
| --- | ---: |
| `rainfall_lag_z` | 0.2076 |
| `humidity_lag_z` | -0.0216 |
| `temperature_lag_z` | 0.0826 |
| `idhm_z` | 0.5238 |
| `log_cases_lag_z` | 0.9367 |
| `neighbor_log_cases_lag_z` | 0.2771 |

Current full-data mobility effects:

| Model | Term | Posterior mean | 95% interval summary |
| --- | --- | ---: | --- |
| S4-road | `road_connectivity_log_z` | 0.2089 | Crossed zero. |
| S5-air | `air_passenger_week_lag_log_z` | 0.0363 | Positive in full-data fit; nearly crossed zero in train fit. |

Interpretation:

- Lagged own cases are the strongest predictor.
- Neighboring lagged cases add meaningful spatial transmission signal.
- Rainfall remains positive and interpretable after accounting for lagged cases and spatial structure.
- Road and air mobility provide only small predictive changes and are best treated as sensitivity/extension models.

## Repository Structure

```text
.
|-- README.md
|-- base_model/
|   |-- base_model.py
|   |-- base_model_covariates.py
|   |-- base_model_lag_weather.py
|   |-- base_model_interpolation.py
|   |-- base_model_lag_cases.py
|   |-- base_model_lag_cases_weather.py
|   `-- base_model_lag_cases_weather_interpolation.py
|-- spatial_R/
|   |-- spatial_inla_model_s1.R
|   |-- spatial_inla_model_s2.R
|   |-- spatial_inla_model_s3.R
|   |-- spatial_inla_model_s4_road.R
|   `-- spatial_inla_model_s5_air.R
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
| `base_model/base_model.py` | M0 true baseline with no covariates. |
| `base_model/base_model_covariates.py` | M1 same-week weather and IDHM covariate model. |
| `base_model/base_model_lag_weather.py` | M2 lagged-weather covariate model. |
| `base_model/base_model_interpolation.py` | M3 same-week covariate model with humidity/temperature interpolation. |
| `base_model/base_model_lag_cases.py` | M4 same-week covariates plus lagged log cases. |
| `base_model/base_model_lag_cases_weather.py` | M5 lagged weather plus lagged log cases. |
| `base_model/base_model_lag_cases_weather_interpolation.py` | M6 lagged weather plus lagged log cases with humidity/temperature interpolation. |
| `spatial_R/spatial_inla_model_s1.R` | S1 spatial BYM2 extension of M5. |
| `spatial_R/spatial_inla_model_s2.R` | S2 BYM2 + adjacency-based neighboring lagged cases. |
| `spatial_R/spatial_inla_model_s3.R` | S3 BYM2 + distance-weighted neighboring lagged cases. |
| `spatial_R/spatial_inla_model_s4_road.R` | S4-road BYM2 + adjacency lag + road connectivity. |
| `spatial_R/spatial_inla_model_s5_air.R` | S5-air BYM2 + adjacency lag + previous-month air mobility. |

## Python Quick Start

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy pandas pyarrow pymc arviz matplotlib scikit-learn
```

Run the true baseline:

```bash
python base_model/base_model.py
```

Run M5:

```bash
python base_model/base_model_lag_cases_weather.py
```

Run M6:

```bash
python base_model/base_model_lag_cases_weather_interpolation.py
```

## R-INLA Quick Start

From the project folder:

```bash
cd /Users/ethanmassey/VS_Code_Test/Hierarchical-Bayesian-Model
```

Run the recommended main spatial model:

```bash
Rscript spatial_R/spatial_inla_model_s2.R
```

Run S5-air:

```bash
Rscript spatial_R/spatial_inla_model_s5_air.R
```

Run only dataframe and leakage checks without fitting the model:

```bash
INLA_RUN_MODEL=0 Rscript spatial_R/spatial_inla_model_s2.R
```

The INLA scripts may print messages such as:

```text
sh: /bin/kstat: No such file or directory
```

These are usually harmless Mac/R environment messages and do not necessarily indicate model failure.

## Configuration

The main settings are defined near the top of each model script.

Common Python settings:

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
| `RUN_TRAIN_TEST_EVALUATION` | Whether to run held-out train/test evaluation. |

Common R-INLA settings:

| Setting | Description |
| --- | --- |
| `CASE_LAG_WEEKS` | Own-case and neighbor-case lag length. |
| `WEATHER_LAG_WEEKS` | Weather lag length. |
| `TRAIN_START_YEAR`, `TRAIN_END_YEAR` | Training years. |
| `TEST_START_YEAR`, `TEST_END_YEAR` | Held-out test years. |
| `INLA_RUN_MODEL` | Set to `0` to build/check data without fitting. |

For fast PyMC smoke tests, reduce:

```python
DRAWS = 100
TUNE = 200
CHAINS = 2
CORES = 2
MAKE_PLOTS = False
```

For final research results, use larger sampling settings and report convergence diagnostics.

## Publication Notes

For a research paper, the strongest current structure is:

```text
M5 as the main non-spatial model
S2 as the main spatial model
S3-S5 as sensitivity or extension models
```

Recommended tables and figures:

- Model comparison table for M0-M6 and S1-S5.
- Posterior coefficient table for the final model.
- Observed-vs-predicted plot for the held-out year.
- Residuals by year and municipality.
- Map of BYM2 spatial random effects.
- Map of prediction error or residual burden.
- Sensitivity table showing whether rainfall changes under S3-S5.

The current WAPE-based accuracy around 53% should be described as moderate predictive performance, not classification accuracy. The model is best framed as an interpretable spatial epidemiological model rather than a purely optimized machine-learning forecaster.

## Reproducibility Notes

- Python scripts use `RANDOM_SEED = 42`.
- Exact Bayesian sampling results may vary across machines and package versions.
- Rows included in each model may differ if interpolation fills values that another model drops.
- When comparing predictive metrics, confirm whether models are evaluated on the same held-out rows.
- Generated CSV outputs are disabled by default in most Python model scripts with `SAVE_OUTPUTS = False`.
- INLA output may vary slightly across R-INLA versions.

## License

No license file is currently included. Add a license before distributing or reusing this project publicly.
