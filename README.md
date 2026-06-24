# Hierarchical Bayesian Dengue Modeling in Rio de Janeiro

The main research question is whether rainfall and related climate conditions help explain dengue cases, and whether those effects vary across space. The current strongest non-spatial model is M5. The current main spatial model is S2, with S6 used to test the key rainfall-by-region research question.

## Project Status

The models are organized into 3 layers:

1. **M0-M6:** non-spatial hierarchical Bayesian models fit in Python/PyMC.
2. **S1-S5:** spatial and mobility extensions fit in R-INLA.
3. **S6-S7:** region-varying climate-effect models fit in R-INLA.

Current interpretation:

- **M5** strongest non-spatial model because lagged cases add a major predictive signal.
- **S2** cleanest main spatial model because it adds BYM2 spatial structure and adjacency-based neighboring lagged cases.
- **S6** most important model for our question because it tests whether rainfall effects differ by region
- **S7** useful sensitivity model that tests whether temperature also varies by region.
- **S3-S5** are spatial or mobility sensitivity models. They are useful, but they have not replaced S2.

## Data

Primary modeling file:

```text
data/complete_combined_datasets.csv
```

Support files:

| File | Role |
| --- | --- |
| `data/adjacency_matrix_correct.parquet` | Municipality adjacency graph for BYM2 spatial structure and adjacency-based neighboring case lag. |
| `data/RJ.json` | GeoJSON boundary file for Rio de Janeiro municipality maps. |
| `data/municipios.csv` | Municipality metadata and IBGE lookup support. |
| `data/hub_pop_density.csv` | IBGE, mesoregion, and population-density support data. |
| `data/fluvi_road_ibge.parquet` | Road/fluvial connectivity support used in S4. |
| `data/aero_anac_2017_2023.parquet` | Air passenger mobility support used in S5. |

Expected core columns in `complete_combined_datasets.csv`:

| Column | Meaning |
| --- | --- |
| `municipio` | Municipality name. |
| `year` | Year. |
| `week` | Epidemiological/ISO week. |
| `cases` | Weekly dengue case count. |
| `rainfall` | Weekly rainfall covariate. |
| `humidity` | Weekly humidity covariate. |
| `temperature` | Weekly temperature covariate. |
| `idhm` | Municipal Human Development Index. |

Project layout:

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
|   |-- spatial_inla_model_s5_air.R
|   |-- spatial_inla_model_s6_rainfall_region.R
|   |-- spatial_inla_model_s7_temperature_region.R
|   |-- map_s2_unexplained_effects.R
|   |-- map_s6_rainfall_region_effect.R
|   |-- map_s7_temperature_region_effect.R
|   `-- s2_morans_i_diagnostic.R
|-- outputs/
|   |-- s2_unexplained_spatial_effects.csv
|   |-- s2_unexplained_spatial_relative_risk_map.png
|   |-- s2_morans_i_diagnostic.csv
|   |-- s6_rainfall_region_effects.csv
|   |-- s6_rainfall_region_effect_map.png
|   |-- s7_temperature_region_effects.csv
|   |-- s7_temperature_region_effect_map.png
|   `-- run_logs/
`-- data/
```

## Model Lineup

### Non-Spatial Models

| Model | Script | Model definition | Main purpose |
| --- | --- | --- | --- |
| M0 | `base_model/base_model.py` | Null hierarchical model with no covariates. | True baseline. |
| M1 | `base_model/base_model_covariates.py` | M0 + same-week rainfall, humidity, temperature, and IDHM. | Tests whether basic covariates help before adding lags. |
| M2 | `base_model/base_model_lag_weather.py` | M0 + lagged rainfall, humidity, temperature, and IDHM. | Tests delayed weather effects. |
| M3 | `base_model/base_model_interpolation.py` | M1 + leakage-free interpolation for humidity and temperature. | Tests whether filling missing climate values helps. |
| M4 | `base_model/base_model_lag_cases.py` | M1 + lagged log cases. | Tests temporal persistence in dengue cases. |
| M5 | `base_model/base_model_lag_cases_weather.py` | M0 + lagged weather + IDHM + lagged log cases. | Main non-spatial model. |
| M6 | `base_model/base_model_lag_cases_weather_interpolation.py` | M5 + leakage-free interpolation. | Tests interpolation after adding temporal structure. |

### Spatial and Mobility Models

| Model | Script | Model definition | Main purpose |
| --- | --- | --- | --- |
| S1 | `spatial_R/spatial_inla_model_s1.R` | M5 + BYM2 spatial random effect. | Tests unexplained spatial risk. |
| S2 | `spatial_R/spatial_inla_model_s2.R` | S1 + adjacency-based neighboring lagged cases. | Main spatial model. |
| S3 | `spatial_R/spatial_inla_model_s3.R` | S1 + distance-weighted neighboring lagged cases. | Sensitivity test for distance-based spread. |
| S4 | `spatial_R/spatial_inla_model_s4_road.R` | S2 + road connectivity. | Tests road mobility. |
| S5 | `spatial_R/spatial_inla_model_s5_air.R` | S2 + previous-month air passenger mobility. | Tests air mobility without future leakage. |
| S6 | `spatial_R/spatial_inla_model_s6_rainfall_region.R` | S2 + rainfall-by-IBGE-mesoregion interaction. | Main rainfall heterogeneity model. |
| S7 | `spatial_R/spatial_inla_model_s7_temperature_region.R` | S2 + temperature-by-IBGE-mesoregion interaction. | Climate sensitivity model for temperature. |


### Model Comparison

| Model | Test MAE | Test RMSE | Test WAPE | Accuracy % | Test R2 | DIC | WAIC | Current interpretation |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| M0 | 18.8 | 96.95 | 0.88 | 12.47% | 0.19 | N/A | N/A | True null baseline; very weak prediction by design. |
| M1 | 19.14 | 97.88 | 0.875 | 12.49% | 0.1864 | N/A | N/A | Same-week covariates barely improved over M0. |
| M2 | 19.07 | 97.18 | 0.8692 | 13.07% | 0.2 | N/A | N/A | Lagged weather model; useful comparison but not yet saved in outputs. |
| M3 | 18.95 | 96.78 | 0.869 | 13.03% | 0.199 | N/A | N/A | Interpolation sensitivity model; interpolation did not clearly help in later models. |
| M4 | 4.7 | 12.94 | 0.5 | 49.5% | 0.825 | N/A | N/A | Lagged cases substantially improved prediction compared with weather-only models. |
| M5 | 4.3 | 11.7 | 0.49 | 50.53% | 0.84 | N/A | N/A | Strongest non-spatial model; main non-spatial comparison point. |
| M6 | 4.5 | 12.5 | 0.52 | 48% | 0.818 | N/A | N/A | Interpolation did not improve the M5 structure in the latest run. |
| S1 | 4.1906 | 11.3094 | 0.4839 | 51.6059% | 0.8511 | 93995.47 | 94367.80 | BYM2 spatial structure improved over M5. |
| S2 | 4.0437 | 10.6163 | 0.4670 | 53.3025% | 0.8688 | 93818.78 | 93989.47 | Main spatial model. |
| S3 | 4.1939 | 11.3304 | 0.4843 | 51.5682% | 0.8506 | 93835.93 | 94015.49 | Distance-weighted spatial lag did not beat adjacency-based S2. |
| S4-road | 4.0625 | 10.6889 | 0.4692 | 53.0849% | 0.8670 | 93819.11 | 93979.10 | Road effect was small; interval crossed zero. |
| S5-air | 4.0314 | 10.5433 | 0.4656 | 53.4441% | 0.8706 | 93818.31 | 93980.68 | Air effect was small; not a major improvement over S2. |
| S6 | 4.0198 | 10.5809 | 0.4642 | 53.5787% | 0.8697 | 93722.40 | 93880.28 | Rainfall effect varies by region; important for the research question. |
| S7 | 4.1713 | 11.0831 | 0.4817 | 51.8289% | 0.8570 | 93729.90 | 93914.20 | Temperature effect varies in some regions, but this is a secondary sensitivity model. |

WAPE-based accuracy:

```text
accuracy_pct = 100 * (1 - sum(abs(actual - predicted)) / sum(actual))
```

### Main S2 Fixed Effects

S2 is the clean main spatial model

| Term | Posterior mean | Interpretation |
| --- | ---: | --- |
| `rainfall_lag_z` | 0.2085 | Positive rainfall association after a weather lag. |
| `humidity_lag_z` | -0.0213 | Very small negative association. |
| `temperature_lag_z` | 0.0833 | Positive but smaller than rainfall. |
| `idhm_z` | 0.5238 | Higher IDHM municipalities have higher expected reported cases in this model. |
| `log_cases_lag_z` | 0.9368 | Strong temporal persistence; strongest predictor. |
| `neighbor_log_cases_lag_z` | 0.2772 | Neighboring recent cases add spatial transmission pressure. |

### S6 Rainfall-by-Region Effects

S6 estimates one rainfall effect for each official IBGE mesoregion. 

| Region | Rainfall effect | Rainfall relative risk | 95% RR interval | Interpretation |
| --- | ---: | ---: | --- | --- |
| Baixadas | -0.103 | 0.902 | 0.785 to 1.035 | Near-neutral to slightly negative; interval overlaps 1. |
| Centro Fluminense | 0.258 | 1.294 | 1.143 to 1.463 | Clear positive rainfall association. |
| Metropolitana do Rio de Janeiro | 0.146 | 1.157 | 1.086 to 1.231 | Positive reference-region rainfall effect. |
| Noroeste Fluminense | 0.284 | 1.329 | 1.175 to 1.500 | Strongest rainfall association. |
| Norte Fluminense | 0.176 | 1.192 | 1.028 to 1.381 | Positive rainfall association. |
| Sul Fluminense | 0.154 | 1.166 | 1.042 to 1.304 | Positive rainfall association. |

### S7 Temperature-by-Region Effects

S7 repeats the region-interaction logic for temperature. 

| Region | Temperature effect | Temperature relative risk | 95% RR interval | Interpretation |
| --- | ---: | ---: | --- | --- |
| Baixadas | 0.005 | 1.005 | 0.897 to 1.126 | Essentially neutral. |
| Centro Fluminense | 0.213 | 1.237 | 1.103 to 1.387 | Positive temperature association. |
| Metropolitana do Rio de Janeiro | 0.010 | 1.010 | 0.963 to 1.058 | Essentially neutral reference effect. |
| Noroeste Fluminense | 0.273 | 1.314 | 1.176 to 1.469 | Strongest temperature association. |
| Norte Fluminense | 0.111 | 1.118 | 0.994 to 1.256 | Positive but interval nearly touches 1. |
| Sul Fluminense | 0.107 | 1.113 | 1.011 to 1.225 | Positive temperature association. |

S7 is not the central model; it is included as a sensitivity model.

### Spatial Diagnostic Results

S2 residual spatial autocorrelation was evaluated with Moran's I on municipality-level standardized residuals from the negative-binomial model

| Variable | Moran's I | Permutation p-value | Interpretation |
| --- | ---: | ---: | --- |
| `municipality_aggregated_pearson_residual` | 0.1092 | 0.010 | Primary diagnostic; suggests modest remaining positive residual spatial autocorrelation. |
| `municipality_mean_pearson_residual` | 0.1680 | 0.001 | Sensitivity check; also suggests remaining positive residual spatial autocorrelation. |
| `municipality_mean_deviance_residual` | 0.0254 | 0.500 | Sensitivity check; no significant residual spatial autocorrelation detected. |

The Pearson-residual diagnostics suggest that S2 improves the spatial structure but does not remove every spatial pattern in the residuals. The deviance-residual sensitivity check is not significant, so the evidence is not one-dimensional. 

## Methodology

### Why a Negative-Binomial Model?

Dengue cases are nonnegative counts. A Poisson model often assumes too little variance for epidemiological case data, where outbreaks create bursts and overdispersion. The negative-binomial likelihood is more flexible:

```text
cases_it ~ NegativeBinomial(mu_it, alpha)
```

where `mu_it` is the expected case count for municipality `i` at week `t`, and `alpha` controls overdispersion.

### Log Link and Coefficient Interpretation

The models use a log link:

```text
log(mu_it) = intercept + predictors
```

This keeps predicted case counts positive. It also makes coefficients interpretable as multiplicative changes in expected cases. If a coefficient is `0.20`, then:

```text
exp(0.20) = 1.22
```

### Hierarchical Structure

The non-spatial models include repeated observations for each municipality and week. The basic structure is:

```text
log(mu_it) = intercept
             + municipality effect
             + week-of-year effect
             + year effect
             + covariate effects
```

### Spatial Structure

```text
spatial risk = structured adjacency effect + unstructured municipality effect
```

S2 also includes an adjacency-based lagged neighboring case term:

```text
neighbor_cases_lag[i,t] =
  mean cases among adjacent municipalities at t - CASE_LAG_WEEKS

neighbor_log_cases_lag[i,t] =
  log1p(neighbor_cases_lag[i,t])
```

### Lags

```text
CASE_LAG_WEEKS = 4
WEATHER_LAG_WEEKS = 6
```

Lagged cases:

```text
log_cases_lag = log1p(cases exactly CASE_LAG_WEEKS earlier)
```

Lagged weather:

```text
rainfall_lag    = rainfall exactly WEATHER_LAG_WEEKS earlier
humidity_lag    = humidity exactly WEATHER_LAG_WEEKS earlier
temperature_lag = temperature exactly WEATHER_LAG_WEEKS earlier
```

### Interpolation

Rainfall and IDHM were treated as complete in the current dataset. Missing humidity and temperature values were handled only in the interpolation models.

Interpolation models use municipality-level linear temporal interpolation for internal gaps:

```text
x(t) = x(t0) + ((t - t0) / (t1 - t0)) * (x(t1) - x(t0))
```

where `t0` and `t1` are observed dates around the missing value. The scripts limit interpolation to reasonably short internal gaps.

### Leakage Prevention

The project is designed to avoid using future information in prediction.

Key checks:

- Own-case lags must come from dates strictly before the current row.
- Neighbor-case lags must come from dates strictly before the current row.
- Held-out 2023 test rows are dropped if their own-case or neighbor-case lag would come from inside the test period.
- Standardization parameters for train/test evaluation are fit on training data only.
- S5 air mobility uses previous-month mobility, not future annual totals.
- Interpolation models are designed to avoid using held-out test values while training.

## Outputs and Figures

| Output | Created by | Meaning |
| --- | --- | --- |
| `outputs/s2_unexplained_spatial_effects.csv` | `spatial_R/map_s2_unexplained_effects.R` | Municipality-level BYM2 residual spatial effects and relative risks. |
| `outputs/s2_unexplained_spatial_relative_risk_map.png` | `spatial_R/map_s2_unexplained_effects.R` | Map of unexplained residual spatial relative risk. |
| `outputs/s2_morans_i_diagnostic.csv` | `spatial_R/s2_morans_i_diagnostic.R` | Moran's I diagnostic using standardized negative-binomial residuals. |
| `outputs/s2_standardized_residuals_by_municipio.csv` | `spatial_R/s2_morans_i_diagnostic.R` | Municipality-level Pearson and deviance residual summaries used for Moran's I. |
| `outputs/s6_rainfall_region_effects.csv` | `spatial_R/map_s6_rainfall_region_effect.R` | Region-specific rainfall effects and relative risks. |
| `outputs/s6_rainfall_region_effect_map.png` | `spatial_R/map_s6_rainfall_region_effect.R` | Map of rainfall relative risk by IBGE mesoregion. |
| `outputs/s7_temperature_region_effects.csv` | `spatial_R/map_s7_temperature_region_effect.R` | Region-specific temperature effects and relative risks. |
| `outputs/s7_temperature_region_effect_map.png` | `spatial_R/map_s7_temperature_region_effect.R` | Map of temperature relative risk by IBGE mesoregion. |
| `outputs/run_logs/` | S1-S7 model runs | Saved console logs used to update the model-results table. |

## How to Run

From the project folder:

```bash
cd /Users/ethanmassey/VS_Code_Test/Hierarchical-Bayesian-Model
```

Run the main non-spatial model:

```bash
python base_model/base_model_lag_cases_weather.py
```

Run the main spatial model:

```bash
Rscript spatial_R/spatial_inla_model_s2.R
```

Run the rainfall-by-region model:

```bash
Rscript spatial_R/spatial_inla_model_s6_rainfall_region.R
```

Create the rainfall region map:

```bash
Rscript spatial_R/map_s6_rainfall_region_effect.R
```

Run the temperature-by-region sensitivity model:

```bash
Rscript spatial_R/spatial_inla_model_s7_temperature_region.R
```

Create the temperature region map:

```bash
Rscript spatial_R/map_s7_temperature_region_effect.R
```

Run the S2 residual spatial diagnostic:

```bash
Rscript spatial_R/s2_morans_i_diagnostic.R
```

## Configuration

Common settings are near the top of the scripts.

| Setting | Meaning |
| --- | --- |
| `DATA_START_YEAR`, `DATA_END_YEAR` | Full modeling period. |
| `TRAIN_START_YEAR`, `TRAIN_END_YEAR` | Training period. |
| `TEST_START_YEAR`, `TEST_END_YEAR` | Held-out test period. |
| `CASE_LAG_WEEKS` | Number of weeks used for own-case and neighbor-case lags. |
| `WEATHER_LAG_WEEKS` | Number of weeks used for lagged weather. |
| `INTERPOLATION_LIMIT_WEEKS` | Maximum gap length for interpolation models. |
| `RUN_TRAIN_TEST_EVALUATION` | Whether to run train/test evaluation. |
| `SAVE_OUTPUTS` | Whether to save generated metrics and predictions. |
| `INLA_NUM_THREADS` | Threading setting for R-INLA. |

## Recommended Paper Framing

A good current paper structure for publication is:

1. Start with **M0-M5** to show that lagged cases are essential and that weather-only covariates are not enough.
2. Use **M5** as the main non-spatial benchmark.
3. Use **S2** as the main spatial model because it combines BYM2 residual spatial risk with observed neighboring dengue pressure.
4. Use **S6** to answer the rainfall research question: rainfall effects vary across IBGE mesoregions.
5. Use **S7** as a sensitivity model showing whether temperature has similar spatial heterogeneity.
6. Treat **S3-S5** as sensitivity or extension models rather than the main result.

Suggested tables and figures:

- Model comparison table for M0-M6 and S1-S7, with WAPE accuracy and R2 where available.
- Posterior coefficient table for S2.
- Region-specific rainfall table from S6.
- Rainfall-by-region map from S6.
- Residual spatial relative-risk map from S2.
- Moran's I diagnostic table for residual spatial autocorrelation.
- Sensitivity table for temperature, road mobility, and air mobility.

## Reproducibility Notes

- Python scripts use a fixed random seed where configured.
- Exact MCMC results may vary slightly across runs and package versions.
- INLA results may vary slightly across R-INLA versions.
- Models may use different row sets if interpolation fills rows that non-interpolation models drop.
- When comparing predictive performance, confirm whether the held-out test rows are identical.
- Several scripts default to `SAVE_OUTPUTS = FALSE`, so not every run leaves a CSV behind.

