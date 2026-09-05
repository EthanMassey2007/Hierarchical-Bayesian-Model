# Hierarchical Bayesian Dengue Modeling in Rio de Janeiro

M5 is the strongest non-spatial benchmark. S2 is the main fixed-rainfall spatial model. S6 and S8 test spatial rainfall variation, while S10 and S11 extend the rainfall effect across time and space-time.

## Project Status

The models are organized into 3 layers:

1. **M0-M6:** non-spatial hierarchical Bayesian models fit in Python/PyMC.
2. **S1-S5:** spatial and mobility extensions fit in R-INLA.
3. **S6-S11:** spatial, temporal, and spatiotemporal climate-effect models fit in R-INLA.

Current interpretation:

- **M5** strongest non-spatial model because lagged cases add a major predictive signal.
- **S2** cleanest main spatial model because it adds BYM2 spatial structure and adjacency-based neighboring lagged cases.
- **S6** tests whether rainfall effects differ by region.
- **S8** tests whether rainfall effects differ by municipality.
- **S10** tests whether rainfall effects vary over time.
- **S11** tests an additive municipality-plus-time rainfall effect.
- **S7** useful sensitivity model that tests whether temperature also varies by region.
- **S9** exploratory municipality-level temperature random-slope map.
- **S3-S5** are spatial or mobility sensitivity models. They are useful, but they have not replaced S2.

## Data

Primary modeling file:

```text
data/complete_combined_datasets.csv
```

Detailed data provenance and citation information is documented in
[`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md), with BibTeX entries in
[`references.bib`](references.bib).

Reproducibility setup and check commands are documented in
[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md). A lightweight `Makefile` provides
syntax checks and common model/map targets.

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
|-- DATA_AVAILABILITY.md
|-- REPRODUCIBILITY.md
|-- references.bib
|-- requirements.txt
|-- Makefile
|-- models/
|   |-- pymc/
|   |   |-- base_model.py
|   |   |-- base_model_covariates.py
|   |   |-- base_model_lag_weather.py
|   |   |-- base_model_interpolation.py
|   |   |-- base_model_lag_cases.py
|   |   |-- base_model_lag_cases_weather.py
|   |   `-- base_model_lag_cases_weather_interpolation.py
|   `-- r_inla/
|       |-- baseline/
|       |   |-- base_model_r_m0.R
|       |   |-- base_model_r_m1_covariates.R
|       |   |-- base_model_r_m2_lag_weather.R
|       |   |-- base_model_r_m3_interpolation.R
|       |   |-- base_model_r_m4_lag_cases.R
|       |   `-- base_model_r_m5_lag_weather_cases.R
|       `-- spatial/
|           |-- spatial_inla_model_s1.R
|           |-- spatial_inla_model_s2.R
|           |-- spatial_inla_model_s3.R
|           |-- spatial_inla_model_s4_road.R
|           |-- spatial_inla_model_s5_air.R
|           |-- spatial_inla_model_s6_rainfall_region.R
|           |-- spatial_inla_model_s7_temperature_region.R
|           |-- spatial_inla_model_s8_rainfall_municipality.R
|           |-- spatial_inla_model_s9_temperature_municipality.R
|           |-- spatial_inla_model_s10_rainfall_time.R
|           |-- spatial_inla_model_s11_rainfall_spacetime.R
|           |-- lag_sensitivity_climate_s2.R
|           |-- map_s2_unexplained_effects.R
|           |-- map_s2_spacetime_weather_residual_animations.R
|           |-- map_s6_rainfall_region_effect.R
|           |-- map_s7_temperature_region_effect.R
|           |-- s2_morans_i_diagnostic.R
|           `-- install_packages.R
|-- scripts/
|   |-- analysis/
|   |   |-- run_all_models_collect_results.R
|   |   |-- correlation_matrix.py
|   |   `-- unexplained_effects_map.py
|   `-- figures/
|       |-- plot_weekly_dengue_cases.py
|       |-- plot_municipal_mean_annual_incidence.py
|       |-- plot_model_comparison_figures.py
|       |-- plot_observed_predicted_s6_s11.py
|       `-- plot_s11_rainfall_effect_heterogeneity.py
|-- outputs/
|   |-- s2_unexplained_spatial_effects.csv
|   |-- s2_unexplained_spatial_relative_risk_map.png
|   |-- s2_morans_i_diagnostic.csv
|   |-- s6_rainfall_region_effects.csv
|   |-- s6_rainfall_region_effect_map.png
|   |-- s7_temperature_region_effects.csv
|   |-- s7_temperature_region_effect_map.png
|   |-- s8_rainfall_municipality_effects.csv
|   |-- s8_rainfall_municipality_effect_map.png
|   |-- s9_temperature_municipality_effects.csv
|   |-- s9_temperature_municipality_effect_map.png
|   |-- s10_rainfall_time_effects.csv
|   |-- s10_rainfall_time_effect_plot.png
|   |-- s10_rainfall_time_model_criteria.csv
|   |-- s10_rainfall_time_train_test_metrics.csv
|   |-- s11_rainfall_spacetime_municipality_effects.csv
|   |-- s11_rainfall_spacetime_time_effects.csv
|   |-- s11_rainfall_spacetime_combined_effects.csv
|   |-- s11_rainfall_spacetime_municipality_effect_map.png
|   |-- s11_rainfall_spacetime_time_effect_plot.png
|   |-- s11_rainfall_spacetime_model_criteria.csv
|   |-- s11_rainfall_spacetime_train_test_metrics.csv
|   |-- correlation_matrix_*
|   |-- s2_spacetime_rainfall_lag_animation_2023.gif
|   |-- s2_spacetime_temperature_lag_animation_2023.gif
|   |-- s2_spacetime_residual_animation_2023.gif
|   |-- descriptive_figures/
|   |-- model_comparison_figures/
|   |-- observed_predicted_figures/
|   |-- rainfall_effect_heterogeneity_figures/
|   |-- misc/
|   `-- run_logs/
`-- data/
```

## Python and R Ecosystem

The project has two modeling layers:

```text
data/complete_combined_datasets.csv
        |
        |-- Python/PyMC models
        |       `-- M0-M6: non-spatial baseline and feature tests
        |
        `-- R-INLA models
                `-- S1-S11: spatial models, mobility tests, climate-effect maps, diagnostics
```

- **Python/PyMC:** used for early non-spatial model development and sensitivity checks.
- **R-INLA non-spatial baselines:** used for fast final baseline comparison against the spatial rainfall models.
- **R-INLA:** used for spatial models because BYM2 spatial models are much faster in INLA than PyMC MCMC.
- **Shared input:** both ecosystems read `data/complete_combined_datasets.csv`.
- **No dependency on Python outputs:** the R scripts rebuild their own model data from the combined dataset.
- **Main workflow:** M5 becomes the non-spatial benchmark, S2 becomes the fixed-rainfall spatial benchmark, then S6/S8/S10/S11 test whether rainfall effects vary across space and time.

### Shared Data Flow

Most scripts follow the same basic workflow:

```text
read combined dataset
clean names
create dates
create model features
drop unusable rows
split train/test
standardize covariates
fit negative-binomial model
evaluate predictions
```

- **M0:** no covariates.
- **M5:** lagged weather, IDHM, and lagged own cases.
- **S2:** M5-style features plus lagged neighboring cases.
- **S6:** S2 plus rainfall effects that vary by IBGE mesoregion.
- **S8:** S2 plus rainfall effects that vary by municipality.
- **S10:** S2 plus rainfall effects that vary over ordered weeks.
- **S11:** S2 plus additive municipality and time rainfall effects.
- **S7:** S2 plus temperature effects that vary by IBGE mesoregion.

### R-INLA Baseline Layer

The final fast baseline ladder can be run in R-INLA:

| File | Model | What it does |
| --- | --- | --- |
| `base_model_r_m0.R` | R-M0 | Baseline with municipality, week, and year random effects; no fixed covariates. |
| `base_model_r_m1_covariates.R` | R-M1 | Adds same-week rainfall, humidity, temperature, and IDHM. |
| `base_model_r_m2_lag_weather.R` | R-M2 | Adds 12-week lagged rainfall, humidity, temperature, and IDHM. |
| `base_model_r_m3_interpolation.R` | R-M3 | Adds leakage-free humidity and temperature interpolation to the same-week covariate model. |
| `base_model_r_m4_lag_cases.R` | R-M4 | Adds four-week log lagged cases to the same-week covariate model. |
| `base_model_r_m5_lag_weather_cases.R` | R-M5 | Main non-spatial benchmark with lagged weather and log lagged cases. |

These R baseline files save:

| Output | Meaning |
| --- | --- |
| `r_m*_model_criteria.csv` | DIC, WAIC, LPML, and mean log CPO. |
| `r_m*_fixed_effects.csv` | INLA fixed-effect summaries and uncertainty intervals. |
| `r_m*_train_test_metrics.csv` | MAE, RMSE, WAPE, accuracy percentage, and R2. |

All R-INLA model scripts can be run and summarized with:

```bash
make all-results
```

This writes one log per model in `outputs/run_logs/`, then combines model criteria and held-out prediction metrics into `outputs/all_model_results_table.csv`. To rebuild the combined result CSVs from existing model outputs without refitting, run:

```bash
make collect-results
```

### Python/PyMC Layer

The Python scripts live in `base_model/`. They were used to build and validate the original non-spatial comparison ladder. They are slower because they use MCMC, so the R-INLA versions are preferred for final repeated model comparison.

| File | Model | What it does |
| --- | --- | --- |
| `base_model/base_model.py` | M0 | True baseline; no covariates. |
| `base_model/base_model_covariates.py` | M1 | Adds same-week rainfall, humidity, temperature, and IDHM. |
| `base_model/base_model_lag_weather.py` | M2 | Uses lagged weather instead of same-week weather. |
| `base_model/base_model_interpolation.py` | M3 | Adds leakage-free interpolation for humidity and temperature. |
| `base_model/base_model_lag_cases.py` | M4 | Adds lagged log cases. |
| `base_model/base_model_lag_cases_weather.py` | M5 | Main non-spatial model. |
| `base_model/base_model_lag_cases_weather_interpolation.py` | M6 | M5 plus interpolation. |

Common Python script structure:

| Code section | What it does |
| --- | --- |
| Settings | Years, lags, sampling settings, save options. |
| Paths | Finds `data/` and `outputs/`. |
| Cleaning | Cleans columns, municipality names, and dates. |
| Feature building | Adds covariates, lags, or interpolation. |
| `build_model_dataframe()` | Creates the final modeling dataframe. |
| `prepare_arrays()` | Converts data into PyMC-ready arrays. |
| `fit_model()` | Fits the negative-binomial Bayesian model. |
| `posterior_expected_cases()` | Creates expected case predictions. |
| `run_train_test_evaluation()` | Fits on train data and evaluates test data. |
| `main()` | Runs the script. |

Basic Python model shape:

```text
cases_it ~ NegativeBinomial(mu_it, alpha)

log(mu_it) =
  intercept
  + municipality random effect
  + week-of-year random effect
  + year random effect
  + model-specific covariates
```

- `municipality random effect`: stable differences between municipalities.
- `week-of-year random effect`: seasonality.
- `year random effect`: year-to-year shifts.
- `model-specific covariates`: rainfall, temperature, IDHM, lags, or interpolation depending on the model.

### R-INLA Spatial Layer

The R scripts live in `spatial_R/`. They add spatial structure on top of the M5-style model.

| File | Model | What it does |
| --- | --- | --- |
| `spatial_R/spatial_inla_model_s1.R` | S1 | Adds BYM2 spatial random effects. |
| `spatial_R/spatial_inla_model_s2.R` | S2 | Main spatial model; adds adjacency-based lagged neighboring cases. |
| `spatial_R/spatial_inla_model_s3.R` | S3 | Tests distance-weighted neighboring cases. |
| `spatial_R/spatial_inla_model_s4_road.R` | S4 | Tests road/fluvial mobility. |
| `spatial_R/spatial_inla_model_s5_air.R` | S5 | Tests lagged air passenger mobility. |
| `spatial_R/spatial_inla_model_s6_rainfall_region.R` | S6 | Tests rainfall effects by region. |
| `spatial_R/spatial_inla_model_s7_temperature_region.R` | S7 | Tests temperature effects by region. |
| `spatial_R/spatial_inla_model_s8_rainfall_municipality.R` | S8 | Exploratory municipality-level rainfall random slopes. |
| `spatial_R/spatial_inla_model_s9_temperature_municipality.R` | S9 | Exploratory municipality-level temperature random slopes. |
| `spatial_R/spatial_inla_model_s10_rainfall_time.R` | S10 | Tests time-varying rainfall random slopes. |
| `spatial_R/spatial_inla_model_s11_rainfall_spacetime.R` | S11 | Tests additive municipality-plus-time rainfall random slopes. |

Main S2 model terms:

```text
cases_it ~ NegativeBinomial(mu_it)

log(mu_it) =
  intercept
  + rainfall_lag_z
  + humidity_lag_z
  + temperature_lag_z
  + idhm_z
  + log_cases_lag_z
  + neighbor_log_cases_lag_z
  + week random effect
  + year random effect
  + BYM2 spatial random effect
```

- **BYM2 structured effect:** uses `data/adjacency_matrix_correct.parquet`.
- **BYM2 unstructured effect:** captures municipality-specific noise.
- **Neighbor lag:** uses adjacent municipalities' dengue cases from the lag period.
- **S6/S7 region effects:** use IBGE mesoregion information from `data/hub_pop_density.csv`.
- **S8/S9 municipality effects:** use partially pooled municipality-level random slopes.
- **S10 rainfall time effect:** uses an RW1 time-varying rainfall slope across ordered weeks.
- **S11 rainfall space-time effect:** combines a municipality rainfall random slope and an RW1 time-varying rainfall slope.

### Inferential Rainfall Model Ladder

The rainfall-focused model ladder is:

| Model | Rainfall specification | Main question |
| --- | --- | --- |
| M5 | Fixed rainfall effect, no spatial INLA structure | Non-spatial benchmark. |
| S2 | Fixed rainfall effect with BYM2 and neighboring lagged cases | Spatial benchmark with conventional rainfall effect. |
| S6 | Rainfall effect varies by IBGE mesoregion | Is rainfall sensitivity region-specific? |
| S8 | Rainfall effect varies by municipality | Is there finer local rainfall heterogeneity? |
| S10 | Rainfall effect varies over ordered weeks | Does rainfall sensitivity change over time? |
| S11 | Rainfall effect varies additively by municipality and time | Does evidence support a spatiotemporally varying rainfall effect? |

Primary inferential outputs:

| Output | Meaning |
| --- | --- |
| rainfall coefficient / relative risk | Direction and size of the rainfall association. |
| posterior interval | Uncertainty in the rainfall association. |
| DIC / WAIC / CPO | Evidence and fit comparison across rainfall specifications. |
| rainfall maps / time plots | Where and when rainfall sensitivity is stronger or weaker. |
| MAE / RMSE / WAPE | Secondary held-out predictive validation. |

### R Diagnostics and Map Scripts

These scripts turn fitted spatial models into maps and diagnostics.

| File | Depends on | Output | Purpose |
| --- | --- | --- | --- |
| `spatial_R/map_s2_unexplained_effects.R` | S2 | `outputs/s2_unexplained_spatial_effects.csv`, `outputs/s2_unexplained_spatial_relative_risk_map.png` | Maps unexplained spatial relative risk. |
| `spatial_R/s2_morans_i_diagnostic.R` | S2 | `outputs/s2_morans_i_diagnostic.csv`, `outputs/s2_standardized_residuals_by_municipio.csv` | Runs Moran's I on standardized negative-binomial residuals. |
| `spatial_R/map_s6_rainfall_region_effect.R` | S6 | `outputs/s6_rainfall_region_effects.csv`, `outputs/s6_rainfall_region_effect_map.png` | Maps rainfall relative risk by region. |
| `spatial_R/map_s7_temperature_region_effect.R` | S7 | `outputs/s7_temperature_region_effects.csv`, `outputs/s7_temperature_region_effect_map.png` | Maps temperature relative risk by region. |
| `spatial_R/map_s2_spacetime_weather_residual_animations.R` | S2 | `outputs/s2_spacetime_*_animation_2023.gif`, `outputs/s2_spacetime_weather_residuals_by_municipio_week.csv` | Animates 2023 lagged rainfall, lagged temperature, and S2 residuals. |

Some R helper scripts source S2 with:

```text
INLA_RUN_MODEL=0
```

That lets them reuse S2 functions without immediately running the full S2 model when the file is imported.

### How the Files Work Together

Model path:

```text
M0 -> M1 -> M2 -> M3 -> M4 -> M5 -> M6
                         |
                         `-- M5 becomes the non-spatial benchmark

M5 logic -> S1 -> S2 -> S3/S4/S5
                 |
                 |-- S6: rainfall varies by region
                 |-- S7: temperature varies by region
                 |-- S8: rainfall varies by municipality
                 |-- S9: temperature varies by municipality
                 |-- S10: rainfall varies over time
                 |-- S11: rainfall varies by municipality and time
                 |-- S2 weekly rainfall/temperature/residual GIFs
                 |-- S2 unexplained-effect map
                 `-- S2 Moran's I residual diagnostic
```

How to think about the models:

- **M0-M6:** build the baseline story.
- **M5:** best non-spatial benchmark.
- **S1:** asks whether spatial random effects help.
- **S2:** main spatial model.
- **S3:** checks whether distance weighting beats adjacency.
- **S4-S5:** check mobility.
- **S6:** main rainfall heterogeneity model.
- **S7:** temperature sensitivity model.
- **S8-S9:** exploratory municipality-level climate-effect maps.
- **S10-S11:** temporal and additive spatiotemporal rainfall-effect models for the central rainfall hypothesis.
- **Maps/diagnostics:** explain spatial effects after modeling.

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
| S8 | `spatial_R/spatial_inla_model_s8_rainfall_municipality.R` | S2 + municipality-level rainfall random slope. | Exploratory fine-scale rainfall heterogeneity map. |
| S9 | `spatial_R/spatial_inla_model_s9_temperature_municipality.R` | S2 + municipality-level temperature random slope. | Exploratory fine-scale temperature heterogeneity map. |
| S10 | `spatial_R/spatial_inla_model_s10_rainfall_time.R` | S2 + RW1 time-varying rainfall random slope. | Tests whether rainfall sensitivity changes across ordered weeks. |
| S11 | `spatial_R/spatial_inla_model_s11_rainfall_spacetime.R` | S2 + municipality rainfall random slope + RW1 time-varying rainfall random slope. | Tests additive spatiotemporal rainfall heterogeneity. |

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
| S10 | 4.1318 | 10.7383 | 0.4771 | 52.2855% | 0.8658 | 93057.94 | 93434.95 | Time-varying rainfall effect improves information criteria but not held-out prediction over S2/S6. |
| S11 | 4.6196 | 12.3057 | 0.5335 | 46.6514% | 0.8237 | 92566.65 | 92906.11 | Best information criteria among rainfall-effect models, but weaker held-out prediction; interpret as inferential rather than predictive evidence. |

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
Region-specific intervals use an approximate normal interval for the main climate effect plus the region interaction. This avoids adding marginal quantiles directly.

| Region | Rainfall effect | Rainfall relative risk | 95% RR interval | Interpretation |
| --- | ---: | ---: | --- | --- |
| Baixadas | -0.103 | 0.902 | 0.818 to 0.996 | Slight negative rainfall association under the approximate interval. |
| Centro Fluminense | 0.258 | 1.294 | 1.186 to 1.412 | Clear positive rainfall association. |
| Metropolitana do Rio de Janeiro | 0.146 | 1.157 | 1.087 to 1.232 | Positive reference-region rainfall effect. |
| Noroeste Fluminense | 0.284 | 1.329 | 1.219 to 1.448 | Strongest rainfall association. |
| Norte Fluminense | 0.176 | 1.192 | 1.073 to 1.325 | Positive rainfall association. |
| Sul Fluminense | 0.154 | 1.166 | 1.077 to 1.264 | Positive rainfall association. |

### S7 Temperature-by-Region Effects

S7 repeats the region-interaction logic for temperature.
Region-specific intervals use the same approximate normal interval logic as S6.

| Region | Temperature effect | Temperature relative risk | 95% RR interval | Interpretation |
| --- | ---: | ---: | --- | --- |
| Baixadas | 0.005 | 1.005 | 0.926 to 1.091 | Essentially neutral. |
| Centro Fluminense | 0.213 | 1.237 | 1.140 to 1.343 | Positive temperature association. |
| Metropolitana do Rio de Janeiro | 0.010 | 1.010 | 0.963 to 1.058 | Essentially neutral reference effect. |
| Noroeste Fluminense | 0.273 | 1.314 | 1.214 to 1.423 | Strongest temperature association. |
| Norte Fluminense | 0.111 | 1.118 | 1.028 to 1.216 | Positive temperature association. |
| Sul Fluminense | 0.107 | 1.113 | 1.040 to 1.191 | Positive temperature association. |

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
WEATHER_LAG_WEEKS = 12
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
| `outputs/s8_rainfall_municipality_effects.csv` | `spatial_R/spatial_inla_model_s8_rainfall_municipality.R` | Municipality-specific rainfall random-slope summaries. |
| `outputs/s8_rainfall_municipality_effect_map.png` | `spatial_R/spatial_inla_model_s8_rainfall_municipality.R` | Exploratory municipality-level rainfall relative-risk map. |
| `outputs/s9_temperature_municipality_effects.csv` | `spatial_R/spatial_inla_model_s9_temperature_municipality.R` | Municipality-specific temperature random-slope summaries. |
| `outputs/s9_temperature_municipality_effect_map.png` | `spatial_R/spatial_inla_model_s9_temperature_municipality.R` | Exploratory municipality-level temperature relative-risk map. |
| `outputs/s2_spacetime_weather_residuals_by_municipio_week.csv` | `spatial_R/map_s2_spacetime_weather_residual_animations.R` | Weekly S2 fitted means and standardized Pearson residuals. |
| `outputs/s2_spacetime_rainfall_lag_animation_2023.gif` | `spatial_R/map_s2_spacetime_weather_residual_animations.R` | Animated 2023 map of lagged rainfall used by S2. |
| `outputs/s2_spacetime_temperature_lag_animation_2023.gif` | `spatial_R/map_s2_spacetime_weather_residual_animations.R` | Animated 2023 map of lagged temperature used by S2. |
| `outputs/s2_spacetime_residual_animation_2023.gif` | `spatial_R/map_s2_spacetime_weather_residual_animations.R` | Animated 2023 map of full-data S2 standardized Pearson residuals. |
| `outputs/s10_rainfall_time_effects.csv` | `spatial_R/spatial_inla_model_s10_rainfall_time.R` | Time-varying rainfall effect summaries. |
| `outputs/s10_rainfall_time_effect_plot.png` | `spatial_R/spatial_inla_model_s10_rainfall_time.R` | Plot of the RW1 rainfall time effect. |
| `outputs/s10_rainfall_time_model_criteria.csv` | `spatial_R/spatial_inla_model_s10_rainfall_time.R` | DIC and WAIC for S10. |
| `outputs/s10_rainfall_time_train_test_metrics.csv` | `spatial_R/spatial_inla_model_s10_rainfall_time.R` | Train/test metrics for S10. |
| `outputs/s11_rainfall_spacetime_municipality_effects.csv` | `spatial_R/spatial_inla_model_s11_rainfall_spacetime.R` | Municipality-level rainfall random-slope summaries from S11. |
| `outputs/s11_rainfall_spacetime_time_effects.csv` | `spatial_R/spatial_inla_model_s11_rainfall_spacetime.R` | Time-varying rainfall random-slope summaries from S11. |
| `outputs/s11_rainfall_spacetime_combined_effects.csv` | `spatial_R/spatial_inla_model_s11_rainfall_spacetime.R` | Additive municipality-by-time rainfall-effect summaries. |
| `outputs/s11_rainfall_spacetime_municipality_effect_map.png` | `spatial_R/spatial_inla_model_s11_rainfall_spacetime.R` | Map of S11 municipality rainfall effects. |
| `outputs/s11_rainfall_spacetime_time_effect_plot.png` | `spatial_R/spatial_inla_model_s11_rainfall_spacetime.R` | Plot of the S11 rainfall time effect. |
| `outputs/s11_rainfall_spacetime_model_criteria.csv` | `spatial_R/spatial_inla_model_s11_rainfall_spacetime.R` | DIC and WAIC for S11. |
| `outputs/s11_rainfall_spacetime_train_test_metrics.csv` | `spatial_R/spatial_inla_model_s11_rainfall_spacetime.R` | Train/test metrics for S11. |
| `outputs/all_model_criteria.csv` | `run_all_models_collect_results.R` | Combined DIC, WAIC, LPML, and related criteria from all R-INLA models. |
| `outputs/all_model_train_test_metrics.csv` | `run_all_models_collect_results.R` | Combined train/test MAE, RMSE, WAPE, accuracy percentage, and R2 from all R-INLA models. |
| `outputs/all_model_results_table.csv` | `run_all_models_collect_results.R` | Compact model-comparison table combining WAIC/DIC with held-out test metrics and ranks. |
| `outputs/all_model_run_status.csv` | `run_all_models_collect_results.R` | Status file showing which model scripts completed or failed during a full run. |
| `outputs/correlation_matrix_*` | `correlation_matrix.py` and prior correlation-matrix runs | Covariate correlation matrices, heatmaps, and metadata. |
| `outputs/r_m*_model_criteria.csv`, `outputs/r_m*_fixed_effects.csv`, `outputs/r_m*_train_test_metrics.csv` | `base_model_r_m*.R` | Fast R-INLA baseline summaries for R-M0 through R-M5. |
| `outputs/run_logs/` | R-INLA model runs | Saved console logs used to check full model-run status. |

Map interpretation note:

- Static S2/S6/S7/S8/S9/S11 maps are full-data explanatory maps, not held-out prediction maps.
- The animated residual GIF is a full-data residual diagnostic. It shows weeks and places where S2 underpredicted or overpredicted after accounting for the fitted model.
- Use train/test metrics for predictive claims. Use maps for interpretation, diagnostics, and spatial pattern description.

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

Create the weekly rainfall, temperature, and unexplained-residual GIFs:

```bash
Rscript spatial_R/map_s2_spacetime_weather_residual_animations.R
```

Run the exploratory municipality-level climate-effect maps:

```bash
Rscript spatial_R/spatial_inla_model_s8_rainfall_municipality.R
Rscript spatial_R/spatial_inla_model_s9_temperature_municipality.R
```

Run the time-varying and spatiotemporal rainfall-effect models:

```bash
Rscript spatial_R/spatial_inla_model_s10_rainfall_time.R
Rscript spatial_R/spatial_inla_model_s11_rainfall_spacetime.R
```

Run the fast R-INLA non-spatial baseline ladder:

```bash
Rscript base_model_r_m0.R
Rscript base_model_r_m1_covariates.R
Rscript base_model_r_m4_lag_cases.R
Rscript base_model_r_m5_lag_weather_cases.R
```

Create the final-model covariate correlation matrix:

```bash
python correlation_matrix.py
```

Run lightweight syntax checks:

```bash
make check
```

Convenience Make targets are available for the most common scripts:

```bash
make m0
make m5
make s1
make s2
make s6
make s7
make diagnostics
make maps
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

R package note:

- The GIF script requires `gganimate`, `gifski`, and `scales` in addition to the core R-INLA mapping packages.
- If those packages are missing, install them with:

```bash
Rscript -e 'install.packages(c("gganimate", "gifski", "scales"), repos="https://cloud.r-project.org")'
```

## Recommended Paper Framing

A good current paper structure for publication is:

1. Start with **M0-M5** to show that lagged cases are essential and that weather-only covariates are not enough.
2. Use **M5** as the main non-spatial benchmark.
3. Use **S2** as the main spatial model because it combines BYM2 residual spatial risk with observed neighboring dengue pressure.
4. Use **S6** to answer the spatial rainfall question: rainfall effects vary across IBGE mesoregions.
5. Use **S10** and **S11** to test whether rainfall sensitivity changes over time and across space-time.
6. Use **S7** as a sensitivity model showing whether temperature has similar spatial heterogeneity.
7. Treat **S3-S5** as sensitivity or extension models rather than the main result.

Suggested tables and figures:

- Model comparison table for M0-M6 and S1-S11, with WAPE accuracy and R2 where available.
- Posterior coefficient table for S2.
- Region-specific rainfall table from S6.
- Rainfall-by-region map from S6.
- Time-varying rainfall-effect plot from S10.
- Municipality and time rainfall-effect summaries from S11.
- Residual spatial relative-risk map from S2.
- Moran's I diagnostic table for residual spatial autocorrelation.
- Sensitivity table for temperature, road mobility, and air mobility.
- Optional supplementary GIFs for 2023 lagged rainfall, lagged temperature, and S2 weekly residuals.
- Optional supplementary S8/S9 municipality-level effect maps, clearly labeled as exploratory.

## Reproducibility Notes

- Python scripts use a fixed random seed where configured.
- Exact MCMC results may vary slightly across runs and package versions.
- INLA results may vary slightly across R-INLA versions.
- Models may use different row sets if interpolation fills rows that non-interpolation models drop.
- When comparing predictive performance, confirm whether the held-out test rows are identical.
- Several scripts default to `SAVE_OUTPUTS = FALSE`, so not every run leaves a CSV behind.
- Before journal submission or public archival, confirm the redistribution status of all files under `data/`; `DATA_AVAILABILITY.md` currently flags `data/RJ.json` as still needing confirmed provenance.
- Use `references.bib` for source citations and add a repository-level software citation or DOI if the code is archived on Zenodo, OSF, or another preservation service.
