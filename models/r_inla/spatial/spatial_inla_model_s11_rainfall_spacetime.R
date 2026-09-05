# =========================================================
# S11: S2 + spatiotemporally varying rainfall effect using R-INLA
# =========================================================
# Uses the S2 model/data pipeline and adds additive municipality-level
# and time-varying random slopes for lagged rainfall. This estimates
# whether rainfall sensitivity varies across both space and time.
#
# Run from the project root:
#   Rscript models/r_inla/spatial/spatial_inla_model_s11_rainfall_spacetime.R
#
# Outputs:
#   outputs/s11_rainfall_spacetime_municipality_effects.csv
#   outputs/s11_rainfall_spacetime_time_effects.csv
#   outputs/s11_rainfall_spacetime_combined_effects.csv
#   outputs/s11_rainfall_spacetime_municipality_effect_map.png
#   outputs/s11_rainfall_spacetime_time_effect_plot.png

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
  library(sf)
  library(ggplot2)
})


# =========================================================
# Paths
# =========================================================
script_arg <- commandArgs(trailingOnly = FALSE)
script_file_arg <- script_arg[grepl("^--file=", script_arg)]
if (length(script_file_arg) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  SCRIPT_DIR <- getwd()
}

PROJECT_DIR_OVERRIDE <- Sys.getenv("HBM_PROJECT_DIR", "")
if (nzchar(PROJECT_DIR_OVERRIDE)) {
  PROJECT_DIR <- normalizePath(PROJECT_DIR_OVERRIDE)
} else if (dir.exists(file.path(getwd(), "data")) && dir.exists(file.path(getwd(), "models", "r_inla", "spatial"))) {
  PROJECT_DIR <- normalizePath(getwd())
} else if (dir.exists(file.path(SCRIPT_DIR, "data"))) {
  PROJECT_DIR <- normalizePath(SCRIPT_DIR)
} else {
  PROJECT_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
}

DATA_DIR_PROJECT <- file.path(PROJECT_DIR, "data")
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
S2_SCRIPT <- file.path(PROJECT_DIR, "models", "r_inla", "spatial", "spatial_inla_model_s2.R")
RJ_GEOJSON <- file.path(DATA_DIR_PROJECT, "RJ.json")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

EFFECTS_CSV <- file.path(OUTPUT_DIR, "s11_rainfall_spacetime_municipality_effects.csv")
TIME_EFFECTS_CSV <- file.path(OUTPUT_DIR, "s11_rainfall_spacetime_time_effects.csv")
COMBINED_EFFECTS_CSV <- file.path(OUTPUT_DIR, "s11_rainfall_spacetime_combined_effects.csv")
FIXED_EFFECTS_CSV <- file.path(OUTPUT_DIR, "s11_rainfall_spacetime_fixed_effects.csv")
RAINFALL_AVERAGE_EFFECT_CSV <- file.path(OUTPUT_DIR, "s11_rainfall_spacetime_average_rainfall_effect.csv")
CRITERIA_CSV <- file.path(OUTPUT_DIR, "s11_rainfall_spacetime_model_criteria.csv")
METRICS_CSV <- file.path(OUTPUT_DIR, "s11_rainfall_spacetime_train_test_metrics.csv")
FULL_PREDICTIONS_CSV <- file.path(OUTPUT_DIR, "s11_rainfall_spacetime_full_predictions.csv")
TEST_PREDICTIONS_CSV <- file.path(OUTPUT_DIR, "s11_rainfall_spacetime_test_predictions.csv")
MAP_PNG <- file.path(OUTPUT_DIR, "s11_rainfall_spacetime_municipality_effect_map.png")
TIME_PNG <- file.path(OUTPUT_DIR, "s11_rainfall_spacetime_time_effect_plot.png")


# =========================================================
# Load S2 functions without running S2 main()
# =========================================================
Sys.setenv(INLA_RUN_MODEL = "0")
source(S2_SCRIPT)

# Ensure sourced S2 globals resolve to project-root data/.
BASE_DIR <- PROJECT_DIR
DATA_DIR <- DATA_DIR_PROJECT
OUTPUT_DIR <- file.path(BASE_DIR, "outputs")
COMBINED_FILE <- file.path(DATA_DIR, "complete_combined_datasets.csv")
MUNICIPIOS_FILE <- file.path(DATA_DIR, "municipios.csv")
HUB_FILE <- file.path(DATA_DIR, "hub_pop_density.csv")
ADJACENCY_FILE <- file.path(DATA_DIR, "adjacency_matrix_correct.parquet")
GRAPH_FILE <- file.path(tempdir(), "rj_municipality_inla_s11.graph")


# =========================================================
# S11 helpers
# =========================================================
add_spacetime_slope_indices <- function(df) {
  slope_lookup <- unique(df[, .(municipio, ibge_code)])
  setorder(slope_lookup, municipio)
  slope_lookup[, rainfall_municipality_slope_idx := .I]

  df <- merge(df, slope_lookup, by = c("municipio", "ibge_code"), all.x = TRUE)

  time_lookup <- unique(df[, .(date, year, week)])
  setorder(time_lookup, date)
  time_lookup[, rainfall_time_idx := .I]

  df <- merge(df, time_lookup, by = c("date", "year", "week"), all.x = TRUE)
  setorder(df, municipio, date)

  list(df = df, municipality_lookup = slope_lookup, time_lookup = time_lookup)
}

build_s11_formula <- function(graph_file) {
  z_covariates <- paste0(BASE_COVARIATES, "_z")
  formula <- as.formula(paste(
    "cases ~ 1 +",
    paste(z_covariates, collapse = " + "),
    "+ f(rainfall_municipality_slope_idx, rainfall_lag_z, model = 'iid')",
    "+ f(rainfall_time_idx, rainfall_lag_z, model = 'rw1', scale.model = TRUE, constr = TRUE)",
    "+ f(week_idx, model = 'iid')",
    "+ f(year_idx, model = 'iid')",
    "+ f(spatial_idx, model = 'bym2', graph = graph_file, scale.model = TRUE)",
    "+ f(isolated_idx, model = 'iid')"
  ))
  environment(formula) <- environment()
  formula
}

fit_s11_inla <- function(model_dt, graph_file) {
  formula <- build_s11_formula(graph_file)
  fit <- inla(
    formula,
    family = "nbinomial",
    data = model_dt,
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    num.threads = INLA_NUM_THREADS,
    verbose = FALSE
  )
  fit$.model_data <- model_dt
  fit$.formula <- formula
  fit
}

predict_s11_mean <- function(fit, new_dt) {
  pred_dt <- copy(new_dt)
  pred_dt[, cases := NA_integer_]

  fit_dt <- rbindlist(
    list(fit$.model_data, pred_dt),
    use.names = TRUE,
    fill = TRUE
  )

  result <- inla(
    fit$.formula,
    family = "nbinomial",
    data = fit_dt,
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
    num.threads = INLA_NUM_THREADS,
    verbose = FALSE
  )

  pred_rows <- (nrow(fit$.model_data) + 1L):nrow(fit_dt)
  result$summary.fitted.values$mean[pred_rows]
}

extract_rainfall_fixed <- function(fit) {
  fixed <- as.data.table(fit$summary.fixed, keep.rownames = "term")
  rainfall <- fixed[term == "rainfall_lag_z"]
  if (nrow(rainfall) != 1) {
    stop("Could not find rainfall_lag_z in fixed effects.")
  }
  rainfall
}

save_fixed_effect_outputs <- function(fit, model_id, fixed_csv, rainfall_csv) {
  fixed <- as.data.table(fit$summary.fixed, keep.rownames = "term")
  fixed[, model := model_id]
  setcolorder(fixed, c("model", "term", setdiff(names(fixed), c("model", "term"))))
  fwrite(fixed, fixed_csv)

  rainfall <- copy(fixed[term == "rainfall_lag_z"])
  if (nrow(rainfall) != 1) {
    warning("Could not find rainfall_lag_z in fixed effects; rainfall average-effect CSV was not written.")
    return(invisible(NULL))
  }

  rainfall[, `:=`(
    relative_risk_mean = exp(mean),
    relative_risk_q025 = exp(`0.025quant`),
    relative_risk_q975 = exp(`0.975quant`),
    interpretation = fifelse(
      `0.025quant` > 0,
      "positive",
      fifelse(`0.975quant` < 0, "negative", "not clearly different from null")
    )
  )]
  fwrite(rainfall, rainfall_csv)
  invisible(rainfall)
}

summarize_municipality_rainfall_effects <- function(fit, municipality_lookup) {
  rainfall <- extract_rainfall_fixed(fit)

  random <- as.data.table(fit$summary.random$rainfall_municipality_slope_idx)
  random[, rainfall_municipality_slope_idx := as.integer(ID)]
  effects <- merge(municipality_lookup, random, by = "rainfall_municipality_slope_idx", all.x = TRUE)

  effects[, rainfall_effect_mean := rainfall$mean + mean]
  effects[, rainfall_effect_sd := sqrt(rainfall$sd^2 + sd^2)]
  effects[, rainfall_effect_q025 := rainfall_effect_mean - 1.96 * rainfall_effect_sd]
  effects[, rainfall_effect_q975 := rainfall_effect_mean + 1.96 * rainfall_effect_sd]
  effects[, rainfall_relative_risk := exp(rainfall_effect_mean)]
  effects[, rainfall_relative_risk_q025 := exp(rainfall_effect_q025)]
  effects[, rainfall_relative_risk_q975 := exp(rainfall_effect_q975)]

  effects <- effects[, .(
    municipio,
    ibge_code,
    rainfall_municipality_slope_idx,
    rainfall_effect_mean,
    rainfall_effect_sd,
    rainfall_effect_q025,
    rainfall_effect_q975,
    rainfall_relative_risk,
    rainfall_relative_risk_q025,
    rainfall_relative_risk_q975
  )]
  setorder(effects, municipio)
  effects
}

summarize_time_rainfall_effects <- function(fit, time_lookup) {
  rainfall <- extract_rainfall_fixed(fit)

  random <- as.data.table(fit$summary.random$rainfall_time_idx)
  random[, rainfall_time_idx := as.integer(ID)]
  effects <- merge(time_lookup, random, by = "rainfall_time_idx", all.x = TRUE)

  effects[, rainfall_effect_mean := rainfall$mean + mean]
  effects[, rainfall_effect_sd := sqrt(rainfall$sd^2 + sd^2)]
  effects[, rainfall_effect_q025 := rainfall_effect_mean - 1.96 * rainfall_effect_sd]
  effects[, rainfall_effect_q975 := rainfall_effect_mean + 1.96 * rainfall_effect_sd]
  effects[, rainfall_relative_risk := exp(rainfall_effect_mean)]
  effects[, rainfall_relative_risk_q025 := exp(rainfall_effect_q025)]
  effects[, rainfall_relative_risk_q975 := exp(rainfall_effect_q975)]

  effects <- effects[, .(
    date,
    year,
    week,
    rainfall_time_idx,
    rainfall_effect_mean,
    rainfall_effect_sd,
    rainfall_effect_q025,
    rainfall_effect_q975,
    rainfall_relative_risk,
    rainfall_relative_risk_q025,
    rainfall_relative_risk_q975
  )]
  setorder(effects, date)
  effects
}

summarize_combined_spacetime_effects <- function(fit, municipality_lookup, time_lookup) {
  rainfall <- extract_rainfall_fixed(fit)

  municipality_random <- as.data.table(fit$summary.random$rainfall_municipality_slope_idx)
  municipality_random[, rainfall_municipality_slope_idx := as.integer(ID)]
  municipality_random <- merge(
    municipality_lookup,
    municipality_random[, .(rainfall_municipality_slope_idx, municipality_mean = mean, municipality_sd = sd)],
    by = "rainfall_municipality_slope_idx",
    all.x = TRUE
  )

  time_random <- as.data.table(fit$summary.random$rainfall_time_idx)
  time_random[, rainfall_time_idx := as.integer(ID)]
  time_random <- merge(
    time_lookup,
    time_random[, .(rainfall_time_idx, time_mean = mean, time_sd = sd)],
    by = "rainfall_time_idx",
    all.x = TRUE
  )

  combined <- CJ(
    rainfall_municipality_slope_idx = municipality_random$rainfall_municipality_slope_idx,
    rainfall_time_idx = time_random$rainfall_time_idx
  )
  combined <- merge(combined, municipality_random, by = "rainfall_municipality_slope_idx", all.x = TRUE)
  combined <- merge(combined, time_random, by = "rainfall_time_idx", all.x = TRUE)

  combined[, rainfall_effect_mean := rainfall$mean + municipality_mean + time_mean]
  combined[, rainfall_effect_sd := sqrt(rainfall$sd^2 + municipality_sd^2 + time_sd^2)]
  combined[, rainfall_effect_q025 := rainfall_effect_mean - 1.96 * rainfall_effect_sd]
  combined[, rainfall_effect_q975 := rainfall_effect_mean + 1.96 * rainfall_effect_sd]
  combined[, rainfall_relative_risk := exp(rainfall_effect_mean)]
  combined[, rainfall_relative_risk_q025 := exp(rainfall_effect_q025)]
  combined[, rainfall_relative_risk_q975 := exp(rainfall_effect_q975)]

  combined <- combined[, .(
    municipio,
    ibge_code,
    date,
    year,
    week,
    rainfall_effect_mean,
    rainfall_effect_sd,
    rainfall_effect_q025,
    rainfall_effect_q975,
    rainfall_relative_risk,
    rainfall_relative_risk_q025,
    rainfall_relative_risk_q975
  )]
  setorder(combined, municipio, date)
  combined
}

plot_rainfall_municipality_map <- function(effects) {
  rj <- st_read(RJ_GEOJSON, quiet = TRUE)
  rj$ibge_code <- as.integer(rj$GEOCODIGO)

  map_sf <- merge(rj, effects, by = "ibge_code", all.x = TRUE, sort = FALSE)

  missing_geo_effects <- map_sf[is.na(map_sf$rainfall_relative_risk), ]
  if (nrow(missing_geo_effects) > 0) {
    cat("Warning: GeoJSON contains municipalities outside the S11 model set:", nrow(missing_geo_effects), "\n")
  }

  p <- ggplot(map_sf) +
    geom_sf(aes(fill = rainfall_relative_risk), color = "grey70", linewidth = 0.12) +
    scale_fill_gradient2(
      low = "#2c7bb6",
      mid = "white",
      high = "#d7191c",
      midpoint = 1,
      na.value = "grey92",
      name = "Rainfall relative\nrisk"
    ) +
    labs(
      title = "S11 Spatiotemporal Rainfall Effect",
      subtitle = "Relative change in expected dengue cases for a one-SD increase in lagged rainfall",
      caption = "Municipality component from additive municipality + RW1 time rainfall slopes. Values above 1 indicate higher expected cases."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9),
      legend.position = "right"
    )

  ggsave(MAP_PNG, p, width = 9.5, height = 7, dpi = 300)
  p
}

plot_rainfall_time_effect <- function(effects) {
  p <- ggplot(effects, aes(x = date, y = rainfall_relative_risk)) +
    geom_hline(yintercept = 1, linewidth = 0.35, color = "grey45") +
    geom_ribbon(
      aes(ymin = rainfall_relative_risk_q025, ymax = rainfall_relative_risk_q975),
      fill = "#9ecae1",
      alpha = 0.35
    ) +
    geom_line(color = "#08519c", linewidth = 0.65) +
    labs(
      title = "S11 Time Component of Rainfall Effect",
      subtitle = "Relative change in expected dengue cases for a one-SD increase in lagged rainfall",
      x = "Week",
      y = "Rainfall relative risk",
      caption = "RW1 time component from additive spatiotemporal rainfall-slope model."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9)
    )

  ggsave(TIME_PNG, p, width = 10, height = 5.5, dpi = 300)
  p
}

prepare_s11_dataframe <- function(df) {
  week_levels <- sort(unique(df$week))
  year_levels <- sort(unique(df$year))

  df[, week_idx := match(week, week_levels)]
  df[, year_idx := match(year, year_levels)]

  spatial_lookup <- unique(df[, .(municipio, ibge_code)])
  setorder(spatial_lookup, municipio)

  graph_info <- write_inla_graph(spatial_lookup$ibge_code, GRAPH_FILE)
  spatial_lookup <- merge(
    spatial_lookup,
    graph_info$lookup,
    by = "ibge_code",
    all.x = TRUE,
    sort = FALSE
  )

  spatial_lookup[, isolated_idx := NA_integer_]
  isolated_rows <- which(spatial_lookup$is_isolated)
  if (length(isolated_rows) > 0) {
    spatial_lookup[isolated_rows, isolated_idx := seq_along(isolated_rows)]
  }

  df <- merge(df, spatial_lookup, by = c("municipio", "ibge_code"), all.x = TRUE)
  setorder(df, municipio, date)

  slope_objects <- add_spacetime_slope_indices(df)
  list(
    df = slope_objects$df,
    municipality_lookup = slope_objects$municipality_lookup,
    time_lookup = slope_objects$time_lookup
  )
}


# =========================================================
# Main
# =========================================================
main <- function() {
  df <- build_model_dataframe()
  prepared <- prepare_s11_dataframe(df)
  df <- prepared$df
  municipality_lookup <- prepared$municipality_lookup
  time_lookup <- prepared$time_lookup

  full_dt <- standardize_full(copy(df), BASE_COVARIATES)
  full_fit <- fit_s11_inla(full_dt, GRAPH_FILE)

  cat("\nS11 full-data fixed effects:\n")
  print(full_fit$summary.fixed)
  average_rainfall_effect <- save_fixed_effect_outputs(
    full_fit,
    "S11",
    FIXED_EFFECTS_CSV,
    RAINFALL_AVERAGE_EFFECT_CSV
  )

  cat("\nS11 full-data model criteria:\n")
  criteria <- data.table(
    model = "S11",
    rainfall_specification = "additive municipality + RW1 time rainfall slopes",
    dic = full_fit$dic$dic,
    waic = full_fit$waic$waic
  )
  print(criteria)
  fwrite(criteria, CRITERIA_CSV)

  municipality_effects <- summarize_municipality_rainfall_effects(full_fit, municipality_lookup)
  time_effects <- summarize_time_rainfall_effects(full_fit, time_lookup)
  combined_effects <- summarize_combined_spacetime_effects(full_fit, municipality_lookup, time_lookup)
  fwrite(municipality_effects, EFFECTS_CSV)
  fwrite(time_effects, TIME_EFFECTS_CSV)
  fwrite(combined_effects, COMBINED_EFFECTS_CSV)
  full_predictions <- copy(full_dt[, .(municipio, year, week, date, cases)])
  full_predictions[, predicted_cases := full_fit$summary.fitted.values$mean]
  fwrite(full_predictions, FULL_PREDICTIONS_CSV)
  plot_rainfall_municipality_map(municipality_effects)
  plot_rainfall_time_effect(time_effects)

  cat("\nS11 spatiotemporal rainfall outputs written:\n")
  cat("Municipality CSV:", EFFECTS_CSV, "\n")
  cat("Time CSV:", TIME_EFFECTS_CSV, "\n")
  cat("Combined municipality-week CSV:", COMBINED_EFFECTS_CSV, "\n")
  cat("Fixed effects:", FIXED_EFFECTS_CSV, "\n")
  cat("Average rainfall effect:", RAINFALL_AVERAGE_EFFECT_CSV, "\n")
  cat("Criteria:", CRITERIA_CSV, "\n")
  cat("Full-data fitted predictions:", FULL_PREDICTIONS_CSV, "\n")
  cat("Municipality map:", MAP_PNG, "\n")
  cat("Time plot:", TIME_PNG, "\n")

  if (RUN_TRAIN_TEST_EVALUATION) {
    train_dt <- df[year >= TRAIN_START_YEAR & year <= TRAIN_END_YEAR]
    test_dt <- df[year >= TEST_START_YEAR & year <= TEST_END_YEAR]

    test_start_date <- min(test_dt$date)
    test_rows_before_lag_filter <- nrow(test_dt)
    test_dt <- test_dt[
      cases_lag_source_date < test_start_date &
        neighbor_lag_source_date < test_start_date
    ]
    dropped_test_lag_rows <- test_rows_before_lag_filter - nrow(test_dt)

    scaled <- standardize_train_test(copy(train_dt), copy(test_dt), BASE_COVARIATES)
    train_dt <- scaled$train
    test_dt <- scaled$test

    train_fit <- fit_s11_inla(train_dt, GRAPH_FILE)
    train_pred <- train_fit$summary.fitted.values$mean
    test_pred <- predict_s11_mean(train_fit, test_dt)

    train_metrics <- compute_metrics(train_dt$cases, train_pred)
    train_metrics[, split := "train"]
    test_metrics <- compute_metrics(test_dt$cases, test_pred)
    test_metrics[, split := "test"]
    metrics <- rbindlist(list(train_metrics, test_metrics), use.names = TRUE)
    metrics[, model := "S11"]
    setcolorder(metrics, c("model", "split", "mae", "rmse", "wape", "accuracy_pct", "r2"))

    cat("\nTrain/test evaluation split:\n")
    cat("Train rows:", nrow(train_dt), "\n")
    cat("Test rows:", nrow(test_dt), "\n")
    cat("Dropped test rows with own/neighbor lagged cases inside test period:", dropped_test_lag_rows, "\n")
    cat("No-leakage policy: scaler fit on training data only.\n")
    cat("No-leakage policy: test own-case and neighbor-case lags must come from before test start.\n")

    cat("\nS11 train/test metrics:\n")
    print(metrics)
    fwrite(metrics, METRICS_CSV)
    output <- copy(test_dt[, .(municipio, year, week, date, cases)])
    output[, predicted_cases := test_pred]
    fwrite(output, TEST_PREDICTIONS_CSV)
    cat("Metrics:", METRICS_CSV, "\n")
    cat("Test predictions:", TEST_PREDICTIONS_CSV, "\n")
  }

  cat("\nS11 spatiotemporal rainfall effects:\n")
  print(municipality_effects)
  print(time_effects)
}

main()
