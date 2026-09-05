# =========================================================
# S10: S2 + time-varying rainfall effect using R-INLA
# =========================================================
# Uses the S2 model/data pipeline and adds a time-varying random slope
# for lagged rainfall. This estimates when rainfall sensitivity changes
# over the study period while partially pooling adjacent weeks.
#
# Run from the project root:
#   Rscript models/r_inla/spatial/spatial_inla_model_s10_rainfall_time.R
#
# Outputs:
#   outputs/s10_rainfall_time_effects.csv
#   outputs/s10_rainfall_time_effect_plot.png

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
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
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

EFFECTS_CSV <- file.path(OUTPUT_DIR, "s10_rainfall_time_effects.csv")
FIXED_EFFECTS_CSV <- file.path(OUTPUT_DIR, "s10_rainfall_time_fixed_effects.csv")
RAINFALL_AVERAGE_EFFECT_CSV <- file.path(OUTPUT_DIR, "s10_rainfall_time_average_rainfall_effect.csv")
CRITERIA_CSV <- file.path(OUTPUT_DIR, "s10_rainfall_time_model_criteria.csv")
METRICS_CSV <- file.path(OUTPUT_DIR, "s10_rainfall_time_train_test_metrics.csv")
EFFECT_PNG <- file.path(OUTPUT_DIR, "s10_rainfall_time_effect_plot.png")


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
GRAPH_FILE <- file.path(tempdir(), "rj_municipality_inla_s10.graph")


# =========================================================
# S10 helpers
# =========================================================
add_time_slope_index <- function(df) {
  slope_lookup <- unique(df[, .(date, year, week)])
  setorder(slope_lookup, date)
  slope_lookup[, rainfall_time_idx := .I]

  df <- merge(df, slope_lookup, by = c("date", "year", "week"), all.x = TRUE)
  setorder(df, municipio, date)

  list(df = df, slope_lookup = slope_lookup)
}

build_s10_formula <- function(graph_file) {
  z_covariates <- paste0(BASE_COVARIATES, "_z")
  formula <- as.formula(paste(
    "cases ~ 1 +",
    paste(z_covariates, collapse = " + "),
    "+ f(rainfall_time_idx, rainfall_lag_z, model = 'rw1', scale.model = TRUE, constr = TRUE)",
    "+ f(week_idx, model = 'iid')",
    "+ f(year_idx, model = 'iid')",
    "+ f(spatial_idx, model = 'bym2', graph = graph_file, scale.model = TRUE)",
    "+ f(isolated_idx, model = 'iid')"
  ))
  environment(formula) <- environment()
  formula
}

fit_s10_inla <- function(model_dt, graph_file) {
  formula <- build_s10_formula(graph_file)
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

predict_s10_mean <- function(fit, new_dt) {
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

summarize_time_rainfall_effects <- function(fit, slope_lookup) {
  fixed <- as.data.table(fit$summary.fixed, keep.rownames = "term")
  rainfall <- fixed[term == "rainfall_lag_z"]
  if (nrow(rainfall) != 1) {
    stop("Could not find rainfall_lag_z in fixed effects.")
  }

  random <- as.data.table(fit$summary.random$rainfall_time_idx)
  random[, rainfall_time_idx := as.integer(ID)]
  effects <- merge(slope_lookup, random, by = "rainfall_time_idx", all.x = TRUE)

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
      title = "S10 Time-Varying Rainfall Effect",
      subtitle = "Relative change in expected dengue cases for a one-SD increase in lagged rainfall",
      x = "Week",
      y = "Rainfall relative risk",
      caption = "RW1 time-varying rainfall slope. Values above 1 indicate higher expected cases."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9)
    )

  ggsave(EFFECT_PNG, p, width = 10, height = 5.5, dpi = 300)
  p
}

prepare_s10_dataframe <- function(df) {
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

  slope_objects <- add_time_slope_index(df)
  list(
    df = slope_objects$df,
    slope_lookup = slope_objects$slope_lookup
  )
}


# =========================================================
# Main
# =========================================================
main <- function() {
  df <- build_model_dataframe()
  prepared <- prepare_s10_dataframe(df)
  df <- prepared$df
  slope_lookup <- prepared$slope_lookup

  full_dt <- standardize_full(copy(df), BASE_COVARIATES)
  full_fit <- fit_s10_inla(full_dt, GRAPH_FILE)

  cat("\nS10 full-data fixed effects:\n")
  print(full_fit$summary.fixed)
  average_rainfall_effect <- save_fixed_effect_outputs(
    full_fit,
    "S10",
    FIXED_EFFECTS_CSV,
    RAINFALL_AVERAGE_EFFECT_CSV
  )

  cat("\nS10 full-data model criteria:\n")
  criteria <- data.table(
    model = "S10",
    rainfall_specification = "time-varying RW1 rainfall slope",
    dic = full_fit$dic$dic,
    waic = full_fit$waic$waic
  )
  print(criteria)
  fwrite(criteria, CRITERIA_CSV)

  full_effects <- summarize_time_rainfall_effects(full_fit, slope_lookup)
  fwrite(full_effects, EFFECTS_CSV)
  plot_rainfall_time_effect(full_effects)

  cat("\nS10 time-varying rainfall outputs written:\n")
  cat("CSV:", EFFECTS_CSV, "\n")
  cat("Fixed effects:", FIXED_EFFECTS_CSV, "\n")
  cat("Average rainfall effect:", RAINFALL_AVERAGE_EFFECT_CSV, "\n")
  cat("Criteria:", CRITERIA_CSV, "\n")
  cat("Plot:", EFFECT_PNG, "\n")

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

    train_fit <- fit_s10_inla(train_dt, GRAPH_FILE)
    train_pred <- train_fit$summary.fitted.values$mean
    test_pred <- predict_s10_mean(train_fit, test_dt)

    train_metrics <- compute_metrics(train_dt$cases, train_pred)
    train_metrics[, split := "train"]
    test_metrics <- compute_metrics(test_dt$cases, test_pred)
    test_metrics[, split := "test"]
    metrics <- rbindlist(list(train_metrics, test_metrics), use.names = TRUE)
    metrics[, model := "S10"]
    setcolorder(metrics, c("model", "split", "mae", "rmse", "wape", "accuracy_pct", "r2"))

    cat("\nTrain/test evaluation split:\n")
    cat("Train rows:", nrow(train_dt), "\n")
    cat("Test rows:", nrow(test_dt), "\n")
    cat("Dropped test rows with own/neighbor lagged cases inside test period:", dropped_test_lag_rows, "\n")
    cat("No-leakage policy: scaler fit on training data only.\n")
    cat("No-leakage policy: test own-case and neighbor-case lags must come from before test start.\n")

    cat("\nS10 train/test metrics:\n")
    print(metrics)
    fwrite(metrics, METRICS_CSV)
    cat("Metrics:", METRICS_CSV, "\n")
  }

  cat("\nS10 time-specific rainfall effects:\n")
  print(full_effects)
}

main()
