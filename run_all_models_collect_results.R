# =========================================================
# Run all R-INLA model files and collect model-comparison results
# =========================================================
# Default behavior runs each model script, then combines model criteria and
# held-out train/test metrics into publication-ready CSV files. Set
# RUN_MODEL_SCRIPTS=0 to collect existing outputs without refitting.

suppressPackageStartupMessages({
  library(data.table)
})


# =========================================================
# Paths and options
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
  BASE_DIR <- normalizePath(PROJECT_DIR_OVERRIDE)
} else if (dir.exists(file.path(getwd(), "data"))) {
  BASE_DIR <- normalizePath(getwd())
} else if (dir.exists(file.path(SCRIPT_DIR, "data"))) {
  BASE_DIR <- normalizePath(SCRIPT_DIR)
} else {
  BASE_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
}

OUTPUT_DIR <- file.path(BASE_DIR, "outputs")
LOG_DIR <- file.path(OUTPUT_DIR, "run_logs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

RUN_MODEL_SCRIPTS <- Sys.getenv("RUN_MODEL_SCRIPTS", "1") != "0"
RSCRIPT_BIN <- Sys.getenv("RSCRIPT", "Rscript")

COMBINED_CRITERIA_CSV <- file.path(OUTPUT_DIR, "all_model_criteria.csv")
COMBINED_METRICS_CSV <- file.path(OUTPUT_DIR, "all_model_train_test_metrics.csv")
RESULTS_TABLE_CSV <- file.path(OUTPUT_DIR, "all_model_results_table.csv")
RUN_STATUS_CSV <- file.path(OUTPUT_DIR, "all_model_run_status.csv")


# =========================================================
# Model registry
# =========================================================
model_registry <- data.table(
  model = c(
    "R_M0", "R_M1", "R_M2", "R_M3", "R_M4", "R_M5",
    "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11"
  ),
  model_group = c(
    rep("non_spatial", 6),
    rep("spatial", 11)
  ),
  script = c(
    "base_model_r_m0.R",
    "base_model_r_m1_covariates.R",
    "base_model_r_m2_lag_weather.R",
    "base_model_r_m3_interpolation.R",
    "base_model_r_m4_lag_cases.R",
    "base_model_r_m5_lag_weather_cases.R",
    file.path("spatial_R", "spatial_inla_model_s1.R"),
    file.path("spatial_R", "spatial_inla_model_s2.R"),
    file.path("spatial_R", "spatial_inla_model_s3.R"),
    file.path("spatial_R", "spatial_inla_model_s4_road.R"),
    file.path("spatial_R", "spatial_inla_model_s5_air.R"),
    file.path("spatial_R", "spatial_inla_model_s6_rainfall_region.R"),
    file.path("spatial_R", "spatial_inla_model_s7_temperature_region.R"),
    file.path("spatial_R", "spatial_inla_model_s8_rainfall_municipality.R"),
    file.path("spatial_R", "spatial_inla_model_s9_temperature_municipality.R"),
    file.path("spatial_R", "spatial_inla_model_s10_rainfall_time.R"),
    file.path("spatial_R", "spatial_inla_model_s11_rainfall_spacetime.R")
  ),
  criteria_file = c(
    "r_m0_model_criteria.csv",
    "r_m1_model_criteria.csv",
    "r_m2_model_criteria.csv",
    "r_m3_model_criteria.csv",
    "r_m4_model_criteria.csv",
    "r_m5_model_criteria.csv",
    "spatial_inla_s1_model_criteria.csv",
    "spatial_inla_s2_model_criteria.csv",
    "spatial_inla_s3_model_criteria.csv",
    "spatial_inla_s4_road_model_criteria.csv",
    "spatial_inla_s5_air_model_criteria.csv",
    "spatial_inla_s6_rainfall_region_model_criteria.csv",
    "spatial_inla_s7_temperature_region_model_criteria.csv",
    "s8_rainfall_municipality_model_criteria.csv",
    "s9_temperature_municipality_model_criteria.csv",
    "s10_rainfall_time_model_criteria.csv",
    "s11_rainfall_spacetime_model_criteria.csv"
  ),
  metrics_file = c(
    "r_m0_train_test_metrics.csv",
    "r_m1_train_test_metrics.csv",
    "r_m2_train_test_metrics.csv",
    "r_m3_train_test_metrics.csv",
    "r_m4_train_test_metrics.csv",
    "r_m5_train_test_metrics.csv",
    "spatial_inla_s1_train_test_metrics.csv",
    "spatial_inla_s2_train_test_metrics.csv",
    "spatial_inla_s3_train_test_metrics.csv",
    "spatial_inla_s4_road_train_test_metrics.csv",
    "spatial_inla_s5_air_train_test_metrics.csv",
    "spatial_inla_s6_rainfall_region_train_test_metrics.csv",
    "spatial_inla_s7_temperature_region_train_test_metrics.csv",
    "s8_rainfall_municipality_train_test_metrics.csv",
    "s9_temperature_municipality_train_test_metrics.csv",
    "s10_rainfall_time_train_test_metrics.csv",
    "s11_rainfall_spacetime_train_test_metrics.csv"
  )
)

model_registry[, script_path := file.path(BASE_DIR, script)]
model_registry[, criteria_path := file.path(OUTPUT_DIR, criteria_file)]
model_registry[, metrics_path := file.path(OUTPUT_DIR, metrics_file)]
model_registry[, log_path := file.path(LOG_DIR, paste0(model, ".log"))]


# =========================================================
# Running and collecting helpers
# =========================================================
run_model_script <- function(row) {
  model <- row$model
  script_path <- row$script_path
  log_path <- row$log_path

  if (!file.exists(script_path)) {
    return(data.table(
      model = model,
      script = row$script,
      status = "missing_script",
      exit_status = NA_integer_,
      log_path = log_path
    ))
  }

  cat("\n=========================================================\n")
  cat("Running", model, "\n")
  cat("Script:", script_path, "\n")
  cat("Log:", log_path, "\n")
  cat("=========================================================\n")

  start_time <- Sys.time()
  output <- system2(
    RSCRIPT_BIN,
    args = script_path,
    stdout = TRUE,
    stderr = TRUE
  )
  exit_status <- attr(output, "status")
  if (is.null(exit_status)) {
    exit_status <- 0L
  }
  writeLines(output, log_path)

  data.table(
    model = model,
    script = row$script,
    status = if (exit_status == 0L) "ok" else "failed",
    exit_status = as.integer(exit_status),
    started_at = as.character(start_time),
    finished_at = as.character(Sys.time()),
    log_path = log_path
  )
}

read_optional_csv <- function(path, model, kind) {
  if (!file.exists(path)) {
    return(data.table(model = model, missing_file = basename(path), result_type = kind))
  }
  dt <- fread(path)
  if (nrow(dt) == 0) {
    return(data.table(model = model, missing_file = basename(path), result_type = kind))
  }
  if (!("model" %in% names(dt))) {
    dt[, model := model]
    setcolorder(dt, c("model", setdiff(names(dt), "model")))
  }
  dt[, result_type := kind]
  dt[, source_file := basename(path)]
  dt
}

collect_results <- function() {
  criteria <- rbindlist(
    lapply(seq_len(nrow(model_registry)), function(i) {
      row <- model_registry[i]
      read_optional_csv(row$criteria_path, row$model, "criteria")
    }),
    fill = TRUE
  )

  metrics <- rbindlist(
    lapply(seq_len(nrow(model_registry)), function(i) {
      row <- model_registry[i]
      read_optional_csv(row$metrics_path, row$model, "metrics")
    }),
    fill = TRUE
  )

  criteria <- merge(
    model_registry[, .(model, model_group)],
    criteria,
    by = "model",
    all.y = TRUE,
    sort = FALSE
  )
  metrics <- merge(
    model_registry[, .(model, model_group)],
    metrics,
    by = "model",
    all.y = TRUE,
    sort = FALSE
  )

  if (!("split" %in% names(metrics))) {
    metrics[, split := NA_character_]
  }
  test_metrics <- metrics[split == "test"]
  needed_criteria_cols <- c("dic", "waic")
  for (col in needed_criteria_cols) {
    if (!(col %in% names(criteria))) {
      criteria[, (col) := NA_real_]
    }
  }
  needed_metric_cols <- c("mae", "rmse", "wape", "accuracy_pct", "r2")
  for (col in needed_metric_cols) {
    if (!(col %in% names(test_metrics))) {
      test_metrics[, (col) := NA_real_]
    }
  }

  table <- merge(
    criteria[, .(model, model_group, dic, waic)],
    test_metrics[, .(model, mae, rmse, wape, accuracy_pct, r2)],
    by = "model",
    all = TRUE,
    sort = FALSE
  )
  table <- merge(
    model_registry[, .(model, model_group, script)],
    table,
    by = c("model", "model_group"),
    all.x = TRUE,
    sort = FALSE
  )

  if ("waic" %in% names(table)) {
    table[, waic_rank := frank(waic, ties.method = "min", na.last = "keep")]
  }
  if ("wape" %in% names(table)) {
    table[, wape_rank := frank(wape, ties.method = "min", na.last = "keep")]
  }
  if ("rmse" %in% names(table)) {
    table[, rmse_rank := frank(rmse, ties.method = "min", na.last = "keep")]
  }

  fwrite(criteria, COMBINED_CRITERIA_CSV)
  fwrite(metrics, COMBINED_METRICS_CSV)
  fwrite(table, RESULTS_TABLE_CSV)

  cat("\nCombined outputs written:\n")
  cat("Criteria:", COMBINED_CRITERIA_CSV, "\n")
  cat("Metrics:", COMBINED_METRICS_CSV, "\n")
  cat("Results table:", RESULTS_TABLE_CSV, "\n")

  invisible(list(criteria = criteria, metrics = metrics, table = table))
}


# =========================================================
# Main
# =========================================================
main <- function() {
  cat("Project:", BASE_DIR, "\n")
  cat("Run model scripts:", RUN_MODEL_SCRIPTS, "\n")

  if (RUN_MODEL_SCRIPTS) {
    run_status <- rbindlist(
      lapply(seq_len(nrow(model_registry)), function(i) run_model_script(model_registry[i])),
      fill = TRUE
    )
    fwrite(run_status, RUN_STATUS_CSV)
    cat("\nRun status:", RUN_STATUS_CSV, "\n")

    failed <- run_status[status != "ok"]
    if (nrow(failed) > 0) {
      cat("\nSome models failed. Collecting available outputs anyway.\n")
      print(failed[, .(model, status, exit_status, log_path)])
    }
  }

  collect_results()
}

main()
