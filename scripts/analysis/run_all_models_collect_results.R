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
  candidates <- normalizePath(
    file.path(SCRIPT_DIR, c("..", "../..", "../../..")),
    mustWork = FALSE
  )
  matches <- candidates[dir.exists(file.path(candidates, "data"))]
  if (length(matches) == 0) {
    stop("Could not find project root. Run from the project root or set HBM_PROJECT_DIR.")
  }
  BASE_DIR <- matches[1]
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
    file.path("models", "r_inla", "baseline", "base_model_r_m0.R"),
    file.path("models", "r_inla", "baseline", "base_model_r_m1_covariates.R"),
    file.path("models", "r_inla", "baseline", "base_model_r_m2_lag_weather.R"),
    file.path("models", "r_inla", "baseline", "base_model_r_m3_interpolation.R"),
    file.path("models", "r_inla", "baseline", "base_model_r_m4_lag_cases.R"),
    file.path("models", "r_inla", "baseline", "base_model_r_m5_lag_weather_cases.R"),
    file.path("models", "r_inla", "spatial", "spatial_inla_model_s1.R"),
    file.path("models", "r_inla", "spatial", "spatial_inla_model_s2.R"),
    file.path("models", "r_inla", "spatial", "spatial_inla_model_s3.R"),
    file.path("models", "r_inla", "spatial", "spatial_inla_model_s4_road.R"),
    file.path("models", "r_inla", "spatial", "spatial_inla_model_s5_air.R"),
    file.path("models", "r_inla", "spatial", "spatial_inla_model_s6_rainfall_region.R"),
    file.path("models", "r_inla", "spatial", "spatial_inla_model_s7_temperature_region.R"),
    file.path("models", "r_inla", "spatial", "spatial_inla_model_s8_rainfall_municipality.R"),
    file.path("models", "r_inla", "spatial", "spatial_inla_model_s9_temperature_municipality.R"),
    file.path("models", "r_inla", "spatial", "spatial_inla_model_s10_rainfall_time.R"),
    file.path("models", "r_inla", "spatial", "spatial_inla_model_s11_rainfall_spacetime.R")
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
    args = shQuote(script_path),
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

source(file.path(BASE_DIR, "scripts", "analysis", "common_evaluation.R"))

input_fingerprints <- function() {
  files <- sort(unique(c(model_registry$script_path,
    file.path(BASE_DIR, "scripts", "analysis", c("common_evaluation.R", "run_all_models_collect_results.R")),
    list.files(file.path(BASE_DIR, "data"), recursive = TRUE, full.names = TRUE))))
  files <- files[!dir.exists(files)]
  data.table(path = files, md5 = unname(tools::md5sum(files)))
}

collect_results <- function(run_id) {
  run_dir <- file.path(OUTPUT_DIR, "evaluation_runs", run_id)
  provenance_file <- file.path(run_dir, "input_fingerprints.csv")
  if (!file.exists(provenance_file)) stop("Missing run provenance; rerun all models.")
  if (!identical(fread(provenance_file), input_fingerprints())) {
    stop("Code or input data changed since this run began; do not combine stale outputs.")
  }
  manifests <- list()
  predictions <- list()
  for (id in model_registry$model) {
    manifest_file <- file.path(run_dir, paste0(id, "_manifest.csv"))
    prediction_file <- file.path(run_dir, paste0(id, "_predictions.csv"))
    if (!file.exists(manifest_file) || !file.exists(prediction_file)) stop("Missing evaluation bundle: ", id)
    manifest <- fread(manifest_file)
    if (nrow(manifest) != 1L || manifest$model != id || manifest$run_id != run_id ||
        manifest$protocol != EVALUATION_PROTOCOL || manifest$train_start_year != 2017 ||
        manifest$train_end_year != 2022) stop("Incompatible evaluation manifest: ", id)
    if (manifest$prediction_md5 != unname(tools::md5sum(prediction_file))) stop("Prediction file changed: ", id)
    pred <- fread(prediction_file, colClasses = list(character = c("municipio", "date")))
    validate_observations(pred, TRUE)
    if (nrow(pred) != manifest$eligible_test_n || observation_signature(pred) != manifest$eligible_signature) {
      stop("Prediction sample does not match manifest: ", id)
    }
    if (any(as.Date(pred$date) < EVALUATION_START | as.Date(pred$date) > EVALUATION_END)) {
      stop("Evaluation dates outside common protocol: ", id)
    }
    manifests[[id]] <- manifest
    predictions[[id]] <- pred
  }
  manifest <- rbindlist(manifests)
  if (any(!is.finite(manifest$dic)) || any(!is.finite(manifest$waic))) stop("Invalid model criteria.")
  common <- score_common_predictions(predictions)
  metrics <- merge(model_registry[, .(model, model_group)], common$metrics, by = "model", sort = FALSE)
  metrics[, run_id := run_id]
  criteria <- manifest[, .(model, dic, waic, fit_n, fit_start, fit_end, fit_signature)]
  # DIC/WAIC may only be ranked within identical fitted-response cohorts.
  criteria[, criteria_group := paste0("fit_", match(fit_signature, unique(fit_signature)))]
  criteria[, waic_rank_within_fit_sample := frank(waic, ties.method = "min"), by = fit_signature]
  criteria[, dic_rank_within_fit_sample := frank(dic, ties.method = "min"), by = fit_signature]
  table <- merge(metrics, criteria, by = "model", sort = FALSE)
  table <- table[match(model_registry$model, model)]
  table[, wape_rank := frank(wape, ties.method = "min", na.last = "keep")]
  table[, rmse_rank := frank(rmse, ties.method = "min", na.last = "keep")]
  fwrite(common$sample, file.path(run_dir, "common_test_observations.csv"))
  fwrite(common$audit, file.path(run_dir, "common_sample_audit.csv"))
  fwrite(common$predictions, file.path(run_dir, "common_test_predictions.csv"))
  fwrite(criteria, COMBINED_CRITERIA_CSV)
  fwrite(metrics, COMBINED_METRICS_CSV)
  fwrite(table, RESULTS_TABLE_CSV)
  fwrite(table, file.path(run_dir, "results_table.csv"))
  writeLines(run_id, file.path(OUTPUT_DIR, "current_evaluation_run.txt"))
  fwrite(data.table(run_id = run_id, status = "complete"), file.path(OUTPUT_DIR, "evaluation_status.csv"))
  cat("\nCommon held-out observations:", nrow(common$sample), "\n")
  print(common$audit)
  print(table[, .(model, test_n, mae, rmse, wape, r2, criteria_group)])
  cat("\nResults:", RESULTS_TABLE_CSV, "\n")
  invisible(table)
}

main <- function() {
  cat("Project:", BASE_DIR, "\n")
  if (RUN_MODEL_SCRIPTS) {
    run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid())
    Sys.setenv(HBM_EVALUATION_RUN = run_id, HBM_PROJECT_DIR = BASE_DIR)
    run_dir <- file.path(OUTPUT_DIR, "evaluation_runs", run_id)
    dir.create(run_dir, recursive = TRUE)
    archive <- file.path(run_dir, "previous_aggregate_outputs")
    dir.create(archive)
    old <- c(COMBINED_CRITERIA_CSV, COMBINED_METRICS_CSV, RESULTS_TABLE_CSV, RUN_STATUS_CSV)
    file.copy(old[file.exists(old)], archive)
    fwrite(input_fingerprints(), file.path(run_dir, "input_fingerprints.csv"))
    fwrite(data.table(run_id = run_id, status = "running"), file.path(OUTPUT_DIR, "evaluation_status.csv"))
    succeeded <- FALSE
    on.exit({
      if (!succeeded) fwrite(data.table(run_id = run_id, status = "failed"), file.path(OUTPUT_DIR, "evaluation_status.csv"))
    })
    for (i in seq_len(nrow(model_registry))) {
      status <- run_model_script(model_registry[i])
      fwrite(status, file.path(run_dir, "run_status.csv"), append = i > 1L)
      if (status$status != "ok") stop("Model failed: ", status$model, ". See ", status$log_path)
    }
    file.copy(file.path(run_dir, "run_status.csv"), RUN_STATUS_CSV, overwrite = TRUE)
    collect_results(run_id)
    succeeded <- TRUE
  } else {
    status_file <- file.path(OUTPUT_DIR, "evaluation_status.csv")
    if (!file.exists(status_file) || fread(status_file)$status != "complete") {
      stop("Latest evaluation is incomplete or failed. Run make all-results.")
    }
    pointer <- file.path(OUTPUT_DIR, "current_evaluation_run.txt")
    if (!file.exists(pointer)) stop("No verified common-sample run. Run make all-results first.")
    run_id <- readLines(pointer, warn = FALSE)[1]
    if (!grepl("^[A-Za-z0-9_-]+$", run_id)) stop("Invalid saved run ID.")
    collect_results(run_id)
  }
}

main()
