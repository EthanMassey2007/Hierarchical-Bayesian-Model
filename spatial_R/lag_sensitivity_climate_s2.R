# =========================================================
# Climate lag sensitivity analysis for S2
# =========================================================
# Refit the same S2 negative-binomial BYM2 model across alternative shared
# climate lags. The dengue own-case and neighbor-case lags stay fixed at 4
# weeks so the comparison isolates the climate lag choice. Candidate climate
# lags are compared on a common analytic sample and with rolling-origin yearly
# validation, which is more suitable for publication than a single small
# holdout period.

suppressPackageStartupMessages({
  library(data.table)
})


# =========================================================
# Configuration
# =========================================================
MODEL_ID <- "s2_climate_lag_sensitivity"
CASE_LAG_FIXED_WEEKS <- 4L
DEFAULT_CLIMATE_LAGS <- 4:14
DEFAULT_VALIDATION_YEARS <- 2020:2023
ENFORCE_COMMON_ROWS <- TRUE
PRIMARY_SELECTION_METRIC <- "cv_wape"
VALIDATION_POLICY <- "rolling_origin_one_step_ahead"

parse_lag_values <- function(x, default, label) {
  if (is.na(x) || !nzchar(trimws(x))) {
    return(default)
  }
  values <- as.integer(trimws(unlist(strsplit(x, ","))))
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    stop(sprintf("%s did not contain any integer values.", label))
  }
  sort(unique(values))
}

CLIMATE_LAGS <- parse_lag_values(
  Sys.getenv("CLIMATE_LAGS", unset = NA_character_),
  DEFAULT_CLIMATE_LAGS,
  "CLIMATE_LAGS"
)
VALIDATION_YEARS <- parse_lag_values(
  Sys.getenv("VALIDATION_YEARS", unset = NA_character_),
  DEFAULT_VALIDATION_YEARS,
  "VALIDATION_YEARS"
)


# =========================================================
# Load S2 functions without running S2 main()
# =========================================================
script_arg <- commandArgs(trailingOnly = FALSE)
script_file_arg <- script_arg[grepl("^--file=", script_arg)]
if (length(script_file_arg) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  SCRIPT_DIR <- getwd()
}
S2_SCRIPT <- file.path(SCRIPT_DIR, "spatial_inla_model_s2.R")
if (!file.exists(S2_SCRIPT)) {
  S2_SCRIPT <- file.path(SCRIPT_DIR, "spatial_R", "spatial_inla_model_s2.R")
}
if (!file.exists(S2_SCRIPT)) {
  stop("Could not locate spatial_inla_model_s2.R.")
}

old_inla_run_model <- Sys.getenv("INLA_RUN_MODEL", unset = NA_character_)
Sys.setenv(INLA_RUN_MODEL = "0")
source(S2_SCRIPT)
if (is.na(old_inla_run_model)) {
  Sys.unsetenv("INLA_RUN_MODEL")
} else {
  Sys.setenv(INLA_RUN_MODEL = old_inla_run_model)
}


# =========================================================
# Sensitivity helpers
# =========================================================
extract_lpml <- function(fit) {
  cpo_values <- fit$cpo$cpo
  cpo_values <- cpo_values[is.finite(cpo_values) & cpo_values > 0]
  if (length(cpo_values) == 0) {
    return(NA_real_)
  }
  sum(log(cpo_values))
}

make_row_key <- function(dt) {
  paste(dt$municipio, dt$date, sep = "|")
}

build_common_row_keys <- function(climate_lags) {
  if (!ENFORCE_COMMON_ROWS) {
    return(NULL)
  }

  cat("\nBuilding common municipality-week sample across candidate lags...\n")
  key_sets <- lapply(climate_lags, function(lag_value) {
    WEATHER_LAG_WEEKS <<- as.integer(lag_value)
    CASE_LAG_WEEKS <<- CASE_LAG_FIXED_WEEKS
    lag_dt <- build_model_dataframe()
    unique(make_row_key(lag_dt))
  })

  common_keys <- Reduce(intersect, key_sets)
  if (length(common_keys) == 0) {
    stop("No common rows are available across the requested climate lags.")
  }

  cat("Common rows available to every lag:", length(common_keys), "\n")
  common_keys
}

prepare_spatial_data <- function(df) {
  week_levels <- sort(unique(df$week))
  year_levels <- sort(unique(df$year))

  df[, week_idx := match(week, week_levels)]
  df[, year_idx := match(year, year_levels)]

  spatial_lookup <- unique(df[, .(municipio, ibge_code)])
  setorder(spatial_lookup, municipio)

  graph_info <- write_inla_graph(spatial_lookup$ibge_code)
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
  df
}

fit_validation_fold <- function(df, validation_year) {
  train_dt <- df[year < validation_year]
  test_dt <- df[year == validation_year]

  if (nrow(train_dt) == 0 || nrow(test_dt) == 0) {
    stop(sprintf("Validation year %s does not have both train and test rows.", validation_year))
  }

  scaled <- standardize_train_test(copy(train_dt), copy(test_dt), BASE_COVARIATES)
  train_dt <- scaled$train
  test_dt <- scaled$test

  fold_fit <- fit_spatial_inla(train_dt, GRAPH_FILE)
  fold_fit$.model_data <- train_dt
  fold_fit$.formula <- build_inla_formula(GRAPH_FILE)

  train_pred <- fold_fit$summary.fitted.values$mean
  test_pred <- predict_inla_mean(fold_fit, test_dt)

  train_metrics <- compute_metrics(train_dt$cases, train_pred)
  test_metrics <- compute_metrics(test_dt$cases, test_pred)
  test_abs_error <- abs(as.numeric(test_dt$cases) - as.numeric(test_pred))
  test_sq_error <- (as.numeric(test_dt$cases) - as.numeric(test_pred))^2

  data.table(
    validation_year = as.integer(validation_year),
    train_year_start = min(train_dt$year),
    train_year_end = max(train_dt$year),
    train_rows = nrow(train_dt),
    test_rows = nrow(test_dt),
    test_cases_sum = sum(as.numeric(test_dt$cases)),
    test_abs_error_sum = sum(test_abs_error),
    test_sq_error_sum = sum(test_sq_error),
    train_mae = train_metrics$mae,
    train_rmse = train_metrics$rmse,
    train_wape = train_metrics$wape,
    train_accuracy_pct = train_metrics$accuracy_pct,
    train_r2 = train_metrics$r2,
    test_mae = test_metrics$mae,
    test_rmse = test_metrics$rmse,
    test_wape = test_metrics$wape,
    test_accuracy_pct = test_metrics$accuracy_pct,
    test_r2 = test_metrics$r2
  )
}

summarize_validation_folds <- function(folds) {
  folds[, .(
    cv_folds = .N,
    cv_train_rows_min = min(train_rows),
    cv_train_rows_max = max(train_rows),
    cv_test_rows_total = sum(test_rows),
    cv_test_cases_total = sum(test_cases_sum),
    cv_mae = sum(test_abs_error_sum) / sum(test_rows),
    cv_rmse = sqrt(sum(test_sq_error_sum) / sum(test_rows)),
    cv_wape = sum(test_abs_error_sum) / max(sum(test_cases_sum), 1e-9),
    cv_accuracy_pct = max(0, 100 * (1 - sum(test_abs_error_sum) / max(sum(test_cases_sum), 1e-9))),
    cv_r2 = weighted.mean(test_r2, test_rows)
  )]
}

fit_one_climate_lag <- function(climate_lag_weeks, common_row_keys = NULL) {
  cat("\n=========================================================\n")
  cat("Fitting S2 climate lag sensitivity model\n")
  cat("Climate lag:", climate_lag_weeks, "weeks\n")
  cat("Case lag:", CASE_LAG_FIXED_WEEKS, "weeks\n")
  cat("=========================================================\n")

  WEATHER_LAG_WEEKS <<- as.integer(climate_lag_weeks)
  CASE_LAG_WEEKS <<- CASE_LAG_FIXED_WEEKS
  SAVE_OUTPUTS <<- FALSE
  RUN_TRAIN_TEST_EVALUATION <<- TRUE

  df <- build_model_dataframe()
  n_available_rows <- nrow(df)
  if (!is.null(common_row_keys)) {
    df <- df[make_row_key(df) %in% common_row_keys]
    setorder(df, municipio, date)
  }
  df <- prepare_spatial_data(df)

  full_dt <- standardize_full(copy(df), BASE_COVARIATES)
  full_fit <- fit_spatial_inla(full_dt, GRAPH_FILE)

  fold_results <- rbindlist(lapply(VALIDATION_YEARS, function(validation_year) {
    fit_validation_fold(df, validation_year)
  }))
  fold_results[, climate_lag_weeks := as.integer(climate_lag_weeks)]
  setcolorder(fold_results, c("climate_lag_weeks", setdiff(names(fold_results), "climate_lag_weeks")))

  cv_summary <- summarize_validation_folds(fold_results)

  summary <- data.table(
    model = MODEL_ID,
    climate_lag_weeks = as.integer(climate_lag_weeks),
    case_lag_weeks = CASE_LAG_FIXED_WEEKS,
    validation_years = paste(VALIDATION_YEARS, collapse = ","),
    primary_selection_metric = PRIMARY_SELECTION_METRIC,
    validation_policy = VALIDATION_POLICY,
    common_row_filter = !is.null(common_row_keys),
    n_available_rows_before_common_filter = n_available_rows,
    n_model_rows = nrow(df),
    dic = full_fit$dic$dic,
    waic = full_fit$waic$waic,
    lpml = extract_lpml(full_fit)
  )

  list(summary = cbind(summary, cv_summary), folds = fold_results)
}


# =========================================================
# Main
# =========================================================
main <- function() {
  cat("Candidate shared climate lags:", paste(CLIMATE_LAGS, collapse = ", "), "\n")
  cat("Fixed own-case and neighbor-case lag:", CASE_LAG_FIXED_WEEKS, "weeks\n")
  cat("Validation years:", paste(VALIDATION_YEARS, collapse = ", "), "\n")
  cat("Primary selection metric:", PRIMARY_SELECTION_METRIC, "\n")
  cat("Common-row comparison:", ENFORCE_COMMON_ROWS, "\n")

  common_row_keys <- build_common_row_keys(CLIMATE_LAGS)

  lag_outputs <- lapply(CLIMATE_LAGS, function(lag_value) {
    tryCatch(
      fit_one_climate_lag(lag_value, common_row_keys),
      error = function(e) {
        list(
          summary = data.table(
            model = MODEL_ID,
            climate_lag_weeks = as.integer(lag_value),
            case_lag_weeks = CASE_LAG_FIXED_WEEKS,
            validation_years = paste(VALIDATION_YEARS, collapse = ","),
            primary_selection_metric = PRIMARY_SELECTION_METRIC,
            validation_policy = VALIDATION_POLICY,
            common_row_filter = ENFORCE_COMMON_ROWS,
            error = conditionMessage(e)
          ),
          folds = data.table(
            climate_lag_weeks = as.integer(lag_value),
            error = conditionMessage(e)
          )
        )
      }
    )
  })

  results <- rbindlist(
    lapply(lag_outputs, function(x) x$summary),
    fill = TRUE
  )
  fold_results <- rbindlist(
    lapply(lag_outputs, function(x) x$folds),
    fill = TRUE
  )

  if (!("error" %in% names(results))) {
    results[, error := NA_character_]
  }
  if (!("error" %in% names(fold_results))) {
    fold_results[, error := NA_character_]
  }

  successful <- results[is.na(error)]
  if (nrow(successful) > 0) {
    successful[, waic_rank := frank(waic, ties.method = "min")]
    successful[, dic_rank := frank(dic, ties.method = "min")]
    successful[, cv_rmse_rank := frank(cv_rmse, ties.method = "min")]
    successful[, cv_wape_rank := frank(cv_wape, ties.method = "min")]
    successful[, best_by_waic := waic_rank == 1]
    successful[, best_by_dic := dic_rank == 1]
    successful[, best_by_cv_rmse := cv_rmse_rank == 1]
    successful[, best_by_cv_wape := cv_wape_rank == 1]

    results <- merge(
      results,
      successful[, .(
        climate_lag_weeks,
        waic_rank,
        dic_rank,
        cv_rmse_rank,
        cv_wape_rank,
        best_by_waic,
        best_by_dic,
        best_by_cv_rmse,
        best_by_cv_wape
      )],
      by = "climate_lag_weeks",
      all.x = TRUE,
      sort = FALSE
    )
    setorder(results, cv_wape_rank, cv_rmse_rank, waic_rank, climate_lag_weeks)
  }

  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  summary_file <- file.path(OUTPUT_DIR, "s2_climate_lag_sensitivity.csv")
  folds_file <- file.path(OUTPUT_DIR, "s2_climate_lag_sensitivity_folds.csv")
  table_file <- file.path(OUTPUT_DIR, "s2_climate_lag_sensitivity_table.csv")
  fwrite(results, summary_file)
  fwrite(fold_results, folds_file)

  table_cols <- c(
    "climate_lag_weeks",
    "case_lag_weeks",
    "n_model_rows",
    "cv_folds",
    "cv_test_rows_total",
    "cv_wape",
    "cv_accuracy_pct",
    "cv_rmse",
    "waic",
    "dic",
    "cv_wape_rank",
    "cv_rmse_rank",
    "best_by_cv_wape"
  )
  available_table_cols <- intersect(table_cols, names(results))
  manuscript_table <- copy(results[, ..available_table_cols])
  fwrite(manuscript_table, table_file)

  cat("\nClimate lag sensitivity summary:\n")
  print(results)
  cat("\nWrote:", summary_file, "\n")
  cat("Wrote:", folds_file, "\n")
  cat("Wrote:", table_file, "\n")

  if (nrow(successful) > 0) {
    best <- successful[which.min(cv_wape)]
    cat(
      "\nBest shared climate lag by cross-validated WAPE:",
      best$climate_lag_weeks,
      "weeks\n"
    )
  }
}

main()
