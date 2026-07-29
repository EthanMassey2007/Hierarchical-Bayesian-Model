# =========================================================
# R-M0: non-spatial baseline using R-INLA
# =========================================================
# Converted from base_model.py.
# Uses no fixed covariates.

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
})

MODEL_ID <- "R_M0"
MODEL_DESCRIPTION <- "base only"

DATA_START_YEAR <- 2017
DATA_END_YEAR <- 2023
TRAIN_START_YEAR <- 2017
TRAIN_END_YEAR <- 2022
TEST_START_YEAR <- 2023
TEST_END_YEAR <- 2023

CASE_LAG_WEEKS <- 4
WEATHER_LAG_WEEKS <- 6
INLA_NUM_THREADS <- "4:1"
inla.setOption(num.threads = INLA_NUM_THREADS)

USE_SAME_WEEK_WEATHER <- FALSE
USE_WEATHER_LAGS <- FALSE
USE_CASE_LAG <- FALSE

BASE_COVARIATES <- character()

WEATHER_COVARIATES <- c("rainfall", "humidity", "temperature")

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

DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

COMBINED_FILE <- file.path(DATA_DIR, "complete_combined_datasets.csv")
CRITERIA_CSV <- file.path(OUTPUT_DIR, "r_m0_model_criteria.csv")
METRICS_CSV <- file.path(OUTPUT_DIR, "r_m0_train_test_metrics.csv")
FIXED_EFFECTS_CSV <- file.path(OUTPUT_DIR, "r_m0_fixed_effects.csv")

clean_columns <- function(dt) {
  setnames(dt, trimws(tolower(names(dt))))
  dt
}

normalize_name <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x <- trimws(tolower(x))
  x <- gsub("/[a-z]{2}$", "", x)
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

iso_week_to_date <- function(year, week) {
  jan4 <- as.Date(sprintf("%d-01-04", year))
  jan4_weekday <- as.integer(format(jan4, "%u"))
  week1_monday <- jan4 - (jan4_weekday - 1)
  week1_monday + (as.integer(week) - 1L) * 7L
}

compute_metrics <- function(y_true, y_pred) {
  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)
  mae <- mean(abs(y_true - y_pred))
  rmse <- sqrt(mean((y_true - y_pred)^2))
  wape <- sum(abs(y_true - y_pred)) / max(sum(abs(y_true)), 1e-9)
  accuracy_pct <- max(0, 100 * (1 - wape))
  sst <- sum((y_true - mean(y_true))^2)
  sse <- sum((y_true - y_pred)^2)
  r2 <- if (sst > 0) 1 - sse / sst else NA_real_
  data.table(mae = mae, rmse = rmse, wape = wape, accuracy_pct = accuracy_pct, r2 = r2)
}

standardize_train_test <- function(train_dt, test_dt, covariates) {
  if (length(covariates) == 0) {
    return(list(train = train_dt, test = test_dt))
  }
  means <- train_dt[, lapply(.SD, mean), .SDcols = covariates]
  sds <- train_dt[, lapply(.SD, sd), .SDcols = covariates]
  for (col in covariates) {
    mu <- means[[col]]
    sigma <- sds[[col]]
    if (is.na(sigma) || sigma == 0) {
      stop(sprintf("Cannot standardize %s because training sd is zero/NA.", col))
    }
    train_dt[, paste0(col, "_z") := (get(col) - mu) / sigma]
    test_dt[, paste0(col, "_z") := (get(col) - mu) / sigma]
  }
  list(train = train_dt, test = test_dt)
}

standardize_full <- function(dt, covariates) {
  if (length(covariates) == 0) {
    return(dt)
  }
  for (col in covariates) {
    mu <- mean(dt[[col]])
    sigma <- sd(dt[[col]])
    if (is.na(sigma) || sigma == 0) {
      stop(sprintf("Cannot standardize %s because sd is zero/NA.", col))
    }
    dt[, paste0(col, "_z") := (get(col) - mu) / sigma]
  }
  dt
}

build_model_dataframe <- function() {
  df <- fread(COMBINED_FILE)
  df <- clean_columns(df)

  required <- c("municipio", "year", "week", "cases")
  if (USE_SAME_WEEK_WEATHER || USE_WEATHER_LAGS) {
    required <- c(required, WEATHER_COVARIATES, "idhm")
  }
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  df[, municipio := normalize_name(municipio)]
  numeric_cols <- intersect(c("year", "week", "cases", "idhm", WEATHER_COVARIATES), names(df))
  for (col in numeric_cols) {
    df[, (col) := as.numeric(get(col))]
  }

  original_rows <- nrow(df)
  df <- df[!is.na(municipio) & !is.na(year) & !is.na(week) & !is.na(cases)]
  df[, year := as.integer(year)]
  df[, week := as.integer(week)]
  df <- df[year >= DATA_START_YEAR & year <= DATA_END_YEAR]
  df[, date := iso_week_to_date(year, week)]
  df <- df[!is.na(date)]
  setorder(df, municipio, date)

  if (USE_WEATHER_LAGS) {
    weather_lookup <- df[, c("municipio", "date", WEATHER_COVARIATES), with = FALSE]
    setnames(weather_lookup, WEATHER_COVARIATES, paste0(WEATHER_COVARIATES, "_lag"))
    weather_lookup[, date := date + WEATHER_LAG_WEEKS * 7L]
    weather_lookup[, weather_lag_source_date := date - WEATHER_LAG_WEEKS * 7L]
    df <- merge(df, weather_lookup, by = c("municipio", "date"), all.x = TRUE)
    bad_weather_lag <- df[
      !is.na(weather_lag_source_date) &
        (weather_lag_source_date != date - WEATHER_LAG_WEEKS * 7L |
           weather_lag_source_date >= date)
    ]
    if (nrow(bad_weather_lag) > 0) {
      stop("Weather lag leakage check failed.")
    }
  }

  if (USE_CASE_LAG) {
    case_lookup <- df[, .(municipio, date, cases_lag = cases)]
    case_lookup[, date := date + CASE_LAG_WEEKS * 7L]
    case_lookup[, cases_lag_source_date := date - CASE_LAG_WEEKS * 7L]
    df <- merge(df, case_lookup, by = c("municipio", "date"), all.x = TRUE)
    bad_case_lag <- df[
      !is.na(cases_lag) &
        (cases_lag_source_date != date - CASE_LAG_WEEKS * 7L |
           cases_lag_source_date >= date)
    ]
    if (nrow(bad_case_lag) > 0) {
      stop("Case lag leakage check failed.")
    }
    df[, log_cases_lag := log1p(cases_lag)]
  }

  df[, cases := pmax(as.integer(round(cases)), 0L)]
  rows_before_drop <- nrow(df)
  if (length(BASE_COVARIATES) > 0) {
    df <- df[complete.cases(df[, ..BASE_COVARIATES])]
  }

  df[, municipio_idx := as.integer(factor(municipio))]
  df[, week_idx := as.integer(factor(week, levels = sort(unique(week))))]
  df[, year_idx := as.integer(factor(year, levels = sort(unique(year))))]
  setorder(df, municipio, date)

  cat(MODEL_ID, "dataframe built\n")
  cat("Original rows:", original_rows, "\n")
  cat("Model rows:", nrow(df), "\n")
  cat("Dropped missing model covariates:", rows_before_drop - nrow(df), "\n")
  cat("Municipios:", uniqueN(df$municipio), "\n")
  cat("Years:", paste(sort(unique(df$year)), collapse = ", "), "\n")
  df
}

build_formula <- function() {
  z_covariates <- if (length(BASE_COVARIATES) > 0) paste0(BASE_COVARIATES, "_z") else character()
  fixed_terms <- if (length(z_covariates) > 0) paste(z_covariates, collapse = " + ") else "1"
  formula_text <- paste(
    "cases ~",
    fixed_terms,
    "+ f(municipio_idx, model = 'iid')",
    "+ f(week_idx, model = 'iid')",
    "+ f(year_idx, model = 'iid')"
  )
  as.formula(formula_text)
}

fit_inla_model <- function(model_dt) {
  formula <- build_formula()
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

predict_inla_mean <- function(fit, new_dt) {
  pred_dt <- copy(new_dt)
  pred_dt[, cases := NA_integer_]
  fit_dt <- rbindlist(list(fit$.model_data, pred_dt), use.names = TRUE, fill = TRUE)
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

criteria_table <- function(fit) {
  cpo_values <- fit$cpo$cpo
  cpo_values <- cpo_values[is.finite(cpo_values) & cpo_values > 0]
  data.table(
    model = MODEL_ID,
    description = MODEL_DESCRIPTION,
    dic = fit$dic$dic,
    waic = fit$waic$waic,
    lpml = sum(log(cpo_values)),
    mean_log_cpo = mean(log(cpo_values))
  )
}

main <- function() {
  df <- build_model_dataframe()
  full_dt <- standardize_full(copy(df), BASE_COVARIATES)
  full_fit <- fit_inla_model(full_dt)

  fixed_effects <- as.data.table(full_fit$summary.fixed, keep.rownames = "term")
  criteria <- criteria_table(full_fit)

  cat("\n", MODEL_ID, " full-data fixed effects:\n", sep = "")
  print(fixed_effects)
  cat("\n", MODEL_ID, " full-data criteria:\n", sep = "")
  print(criteria)
  fwrite(fixed_effects, FIXED_EFFECTS_CSV)
  fwrite(criteria, CRITERIA_CSV)

  train_dt <- df[year >= TRAIN_START_YEAR & year <= TRAIN_END_YEAR]
  test_dt <- df[year >= TEST_START_YEAR & year <= TEST_END_YEAR]
  dropped_test_lag_rows <- 0L
  if (USE_CASE_LAG) {
    test_start_date <- min(test_dt$date)
    rows_before <- nrow(test_dt)
    test_dt <- test_dt[cases_lag_source_date < test_start_date]
    dropped_test_lag_rows <- rows_before - nrow(test_dt)
  }

  scaled <- standardize_train_test(copy(train_dt), copy(test_dt), BASE_COVARIATES)
  train_dt <- scaled$train
  test_dt <- scaled$test
  train_fit <- fit_inla_model(train_dt)
  train_pred <- train_fit$summary.fitted.values$mean
  test_pred <- predict_inla_mean(train_fit, test_dt)

  train_metrics <- compute_metrics(train_dt$cases, train_pred)
  train_metrics[, split := "train"]
  test_metrics <- compute_metrics(test_dt$cases, test_pred)
  test_metrics[, split := "test"]
  metrics <- rbindlist(list(train_metrics, test_metrics), use.names = TRUE)
  setcolorder(metrics, c("split", "mae", "rmse", "wape", "accuracy_pct", "r2"))

  cat("\nTrain/test evaluation split:\n")
  cat("Train rows:", nrow(train_dt), "\n")
  cat("Test rows:", nrow(test_dt), "\n")
  cat("Dropped test rows with lagged cases inside test period:", dropped_test_lag_rows, "\n")
  cat("No-leakage policy: scaler fit on training data only.\n")
  if (USE_CASE_LAG) {
    cat("No-leakage policy: test case lags must come from before test start.\n")
  }
  cat("\n", MODEL_ID, " train/test metrics:\n", sep = "")
  print(metrics)
  fwrite(metrics, METRICS_CSV)

  cat("\nSaved outputs:\n")
  cat("Criteria:", CRITERIA_CSV, "\n")
  cat("Fixed effects:", FIXED_EFFECTS_CSV, "\n")
  cat("Metrics:", METRICS_CSV, "\n")
}

main()
