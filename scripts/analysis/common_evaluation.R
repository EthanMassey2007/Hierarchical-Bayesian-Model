# Shared fixed-window evaluation contract for the manuscript's R-INLA ladder.
# Full-data descriptive fits remain separate. Only held-out predictions are scored.
EVALUATION_START <- as.Date("2023-01-02")
EVALUATION_END <- as.Date("2023-01-23")
EVALUATION_PROTOCOL <- "common_observations_2023_january_v1"

restrict_evaluation_window <- function(dt) {
  result <- copy(dt)[as.Date(date) >= EVALUATION_START & as.Date(date) <= EVALUATION_END]
  if (!nrow(result)) stop("No eligible observations in the fixed evaluation window.")
  result
}

validate_observations <- function(dt, predictions = FALSE) {
  required <- c("municipio", "date", "cases", if (predictions) "predicted_cases")
  if (!all(required %in% names(dt))) stop("Missing observation/prediction columns.")
  if (!nrow(dt)) stop("Empty observation set.")
  if (anyNA(dt[, ..required])) stop("Missing observation/prediction values.")
  if (any(!nzchar(dt$municipio)) || anyNA(as.Date(dt$date))) stop("Invalid observation keys.")
  if (anyDuplicated(dt[, .(municipio, date)])) stop("Duplicate municipality/date keys.")
  if (any(!is.finite(dt$cases) | dt$cases < 0)) stop("Invalid observed counts.")
  if (predictions && any(!is.finite(dt$predicted_cases) | dt$predicted_cases < 0)) {
    stop("Invalid predictions; do not silently remove failed predictions.")
  }
  invisible(TRUE)
}

observation_signature <- function(dt) {
  validate_observations(dt)
  canonical <- copy(dt[, .(municipio, date = as.character(as.Date(date)), cases)])
  setorder(canonical, municipio, date)
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  fwrite(canonical, path)
  unname(tools::md5sum(path))
}

write_evaluation_bundle <- function(model_id, train_dt, test_dt, predictions,
                                    full_dt, output_dir, criteria,
                                    train_start_year, train_end_year) {
  if (train_start_year != 2017 || train_end_year != 2022) stop("Unexpected training years.")
  if (length(predictions) != nrow(test_dt)) stop("Prediction length mismatch.")
  pred <- copy(test_dt[, .(municipio, date = as.character(as.Date(date)), cases)])
  pred[, predicted_cases := as.numeric(predictions)]
  validate_observations(pred, TRUE)
  if (any(as.Date(pred$date) < EVALUATION_START | as.Date(pred$date) > EVALUATION_END)) {
    stop("Predictions outside the declared evaluation window.")
  }
  if (any(as.Date(train_dt$date) >= EVALUATION_START)) stop("Training overlaps evaluation.")
  for (col in intersect(c("cases_lag_source_date", "neighbor_lag_source_date", "distance_lag_source_date"), names(test_dt))) {
    if (anyNA(test_dt[[col]]) || any(as.Date(test_dt[[col]]) >= EVALUATION_START)) {
      stop("Case lag is unavailable at the evaluation origin: ", col)
    }
  }
  run_id <- Sys.getenv("HBM_EVALUATION_RUN", "standalone")
  if (!grepl("^[A-Za-z0-9_-]+$", run_id)) stop("Invalid evaluation run ID.")
  bundle_dir <- file.path(output_dir, "evaluation_runs", run_id)
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  prediction_file <- file.path(bundle_dir, paste0(model_id, "_predictions.csv"))
  fwrite(pred, prediction_file)
  manifest <- data.table(
    model = model_id, run_id = run_id, protocol = EVALUATION_PROTOCOL,
    train_start_year = train_start_year, train_end_year = train_end_year,
    train_n = nrow(train_dt), train_signature = observation_signature(train_dt),
    fit_n = nrow(full_dt), fit_signature = observation_signature(full_dt),
    fit_start = as.character(min(as.Date(full_dt$date))),
    fit_end = as.character(max(as.Date(full_dt$date))),
    eligible_test_n = nrow(pred), eligible_signature = observation_signature(pred),
    prediction_md5 = unname(tools::md5sum(prediction_file)),
    dic = criteria$dic[1], waic = criteria$waic[1]
  )
  fwrite(manifest, file.path(bundle_dir, paste0(model_id, "_manifest.csv")))
}

score_predictions <- function(dt) {
  validate_observations(dt, TRUE)
  error <- dt$cases - dt$predicted_cases
  total <- sum(abs(dt$cases))
  sst <- sum((dt$cases - mean(dt$cases))^2)
  data.table(mae = mean(abs(error)), rmse = sqrt(mean(error^2)),
             wape = if (total > 0) sum(abs(error)) / total else NA_real_,
             r2 = if (sst > 0) 1 - sum(error^2) / sst else NA_real_)
}

score_common_predictions <- function(predictions) {
  if (is.null(names(predictions)) || anyDuplicated(names(predictions))) stop("Models must have unique names.")
  lapply(predictions, validate_observations, predictions = TRUE)
  keys <- Reduce(function(x, y) merge(x, y, by = c("municipio", "date")),
                 lapply(predictions, function(x) x[, .(municipio, date)]))
  if (!nrow(keys)) stop("No common evaluation observations.")
  setorder(keys, municipio, date)
  aligned <- lapply(predictions, function(x) merge(keys, x, by = c("municipio", "date"), sort = TRUE))
  reference <- aligned[[1]]$cases
  if (!all(vapply(aligned, function(x) identical(x$cases, reference), logical(1)))) {
    stop("Observed counts disagree between models on shared keys.")
  }
  sample <- aligned[[1]][, .(municipio, date, cases)]
  signature <- observation_signature(sample)
  metrics <- rbindlist(lapply(names(aligned), function(id) {
    result <- score_predictions(aligned[[id]])
    result[, `:=`(model = id, split = "test", evaluation_protocol = EVALUATION_PROTOCOL,
                   evaluation_signature = signature, test_n = nrow(sample),
                   test_start = min(sample$date), test_end = max(sample$date))]
    result
  }))
  audit <- data.table(model = names(predictions),
                      eligible_test_n = vapply(predictions, nrow, integer(1)),
                      common_test_n = nrow(sample))
  audit[, excluded_from_common := eligible_test_n - common_test_n]
  long <- rbindlist(lapply(names(aligned), function(id) {
    value <- copy(aligned[[id]])
    value[, model := id]
    value
  }))
  list(metrics = metrics, sample = sample, audit = audit, predictions = long)
}

# Per-model diagnostics retain their native row sets; make that scope explicit.
annotate_model_metrics <- function(metrics, train_dt, test_dt) {
  result <- copy(metrics)
  for (part in c("train", "test")) {
    dt <- if (part == "train") train_dt else test_dt
    result[split == part, `:=`(
      evaluation_scope = if (part == "train") "training_fit" else "model_eligible_not_cross_model",
      n_observations = nrow(dt), sample_start = as.character(min(as.Date(dt$date))),
      sample_end = as.character(max(as.Date(dt$date))), sample_signature = observation_signature(dt)
    )]
  }
  result
}
