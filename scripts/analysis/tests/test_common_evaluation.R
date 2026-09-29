suppressPackageStartupMessages(library(data.table))
source("scripts/analysis/common_evaluation.R")
expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
a <- data.table(municipio = c("a", "b", "c"), date = rep("2023-01-02", 3),
                cases = c(2, 4, 100), predicted_cases = c(1, 6, 0))
b <- data.table(municipio = c("b", "a", "d"), date = rep("2023-01-02", 3),
                cases = c(4, 2, 200), predicted_cases = c(4, 2, 0))
result <- score_common_predictions(list(M0 = a, S11 = b))
stopifnot(nrow(result$sample) == 2L, all(result$metrics$test_n == 2L),
          result$metrics[model == "M0", mae] == 1.5,
          result$metrics[model == "M0", wape] == 0.5,
          result$metrics[model == "M0", rmse] == sqrt(2.5),
          result$metrics[model == "M0", r2] == -1.5,
          result$metrics[model == "S11", rmse] == 0,
          all(result$audit$excluded_from_common == 1L))
stopifnot(identical(observation_signature(a), observation_signature(a[3:1])))
expect_error(score_common_predictions(list(A = rbind(a, a[1]), B = b)), "Duplicate")
bad <- copy(b); bad[1, cases := 99]
expect_error(score_common_predictions(list(A = a, B = bad)), "disagree")
bad <- copy(b); bad[1, predicted_cases := NA_real_]
expect_error(score_common_predictions(list(A = a, B = bad)), "Missing")
bad <- copy(b); bad[1, predicted_cases := Inf]
expect_error(score_common_predictions(list(A = a, B = bad)), "Invalid predictions")
expect_error(score_common_predictions(list(A = a, B = b[3])), "No common")
expect_error(restrict_evaluation_window(data.table(date = as.Date("2023-02-01"))), "No eligible")
zeros <- copy(a[1]); zeros[, `:=`(cases = 0, predicted_cases = 0)]
stopifnot(is.na(score_predictions(zeros)$wape), is.na(score_predictions(zeros)$r2))
train <- copy(a); train[, date := "2022-12-05"]
test <- copy(a); test[, cases_lag_source_date := as.Date("2022-12-05")]
folder <- tempfile(); dir.create(folder)
write_evaluation_bundle("M0", train, test, test$predicted_cases, rbind(train, a),
                        folder, data.table(dic=1,waic=2), 2017, 2022)
manifest <- fread(file.path(folder, "evaluation_runs/standalone/M0_manifest.csv"))
stopifnot(manifest$eligible_test_n == 3L, manifest$fit_n == 6L)
test[1, cases_lag_source_date := as.Date("2023-01-02")]
expect_error(write_evaluation_bundle("M0", train, test, test$predicted_cases,
  rbind(train,a), folder, data.table(dic=1,waic=2), 2017, 2022), "unavailable")
cat("Common-sample scoring, invalid input, provenance, and leakage checks passed.\n")

# Integration: the collector must publish aligned scores and reject stale bundles.
collector <- new.env(parent = globalenv())
expressions <- parse("scripts/analysis/run_all_models_collect_results.R")
for (expr in expressions[-length(expressions)]) eval(expr, collector)
collector$OUTPUT_DIR <- tempfile(); dir.create(collector$OUTPUT_DIR)
collector$COMBINED_CRITERIA_CSV <- file.path(collector$OUTPUT_DIR, "criteria.csv")
collector$COMBINED_METRICS_CSV <- file.path(collector$OUTPUT_DIR, "metrics.csv")
collector$RESULTS_TABLE_CSV <- file.path(collector$OUTPUT_DIR, "results.csv")
collector$model_registry <- data.table(model = c("A", "B"), model_group = "fixture")
collector$input_fingerprints <- function() data.table(path = "fixture", md5 = "fixture")
Sys.setenv(HBM_EVALUATION_RUN = "test_run")
write_evaluation_bundle("A", train, a, a$predicted_cases, rbind(train,a),
  collector$OUTPUT_DIR, data.table(dic=1,waic=2), 2017, 2022)
write_evaluation_bundle("B", train, b, b$predicted_cases, rbind(train,a),
  collector$OUTPUT_DIR, data.table(dic=2,waic=3), 2017, 2022)
run_dir <- file.path(collector$OUTPUT_DIR, "evaluation_runs/test_run")
fwrite(collector$input_fingerprints(), file.path(run_dir, "input_fingerprints.csv"))
collected <- collector$collect_results("test_run")
stopifnot(all(collected$test_n == 2L), length(unique(collected$criteria_group)) == 1L,
          !"accuracy_pct" %in% names(collected), !"waic_rank" %in% names(collected))
original_table <- unname(tools::md5sum(collector$RESULTS_TABLE_CSV))
expect_error(collector$collect_results("missing"), "Missing run provenance")
writeLines("altered", file.path(run_dir, "A_predictions.csv"))
expect_error(collector$collect_results("test_run"), "Prediction file changed")
stopifnot(unname(tools::md5sum(collector$RESULTS_TABLE_CSV)) == original_table)
Sys.unsetenv("HBM_EVALUATION_RUN")
cat("Collector integration and stale-output rejection checks passed.\n")
