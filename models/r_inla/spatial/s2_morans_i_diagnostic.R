# =========================================================
# S2 Moran's I diagnostic on standardized NB residuals
# =========================================================
# Fits the S2 negative-binomial INLA model, computes municipality-aggregated
# standardized Pearson and deviance residuals, and tests those residuals for
# remaining spatial autocorrelation with Moran's I.
#
# Run from the project root:
#   Rscript models/r_inla/spatial/s2_morans_i_diagnostic.R
#
# Outputs:
#   outputs/s2_morans_i_diagnostic.csv
#   outputs/s2_standardized_residuals_by_municipio.csv

suppressPackageStartupMessages({
  library(INLA)
  library(arrow)
  library(data.table)
  library(Matrix)
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
} else if (dir.exists(file.path(getwd(), "data"))) {
  PROJECT_DIR <- normalizePath(getwd())
} else if (dir.exists(file.path(SCRIPT_DIR, "data"))) {
  PROJECT_DIR <- normalizePath(SCRIPT_DIR)
} else {
  candidates <- normalizePath(file.path(SCRIPT_DIR, c("..", "../..", "../../..")), mustWork = FALSE)
  matches <- candidates[dir.exists(file.path(candidates, "data"))]
  if (length(matches) == 0) {
    stop("Could not find project root. Run from the project root or set HBM_PROJECT_DIR.")
  }
  PROJECT_DIR <- matches[1]
}

DATA_DIR <- file.path(PROJECT_DIR, "data")
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
S2_SCRIPT <- file.path(PROJECT_DIR, "models", "r_inla", "spatial", "spatial_inla_model_s2.R")
ADJACENCY_FILE_DIAG <- file.path(DATA_DIR, "adjacency_matrix_correct.parquet")
MORANS_OUTPUT <- file.path(OUTPUT_DIR, "s2_morans_i_diagnostic.csv")
RESIDUALS_OUTPUT <- file.path(OUTPUT_DIR, "s2_standardized_residuals_by_municipio.csv")
GRAPH_FILE_DIAG <- file.path(tempdir(), "rj_municipality_inla_residual_diagnostic.graph")

PERMUTATIONS <- 999L
RANDOM_SEED <- 42L

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# =========================================================
# Load S2 functions without running S2 main()
# =========================================================
Sys.setenv(INLA_RUN_MODEL = "0")
source(S2_SCRIPT)

# Keep sourced S2 globals pointed at the project root.
BASE_DIR <- PROJECT_DIR
DATA_DIR <- file.path(BASE_DIR, "data")
COMBINED_FILE <- file.path(DATA_DIR, "complete_combined_datasets.csv")
MUNICIPIOS_FILE <- file.path(DATA_DIR, "municipios.csv")
HUB_FILE <- file.path(DATA_DIR, "hub_pop_density.csv")
ADJACENCY_FILE <- file.path(DATA_DIR, "adjacency_matrix_correct.parquet")
GRAPH_FILE <- GRAPH_FILE_DIAG


# =========================================================
# Helpers
# =========================================================
read_adjacency_for_codes <- function(ibge_codes) {
  adj <- as.data.table(read_parquet(ADJACENCY_FILE_DIAG))

  if (!("co_muni_ori" %in% names(adj))) {
    stop("adjacency_matrix_correct.parquet must contain co_muni_ori.")
  }

  adj[, co_muni_ori := as.integer(co_muni_ori)]
  ibge_codes <- as.integer(ibge_codes)
  code_cols <- as.character(ibge_codes)

  missing_rows <- setdiff(ibge_codes, adj$co_muni_ori)
  missing_cols <- setdiff(code_cols, names(adj))
  if (length(missing_rows) > 0 || length(missing_cols) > 0) {
    stop("Adjacency matrix is missing one or more model municipalities.")
  }

  adj <- adj[match(ibge_codes, co_muni_ori)]
  w <- as.matrix(adj[, ..code_cols])
  storage.mode(w) <- "numeric"
  w[is.na(w)] <- 0
  w[w != 0] <- 1
  diag(w) <- 0

  # Match the S2 BYM2 graph construction.
  w <- ((w + t(w)) > 0) * 1
  rownames(w) <- code_cols
  colnames(w) <- code_cols
  w
}

morans_i <- function(x, w) {
  x <- as.numeric(x)
  if (length(x) != nrow(w)) {
    stop("Length of x must equal number of adjacency rows.")
  }

  valid <- is.finite(x)
  if (!all(valid)) {
    x <- x[valid]
    w <- w[valid, valid, drop = FALSE]
  }

  n <- length(x)
  x_centered <- x - mean(x)
  s0 <- sum(w)
  denom <- sum(x_centered^2)

  if (n < 3 || s0 == 0 || denom == 0) {
    return(NA_real_)
  }

  as.numeric((n / s0) * (sum(w * outer(x_centered, x_centered)) / denom))
}

permutation_test <- function(x, w, permutations = PERMUTATIONS) {
  observed <- morans_i(x, w)
  if (is.na(observed)) {
    return(list(
      observed = NA_real_,
      expected_random = NA_real_,
      permuted_sd = NA_real_,
      p_value = NA_real_
    ))
  }

  permuted <- numeric(permutations)
  for (i in seq_len(permutations)) {
    permuted[i] <- morans_i(sample(x), w)
  }

  # Two-sided permutation p-value with +1 correction.
  p_value <- (sum(abs(permuted) >= abs(observed)) + 1) / (permutations + 1)

  list(
    observed = observed,
    expected_random = mean(permuted),
    permuted_sd = sd(permuted),
    p_value = p_value
  )
}

extract_nb_size <- function(fit) {
  hyper <- as.data.table(fit$summary.hyperpar, keep.rownames = "term")
  size_rows <- hyper[grepl("size", term, ignore.case = TRUE)]

  if (nrow(size_rows) == 0) {
    stop(sprintf(
      "Could not identify the negative-binomial size parameter. Hyperparameters found: %s",
      paste(hyper$term, collapse = ", ")
    ))
  }

  size <- as.numeric(size_rows$mean[1])
  if (!is.finite(size) || size <= 0) {
    stop(sprintf("Invalid negative-binomial size parameter: %s", size))
  }

  size
}

compute_standardized_residuals <- function(model_dt, fit) {
  n <- nrow(model_dt)
  mu <- as.numeric(fit$summary.fitted.values$mean[seq_len(n)])
  y <- as.numeric(model_dt$cases)
  size <- extract_nb_size(fit)

  if (length(mu) != n || any(!is.finite(mu)) || any(mu <= 0)) {
    stop("Invalid fitted means from S2 fit.")
  }

  variance <- mu + (mu^2 / size)
  pearson <- (y - mu) / sqrt(variance)

  # Negative-binomial deviance residual using Var(Y)=mu+mu^2/size.
  y_log_term <- ifelse(y == 0, 0, y * log(y / mu))
  deviance_component <- 2 * (
    y_log_term - (y + size) * log((y + size) / (mu + size))
  )
  deviance_component <- pmax(deviance_component, 0)
  deviance <- sign(y - mu) * sqrt(deviance_component)

  residual_dt <- copy(model_dt[, .(municipio, ibge_code, year, week, date, cases)])
  residual_dt[, fitted_mu := mu]
  residual_dt[, nb_size := size]
  residual_dt[, nb_variance := variance]
  residual_dt[, pearson_residual := pearson]
  residual_dt[, deviance_residual := deviance]

  residual_dt
}

aggregate_residuals_by_municipio <- function(residual_dt) {
  residual_dt[, .(
    n_observations = .N,
    observed_cases = sum(cases),
    fitted_cases = sum(fitted_mu),
    pearson_residual = sum(cases - fitted_mu) / sqrt(sum(nb_variance)),
    mean_pearson_residual = mean(pearson_residual),
    mean_deviance_residual = mean(deviance_residual),
    median_deviance_residual = median(deviance_residual)
  ), by = .(municipio, ibge_code)]
}

build_s2_full_fit <- function() {
  df <- build_model_dataframe()

  week_levels <- sort(unique(df$week))
  year_levels <- sort(unique(df$year))

  df[, week_idx := match(week, week_levels)]
  df[, year_idx := match(year, year_levels)]

  spatial_lookup <- unique(df[, .(municipio, ibge_code)])
  setorder(spatial_lookup, municipio)

  graph_info <- write_inla_graph(spatial_lookup$ibge_code, GRAPH_FILE_DIAG)
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

  full_dt <- standardize_full(copy(df), BASE_COVARIATES)
  full_fit <- fit_spatial_inla(full_dt, GRAPH_FILE_DIAG)

  list(data = full_dt, fit = full_fit)
}

run_morans_diagnostics <- function() {
  fit_objects <- build_s2_full_fit()
  residual_dt <- compute_standardized_residuals(fit_objects$data, fit_objects$fit)
  residual_by_municipio <- aggregate_residuals_by_municipio(residual_dt)
  setorder(residual_by_municipio, ibge_code)

  w <- read_adjacency_for_codes(residual_by_municipio$ibge_code)

  set.seed(RANDOM_SEED)

  pearson_result <- permutation_test(residual_by_municipio$pearson_residual, w)
  mean_pearson_result <- permutation_test(residual_by_municipio$mean_pearson_residual, w)
  deviance_result <- permutation_test(residual_by_municipio$mean_deviance_residual, w)

  diagnostics <- rbindlist(list(
    data.table(
      variable = "municipality_aggregated_pearson_residual",
      morans_i = pearson_result$observed,
      expected_random = pearson_result$expected_random,
      permuted_sd = pearson_result$permuted_sd,
      permutation_p_value = pearson_result$p_value,
      permutations = PERMUTATIONS,
      interpretation = "Primary diagnostic: Moran's I on municipality-aggregated Pearson residuals from the negative-binomial S2 model."
    ),
    data.table(
      variable = "municipality_mean_pearson_residual",
      morans_i = mean_pearson_result$observed,
      expected_random = mean_pearson_result$expected_random,
      permuted_sd = mean_pearson_result$permuted_sd,
      permutation_p_value = mean_pearson_result$p_value,
      permutations = PERMUTATIONS,
      interpretation = "Sensitivity diagnostic: Moran's I on municipality mean Pearson residuals."
    ),
    data.table(
      variable = "municipality_mean_deviance_residual",
      morans_i = deviance_result$observed,
      expected_random = deviance_result$expected_random,
      permuted_sd = deviance_result$permuted_sd,
      permutation_p_value = deviance_result$p_value,
      permutations = PERMUTATIONS,
      interpretation = "Sensitivity diagnostic: Moran's I on municipality mean deviance residuals."
    )
  ), use.names = TRUE)

  fwrite(residual_by_municipio, RESIDUALS_OUTPUT)
  fwrite(diagnostics, MORANS_OUTPUT)

  cat("\nS2 standardized-residual Moran's I diagnostic written:\n")
  cat("Residuals CSV:", RESIDUALS_OUTPUT, "\n")
  cat("Diagnostic CSV:", MORANS_OUTPUT, "\n\n")
  print(diagnostics)

  cat("\nHow to read this:\n")
  cat("Moran's I > 0 means neighboring municipalities tend to have similar standardized residuals.\n")
  cat("Moran's I near 0 means little residual spatial autocorrelation remains.\n")
  cat("Small permutation p-values suggest residual spatial autocorrelation remains after S2.\n")
}


run_morans_diagnostics()
