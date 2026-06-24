# =========================================================
# S2 Moran's I diagnostic for unexplained spatial effects
# =========================================================
# Uses the S2 unexplained spatial effects CSV and the model adjacency matrix
# to test whether residual spatial effects remain spatially autocorrelated.
#
# Run from the project root:
#   Rscript spatial_R/s2_morans_i_diagnostic.R
#
# Inputs:
#   outputs/s2_unexplained_spatial_effects.csv
#   data/adjacency_matrix_correct.parquet
#
# Output:
#   outputs/s2_morans_i_diagnostic.csv

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
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

PROJECT_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
DATA_DIR <- file.path(PROJECT_DIR, "data")
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")

EFFECTS_CSV <- file.path(OUTPUT_DIR, "s2_unexplained_spatial_effects.csv")
ADJACENCY_FILE <- file.path(DATA_DIR, "adjacency_matrix_correct.parquet")
MORANS_OUTPUT <- file.path(OUTPUT_DIR, "s2_morans_i_diagnostic.csv")

PERMUTATIONS <- 999L
RANDOM_SEED <- 42L


# =========================================================
# Helpers
# =========================================================
read_effects <- function() {
  if (!file.exists(EFFECTS_CSV)) {
    stop(sprintf(
      "Missing %s. Run spatial_R/map_s2_unexplained_effects.R first.",
      EFFECTS_CSV
    ))
  }

  effects <- fread(EFFECTS_CSV)
  required <- c("ibge_code", "spatial_effect_mean", "residual_spatial_rr")
  missing <- setdiff(required, names(effects))
  if (length(missing) > 0) {
    stop(sprintf("Effects CSV is missing columns: %s", paste(missing, collapse = ", ")))
  }

  effects[, ibge_code := as.integer(ibge_code)]
  effects <- effects[!is.na(ibge_code)]
  effects
}

read_adjacency_for_effects <- function(ibge_codes) {
  adj <- as.data.table(read_parquet(ADJACENCY_FILE))

  if (!("co_muni_ori" %in% names(adj))) {
    stop("adjacency_matrix_correct.parquet must contain co_muni_ori.")
  }

  adj[, co_muni_ori := as.integer(co_muni_ori)]
  ibge_codes <- as.integer(ibge_codes)
  code_cols <- as.character(ibge_codes)

  missing_rows <- setdiff(ibge_codes, adj$co_muni_ori)
  missing_cols <- setdiff(code_cols, names(adj))
  if (length(missing_rows) > 0 || length(missing_cols) > 0) {
    stop("Adjacency matrix is missing one or more S2 effect municipalities.")
  }

  adj <- adj[match(ibge_codes, co_muni_ori)]
  w <- as.matrix(adj[, ..code_cols])
  storage.mode(w) <- "numeric"
  w[is.na(w)] <- 0
  w[w != 0] <- 1
  diag(w) <- 0

  # Keep adjacency symmetric, matching the S2 BYM2 graph construction.
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
    return(list(observed = NA_real_, p_value = NA_real_))
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

run_morans_diagnostics <- function() {
  effects <- read_effects()
  setorder(effects, ibge_code)
  w <- read_adjacency_for_effects(effects$ibge_code)

  set.seed(RANDOM_SEED)

  spatial_effect_result <- permutation_test(effects$spatial_effect_mean, w)
  rr_result <- permutation_test(effects$residual_spatial_rr, w)

  diagnostics <- rbindlist(list(
    data.table(
      variable = "spatial_effect_mean",
      morans_i = spatial_effect_result$observed,
      expected_random = spatial_effect_result$expected_random,
      permuted_sd = spatial_effect_result$permuted_sd,
      permutation_p_value = spatial_effect_result$p_value,
      permutations = PERMUTATIONS,
      interpretation = "Primary diagnostic on the log-scale BYM2 spatial effect."
    ),
    data.table(
      variable = "residual_spatial_rr",
      morans_i = rr_result$observed,
      expected_random = rr_result$expected_random,
      permuted_sd = rr_result$permuted_sd,
      permutation_p_value = rr_result$p_value,
      permutations = PERMUTATIONS,
      interpretation = "Secondary diagnostic on exp(BYM2 effect), the mapped relative-risk scale."
    )
  ), use.names = TRUE)

  fwrite(diagnostics, MORANS_OUTPUT)

  cat("\nS2 Moran's I diagnostic written:\n")
  cat("CSV:", MORANS_OUTPUT, "\n\n")
  print(diagnostics)

  cat("\nHow to read this:\n")
  cat("Moran's I > 0 means neighboring municipalities tend to have similar residual spatial effects.\n")
  cat("Moran's I near 0 means little residual spatial autocorrelation remains.\n")
  cat("Small permutation p-values suggest the spatial pattern is unlikely under random reassignment.\n")
}


run_morans_diagnostics()
