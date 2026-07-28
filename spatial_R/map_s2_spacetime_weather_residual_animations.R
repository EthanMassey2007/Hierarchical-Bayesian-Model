# =========================================================
# Animated S2 weekly maps: rainfall, temperature, residuals
# =========================================================
# Fits S2, computes standardized negative-binomial Pearson residuals, and
# creates weekly animated maps for:
#   1. lagged rainfall used by S2
#   2. lagged temperature used by S2
#   3. unexplained S2 residuals
#
# Run from the project root:
#   Rscript spatial_R/map_s2_spacetime_weather_residual_animations.R
#
# If this file is outside the project folder, run:
#   HBM_PROJECT_DIR=/Users/ethanmassey/VS_Code_Test/Hierarchical-Bayesian-Model Rscript path/to/map_s2_spacetime_weather_residual_animations.R
#
# Outputs:
#   outputs/s2_spacetime_weather_residuals_by_municipio_week.csv
#   outputs/s2_spacetime_rainfall_lag_animation_2023.gif
#   outputs/s2_spacetime_temperature_lag_animation_2023.gif
#   outputs/s2_spacetime_residual_animation_2023.gif

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
  library(sf)
  library(ggplot2)
})

if (!requireNamespace("gganimate", quietly = TRUE)) {
  stop("Package 'gganimate' is required. Install it with: install.packages('gganimate')")
}

if (!requireNamespace("gifski", quietly = TRUE)) {
  stop("Package 'gifski' is required. Install it with: install.packages('gifski')")
}

if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package 'scales' is required. Install it with: install.packages('scales')")
}


# =========================================================
# Configuration
# =========================================================
ANIMATION_YEAR <- 2023
FPS <- 2
FRAME_WIDTH <- 1000
FRAME_HEIGHT <- 800
FRAME_RES <- 120


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
} else if (dir.exists(file.path(getwd(), "data")) && dir.exists(file.path(getwd(), "spatial_R"))) {
  PROJECT_DIR <- normalizePath(getwd())
} else if (dir.exists(file.path(SCRIPT_DIR, "data"))) {
  PROJECT_DIR <- normalizePath(SCRIPT_DIR)
} else {
  PROJECT_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
}

DATA_DIR_PROJECT <- file.path(PROJECT_DIR, "data")
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
S2_SCRIPT <- file.path(PROJECT_DIR, "spatial_R", "spatial_inla_model_s2.R")
RJ_GEOJSON <- file.path(DATA_DIR_PROJECT, "RJ.json")
GRAPH_FILE_ANIMATION <- file.path(tempdir(), "rj_municipality_inla_spacetime_weather_residual.graph")

SPACETIME_CSV <- file.path(OUTPUT_DIR, "s2_spacetime_weather_residuals_by_municipio_week.csv")
RAINFALL_GIF <- file.path(
  OUTPUT_DIR,
  sprintf("s2_spacetime_rainfall_lag_animation_%s.gif", ANIMATION_YEAR)
)
TEMPERATURE_GIF <- file.path(
  OUTPUT_DIR,
  sprintf("s2_spacetime_temperature_lag_animation_%s.gif", ANIMATION_YEAR)
)
RESIDUAL_GIF <- file.path(
  OUTPUT_DIR,
  sprintf("s2_spacetime_residual_animation_%s.gif", ANIMATION_YEAR)
)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# =========================================================
# Load S2 functions without running S2 main()
# =========================================================
Sys.setenv(INLA_RUN_MODEL = "0")
source(S2_SCRIPT)

# Keep sourced S2 globals pointed at the project root.
BASE_DIR <- PROJECT_DIR
DATA_DIR <- DATA_DIR_PROJECT
OUTPUT_DIR <- file.path(BASE_DIR, "outputs")
COMBINED_FILE <- file.path(DATA_DIR, "complete_combined_datasets.csv")
MUNICIPIOS_FILE <- file.path(DATA_DIR, "municipios.csv")
HUB_FILE <- file.path(DATA_DIR, "hub_pop_density.csv")
ADJACENCY_FILE <- file.path(DATA_DIR, "adjacency_matrix_correct.parquet")
GRAPH_FILE <- GRAPH_FILE_ANIMATION


# =========================================================
# Residual helpers
# =========================================================
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

compute_spacetime_table <- function(model_dt, fit) {
  n <- nrow(model_dt)
  mu <- as.numeric(fit$summary.fitted.values$mean[seq_len(n)])
  y <- as.numeric(model_dt$cases)
  size <- extract_nb_size(fit)

  if (length(mu) != n || any(!is.finite(mu)) || any(mu <= 0)) {
    stop("Invalid fitted means from S2 fit.")
  }

  variance <- mu + (mu^2 / size)
  pearson <- (y - mu) / sqrt(variance)

  out <- copy(model_dt[, .(
    municipio,
    ibge_code,
    year,
    week,
    date,
    cases,
    rainfall_lag,
    temperature_lag
  )])
  out[, fitted_mu := mu]
  out[, nb_size := size]
  out[, nb_variance := variance]
  out[, pearson_residual := pearson]
  out[, residual_direction := fifelse(
    pearson_residual > 0,
    "underpredicted",
    fifelse(pearson_residual < 0, "overpredicted", "matched")
  )]

  out
}

build_s2_full_fit <- function() {
  df <- build_model_dataframe()

  week_levels <- sort(unique(df$week))
  year_levels <- sort(unique(df$year))

  df[, week_idx := match(week, week_levels)]
  df[, year_idx := match(year, year_levels)]

  spatial_lookup <- unique(df[, .(municipio, ibge_code)])
  setorder(spatial_lookup, municipio)

  graph_info <- write_inla_graph(spatial_lookup$ibge_code, GRAPH_FILE_ANIMATION)
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
  full_fit <- fit_spatial_inla(full_dt, GRAPH_FILE_ANIMATION)

  list(data = full_dt, fit = full_fit)
}


# =========================================================
# Animation helpers
# =========================================================
build_map_sf <- function(spacetime_dt) {
  map_year <- spacetime_dt[year == ANIMATION_YEAR]
  if (nrow(map_year) == 0) {
    stop(sprintf("No rows found for ANIMATION_YEAR = %s.", ANIMATION_YEAR))
  }

  rj <- st_read(RJ_GEOJSON, quiet = TRUE)
  rj$ibge_code <- as.integer(rj$GEOCODIGO)

  map_sf <- merge(rj, map_year, by = "ibge_code", all.x = FALSE, sort = FALSE)
  map_sf$week_label <- sprintf("Year %s - Week %02d", map_sf$year, map_sf$week)
  map_sf
}

save_fill_animation <- function(map_sf, fill_col, title, subtitle, legend_name, output_file, midpoint = NULL) {
  values <- map_sf[[fill_col]]
  if (all(is.na(values))) {
    stop(sprintf("All values are NA for %s.", fill_col))
  }

  if (is.null(midpoint)) {
    fill_scale <- scale_fill_viridis_c(
      option = "C",
      na.value = "grey92",
      name = legend_name
    )
  } else {
    limit <- quantile(abs(values), 0.98, na.rm = TRUE)
    limit <- max(2, as.numeric(limit))
    fill_scale <- scale_fill_gradient2(
      low = "#2c7bb6",
      mid = "white",
      high = "#d7191c",
      midpoint = midpoint,
      limits = c(-limit, limit),
      oob = scales::squish,
      na.value = "grey92",
      name = legend_name
    )
  }

  p <- ggplot(map_sf) +
    geom_sf(aes(fill = .data[[fill_col]]), color = "grey70", linewidth = 0.10) +
    fill_scale +
    labs(
      title = paste0(title, ": {closest_state}"),
      subtitle = subtitle,
      caption = "Weekly municipality-level panel for Rio de Janeiro."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9),
      legend.position = "right"
    ) +
    gganimate::transition_states(week_label, transition_length = 1, state_length = 2) +
    gganimate::ease_aes("linear")

  animation <- gganimate::animate(
    p,
    fps = FPS,
    width = FRAME_WIDTH,
    height = FRAME_HEIGHT,
    res = FRAME_RES,
    renderer = gganimate::gifski_renderer(output_file)
  )

  invisible(animation)
}


# =========================================================
# Main
# =========================================================
main <- function() {
  fit_objects <- build_s2_full_fit()
  spacetime_dt <- compute_spacetime_table(fit_objects$data, fit_objects$fit)
  fwrite(spacetime_dt, SPACETIME_CSV)

  map_sf <- build_map_sf(spacetime_dt)

  save_fill_animation(
    map_sf = map_sf,
    fill_col = "rainfall_lag",
    title = "Lagged Rainfall Used By S2",
    subtitle = "Rainfall value from the S2 weather lag period",
    legend_name = "Lagged\nrainfall",
    output_file = RAINFALL_GIF
  )

  save_fill_animation(
    map_sf = map_sf,
    fill_col = "temperature_lag",
    title = "Lagged Temperature Used By S2",
    subtitle = "Temperature value from the S2 weather lag period",
    legend_name = "Lagged\ntemperature",
    output_file = TEMPERATURE_GIF
  )

  save_fill_animation(
    map_sf = map_sf,
    fill_col = "pearson_residual",
    title = "S2 Unexplained Space-Time Residuals",
    subtitle = "Red = observed cases above fitted expectation; blue = observed cases below fitted expectation",
    legend_name = "Standardized\nresidual",
    output_file = RESIDUAL_GIF,
    midpoint = 0
  )

  cat("\nS2 space-time animation outputs written:\n")
  cat("CSV:", SPACETIME_CSV, "\n")
  cat("Rainfall GIF:", RAINFALL_GIF, "\n")
  cat("Temperature GIF:", TEMPERATURE_GIF, "\n")
  cat("Unexplained residual GIF:", RESIDUAL_GIF, "\n")
  cat("\nInterpretation:\n")
  cat("Rainfall and temperature GIFs show the lagged climate covariates used by S2.\n")
  cat("Residual GIF shows where S2 underpredicts or overpredicts each week.\n")
  cat("Moving residual clusters suggest space-time structure not fully captured by S2.\n")
}

main()
