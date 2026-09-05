# =========================================================
# Map S2 unexplained spatial effects
# =========================================================
# Fits the S2 R-INLA model, extracts the posterior mean BYM2
# municipality spatial effect, converts it to residual spatial
# relative risk, and maps it with data/RJ.json.

suppressPackageStartupMessages({
  library(INLA)
  library(data.table)
  library(sf)
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
} else if (dir.exists(file.path(getwd(), "data"))) {
  PROJECT_DIR <- normalizePath(getwd())
} else {
  candidates <- normalizePath(
    file.path(SCRIPT_DIR, c("..", "../..", "../../..")),
    mustWork = FALSE
  )
  matches <- candidates[dir.exists(file.path(candidates, "data"))]
  if (length(matches) == 0) {
    stop("Could not find project root. Run from the project root or set HBM_PROJECT_DIR.")
  }
  PROJECT_DIR <- matches[1]
}
DATA_DIR_PROJECT <- file.path(PROJECT_DIR, "data")
S2_SCRIPT <- file.path(PROJECT_DIR, "models", "r_inla", "spatial", "spatial_inla_model_s2.R")
RJ_GEOJSON <- file.path(DATA_DIR_PROJECT, "RJ.json")
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

EFFECTS_CSV <- file.path(OUTPUT_DIR, "s2_unexplained_spatial_effects.csv")
MAP_PNG <- file.path(OUTPUT_DIR, "s2_unexplained_spatial_relative_risk_map.png")


# =========================================================
# Load S2 functions without running S2 main()
# =========================================================
Sys.setenv(INLA_RUN_MODEL = "0")
source(S2_SCRIPT)

# spatial_inla_model_s2.R is now inside models/r_inla/spatial/, while data/ is at
# project root. Reset S2's path globals after sourcing.
BASE_DIR <- PROJECT_DIR
DATA_DIR <- DATA_DIR_PROJECT
COMBINED_FILE <- file.path(DATA_DIR, "complete_combined_datasets.csv")
MUNICIPIOS_FILE <- file.path(DATA_DIR, "municipios.csv")
HUB_FILE <- file.path(DATA_DIR, "hub_pop_density.csv")
ADJACENCY_FILE <- file.path(DATA_DIR, "adjacency_matrix_correct.parquet")
GRAPH_FILE <- file.path(tempdir(), "rj_municipality_inla.graph")


# =========================================================
# Fit S2 and extract random spatial effects
# =========================================================
build_s2_full_fit <- function() {
  df <- build_model_dataframe()

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

  full_dt <- standardize_full(copy(df), BASE_COVARIATES)
  full_fit <- fit_spatial_inla(full_dt, GRAPH_FILE)

  list(fit = full_fit, spatial_lookup = spatial_lookup)
}

extract_spatial_effects <- function(fit, spatial_lookup) {
  bym2_summary <- as.data.table(fit$summary.random$spatial_idx)
  bym2_summary[, spatial_idx := as.integer(ID)]
  bym2_summary <- bym2_summary[, .(
    spatial_idx,
    bym2_mean = mean,
    bym2_sd = sd,
    bym2_q025 = `0.025quant`,
    bym2_q975 = `0.975quant`
  )]

  effects <- merge(spatial_lookup, bym2_summary, by = "spatial_idx", all.x = TRUE)

  if (!is.null(fit$summary.random$isolated_idx)) {
    isolated_summary <- as.data.table(fit$summary.random$isolated_idx)
    isolated_summary[, isolated_idx := as.integer(ID)]
    isolated_summary <- isolated_summary[, .(
      isolated_idx,
      isolated_mean = mean,
      isolated_sd = sd,
      isolated_q025 = `0.025quant`,
      isolated_q975 = `0.975quant`
    )]
    effects <- merge(effects, isolated_summary, by = "isolated_idx", all.x = TRUE)
  } else {
    effects[, `:=`(
      isolated_mean = NA_real_,
      isolated_sd = NA_real_,
      isolated_q025 = NA_real_,
      isolated_q975 = NA_real_
    )]
  }

  effects[, spatial_effect_mean := fifelse(is_isolated, isolated_mean, bym2_mean)]
  effects[, spatial_effect_sd := fifelse(is_isolated, isolated_sd, bym2_sd)]
  effects[, spatial_effect_q025 := fifelse(is_isolated, isolated_q025, bym2_q025)]
  effects[, spatial_effect_q975 := fifelse(is_isolated, isolated_q975, bym2_q975)]

  missing_effects <- effects[is.na(spatial_effect_mean)]
  if (nrow(missing_effects) > 0) {
    stop(sprintf(
      "Missing extracted spatial effects for IBGE codes: %s",
      paste(missing_effects$ibge_code, collapse = ", ")
    ))
  }

  effects[, residual_spatial_rr := exp(spatial_effect_mean)]
  effects[, residual_spatial_rr_q025 := exp(spatial_effect_q025)]
  effects[, residual_spatial_rr_q975 := exp(spatial_effect_q975)]

  setorder(effects, municipio)
  effects[, .(
    municipio,
    ibge_code,
    is_isolated,
    spatial_effect_mean,
    spatial_effect_sd,
    spatial_effect_q025,
    spatial_effect_q975,
    residual_spatial_rr,
    residual_spatial_rr_q025,
    residual_spatial_rr_q975
  )]
}


# =========================================================
# Map
# =========================================================
plot_spatial_effect_map <- function(effects) {
  rj <- st_read(RJ_GEOJSON, quiet = TRUE)
  rj$ibge_code <- as.integer(rj$GEOCODIGO)

  map_dt <- merge(
    rj,
    effects,
    by = "ibge_code",
    all.x = TRUE,
    sort = FALSE
  )

  missing_geo <- effects[!(ibge_code %in% map_dt$ibge_code)]
  if (nrow(missing_geo) > 0) {
    stop(sprintf(
      "RJ.json is missing model municipalities: %s",
      paste(missing_geo$ibge_code, collapse = ", ")
    ))
  }

  p <- ggplot(map_dt) +
    geom_sf(aes(fill = residual_spatial_rr), color = "grey70", linewidth = 0.12) +
    scale_fill_gradient2(
      low = "#2c7bb6",
      mid = "white",
      high = "#d7191c",
      midpoint = 1,
      na.value = "grey92",
      name = "Residual spatial\nrelative risk"
    ) +
    labs(
      title = "S2 Residual Spatial Relative Risk",
      subtitle = "exp(posterior mean BYM2 spatial effect); adjusted for lagged weather, IDHM, own-case lag, and neighboring-case lag",
      caption = "Values above 1 indicate higher unexplained spatial dengue risk; values below 1 indicate lower unexplained spatial risk."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  ggsave(MAP_PNG, p, width = 9, height = 7, dpi = 300)
  p
}


# =========================================================
# Main
# =========================================================
main <- function() {
  fit_objects <- build_s2_full_fit()
  effects <- extract_spatial_effects(fit_objects$fit, fit_objects$spatial_lookup)

  fwrite(effects, EFFECTS_CSV)
  plot_spatial_effect_map(effects)

  cat("\nS2 unexplained spatial effects written:\n")
  cat("CSV:", EFFECTS_CSV, "\n")
  cat("Map:", MAP_PNG, "\n")
  cat("\nTop higher unexplained spatial relative risks:\n")
  print(head(effects[order(-residual_spatial_rr)], 10))
  cat("\nTop lower unexplained spatial relative risks:\n")
  print(head(effects[order(residual_spatial_rr)], 10))
}

main()
