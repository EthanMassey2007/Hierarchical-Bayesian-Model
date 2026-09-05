# =========================================================
# Map S7 region-specific temperature effects
# =========================================================
# Fits S7 if needed, extracts the estimated temperature effect by current IBGE
# intermediate geographic region, joins it to data/RJ.json, and maps the
# temperature relative risk.
#
# Run from the project root:
#   Rscript models/r_inla/spatial/map_s7_temperature_region_effect.R
#
# Outputs:
#   outputs/s7_temperature_region_effects.csv
#   outputs/s7_temperature_region_effect_map.png

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

DATA_DIR_PROJECT <- file.path(PROJECT_DIR, "data")
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
S7_SCRIPT <- file.path(PROJECT_DIR, "models", "r_inla", "spatial", "spatial_inla_model_s7_temperature_region.R")
RJ_GEOJSON <- file.path(DATA_DIR_PROJECT, "RJ.json")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

EFFECTS_CSV <- file.path(OUTPUT_DIR, "s7_temperature_region_effects.csv")
MAP_PNG <- file.path(OUTPUT_DIR, "s7_temperature_region_effect_map.png")


# =========================================================
# Load S7 functions without running S7 main()
# =========================================================
Sys.setenv(INLA_RUN_MODEL = "0")
source(S7_SCRIPT)

# Ensure path globals resolve to project-root data/ after sourcing.
BASE_DIR <- PROJECT_DIR
DATA_DIR <- DATA_DIR_PROJECT
OUTPUT_DIR <- file.path(BASE_DIR, "outputs")
COMBINED_FILE <- file.path(DATA_DIR, "complete_combined_datasets.csv")
MUNICIPIOS_FILE <- file.path(DATA_DIR, "municipios.csv")
HUB_FILE <- file.path(DATA_DIR, "hub_pop_density.csv")
REGION_LOOKUP_FILE <- file.path(DATA_DIR, "rj_ibge_intermediate_regions.csv")
ADJACENCY_FILE <- file.path(DATA_DIR, "adjacency_matrix_correct.parquet")
GRAPH_FILE <- file.path(tempdir(), "rj_municipality_inla.graph")


# =========================================================
# Fit S7 and extract region-specific temperature effects
# =========================================================
build_s7_full_fit <- function() {
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
  full_dt <- add_temperature_region_interactions(full_dt)
  full_fit <- fit_spatial_inla(full_dt, GRAPH_FILE)

  list(
    fit = full_fit,
    effects = summarize_region_temperature_effects(full_fit),
    region_lookup = unique(df[, .(ibge_code, municipio, region)])
  )
}


# =========================================================
# Map
# =========================================================
plot_temperature_region_map <- function(effects, region_lookup) {
  rj <- st_read(RJ_GEOJSON, quiet = TRUE)
  rj$ibge_code <- as.integer(rj$GEOCODIGO)

  map_lookup <- merge(region_lookup, effects, by = "region", all.x = TRUE)

  missing_effects <- map_lookup[is.na(temperature_relative_risk)]
  if (nrow(missing_effects) > 0) {
    stop(sprintf(
      "Missing temperature effects for regions: %s",
      paste(unique(missing_effects$region), collapse = ", ")
    ))
  }

  map_sf <- merge(rj, map_lookup, by = "ibge_code", all.x = TRUE, sort = FALSE)

  missing_geo_effects <- map_sf[is.na(map_sf$temperature_relative_risk), ]
  if (nrow(missing_geo_effects) > 0) {
    cat("Warning: GeoJSON contains municipalities outside the S7 model set:", nrow(missing_geo_effects), "\n")
  }

  p <- ggplot(map_sf) +
    geom_sf(aes(fill = temperature_relative_risk), color = "grey70", linewidth = 0.12) +
    scale_fill_gradient2(
      low = "#2c7bb6",
      mid = "white",
      high = "#d7191c",
      midpoint = 1,
      na.value = "grey92",
      name = "Temperature relative\nrisk"
    ) +
    labs(
      title = "S7 Region-Specific Temperature Effect",
      subtitle = "Relative change in expected dengue cases for a one-SD increase in lagged temperature",
      caption = "Temperature-by-IBGE-intermediate-region interactions. Values above 1 indicate higher expected cases."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9),
      legend.position = "right"
    )

  ggsave(MAP_PNG, p, width = 9.5, height = 7, dpi = 300)
  p
}


# =========================================================
# Main
# =========================================================
main <- function() {
  fit_objects <- build_s7_full_fit()
  effects <- fit_objects$effects
  region_lookup <- fit_objects$region_lookup

  fwrite(effects, EFFECTS_CSV)
  plot_temperature_region_map(effects, region_lookup)

  cat("\nS7 temperature-region map outputs written:\n")
  cat("CSV:", EFFECTS_CSV, "\n")
  cat("Map:", MAP_PNG, "\n\n")
  print(effects)
}

main()
