# =========================================================
# S3: M5 + BYM2 spatial effect + distance-weighted lagged cases using R-INLA
# =========================================================
# This script mirrors the non-spatial M5 feature set:
#   rainfall_lag + humidity_lag + temperature_lag + idhm + log_cases_lag
# and adds:
#   1. a BYM2 spatial random effect over connected municipalities,
#   2. an IID effect for municipalities isolated in the adjacency graph,
#   3. distance_log_cases_lag, a distance-decayed spatial lag of cases.

suppressPackageStartupMessages({
  library(INLA)
  library(arrow)
  library(data.table)
  library(Matrix)
})


# =========================================================
# Helpers
# =========================================================
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


# =========================================================
# Configuration
# =========================================================
DATA_START_YEAR <- 2017
DATA_END_YEAR <- 2023

TRAIN_START_YEAR <- 2017
TRAIN_END_YEAR <- 2022
TEST_START_YEAR <- 2023
TEST_END_YEAR <- 2023

CASE_LAG_WEEKS <- 4
WEATHER_LAG_WEEKS <- 6
DISTANCE_DECAY_KM <- 50
DISTANCE_CUTOFF_KM <- 150

RUN_TRAIN_TEST_EVALUATION <- TRUE
SAVE_OUTPUTS <- FALSE
INLA_NUM_THREADS <- "4:1"
inla.setOption(num.threads = INLA_NUM_THREADS)

BASE_COVARIATES <- c(
  "rainfall_lag",
  "humidity_lag",
  "temperature_lag",
  "idhm",
  "log_cases_lag",
  "distance_log_cases_lag"
)

LAGGED_WEATHER_COVARIATES <- c(
  "rainfall",
  "humidity",
  "temperature"
)


# =========================================================
# Paths
# =========================================================
script_arg <- commandArgs(trailingOnly = FALSE)
script_file_arg <- script_arg[grepl("^--file=", script_arg)]
if (length(script_file_arg) > 0) {
  BASE_DIR <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  BASE_DIR <- getwd()
}
DATA_DIR <- file.path(BASE_DIR, "data")

COMBINED_FILE <- file.path(DATA_DIR, "complete_combined_datasets.csv")
MUNICIPIOS_FILE <- file.path(DATA_DIR, "municipios.csv")
HUB_FILE <- file.path(DATA_DIR, "hub_pop_density.csv")
ADJACENCY_FILE <- file.path(DATA_DIR, "adjacency_matrix_correct.parquet")
CONNECTIVITY_FILE <- file.path(DATA_DIR, "fluvi_road_ibge.parquet")

GRAPH_FILE <- file.path(tempdir(), "rj_municipality_inla.graph")


# =========================================================
# General helpers
# =========================================================
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

  data.table(
    mae = mae,
    rmse = rmse,
    wape = wape,
    accuracy_pct = accuracy_pct,
    r2 = r2
  )
}

standardize_train_test <- function(train_dt, test_dt, covariates) {
  means <- train_dt[, lapply(.SD, mean), .SDcols = covariates]
  sds <- train_dt[, lapply(.SD, sd), .SDcols = covariates]

  for (col in covariates) {
    mu <- means[[col]]
    sigma <- sds[[col]]
    if (is.na(sigma) || sigma == 0) {
      stop(sprintf("Cannot standardize %s because its training sd is zero/NA.", col))
    }
    scaled_col <- paste0(col, "_z")
    train_dt[, (scaled_col) := (get(col) - mu) / sigma]
    test_dt[, (scaled_col) := (get(col) - mu) / sigma]
  }

  list(train = train_dt, test = test_dt)
}

standardize_full <- function(dt, covariates) {
  for (col in covariates) {
    mu <- mean(dt[[col]])
    sigma <- sd(dt[[col]])
    if (is.na(sigma) || sigma == 0) {
      stop(sprintf("Cannot standardize %s because its sd is zero/NA.", col))
    }
    dt[, paste0(col, "_z") := (get(col) - mu) / sigma]
  }
  dt
}


# =========================================================
# Municipality lookup
# =========================================================
build_municipio_lookup <- function() {
  municipios <- fread(MUNICIPIOS_FILE)
  municipios <- clean_columns(municipios)

  if (!all(c("city", "ibgeid") %in% names(municipios))) {
    stop("municipios.csv must contain city and ibgeid columns.")
  }

  if ("state" %in% names(municipios)) {
    municipios[, state_code := toupper(as.character(state))]
  } else {
    municipios[, state_code := toupper(sub(".*\\/([A-Za-z]{2})$", "\\1", city))]
  }
  municipios[, municipio := normalize_name(city)]
  municipios[, ibge_code := as.integer(ibgeid)]
  municipios <- municipios[state_code == "RJ" & !is.na(ibge_code)]
  lookup <- unique(municipios[, .(municipio, ibge_code)])
  lookup <- lookup[!duplicated(municipio)]

  hub <- fread(HUB_FILE)
  hub <- clean_columns(hub)
  if (all(c("uf", "nm_municipio", "co_ibge") %in% names(hub))) {
    hub <- hub[toupper(as.character(uf)) == "RJ"]
    hub[, municipio := normalize_name(nm_municipio)]
    hub[, ibge_code := as.integer(co_ibge)]
    hub_lookup <- unique(hub[!is.na(ibge_code), .(municipio, ibge_code)])
    hub_lookup <- hub_lookup[!duplicated(municipio)]
    lookup <- rbindlist(
      list(
        lookup[, source_priority := 1L],
        hub_lookup[, source_priority := 2L]
      ),
      use.names = TRUE
    )
    setorder(lookup, municipio, source_priority)
    lookup <- lookup[!duplicated(municipio)]
    lookup[, source_priority := NULL]
  }

  ambiguous <- lookup[, .N, by = municipio][N > 1]
  if (nrow(ambiguous) > 0) {
    stop("Ambiguous municipality-to-IBGE lookup entries remain after RJ filtering.")
  }

  lookup
}


# =========================================================
# Adjacency helpers
# =========================================================
read_model_adjacency <- function(model_ibge_codes) {
  adj <- as.data.table(read_parquet(ADJACENCY_FILE))

  if (!("co_muni_ori" %in% names(adj))) {
    stop("adjacency_matrix_correct.parquet must contain co_muni_ori.")
  }

  adj[, co_muni_ori := as.integer(co_muni_ori)]
  model_ibge_codes <- as.integer(model_ibge_codes)
  model_code_chr <- as.character(model_ibge_codes)

  missing_rows <- setdiff(model_ibge_codes, adj$co_muni_ori)
  missing_cols <- setdiff(model_code_chr, names(adj))
  if (length(missing_rows) > 0 || length(missing_cols) > 0) {
    stop("Adjacency matrix is missing one or more model municipalities.")
  }

  adj <- adj[match(model_ibge_codes, co_muni_ori)]
  adj_mat <- as.matrix(adj[, ..model_code_chr])
  storage.mode(adj_mat) <- "numeric"
  adj_mat[is.na(adj_mat)] <- 0
  adj_mat[adj_mat != 0] <- 1
  diag(adj_mat) <- 0

  # Keep adjacency symmetric for both BYM2 and the neighbor spatial lag.
  adj_mat <- ((adj_mat + t(adj_mat)) > 0) * 1
  rownames(adj_mat) <- model_code_chr
  colnames(adj_mat) <- model_code_chr

  edge_idx <- which(adj_mat > 0, arr.ind = TRUE)
  edges <- data.table(
    ibge_code = model_ibge_codes[edge_idx[, 1]],
    neighbor_ibge_code = model_ibge_codes[edge_idx[, 2]]
  )

  list(
    matrix = adj_mat,
    edges = edges,
    degree = rowSums(adj_mat)
  )
}


# =========================================================
# Distance-weighted spatial lag helpers
# =========================================================
haversine_distance_km <- function(lat1, lon1, lat2, lon2) {
  radius_km <- 6371.0088
  to_rad <- pi / 180
  phi1 <- lat1 * to_rad
  phi2 <- lat2 * to_rad
  delta_phi <- (lat2 - lat1) * to_rad
  delta_lambda <- (lon2 - lon1) * to_rad

  a <- sin(delta_phi / 2)^2 +
    cos(phi1) * cos(phi2) * sin(delta_lambda / 2)^2
  radius_km * 2 * atan2(sqrt(a), sqrt(1 - a))
}

manual_coordinate_fallback <- function() {
  # These five municipalities are absent from fluvi_road_ibge.parquet.
  # Coordinates are city-level latitude/longitude used only to complete S3 distances.
  data.table(
    ibge_code = c(3300456L, 3302700L, 3302858L, 3304144L, 3305109L),
    municipio = c(
      "belford roxo",
      "marica",
      "mesquita",
      "queimados",
      "sao joao de meriti"
    ),
    lat = c(-22.7589, -22.9651, -22.7674, -22.7200, -22.8058),
    lon = c(-43.4119, -42.9783, -43.4218, -43.5600, -43.3729),
    coord_source = "manual_city_coordinate_fallback"
  )
}

read_municipality_coordinates <- function(model_ibge_codes) {
  conn <- as.data.table(read_parquet(CONNECTIVITY_FILE))

  required_cols <- c(
    "co_muni_ori",
    "co_muni_des",
    "long_muni_ori",
    "lat_muni_ori",
    "long_muni_des",
    "lat_muni_des"
  )
  missing_cols <- setdiff(required_cols, names(conn))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "fluvi_road_ibge.parquet is missing coordinate columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  coords <- rbindlist(
    list(
      conn[, .(
        ibge_code = as.integer(co_muni_ori),
        lat = as.numeric(lat_muni_ori),
        lon = as.numeric(long_muni_ori),
        coord_source = "fluvi_road_origin"
      )],
      conn[, .(
        ibge_code = as.integer(co_muni_des),
        lat = as.numeric(lat_muni_des),
        lon = as.numeric(long_muni_des),
        coord_source = "fluvi_road_destination"
      )]
    ),
    use.names = TRUE
  )
  coords <- coords[!is.na(ibge_code) & !is.na(lat) & !is.na(lon)]
  coords <- coords[
    ,
    .(
      lat = median(lat),
      lon = median(lon),
      coord_source = "fluvi_road_median_coordinate"
    ),
    by = ibge_code
  ]

  model_ibge_codes <- as.integer(model_ibge_codes)
  fallback <- manual_coordinate_fallback()
  missing_before_fallback <- setdiff(model_ibge_codes, coords$ibge_code)
  if (length(missing_before_fallback) > 0) {
    fallback <- fallback[ibge_code %in% missing_before_fallback]
    coords <- rbindlist(list(coords, fallback[, .(ibge_code, lat, lon, coord_source)]), use.names = TRUE)
  }

  missing_after_fallback <- setdiff(model_ibge_codes, coords$ibge_code)
  if (length(missing_after_fallback) > 0) {
    stop(sprintf(
      "Missing coordinates for model IBGE codes: %s",
      paste(missing_after_fallback, collapse = ", ")
    ))
  }

  coords <- coords[ibge_code %in% model_ibge_codes]
  coords <- coords[match(model_ibge_codes, ibge_code)]
  if (any(is.na(coords$lat)) || any(is.na(coords$lon))) {
    stop("Coordinate lookup failed for one or more model municipalities.")
  }

  coords
}

build_distance_weights <- function(model_ibge_codes) {
  coords <- read_municipality_coordinates(model_ibge_codes)

  pairs <- CJ(
    ibge_code = as.integer(model_ibge_codes),
    source_ibge_code = as.integer(model_ibge_codes)
  )
  pairs <- pairs[ibge_code != source_ibge_code]

  origin <- coords[, .(ibge_code, origin_lat = lat, origin_lon = lon)]
  source <- coords[, .(
    source_ibge_code = ibge_code,
    source_lat = lat,
    source_lon = lon
  )]
  pairs <- merge(pairs, origin, by = "ibge_code", all.x = TRUE)
  pairs <- merge(pairs, source, by = "source_ibge_code", all.x = TRUE)

  pairs[, distance_km := haversine_distance_km(
    origin_lat,
    origin_lon,
    source_lat,
    source_lon
  )]
  pairs <- pairs[distance_km > 0 & distance_km <= DISTANCE_CUTOFF_KM]
  pairs[, raw_weight := exp(-distance_km / DISTANCE_DECAY_KM)]
  pairs[, weight := raw_weight / sum(raw_weight), by = ibge_code]

  if (any(is.na(pairs$weight)) || any(!is.finite(pairs$weight))) {
    stop("Distance-weight construction produced invalid weights.")
  }

  list(
    weights = pairs[, .(ibge_code, source_ibge_code, distance_km, weight)],
    coords = coords
  )
}


# =========================================================
# M5 + distance-weighted spatial-lag feature construction
# =========================================================
build_model_dataframe <- function() {
  df <- fread(COMBINED_FILE)
  df <- clean_columns(df)

  required <- c("municipio", "year", "week", "cases", "idhm", LAGGED_WEATHER_COVARIATES)
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  df[, municipio := normalize_name(municipio)]
  numeric_cols <- c("year", "week", "cases", "idhm", LAGGED_WEATHER_COVARIATES)
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

  lookup <- build_municipio_lookup()
  df <- merge(df, lookup, by = "municipio", all.x = TRUE)
  missing_ibge <- sort(unique(df[is.na(ibge_code), municipio]))
  if (length(missing_ibge) > 0) {
    stop(sprintf(
      "Municipios missing IBGE lookup: %s",
      paste(missing_ibge, collapse = ", ")
    ))
  }
  df[, ibge_code := as.integer(ibge_code)]

  setorder(df, municipio, date)

  weather_lookup <- df[, c("municipio", "date", LAGGED_WEATHER_COVARIATES), with = FALSE]
  setnames(
    weather_lookup,
    LAGGED_WEATHER_COVARIATES,
    paste0(LAGGED_WEATHER_COVARIATES, "_lag")
  )
  weather_lookup[, date := date + WEATHER_LAG_WEEKS * 7L]
  weather_lookup[, weather_lag_source_date := date - WEATHER_LAG_WEEKS * 7L]
  df <- merge(df, weather_lookup, by = c("municipio", "date"), all.x = TRUE)

  expected_weather_source_date <- df$date - WEATHER_LAG_WEEKS * 7L
  bad_weather_lag <- df[
    !is.na(weather_lag_source_date) &
      (weather_lag_source_date != expected_weather_source_date |
         weather_lag_source_date >= date)
  ]
  if (nrow(bad_weather_lag) > 0) {
    stop("Weather lag leakage check failed.")
  }

  case_lookup <- df[, .(municipio, date, cases_lag = cases)]
  case_lookup[, date := date + CASE_LAG_WEEKS * 7L]
  case_lookup[, cases_lag_source_date := date - CASE_LAG_WEEKS * 7L]
  df <- merge(df, case_lookup, by = c("municipio", "date"), all.x = TRUE)

  expected_case_source_date <- df$date - CASE_LAG_WEEKS * 7L
  bad_case_lag <- df[
    !is.na(cases_lag) &
      (cases_lag_source_date != expected_case_source_date |
         cases_lag_source_date >= date)
  ]
  if (nrow(bad_case_lag) > 0) {
    stop("Case lag leakage check failed.")
  }

  df[, log_cases_lag := log1p(cases_lag)]

  df[, row_id := .I]
  distance_info <- build_distance_weights(sort(unique(df$ibge_code)))
  distance_weights <- distance_info$weights

  distance_targets <- df[, .(
    row_id,
    ibge_code,
    date,
    distance_lag_source_date = date - CASE_LAG_WEEKS * 7L
  )]
  distance_candidates <- merge(
    distance_targets,
    distance_weights,
    by = "ibge_code",
    all.x = TRUE,
    allow.cartesian = TRUE
  )

  distance_source <- df[, .(
    source_ibge_code = ibge_code,
    distance_lag_source_date = date,
    source_cases = cases
  )]
  distance_candidates <- merge(
    distance_candidates,
    distance_source,
    by = c("source_ibge_code", "distance_lag_source_date"),
    all.x = TRUE
  )

  distance_summary <- distance_candidates[
    !is.na(source_cases),
    .(
      distance_cases_lag = sum(weight * source_cases) / sum(weight),
      distance_lag_n = .N,
      distance_weight_sum = sum(weight)
    ),
    by = row_id
  ]
  df <- merge(df, distance_summary, by = "row_id", all.x = TRUE, sort = FALSE)
  df[is.na(distance_cases_lag), distance_cases_lag := 0]
  df[is.na(distance_lag_n), distance_lag_n := 0L]
  df[is.na(distance_weight_sum), distance_weight_sum := 0]
  df[, distance_lag_source_date := date - CASE_LAG_WEEKS * 7L]
  df[, distance_log_cases_lag := log1p(distance_cases_lag)]

  bad_distance_lag <- df[
    distance_lag_n > 0 &
      (distance_lag_source_date >= date)
  ]
  if (nrow(bad_distance_lag) > 0) {
    stop("Distance-weighted case lag leakage check failed.")
  }

  df[, cases := pmax(as.integer(round(cases)), 0L)]

  rows_before_drop <- nrow(df)
  df <- df[complete.cases(df[, ..BASE_COVARIATES])]

  cat("M5/S3 dataframe built\n")
  cat("Original rows:", original_rows, "\n")
  cat("Model rows:", nrow(df), "\n")
  cat("Dropped missing S3 covariates:", rows_before_drop - nrow(df), "\n")
  cat("Municipios:", uniqueN(df$municipio), "\n")
  cat("Years:", paste(sort(unique(df$year)), collapse = ", "), "\n")
  cat("Distance decay km:", DISTANCE_DECAY_KM, "\n")
  cat("Distance cutoff km:", DISTANCE_CUTOFF_KM, "\n")
  cat("Rows with at least one distance-weighted lag value:", sum(df$distance_lag_n > 0), "\n")
  cat("Manual coordinate fallbacks used:", sum(distance_info$coords$coord_source == "manual_city_coordinate_fallback"), "\n")

  df[, row_id := NULL]
  setorder(df, municipio, date)
  df
}


# =========================================================
# INLA graph construction
# =========================================================
write_inla_graph <- function(model_ibge_codes, graph_file = GRAPH_FILE) {
  adj <- as.data.table(read_parquet(ADJACENCY_FILE))

  if (!("co_muni_ori" %in% names(adj))) {
    stop("adjacency_matrix_correct.parquet must contain co_muni_ori.")
  }

  adj[, co_muni_ori := as.integer(co_muni_ori)]
  model_ibge_codes <- as.integer(model_ibge_codes)
  model_code_chr <- as.character(model_ibge_codes)

  missing_rows <- setdiff(model_ibge_codes, adj$co_muni_ori)
  missing_cols <- setdiff(model_code_chr, names(adj))
  if (length(missing_rows) > 0 || length(missing_cols) > 0) {
    stop("Adjacency matrix is missing one or more model municipalities.")
  }

  adj <- adj[match(model_ibge_codes, co_muni_ori)]
  adj_mat <- as.matrix(adj[, ..model_code_chr])
  storage.mode(adj_mat) <- "numeric"
  adj_mat[is.na(adj_mat)] <- 0
  adj_mat[adj_mat != 0] <- 1
  diag(adj_mat) <- 0

  # Keep the graph symmetric for the BYM2 spatial effect.
  adj_mat <- ((adj_mat + t(adj_mat)) > 0) * 1

  degree <- rowSums(adj_mat)
  connected_mask <- degree > 0
  connected_codes <- model_ibge_codes[connected_mask]
  connected_mat <- adj_mat[connected_mask, connected_mask, drop = FALSE]

  n <- nrow(connected_mat)
  graph_obj <- inla.matrix2graph(Matrix(connected_mat, sparse = TRUE))
  inla.write.graph(graph_obj, filename = graph_file, mode = "ascii")
  graph <- inla.read.graph(graph_file)

  cat("INLA graph written:", graph_file, "\n")
  cat("Connected graph nodes:", n, "\n")
  cat("Connected graph edges:", sum(connected_mat) / 2, "\n")
  cat("Isolated municipalities:", sum(!connected_mask), "\n")

  list(
    graph = graph,
    graph_file = graph_file,
    lookup = data.table(
      ibge_code = model_ibge_codes,
      degree = degree,
      is_isolated = !connected_mask,
      spatial_idx = ifelse(
        connected_mask,
        match(model_ibge_codes, connected_codes),
        NA_integer_
      )
    )
  )
}


# =========================================================
# Model fitting
# =========================================================
build_inla_formula <- function(graph_file) {
  z_covariates <- paste0(BASE_COVARIATES, "_z")
  formula <- as.formula(paste(
    "cases ~ 1 +",
    paste(z_covariates, collapse = " + "),
    "+ f(week_idx, model = 'iid')",
    "+ f(year_idx, model = 'iid')",
    "+ f(spatial_idx, model = 'bym2', graph = graph_file, scale.model = TRUE)",
    "+ f(isolated_idx, model = 'iid')"
  ))
  environment(formula) <- environment()
  formula
}

fit_spatial_inla <- function(train_dt, graph_file) {
  formula <- build_inla_formula(graph_file)

  inla(
    formula,
    family = "nbinomial",
    data = train_dt,
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    num.threads = INLA_NUM_THREADS,
    verbose = FALSE
  )
}

predict_inla_mean <- function(fit, new_dt) {
  z_covariates <- paste0(BASE_COVARIATES, "_z")
  pred_dt <- copy(new_dt)
  pred_dt[, cases := NA_integer_]

  # INLA prediction is done by appending rows with missing response.
  fit_dt <- rbindlist(
    list(fit$.model_data, pred_dt),
    use.names = TRUE,
    fill = TRUE
  )

  formula <- fit$.formula
  result <- inla(
    formula,
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


# =========================================================
# Main
# =========================================================
main <- function() {
  df <- build_model_dataframe()

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

  full_dt <- standardize_full(copy(df), BASE_COVARIATES)
  full_fit <- fit_spatial_inla(full_dt, GRAPH_FILE)
  full_fit$.model_data <- full_dt
  full_fit$.formula <- build_inla_formula(GRAPH_FILE)

  cat("\nS3 full-data fixed effects:\n")
  print(full_fit$summary.fixed)

  cat("\nS3 full-data model criteria:\n")
  print(data.table(
    dic = full_fit$dic$dic,
    waic = full_fit$waic$waic
  ))

  if (RUN_TRAIN_TEST_EVALUATION) {
    train_dt <- df[year >= TRAIN_START_YEAR & year <= TRAIN_END_YEAR]
    test_dt <- df[year >= TEST_START_YEAR & year <= TEST_END_YEAR]

    test_start_date <- min(test_dt$date)
    test_rows_before_lag_filter <- nrow(test_dt)
    test_dt <- test_dt[
      cases_lag_source_date < test_start_date &
        distance_lag_source_date < test_start_date
    ]
    dropped_test_lag_rows <- test_rows_before_lag_filter - nrow(test_dt)

    scaled <- standardize_train_test(copy(train_dt), copy(test_dt), BASE_COVARIATES)
    train_dt <- scaled$train
    test_dt <- scaled$test

    train_fit <- fit_spatial_inla(train_dt, GRAPH_FILE)
    train_fit$.model_data <- train_dt
    train_fit$.formula <- build_inla_formula(GRAPH_FILE)

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
    cat("Dropped test rows with own/distance-weighted lagged cases inside test period:", dropped_test_lag_rows, "\n")
    cat("No-leakage policy: scaler fit on training data only.\n")
    cat("No-leakage policy: test own-case and distance-weighted case lags must come from before test start.\n")

    cat("\nS3 train/test metrics:\n")
    print(metrics)

    cat("\nS3 train-fit fixed effects:\n")
    print(train_fit$summary.fixed)

    if (SAVE_OUTPUTS) {
      fwrite(metrics, file.path(BASE_DIR, "spatial_inla_s3_train_test_metrics.csv"))
      output <- copy(test_dt[, .(municipio, year, week, date, cases)])
      output[, predicted_cases := test_pred]
      fwrite(output, file.path(BASE_DIR, "spatial_inla_s3_test_predictions.csv"))
    }
  }
}

if (Sys.getenv("INLA_RUN_MODEL", "1") == "1") {
  main()
}
