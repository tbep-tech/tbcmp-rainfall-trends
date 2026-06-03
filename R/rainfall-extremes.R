###############################################################################
# Tampa Bay Coastal Master Plan — Rainfall Extremes Trend Analysis
# Replicates key analyses from:
#   Wright, D.B., Bosma, C.D., & Lopez-Cantu, T. (2019).
#   "U.S. Hydrologic Design Standards Insufficient Due to Large Increases in
#    Frequency of Rainfall Extremes." Geophysical Research Letters, 46, 8144-8153.
#   https://doi.org/10.1029/2019GL083235
#
# Scope: All GHCN-Daily stations within the Tampa Bay Coastal Master Plan project
# area. Study period: 1920–2025 (extending the paper's primary analysis window)
#
# Analyses replicated:
#   1. GHCN station retrieval and spatial filter to TBCMP county areas
#   2. Exceedance counting relative to NOAA Atlas 14 IDF estimates (Section 3.1)
#   3. Negative binomial regression trend analysis (Section 3.1)
#   4. Rainstorm cluster identification (Section 3.2 / Section 2.4)
#   5. Design-vs-observed ARI comparison (Section 3.3 / Section 2.5)
#   6. Publication-quality figures for all analyses
#
# Dependencies (install once before sourcing):
#   install.packages(c(
#     "tidyverse", "sf", "terra", "rnoaa", "MASS", "pscl",
#     "zoo", "lubridate", "ggplot2", "patchwork", "scales",
#     "httr", "jsonlite", "tigris", "units"
#   ))
#
# Data sources fetched automatically at runtime:
#   • GHCN-Daily via rnoaa::ghcnd_*()
#   • NOAA Atlas 14 IDF values via NOAA Precipitation Frequency Data Server
#     (PFDS) REST API  <https://hdsc.nws.noaa.gov/pfds/>
#   • Tampa Bay watershed boundary via the USGS StreamStats or NHD
#     (a local GeoJSON fallback is provided if API is unavailable)
###############################################################################


# ── 0. Setup ─────────────────────────────────────────────────────────────────

library(tidyverse)
library(sf)
library(terra)
library(rnoaa)        # GHCN-Daily access
library(readnoaa)     # GHCN-Daily access, updated package
library(MASS)         # glm.nb (negative binomial)
library(pscl)         # zero-inflated models if needed
library(zoo)          # rollmean for moving-average smoother
library(lubridate)
library(ggplot2)
library(patchwork)
library(scales)
library(httr)
library(jsonlite)
library(tigris)       # Florida county shapefiles
library(units)

options(tigris_use_cache = TRUE)

# Set your NOAA API token (free at https://www.ncdc.noaa.gov/cdo-web/token)
# Replace the string below or set the environment variable before sourcing.
# Sys.setenv(NOAA_KEY = "YOUR_TOKEN_HERE")
noaa_key <- Sys.getenv("noaakey")
if (nchar(noaa_key) == 0) {
  message(
    "⚠  NOAA_KEY environment variable not set.\n",
    "   Register at https://www.ncdc.noaa.gov/cdo-web/token and run:\n",
    "   Sys.setenv(NOAA_KEY = 'your_token')\n",
    "   before sourcing this script."
  )
}

set.seed(2026)  # reproducibility

STUDY_START <- 1900L
STUDY_END   <- 2025L
MAX_MISSING_YEARS <- 25L          # stations with > 10 missing years excluded
MISSING_DAY_THRESHOLD <- 90L     # days missing → year treated as missing
CLUSTER_DIST_KM <- 25           # km radius for storm clustering
CLUSTER_DAY_WIN <- 2             # ±days for storm clustering


# ── 1. Tampa Bay Coastal Master Plan Boundary ──────────────────────────────────────────
# The Tampa Bay Coastal Master Plan project area spans Hillsborough, Pinellas,
# Pasco, Manatee, Sarasota, Citrus and Hernando counties.
# We use HUC-6 polygons from the USGS NHD as the authoritative boundary.
# A county-union fallback is provided in case the NHD service is unavailable.

fetch_watershed_boundary <- function() {
  message("Fetching Tampa Bay watershed (HUC-6 031002) from USGS NHD …")

  # USGS National Map / WBD GeoJSON endpoint — Layer 3 = 6-digit HU (Basin)
  url <- paste0(
    "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/3/query",
    "?where=huc6%3D%27031002%27",
    "&outFields=*&f=geojson"
  )

  resp <- tryCatch(GET(url, timeout(30)), error = function(e) NULL)

  if (!is.null(resp) && status_code(resp) == 200) {
    gj  <- content(resp, as = "text", encoding = "UTF-8")
    shp <- st_read(gj, quiet = TRUE)
    message("  \u2713 HUC-6 boundary loaded from USGS NHD.")
    return(st_transform(shp, 4326))
  }

  # ── Fallback: union of core Tampa Bay counties ────────────────────────────
  message("  NHD unavailable — falling back to county union.")
  fl_counties <- tigris::counties(state = "FL", cb = TRUE, year = 2020,
                                  progress_bar = FALSE) |>
    st_transform(4326)

  tb_counties <- c("Hillsborough", "Pinellas", "Pasco",
                   "Manatee", "Sarasota", "Citrus", "Hernando")

  watershed <- fl_counties |>
    filter(NAME %in% tb_counties) |>
    st_union() |>
    st_as_sf() |>
    mutate(name = "Tampa Bay Coastal Master Plan (county approx.)")

  return(watershed)

  saveRDS(watershed, "./data-raw/tbcmp_counties.rds")
}

usgs_watershed_sf <- fetch_watershed_boundary()
watershed_sf <- tbcmp_counties |>
                st_transform(4326)

# Bounding box for station search (add 0.1° buffer)
bb <- st_bbox(watershed_sf)
buffer <- 0.25
bbx <- c(xmin = bb["xmin"] - buffer,
         ymin = bb["ymin"] - buffer,
         xmax = bb["xmax"] + buffer,
         ymax = bb["ymax"] + buffer)

# ── 2. Fetch GHCN-Daily Station Inventory ────────────────────────────────────

message("Searching GHCN-Daily stations in Tampa Bay region …")

stations <- ghcnd_stations()

ghcn_stations_raw <- stations |>
  filter(
    latitude  >= bb["ymin"] - 0.25,
    latitude  <= bb["ymax"] + 0.25,
    longitude >= bb["xmin"] - 0.25,
    longitude <= bb["xmax"] + 0.25,
    element == "PRCP",
    last_year  >= STUDY_END - 10,
    first_year <= STUDY_START + 50   # allow stations starting a bit late
  )

#ghcn_stations_raw <- noaa_stations(bbox = c(27.10299, -83.10068, 28.65894, -81.63309),
#                                   limit = 1000L,
#                                   cache = FALSE)

# Spatial filter: keep only stations inside the watershed polygon
stations_sf <- st_as_sf(ghcn_stations_raw,
                        coords = c("longitude", "latitude"),
                        crs = 4326)

#stations_in_ws <- stations_sf[st_within(stations_sf, watershed_sf,
#                                        sparse = FALSE)[, 1], ]

stations_in_ws <- stations_sf #Keep all from the BB

list_nws_in_ws <- as.list(stations_in_ws$id)

message(sprintf("  Found %d candidate stations inside watershed boundary.",
                nrow(stations_in_ws)))


# ── 3. Download Daily Precipitation Records ───────────────────────────────────
# readnoaa caches to ~/.readnoaa by default; repeated runs reuse local copies.

message("Downloading GHCN-Daily PRCP records (this may take several minutes) …")

download_station_prcp <- function(station_id) {
  tryCatch({
    dat <- noaa_daily(station = station_id,
                      start_date = paste0(STUDY_START, "-01-01"),
                      end_date =   paste0(STUDY_END,   "-12-31"),
                      datatypes = c("PRCP"),
                      units = "metric",
                      cache = TRUE) |>
            dplyr::select(date, prcp) |>
      mutate(
        station_id = station_id,
        year       = year(date)
      )
    return(dat)
  }, error = function(e) {
    message("  ⚠  Could not download ", station_id, ": ", conditionMessage(e))
    return(NULL)
  })
}

all_prcp <- map(stations_in_ws$id, download_station_prcp) |>
  compact() |>
  bind_rows()

message(sprintf("  Downloaded records for %d stations.",
                n_distinct(all_prcp$station_id)))


# ── 4. Quality Control: Remove Stations with Excessive Missing Data ───────────

assess_completeness <- function(dat) {
  dat |>
    group_by(station_id, year) |>
    summarise(n_obs   = sum(!is.na(prcp)),
              missing = 365 - n_obs,           # approximate
              .groups = "drop") |>
    mutate(year_missing = missing > MISSING_DAY_THRESHOLD) |>
    group_by(station_id) |>
    summarise(n_missing_years = sum(year_missing), .groups = "drop") |>
    filter(n_missing_years <= MAX_MISSING_YEARS)
}

good_stations <- assess_completeness(all_prcp)
prcp_qc <- all_prcp |> semi_join(good_stations, by = "station_id")

# Keep station metadata for good stations
stations_final <- stations_in_ws |>
  bind_cols(
    st_coordinates(stations_in_ws) |>
      as.data.frame() |>
      rename(longitude = X, latitude = Y)) |>
  filter(id %in% good_stations$station_id) |>
  dplyr::select(station_id = id, name, latitude, longitude) |>
  st_drop_geometry()


message(sprintf(
  "  %d stations pass completeness filter (<= %d missing years).",
  nrow(stations_final), MAX_MISSING_YEARS
))


# ── 5. Compute Multi-Day Running Totals (24-hr and 7-day) ─────────────────────
# Paper uses 24-hr and up to 7-day durations. We compute 1-day (calendar)
# and 7-day running sums. Currently, this workflow neglects the potential for
# multiple exceedances in a given year when comparing to ATlas 14 IDF values.

compute_running_totals <- function(dat, durations_days = c(1, 2, 3, 4, 7)) {
  dat |>
    arrange(station_id, date) |>
    group_by(station_id) |>
    mutate(across(
      prcp,
      list(
        d1 = ~ .x,
        d2 = ~ zoo::rollsum(.x, 2, fill = NA, align = "right"),
        d3 = ~ zoo::rollsum(.x, 3, fill = NA, align = "right"),
        d4 = ~ zoo::rollsum(.x, 4, fill = NA, align = "right"),
        d7 = ~ zoo::rollsum(.x, 7, fill = NA, align = "right")
      ),
      .names = "{.fn}"
    )) |>
    ungroup() |>
    # Non-overlapping: keep only one value per duration window per year
    # (simplified: take annual maximum for each duration)
    group_by(station_id, year) |>
    summarise(
      ann_max_1d = max(d1, na.rm = TRUE),
      ann_max_2d = max(d2, na.rm = TRUE),
      ann_max_3d = max(d3, na.rm = TRUE),
      ann_max_4d = max(d4, na.rm = TRUE),
      ann_max_7d = max(d7, na.rm = TRUE),
      .groups = "drop"
    )
}

annual_maxima <- compute_running_totals(prcp_qc)

message("  Annual maxima computed for all duration windows.")

# ── 6. Fetch NOAA Atlas 14 IDF Thresholds ─────────────────────────────────────
#
# Correct endpoint: NOAA Precipitation Frequency Data Server (PFDS) "fe_text"
# CSV service.  Each call returns a plain-text CSV for a single lat/lon point
# covering all standard durations (5-min through 60-day) and all standard ARIs
# (1-yr through 1000-yr).  Values are returned in INCHES; we convert to mm.
#
# Documented URL pattern (CONUS):
#   https://hdsc.nws.noaa.gov/cgi-bin/hdsc/new/fe_text_mean.csv
#     ?lat=<lat>&lon=<lon>&type=pf&data=depth&units=us&series=pds
#
# Response format (first ~30 rows are metadata comment lines beginning with
# "by:" or blank; the data block starts with a header row):
#
#   "","1-yr","2-yr","5-yr","10-yr","25-yr","50-yr","100-yr","200-yr",
#     "500-yr","1000-yr"
#   "5-min",0.33,0.40, ...
#   "10-min", ...
#   ...
#   "24-hr",3.10,3.68, ...
#   "48-hr",3.72,4.42, ...
#   ...
#
#
# We query once per UNIQUE (lat, lon) pair rounded to 4 decimal places,
# cache results in a named list, and apply retries with exponential back-off
# to handle transient NOAA server errors gracefully.
#
# Column naming convention used downstream:
#   dur_24hr_ari2, dur_24hr_ari10, dur_24hr_ari100, dur_24hr_ari500
#   dur_48hr_ari2, dur_48hr_ari10, dur_48hr_ari100, dur_48hr_ari500
# These names are constructed and consumed consistently inside this section
# and in identify_exceedances() below.

# Map from the text labels in the PFDS CSV to the short duration tags we use
PFDS_DUR_MAP <- c(
  "24-hr:" = "24hr",
  "2-day:" = "48hr"
)

# Return-period column headers as they appear in the PFDS CSV response
PFDS_ARI_LABELS <- c(
  "1"    =    1L,
  "2"    =    2L,
  "5"    =    5L,
  "10"   =   10L,
  "25"   =   25L,
  "50"   =   50L,
  "100"  =  100L,
  "200"  =  200L,
  "500"  =  500L,
  "1000" = 1000L
)

#' Fetch Atlas 14 IDF depth estimates for a single location.
#'
#' @param lat  Numeric latitude  (decimal degrees, WGS84)
#' @param lon  Numeric longitude (decimal degrees, WGS84; negative = West)
#' @param max_tries  Integer number of HTTP attempts before giving up
#' @return A tidy data frame with columns: duration_tag, ari_yr, idf_mm
#'         or NULL on complete failure.

fetch_atlas14_idf <- function(lat, lon, max_tries = 4L) {

  # NOAA PFDS "fe_text_mean" endpoint — returns plain-text CSV, CONUS only
  base_url <- "https://hdsc.nws.noaa.gov/cgi-bin/hdsc/new/fe_text_mean.csv"

  params <- list(
    lat    = sprintf("%.4f", lat),
    lon    = sprintf("%.4f", lon),  # lon is already negative for West
    type   = "pf",
    data   = "depth",
    units  = "us",      # return values in INCHES (only valid option for PFDS)
    series = "pds"      # partial-duration series  (matches Atlas 14 convention)
  )

  resp <- NULL
  for (attempt in seq_len(max_tries)) {
    resp <- tryCatch(
      httr::GET(base_url, query = params,
                httr::timeout(30),
                httr::user_agent("R/rainfall-extremes-study")),
      error = function(e) NULL
    )
    if (!is.null(resp) && httr::status_code(resp) == 200L) break
    wait_sec <- 2 ^ attempt          # exponential back-off: 2, 4, 8 sec
    message(sprintf(
      "    attempt %d failed (status %s) — retrying in %ds …",
      attempt,
      if (is.null(resp)) "connection error" else httr::status_code(resp),
      wait_sec
    ))
    Sys.sleep(wait_sec)
  }

  if (is.null(resp) || httr::status_code(resp) != 200L) return(NULL)

  raw_text <- httr::content(resp, as = "text", encoding = "UTF-8")

  # ── Parse the CSV response ──────────────────────────────────────────────────
  # The PFDS response contains metadata comment lines before the data table.
  # Comment lines start with "by:", contain colons, or are blank.
  # The actual CSV data block begins with a header row whose first token is
  # empty (row label column) followed by ARI labels like "1-yr","2-yr", …
  # We split on newlines, drop comment/blank lines, and read the remainder.

  lines <- strsplit(raw_text, "\n")[[1]]
  lines <- lines[nchar(trimws(lines)) > 0]

  # Identify the header row: it contains "1-yr" (first ARI column)
  header_idx <- which(grepl(", 1,2,5", lines, fixed = TRUE))

  if (length(header_idx) == 0) {
    message(sprintf("    ⚠  PFDS response for (%.4f, %.4f) has no recognisable header.",
                    lat, lon))
    return(NULL)
  }
  header_idx <- header_idx[1]

  # Read from header row onward as a CSV string
  data_block <- paste(lines[header_idx:length(lines)], collapse = "\n")
  data_block <- data_block[nchar(data_block) > 0]

  pfds_df <- tryCatch(
    read.csv(text = data_block,
             header      = TRUE,
             check.names = FALSE,
             stringsAsFactors = FALSE),
    error = function(e) {
      message(sprintf("    ⚠  CSV parse error for (%.4f, %.4f): %s",
                      lat, lon, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(pfds_df) || nrow(pfds_df) == 0) return(NULL)

  # First column contains duration labels (e.g. "5-min", "24-hr").
  # Rename it for clarity.
  names(pfds_df)[1] <- "duration_label"

  # Keep only durations needed for this study
  pfds_df <- pfds_df[pfds_df$duration_label %in% names(PFDS_DUR_MAP), ,
                     drop = FALSE]
  if (nrow(pfds_df) == 0) return(NULL)

  # ── Reshape to tidy long format ─────────────────────────────────────────────
  # Identify which ARI columns are present in this response
  ari_cols_present <- intersect(names(PFDS_ARI_LABELS), names(pfds_df))
  if (length(ari_cols_present) == 0) return(NULL)

  tidy_df <- pfds_df |>
    dplyr::select(duration_label, dplyr::all_of(ari_cols_present)) |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(ari_cols_present),
      names_to  = "ari_label",
      values_to = "idf_in"           # values are in INCHES at this point
    ) |>
    dplyr::mutate(
      duration_tag = PFDS_DUR_MAP[duration_label],   # "24hr" or "48hr"
      ari_yr       = PFDS_ARI_LABELS[ari_label],     # integer year
      idf_mm       = as.numeric(idf_in) * 25.4       # inches → millimetres
    ) |>
    dplyr::filter(!is.na(idf_mm), !is.na(ari_yr)) |>
    dplyr::select(duration_tag, ari_yr, idf_mm)

  return(tidy_df)
}


# ── Fetch IDF for every QC-passing station, with caching ─────────────────────
# We round coordinates to 3 decimal places (~100 m) before keying the cache,
# since PFDS grid resolution is ~2.5 arcmin.  This avoids redundant calls for
# co-located or near-co-located stations.

message("Fetching Atlas 14 IDF values for each station via NOAA PFDS …")
message("  (Each call may take 1–5 s; a progress counter is shown every 5 stations.)")

coord_cache <- list()   # key = "lat_lon" → tidy IDF data frame

idf_list <- vector("list", nrow(stations_final))

for (i in seq_len(nrow(stations_final))) {
  sid <- stations_final$station_id[i]
  lat <- round(stations_final$latitude[i],  3)
  lon <- round(stations_final$longitude[i], 3)
  cache_key <- paste(lat, lon, sep = "_")

  if (!is.null(coord_cache[[cache_key]])) {
    # Re-use cached result for this coordinate
    idf_list[[i]] <- coord_cache[[cache_key]] |>
      dplyr::mutate(station_id = sid)
  } else {
    result <- fetch_atlas14_idf(lat, lon)
    coord_cache[[cache_key]] <- result       # cache even if NULL
    if (!is.null(result)) {
      idf_list[[i]] <- result |> dplyr::mutate(station_id = sid)
    }
    Sys.sleep(0.6)   # ~1.6 req/s — well within NOAA's tolerance
  }

  if (i %% 5 == 0 || i == nrow(stations_final)) {
    message(sprintf("  … %d / %d stations processed", i, nrow(stations_final)))
  }
}

idf_all <- purrr::compact(idf_list) |> dplyr::bind_rows() |>
           mutate(idf_mm = case_when(duration_tag == "24hr" ~ idf_mm/1.13,
                                     duration_tag == "48hr" ~ idf_mm/1.04))

if (nrow(idf_all) == 0) {
  stop(
    "No Atlas 14 IDF values were retrieved.\n",
    "Check your internet connection and that the NOAA PFDS service is reachable:\n",
    "  https://hdsc.nws.noaa.gov/cgi-bin/hdsc/new/fe_text_mean.csv",
    "    ?lat=27.9506&lon=-82.4572&type=pf&data=depth&units=us&series=pds"
  )
}

message(sprintf("  ✓ Atlas 14 IDF values retrieved for %d of %d stations.",
                dplyr::n_distinct(idf_all$station_id), nrow(stations_final)))

# ── Pivot to wide format for downstream exceedance lookup ─────────────────────
# Column names follow the pattern:  dur_<tag>_ari<N>
#   e.g. dur_24hr_ari100,  dur_48hr_ari2
# This scheme is explicit and avoids the name-collision risk of the original
# names_glue approach when duration tags contain non-alphanumeric characters.

idf_wide <- idf_all |>
  dplyr::filter(duration_tag %in% c("24hr", "48hr"),
                ari_yr       %in% c(1L, 2L, 5L, 10L, 25L, 50L, 100L, 200L)) |>
  dplyr::mutate(col_name = paste0("dur_", duration_tag, "_ari", ari_yr)) |>
  dplyr::select(station_id, col_name, idf_mm) |>
  tidyr::pivot_wider(
    id_cols     = station_id,
    names_from  = col_name,
    values_from = idf_mm
  )

# ── 7. Identify Atlas 14 Exceedances (§2.3 / §3.1) ───────────────────────────
# For each station-year, check whether the adjusted annual maximum
# exceeds the Atlas 14 threshold for each ARI of interest.

TARGET_ARIS <- c(1, 2, 5, 10, 25, 50, 100, 200)

identify_exceedances <- function(annual_max_df, idf_wide_df) {
  joined <- annual_max_df |>
    dplyr::left_join(idf_wide_df, by = "station_id")

  exceedance_rows <- purrr::map_dfr(TARGET_ARIS, function(ari) {
    # Column names follow the dur_24hr_ariN scheme established in section 6
    col_24h <- paste0("dur_24hr_ari", ari)

    if (!col_24h %in% names(joined)) {
      message(sprintf("    \u26a0  Column '%s' not found in idf_wide -- skipping ARI %d.",
                      col_24h, ari))
      return(dplyr::tibble())
    }

    joined |>
      dplyr::filter(!is.na(.data[[col_24h]])) |>
      dplyr::mutate(
        ari_design  = ari,
        duration    = "24hr",
        threshold   = .data[[col_24h]],    # mm, converted from inches in section 6
        exceeds     = ann_max_1d >= threshold  # Observed 1-day max rainfall in a given year
      ) |>
      dplyr::select(station_id, year, ari_design, duration, threshold,
                    obs_rainfall = ann_max_1d, exceeds)
  })

  return(exceedance_rows)
}

exceedances <- identify_exceedances(annual_maxima, idf_wide)
message(sprintf("  Exceedance table built: %d station-year rows.", nrow(exceedances)))


# ── 8. Aggregate Annual Exceedance Counts (CONUS-style pooling) ───────────────
# Sum exceedances across all stations for each year × ARI combination.
# This is the regional aggregation that boosts signal-to-noise (§2.3).

annual_counts <- exceedances |>
  filter(exceeds) |>
  group_by(year, ari_design, duration) |>
  summarise(n_exceedances = n(), .groups = "drop")

# Fill years with zero counts
full_grid <- expand.grid(
  year       = STUDY_START:STUDY_END,
  ari_design = TARGET_ARIS,
  duration   = "24hr"
)

annual_counts_full <- full_grid |>
  left_join(annual_counts, by = c("year", "ari_design", "duration")) |>
  replace_na(list(n_exceedances = 0L))


# ── 9. Negative Binomial Regression Trend Analysis (§2.3 / §3.1) ─────────────
# Paper uses negative binomial regression (Hilbe 2011) because exceedance
# counts are overdispersed.  We regress annual counts on calendar year.

fit_nb_trend <- function(counts_df, ari_val) {
  dat <- counts_df |>
    filter(ari_design == ari_val) |>
    mutate(year_centered = year - mean(year))  # center for numerical stability

  if (sum(dat$n_exceedances) == 0) {
    return(tibble(ari = ari_val, pct_change = NA_real_,
                  p_value = NA_real_, converged = FALSE))
  }

  mod <- tryCatch(
    MASS::glm.nb(n_exceedances ~ year_centered, data = dat),
    error = function(e) NULL
  )

  if (is.null(mod)) {
    return(tibble(ari = ari_val, pct_change = NA_real_,
                  p_value = NA_real_, converged = FALSE))
  }

  coef_yr  <- coef(mod)["year_centered"]
  se_yr    <- sqrt(vcov(mod)["year_centered", "year_centered"])
  p_val    <- 2 * pnorm(-abs(coef_yr / se_yr))

  # Percent change over the full study period (125 years)
  n_years   <- STUDY_END - STUDY_START
  pct_chg   <- (exp(coef_yr * n_years) - 1) * 100

  tibble(
    ari        = ari_val,
    coef_year  = coef_yr,
    pct_change = pct_chg,
    p_value    = p_val,
    converged  = mod$converged,
    model      = list(mod),
    data       = list(dat)
  )
}

nb_results <- map_dfr(TARGET_ARIS, ~ fit_nb_trend(annual_counts_full, .x))

message("\n── Negative Binomial Regression Results ──────────────────────────────")
nb_results |>
  dplyr::select(ari, pct_change, p_value, converged) |>
  mutate(
    pct_change = round(pct_change, 1),
    p_value    = signif(p_value, 3),
    sig        = case_when(
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.10  ~ ".",
      TRUE            ~ ""
    )
  ) |>
  print()


# ── 10. Storm Cluster Identification (§2.4 / §3.2) ───────────────────────────
# A storm "cluster" groups exceedances at multiple stations that occurred
# within CLUSTER_DIST_KM km and CLUSTER_DAY_WIN days of each other.
# This avoids conflating storm frequency with storm spatial extent.

identify_clusters <- function(exceedance_df, stations_df,
                              ari_target = 100,
                              dist_km    = CLUSTER_DIST_KM,
                              day_win    = CLUSTER_DAY_WIN) {
  # Work with individual daily exceedances (not annual maxima)
  # For cluster purposes, use the raw daily exceedance dates
  # We re-derive from full daily data
  exc_daily <- exceedance_df |>
    filter(ari_design == ari_target, exceeds) |>
    left_join(
      stations_df |> dplyr::select(station_id, latitude, longitude),
      by = "station_id"
    )

  if (nrow(exc_daily) == 0) {
    message(sprintf("  No %d-year exceedances found for clustering.", ari_target))
    return(tibble())
  }

  # Process each year independently
  cluster_results <- map_dfr(STUDY_START:STUDY_END, function(yr) {
    yr_exc <- exc_daily |> filter(year == yr)
    if (nrow(yr_exc) == 0) return(tibble(year = yr, n_clusters = 0L))

    # Distance matrix (Haversine via sf)
    pts <- st_as_sf(yr_exc,
                    coords = c("longitude", "latitude"), crs = 4326) |>
      st_transform(32617)  # UTM zone 17N (Florida)

    dist_m <- as.matrix(st_distance(pts))

    # Greedy recursive clustering
    assigned    <- rep(0L, nrow(yr_exc))
    cluster_id  <- 1L

    for (i in seq_len(nrow(yr_exc))) {
      if (assigned[i] != 0) next
      assigned[i] <- cluster_id
      # Find all un-assigned points within distance threshold
      within_dist <- which(dist_m[i, ] <= units::set_units(dist_km * 1000, "m") & assigned == 0)
      if (length(within_dist) > 0) assigned[within_dist] <- cluster_id
      cluster_id <- cluster_id + 1L
    }

    tibble(year = yr, n_clusters = max(assigned, na.rm = TRUE))
  })

  return(cluster_results)
}

message("Identifying storm clusters for 10-year exceedances …")
clusters_10yr <- identify_clusters(exceedances, stations_final, ari_target = 10)

message("Identifying storm clusters for 25-year exceedances …")
clusters_25yr <- identify_clusters(exceedances, stations_final, ari_target = 25)

message("Identifying storm clusters for 50-year exceedances …")
clusters_100yr  <- identify_clusters(exceedances, stations_final, ari_target = 100)

cluster_trends <- bind_rows(
  clusters_10yr |> mutate(ari = "10-year"),
  clusters_25yr |> mutate(ari = "25-year"),
  clusters_100yr  |> mutate(ari = "100-year")
) |>
  group_by(ari) |>
  mutate(year_centered = year - mean(year))

fit_cluster_trend <- function(cluster_df, ari_label) {
  dat <- cluster_df |> filter(ari == ari_label)
  mod <- tryCatch(
    MASS::glm.nb(n_clusters ~ year_centered, data = dat),
    error = function(e) NULL
  )
  if (is.null(mod)) return(NULL)
  pct_chg <- (exp(coef(mod)["year_centered"] * (STUDY_END - STUDY_START)) - 1) * 100
  p_val   <- summary(mod)$coefficients["year_centered", "Pr(>|z|)"]
  list(model    = mod,
       pct_chg  = pct_chg,
       p_val    = p_val,
       ari_label = ari_label)
}

trend_10 <- fit_cluster_trend(cluster_trends, "10-year")
trend_25 <- fit_cluster_trend(cluster_trends, "25-year")
trend_100  <- fit_cluster_trend(cluster_trends, "100-year")

message("\n── Cluster Trend Results ──────────────────────────────────────────────")
for (tr in list(trend_10, trend_25, trend_100)) {
  if (!is.null(tr)) {
    message(sprintf("  %s clusters: %+.0f%% change (p = %.3f)",
                    tr$ari_label, tr$pct_chg, tr$p_val))
  }
}


# ── 11. Design vs. Observed ARI (§2.5 / §3.3) ────────────────────────────────
# "Observed ARI" = 1 / (annual exceedance fraction across stations)
# When observed ARI < design ARI, infrastructure is underperforming.

n_stations_per_year <- prcp_qc |>
  group_by(year) |>
  summarise(n_active = n_distinct(station_id), .groups = "drop")

observed_ari_100 <- exceedances |>
  filter(ari_design == 100) |>
  group_by(year) |>
  summarise(n_exc = sum(exceeds), .groups = "drop") |>
  left_join(n_stations_per_year, by = "year") |>
  mutate(
    exc_fraction  = n_exc / n_active,
    observed_ari  = ifelse(exc_fraction > 0, ((1 / exc_fraction)/100), 0)
  )

# 5-year centered moving average smoother (7-yr used in paper's approach §2.5)
obs_ari_smooth <- observed_ari_100 |>
  arrange(year) |>
  mutate(
    ari_smooth = zoo::rollmean(observed_ari, k = 5,
                               fill = NA, align = "center")
  )


# ── 12. Figures ───────────────────────────────────────────────────────────────
# Figure A — Station map (analogue of Figure 1b)
# Figure B — Annual 50-year exceedance counts + NB regression (Figure 1c)
# Figure C — Trend summary by ARI (Figure 2)
# Figure D — Cluster trends (Figure 3c/3d)
# Figure E — Observed vs. design ARI time series (Figure 4a)

theme_wb <- function(...) {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92"),
      legend.position  = "bottom",
      ...
    )
}

# ── Fig A: Station map ────────────────────────────────────────────────────────
exceedance_totals <- exceedances |>
  filter(ari_design == 100, exceeds) |>
  count(station_id, name = "n_100yr_exc")

station_plot_df <- stations_final |>
  left_join(exceedance_totals, by = "station_id") |>
  mutate(
    n_100yr_exc = replace_na(n_100yr_exc, 0),
    exc_cat     = cut(n_100yr_exc,
                      breaks = c(-Inf, 0, 1, 2, Inf),
                      labels = c("0", "1", "2", "≥3"))
  )

fl_outline <- tigris::states(cb = TRUE, year = 2020, progress_bar = FALSE) |>
  filter(STUSPS == "FL") |>
  st_transform(4326)

fig_A <- ggplot() +
  geom_sf(data = fl_outline,       fill = "grey96", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = watershed_sf,     fill = "#cce5ff", colour = "#2171b5",
          linewidth = 0.7, alpha = 0.4) +
  geom_point(
    data = station_plot_df,
    aes(x = longitude, y = latitude, fill = exc_cat, size = exc_cat),
    shape = 21, colour = "black", stroke = 0.4
  ) +
  scale_fill_manual(
    values = c("0" = "grey70", "1" = "#fc8d59", "2" = "#d7301f", "≥3" = "#7f0000"),
    name   = "# of 1-day\n100-yr exceedances"
  ) +
  scale_size_manual(
    values = c("0" = 2, "1" = 3, "2" = 4, "≥3" = 5),
    name   = "# of 1-day\n100-yr exceedances"
  ) +
  coord_sf(xlim = c(bb["xmin"] - 0.5, bb["xmax"] + 0.5),
           ylim = c(bb["ymin"] - 0.5, bb["ymax"] + 0.5)) +
  labs(
    title    = "GHCN Stations — TBCMP Project Area",
    subtitle = sprintf("Total 100-year, 24-hr exceedances (1900–2025); N = %d stations",
                       nrow(stations_final)),
    x = "Longitude", y = "Latitude"
  ) +
  theme_wb()

# ── Fig B: Annual exceedance counts + regression ──────────────────────────────
nb_100 <- nb_results |> filter(ari == 100)

if (!is.null(nb_100$model[[1]])) {
  pred_df <- tibble(
    year          = STUDY_START:STUDY_END,
    year_centered = (STUDY_START:STUDY_END) - mean(STUDY_START:STUDY_END)
  ) |>
    (\(d) mutate(d,
      fit = predict(nb_100$model[[1]], newdata = d, type = "response"),
      se  = predict(nb_100$model[[1]], newdata = d, type = "response",
                          se.fit = TRUE)$se.fit,
      lwr = fit - 1.96 * se,
      upr = fit + 1.96 * se
    ))()
} else {
  pred_df <- tibble(year = integer(), fit = numeric(),
                    lwr = numeric(), upr = numeric())
}

count_100 <- annual_counts_full |> filter(ari_design == 100)
label_txt <- if (!is.na(nb_100$pct_change)) {
  sprintf("+%d%% increase, p = %.3f",
          round(nb_100$pct_change),
          nb_100$p_value)
} else "Insufficient data"

fig_B <- ggplot() +
  geom_ribbon(data = pred_df,
              aes(x = year, ymin = pmax(lwr, 0), ymax = upr),
              fill = "#2171b5", alpha = 0.25) +
  geom_line(data  = pred_df,
            aes(x = year, y = fit),
            colour = "#2171b5", linewidth = 1.1) +
  geom_line(data   = count_100,
            aes(x  = year, y = n_exceedances),
            colour = "black", linewidth = 0.5) +
  annotate("text", x = STUDY_START + 2, y = max(count_100$n_exceedances) * 0.95,
           label = label_txt, hjust = 0, size = 3.5, colour = "#2171b5") +
  scale_x_continuous(breaks = seq(1900, 2030, 10)) +
  labs(
    title    = "Trend in 100-year, 24-hr Rainfall Exceedances",
    subtitle = "TBCMP Project Area (all GHCN stations pooled)",
    x = "Year", y = "Number of 100-year exceedances"
  ) +
  theme_wb()

# ── Fig C: Trend summary heatmap by ARI ───────────────────────────────────────
trend_summary <- nb_results |>
  dplyr::select(ari, pct_change, p_value) |>
  mutate(
    pct_change_lab = ifelse(!is.na(pct_change),
                            sprintf("%+.0f%%", pct_change), "NA"),
    sig_lab = case_when(
      p_value < 0.01 ~ "p < 0.01",
      p_value < 0.05 ~ "p < 0.05",
      p_value < 0.10 ~ "p < 0.10",
      TRUE           ~ "p ≥ 0.10"
    ),
    ari_label = paste0(ari, "-year")
  )

fig_C <- ggplot(trend_summary,
                aes(x = "24h", y = fct_reorder(ari_label, ari))) +
  geom_tile(aes(fill = sig_lab), colour = "white", linewidth = 1.5) +
  geom_text(aes(label = pct_change_lab), size = 5, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "p < 0.01" = "#d73027",
      "p < 0.05" = "#f46d43",
      "p < 0.10" = "#fdae61",
      "p ≥ 0.10" = "#e0f3f8"
    ),
    name = "Significance"
  ) +
  labs(
    title    = "% Change in Exceedance Frequency (1900–2025)",
    subtitle = "TBCMP Project Area + | 24-hr duration",
    x        = "Duration",
    y        = "Design ARI"
  ) +
  theme_wb(axis.text.x = element_text(angle = 0))

# ── Fig D: Cluster trends ─────────────────────────────────────────────────────
pred_clusters <- function(tr_obj, cluster_df, ari_label) {
  if (is.null(tr_obj)) return(NULL)
  dat <- cluster_df |> filter(ari == ari_label)
  dat |> mutate(fit = predict(tr_obj$model, newdata = dat, type = "response"))
}
clust_pred_10 <- pred_clusters(trend_10, cluster_trends, "10-year")
clust_pred_25 <- pred_clusters(trend_25, cluster_trends, "25-year")
clust_pred_100  <- pred_clusters(trend_100,  cluster_trends, "100-year")

fig_D_data <- bind_rows(
  cluster_trends |> filter(ari == "10-year")  |> mutate(type = "Observed (10-yr)"),
  cluster_trends |> filter(ari == "25-year") |> mutate(type = "Observed (25-yr)"),
  cluster_trends |> filter(ari == "100-year") |> mutate(type = "Observed (100-yr)"),
)

label_100  <- if (!is.null(trend_100))
  sprintf("%+.0f%%, p = %.3f", trend_100$pct_chg,  trend_100$p_val) else ""
label_25 <- if (!is.null(trend_25))
  sprintf("%+.0f%%, p = %.3f", trend_25$pct_chg, trend_25$p_val) else ""
label_10 <- if (!is.null(trend_10))
  sprintf("%+.0f%%, p = %.3f", trend_10$pct_chg, trend_10$p_val) else ""

fig_D <- ggplot(fig_D_data, aes(x = year, y = n_clusters, colour = ari)) +
  geom_line(linewidth = 0.5, alpha = 0.7) +
  {if (!is.null(clust_pred_10))
    geom_line(data = clust_pred_10 |> mutate(ari = "10-year"),
              aes(y = fit), linewidth = 1.2, linetype = "dotted")} +
  {if (!is.null(clust_pred_25))
    geom_line(data = clust_pred_25 |> mutate(ari = "25-year"),
              aes(y = fit), linewidth = 1.2, linetype = "dashed")} +
  {if (!is.null(clust_pred_100))
    geom_line(data = clust_pred_100 |> mutate(ari = "100-year"),
              aes(y = fit), linewidth = 1.2)} +
  annotate("text", x = STUDY_END - 2, y = max(fig_D_data$n_clusters, na.rm = TRUE),
           label = paste("10-year:", label_10),
           hjust = 1, size = 3, colour = "grey") +
  annotate("text", x = STUDY_END - 2, y = max(fig_D_data$n_clusters, na.rm = TRUE) * 0.8,
           label = paste("25-year:", label_25),
           hjust = 1, size = 3, colour = "#5e3c99") +
  annotate("text", x = STUDY_END - 2, y = max(fig_D_data$n_clusters, na.rm = TRUE) * 0.6,
           label = paste("100-year:", label_100),
           hjust = 1, size = 3, colour = "#e66101") +
  scale_colour_manual(values = c("10-year" = "grey", "25-year" = "#5e3c99", "100-year" = "#e66101")) +
  scale_x_continuous(breaks = seq(1900, 2030, 10)) +
  labs(
    title    = "Rainstorm Cluster Trends — TBCMP Project Area",
    subtitle = "Clusters within 25 km / 1 day",
    x = "Year", y = "Number of clusters",
    colour = "Design ARI"
  ) +
  theme_wb()

# ── Fig E: Observed vs. design ARI ───────────────────────────────────────────
fig_E <- ggplot(obs_ari_smooth, aes(x = year)) +
  geom_hline(yintercept = 0.01, colour = "black", linewidth = 0.8,
             linetype = "dashed") +
  geom_line(aes(y = observed_ari), colour = "grey60", linewidth = 0.5, na.rm = TRUE) +
  geom_line(aes(y = ari_smooth),   colour = "#1a9641", linewidth = 1.3, na.rm = TRUE) +
  annotate("text", x = STUDY_START + 1, y = 105,
           label = "Design ARI = 100 years", hjust = 0, size = 3.2) +
  scale_y_continuous(
    name   = "Exceedance probability",
    limits = c(0, 0.1),
    sec.axis = sec_axis(~ . / 100, name = "Exceedance probability",
                        labels = scales::label_percent(accuracy = 0.01))
  ) +
  scale_x_continuous(breaks = seq(1900, 2030, 10)) +
  labs(
    title    = "Observed vs. Design ARI — 100-year, 24-hr Event",
    subtitle = "TBCMP Project Area (5-yr moving average in green)",
    x = "Year"
  ) +
  theme_wb()

# ── Assemble and save ─────────────────────────────────────────────────────────
combined_plot <- (fig_A | fig_B) /
  (fig_C | fig_D) /
  fig_E +
  plot_annotation(
    title   = "Tampa Bay Coastal Mster Plan — Rainfall Extremes Trend Analysis",
    subtitle = "Methods follow Wright et al. (2019), GRL 46:8144–8153",
    theme   = theme(plot.title    = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 10))
  )

output_dir <- "tbcmp_rainfall_output"
dir.create(output_dir, showWarnings = FALSE)

ggsave(file.path(output_dir, "fig_combined.png"),
       combined_plot, width = 18, height = 20, dpi = 150)

ggsave(file.path(output_dir, "fig_A_station_map.png"),
       fig_A, width = 8, height = 6, dpi = 150)
ggsave(file.path(output_dir, "fig_B_exceedance_trend.png"),
       fig_B, width = 8, height = 5, dpi = 150)
ggsave(file.path(output_dir, "fig_C_trend_heatmap.png"),
       fig_C, width = 5, height = 5, dpi = 150)
ggsave(file.path(output_dir, "fig_D_clusters.png"),
       fig_D, width = 8, height = 5, dpi = 150)
ggsave(file.path(output_dir, "fig_E_ari_comparison.png"),
       fig_E, width = 8, height = 5, dpi = 150)

message(sprintf("\n✓ Figures saved to '%s/'.", output_dir))


# ── 13. Export Summary Tables ─────────────────────────────────────────────────

write_csv(stations_final,
          file.path(output_dir, "stations_final.csv"))

write_csv(annual_counts_full,
          file.path(output_dir, "annual_exceedance_counts.csv"))

write_csv(
  nb_results |> dplyr::select(ari, pct_change, p_value, converged),
  file.path(output_dir, "nb_regression_results.csv")
)

write_csv(
  cluster_trends |> dplyr::select(year, ari, n_clusters),
  file.path(output_dir, "cluster_counts.csv")
)

write_csv(
  obs_ari_smooth |> dplyr::select(year, n_exc, n_active, exc_fraction,
                           observed_ari, ari_smooth),
  file.path(output_dir, "observed_ari.csv")
)

message("✓ Summary tables exported.")
message("\n═══ Analysis complete. ═══════════════════════════════════════════════")


# ── 14. Session / Reproducibility Info ────────────────────────────────────────
# Printed at end so the analyst can include it in supplementary material.
message("\nSession info:")
print(sessionInfo())

