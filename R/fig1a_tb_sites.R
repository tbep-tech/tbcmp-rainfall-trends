###############################################################################
# Figure 1a — Daily Rainfall with IDF Thresholds and Exceedance Markers
# Tampa Bay Watershed GHCN Stations, 1950–2017
#
# Reproduces the structure of Figure 1a in:
#   Wright, D.B., Bosma, C.D., & Lopez-Cantu, T. (2019).
#   GRL 46, 8144–8153. https://doi.org/10.1029/2019GL083235
#
# The original Figure 1a shows, for a single GHCN station (Chapel Hill, NC):
#   • Vertical bars  — calendar 1-day rainfall (mm) for every year 1950–2017
#   • Three horizontal reference lines — Atlas 14 IDF thresholds for the
#     10-year, 100-year, and 1000-year, 24-hr events
#   • Coloured exceedance markers — bars that cross each threshold are
#     highlighted in the threshold's colour
#   • Annotation — for each threshold: number of observed exceedances and
#     the binomial probability of seeing that many or more in a stationary
#     climate (p(k ≥ observed | n=68, p=1/ARI))
#
# This script replicates that figure for every station present in idf_all,
# producing one panel per station.  It assumes the parent script
# (tampa_bay_rainfall_extremes.R) has already been sourced so that the
# following objects exist in the global environment:
#
#   prcp_qc        — daily PRCP records (station_id, date, year, prcp)
#   idf_all        — tidy Atlas 14 thresholds (station_id, duration_tag,
#                    ari_yr, idf_mm)
#   stations_final — station metadata (station_id, name, latitude, longitude)
#   STUDY_START, STUDY_END — integer scalar (1950, 2017)
#
# Outputs:
#   tampa_bay_rainfall_output/fig1a_<station_id>.png  — one file per station
#   tampa_bay_rainfall_output/fig1a_all_stations.pdf  — all stations in one PDF
###############################################################################


# ── 0. Prerequisites ──────────────────────────────────────────────────────────

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)

# Verify that required upstream objects exist
required_objects <- c("prcp_qc", "idf_all", "stations_final",
                      "STUDY_START", "STUDY_END")
missing_objects  <- required_objects[!sapply(required_objects, exists)]
if (length(missing_objects) > 0) {
  stop(
    "The following objects are missing from the environment.\n",
    "Please source tampa_bay_rainfall_extremes.R first:\n  ",
    paste(missing_objects, collapse = ", ")
  )
}

output_dir <- "tbcmp_rainfall_output"
dir.create(output_dir, showWarnings = FALSE)

# Study period length — used in the binomial probability calculation
N_YEARS <- STUDY_END - STUDY_START + 1L   # 126 years (1900–2025 inclusive)

# ARI values to draw as threshold lines, similar to Figure 1a
PLOT_ARIS <- c(10L, 25L, 100L)

# Colours for each threshold line/marker — match the paper's red palette
ARI_COLOURS <- c(
  "10"   = "#fc8d59",   # orange-red  (10-year)
  "25"  = "#d7191c",   # red         (25-year)
  "100" = "#7f0000"    # dark red    (100-year)
)

ARI_LABELS <- c(
  "10"   = "10-year",
  "25"  = "25-year",
  "100" = "100-year"
)


# ── 1. Helper: binomial exceedance probability ────────────────────────────────
# p(X >= k | n, p) where p = 1/ARI.  Matches the annotation in Figure 1a.
# The paper quotes p(exactly k exceedances), which is pbinom(k, n, 1/ARI, FALSE)
# for the single-event probability, but the annotation text says
# "p(N exceedances)" — meaning the exact binomial probability of observing
# exactly that count given a stationary 1/ARI probability each year.
binomial_exact_prob <- function(k, n, ari) {
  dbinom(k, size = n, prob = 1 / ari)
}


# ── 2. Helper: extract 24-hr IDF thresholds for one station ──────────────────
# idf_all contains all ARIs; we filter to the 24-hr duration and the three
# PLOT_ARIS.  The 1000-year ARI is stored as ari_yr = 1000.
get_thresholds <- function(station_id_val, idf_all_df, aris = PLOT_ARIS) {
  idf_all_df |>
    dplyr::filter(
      station_id   == station_id_val,
      duration_tag == "24hr",
      ari_yr       %in% aris
    ) |>
    dplyr::select(ari_yr, idf_mm) |>
    dplyr::mutate(
      ari_chr   = as.character(ari_yr),
      label     = ARI_LABELS[ari_chr],
      colour    = ARI_COLOURS[ari_chr]
    )
}


# ── 3. Helper: build the annual-maximum 1-day series for one station ──────────
# Figure 1a plots the ANNUAL MAXIMUM calendar 1-day rainfall (not every day).
# The bar for each year represents the single largest daily total that year.
# After retrieving the max, we apply the Atlas 14 adjustment factor (/1.13)
# to convert calendar-day to 24-hr equivalent — matching how thresholds are
# applied in the main analysis (§2.3 of the paper / §5 of the main script).
get_annual_max_series <- function(station_id_val, prcp_qc_df) {
  prcp_qc_df |>
    dplyr::filter(station_id == station_id_val, !is.na(prcp)) |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      ann_max_raw = max(prcp, na.rm = TRUE),
      .groups     = "drop"
    ) |>
    dplyr::mutate(
      # Adjustment factor converts calendar 1-day max to 24-hr equivalent
      ann_max_24hr = ann_max_raw / 1.13
    ) |>
    # Ensure the full study period is present (years with no data → NA bar)
    dplyr::right_join(
      tibble::tibble(year = STUDY_START:STUDY_END),
      by = "year"
    ) |>
    dplyr::arrange(year)
}


# ── 4. Helper: annotate exceedances and compute binomial probabilities ─────────
annotate_exceedances <- function(annual_series, thresholds_df) {
  # For each ARI threshold, flag years where ann_max_24hr exceeds the threshold
  # and compute the exact binomial probability of that exceedance count.
  purrr::map_dfr(seq_len(nrow(thresholds_df)), function(i) {
    thr_row  <- thresholds_df[i, ]
    ari_val  <- thr_row$ari_yr
    thr_mm   <- thr_row$idf_mm

    exc_years <- annual_series |>
      dplyr::filter(!is.na(ann_max_24hr), ann_max_24hr >= thr_mm)

    k      <- nrow(exc_years)
    p_binom <- binomial_exact_prob(k, N_YEARS, ari_val)

    tibble::tibble(
      ari_yr     = ari_val,
      ari_chr    = as.character(ari_val),
      threshold  = thr_mm,
      n_exc      = k,
      p_binom    = p_binom,
      colour     = thr_row$colour,
      label      = thr_row$label,
      # Annotation string matching the paper's format:
      # "100-year: 1 exceedance; p(1 exceedance) = 0.347"
      annot_text = sprintf(
        "%s: %d exceedance%s; p(%d) = %.3f",
        thr_row$label, k, if (k == 1) "" else "s", k, p_binom
      )
    )
  })
}


# ── 5. Core plotting function ─────────────────────────────────────────────────
# Produces one ggplot object that faithfully replicates Figure 1a for a
# single station.
make_fig1a <- function(station_id_val,
                       prcp_qc_df   = prcp_qc,
                       idf_all_df   = idf_all,
                       stations_df  = stations_final) {

  # Station display name
  stn_meta <- stations_df |>
    dplyr::filter(station_id == station_id_val)
  stn_name <- if (nrow(stn_meta) > 0 && !is.na(stn_meta$name[1])) {
    stn_meta$name[1]
  } else {
    station_id_val
  }
  stn_coords <- if (nrow(stn_meta) > 0) {
    sprintf("%.4f°N, %.4f°W", stn_meta$latitude[1], abs(stn_meta$longitude[1]))
  } else ""

  # Build data layers
  ann_series  <- get_annual_max_series(station_id_val, prcp_qc_df)
  thresholds  <- get_thresholds(station_id_val, idf_all_df)

  # If this station has no IDF data, skip gracefully
  if (nrow(thresholds) == 0) {
    message(sprintf("  ⚠  No IDF thresholds found for %s — skipping.", station_id_val))
    return(NULL)
  }

  exc_summary <- annotate_exceedances(ann_series, thresholds)

  # Mark individual bar colours: a year gets the darkest threshold it exceeds
  # (matching the paper, where the bar colour reflects the highest exceedance)
  bar_colours <- ann_series |>
    dplyr::mutate(bar_col = "steelblue3") |>           # default colour
    dplyr::left_join(
      # For each year, find the highest-ARI threshold exceeded
      purrr::map_dfr(seq_len(nrow(thresholds)), function(i) {
        thr_row <- thresholds[i, ]
        ann_series |>
          dplyr::filter(!is.na(ann_max_24hr),
                        ann_max_24hr >= thr_row$idf_mm) |>
          dplyr::mutate(
            ari_exceeded = thr_row$ari_yr,
            exc_colour   = thr_row$colour
          ) |>
          dplyr::select(year, ari_exceeded, exc_colour)
      }) |>
        dplyr::group_by(year) |>
        dplyr::slice_max(ari_exceeded, n = 1, with_ties = FALSE) |>
        dplyr::ungroup(),
      by = "year"
    ) |>
    dplyr::mutate(
      bar_col = dplyr::if_else(!is.na(exc_colour), exc_colour, bar_col)
    )

  # y-axis ceiling: a little above the highest threshold or the max observation
  y_max <- max(
    max(ann_series$ann_max_24hr, na.rm = TRUE),
    max(thresholds$idf_mm)
  ) * 1.08

  # Build annotation label dataframe — stacked inside the plot area near
  # the right edge, one line per threshold (top-down = highest ARI first),
  # matching Figure 1a layout
  n_thr <- nrow(exc_summary)
  annot_df <- exc_summary |>
    dplyr::arrange(dplyr::desc(ari_yr)) |>
    dplyr::mutate(
      x_pos = STUDY_END - 0.5,
      # Place labels just above each threshold line, with a small upward offset
      # so they don't overlap the line itself
      y_pos = threshold + (y_max * 0.025)
    )

  # ── ggplot assembly ──────────────────────────────────────────────────────────
  p <- ggplot() +

    # ── Rainfall bars ──────────────────────────────────────────────────────────
    geom_col(
      data    = bar_colours,
      mapping = aes(x = year, y = ann_max_24hr, fill = bar_col),
      width   = 0.85,
      colour  = NA
    ) +
    scale_fill_identity() +

    # ── IDF threshold lines ───────────────────────────────────────────────────
    # Draw in order lowest → highest so higher thresholds render on top
    purrr::map(
      seq_len(nrow(thresholds)),
      function(i) {
        thr_row <- thresholds[i, ]
        geom_hline(
          yintercept = thr_row$idf_mm,
          colour     = thr_row$colour,
          linewidth  = 0.9,
          linetype   = "solid"
        )
      }
    ) +

    # ── Threshold labels (right-aligned, above each line) ─────────────────────
    geom_text(
      data    = annot_df,
      mapping = aes(x = x_pos, y = y_pos, label = annot_text, colour = colour),
      hjust   = 1,
      vjust   = 0,
      size    = 2.9,
      fontface = "plain",
      show.legend = FALSE
    ) +
    scale_colour_identity() +

    # ── Scales and labels ─────────────────────────────────────────────────────
    scale_x_continuous(
      breaks = seq(1950, 2020, 10),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.02)),
      limits = c(0, y_max),
      labels = scales::label_number(accuracy = 1)
    ) +
    labs(
      title    = sprintf("1-day rainfall: %s", stn_name),
      subtitle = sprintf(
        "Station %s  |  %s  |  %d–%d  |  N = %d station-years",
        station_id_val, stn_coords, STUDY_START, STUDY_END, N_YEARS
      ),
      x = "Year",
      y = "1-Day Rainfall (mm), 24-hr adjusted"
    ) +

    # ── Theme — matches the paper's clean, minimal style ─────────────────────
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "grey40"),
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(size = 9, colour = "grey30"),
      axis.title         = element_text(size = 10),
      axis.text          = element_text(size = 9),
      plot.margin        = margin(8, 12, 6, 8)
    )

  return(p)
}


# ── 6. Generate plots for all stations in idf_all ────────────────────────────
# idf_all contains one row per station × duration_tag × ari_yr.
# We take the unique set of station IDs that have at least one 24-hr threshold.

stations_with_idf <- idf_all |>
  dplyr::filter(duration_tag == "24hr") |>
  dplyr::pull(station_id) |>
  unique() |>
  sort()

message(sprintf(
  "\nBuilding Figure 1a panels for %d stations with Atlas 14 IDF data …",
  length(stations_with_idf)
))

plot_list <- vector("list", length(stations_with_idf))
names(plot_list) <- stations_with_idf

for (i in seq_along(stations_with_idf)) {
  sid <- stations_with_idf[i]
  message(sprintf("  [%d/%d] %s", i, length(stations_with_idf), sid))

  plot_list[[sid]] <- tryCatch(
    make_fig1a(sid),
    error = function(e) {
      message(sprintf("    ✗ Error: %s", conditionMessage(e)))
      NULL
    }
  )
}

# Drop any stations that returned NULL
plot_list <- purrr::compact(plot_list)
message(sprintf("  ✓ %d panels produced.", length(plot_list)))


# ── 7. Save individual PNGs ───────────────────────────────────────────────────
message("\nSaving individual station PNGs …")

purrr::iwalk(plot_list, function(p, sid) {
  fname <- file.path(output_dir,
                     paste0("fig1a_", gsub(":", "_", sid), ".png"))
  ggplot2::ggsave(fname, p, width = 10, height = 4, dpi = 150)
})

message(sprintf("  ✓ Individual PNGs saved to '%s/'.", output_dir))


# ── 8. Save all-station PDF (one page per station) ────────────────────────────
# A single PDF is convenient for reviewing all stations at once.

message("Saving combined PDF (one page per station) …")

pdf_path <- file.path(output_dir, "fig1a_all_stations.pdf")
pdf(pdf_path, width = 11, height = 4.5, onefile = TRUE)
for (p in plot_list) print(p)
dev.off()

message(sprintf("  ✓ Combined PDF saved: %s", pdf_path))


# ── 9. Optional multi-panel overview (up to 6 stations per page) ─────────────
# A compact grid view: useful for quick visual comparison across stations.
# Uses patchwork; limited to the first 24 stations to keep file size reasonable.

if (length(plot_list) >= 2) {

  n_overview <- min(length(plot_list), 24L)
  overview_plots <- plot_list[seq_len(n_overview)]

  # Reduce text size for the small multi-panel layout
  overview_plots <- purrr::map(overview_plots, function(p) {
    p + theme(
      plot.title    = element_text(size = 8,  face = "bold"),
      plot.subtitle = element_text(size = 6.5),
      axis.title    = element_text(size = 7.5),
      axis.text     = element_text(size = 7),
      plot.margin   = margin(4, 6, 3, 4)
    )
  })

  # 2-column layout — each row is one station
  n_cols   <- 2L
  n_rows   <- ceiling(n_overview / n_cols)
  combined <- purrr::reduce(overview_plots, `+`) +
    patchwork::plot_layout(ncol = n_cols) +
    patchwork::plot_annotation(
      title    = "Figure 1a — Daily Rainfall with Atlas 14 IDF Thresholds",
      subtitle = sprintf(
        "Tampa Bay Watershed GHCN stations | 1950–2017 | first %d of %d stations shown",
        n_overview, length(plot_list)
      ),
      theme = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10)
      )
    )

  overview_path <- file.path(output_dir, "fig1a_overview_grid.png")
  ggplot2::ggsave(
    overview_path, combined,
    width  = 18,
    height = n_rows * 4.2,
    dpi    = 130,
    limitsize = FALSE
  )
  message(sprintf("  ✓ Overview grid saved: %s", overview_path))
}

message("\n═══ fig1a_tampa_bay.R complete. ════════════════════════════════════════")
