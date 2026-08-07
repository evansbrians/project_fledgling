# overview -----------------------------------------------------------------

# 5_localize_windows.R  --  step 5: localization
#
# Stage 2 of the localization workflow: localize each detection window.
#
# Input is the Stage 1 output (detection_windows.parquet): windows already
# subset, gated, and carrying node_lon / node_lat. This stage does NO filtering
# -- it simply localizes every window with the best model:
#   conditional drop of rssi <= curve asymptote, top-8 nodes, multilaterate
#   weighted by 1 / (rmse_node^2 * dist_est^2). Each fix gets a confidence
#   (predicted distance to nearest node) and an expected error.
#
# Calibration: Blair BOBO / EAME tags at 0.5 m, 150 s windowing (per species).
# Work is CHECKPOINTED per tag (one parquet per tag), so it can be stopped and
# resumed; a combined localizations.parquet is assembled at the end.
#
# Output fields: tag_id, species, window_id, window_start_time,
#                window_end_time, lon, lat, n_nodes, conf_min_dist, exp_err_m,
#                weight
#
# Requires: tidyverse, arrow, lubridate, gslnls, minpack.lm
# (fs, cli, glue, rlang via ::)

# setup --------------------------------------------------------------------

library(tidyverse)
library(arrow)
library(lubridate)
library(oxbowR)
library(fs)
library(cli)
library(glue)

# configuration ------------------------------------------------------------

config <-
  list(
    windows_path = "data/localization/detection_windows.parquet",
    tags_csv = "data/tags_with_shortened_ids.csv",
    cal_csv = "data/calibration/cal_rssi_v_dist_blair.csv",
    out_dir = "data/localization/localizations_output",
    combined_out = "data/localization/localizations.parquet",
    species_keep = c("BOBO", "EAME"),
    cal_height = 0.5,
    window_s = 150,
    top_k = 8,
    max_dist = 1000
  )

set.seed(2358)
dir_create(config$out_dir)

# run ----------------------------------------------------------------------

windows_dataset <- arrow::open_dataset(config$windows_path)

present_tags <-
  windows_dataset %>%
  distinct(tag_id) %>%
  collect()

tags <-
  read_csv(config$tags_csv, show_col_types = FALSE) %>%
  filter(
    species %in% config$species_keep
  ) %>%
  semi_join(present_tags, by = "tag_id")

calibration <-
  read_csv(config$cal_csv, show_col_types = FALSE) %>%
  mutate(
    rssi = as.numeric(rssi)
  )

curves <-
  config$species_keep %>%
  set_names() %>%
  map(
    \(.species) {
      fit_species_curve(
        calibration, .species,
        config$cal_height, config$window_s, config$max_dist
      )
    }
  )

localizations <-
  tags %>%
  select(tag_id, species) %>%
  pmap(
    \(tag_id, species) {
      localize_tag(
        tag_id, species,
        curves, windows_dataset, config$out_dir,
        config$top_k, config$max_dist
      )
    }
  ) %>%
  bind_rows()

arrow::write_parquet(localizations, config$combined_out)

localizations %>%
  summarize(
    n_fixes = n(),
    n_days = n_distinct(as_date(window_start_time)),
    median_conf_m = median(conf_min_dist),
    .by = c(tag_id, species)
  ) %>%
  write_csv(glue("{config$out_dir}/localization_summary.csv"))
