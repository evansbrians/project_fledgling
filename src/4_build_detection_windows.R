# overview -----------------------------------------------------------------

# 4_build_detection_windows.R  --  step 4: filtering
#
# Stage 1 of the localization workflow: subset + window the raw detections.
#
# For each tag, detections are grouped into fixed time windows (default 150 s).
# Each (tag, window, node) is reduced to its MAX rssi (our localization
# measure). Node positions (monthly medians from node_health, to track node
# drift) are attached, and each window is gated to the nodes within `gate_m` of
# its strongest-rssi node. Only windows that could then be localized are kept:
# at least `min_nodes` distinct gated nodes and at least one above the noise
# floor. Every surviving window gets a unique `window_id`.
#
# Doing the gate here (rather than at localization) gives the true count of
# localizable windows per tag and leaves the localization stage as pure
# localization -- all filtering happens once, in this script.
#
# The heavy lifting runs inside DuckDB (the `raw` table has ~290M rows), so the
# whole build is one query. Output is either a flat parquet file or a DuckDB
# view / table -- set by `output_type`. All filtering thresholds are function
# arguments, so they are easy to relax if a tag returns too few windows.
#
# Output fields: tag_id, window_id, window_start_time, window_end_time,
#                node_id, rssi, node_lon, node_lat
#
# Requires: tidyverse, DBI, duckdb, glue

# setup --------------------------------------------------------------------

library(tidyverse)
library(DBI)
library(duckdb)
library(glue)

# custom functions ---------------------------------------------------------

library(oxbowR)

# configuration ------------------------------------------------------------

config <-
  list(
    db_path = "data/oxbow.duckdb",
    tags_csv = "data/tags_with_shortened_ids.csv",
    species_keep = c("BOBO", "EAME"),
    window_s = 150,
    min_nodes = 3,
    rssi_floor = -104,
    min_beacons = 1,
    gate_m = 1000,
    output_type = "parquet",
    dest = "data/localization/detection_windows.parquet"
  )

# run ----------------------------------------------------------------------

tags <-
  read_csv(config$tags_csv, show_col_types = FALSE) %>%
  filter(
    species %in% config$species_keep
  )

con <-
  dbConnect(
    duckdb(),
    dbdir = config$db_path,
    read_only = config$output_type == "parquet"
  )

build_detection_windows(
  .con = con,
  .tag_ids = tags$tag_id,
  .window_s = config$window_s,
  .min_nodes = config$min_nodes,
  .rssi_floor = config$rssi_floor,
  .min_beacons = config$min_beacons,
  .gate_m = config$gate_m,
  .output_type = config$output_type,
  .dest = config$dest
)

dbDisconnect(con, shutdown = TRUE)
