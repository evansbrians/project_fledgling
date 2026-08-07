# 8_detect_movement_states.R  --  step 8: movement states

# Stationary/moving classification from node RSSI variance, via a two-state
# hierarchical HMM. Runs after step 5 (localize_windows.R); independent of the KDE
# steps.

library(cmdstanr)
library(posterior)
library(tidyverse)
library(cli)

library(oxbowR)

# Detections ---------------------------------------------------------------

# Adjust the table and column names to match oxbow.duckdb. The model needs
# raw node-level detections (one row per tag-node-detection), not the
# localizations, since the variance is computed within each node.

con <-
  DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = "/Volumes/ssd980/oxbow_2026/oxbow.duckdb",
    read_only = TRUE
  )

detections <-
  tbl(con, "detections") %>%
  select(
    tag_id,
    node_id,
    datetime = time,
    rssi
  ) %>%
  collect()

DBI::dbDisconnect(con)

# Field test ---------------------------------------------------------------

# Fit the technician walk-and-set-down trial first, where the true state is
# known, to check that the two fitted emission means separate where they
# should and to choose window_minutes.

field_windows <-
  field_test_detections %>%
  build_rssi_windows(
    window_minutes = 60,
    tz = "America/New_York"
  )

field_model <- fit_rssi_state_hmm(field_windows)

field_states <- extract_states(field_model)

# Confusion against the hand-scored truth (truth: tag_id, datetime, state):

field_states %>%
  left_join(
    field_truth,
    by = join_by(tag_id, datetime),
    suffix = c("_fit", "_true")
  ) %>%
  count(state_true, state_fit)

# Full run -----------------------------------------------------------------

windows <-
  detections %>%
  build_rssi_windows(
    window_minutes = 60,
    step_minutes = 60,
    min_detections = 10,
    ref_hours = c(22, 3),
    tz = "America/New_York"
  )

model <- fit_rssi_state_hmm(windows)

model$fit$summary(
  variables =
    c(
      "mu_stat",
      "log_delta",
      "sigma",
      "tau_stat",
      "tau_delta",
      "p_stay_stationary",
      "p_stay_moving"
    )
)

states <- extract_states(model)

# Inspect --------------------------------------------------------------------

# Broad pass over one tag, then a narrower interval:

plot_rssi_states(
  model,
  states,
  tag = "2A661E33"
)

plot_rssi_states(
  model,
  states,
  tag = "2A661E33",
  date_range =
    c(
      ymd_hms("2026-05-14 00:00:00"),
      ymd_hms("2026-05-16 00:00:00")
    )
)

# Write ----------------------------------------------------------------------

arrow::write_parquet(
  states,
  "data/tag_states.parquet"
)

# Optionally carry the states onto the localizations:

localizations %>%
  assign_states(states) %>%
  arrow::write_parquet("data/localizations_states.parquet")
