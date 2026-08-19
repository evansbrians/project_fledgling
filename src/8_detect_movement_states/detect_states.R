# Detect movement states

# Stationary/moving states are classified from node RSSI variance with a
# two-state HMM. The stationary state is defined per species from the
# late-night (22:00-03:00 local) windows of all tags; each tag is then fit
# independently and written to a csv in data/movement_states.

# Note: This runs on the *original* detection data (the raw table),
# independent of the localization workflow.

# setup --------------------------------------------------------------------

library(cmdstanr)
library(fs)
library(glue)
library(tidyverse)

# Get custom functions used:

source("src/8_detect_movement_states/rssi_state_hmm_src.R")

# Connect to the database:

con <-
  DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = "data/oxbow.duckdb",
    read_only = TRUE
  )

# Deployment tags:

deployed_tags <-
  tbl(con, "deployed_tags") %>%

  # Subset to species of interest:

  filter(
    species %in% c("BOBO", "EAME")
  ) %>%

  # Pass records to memory:

  collect()

# Config file defines the configuration of the runs:

config <-
  list(
    out_dir = "data/movement_states/movement_states_raw",
    window_minutes = 60,
    step_minutes = 60,
    min_detections = 10,
    ref_hours = c(22, 3),
    utc_night_hours = c(1, 9),
    tz = "America/New_York"
  )

# Create the write directory:

dir_create(config$out_dir)

# Compile the Stan model (it's shared across all tags):

stan_model <- cmdstan_model("src/8_detect_movement_states/rssi_state_hmm.stan")

# reference state ----------------------------------------------------------

# Define the stationary state per species from the late-night windows of
# all tags:

references <-
  c("BOBO", "EAME") %>%
  set_names() %>%
  map(
    \(.species) {
      species_tags <-
        deployed_tags %>%
        filter(species == .species) %>%
        pull(tag_id)

      collect_detections(
        con,
        species_tags,
        utc_hours = config$utc_night_hours
      ) %>%
        build_windows() %>%
        summarize_reference()
    }
  )

# states -------------------------------------------------------------------

# Fit each tag and write its states to csv. Finished tags are skipped;
# failed models are reported and skipped:

deployed_tags %>%
  distinct(tag_id, species) %>%
  pwalk(
    \(tag_id, species) {
      out_path <- glue("{config$out_dir}/tag_states_{tag_id}.csv")

      if (file_exists(out_path)) {
        return(NULL)
      }

      rlang::try_fetch(
        {
          tag_states <-
            collect_detections(con, tag_id) %>%
            build_windows() %>%
            detect_tag_states(
              stan_model = stan_model,
              reference = references[[species]],
              tz = config$tz
            )

          # Write datetime as local ISO 8601 with its UTC offset:

          tag_states %>%
            mutate(
              datetime =
                format_ISO8601(datetime, usetz = TRUE)
            ) %>%
            write_csv(out_path)
        },
        error = function(cnd) {
          cli::cli_warn(
            "States failed for tag {.val {tag_id}}.",
            parent = cnd
          )
        }
      )
    }
  )

# Combine the tag csvs, with datetime displayed in local time:

states <-
  dir_ls(config$out_dir, glob = "*.csv") %>%
  read_csv(show_col_types = FALSE) %>%
  mutate(
    datetime =
      with_tz(datetime, tzone = config$tz)
  )

# inspect ------------------------------------------------------------------

# Grab a tag with a known fate:

plot_states_all <-
  plot_rssi_states(
    states,
    tag = "192D2D4B",
    tz = config$tz
  )

plot_rssi_states(
  states,
  tag = "192D2D4B",
  tz = config$tz,
  date_range =
    c(
      ymd_hms("2023-06-24 00:00:00", tz = config$tz),
      ymd_hms("2023-06-27 00:00:00", tz = config$tz)
    )
)

# Disconnect from the database:

DBI::dbDisconnect(con)

# write --------------------------------------------------------------------

# Save the combined copy:

states %>%
  write_csv(
    "data/movement_states/tag_movement_states.csv"
  )
