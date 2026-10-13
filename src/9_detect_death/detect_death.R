# Detect death

# Exact change-point posterior: each bird is alive until an unknown day,
# destroyed at record end (predation/haying, with an aerial-predation
# signature), or alive at end (departed). Daily evidence: the daytime
# movement category (moving day, a five-hour stationary run past the
# nest phase, or neither), detection intensity vs the bird's own baseline
# on every day, silence between and after detections, tag recovery, and
# BOBO pings past October 31. Emission rates are fixed constants
# calibrated from the rule-labeled tag-days (see project notes); no MCMC
# is needed.

# setup ---------------------------------------------------------------------

library(fs)
library(tidyverse)

# Get custom functions used:

source("src/9_detect_death/detect_death_src.R")

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
    out_dir = "data/detect_death",
    states_path = "data/movement_states/tag_movement_states_daytime.csv",
    nests_path = "data/ctt_study_nestlocations_fledgedates.csv",
    tz = "America/New_York",
    pars =
      list(
        min_win = 5,
        still_hours = 5,
        nest_days = 7,
        alive =
          list(
            move = 0.875,
            move_sparse = 0.65,
            move_nest = 0.5,
            still = 0.03,
            ratio = c(0, 0.5)
          ),
        dead =
          list(
            move = 0.02,
            move_sparse = 0.02,
            move_nest = 0.02,
            still = 0.8,
            ratio =
              list(
                w = 0.5,
                dense = c(-0.15, 0.8),
                sparse = c(-3.6, 2)
              )
          ),
        nu = 5,
        p_sil_alive = 0.12,
        p_sil_dead = 0.32,
        p_lost_dead = 0.5,
        late_alive_p = 0.1,
        w_die = 0.45,
        w_destroy = 0.15,
        w_alive = 0.4,
        depart_mid = 206,
        depart_scale = 10,
        depart_min = 0.02,
        depart_max = 0.5,
        depart_min_age = 21,
        cold_stop_ratio = 0.1,
        cold_stop_depart = 0.2,
        p_recover_alive = 0.02,
        min_destroy_silence = 3,
        aerial_nodes_mult = 2,
        aerial_min_nodes = 20,
        aerial_max_age = 14,
        aerial_boost = 6
      )
  )

# Create the write directory:

dir_create(config$out_dir)

# daily evidence -------------------------------------------------------------

# Detections and nodes reached per tag-day:

node_days <- collect_daily_node_stats(con, deployed_tags$tag_id)

DBI::dbDisconnect(con)

# Daytime movement states:

states_daytime <-
  read_csv(config$states_path, show_col_types = FALSE)

# Moving windows and the longest stationary run per tag-day:

move_days <-
  states_daytime %>%
  arrange(tag_id, datetime) %>%

  # Number consecutive stationary windows:

  mutate(
    run =
      cumsum(
        state != "stationary" |
          replace_na(datetime - lag(datetime) > dhours(1), TRUE)
      ),
    .by = c(tag_id, date)
  ) %>%
  mutate(
    run_length = sum(state == "stationary"),
    .by = c(tag_id, date, run)
  ) %>%
  summarize(
    n_win = n(),
    n_move = sum(state == "moving"),
    still_run = max(run_length),
    .by = c(tag_id, date)
  )

# Fledge dates:

nests <-
  read_csv(config$nests_path, show_col_types = FALSE) %>%
  select(tag_id, fledge_date = est_fledge_date)

# One modeled row per observed tag-day:

death_days <-
  build_death_days(node_days, move_days, deployed_tags, nests)

# posterior ------------------------------------------------------------------

# Exact death posterior per tag:

death_posteriors <-
  death_days %>%
  distinct(tag_id) %>%
  pull(tag_id) %>%
  set_names() %>%
  map(
    \(.tag_id) {
      death_posterior(
        filter(death_days, tag_id == .tag_id),
        pars = config$pars
      )
    }
  )

# Death-day probability per tag-day:

death_days <-
  death_days %>%
  mutate(
    p_death_day =
      death_posteriors %>%
      map("p_death_day") %>%
      list_c()
  )

# Point estimates per tag:

death_estimates <-
  death_posteriors %>%
  map(
    \(.posterior) {
      tibble(
        death_date = .posterior$death_date,
        p_death_date = .posterior$p_death_date,
        p_died = .posterior$p_died,
        kind = .posterior$kind,
        aerial = .posterior$aerial
      )
    }
  ) %>%
  list_rbind(names_to = "tag_id") %>%

  # First stationary window times the death:

  left_join(
    states_daytime %>%
      filter(state == "stationary") %>%
      summarize(
        death_datetime = min(datetime),
        .by = c(tag_id, date)
      ),
    by = join_by(tag_id, death_date == date)
  ) %>%
  left_join(
    deployed_tags %>%
      distinct(tag_id, species, recover_time),
    by = "tag_id"
  ) %>%
  select(
    tag_id,
    species,
    death_date,
    death_datetime,
    p_death_date,
    p_died,
    kind,
    aerial,
    recover_time
  )

# write ----------------------------------------------------------------------

# Local ISO 8601 with UTC offset:

death_estimates %>%
  mutate(
    death_datetime =
      death_datetime %>%
      with_tz(tzone = config$tz) %>%
      format_ISO8601(usetz = TRUE),
    recover_time =
      recover_time %>%
      with_tz(tzone = config$tz) %>%
      format_ISO8601(usetz = TRUE)
  ) %>%
  write_csv("data/detect_death/tag_death_estimates.csv")

# site data ------------------------------------------------------------------

# Daytime states, in local plot time:

site_states <-
  states_daytime %>%
  arrange(tag_id, datetime) %>%
  mutate(
    datetime =
      datetime %>%
      with_tz(tzone = config$tz) %>%
      format("%Y-%m-%d %H:%M"),
    y = round(y, 3),
    confidence = round(confidence, 3)
  ) %>%
  select(
    tag_id,
    datetime,
    y,
    state,
    confidence
  )

# Death estimates and recoveries, in local plot time:

site_tags <-
  death_estimates %>%
  mutate(
    death_datetime =
      death_datetime %>%
      with_tz(tzone = config$tz) %>%
      format("%Y-%m-%d %H:%M"),
    recover_time =
      recover_time %>%
      with_tz(tzone = config$tz) %>%
      format("%Y-%m-%d %H:%M"),
    confidence = round(p_death_date, 3),
    p_died = round(p_died, 3)
  )

# Death-day posterior, for the site:

site_posterior <-
  death_days %>%
  filter(p_death_day >= 0.001) %>%
  mutate(p = round(p_death_day, 4)) %>%
  select(
    tag_id,
    date,
    p
  )

# Bundle as a JS data file:

str_c(
  "window.STATES_DATA = ",
  jsonlite::toJSON(
    list(
      tags = site_tags,
      states = site_states,
      posterior = site_posterior
    ),
    na = "null"
  ),
  ";"
) %>%
  write_lines("brian_sandbox/states_data.js")
