# Oxbow sun times

# Dawn, sunrise, sunset, and dusk at the Oxbow site for every day the tag
# movement states span (one row per local date).

# setup --------------------------------------------------------------------

library(suncalc)
library(tidyverse)

# Get custom functions used:

source("src/8_detect_movement_states/rssi_state_hmm_src.R")

# Get state data:

movement_states <-
  read_csv("data/movement_states/tag_movement_states.csv")

# sun times ----------------------------------------------------------------

sun_times <-
  movement_states %>%
  
  # One row per local date represented in the states:
  
  distinct(
    date =
      as_date(datetime, tz = "America/New_York")
  ) %>%
  arrange(date) %>%
  
  # Add each sun event of that date:
  
  mutate(
    dawn = get_sun_times(date, "dawn"),
    sunrise = get_sun_times(date, "sunrise"),
    sunset = get_sun_times(date, "sunset"),
    dusk = get_sun_times(date, "dusk")
  )

# subset movement state data to daytime records ---------------------------

movement_states_daytime <-
  movement_states %>% 
  mutate(
    date = as_date(datetime)
  ) %>% 
  left_join(sun_times, by = "date") %>% 
  filter(
    between(
      datetime,
      dawn,
      sunset
    )
  )

# inspect -----------------------------------------------------------------

# Day and night:

movement_states %>% 
  filter(tag_id == "192D2D4B") %>% 
  ggplot() +
  aes(
    x = datetime,
    y = y,
    color = state
  ) +
  geom_line(
    colour = "grey70",
    linewidth = 0.3
  ) +
  geom_point(
    size = 1.4
  ) +
  scale_colour_manual(
    values =
      c(
        stationary = "#2c7fb8",
        moving = "#d95f0e"
      )
  ) +
  scale_alpha_manual(
    values =
      c(
        `TRUE` = 0.4,
        `FALSE` = 1
      ),
    guide = "none"
  ) +
  labs(
    y = "log pooled within-node RSSI variance",
    colour = "State",
    title = "192D2D4B"
  ) +
  theme_minimal()

# Subset to just daytime:

movement_states_daytime %>% 
  filter(tag_id == "192D2D4B") %>% 
  ggplot() +
  aes(
    x = datetime,
    y = y,
    color = state
  ) +
  geom_line(
    colour = "grey70",
    linewidth = 0.3
  ) +
  geom_point(
    size = 1.4
  ) +
  scale_colour_manual(
    values =
      c(
        stationary = "#2c7fb8",
        moving = "#d95f0e"
      )
  ) +
  scale_alpha_manual(
    values =
      c(
        `TRUE` = 0.4,
        `FALSE` = 1
      ),
    guide = "none"
  ) +
  labs(
    y = "log pooled within-node RSSI variance",
    colour = "State",
    title = "192D2D4B"
  ) +
  theme_minimal()

# write -------------------------------------------------------------------

movement_states_daytime %>% 
  write_csv(
    "data/movement_states/tag_movement_states_daytime.csv"
  )


