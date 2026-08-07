
# Calibration using Blair's data

# setup -------------------------------------------------------------------

library(celltracktech)
library(DBI)
library(duckdb)
library(duckplyr)
library(sf)
library(tmap)
library(tidyverse)

tmap_mode("view")

# Custom functions:

library(oxbowR)

# Connect to database:

oxbow_db <- 
  dbConnect(
    duckdb(), 
    dbdir = "./data/oxbow.duckdb", 
    read_only = FALSE
  )

# Basemap for mapping:

basemap <- 
  leaflet() %>% 
  addProviderTiles(
    provider = providers$Esri.WorldImagery,
    options = tileOptions(maxZoom = 21)
  )

# Blair tags:

calibration_tags <-
  read_rds("data/calibration/calibration_tags_blair.rds")

# process bad elf ---------------------------------------------------------

# Bad elf file:

bad_elf <- 
  read_rds("data/calibration/bad_elf_blair.rds")

# Convert to a shapefile (for checking):

bad_elf_pts <-
  bad_elf %>% 
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326
  )

# Have a look:

basemap %>% 
  addPolylines(
    data = 
      bad_elf_pts %>% 
      summarize(do_union = FALSE) %>% 
      st_cast("LINESTRING"),
    weight = 2.75,
    opacity = 0.7,
    dashArray = "2, 5",
    color = "#ffff00"
  )

# get node health from database and process -------------------------------

# Nodes are sometimes back in the office or in the field house, so load the
# extent:

oxbow_bbox <-
  st_read("data/oxbow_extent.geojson") %>% 
  st_bbox()

# Get node health files:

node_locations <- 
  tbl(oxbow_db, "node_health") %>% 
  filter(
    
    # Node health files were incomplete, so I will get node health files for 1
    # week before and after the start times:
    
    time >= !!min(bad_elf$datetime) - days(7),
    time <= !!max(bad_elf$datetime) + days(7),
    longitude > !!oxbow_bbox[["xmin"]],
    longitude < !!oxbow_bbox[["xmax"]],
    latitude  > !!oxbow_bbox[["ymin"]],
    latitude  < !!oxbow_bbox[["ymax"]]
  ) %>% 
  select(node_id, longitude:latitude) %>% 
  collect() %>% 
  
  # Get node locations:
  
  calculate_node_locations() %>% 
  select(
    node_id,
    lon = median_lon,
    lat = median_lat
  )

# Have a look:

basemap %>% 
  addPolylines(
    data = 
      bad_elf_pts %>% 
      summarize(do_union = FALSE) %>% 
      st_cast("LINESTRING"),
    weight = 2.75,
    opacity = 0.7,
    dashArray = "2, 5",
    color = "#ffff00"
  ) %>% 
  addCircleMarkers(
    data = node_locations, 
    lng = ~ lon, 
    lat = ~ lat, 
    radius = 5, 
    color = "#000000", 
    opacity = 0.9,
    weight = 2,
    fillColor = "#00ffff", 
    fillOpacity = 0.9, 
    label = ~ node_id
  )

# get detections from database and process --------------------------------

detections <- 
  tbl(oxbow_db, "raw") %>% 
  filter(
    
    # Query based on start and end times:
    
    time >= !!min(bad_elf$datetime),
    time <= !!max(bad_elf$datetime),
    
    # Query based on the tags used for calibration:
    
    tag_id %in% calibration_tags$tag_id
  ) %>% 
  select(
    time,
    tag_id:node_id,
    rssi = tag_rssi
  ) %>% 
  collect()

# Disconnect from the database:

DBI::dbDisconnect(oxbow_db)

# rssi vs. distance -------------------------------------------------------

rssi_v_dist <- 
  detections %>% 
  mutate(
    time = round_date(time, "1 secs")
  ) %>% 
  arrange(time) %>% 
  
  # Add the node locations for a given detection:
  
  inner_join(
    node_locations %>% 
      select(
        node_id,
        node_lon = lon,
        node_lat = lat
      ),
    by = "node_id"
  ) %>% 
  
  # Add bad elf locations for a given detection:
  
  inner_join(
    bad_elf %>% 
      mutate(
        time = round_date(datetime, "1 secs")
      ) %>% 
      select(
        time = datetime,
        gps_lon = longitude,
        gps_lat = latitude
      ),
    by = "time"
  ) %>% 
  
  # Calculate the distance between the node and gps for each detection (note:
  # haversine is a celltrack function, I adjusted it so that it's more than
  # twice as fast and uses ~15% less memory):
  
  mutate(
    distance = 
      haversine(
        node_lat, 
        node_lon, 
        gps_lat, 
        gps_lon
      )
  ) %>% 
  
  # Add species:
  
  left_join(
    calibration_tags,
    by = "tag_id"
  )

# Have a look too make sure the Bad Elf and node times are aligned (due mainly
# to the time zone question):

rssi_v_dist %>% 
  ggplot() +
  aes(
    x = distance,
    y = rssi,
    color = node_id
  ) +
  facet_wrap(height_m ~ tag_type) +
  geom_point()

# For an extra double-check, it might be good to go individually through the
# nodes (in case a node clock is off):

rssi_v_dist %>% 
  split(.$node_id) %>% 
  map(
    ~ .x %>% 
      ggplot() +
      aes(
        x = distance,
        y = rssi,
        color = node_id
      ) +
      facet_wrap(~ tag_type) +
      geom_point()
  )

# All looks reasonable!

# fit calibration model ---------------------------------------------------

calibration_fits <-
  c(
    bobo = "BOBO",
    eame = "EAME"
  ) %>% 
  map(
    \(.species) {
      rssi_v_dist_data <-
        rssi_v_dist %>% 
        filter(height_m == 2) %>% 
        mutate(
          rssi = as.numeric(rssi),
          .keep = "unused"
        ) %>% 
        filter(tag_type == .species)
      
      fit <- 
        gslnls::gsl_nls(
          rssi ~ a - b * exp(-c * distance),
          rssi_v_dist_data,
          start = 
            list(
              a = -105, 
              b = -60, 
              c = 0.17
            ),
          control = 
            gslnls::gsl_nls_control(
              maxiter = 1000,
              scale = "levenberg",
              solver = "svd"
            )
        )
      list(
        species = .species,
        data = rssi_v_dist_data,
        fit = fit,
        coefs = coef(fit)
      )
    }
  )

# calibration error -------------------------------------------------------

# Get calibration error:

calibration_error <-
  calibration_fits %>% 
  keep(~ !is.null(.x)) %>% 
  map(
    ~ get_calibration_error(.x)
  ) %>% 
  bind_rows()

# Overall calibration error by species:

calibration_error_summary <-
  calibration_error %>% 
  summarize_distance_error(.by = species)

# Calibration error by node and species:

calibration_error_by_node <-
  calibration_error %>% 
  summarize_distance_error(
    .by = c(species, node_id)
  )

# Calibration error by species and distance:

calibration_error_by_distance <-
  calibration_error %>% 
  mutate(
    distance_bin = cut_width(distance, width = 10, boundary = 0)
  ) %>% 
  summarize_distance_error(
    .by = c(species, distance_bin)
  )

# Calibration error by species, distance, and node:

calibration_error_by_node_distance <- 
  calibration_error %>% 
  mutate(
    distance_bin = cut_width(distance, width = 10, boundary = 0)
  ) %>% 
  summarize_distance_error(
    .by = c(species, node_id, distance_bin)
  )

# exploring calibration error ---------------------------------------------

# Modeling error:

error_fits <-
  c(
    "species",
    "species + node_id",
    "species + rssi",
    "species + node_id + rssi",
    "species + node_id * rssi"
  ) %>% 
  str_c("distance_error ~ ", .) %>% 
  set_names() %>% 
  map(
    ~ glm(
      as.formula(.x),
      data = calibration_error
    )
  ) 

AICcmodavg::aictab(error_fits)

# Distance error by node and species:

calibration_error %>% 
  mutate(
    mean_dist_error = mean(distance_error, na.rm = TRUE),
    .by = node_id
  ) %>% 
  mutate(
    node_id = fct_reorder(node_id, mean_dist_error)
  ) %>% 
  ggplot() +
  aes(
    x = distance_error,
    y = node_id,
    fill = species
  ) +
  geom_boxplot(outliers = FALSE) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  )

# Is there evidence of node-specific error?

calibration_error_by_node %>% 
  mutate(
    node_id = fct_reorder(node_id, bias_m)
  ) %>% 
  ggplot() +
  aes(
    x = bias_m,
    y = node_id
  ) +
  geom_boxplot() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  )

# Is there evidence of distance-specific error?

calibration_error_by_distance %>% 
  ggplot() +
  aes(
    x = mean_distance_m,
    y = rmse_m,
    color = species
  ) +
  geom_point() +
  geom_line() +
  labs(
    x = "Known distance to node (m)",
    y = "Distance RMSE (m)"
  )

# Is there evidence of distance-specific error by node?

calibration_error_by_node_distance %>% 
  filter(species == "BOBO") %>%
  ggplot() +
  aes(
    x = bias_m,
    y = rmse_m,
    color = node_id
  ) +
  geom_point() +
  geom_line() +
  facet_wrap(~ species)
labs(
  x = "Known distance to node (m)",
  y = "Distance RMSE (m)"
)

# node specific calibration? ----------------------------------------------

# Fit model by tag type and node:

calibration_fits_node_specific <-
  c(
    bobo = "BOBO",
    eame = "EAME"
  ) %>% 
  map(
    \(.species) {
      
      rssi_v_dist_list <- 
        rssi_v_dist %>% 
        filter(tag_type == .species) %>% 
        filter(distance < 1000) %>% 
        filter(
          n() > 200,
          .by = node_id
        ) %>%
        mutate(rssi = as.numeric(rssi)) %>% 
        split(.$node_id)
      
      rssi_v_dist_fits <-
        rssi_v_dist_list %>% 
        map(
          \(.node_data) {
            fit <-
              gslnls::gsl_nls(
                fn = rssi ~ a - b * exp(-c * distance),
                data = .node_data,
                start = 
                  list(
                    a = -105, 
                    b = -60, 
                    c = 0.17
                  ),
                control = 
                  gslnls::gsl_nls_control(
                    maxiter = 1000,
                    scale = "levenberg",
                    solver = "svd",
                    warnOnly = TRUE
                  )
              )
            list(
              data = .node_data,
              fit = fit,
              coefs = coef(fit)
            )
          }
        )
    }
  )

# Get calibration error:

calibration_error_node_specific <- 
  calibration_fits_node_specific %>% 
  map(
    \(.species_list) {
      .species_list %>% 
        map(
          \(.node_list) {
            get_calibration_error(.node_list)
          }
        ) %>% 
        bind_rows()
    }
  ) %>% 
  bind_rows()

# Overall calibration error by species and node:

calibration_error_summary_node <-
  calibration_error_node_specific %>% 
  summarize_distance_error(
    .by = c(species)
  )

# Calibration error by species and node:

calibration_error_node_specific %>% 
  summarize_distance_error(
    .by = c(species, node_id)
  ) %>% 
  arrange(node_id) %>% 
  print(n = Inf)

# Calibration error by species, distance, and node:

calibration_error_node_specific %>% 
  mutate(
    distance_bin = cut_width(distance, width = 10, boundary = 0)
  ) %>% 
  summarize_distance_error(
    .by = c(species, distance_bin, node_id)
  ) %>% 
  arrange(distance_bin)

# explore which calibration model/usage works best for localization -------

# Note: While it's more honest to use cross-validation (given that detections
# below contribute to both the distance estimation and localization), during
# this exploration phase I'm going to keep it simple!

# Make an events file:

set.seed(2358)

detection_events_start <- 
  rssi_v_dist %>% 
  
  # Add coefficients:
  
  left_join(
    calibration_error %>% 
      select(species, a:c) %>% 
      distinct(),
    by = join_by(tag_type == species)
  ) %>% 
  
  # Add grouping:
  
  group_by(time, tag_type, tag_id) %>% 
  
  # Filter to minimum detection criteria:
  
  filter(
    n() >= 3,
    !all(rssi < -104)
  ) %>% 
  ungroup()

# Get detection events:

detection_events <- 
  lst(
    BOBO = "BOBO",
    EAME = "EAME"
  ) %>% 
  map(
    ~ detection_events_start %>% 
      filter(tag_type == .x) %>% 
      group_by(time) %>% 
      group_split() %>% 
      sample(size = 100)
  ) %>% 
  bind_rows()

## method 1: global coefficient score with no weighting -------------------

### bobo ------------------------------------------------------------------

max_distance <- 450

# Predict distances:

events_with_distance <- 
  detection_events %>% 
  filter(tag_type == "BOBO") %>% 
  mutate(
    rssi_ceiling = a - b,
    rssi_asymptote = a,
    log_arg = (rssi - rssi_asymptote) / (-b),
    dist_raw = -log(log_arg) / c,
    dist_estimate =
      case_when(
        rssi >= rssi_ceiling ~ 1,
        rssi <= rssi_asymptote ~ max_distance,
        .default =
          pmax(dist_raw, 1) %>%
          pmin(max_distance)
      )
  ) %>% 
  select(!tag_type:dist_raw)

# Fit models:

events_fit_global <-
  events_with_distance %>%
  split(.$time) %>% 
  map(
    \(.x) {
      
      # Start with the closest node (that we have the most confidence in):
      
      start_node <-
        slice_max(
          .x,
          order_by = rssi,
          with_ties = FALSE
        )
      
      tryCatch(
        minpack.lm::nlsLM(
          dist_estimate ~
            haversine(
              lat1 = node_lat,
              lon1 = node_lon,
              lat2 = ml_lat,
              lon2 = ml_lon
            ),
          data = .x,
          start =
            list(
              ml_lat = start_node$node_lat,
              ml_lon = start_node$node_lon
            ),
          control =
            minpack.lm::nls.lm.control(
              maxiter = 200,
              ftol = .Machine$double.eps,
              ptol = .Machine$double.eps
            )
        ),
        error = function(e) NULL,
        warning = function(w) NULL
      ) %>% 
        broom::tidy() %>% 
        select(term:estimate) %>% 
        pivot_wider(
          names_from = term,
          values_from = estimate
        ) %>% 
        bind_cols(.x, .) %>% 
        mutate(
          ml_dist = 
            haversine(
              lon1 = ml_lon,
              lat1 = ml_lat,
              lon2 = gps_lon,
              lat2 = gps_lat
            )
        )
    }
  ) %>% 
  bind_rows()

events_fit_global %>% 
  get_summary_stats(ml_dist)

## method 2: weight by node error -----------------------------------------

events_with_node_weights <-
  events_with_distance %>% 
  left_join(
    calibration_error_by_node %>% 
      filter(species == "BOBO") %>% 
      select(node_id, node_weight = rmse_m) %>% 
      mutate(
        node_weight =
          node_weight %>% 
          replace_na(
            median(., na.rm = TRUE)
          ) %>% 
          { 1 / .^2}
      ),
    by = "node_id"
  )

# Fit models:

events_fit_node_weighted <-
  events_with_node_weights %>% 
  split(.$time) %>% 
  map(
    \(.x) {
      
      # Start with the closest node (that we have the most confidence in):
      
      start_node <-
        slice_max(
          .x,
          order_by = rssi,
          with_ties = FALSE
        )
      
      mod <- 
        tryCatch(
          minpack.lm::nlsLM(
            dist_estimate ~
              haversine(
                lat1 = node_lat,
                lon1 = node_lon,
                lat2 = ml_lat,
                lon2 = ml_lon
              ),
            data = .x,
            weights = node_weight,
            start =
              list(
                ml_lat = start_node$node_lat,
                ml_lon = start_node$node_lon
              ),
            control =
              minpack.lm::nls.lm.control(
                maxiter = 200,
                ftol = .Machine$double.eps,
                ptol = .Machine$double.eps
              )
          ),
          error = function(e) NULL,
          warning = function(w) NULL
        )
      
      if (!is.null(mod)) {
        mod %>% 
        broom::tidy() %>% 
          select(term:estimate) %>% 
          pivot_wider(
            names_from = term,
            values_from = estimate
          ) %>% 
          bind_cols(.x, .) %>% 
          mutate(
            ml_dist = 
              haversine(
                lon1 = ml_lon,
                lat1 = ml_lat,
                lon2 = gps_lon,
                lat2 = gps_lat
              )
          )
      } else {
        NULL
      }
    }
  ) %>% 
  bind_rows()

lst(
  global = events_fit_global,
  node_weighted = events_fit_node_weighted
) %>% 
  imap(
    \(.x, .name) {
      get_summary_stats(.x, ml_dist) %>% 
        mutate(
          method = .name,
          .before = 1
        )
    }
  ) %>% 
  bind_rows()

# This showed a slight improvement over the previous, but not a huge one.

### method 3: weight by node and mean distance error ----------------------

events_with_distance_distance_weights <-
  events_with_distance %>% 
  bind_rows() %>% 
  mutate(
    distance_bin = 
      cut_width(
        dist_estimate, 
        width = 10, 
        boundary = 0
      )
  ) %>% 
  left_join(
    calibration_error_by_node_distance %>% 
      filter(species == "BOBO") %>% 
      select(node_id, distance_bin, node_weight = rmse_m) %>% 
      mutate(
        node_weight =
          node_weight %>% 
          replace_na(
            median(., na.rm = TRUE)
          ) %>% 
          { 1 / .^2}
      ),
    by = c("node_id", "distance_bin")
  )

# Fit models:

events_fit_node_distance_weighted <-
  events_with_distance_distance_weights %>% 
  split(.$time) %>% 
  map(
    \(.x) {
      
      # Start with the closest node (that we have the most confidence in):
      
      start_node <-
        slice_max(
          .x,
          order_by = rssi,
          with_ties = FALSE
        )
      
      tryCatch(
        minpack.lm::nlsLM(
          dist_estimate ~
            haversine(
              lat1 = node_lat,
              lon1 = node_lon,
              lat2 = ml_lat,
              lon2 = ml_lon
            ),
          data = .x,
          weights = node_weight,
          start =
            list(
              ml_lat = start_node$node_lat,
              ml_lon = start_node$node_lon
            ),
          control =
            minpack.lm::nls.lm.control(
              maxiter = 200,
              ftol = .Machine$double.eps,
              ptol = .Machine$double.eps
            )
        ),
        error = function(e) NULL,
        warning = function(w) NULL
      ) %>% 
        broom::tidy() %>% 
        select(term:estimate) %>% 
        pivot_wider(
          names_from = term,
          values_from = estimate
        ) %>% 
        bind_cols(.x, .) %>% 
        mutate(
          ml_dist = 
            haversine(
              lon1 = ml_lon,
              lat1 = ml_lat,
              lon2 = gps_lon,
              lat2 = gps_lat
            )
        )
    }
  ) %>% 
  bind_rows()

lst(
  global = events_fit_global,
  node_weighted = events_fit_node_weighted,
  distance_weighted = events_fit_node_distance_weighted
) %>% 
  imap(
    \(.x, .name) {
      get_summary_stats(.x, ml_dist) %>% 
        mutate(
          method = .name,
          .before = 1
        )
    }
  ) %>% 
  bind_rows()

# This was a big improvement over the previous! Not where we want to be (max,
# sd, and se are actually worse), but it's starting to look better.

### method 4: offset by mean absolute distances ----------------------------

events_with_distance_bias_m <-
  events_with_distance %>% 
  bind_rows() %>% 
  mutate(
    distance_bin = 
      cut_width(
        dist_estimate, 
        width = 10, 
        boundary = 0
      )
  ) %>% 
  left_join(
    calibration_error_by_node_distance %>% 
      filter(species == "BOBO") %>% 
      select(node_id, distance_bin, bias_m),
    by = c("node_id", "distance_bin")
  ) %>% 
  mutate(
    dist_estimate =
      dist_estimate + bias_m
  )

# Fit models:

events_fit_node_distance_bias_m <-
  events_with_distance_bias_m %>% 
  split(.$time) %>% 
  map(
    \(.x) {
      
      # Start with the closest node (that we have the most confidence in):
      
      start_node <-
        slice_max(
          .x,
          order_by = rssi,
          with_ties = FALSE
        )
      
      tryCatch(
        minpack.lm::nlsLM(
          dist_estimate ~
            haversine(
              lat1 = node_lat,
              lon1 = node_lon,
              lat2 = ml_lat,
              lon2 = ml_lon
            ),
          data = .x,
          start =
            list(
              ml_lat = start_node$node_lat,
              ml_lon = start_node$node_lon
            ),
          control =
            minpack.lm::nls.lm.control(
              maxiter = 200,
              ftol = .Machine$double.eps,
              ptol = .Machine$double.eps
            )
        ),
        error = function(e) NULL,
        warning = function(w) NULL
      ) %>% 
        broom::tidy() %>% 
        select(term:estimate) %>% 
        pivot_wider(
          names_from = term,
          values_from = estimate
        ) %>% 
        bind_cols(.x, .) %>% 
        mutate(
          ml_dist = 
            haversine(
              lon1 = ml_lon,
              lat1 = ml_lat,
              lon2 = gps_lon,
              lat2 = gps_lat
            )
        )
    }
  ) %>% 
  bind_rows()

lst(
  global = events_fit_global,
  node_weighted = events_fit_node_weighted,
  distance_weighted = events_fit_node_distance_weighted,
  bias_adjusted = events_fit_node_distance_bias_m
) %>% 
  imap(
    \(.x, .name) {
      get_summary_stats(.x, ml_dist) %>% 
        mutate(
          method = .name,
          .before = 1
        )
    }
  ) %>% 
  bind_rows()

# Lowest SD, but not so great

### method 5: offset by mean absolute distances and bias ------------------

events_with_distance_bias_m_weights <-
  events_with_distance %>% 
  bind_rows() %>% 
  mutate(
    distance_bin = 
      cut_width(
        dist_estimate, 
        width = 10, 
        boundary = 0
      )
  ) %>% 
  left_join(
    calibration_error_by_node_distance %>% 
      filter(species == "BOBO") %>% 
      select(node_id, distance_bin, bias_m, node_weight = rmse_m) %>% 
      mutate(
        node_weight =
          node_weight %>% 
          replace_na(
            median(., na.rm = TRUE)
          ) %>% 
          { 1 / .^2}
      ),
    by = c("node_id", "distance_bin")
  ) %>% 
  mutate(
    dist_estimate =
      dist_estimate + bias_m
  )

# Fit models:

events_fit_node_distance_bias_m_weighted <-
  events_with_distance_bias_m_weights %>% 
  split(.$time) %>% 
  map(
    \(.x) {
      
      # Start with the closest node (that we have the most confidence in):
      
      start_node <-
        slice_max(
          .x,
          order_by = rssi,
          with_ties = FALSE
        )
      
      tryCatch(
        minpack.lm::nlsLM(
          dist_estimate ~
            haversine(
              lat1 = node_lat,
              lon1 = node_lon,
              lat2 = ml_lat,
              lon2 = ml_lon
            ),
          data = .x,
          weights = node_weight,
          start =
            list(
              ml_lat = start_node$node_lat,
              ml_lon = start_node$node_lon
            ),
          control =
            minpack.lm::nls.lm.control(
              maxiter = 200,
              ftol = .Machine$double.eps,
              ptol = .Machine$double.eps
            )
        ),
        error = function(e) NULL,
        warning = function(w) NULL
      ) %>% 
        broom::tidy() %>% 
        select(term:estimate) %>% 
        pivot_wider(
          names_from = term,
          values_from = estimate
        ) %>% 
        bind_cols(.x, .) %>% 
        mutate(
          ml_dist = 
            haversine(
              lon1 = ml_lon,
              lat1 = ml_lat,
              lon2 = gps_lon,
              lat2 = gps_lat
            )
        )
    }
  ) %>% 
  bind_rows()

lst(
  global = events_fit_global,
  node_weighted = events_fit_node_weighted,
  distance_weighted = events_fit_node_distance_weighted,
  bias_adjusted = events_fit_node_distance_bias_m,
  bias_adjusted_weighted = events_fit_node_distance_bias_m_weighted
) %>% 
  imap(
    \(.x, .name) {
      get_summary_stats(.x, ml_dist) %>% 
        mutate(
          method = .name,
          .before = 1
        )
    }
  ) %>% 
  bind_rows()

# method 6: node specific -------------------------------------------------

events_node_specific_df <-
  bind_rows(events_with_distance) %>% 
  left_join(
    
    # Get nodes-specific coefficients:
    
    calibration_fits_node_specific %>% 
      pluck("bobo") %>% 
      imap(
        ~ bind_rows(.x[["coefs"]]) %>% 
          mutate(
            node_id = .y,
            .before = 1
          )
      ) %>% 
      bind_rows(),
    by = "node_id"
  ) %>% 
  left_join(
    calibration_error_by_node %>% 
      filter(species == "BOBO") %>% 
      select(node_id, node_weight = rmse_m) %>% 
      mutate(
        node_weight =
          node_weight %>% 
          replace_na(
            median(., na.rm = TRUE)
          ) %>% 
          { 1 / .^2}
      ),
    by = "node_id"
  ) %>% 
  mutate(
    rssi_ceiling = a - b,
    rssi_asymptote = a,
    log_arg = (rssi - rssi_asymptote) / (-b),
    dist_raw = -log(log_arg) / c,
    dist_estimate =
      case_when(
        rssi >= rssi_ceiling ~ 1,
        rssi <= rssi_asymptote ~ max_distance,
        .default =
          pmax(dist_raw, 1) %>%
          pmin(max_distance)
      )
  )

events_fit_node_distance_by_node <-
  events_node_specific_df %>% 
  split(.$time) %>% 
  map(
    \(.x) {
      
      # Start with the closest node (that we have the most confidence in):
      
      start_node <-
        slice_max(
          .x,
          order_by = rssi,
          with_ties = FALSE
        )
      
      tryCatch(
        minpack.lm::nlsLM(
          dist_estimate ~
            haversine(
              lat1 = node_lat,
              lon1 = node_lon,
              lat2 = ml_lat,
              lon2 = ml_lon
            ),
          data = .x,
          weights = node_weight,
          start =
            list(
              ml_lat = start_node$node_lat,
              ml_lon = start_node$node_lon
            ),
          control =
            minpack.lm::nls.lm.control(
              maxiter = 200,
              ftol = .Machine$double.eps,
              ptol = .Machine$double.eps
            )
        ),
        error = function(e) NULL,
        warning = function(w) NULL
      ) %>% 
        broom::tidy() %>% 
        select(term:estimate) %>% 
        pivot_wider(
          names_from = term,
          values_from = estimate
        ) %>% 
        bind_cols(.x, .) %>% 
        mutate(
          ml_dist = 
            haversine(
              lon1 = ml_lon,
              lat1 = ml_lat,
              lon2 = gps_lon,
              lat2 = gps_lat
            )
        )
    }
  ) %>% 
  bind_rows()

lst(
  global = events_fit_global,
  node_weighted = events_fit_node_weighted,
  distance_weighted = events_fit_node_distance_weighted,
  bias_adjusted = events_fit_node_distance_bias_m,
  bias_adjusted_weighted = events_fit_node_distance_bias_m_weighted,
  by_node = events_fit_node_distance_by_node
) %>% 
  imap(
    \(.x, .name) {
      get_summary_stats(.x, ml_dist) %>% 
        mutate(
          method = .name,
          .before = 1
        )
    }
  ) %>% 
  bind_rows()

# write to file -----------------------------------------------------------

# RSSI-to-distance calibration table (read by steps 3 and 5).

write_csv(
  rssi_v_dist %>%
    select(
      time, tag_id, tag_type, height_m, node_id,
      node_lon, node_lat, gps_lon, gps_lat, rssi
    ),
  "data/calibration/cal_rssi_v_dist_blair.csv"
)
