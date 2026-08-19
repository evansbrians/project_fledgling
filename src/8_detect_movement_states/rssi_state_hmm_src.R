# rssi_state_hmm.R ---------------------------------------------------------

# Functions for stationary/moving state detection from node RSSI variance.
# Companion to rssi_state_hmm.stan; workflow in detect_states.R.

library(cmdstanr)
library(posterior)
library(tidyverse)

# database -----------------------------------------------------------------

# Query a set of tags' detections, subset to each tag's deployment period
# before extraction (all times UTC):

collect_detections <-
  function(con,
           tag_ids,
           utc_hours = NULL) {
    deployments <-
      tbl(con, "deployed_tags") %>%

      # Subset to the requested tags:

      filter(tag_id %in% tag_ids) %>%

      # Deployments end at recovery or, if unrecovered, at year's end:

      mutate(
        end_time =
          coalesce(
            recover_time,
            make_date(year(deploy_time) + 1L, 1L, 1L)
          )
      ) %>%
      select(
        tag_id,
        deploy_time,
        end_time
      )

    detections <-
      tbl(con, "raw") %>%

      # Subset to detections within each tag's deployment:

      inner_join(
        deployments,
        by = join_by(tag_id)
      ) %>%
      filter(
        time > deploy_time,
        time < end_time
      )

    # Optionally subset to late-night hours (UTC):

    if (!is.null(utc_hours)) {
      detections <-
        detections %>%
        filter(
          hour(time) >= !!utc_hours[1],
          hour(time) < !!utc_hours[2]
        )
    }

    detections %>%
      select(
        tag_id,
        node_id,
        datetime = time,
        rssi = tag_rssi
      ) %>%

      # Pass records to memory:

      collect()
  }

# time windows --------------------------------------------------------------

# Window a set of detections with the shared settings in config:

build_windows <-
  function(detections) {
    detections %>%
      build_rssi_windows(
        window_minutes = config$window_minutes,
        step_minutes = config$step_minutes,
        min_detections = config$min_detections,
        ref_hours = config$ref_hours,
        tz = config$tz
      )
  }

# TRUE for clock hours inside the reference period (may wrap midnight):

is_night_hour <-
  function(clock_hour,
           ref_hours = c(22, 3)) {
    if (ref_hours[1] < ref_hours[2]) {
      clock_hour >= ref_hours[1] & clock_hour < ref_hours[2]
    } else {
      clock_hour >= ref_hours[1] | clock_hour < ref_hours[2]
    }
  }

# Summarize detections into per tag-window RSSI variance statistics
# (step_minutes < window_minutes gives sliding windows):

build_rssi_windows <-
  function(detections,
           window_minutes = 60,
           step_minutes = window_minutes,
           min_detections = 10,
           min_nodes = 1,
           ref_hours = c(22, 3),
           tz = "UTC") {
    if (nrow(detections) == 0) {
      cli::cli_abort("{.arg detections} has no rows.")
    }

    if (window_minutes %% step_minutes != 0) {
      cli::cli_abort(
        "{.arg window_minutes} must be a whole multiple of {.arg step_minutes}."
      )
    }

    n_offsets <- window_minutes / step_minutes

    step_secs <- step_minutes * 60

    origin <- min(detections$datetime)

    # Detection count, sum, and sum of squares per tag, node, and bin:

    bin_stats <-
      detections %>%
      mutate(
        elapsed =
          as.numeric(datetime - origin, units = "secs"),
        bin = floor(elapsed / step_secs)
      ) %>%
      summarize(
        n = n(),
        sum_rssi = sum(rssi),
        sum_sq = sum(rssi^2),
        .by = c(tag_id, node_id, bin)
      )

    # Add each bin to every window that contains it:

    node_stats <-
      bin_stats %>%
      cross_join(
        tibble(offset = 0:(n_offsets - 1))
      ) %>%
      mutate(window = bin - offset) %>%
      summarize(
        n = sum(n),
        sum_rssi = sum(sum_rssi),
        sum_sq = sum(sum_sq),
        .by = c(tag_id, node_id, window)
      )

    node_stats %>%
      filter(n >= 2) %>%

      # Pool within-node variances, weighted by degrees of freedom:

      mutate(
        node_mean = sum_rssi / n,
        node_var = (sum_sq - n * node_mean^2) / (n - 1)
      ) %>%
      summarize(
        n_detections = sum(n),
        n_nodes = n(),
        pooled_var = sum((n - 1) * node_var) / sum(n - 1),
        between_node_var = var(node_mean),
        .by = c(tag_id, window)
      ) %>%
      filter(
        n_detections >= min_detections,
        n_nodes >= min_nodes
      ) %>%

      # Flag reference windows (fully inside the local night period):

      mutate(
        datetime = origin + window * step_secs,
        local_start =
          with_tz(datetime, tzone = tz),
        local_end = local_start + minutes(window_minutes) - seconds(1),
        ref_start = is_night_hour(hour(local_start), ref_hours),
        ref_end = is_night_hour(hour(local_end), ref_hours),
        is_ref = ref_start & ref_end,
        y = log(pooled_var + 0.01)
      ) %>%
      select(!c(local_start, local_end, ref_start, ref_end)) %>%
      arrange(tag_id, datetime)
  }

# Split each tag's windows into gap-free runs (a detection gap is not a
# state transition); runs shorter than min_length are dropped:

add_segments <-
  function(windows,
           min_length = 3) {
    windows %>%
      arrange(tag_id, datetime) %>%
      mutate(
        new_run = replace_na(window - lag(window) != 1, TRUE),
        .by = tag_id
      ) %>%
      mutate(segment = cumsum(new_run)) %>%
      mutate(
        run_length = n(),
        .by = segment
      ) %>%
      filter(run_length >= min_length) %>%
      select(!c(new_run, run_length))
  }

# reference summary ---------------------------------------------------------

# Location and scale of the night windows across the tag population -- the
# prior that defines the stationary state:

summarize_reference <-
  function(windows) {
    reference_y <- windows$y[windows$is_ref]

    n_reference <- length(reference_y)

    if (n_reference < 30) {
      cli::cli_warn(
        "Only {n_reference} reference windows: weak stationary prior."
      )
    }

    list(
      mu_ref_mean = median(reference_y),
      mu_ref_sd = max(mad(reference_y), 0.25)
    )
  }

# fitting -------------------------------------------------------------------

# Fit the HMM to one tag's segmented windows (ref_leak = prior probability
# of moving at night; nu = Student-t df of the emissions):

fit_rssi_state_hmm <-
  function(tag_windows,
           stan_model,
           reference,
           ref_leak = 0.05,
           nu = 7,
           chains = 4,
           iter_warmup = 1000,
           iter_sampling = 1000,
           ...) {
    tag_windows <-
      tag_windows %>%
      mutate(row_id = row_number())

    segments <-
      tag_windows %>%
      summarize(
        seg_start = min(row_id),
        seg_len = n(),
        .by = segment
      )

    stan_data <-
      list(
        n_win = nrow(tag_windows),
        n_seg = nrow(segments),
        y = tag_windows$y,
        is_ref = as.integer(tag_windows$is_ref),
        seg_start = segments$seg_start,
        seg_len = segments$seg_len,
        mu_ref_mean = reference$mu_ref_mean,
        mu_ref_sd = reference$mu_ref_sd,
        delta_mean = log(1.5),
        delta_sd = 0.75,
        alpha_gamma = rbind(c(9, 1), c(3, 7)),
        ref_leak = ref_leak,
        nu = nu,
        save_states = 0L
      )

    fit <-
      stan_model$sample(
        data = stan_data,
        chains = chains,
        parallel_chains = chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        ...
      )

    list(
      stan_model = stan_model,
      fit = fit,
      windows = tag_windows,
      stan_data = stan_data
    )
  }

# state extraction ----------------------------------------------------------

# Posterior state probabilities, via a generated-quantities pass over
# thinned draws. Returns the window statistics with state and confidence
# added (datetime in tz):

extract_states <-
  function(model,
           threshold = 0.5,
           n_draws = 200,
           tz = "UTC") {
    param_draws <-
      model$fit$draws(
        variables =
          c(
            "mu_stat",
            "log_delta",
            "sigma",
            "gamma_rows",
            "rho"
          )
      )

    n_available <- ndraws(param_draws)

    thin_by <- max(floor(n_available / n_draws), 1)

    state_data <- model$stan_data

    state_data$save_states <- 1L

    state_fit <-
      model$stan_model$generate_quantities(
        fitted_params =
          thin_draws(param_draws, thin = thin_by),
        data = state_data
      )

    p_move_mean <-
      state_fit$draws("p_moving") %>%
      summarize_draws(mean) %>%
      pull(mean)

    model$windows %>%

      # Classify each window; confidence = probability of that state:

      mutate(
        datetime =
          with_tz(datetime, tzone = tz),
        p_moving = p_move_mean,
        state =
          if_else(
            p_moving > threshold,
            "moving",
            "stationary"
          ),
        confidence =
          if_else(
            state == "moving",
            p_moving,
            1 - p_moving
          )
      ) %>%
      select(
        tag_id,
        datetime,
        window,
        n_detections,
        n_nodes,
        pooled_var,
        between_node_var,
        is_ref,
        y,
        state,
        confidence
      )
  }

# detection -----------------------------------------------------------------

# Segment, fit, and classify one tag (mapped across tags in
# detect_states.R):

detect_tag_states <-
  function(tag_windows,
           stan_model,
           reference,
           threshold = 0.5,
           tz = "UTC",
           ...) {
    segmented <-
      tag_windows %>%
      add_segments()

    if (nrow(segmented) == 0) {
      cli::cli_abort("No gap-free window runs to fit.")
    }

    fit_rssi_state_hmm(
      segmented,
      stan_model = stan_model,
      reference = reference,
      ...
    ) %>%
      extract_states(
        threshold = threshold,
        tz = tz
      )
  }

# Rolling join of window states onto detections or localizations:

assign_states <-
  function(target,
           states) {
    target %>%
      left_join(
        states,
        by = join_by(tag_id, closest(datetime >= datetime))
      )
  }

# visualization -------------------------------------------------------------

# Log RSSI variance over time, colored by state; reference windows dimmed:

plot_rssi_states <-
  function(states,
           tag,
           date_range = NULL,
           tz = "UTC") {
    plot_data <-
      states %>%
      filter(tag_id == tag) %>%
      mutate(
        datetime =
          with_tz(datetime, tzone = tz),
        state =
          factor(
            state,
            levels = c("stationary", "moving")
          )
      )

    if (!is.null(date_range)) {
      plot_data <-
        plot_data %>%
        filter(
          datetime >= date_range[1],
          datetime <= date_range[2]
        )
    }

    plot_data %>%
      ggplot() +
      aes(
        x = datetime,
        y = y
      ) +
      geom_line(
        color = "grey70",
        linewidth = 0.3
      ) +
      geom_point(
        aes(
          color = state,
          alpha = is_ref
        ),
        size = 1.4
      ) +
      scale_color_manual(
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
        x = NULL,
        y = "log pooled within-node RSSI variance",
        color = "State",
        title = tag
      ) +
      theme_minimal()
  }

# my sun times ------------------------------------------------------------

get_sun_times <-
  function(.date, .sun) {
    getSunlightTimes(
      date = as_date(.date),
      lat = 38.89,
      lon = -78.16,
      keep = .sun,
      tz = "America/New_York"
    ) %>% 
      pull(.sun)
  }
