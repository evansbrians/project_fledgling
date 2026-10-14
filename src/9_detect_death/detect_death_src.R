# detect_death_src.R --------------------------------------------------------

# Functions for the death change-point model (v4, exact posterior).
# Workflow in detect_death.R; channel calibration in the project notes.

library(tidyverse)

# daily evidence -------------------------------------------------------------

# Detections and nodes reached per tag-day, within deployment (UTC):

collect_daily_node_stats <-
  function(con, tag_ids) {
    deployments <-
      tbl(con, "deployed_tags") %>%

      # Subset to the requested tags:

      filter(tag_id %in% tag_ids) %>%

      # Deployments end at recovery or year's end:

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

    tbl(con, "raw") %>%
      inner_join(
        deployments,
        by = join_by(tag_id)
      ) %>%
      filter(
        time > deploy_time,
        time < end_time
      ) %>%

      # UTC day number (tz-cast-free):

      mutate(
        date = floor(epoch(time) / 86400)
      ) %>%

      # Detections per tag, day, and node:

      count(
        tag_id,
        date,
        node_id
      ) %>%

      # Nodes reached per tag-day:

      summarize(
        n_det =
          sum(n, na.rm = TRUE),
        n_nodes = n(),
        .by = c(tag_id, date)
      ) %>%
      collect() %>%

      # Day number to Date:

      mutate(date = as_date(date))
  }

# One modeled row per observed tag-day:

build_death_days <-
  function(node_days,
           move_days,
           deployed_tags,
           nests) {
    deployments <-
      deployed_tags %>%
      mutate(
        end_time =
          coalesce(
            recover_time,
            ceiling_date(deploy_time, unit = "year")
          ),
        deploy_day = floor(as.numeric(deploy_time) / 86400),
        end_day = floor(as.numeric(end_time) / 86400),
        recovered = !is.na(recover_time)
      ) %>%
      distinct(
        tag_id,
        species,
        deploy_day,
        end_day,
        recovered
      )

    node_days %>%

      # Drop corrupt detection dates:

      filter(is.finite(date)) %>%
      left_join(
        move_days,
        by = join_by(tag_id, date)
      ) %>%
      inner_join(
        deployments,
        by = join_by(tag_id)
      ) %>%
      left_join(
        nests,
        by = join_by(tag_id)
      ) %>%
      mutate(
        n_win = replace_na(n_win, 0L),
        n_move = replace_na(n_move, 0L),
        still_run = replace_na(still_run, 0L),
        move_day = n_move > n_win / 2,
        day_num = as.numeric(date),
        days_post_fledge = as.numeric(date - fledge_date)
      ) %>%

      # Keep the deployment containing each day:

      filter(
        day_num >= deploy_day,
        day_num <= end_day
      ) %>%
      arrange(tag_id, date) %>%

      # Alive baselines from moving days; movement seen before each day:

      mutate(
        base_det = median(n_det[move_day]),
        node_base = median(n_nodes[move_day]),
        moved_before =
          lag(
            cumany(move_day),
            default = FALSE
          ),
        gap_after =
          replace_na(lead(day_num) - day_num - 1, 0),
        .by = tag_id
      ) %>%
      mutate(
        base_det = coalesce(base_det, median(n_det)),
        node_base = coalesce(node_base, median(n_nodes)),
        .by = tag_id
      ) %>%
      mutate(
        log_ratio = log(n_det / base_det),
        deploy_age = day_num - deploy_day,
        is_late =
          species == "BOBO" & month(date) >= 11
      )
  }

# channels -------------------------------------------------------------------

# Student-t log density with location and scale:

log_dt <-
  function(x, par, nu) {
    dt(
      (x - par[1]) / par[2],
      df = nu,
      log = TRUE
    ) - log(par[2])
  }

# Daily log likelihood of the detection ratio in one state:

ratio_loglik <-
  function(log_ratio, par, nu) {
    if (is.list(par)) {
      log(
        par$w * exp(log_dt(log_ratio, par$dense, nu)) +
          (1 - par$w) * exp(log_dt(log_ratio, par$sparse, nu))
      )
    } else {
      log_dt(log_ratio, par, nu)
    }
  }

# Daily log likelihood of the daytime movement category in one state:
# moving / still run / neither on full days, moving / not on sparse days,
# nest-phase rates until the bird has moved.

movement_loglik <-
  function(tag_days, q, pars) {
    full <- tag_days$n_win >= pars$min_win

    nest_phase <-
      tag_days$days_post_fledge <= pars$nest_days & !tag_days$moved_before

    still <-
      full & !nest_phase & tag_days$still_run >= pars$still_hours

    q_move <-
      case_when(
        nest_phase ~ q$move_nest,
        full ~ q$move,
        .default = q$move_sparse
      )

    q_still <-
      if_else(full & !nest_phase, q$still, 0)

    case_when(
      tag_days$n_win == 0 ~ 0,
      tag_days$move_day ~ log(q_move),
      still ~ log(q_still),
      .default = log(1 - q_move - q_still)
    )
  }

# posterior ------------------------------------------------------------------

# Exact death posterior for one tag. Hypotheses: death on an observed day,
# tag destroyed at record end (predation/haying), or alive (departed):

death_posterior <-
  function(tag_days, pars) {
    n <- nrow(tag_days)

    # Daily log likelihood, alive state, with the silence that follows:

    la <-
      movement_loglik(tag_days, pars$alive, pars) +
      ratio_loglik(tag_days$log_ratio, pars$alive$ratio, pars$nu) +
      log(1 - pars$p_sil_alive) +
      tag_days$gap_after * log(pars$p_sil_alive) +
      if_else(tag_days$is_late, log(pars$late_alive_p), 0)

    # Daily log likelihood, dead state:

    ld <-
      movement_loglik(tag_days, pars$dead, pars) +
      ratio_loglik(tag_days$log_ratio, pars$dead$ratio, pars$nu) +
      log(1 - pars$p_sil_dead) +
      tag_days$gap_after * log(pars$p_sil_dead)

    ca <- cumsum(la)

    cd <- cumsum(ld)

    # Terminal silence after death: per-day silence or a lost signal:

    s_tail <- tag_days$end_day[n] - tag_days$day_num[n]

    silence_dead <-
      if_else(
        s_tail == 0,
        0,
        max(s_tail * log(pars$p_sil_dead), log(pars$p_lost_dead))
      )

    # Death on day k: alive before, dead after, then silence:

    lp_die <-
      log(pars$w_die) - log(n) +
      lag(ca, default = 0) +
      cd[n] - lag(cd, default = 0) +
      silence_dead

    # Node spike soon after deployment marks aerial predation:

    tail_days <- slice_tail(tag_days, n = 2)

    aerial <-
      any(
        tail_days$deploy_age <= pars$aerial_max_age &
          tail_days$n_nodes > pars$aerial_nodes_mult * tail_days$node_base &
          tail_days$n_nodes >= pars$aerial_min_nodes
      )

    # Destroyed at record end: alive throughout, then silent:

    lp_destroy <-
      if_else(
        s_tail >= pars$min_destroy_silence,
        log(pars$w_destroy) + ca[n] +
          if_else(aerial, log(pars$aerial_boost), 0),
        -Inf
      )

    # Departure needs the season and an independent bird; a record
    # that stops cold from a dense day is rarely a departure:

    cold_stop <-
      tag_days$log_ratio[n] >= log(pars$cold_stop_ratio)

    p_depart <-
      plogis(
        (yday(tag_days$date[n]) - pars$depart_mid) / pars$depart_scale
      ) %>%
      pmin(pars$depart_max) %>%
      pmax(pars$depart_min) *
      if_else(tag_days$days_post_fledge[n] < pars$depart_min_age, 0.1, 1) *
      if_else(cold_stop, pars$cold_stop_depart, 1)

    silence_alive <-
      if_else(
        s_tail == 0,
        0,
        max(s_tail * log(pars$p_sil_alive), log(p_depart))
      )

    # Alive at end, unlikely once the tag was recovered:

    lp_alive <-
      log(pars$w_alive) + ca[n] + silence_alive +
      if_else(tag_days$recovered[n], log(pars$p_recover_alive), 0)

    # Normalize hypothesis weights:

    lp <- c(lp_die, lp_destroy, lp_alive)

    w <- exp(lp - max(lp))

    w <- w / sum(w)

    # Destroyed mass lands on the last day:

    p_death_day <- w[1:n]

    p_death_day[n] <- p_death_day[n] + w[n + 1]

    k <- which.max(p_death_day)

    list(
      p_death_day = p_death_day,
      p_died = 1 - w[n + 2],
      death_date = tag_days$date[k],
      p_death_date = p_death_day[k],
      kind = c(rep("died", n), "destroyed", "alive")[which.max(w)],
      aerial = aerial
    )
  }
