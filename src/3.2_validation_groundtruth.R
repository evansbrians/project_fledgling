
library(leaflet)
library(sf)
library(tidyverse)
library(oxbowR)

calibration_pts <- 
  list.files(
    "data/localization/validation/individual_points/",
    full.names = TRUE
  ) %>% 
  map(
    ~ st_read(.x, quiet = TRUE) %>% 
      mutate(
        across(
          bearing,
          ~ as.numeric(.x)
        )
      )
  ) %>% 
  bind_rows() %>% 
  filter(point_name != "another_test")

basemap <- 
  leaflet() %>% 
  addProviderTiles(
    provider = providers$Esri.WorldImagery,
    options = tileOptions(maxZoom = 21)
  )

basemap %>%
  addCircleMarkers(
    data = calibration_pts,
    radius = 5,
    color = "#000000",
    opacity = 0.9,
    weight = 2,
    fillColor = "#00ffff",
    fillOpacity = 0.9
  )

# =========================================================================
# Ground-truth localization of the technician's four tags
# =========================================================================

# A technician carried four tags plus a GPS during the 2026-06-26 field
# campaign. It was not recorded whether each tag is an EAME or a BOBO tag,
# so we localize every tag under BOTH the Blair BOBO 0.5 m and Blair EAME
# 0.5 m calibration curves and infer the type as the curve that gives
# consistently lower error against the true GPS point.
#
# For each calibration point we use the EXACT (time_start, time_end) window
# (not a fixed 90 s window): per node we take the MAX rssi over that window,
# gate to nodes within 1000 m of the strongest node, drop rssi <= the curve
# asymptote, keep the top-8 nodes, and multilaterate weighted by
# 1 / (rmse_node^2 * dist_est^2). This mirrors the validated Python in
# groundtruth/validate_localization.py and the reference workflow in
# calibration_project/test_folder/bird_space_use.R.

library(DBI)
library(duckdb)
library(ks)


# ---- parameters ---------------------------------------------------------

field_tags <-
  c(
    "194C662A", # EAME
    "331E551E", # BOBO
    "19195552", # BOBO
    "4C336155"  # EAME
  )

cal_height <- 0.5
gate_m <- 1000
top_k <- 8
max_dist <- 1000
utm_crs <- 32618

set.seed(2358)

# ---- 1. Blair calibration curves (BOBO and EAME) at 0.5 m ---------------

# The windowed max-rssi rssi ~ distance table (built for the reference
# workflow); it carries tag_type and height_m columns:

cal_raw <-
  read_csv("data/calibration/cal_rssi_v_dist_blair.csv") %>%
  mutate(
    rssi = as.numeric(rssi),
    win = as.integer(as.numeric(time)) %/% 90L
  )


blair_curves <-
  c("BOBO", "EAME") %>%
  set_names() %>%
  map(fit_blair_curve)

# ---- 2. field detections + node locations from the database -------------

# The geojson time windows are local (America/New_York); the database stores
# UTC. Convert each window to UTC for the query:

cal_windows <-
  calibration_pts %>%
  st_drop_geometry() %>%
  as_tibble() %>%
  transmute(
    point_name,
    true_lon =
      st_coordinates(calibration_pts)[, 1],
    true_lat =
      st_coordinates(calibration_pts)[, 2],
    t0 =
      force_tzs(
        ymd_hms(time_start),
        tzones = "America/New_York",
        tzone_out = "UTC"
      ),
    t1 =
      force_tzs(
        ymd_hms(time_end),
        tzones = "America/New_York",
        tzone_out = "UTC"
      )
  )

campaign_start <- min(cal_windows$t0) - hours(6)
campaign_end <- max(cal_windows$t1) + hours(6)

oxbow_db <-
  dbConnect(
    duckdb(),
    dbdir = "data/oxbow.duckdb",
    read_only = TRUE
  )

field_detections <-
  tbl(oxbow_db, "raw") %>%
  filter(
    tag_id %in% field_tags,
    time >= !!campaign_start,
    time <= !!campaign_end
  ) %>%
  select(
    time,
    tag_id,
    node_id,
    rssi = tag_rssi
  ) %>%
  collect() %>%
  mutate(rssi = as.numeric(rssi))

# node_health has no June-2026 rows, so fall back to the most recent
# positions (2026-04+, the closest node deployment to the campaign):

node_locations <-
  tbl(oxbow_db, "node_health") %>%
  filter(
    time >= !!as.POSIXct("2026-04-01", tz = "UTC"),
    !is.na(longitude)
  ) %>%
  select(node_id, longitude, latitude) %>%
  collect() %>%
  calculate_node_locations() %>%
  select(
    node_id,
    node_lon = median_lon,
    node_lat = median_lat
  )

dbDisconnect(oxbow_db, shutdown = TRUE)

# ---- 3. localize one calibration point under one curve ------------------

# ---- 4. localize every point x tag x curve ------------------------------

point_error_all <-
  cal_windows %>%
  crossing(tag_id = field_tags) %>%
  crossing(curve_type = c("BOBO", "EAME")) %>%
  pmap(
    function(...) {
      row <- tibble(...)
      localize_point(row, blair_curves[[row$curve_type]])
    }
  ) %>%
  bind_rows()

write_csv(
  point_error_all,
  "data/localization/validation/point_error_all.csv"
)

# ---- 5. infer tag type, per-point best, and per-tag summary -------------

if (nrow(point_error_all) > 0) {
  tag_type_inference <-
    point_error_all %>%
    filter(!is.na(error_m)) %>%
    summarize(
      median_error_m = median(error_m),
      .by = c(tag_id, curve)
    ) %>%
    slice_min(
      median_error_m,
      n = 1,
      by = tag_id,
      with_ties = FALSE
    ) %>%
    rename(inferred_type = curve)

  write_csv(
    tag_type_inference,
    "data/localization/validation/tag_type_inference.csv"
  )

  point_error_best <-
    point_error_all %>%
    inner_join(
      tag_type_inference %>%
        select(tag_id, curve = inferred_type),
      by = c("tag_id", "curve")
    )

  write_csv(
    point_error_best,
    "data/localization/validation/point_error_best.csv"
  )

  tag_summary <-
    point_error_best %>%
    filter(!is.na(error_m)) %>%
    summarize(
      inferred_type = first(curve),
      n_points = n(),
      median_error_m = median(error_m),
      mean_error_m = mean(error_m),
      rmse_m = sqrt(mean(error_m^2)),
      q80_m = quantile(error_m, 0.80),
      max_error_m = max(error_m),
      median_conf_m = median(conf_min_dist),
      median_se_m = median(se_m, na.rm = TRUE),
      .by = tag_id
    )

  write_csv(
    tag_summary,
    "data/localization/validation/tag_summary.csv"
  )

  # ---- 6. per-tag KDE across that tag's fixes ---------------------------

  make_tag_kde <-
    function(.tag_id) {
      fixes <-
        point_error_best %>%
        filter(
          tag_id == .tag_id,
          !is.na(lon)
        )
      if (nrow(fixes) < 3) return(NULL)
      xy <-
        fixes %>%
        st_as_sf(
          coords = c("lon", "lat"),
          crs = 4326
        ) %>%
        st_transform(utm_crs) %>%
        st_coordinates()
      kd <-
        ks::kde(
          x = xy,
          H = Hpi(xy),
          compute.cont = TRUE
        )
      polys <-
        c(80, 90, 95) %>%
        map(
          function(.level) {
            thr <- kd$cont[[str_c(.level, "%")]]
            cl <-
              contourLines(
                kd$eval.points[[1]],
                kd$eval.points[[2]],
                kd$estimate,
                levels = thr
              )
            if (!length(cl)) return(NULL)
            rings <-
              cl %>%
              map(
                function(.p) {
                  m <- cbind(.p$x, .p$y)
                  if (!all(m[1, ] == m[nrow(m), ])) m <- rbind(m, m[1, ])
                  m
                }
              )
            st_sf(
              tag_id = .tag_id,
              level = .level,
              geometry =
                st_sfc(
                  st_multipolygon(list(rings)),
                  crs = utm_crs
                )
            )
          }
        ) %>%
        bind_rows()
      polys %>%
        st_make_valid() %>%
        st_transform(4326)
    }

  tag_kdes <-
    field_tags %>%
    map(make_tag_kde) %>%
    bind_rows()

  if (nrow(tag_kdes) > 0) {
    st_write(
      tag_kdes,
      "data/localization/validation/tag_kdes.geojson",
      delete_dsn = TRUE,
      quiet = TRUE
    )
  }

  # ---- 7. map: true points, predicted points, and KDE polygons ---------

  ground_truth_map <-
    basemap %>%
    addCircleMarkers(
      data = calibration_pts,
      radius = 5,
      color = "#000000",
      weight = 2,
      fillColor = "#00ffff",
      fillOpacity = 0.9,
      label = ~ point_name,
      group = "True points"
    ) %>%
    addCircleMarkers(
      data =
        point_error_best %>%
        filter(!is.na(lon)),
      lng = ~ lon,
      lat = ~ lat,
      radius = 4,
      color = "#ff6b6b",
      stroke = FALSE,
      fillOpacity = 0.7,
      group = "Predicted points",
      popup =
        ~ str_c(
            tag_id, " (", curve, ")<br>error ",
            round(error_m), " m<br>conf ", round(conf_min_dist), " m"
          )
    )

  if (nrow(tag_kdes) > 0) {
    ground_truth_map <-
      ground_truth_map %>%
      addPolygons(
        data = tag_kdes,
        weight = 2,
        color = "#440154",
        fill = FALSE,
        dashArray = "4,4",
        group = "KDE polygons",
        label = ~ str_c(tag_id, " ", level, "% KDE")
      )
  }

  ground_truth_map <-
    ground_truth_map %>%
    addLayersControl(
      overlayGroups =
        c(
          "True points",
          "Predicted points",
          "KDE polygons"
        ),
      options = layersControlOptions(collapsed = FALSE)
    )

  htmlwidgets::saveWidget(
    ground_truth_map,
    "data/localization/validation/ground_truth_map.html"
  )

  message(
    str_c(
      "Localized ", sum(!is.na(point_error_best$error_m)),
      " of ", nrow(point_error_best), " point x tag fixes."
    )
  )
} else {
  # No field detections for 2026-06-26 are present in the database (the raw
  # table ends 2026-04-13 and node_health has no June-2026 rows). The
  # workflow above is validated (see groundtruth/validate_localization.py);
  # rerun once the campaign detections are imported to populate the outputs.

  message(
    "No 2026-06-26 detections found in oxbow.duckdb for the four tags -- ",
    "the field campaign has not yet been imported. Localization outputs are ",
    "empty pending that import."
  )
}
