# overview -----------------------------------------------------------------

# 7_field_use.R  --  step 7: KDE to fields
#
# Attribute each bird-day's space use to the farm's fields, using the daily
# KDE rasters (Stage 3) and the field polygons (data/oxbow_fields.geojson).
#
# For each tag-day we integrate the KDE utilization distribution (UD) over each
# field -- the share of UD mass in a field is our proxy for the proportion of
# that day's use spent there. We report it within two home-range isopleths:
#   - 95% : broad daily use
#   - 50% : core use
# UD mass that falls outside every field goes to an "off_field" category.
#
# confidence = M50 / M95 (a field's UD mass in the 50% core divided by its mass
# in the 95% region): near 1 = the field is central to the core (a confident
# attribution); near 0 = the field only appears in the diffuse 95% tail (likely
# transit or positional spillover). It is a field-day property, reported on both
# level rows. Off-field confidence is NA.
#
# Output (tidy): tag_id, species, date, level, field_name, prop_use, confidence
# Work is CHECKPOINTED per tag.
#
# Requires: tidyverse, sf, arrow
# (terra, exactextractr, fs, cli, glue, rlang via ::)

# setup --------------------------------------------------------------------

library(tidyverse)
library(sf)
library(arrow)

# custom functions ---------------------------------------------------------

library(oxbowR)
library(fs)
library(cli)
library(glue)

# terra, exactextractr, fs, cli, and glue are used via `::` only (not attached)
# so they cannot mask tidyverse functions (e.g. tidyr::extract).

# configuration ------------------------------------------------------------

config <-
  list(
    kde_dir = "results/localization/kde_output",
    daily_areas_path = "results/localization/kde_output/daily_areas.csv",
    fields_path = "data/oxbow_fields.geojson",
    out_dir = "results/localization/field_use",
    utm_crs = 32618,
    levels = c(95, 50),
    min_prop = 0.005
  )

set.seed(2358)

dir_create(config$out_dir)

# fields -------------------------------------------------------------------

fields <-
  config$fields_path %>%
  st_read(quiet = TRUE) %>%
  st_transform(config$utm_crs) %>%
  st_make_valid() %>%
  select(
    field_name = `Field.name`
  )

# run ----------------------------------------------------------------------

daily_areas <-
  read_csv(config$daily_areas_path, show_col_types = FALSE) %>%
  mutate(
    tag_id = as.character(tag_id),
    date = as_date(date)
  )

field_use <-
  daily_areas %>%
  distinct(tag_id) %>%
  pull(tag_id) %>%
  imap(
    \(.tag_id, .index) {
      fielduse_process_tag(
        .tag_id, .index,
        daily_areas, config$kde_dir, fields,
        config$out_dir, config$min_prop
      )
    }
  ) %>%
  bind_rows()

write_csv(
  field_use,
  glue("{config$out_dir}/field_use.csv")
)
