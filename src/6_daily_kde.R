# overview -----------------------------------------------------------------

# 6_daily_kde.R  --  step 6: KDE
#
# Stage 3 of the localization workflow: daily kernel density estimates (KDE)
# of space use for each tag -- both raster surfaces and home-range polygons.
#
# Input is the Stage 2 output (localizations.parquet): one row per localized
# window, carrying lon / lat and a `weight` (1 / expected_err^2). For each tag,
# fixes are grouped by day; every bird-day with at least `kde_min_fix` fixes
# gets a weighted KDE. From each day's KDE we write:
#   - a raster (the utilization-distribution surface) as a GeoTIFF, and
#   - 50/80/95% home-range polygons with their areas (hectares).
# The KDE is computed once per day in UTM 18N; the raster is written in UTM
# (native, undistorted for analysis) and the polygons are returned in lon/lat.
#
# Work is CHECKPOINTED per tag; a master areas table is assembled at the end.
#
# Outputs (in `out_dir`):
#   rasters/<tag>/<tag>_<date>.tif  -- daily KDE surface (UTM 18N)
#   kde_<tag>.geojson               -- daily polygons (tag_id, species, date,
#                                      level, n_fixes, area_ha)
#   daily_areas.csv                 -- the polygons' attributes for all tags
#
# Requires: tidyverse, arrow, sf, ks, lubridate
# (terra, fs, cli, glue, rlang via ::)

# setup --------------------------------------------------------------------

library(tidyverse)
library(arrow)
library(sf)
library(ks)

# custom functions ---------------------------------------------------------

library(oxbowR)
library(fs)
library(cli)
library(glue)

# terra is used via `terra::` prefixes only (NOT attached): attaching it masks
# several tidyverse/base functions (e.g. tidyr::extract, and terra's own
# extract / intersect / union / mask). Namespacing keeps tidyverse intact.

# configuration ------------------------------------------------------------

config <-
  list(
    localizations_path = "data/localization/localizations.parquet",
    out_dir = "results/localization/kde_output",
    utm_crs = 32618,
    kde_levels = c(50, 80, 95),
    kde_min_fix = 10
  )

set.seed(2358)

dir_create(config$out_dir)

# run ----------------------------------------------------------------------

localizations_dataset <- arrow::open_dataset(config$localizations_path)

tags <-
  localizations_dataset %>%
  distinct(tag_id, species) %>%
  collect()

daily_areas <-
  tags %>%
  pmap(
    \(tag_id, species) {
      kde_process_tag(
        tag_id, species,
        localizations_dataset, config$out_dir,
        config$kde_min_fix, config$utm_crs, config$kde_levels
      )
    }
  ) %>%
  bind_rows() %>%
  mutate(tag_id = as.character(tag_id))

write_csv(
  daily_areas,
  glue("{config$out_dir}/daily_areas.csv")
)
