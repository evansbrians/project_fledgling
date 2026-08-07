
# This script is for setting up runs. It:
# - Loads the libraries
# - Reads in my custom functions
# - Reads in data files that are used
# - Establishes connections to the relevant data tables

# setup -------------------------------------------------------------------

library(celltracktech)
library(DBI)
library(duckdb)
library(duckplyr)
library(tidyverse)

load_dot_env(file = ".env")

# Custom functions:

library(oxbowR)

# Coefficients:

coefs <-
  read_rds("calibration_martin/coefs_martin.rds") %>% 
  map(
    ~ as.list(.x)
  )

# Oxbow extent:

oxbow_extent <- 
  st_read("data/oxbow_extent.geojson", quiet = TRUE) %>% 
  st_bbox() %>% 
  as.list()

# Oxbow grid:

oxbow_grid <-
  read_rds("data/oxbow_grid_20m.rds")

# Background layer for mapping:

my_tile_url <- "https://mt2.google.com/vt/lyrs=y&x={x}&y={y}&z={z}"


