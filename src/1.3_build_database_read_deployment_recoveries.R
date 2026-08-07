# metadata ----------------------------------------------------------------

# Title: read_deployment_recoveries
# Author: Brian
# Date created: 2024-05-29
# Referenced by: scripts/update_database.R
# References: NA

# Description: This script reads in and cleans the deployment and recovery 
# spreadsheet from Google Sheets.

# Usage note: To execute this code, you may be required to authenticate 
# your user credentials (especially the first time that you run it).

# setup -------------------------------------------------------------------

library(duckplyr)
library(tidyverse)

# Get custom functions used:

library(oxbowR)

# Connect to database:

oxbow_db <- 
  DBI::dbConnect(
    duckdb::duckdb(), 
    dbdir = "./data/oxbow.duckdb", 
    read_only = FALSE
  )

# get deployments and recoveries ------------------------------------------

tag_deployment_recovery <- 
  
  ## connect to file ------------------------------------------------------

# Create path to online file:

file.path(
  the_goog = "https://docs.google.com/spreadsheets/d",
  sheet_id = "1RG3NORA1mhS3oX6mOmO4Tssd-V0jHf1MvfL606yAyHg"
) %>% 
  
  # Read in google sheet:
  
  googlesheets4::read_sheet() %>% 
  
  ## pre-process data -----------------------------------------------------

mutate(
  tag_id = as.character(tag_id),
  species,
  
  # Define deployment time:
  
  deploy_time = 
    
    # Repair times and fill in missing time values:
    
    time %>% 
    str_sub(-5) %>% 
    str_replace("-", "0") %>% 
    replace_na("06:00") %>% 
    
    # Combine deployment dates and times:
    
    str_c(date, ., sep = " ") %>% 
    
    # Convert character datetime to a true datetime with site tz:
    
    as_datetime(
      format = "%Y-%m-%d %H:%M",
      tz = "America/New_York"
    ) %>% 
    
    # Convert tz to UTC (tz of CTT detections):
    
    force_tzs(
      tzones = "America/New_York",
      tzone_out = "UTC",
    ),
  
  # Define recovery time:
  
  recover_time =
    recovered_on %>% 
    
    # Fill in missing time values:
    
    if_else(
      nchar(.) == 10,
      str_c(., " 23:59"),
      .
    ) %>% 
    
    # Convert character datetime to a true datetime with site tz:
    
    as_datetime(
      format = "%Y-%m-%d %H:%M",
      tz = "America/New_York"
    ) %>% 
    
    # Convert tz to UTC (tz of CTT detections):
    
    force_tzs(
      tzones = "America/New_York",
      tzone_out = "UTC",
    ),
  .keep = "none"
)

# add to database ---------------------------------------------------------

DBI::dbWriteTable(
  oxbow_db,
  "deployed_tags",
  tag_deployment_recovery
)

DBI::dbDisconnect(oxbow_db)

