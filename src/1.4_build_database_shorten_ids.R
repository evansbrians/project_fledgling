
# setup -------------------------------------------------------------------

library(duckplyr)
library(tidyverse)

# Connect to database:

oxbow_db <- 
  DBI::dbConnect(
    duckdb::duckdb(), 
    dbdir = "./data/oxbow.duckdb", 
    read_only = FALSE
  )

# Read tags:

tags <- 
  tbl(oxbow_db, "deployed_tags") %>% 
  collect()

# process data ------------------------------------------------------------

tags_with_shortened_id <- 
  tags %>%
  mutate(
    tag_short_location =
      map_int(
        tag_id,
        \(.x) {
          detect_index(
            seq_len(nchar(.x)
            ),
            \(len) {
              sum(
                str_starts(
                  tags$tag_id, 
                  str_sub(.x, 1, len)
                )
              ) == 1
            } 
          )
        }
      ),
    tag_short = 
      str_sub(
        tag_id,
        end = case_when(
          str_sub(
            tag_id, 
            tag_short_location, 
            tag_short_location
          ) == "0" ~ tag_short_location + 1,
          .default = tag_short_location
        )
      ),
    .after = tag_id
  ) %>%
  select(!tag_short_location)

# end session -------------------------------------------------------------

# Write to file:

tags_with_shortened_id %>% 
  write_csv("data/tags_with_shortened_ids.csv")

# Disconnect from database:

DBI::dbDisconnect(oxbow_db)
