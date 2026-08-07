
# setup -------------------------------------------------------------------

library(celltracktech)

load_dot_env(file = ".env")

my_token <- Sys.getenv("API_KEY")

# connect to Database using DuckDB ----------------------------------------

con <- 
  DBI::dbConnect(
    duckdb::duckdb(), 
    dbdir = "./data/oxbow.duckdb", 
    read_only = FALSE
  )

# load data ---------------------------------------------------------------

get_my_data(
  my_token = my_token,
  outpath = "data/compressed_files", 
  db_name = con, 
  myproject = "Grassland Fledgelings",
  begin = as_date("2022-01-01"),
  end = today(),
  filetypes = 
    c(
      "raw",
      "gps",
      "node_health", 
      "log"
    )
)

# update the database -----------------------------------------------------

# Upload the compressed data (i.e. the ‘.csv.gz’ files) into the DuckDB
# database:

update_db(
  con, 
  "data/compressed_files",
  "Grassland Fledgelings"
)

# List tables in database:

DBI::dbListTables(con)



