# For ground-truthing calibration data

# setup -------------------------------------------------------------------

library(celltracktech)
library(DBI)
library(duckdb)
library(duckplyr)
library(sf)
library(mapedit)
library(leaflet)
library(leaflet.extras)
library(tidyverse)

# Custom functions:

library(oxbowR)
library(glue)

# Connect to database:

oxbow_db <- 
  dbConnect(
    duckdb(), 
    dbdir = "./data/oxbow.duckdb", 
    read_only = FALSE
  )

# Fields:

fields <- 
  st_read("data/oxbow_fields.geojson") %>% 
  janitor::clean_names()

# Basemap for mapping:

basemap <- 
  leaflet() %>% 
  addProviderTiles(
    providers$Esri.WorldImagery,
    group = "Satellite",
    options = tileOptions(maxZoom = 21)
  ) %>%
  addProviderTiles(
    providers$OpenStreetMap,
    group = "Street Map"
  ) %>% 
  addPolygons(
    data = fields,
    color = "#fff",
    weight = 1,
    fillColor = "#00ff00",
    fillOpacity = 0.1,
    label = ~ field_name,
    group = "Fields"
  )

# define buildings --------------------------------------------------------

buildings <-
  basemap %>% 
  addDrawToolbar(
    targetGroup = "patches",
    
    # We don't want the options for drawing new shapes:
    
    polylineOptions = FALSE,
    circleOptions = FALSE,
    rectangleOptions = FALSE,
    markerOptions = FALSE,
    circleMarkerOptions = FALSE,
    
    # We do want editing options:
    
    editOptions = editToolbarOptions()
  ) %>% 
  
  # Get returns for the edited content:
  
  mapedit::editMap() %>% 
  pluck("all")

buildings <-
  read_rds("data/localization/validation/main_buildings.rds") %>% 
  bind_rows(
    read_rds("data/localization/validation/other_buildings.rds")
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
    
    time >= !!as_datetime("2026-04-01"),
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
  ) %>% 
  convert_df_to_pts()

# Disconnect from the database:

DBI::dbDisconnect(oxbow_db)

# Have a look:

basemap %>% 
  addCircleMarkers(
    data = node_locations, 
    radius = 5, 
    color = "#000000", 
    opacity = 0.9,
    weight = 2,
    fillColor = "#00ffff", 
    fillOpacity = 0.9, 
    label = ~ node_id
  )

# spacing criteria --------------------------------------------------------

# We need to ensure that the spacing replicates the spacing during Blair's 
# calibration. The median distance between nodes was 163 m, but that leaves
# too few locations to sample, so I am increasing it to 250 m.

buffered_nodes <-
  seq(150, 225, by = 25) %>% 
  as.list() %>% 
  set_names(
    str_c("d", .)
  ) %>% 
  imap(
    ~ node_locations %>% 
      st_buffer(.x) %>% 
      mutate(
        dist = .y,
        .after = node_id
      ) %>% 
      st_transform(32618)
  )

nodes_3_or_more <-
  buffered_nodes %>% 
  imap(
    ~ {
      st_intersection(.x) %>% 
        st_collection_extract("POLYGON", warn = FALSE) %>% 
        st_make_valid() %>% 
        filter(
          !st_is_empty(geometry),
          n.overlaps >= 3
        ) %>% 
        st_union() %>% 
        st_cast("MULTIPOLYGON") %>% 
        st_sf(
          dist = .y,
          geometry = .
        )
    }
  ) %>% 
  bind_rows() %>% 
  st_transform(4326) %>% 
  st_make_valid()


# subsetting --------------------------------------------------------------

# Remove nodes that are not within 250 meters of the one of the overlap zones:

nodes_dist_subset <-
  node_locations %>% 
  mutate(
    min_dist =
      st_distance(., nodes_3_or_more) %>% 
      apply(1, min) %>% 
      as.numeric()
  ) %>% 
  filter(min_dist < 250)

# Get the convex hull for these points

nodes_mch <-
  nodes_dist_subset %>% 
  st_union() %>% 
  st_convex_hull()

# Subset to regions of the polygons that are within the convex hull:

node_areas <- 
  st_intersection(
    filter(nodes_3_or_more, dist == "d200"),
    nodes_mch
  ) %>% 
  
  # Subset to regions of the polygons that are in fields:
  
  st_intersection(fields) %>% 
  
  # Subset to regions of the polygons that are 50m or more from buildings:
  
  st_difference(
    st_union(
      buildings %>% 
        st_buffer(50)
    )
  ) %>% 
  
  # Subset to regions of the polygons that are 5m or more from nodes (not a
  # fair sampling otherwise):
  
  st_difference(
    st_union(
      node_locations %>% 
        st_buffer(5)
    )
  )

# Have a look:

basemap %>% 
  addPolygons(
    data = node_areas,
    weight = 2,
    color = "#000",
    fillColor = "#0000ff"
  ) %>%
  addPolygons(
    data = buildings,
    weight = 2,
    color = "#000",
    fillColor = "#ff0000",
    fillOpacity = 0.5
  ) %>%
  addCircleMarkers(
    data = nodes_dist_subset, 
    radius = 5, 
    color = "#000000", 
    opacity = 0.9,
    weight = 2,
    fillColor = "#00ffff", 
    fillOpacity = 0.9, 
    label = ~ node_id
  )

# sampling ----------------------------------------------------------------

# Convert the node areas to UTM:

node_areas_utm <- 
  st_transform(node_areas, 32618)

# Get four points that are close to field borders (for setting down tags):

field_border_points <-
  fields %>% 
  
  # Get the field borders (in UTM):
  
  st_transform(32618) %>% 
  st_cast("POLYGON", warn = FALSE) %>% 
  st_cast("LINESTRING", warn = FALSE) %>% 
  convert_line_to_points(.density = 1) %>%
  
  # Subset to field borders that are near the acceptable node areas:
  
  st_filter(
    node_areas_utm,
    .predicate = st_is_within_distance,
    dist = 5
  ) %>% 
  
  # Reorder points so its a random draw (there are not very many points
  # available):
  
  slice_sample(
    n = nrow(.), 
    by = field_name
  ) %>% 
  
  # Reduce to 4 sampled points that are a minimum of 75 meters from one another:
  
  reduce_sampled_points(
    .distance_threshold = 75,
    .n_points = 4
  ) %>% 
  
  # Add a point id for mapping:
  
  mutate(
    point_id = glue("fence_{row_number()}")
  )
  
# Get 12 inner field points:

inner_field_points <-
  node_areas_utm %>% 
  
  # Grab a big random sample of points:
  
  st_sample(
    size = 10000, 
    type = "random"
  ) %>%
  st_sf() %>%
  
  # Subset to where they are at least 75 meters from the fence points:
  
  filter(
    st_distance(
      geometry, st_union(field_border_points)
    ) %>% 
      as.numeric() %>% 
      { . > 75 }
  ) %>% 
  
  # Reduce to 12 random points:
  
  reduce_sampled_points(
    .distance_threshold = 75,
    .n_points = 12
  ) %>% 
  
  # Add a point id for mapping:
  
  mutate(
    point_id = glue("field_{row_number()}")
  )

# Point object is a combination of the fence and field points:

rando_pts <-
  list(
    field_border_points,
    inner_field_points
  ) %>%
  map(
    ~ st_transform(.x, 4326)
  ) %>% 
  bind_rows() %>% 
  mutate(
    point_class = 
      point_id %>% 
      str_extract("[a-z]*") %>% 
      factor()
  )

# Have a look -------------------------------------------------------------

# Randomly sample 

basemap %>% 
  addPolygons(
    data = node_areas_no_buildings,
    weight = 2,
    color = "#000",
    fillColor = "#0000ff",
    group = "Sampling area"
  ) %>%
  addCircleMarkers(
    data = nodes_dist_subset, 
    radius = 5, 
    color = "#000000", 
    opacity = 0.9,
    weight = 2,
    fillColor = "#00ffff", 
    fillOpacity = 0.9, 
    label = ~ glue("node: {node_id}"),
    group = "Nodes"
  ) %>% 
  addCircleMarkers(
    data = rando_pts,
    radius = 4,
    color = "#000",
    weight = 2,
    opacity = 1,
    fillColor = 
      ~ colorFactor(
          palette = c("#f00","#00f"),
          domain = point_class
          )(point_class),
    fillOpacity = 0.9,
    group = "Sampling points"
  ) %>%
  
  # Layer control:
  
  addLayersControl(
    baseGroups = c("Satellite", "Street Map"),
    overlayGroups = 
      c(
        "Sampling area",
        "Nodes",
        "Sampling points",
        "Fields"
      ),
    options = layersControlOptions(collapsed = TRUE)
  )

# write files for leaflet map app -----------------------------------------


nodes_dist_subset %>% 
  st_write("data/localization/validation/proc/nodes.geojson")

node_areas %>% 
  st_write("data/localization/validation/proc/sample_areas.geojson")

rando_pts %>% 
  st_write("data/localization/validation/proc/rando_pts.geojson")


