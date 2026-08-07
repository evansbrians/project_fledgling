# This script generates the data frame `detections_with_distance`, which 
# describes the following variables during Blair's 27 Apr 2023 transect walk:

# * time: The date-time (ISO 8601 datetime object) of a given detection
# * tag_id: The unique identifier of the tag that was detected
# * tag_type: Whether the tag is meant for EAME or BOBO
# * height_m: The height from the tag to the ground, in meters
# * node_id: The unique identifier of the node where the tag was detected
# * node_dist: The distance between the tag and the node, in meters
# * rssi: The radio signal strength index of the detection

# setup -------------------------------------------------------------------

library(tidyverse)
library(sf)

# data starts -------------------------------------------------------------

# GPS data:

bad_elf <- 
  read_csv('data/blair/2023-04-27T15-28-30Z/track_points.csv') %>% 
  select(
    time,
    lon = X,
    lat = Y,
    elevation = ele) %>% 
  mutate(time = as_datetime(time)) %>% 
  arrange(time)

# Transmitters used in this study:

transmitters <-
  readxl::read_excel('data/blair/blairEtransmitternumbers.xlsx') %>% 
  pivot_longer(
    tag_id_eame:tag_id_bobo,
    names_to = "tag_type",
    names_prefix = "tag_id_",
    values_to = "tag_id")

# A node health file from the day of sampling:

node_locations <- 
  read_csv('data/blair/CTT-129D85D701A7-node-health.2023-04-27_145549.csv') %>% 
  group_by(node_id = NodeId) %>% 
  filter(NodeRSSI == max(NodeRSSI)) %>% 
  rename(lon = Longitude, lat = Latitude) %>% 
  summarize(
    across(
      c(lon, lat),
      ~ mean(.x, na.rm = TRUE)))

# Detections:

detections <- 
  read_rds('data/blair/detections_2023_apr.rds') %>% 
  
  # Subset to detections made during the transect walk:
  
  filter(
    between(
      time, 
      min(bad_elf$time), 
      max(bad_elf$time))) %>% 
  
  # Subset to transmitters that were a part of this study and add height and
  # tag type variables:
  
  inner_join(transmitters, by = "tag_id")

# Remove any files that will not be needed in subsequent steps:

rm(transmitters)

# distance between bad elf and nodes --------------------------------------

# Nodes and bad elf, spatial:

spatial_node_and_elf <-
  list(node_locations, bad_elf) %>% 
  map(
    ~ st_as_sf(
      .x,
      coords = c("lon", "lat"),
      crs = 4326)) %>% 
  set_names(
    c("node_locations", "bad_elf"))

# Calculate the distance between each bad elf GPS point and each node:

node_distances <- 
  map_dfr(
    1:nrow(node_locations),
    function(i) {
      
      # Grab a node:
      
      focal_node <-
        spatial_node_and_elf %>% 
        pluck("node_locations") %>% 
        slice(i)
      
      # Calculate the distance to the node and return frame:
      
      bad_elf_dist <-
        spatial_node_and_elf %>% 
        pluck("bad_elf") %>% 
        mutate(
          node_id = focal_node$node_id,
          node_dist = 
            st_distance(
              ., 
              focal_node,
              by_element = TRUE)) %>% 
        as_tibble() %>% 
        select(time, elevation, node_id, node_dist)
    })

# Remove any files that will not be needed in subsequent steps:

rm(bad_elf,
   node_locations, 
   spatial_node_and_elf)

# join detections and node distances --------------------------------------

detections_with_distance <- 
  detections %>% 
  left_join(
    node_distances, 
    by = c("node_id", "time")) %>% 
  mutate(node_dist = as.numeric(node_dist)) %>% 
  arrange(time, height_m, tag_type) %>% 
  select(time,
         tag_id, 
         tag_type,
         height_m, 
         node_id,
         node_dist,
         rssi)
  

# Remove any files that will not be needed in subsequent steps:

rm(detections, node_distances)

# just some examples ------------------------------------------------------

# How does RSSI differ by tag type?

detections_with_distance %>% 
  filter(
    between(node_dist, 75, 100)) %>% 
  ggplot(
    aes(x = tag_type, y = rssi)) +
  geom_boxplot(fill = "#dcdcdc") +
  theme_classic()

# How does RSSI differ by tag height?

detections_with_distance %>% 
  filter(
    between(node_dist, 75, 100)) %>% 
  mutate(height_m = factor(height_m)) %>% 
  ggplot(
    aes(x = height_m, y = rssi)) +
  geom_boxplot(fill = "#dcdcdc") +
  theme_classic()

# How does RSSI differ by tag height and tag type?

detections_with_distance %>% 
  filter(
    between(node_dist, 75, 100)) %>% 
  mutate(height_m = factor(height_m)) %>% 
  ggplot(
    aes(x = height_m, 
        y = rssi,
        fill = tag_type)) +
  geom_boxplot() +
  theme_classic()


