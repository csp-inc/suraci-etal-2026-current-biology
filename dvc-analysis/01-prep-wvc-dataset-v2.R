## ---------------------------
##
## Script name: 01-prep-wvc-dataset.R
##
## Purpose of script: A script for cleaning and 
## processing the Washington Department of Transportation 
## Ungulate Collision data and preparing for analysis 
##
## Author: Justin Suraci
##
## ---------------------------

library(readxl)
library(janitor)
library(tidyverse)
library(lubridate)
library(sf)
library(suncalc)

# ** EXPORT SPATIAL PRODUCTS? **
makeExports = FALSE

# Load relevant datasets -------------------------------------------------
crs <- st_crs(5070)
aoi <- st_read("data/puma-aoi.geojson") # Olympic peninsula AOI
load("data/cams-w-puma-deer-occ.rda", verbose = T) # cams w occ data
cams <- cams %>% 
  st_transform(crs)

#-------------------------------------------------------------------------------
# PREP EVENT-BASED WVC DATASET--------------------------------------------------
#-------------------------------------------------------------------------------
# Load collision data and convert to spatial ------------------------------
wvc <- read_csv("data/20250827Suraci_All_Wildlife_Crashes_Multi_Counties.csv", skip = 1) %>% 
  clean_names() %>%
  st_as_sf(., coords=c("wa_state_plane_south_x", "wa_state_plane_south_y"), crs=2286, na.fail=F) %>%
  filter(!st_is_empty(.)) %>% 
  st_transform(crs) %>% 
  mutate(date_time = as.POSIXct(strptime(paste(date_ymd, x24_hr_time), "%Y-%m-%d %H:%M:%S")),
         time_of_day = case_when(lighting_conditions == "Daylight" ~ "Daylight",
                                    lighting_conditions %in% c("Dawn", "Dusk") ~ "Crepuscular",
                                    lighting_conditions %in% c("Dark-No Street Lights", "Dark-Street Lights On",
                                                              "Dark-Street Lights Off", "Dark - Unknown Lightin") ~ "Dark",
                                    lighting_conditions %in% c("Unknown", "Other") ~ NA),
         road_condition = case_when(road_surface_conditions == "Dry" ~ "Dry",
                                    road_surface_conditions %in% c("Unknown", "Other") ~ NA,
                                    road_surface_conditions != "Dry" &
                                      road_surface_conditions != "Unknown" ~ "Wet-Frozen"),
         species = case_when(first_collision_type_object_struck == "Vehicle Strikes Deer" ~ "Deer",
                             first_collision_type_object_struck == "Vehicle Strikes Elk" ~ "Elk",
                             first_collision_type_object_struck == "Vehicle Strikes All Other Non-Domestic Animal" ~ "Other"),
         injured_dead = total_fatalities + 
           total_serious_injuries + 
           total_minor_injuries + 
           total_possible_injuries
         )

# Do some sun calcs with suncalc ------------------------------------------
wvcLL <- wvc %>% st_transform(4326)
coords <- st_coordinates(wvcLL) %>% as.data.frame()
wvcLL$lat <- coords$Y
wvcLL$lon <- coords$X
wvcLL <- wvcLL %>% select(date_time, lat, lon) %>% 
  rename(date = date_time)
sunPos <- getSunlightPosition(data = wvcLL)
moon <- getMoonIllumination(wvcLL$date)
moonPos <- getMoonPosition(data = wvcLL)
sunIllum <- ifelse(sunPos$altitude > 0, cos(pi/2 - sunPos$altitude), 0)
moonIllum <- ifelse(moonPos$altitude > 0, moon$fraction * cos(pi/2 - moonPos$altitude), 0)
wvc$sunIllum <- sunIllum
wvc$moonIllum <- moonIllum

# Calculate some average times for different diel periods
wvcLL2 <- wvcLL %>% mutate(date = as.Date(date))
sunTimes <- getSunlightTimes(data = wvcLL2)
ave_sunrise <- mean(hour(with_tz(sunTimes$sunrise, "America/Los_Angeles")))
ave_sunset <- mean(hour(with_tz(sunTimes$sunset, "America/Los_Angeles")))
ave_dawn <- mean(hour(with_tz(sunTimes$dawn, "America/Los_Angeles")))
ave_dusk <- mean(hour(with_tz(sunTimes$dusk, "America/Los_Angeles")))

# Remove extraneous columns and assign seasons ------------------------------
keeps <- c("date_ymd", "time", "date_time", "time_of_day", "sunIllum", "moonIllum", 
           "species", "veh_1_type", "most_severe_injury_type", "injured_dead",
           "total_fatalities", "total_serious_injuries", "total_minor_injuries",                                                             
           "total_possible_injuries", "weather", "road_condition", "lighting_conditions")
wvc <- wvc %>%
  dplyr::select(any_of(keeps)) %>%
  dplyr::mutate(month = month(date_time)) %>%
  dplyr::mutate(season = 
                  case_when(month %in% c(12, 1, 2) ~ "winter",
                            month %in% c(3, 4, 5) ~ "spring",
                            month %in% c(6, 7, 8) ~ "summer",
                            month %in% c(9, 10, 11) ~ "fall"))

# Add nearest neighbor camera trap ID and distance ------------------------------
nearest_idxs <- st_nearest_feature(wvc, cams)
nearest_cam <- cams[nearest_idxs, ]
nearest_cam_name <- nearest_cam$CamID
nn_cam_dist <- st_distance(wvc, nearest_cam, by_element = TRUE) %>% as.numeric()
wvc <- wvc %>% mutate(CamID = nearest_cam_name,
                      nn_cam_dist = nn_cam_dist)
# Add some camera attributes
camCovs <- cams %>% select(CamID, "deerOccProb", "pumaOccProb",
                           "probUnoccupied", "deerOnlyProb",
                           "deerPumaProb", "killsite_mod_1km",
                           "killsite_mod_100m", "movement_mod_1km",
                           "movement_mod_100m" ) %>% 
  st_drop_geometry()
wvc <- wvc %>% left_join(camCovs, by = "CamID")

#-------------------------------------------------------------------------------
# PREP GRID-BASED WVC COUNT DATASET---------------------------------------------
#-------------------------------------------------------------------------------
# COMPILE ROAD DATA
# Define road types to include 
rd_types = c("primary", "primary_link", "secondary", "secondary_link", 
             "motorway", "motorway_link", "trunk", "trunk_link", 
             "tertiary", "tertiary_link", "residential", 
             "living_street", "busway")

# Filter OSM roads data (download for WA from https://download.geofabrik.de/north-america/us.html)
wa_rds <- st_read("data/washington-osm-roads/gis_osm_roads_free_1.shp") %>% 
  dplyr::filter(fclass %in% rd_types)
aoiRoads <- st_intersection(wa_rds, aoi) %>% 
  st_transform(st_crs(wvc))

#Export shapefile if desired...
if(makeExports){
  aoiRoads %>% dplyr::select(fclass, name, osm_id) %>%
    st_write(dsn = "data/op-roads", layer = "op-roads",
             driver = "ESRI Shapefile", delete_layer = TRUE)
}

#----------------------------------
# CREATE GRID FOR ANALYSIS
# Define grid cell size
# calculate average distance between cameras as potential filter/cut-off for nn_cam_dist
camDist <- st_distance(cams)
diag(camDist) <- Inf
nnDist <- apply(camDist, 1, min)
mean(nnDist)

# Project AOI and make grid across OP
aoiProj <- aoi %>% st_transform(st_crs(wvc))
opGrid <-  aoiProj %>% 
  st_make_grid(cellsize = 2250,     
               square = TRUE,
               what = "polygons")
opGrid <- st_intersection(st_sf(geometry = opGrid), aoiProj) %>% 
  mutate(cellID = row_number(),
         cellArea = as.numeric(st_area(.))/1e6) %>% 
  select(cellID, cellArea)

#----------------------------------
# ADD ROAD DENSITY
# Intersect and assign grid IDs
lines_split <- st_intersection(aoiRoads, opGrid) %>%
  mutate(length_m = st_length(geometry)) %>%
  st_set_geometry(NULL)

# Summarize by grid cell
lengths_by_cell <- lines_split %>%
  group_by(cellID) %>%
  summarize(road_length_m = sum(as.numeric(length_m)))

# Join to grid
opGrid <- opGrid %>% left_join(lengths_by_cell)
opGrid$road_length_m[is.na(opGrid$road_length_m)] <- 0

# calc road density
opGrid <- opGrid %>% 
  mutate(road_length_km = road_length_m/1000,
         road_density = road_length_km/cellArea)

ggplot() +
  # Polygon layer with color ramp
  geom_sf(data = aoiRoads, color = "black", size = 0.4) +
  geom_sf(data = opGrid, aes(fill = road_density), color = NA) +
  scale_fill_viridis_c(option = "plasma", alpha = 0.8) +  # You can choose from "magma", "viridis", etc.
  theme_minimal() + 
  labs(title = "Road density over roads")

#----------------------------------
# ADD WVC COUNT
secIn5Yrs <- 157788000
fiveYrsAgo <- max(wvc$date_time)-secIn5Yrs
tenYrsAgo <- max(wvc$date_time)-(secIn5Yrs*2)
wvc_cell <- st_join(wvc, opGrid[, "cellID"], left = FALSE) %>% 
  st_drop_geometry() %>%
  group_by(cellID) %>% 
  summarise(wvc_count_ungulate = length(which(species == "Deer" | species == "Elk")),
            wvc_count_ungulate_5yrs = length(which(species %in% c("Deer", "Elk") & date_time >= fiveYrsAgo)),
            wvc_count_ungulate_10yrs = length(which(species %in% c("Deer", "Elk") & date_time >= tenYrsAgo)),
            wvc_count_deer = length(which(species == "Deer")),
            wvc_count_deer_5yrs = length(which(species == "Deer" & date_time >= fiveYrsAgo)),
            wvc_count_deer_10yrs = length(which(species == "Deer" & date_time >= tenYrsAgo)),
            n_injured_dead = sum(injured_dead),
            n_s_injured_dead = sum(total_fatalities + total_serious_injuries),
            injured_dead_5yrs = sum(injured_dead[date_time >= fiveYrsAgo]),
            injured_dead_10yrs = sum(injured_dead[date_time >= tenYrsAgo]))

# Join to grid
opGrid <- opGrid %>% left_join(wvc_cell)
opGrid[is.na(opGrid)] <- 0

ggplot() +
  geom_sf(data = wvc, color = "black", size = 0.4) +
  geom_sf(data = opGrid, aes(fill = n_injured_dead), color = NA) +
  scale_fill_viridis_c(option = "plasma",alpha = 0.7) +  
  theme_minimal() + 
  labs(title = "WVC count over incident locations")

# Create export versions for cov extraction
if(makeExports){
  save(opGrid, file = "data/olympic-peninsula-grid/op-rd-grid.rda")
  exportGrid <- opGrid %>% filter(road_density > 0) %>% 
    select(cellID, cellArea, road_length_km, road_density,
           wvc_count_ungulate, wvc_count_ungulate_5yrs, 
           wvc_count_ungulate_10yrs, wvc_count_deer, wvc_count_deer_5yrs,
           wvc_count_deer_10yrs) %>% 
    rename(rd_km = road_length_km,
           rd_dens = road_density,
           wvcUng = wvc_count_ungulate,
           wvcUng5y = wvc_count_ungulate_5yrs,
           wvcUng10y = wvc_count_ungulate_10yrs,
           wvcDeer = wvc_count_deer,
           wvcDeer5y = wvc_count_deer_5yrs,
           wvcDeer10y = wvc_count_deer_10yrs)
  gridCent <- st_centroid(exportGrid)
  st_write(exportGrid, dsn = "data/olympic-peninsula-grid",
           layer = "op-rd-grid", driver = 'ESRI Shapefile', delete_layer = TRUE)
  st_write(gridCent, dsn = "data/olympic-peninsula-grid",
           layer = "op-rd-centroid", driver = 'ESRI Shapefile', delete_layer = TRUE)
}

#----------------------------------
# ADD CAMERA INFO
cam_in_cell <- st_join(cams, opGrid[, "cellID"], left = FALSE) %>% 
  dplyr::select(CamID, cellID, deerOccProb, pumaOccProb, 
                deerOnlyProb, deerPumaProb,
                probUnoccupied, killsite_mod_100m, 
                killsite_mod_1km, movement_mod_100m,
                movement_mod_1km, developed_pcov_100m,
                developed_pcov_1km) %>% 
  st_drop_geometry() %>% 
  group_by(cellID) %>% 
  summarize(deerOcc = mean(deerOccProb),
            pumaOcc = mean(pumaOccProb),
            unOcc = mean(probUnoccupied),
            deerOnly = mean(deerOnlyProb),
            deerPuma = mean(deerPumaProb),
            pumaBinary = ifelse(pumaOcc > 0.75, 1, 0),
            deerBinary = ifelse(deerOcc > 0.75, 1, 0),
            killMod100m = mean(killsite_mod_100m),
            killMod1km = mean(killsite_mod_1km),
            moveMod100m = mean(movement_mod_100m),
            moveMod1km = mean(movement_mod_1km))

# Join to grid
opGridJoin <- opGrid %>% left_join(cam_in_cell)

# Subset to just relevant cells
wvcGrid <- opGridJoin %>% filter(!is.na(deerOcc), road_density > 0)
# Add coordinates of centroids (for autocorrelation check)
wvcGrid$x <- st_coordinates(st_centroid(wvcGrid))[,1]
wvcGrid$y <- st_coordinates(st_centroid(wvcGrid))[,2]

# Add proportion of each cell developed (based on 90-m pixels w/in cell)
covs <- rast("data/ssf-cov-rasters/puma-ssf-focal-covs-90m.tif")
gtrans <- st_transform(wvcGrid, st_crs(covs)) %>% vect()
devPx <- terra::extract(covs$developed, gtrans, fun = sum, na.rm = TRUE)
numPx <- terra::extract(covs$developed, gtrans, 
                        fun = function(x)length(which(!is.na(x))))
devProp <- devPx[,2]/numPx[,2]
wvcGrid$devProp <- devProp

#----------------------------------
# ADD TRAFFIC VOLUME (LOTS OF MISSING DATA...)
# Get traffic data
traffic <- st_read("data/WSDOT-HistoricTrafficData_2022/Traffic_Sections_2022.shp") %>% 
  st_transform(5070) %>% 
  select(AADT)
# Join traffic volume data from 2022 to wvcGrid 
# ** Results in duplicate cells **
traffGrid <- st_join(wvcGrid, traffic)

# Deal with duplicate cells
cellID <- wvcGrid$cellID
aadt_temp <- numeric(length(cellID))
for(i in 1:length(cellID)){
  temp <- traffGrid[which(traffGrid$cellID==cellID[i]),]
  if(!is.na(temp$AADT[1])) tAADT = mean(temp$AADT) else tAADT = NA
  aadt_temp[i] <- tAADT
}

#Add to main dataset
wvcGrid$AADT2022 <- aadt_temp

#-------------------------------------------------------------------------------
# SAVE & WRITE------------------------------------------------------------------
#-------------------------------------------------------------------------------
save(wvc, wvcGrid, 
     file = "data/wvc-analysis-dataset-20251109.rda")
st_write(wvc,
         "data/washington_DOT_collisions_cleaned.gpkg",
         delete_layer = TRUE)


