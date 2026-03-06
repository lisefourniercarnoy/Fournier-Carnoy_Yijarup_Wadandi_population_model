# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Claude's SW habitat predictions
# Task:    Create a shapefile with grid cells
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    November 2025

# -----------------------------------------------------------------------------

# Status: Latest run: 22/01/2026

# -----------------------------------------------------------------------------

rm(list = ls()) # Clear working environment

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))

library(tidyverse) # for manipulating data
library(ggplot2) # for plots
library(sf) # for manipulating shapefiles
library(terra) # for manipulating shapefiles
library(sp) # for manipulating shapefiles
library(nngeo) # to remove rings in shapefiles
library(raster) # for manipulating rasters
library(RColorBrewer) # for colours on plots


## Files used in this script --------------------------------------------------

file_land         <- "data/input_data/Q_aus_land_high_res_no_estuary.shp"
file_MPA          <- "data/input_data/western-australia_marine-parks-all.shp"
file_MPA_fixed    <- "data/output_data/Q_NTZ_manual_clean_01_B.shp"
file_habitats     <- "data/input_data/wadandi_predicted_habitat.RDS" # Claude's SWC habitat predictions
file_shore_hab    <- "data/output_data/Q_manual_shoreline_habitat.shp" # manually categorised shore habitat
file_bathy        <- "data/input_data/AusBathyTopo__Australia__2024_250m_MSL_cog.tif"

## Set the extent -------------------------------------------------------------

common_crs = 7850 # Using GDA2020 / MGA Zone 50, suitable for WA. projection in metres to allow buffer distance to be calculated.

bbox_wadandi <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -33), crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(common_crs)

bbox_north <- st_bbox(c(xmin = 114.4, ymin = -33, xmax = 116.0, ymax = -31), crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(common_crs)

bbox_whole <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -31), crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(common_crs)

## Map of WA coastline --------------------------------------------------------

wa_map <- st_read(file_land) %>% # takes up to 30 sec. because very detailed file
  st_transform(common_crs) %>% 
  st_make_valid() %>% # fix broken geometries
  st_crop(bbox_whole) # crop to site

plot(wa_map$geometry, col = "#eeeeeeff") # Map of WA, cropped to Wadandi Country

wa_map <- wa_map[!st_is_empty(wa_map),]
st_write(wa_map, "data/output_data/01_B_land.shp", append = F)


## Map of the State + Commonwealth Marine Park --------------------------------

MP <- st_read(file_MPA) %>%
  st_transform(common_crs) %>%
  st_make_valid() %>% # fix broken geometries
  st_crop(bbox_whole) %>% 
  filter(zone %in% c("Sanctuary Zone", 
                     "National Park Zone", 
                     "Special Purpose Zone",
                     "Recreational Use Zone"
                     )
         )
# zone names in the dataframe above are inconsistent, and some zones are open to some fleets but some others. 
# you will need to manually tick which polygons are open to fishing and which aren't
write.csv(MP %>% st_drop_geometry(), "data/output_data/01_B_region_MPAs_to_classify.csv")

MP <- tibble::rowid_to_column(MP, "X")

MP_permissions <- read.csv("data/input_data/YIJARUP - zoning fishing permissions - Sheet1.csv") %>% 
  glimpse()
# here, i have every zoning polygon, and whether it's fishable to each fleet or not.
MP <- MP %>%
  left_join(MP_permissions %>% dplyr::select(X, commercial, boat_rec, shore_rec, source), by = "X")
plot(MP[,c("commercial", "boat_rec", "shore_rec")])

# we'll remove zones that are fishable to all fleets
fleets <- c("commercial", "boat_rec", "shore_rec")
MP <- MP %>% dplyr::filter(if_all(all_of(fleets)) != TRUE) # if all fishing permissions are TRUE, remove that zone

NTZ <- MP %>%
  dplyr::select(!c(colour, source)) %>% 
  st_crop(bbox_whole)
plot(NTZ$geometry)

# Save the No-Take Zones layer
plot(NTZ$geometry, col = colour_palette[4]); plot(wa_map$geometry, col = "#eeeeeeff", add = T)
st_write(NTZ, "data/output_data/01_B_NTZ_to_clean.shp", delete_layer = T)


# after manual fix (in QGIS, with 'Vector Overlay > Difference' with the 01_X_wadandi_NTZ_manually_modified_to_reach_land.shp)
NTZ_clean <- st_read(file_MPA_fixed); plot(NTZ_clean$geometry)


### Habitat Files - RUN AGAIN IF THINGS CHANGE --------------------------------

sand <- readRDS("data/output_data/01_A_bathymetry_habitat_rasters.rds")[["psand.fit"]]
reef <- readRDS("data/output_data/01_A_bathymetry_habitat_rasters.rds")[["preef.fit"]]
seagrass <- readRDS("data/output_data/01_A_bathymetry_habitat_rasters.rds")[["pseagrass.fit"]]

hab <- c(reef, sand, seagrass)
hab <- project(hab, crs(wa_map))
hab <- mask(hab, wa_map, inverse = TRUE); names(hab) <- c("reef", "sand", "seagrass")
plot(hab)

#  Crop the raster
bbox_sf <- st_transform(bbox_whole, crs(hab))
bbox_vect <- vect(bbox_sf)
cropped_ras <- crop(hab, bbox_vect)
plot(cropped_ras, range = c(0, 1))

# make a polygon to use to make grids
hab_polygon <- st_as_sf(as.polygons(app(cropped_ras[[1]], fun = function(x) ifelse(is.na(x), NA, 1)), dissolve = TRUE)) %>% 
  st_transform(common_crs); plot(hab_polygon, col = NA, border = "red")


## Create files to fill shore habitat in QGIS ---------------------------------

# The habitat predictions above don't reach the shore.
# Bathymetry files never reach the shore so no models can be fitted on this area.
# I'll export the strip between shore and habitat predictions to manually fill in QGIS
# I'm using the categorisation beach/beach+rocky/rocky, as in the management plan, page 46, available here: https://www.dbca.wa.gov.au/management/plans/ngari-capes-marine-park

# plot(wa_map$geometry)
# plot(hab_polygon$geometry)
# 
# wa_union <- st_union(wa_map)
# hab_union <- st_union(hab_polygon) %>% st_transform(st_crs(wa_union))
# 
# combined <- st_union(wa_union, hab_union); plot(combined) # area to exclude
# 
# shore_hab <- st_difference(bbox_whole, combined)
# plot(shore_hab, col = "red", main = "Area outside both polygons")
# 
# st_write(shore_hab, "data/output_data/01_B_wadandi_shore_habitat_to_categorise.shp", delete_layer = T)

shore_hab <- st_read(file_shore_hab)
plot(shore_hab)
wa_map_union <- st_union(wa_map)
shore_hab <- st_difference(shore_hab, wa_map_union) %>% st_transform(st_crs(hab)) # crop out the areas that are too high res.
plot(shore_hab)
# merge Claude's predicted habitat to the shore habitat

# sand
sand_shore <- rasterize(shore_hab[shore_hab$habitat == "sand",], hab, field = global(hab[["sand"]], fun = "max", na.rm = TRUE)[1,1], touches = TRUE)
sand_layer <- hab[["sand"]]
sand_combined <- cover(sand_layer, sand_shore)
plot(sand_combined)

# reef
reef_shore <- rasterize(shore_hab[shore_hab$habitat == "reef",], hab, field = global(hab[["reef"]], fun = "max", na.rm = TRUE)[1,1], touches = TRUE)
reef_layer <- hab[["reef"]]
reef_combined <- cover(reef_layer, reef_shore)
plot(reef_combined)

# seagrass (use the shore hab of sand and replace with zero, since we have no seagrass on the shore)
seagrass_shore <- rasterize(shore_hab[shore_hab$habitat == "sand",], hab, field = global(hab[["sand"]], fun = "max", na.rm = TRUE)[1,1], touches = TRUE)
plot(seagrass_shore)
seagrass_shore[[1]][!is.na(seagrass_shore[[1]]), ] <- 0
seagrass_layer <- hab[["seagrass"]]
seagrass_combined <- cover(seagrass_layer, seagrass_shore)
plot(seagrass_combined)

# the areas in sand_combined that are rocky are set to zero and vice versa
sand_mask <- !is.na(sand_shore); plot(sand_mask)
reef_combined_clean <- mask(reef_combined, sand_mask, maskvalues = TRUE, updatevalue = 0); plot(reef_combined_clean)

# combine into one
hab_raster <- c(sand_combined, reef_combined_clean, seagrass_combined)
hab_raster <- crop(hab_raster, st_transform(bbox_whole, crs(hab_raster))); plot(hab_raster)


## Make grid cells ------------------------------------------------------------

# make a polygon to use to make grids
# hab_polygon <- st_as_sf(as.polygons(app(hab_raster[[1]], fun = function(x) ifelse(is.na(x), NA, 1)), dissolve = TRUE)) %>% 
#   st_transform(common_crs); plot(hab_polygon, col = NA, border = "red", lwd = 2)

bathy <- rast("data/input_data/AusBathyTopo__Australia__2024_250m_MSL_cog.tif") %>% crop(st_transform(bbox_whole, 4326)); plot(bathy)
bathy <- ifel(bathy$AusBathyTopo__Australia__2024_250m_MSL_cog < -200, NA, bathy$AusBathyTopo__Australia__2024_250m_MSL_cog); plot(bathy)
temp <- st_as_sf(as.polygons(app(bathy[[1]], fun = function(x) ifelse(is.na(x), NA, 1)), dissolve = TRUE)) %>% 
  st_transform(common_crs); plot(temp, col = NA, border = "red", lwd = 2)


### 1. Shore cells ------------------------------------------------------------

# create a 500m buffer around land (takes a while, so use the already-made one)
# distance <- 100 # Buffer distance in meters
# land_buffer <- st_buffer(st_transform(wa_map, common_crs), dist = distance, nQuadSegs = 100) %>%
#   st_make_valid()
# plot(land_buffer$geometry)
# 
# buf <- st_difference(land_buffer, st_transform(wa_map, common_crs)) %>% # takes a long while because the coastline is detailed.
#   st_make_valid() %>%
#   st_crop(bbox_whole)
# buf2 <- st_union(buf)
# 
# buf_poly <- st_cast(buf, "MULTIPOLYGON", warn = FALSE)
# 
# st_write(
#   buf_poly,
#   "data/output_data/01_B_buffer_empty.shp",
#   delete_layer = TRUE
# )

buf <- st_read("data/output_data/Q_buffer_empty_mainland.shp") %>%  # file saved just above, manually modified in QGIS to keep the mainland buffer (no islands)
  st_crop(bbox_whole)
plot(buf$geometry)

water_polygon <- st_difference(temp, st_union(wa_map)); plot(water_polygon, col = NA, border = "red", lwd = 2)
water_offshore_polygon <- st_difference(water_polygon, st_union(buf)) %>% st_transform(common_crs); plot(water_offshore_polygon, col = NA, border = "red", lwd = 2)


## over north area
small_grd_north <- st_make_grid(bbox_north, cellsize = 4000, square = T) %>% # 2000m x 2000m square, or 4km2
  st_as_sf() %>% 
  st_transform(common_crs); plot(small_grd_north)

small_grd_north <- st_intersection(small_grd_north, buf) %>% 
  st_make_valid() %>%
  st_transform(common_crs) %>% 
  dplyr::select(x)

small_grd_north <- st_intersection(small_grd_north, bbox_north) %>% 
  st_as_sf() %>% 
  st_collection_extract("POLYGON"); plot(small_grd_north)

small_grd_north$type <- "shore_north"


## over wadandi area
small_grd_wadandi <- st_make_grid(bbox_wadandi, cellsize = 2000, square = T) %>% # 2000m x 2000m square, or 4km2
  st_as_sf() %>% 
  st_transform(common_crs); plot(small_grd_wadandi)

small_grd_wadandi <- st_intersection(small_grd_wadandi, buf) %>% 
  st_make_valid() %>%
  st_transform(common_crs) %>% 
  dplyr::select(x)

small_grd_wadandi <- st_intersection(small_grd_wadandi, bbox_wadandi) %>% 
  st_as_sf() %>% 
  st_collection_extract("POLYGON"); plot(small_grd_wadandi)

small_grd_wadandi$type <- "shore_wadandi"


### 2. Offshore cells ---------------------------------------------------------

# large cells over north offshore cells.

## over north cells
hab_north <- crop(hab, st_transform(bbox_north, crs(hab)))
plot(hab_north)

big_grd_north <- st_make_grid(bbox_north, cellsize = 10000, square = T) %>% st_transform(common_crs); plot(big_grd_north) # 10000m x 10000m square, or 100km2

big_grd_north <- st_intersection(big_grd_north, water_offshore_polygon) %>% # crop the grid to the extent of the habitat raster
  st_make_valid() %>%
  st_transform(common_crs)

big_grd_north <- st_intersection(big_grd_north, bbox_north) %>% 
  st_as_sf() %>% 
  st_collection_extract("POLYGON"); plot(big_grd_north)

big_grd_north <- st_difference(big_grd_north, buf) %>% dplyr::select(x); plot(big_grd_north)

big_grd_north$type <- "offshore_north"


## over wadandi cells
hab_wadandi <- crop(hab, st_transform(bbox_wadandi, crs(hab)))
plot(hab_wadandi)

big_grd_wadandi <- st_make_grid(bbox_wadandi, cellsize = 4000, square = T) %>% st_transform(common_crs); plot(big_grd_wadandi) # 4000m x 4000m square, or 16km2

big_grd_wadandi <- st_intersection(big_grd_wadandi, water_offshore_polygon) %>% # crop the grid to the extent of the habitat raster
  st_make_valid() %>%
  st_transform(common_crs)

big_grd_wadandi <- st_intersection(big_grd_wadandi, bbox_wadandi) %>% 
  st_as_sf() %>% 
  st_collection_extract("POLYGON"); plot(big_grd_wadandi)

big_grd_wadandi <- st_as_sf(big_grd_wadandi)
big_grd_wadandi <- st_difference(big_grd_wadandi, buf) %>% dplyr::select(x); plot(big_grd_wadandi)

big_grd_wadandi$type <- "offshore_wadandi"


### 3. Merge all grids together -----------------------------------------------

# merge the big and small grids together
shore_grd <- rbind(small_grd_north, small_grd_wadandi) %>% st_transform(common_crs) %>% st_crop(bbox_whole); plot(shore_grd)
offshore_grd <- rbind(big_grd_north, big_grd_wadandi) %>% st_transform(common_crs) %>% st_crop(bbox_whole); plot(offshore_grd)

full_grd <- rbind(shore_grd, offshore_grd) %>% st_transform(common_crs) %>% st_crop(bbox_whole) %>% st_make_valid(); plot(full_grd)

water_vect <- vect(full_grd) %>% project("EPSG:4326"); plot(water_vect)
water_sf <- full_grd %>% 
  st_make_valid() %>% 
  st_cast("MULTIPOLYGON") %>%  # cast to MULTIPOLYGON to preserve complex geometry
  st_transform(common_crs)

# remove grid cells southeast of black point
water_ll <- st_transform(water_sf, 4326)   # transform to lat and long just to cut out
coords <- st_coordinates(st_centroid(water_ll))
water <- water_ll[c(
  coords[, "X"] <= 115.543734 |
    coords[, "Y"] >= -34.413498), 
]
grd_final <- st_transform(water, common_crs)
plot(grd_final)

# Check the grid
ggplot(data = grd_final) +
  geom_sf(color = colour_palette[2], fill = NA) +
  theme_minimal()

# save the grid and clean up in QGIS (removing shore cells that don't perfectly align with the large grids, at the south extremity)
st_write(grd_final, "data/output_data/01_B_grid_to_cleanup.shp", append = F)
grd_final <- st_read("data/output_data/Q_grid_cleaned_up_01_B.shp")


### 4. Extract habitat values to the grid cells -------------------------------

mean_values <- extract(hab_raster, st_transform(grd_final, crs(hab_raster)), fun = mean, na.rm = TRUE); head(mean_values)
summary(is.na(mean_values)) # check there are no NAs
mean_values[is.na(mean_values)] <- 0 # if there are NAs, they are because one habitat is classified as 100% of the cell area. fill in with zeroes

water_vect <- vect(grd_final) %>% project("EPSG:4326"); plot(water_vect)
water_vect <- cbind(water_vect, mean_values)
water_sf <- st_as_sf(water_vect) %>% 
  st_make_valid() %>% 
  st_transform(common_crs)

plot(water_sf)

saveRDS(water_sf, "data/output_data/01_B_water.rds")

## END ##

