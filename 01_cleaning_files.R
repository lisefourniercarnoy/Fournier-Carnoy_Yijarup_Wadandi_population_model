# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Claude's SW habitat predictions
# Task:    Create a shapefile with grid cells
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    November 2025

# -----------------------------------------------------------------------------

# Status:  Currently running trial grids 

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

file_land         <- "data/input_data/aus_land_low_res.shp"
#file_land         <- "data/input_data/aus_land_high_res_no_freshwater.shp"
#file_land         <- "data/input_data/aus_land_high_res.shp" # for plotting only
file_MPA          <- "data/input_data/Collaborative_Australian_Protected_Area_Database_(CAPAD)_2024_-_Marine_SW_selection.shp"
file_MPA          <- "data/input_data/western-australia_marine-parks-all.shp"
file_habitats     <- "data/input_data/wadandi_predicted_habitat.RDS" # Claude's SWC habitat predictions
file_shore_hab    <- "data/input_data/wadandi_shore_hab_categorised.shp" # manually categorised shore habitat


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
  st_crop(bbox_whole) # crop to Wadandi

wa_map <- st_simplify(wa_map, dTolerance = 100)  # Adjust dTolerance for balance between speed and precision

plot(wa_map$geometry, col = "#eeeeeeff") # Map of WA, cropped to Wadandi Country
wa_map <- wa_map |> dplyr::select(-Shape_Area) # remove Shape_Area which is too large to fit in the shapefile column
st_write(wa_map, "data/output_data/01_wadandi_land.shp", append = F)


## Map of the State + Commonwealth Marine Park --------------------------------

MP <- st_read(file_MPA) %>%
  st_transform(common_crs) %>% # set the projection
  st_make_valid() %>% # fix broken geometries
  st_crop(bbox_whole) %>%  # crop to Wadandi
  filter(zone %in% c("Sanctuary Zone", 
                     "National Park Zone", 
                     "Special Purpose Zone"
  )
  )
MP <- MP %>% filter(zone_type != c("Special Purpose Zone (Mining Exclusion) (IUCN VI)"))
plot(MP$geometry)

NTZ <- MP %>%
  #dplyr::select(GIS_AREA, GAZ_DATE, LATEST_GAZ, COMMENTS, ZONE_TYPE) %>% 
  st_crop(bbox_whole)
plot(NTZ)

# Save the No-Take Zones layer
plot(NTZ$geometry, col = colour_palette[4]); plot(wa_map$geometry, col = "#eeeeeeff", add = T)
st_write(NTZ, "data/output_data/01_wadandi_NTZ.shp", delete_layer = T)


### Habitat Files - RUN AGAIN IF THINGS CHANGE --------------------------------

sand <- rast("data/external_data/sw-network_predicted-habitat_psand.tif")
reef <- rast("data/external_data/sw-network_predicted-habitat_pinverts.tif") +
  rast("data/external_data/sw-network_predicted-habitat_pmacro.tif") + 
  rast("data/external_data/sw-network_predicted-habitat_prock.tif")
seagrass <- rast("data/external_data/sw-network_predicted-habitat_pseagrass.tif")

hab <- c(sand, reef, seagrass); names(hab) <- c("sand", "reef", "seagrass")

#  Crop the raster
bbox_sf <- st_transform(bbox_whole, crs(hab))
bbox_vect <- vect(bbox_sf)
cropped_ras <- crop(hab, bbox_vect)
plot(cropped_ras)

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
# st_write(shore_hab, "data/output_data/01_wadandi_shore_habitat_to_categorise.shp", delete_layer = T)

shore_hab <- st_read(file_shore_hab)
plot(shore_hab)
wa_map_union <- st_union(wa_map)
shore_hab <- st_difference(shore_hab, wa_map_union) %>% st_transform(st_crs(hab)) # crop out the areas that are too high res.
plot(shore_hab)
# merge Claude's predicted habitat to the shore habitat

# sand
sand_shore <- rasterize(shore_hab[shore_hab$habitat == "sand",], hab, field =  global(hab[["sand"]], fun = "max", na.rm = TRUE)[1,1], touches = TRUE)
sand_layer <- hab[["sand"]]
sand_combined <- cover(sand_layer, sand_shore)
plot(sand_combined)

# reef
reef_shore <- rasterize(shore_hab[shore_hab$habitat == "reef",], hab, field =  global(hab[["reef"]], fun = "max", na.rm = TRUE)[1,1], touches = TRUE)
reef_layer <- hab[["reef"]]
reef_combined <- cover(reef_layer, reef_shore)
plot(reef_combined)

# seagrass (use the shore hab of sand and replace with zero, since we have no seagrass on the shore)
seagrass_shore <- rasterize(shore_hab[shore_hab$habitat == "sand",], hab, field =  global(hab[["sand"]], fun = "max", na.rm = TRUE)[1,1], touches = TRUE)
plot(seagrass_shore)
seagrass_shore[[1]][!is.na(seagrass_shore[[1]]), ] <- 0
seagrass_layer <- hab[["seagrass"]]
seagrass_combined <- cover(seagrass_layer, seagrass_shore)
plot(seagrass_combined)

# the areas in sand_combined that are rocky are set to zero and vice versa
sand_mask <- !is.na(sand_shore); plot(sand_mask)
reef_combined_clean <- mask(reef_combined, sand_mask, maskvalues = TRUE, updatevalue = 0); plot(reef_combined_clean)

rocky_mask <- !is.na(rocky_shore); plot(rocky_mask)
sand_combined_clean <- mask(sand_combined, rocky_mask, maskvalues = TRUE, updatevalue = 0); plot(sand_combined_clean)

# combine into one
hab_raster <- c(sand_combined_clean, reef_combined_clean, seagrass_combined)
hab_raster <- crop(hab_raster, st_transform(bbox_whole, 4326)); plot(hab_raster)



## Make grid cells ------------------------------------------------------------

# make a polygon to use to make grids
hab_polygon <- st_as_sf(as.polygons(app(hab_raster[[1]], fun = function(x) ifelse(is.na(x), NA, 1)), dissolve = TRUE)) %>% 
  st_transform(common_crs); plot(hab_polygon, col = NA, border = "red", lwd = 2)

# create a 500m buffer around land
distance <- 500 # Buffer distance in meters (adjust as needed)

land_buffer <- st_buffer(wa_map, dist = distance, nQuadSegs = 100) %>%
  st_make_valid() %>%
  st_transform(common_crs)
plot(land_buffer$geometry)

buf <- st_difference(land_buffer, wa_map) %>% 
  st_make_valid() %>% 
  st_crop(bbox_whole)
plot(buf$geometry)


## 1. grid over wadandi ------------------------------------------------------- 

# make a large grid over habitat
big_grd <- st_make_grid(bbox_wadandi, cellsize = 4500, square = T); plot(big_grd) # 4500m x 4500m square, or 20.25km2

big_grd <- st_intersection(big_grd, hab_polygon) %>% # crop the grid to the extent of the habitat raster
  st_make_valid() %>%
  st_transform(common_crs); plot(big_grd)

big_grd <- st_as_sf(big_grd)
big_grd <- st_difference(big_grd, buf) %>% dplyr::select(x); plot(big_grd)
big_grd <- st_difference(big_grd, wa_map) %>% dplyr::select(x); plot(big_grd)

big_grd$type <- "pelagic"

# make a small grid inside that buffer
small_grd <- st_make_grid(bbox_wadandi, cellsize = 2000, square = T) %>% # 2000m x 2000m square, or 4km2
  st_as_sf() %>% 
  st_transform(common_crs); plot(small_grd)

small_grd <- st_intersection(small_grd, buf) %>% 
  st_make_valid() %>%
  st_transform(common_crs) %>% 
  dplyr::select(x)
plot(small_grd)
small_grd$type <- "shore"

# merge the big and small grids together
full_grd <- rbind(big_grd, small_grd) %>% st_transform(common_crs) %>% st_crop(bbox_wadandi); plot(full_grd)

# transform into a vector for use later
full_grd <- full_grd[st_geometry_type(full_grd) %in% c("POLYGON", "MULTIPOLYGON"), ] %>%  # remove problematic geometries
  st_make_valid()

water_vect <- vect(full_grd) %>% project("EPSG:4326"); plot(water_vect)
water_sf <- st_as_sf(water_vect) %>% 
  st_make_valid() %>% 
  st_cast("POLYGON") %>% 
  st_transform(common_crs)

# remove grid cells southeast of black point
water_ll <- st_transform(water_sf, 4326)   # transform to lat and long just to cut out
coords <- st_coordinates(st_centroid(water_ll))
water <- water_ll[c(
  coords[, "X"] <= 115.543734 |
  coords[, "Y"] >= -34.413498), 
]
grd_wadandi <- st_transform(water, common_crs)
plot(grd_wadandi)

# Check the grid
ggplot(data = grd_wadandi) +
  geom_sf(color = colour_palette[2]) +
  theme_minimal()


## 2. make a grid over north cells --------------------------------------------

# make a polygon to use to make grids
hab_north <- crop(hab, st_transform(bbox_north, crs(hab)))
plot(hab_north)

big_grd_north <- st_make_grid(bbox_north, cellsize = 10000, square = T); plot(big_grd_north) # 10000m x 10000m square, or 100km2

big_grd_north <- st_intersection(big_grd_north, hab_polygon) %>% # crop the grid to the extent of the habitat raster
  st_make_valid() %>%
  st_transform(common_crs); plot(big_grd_north)

big_grd_north <- st_as_sf(big_grd_north)
big_grd_north <- st_difference(big_grd_north, buf) %>% dplyr::select(x); plot(big_grd_north)
big_grd_north <- st_difference(big_grd_north, wa_map) %>% dplyr::select(x); plot(big_grd_north)

big_grd_north$type <- "north"

# make a small grid inside that buffer
small_grd_north <- st_make_grid(bbox_north, cellsize = 4000, square = T) %>% # 4000m x 4000m square, or 16km2
  st_as_sf() %>% 
  st_transform(common_crs); plot(small_grd_north)

small_grd_north <- st_intersection(small_grd_north, buf) %>% 
  st_make_valid() %>%
  st_transform(common_crs) %>% 
  dplyr::select(x)
plot(small_grd_north)
small_grd_north$type <- "shore_north"

# merge the big and small grids together
grd_north <- rbind(big_grd_north, small_grd_north) %>% st_crop(bbox_north) %>% st_transform(common_crs); plot(grd_north)


## 3. put north and wadandi together ------------------------------------------
ggplot() +
  geom_sf(data = grd_north) +
  geom_sf(data = grd_wadandi) +
  geom_sf(data = wa_map)

# clean up the grid object type otherwise it's a mess
grd_all <- rbind(st_sf(grd_north) %>% rename("geometry" = x), st_sf(grd_wadandi)) %>% st_make_valid(); plot(grd_all)
grd_all_polygons <- grd_all[st_geometry_type(grd_all) %in% c("POLYGON", "MULTIPOLYGON"),]
grd_all <- st_cast(grd_all_polygons, "POLYGON"); plot(grd_all)


# 4. Extract habitat values to the grid cells
mean_values <- extract(hab_raster, st_transform(grd_all, crs(hab_raster)), fun = mean, na.rm = TRUE); head(mean_values)
summary(is.na(mean_values)) # check there are no NAs

water_vect <- vect(grd_all) %>% project("EPSG:4326"); plot(water_vect)
water_vect <- cbind(water_vect, mean_values)
water_sf <- st_as_sf(water_vect) %>% 
  st_make_valid() %>% 
  #st_cast("POLYGON") %>% 
  st_transform(common_crs)

plot(water_sf)

saveRDS(water_sf, "data/output_data/01_water.rds")

## END ##
