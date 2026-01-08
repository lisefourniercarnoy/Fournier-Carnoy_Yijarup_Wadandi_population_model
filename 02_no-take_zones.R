# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Spatial layers from 01_Cleaning-files.R
# Task:    Intersect habitat grid with NTZ?
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    August 2024

# -----------------------------------------------------------------------------

# Status:  So far so good, may need adjustments from 01

# -----------------------------------------------------------------------------

library(tidyverse) # to manipulate data
library(ggplot2) # for pretty plots
library(sf) # to monipulate spatial layers

rm(list = ls()) # Clean working environment

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))


## Set projection and extent --------------------------------------------------

common_crs = 4283
bbox <- st_bbox(c(xmin = 114.4, xmax = 116.0, ymin = -35.3, ymax = -32.5), crs = common_crs)

## Files used in this script --------------------------------------------------

file_NTZ    <- "data/output_data/01_B_wadandi_NTZ.shp"
file_water  <- "data/output_data/01_B_water.rds"

## Load files -----------------------------------------------------------------

NTZ <-  st_read(file_NTZ) %>%
  st_transform(common_crs) %>%
  mutate(
    restriction_date = case_when(name == "Rottnest"           ~ "01/07/2007", # add when exactly each NTZ was put in place
                                 name == "Shoalwater Islands" ~ "01/01/2007",
                                 name == "Marmion"            ~ "30/01/1992",
                                 name == "Two Rocks"          ~ "01/07/2018",
                                 name == "Ngari Capes"        ~ "10/04/2019",
                                 name == "Geographe"          ~ "01/07/2018",
                                 name == "South-west Corner"  ~ "01/07/2018",
                                 .default = NULL # others aren't in the final grid
                                 ),
    restriction_date = as.list(restriction_date)
  ) %>% 
  st_make_valid()
plot(NTZ$geometry)

water <- readRDS(file_water) %>%
  st_transform(common_crs) %>%
  st_make_valid()
plot(water)

## make the cockburn sound pink snapper closure -------------------------------

test <- matrix(c(
  115.749655, -32.073611, # NE point
  115.608333, -32.073611, # NW point
  115.652778, -32.380851, # SW point
  115.724430, -32.380851, # becher point
  115.849846, -32.313575, # close polygon back to freo
  115.841356, -32.167078  # close polygon back to freo
), 
ncol = 2, byrow = TRUE,
dimnames = list(NULL, c("lon", "lat")))

pts <- st_as_sf(
  data.frame(lon = test[,1], lat = test[,2]),
  coords = c("lon", "lat"),
  crs = 4326
)
cs_closure <- st_convex_hull(st_union(pts))
cs_closure_proj <- st_as_sf(st_transform(cs_closure, st_crs(water)))

sf::sf_use_s2(FALSE)
cs_closure_proj <- st_intersection(cs_closure_proj, st_make_valid(st_union(water)))

plot(cs_closure_proj, col = NA, border = "red", lwd = 2); plot(water$geometry, add = T)

cs_closure_proj <- cs_closure_proj %>% 
  mutate(
    name = "Cockburn Sound temporal closure",
    restriction_date = list(c(2000, 2005, 2023))
  ) %>% 
  rename(
    geometry = x
  )

# 2000: 6 weeks from 15sep 31oct, source: https://www.wa.gov.au/government/media-statements/Court%20Coalition%20Government/Cockburn-Sound-spawning-closure-gets-go-ahead-20000830
# 2005: 4 months from 01oct 15dec source: https://www.wa.gov.au/government/media-statements/Gallop%20Labor%20Government/Spawning-closures-protect-Perth%27s-pink-snapper-20050819
# 2023: 01aug to 31jan, source: https://www.wa.gov.au/government/announcements/recreational-demersal-fishing-closure-aid-stock-recovery

# closure file
NTZ <- bind_rows(cs_closure_proj, NTZ)
plot(NTZ)


## Adjusting the grid to account for NTZ --------------------------------------

# Make a grid for areas that are not fished
NTZarea <- st_intersection(NTZ, water) %>% 
  st_make_valid() %>%
  st_transform(common_crs)
plot(NTZarea$geometry)

# Make a grid for areas that are fished
NTZ_union <- st_union(NTZarea) %>%
  st_make_valid() %>%
  st_transform(common_crs)
plot(NTZ_union)

Fished_area <- st_difference(water, NTZ_union) %>% 
  st_make_valid() %>%
  st_transform(common_crs)
plot(Fished_area$geometry)

# Put this back together with the fished area to create "water" again
df_Fished_area <- st_sf(Fished_area) %>%
  mutate(status = "Fished",
         restriction_date = NA)
plot(df_Fished_area)

df_NTZ_union <- st_sf(NTZarea) %>%
  mutate(status = ifelse(zone_type == "Special Purpose Zone (Shore-based Activities) (IUCN VI)", # for cells in this shore-activities-only zones...
                         "NTZ_for_boat_only", "NTZ_for_shore_and_boat"), # ... boat fishing is not allowed
         status = ifelse(name == "Cockburn Sound temporal closure", 
                         "temporal_closure", status)) %>%
  dplyr::select(names(df_Fished_area))

names(df_Fished_area); names(df_NTZ_union)

water <- rbind(df_NTZ_union, df_Fished_area)
plot(water$geometry)


# Check that the NTZs are where you expect them to be
ggplot(water) +
  geom_sf(aes(fill = status)) +
  theme_void() +
  scale_fill_manual(values = c(colour_palette[4], colour_palette[5], colour_palette[6], colour_palette[2]))

# give a new cell ID to all cells, because a few were cut in two in the process
water$ID <- 1:nrow(water)


## find cockburn sound --------------------------------------------------------

# make a rough polygon of cockburn + warnbro sound
test <- matrix(c(
  115.742245, -32.050302,  # freo
  115.635004, -32.066268,  # straggler_rock
  115.838059, -32.365724,   # point_kennedy
  115.659383, -32.373515,    # warnbro
  115.742245, -32.050302   # close polygon back to freo
), 
ncol = 2, byrow = TRUE,
dimnames = list(NULL, c("lon", "lat")))

pts <- st_as_sf(
  data.frame(lon = test[,1], lat = test[,2]),
  coords = c("lon", "lat"),
  crs = 4326
)
cs <- st_convex_hull(st_union(pts))
cs_proj <- st_transform(cs, st_crs(water))
plot(water$geometry); plot(cs_proj, col = NA, border = "red", lwd = 2, add = TRUE)

# find the cell centroid that are inside the cockburn polygon
plot(st_centroid(water))
inside <- st_within(st_centroid(water), cs_proj)
cs_cell_id <- which(lengths(inside) > 0)

intersections <- st_within(st_centroid(water$geometry), cs_proj)
intersects_logical <- lengths(intersections) > 0

cockburn_cells <- water[intersects_logical, ]
plot(cockburn_cells)
cs_cell_id <- cockburn_cells$ID # find the cell IDs


ggplot() + # plot check BEFORE adding cockburn type
  geom_sf(data = water, aes(fill = type), col = NA) +
  scale_fill_manual(values = colour_palette) +
  theme_minimal()

# modify cockburn cells and change their type to identify them.
water$type <- case_when(
  water$ID %in% cs_cell_id & water$type == "shore_north" ~ "shore_north_cockburn_warnbro",
  water$ID %in% cs_cell_id ~ "north_cockburn_warnbro",
  .default = water$type
  )

ggplot() + # plot check AFTER adding cockburn type
  geom_sf(data = water, aes(fill = type), col = NA) +
  scale_fill_manual(values = colour_palette) +
  theme_minimal()

# Calculating grid cell area and removing cells < 1m2
water <- water %>%
  mutate(cell_area = st_area(geometry),
         ID = row_number()) %>% 
  filter(as.numeric(cell_area)>1)
water <- st_make_valid(water) %>% 
  st_as_sf()


# identify important cells

# NTZ cells
no_take_list_boat_shore <- water %>%
  filter(status == "NTZ_for_boat_only" | status == "NTZ_for_shore_and_boat") %>%
  st_drop_geometry() %>%
  pull(ID) %>%
  as.numeric()

no_take_list_shore <- water %>%
  filter(status == "NTZ_for_shore_and_boat") %>%
  st_drop_geometry() %>%
  pull(ID) %>%
  as.numeric()


# cockburn sound + warnbro sound cells
cs_list <- water$ID[water$type %in% c("north_cockburn_warnbro", "shore_north_cockburn_warnbro")]
 

## Save files for next step ---------------------------------------------------

saveRDS(no_take_list_boat_shore, file = "data/output_data/02_no_take_list_boat_shore.rds")
saveRDS(no_take_list_shore, file = "data/output_data/02_no_take_list_shore.rds")
saveRDS(cs_cell_id, "data/output_data/02_cockburn_cell_id.rds")
saveRDS(st_as_sf(water), file = "data/output_data/02_watergrid.rds")

### END ###
