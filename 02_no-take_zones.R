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

file_NTZ    <- "data/output_data/Q_NTZ_manual_clean_01_B.shp"
file_water  <- "data/output_data/01_B_water.rds"

fleets <- c("commercial", "boat_rec", "shore_rec") # these are the fleets from 01_B "MP_permissions" that describe fishing per marine park zone


## Load files -----------------------------------------------------------------

NTZ <-  st_read(file_NTZ) %>%
  st_transform(common_crs) %>%
  mutate(
    SC_restriction_date = case_when(name == "Rottnest"           ~ "01/07/2007", # add when exactly each NTZ was put in place
                                    name == "Shoalwater Islands" ~ "01/01/2007",
                                    name == "Marmion"            ~ "30/01/1992",
                                    name == "Two Rocks"          ~ "01/07/2018",
                                    name == "Ngari Capes"        ~ "10/04/2019",
                                    name == "Geographe"          ~ "01/07/2018",
                                    name == "South-west Corner"  ~ "01/07/2018",
                                    .default = NULL # others aren't in the final grid
    ),
    SC_restriction_date = as.list(SC_restriction_date)
  ) %>% 
  st_make_valid()
plot(NTZ$geometry)

water <- readRDS(file_water) %>%
  st_transform(common_crs) %>%
  st_make_valid()
plot(water)


## Make a temporal closure area -----------------------------------------------

# TC = temporal closure
TC_poly <- matrix(c(
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
  data.frame(lon = TC_poly[,1], lat = TC_poly[,2]),
  coords = c("lon", "lat"),
  crs = 4326
)
TC <- st_convex_hull(st_union(pts))
TC <- st_as_sf(st_transform(TC, st_crs(water)))

sf::sf_use_s2(FALSE)
TC <- st_intersection(TC, st_make_valid(st_union(water)))

plot(TC, col = NA, border = "red", lwd = 2); plot(water$geometry, add = T)

# below is a list and description of temporal restrictions. i've made it this way because it's easier to see and modify in the code, in case there are more temporal closures or more complex ones.
temporal_closure_details = list(
#   YEAR,   MONTHS,           PERCENT OF EACH MONTH THAT IS FISHED
  c("2000", "9-10",           "0.5-0"),
  c("2005", "10-11-12",       "0-0-0.5"),
  c("2023", "8-9-10-11-12-1", "0-0-0-0-0-0")
)

TC <- TC %>% 
  mutate(
    name = "Cockburn Sound temporal closure", # WHATEVER NAME IS HERE SHOULD HAVE THE EXACT TERMS 'temporal closure' IN IT OR THE REST WONT WORK
    TC_restriction_date = list(as.integer(vapply(temporal_closure_details, `[`, "", 1))),
    TC_restriction_months = list(vapply(temporal_closure_details, `[`, "", 2)),
    TC_restriction_date_perc_fished = list(vapply(temporal_closure_details, `[`, "", 3))
  ) %>% 
  rename(
    geometry = x
  )

# 2000: 15sep - 31oct, source: https://www.wa.gov.au/government/media-statements/Court%20Coalition%20Government/Cockburn-Sound-spawning-closure-gets-go-ahead-20000830
# 2005: 01oct - 15dec source: https://www.wa.gov.au/government/media-statements/Gallop%20Labor%20Government/Spawning-closures-protect-Perth%27s-pink-snapper-20050819
# 2023: 01aug - 31jan, source: https://www.wa.gov.au/government/announcements/recreational-demersal-fishing-closure-aid-stock-recovery

# add the closures in the NTZ dataframe
NTZ <- bind_rows(TC, NTZ)
plot(NTZ[, !sapply(NTZ, is.list)])


## Adjusting the grid to account for temporal/spatial closures ----------------

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

fished_area <- st_difference(water, NTZ_union) %>% 
  st_make_valid() %>%
  st_transform(common_crs)
plot(fished_area$geometry)

# Put this back together with the fished area to create "water" again
fished_area <- st_sf(fished_area) %>%
  mutate(status = "Fished",
         SC_restriction_date = NA,
         TC_restriction_date = NA,
         TC_restriction_months = NA,
         TC_restriction_date_perc_fished = NA)
plot(fished_area$geometry)
fished_area[,fleets] <- NA

NTZ_union <- st_sf(NTZarea) %>%
  mutate(status = case_when(grepl("temporal closure", name, fixed = TRUE) ~ "TC",
                            .default = "SC")) %>% 
  dplyr::select(names(fished_area), all_of(fleets))

water <- rbind(NTZ_union, fished_area)
plot(water$geometry)


# Check that the NTZs are where you expect them to be
ggplot(water) +
  geom_sf(aes(fill = status)) +
  theme_void() +
  scale_fill_manual(values = c(colour_palette[4], colour_palette[5], colour_palette[6], colour_palette[2]))

# give a new cell ID to all cells, because a few were cut in two in the process
water$ID <- 1:nrow(water)


## Identify spawning ground ---------------------------------------------------

# sg = spawning ground

# make a rough polygon of cockburn + warnbro sound
sg_poly <- matrix(c(
  115.742245, -32.050302,  # freo
  115.635004, -32.066268,  # straggler_rock
  115.838059, -32.365724,   # point_kennedy
  115.659383, -32.373515,    # warnbro
  115.742245, -32.050302   # close polygon back to freo
), 
ncol = 2, byrow = TRUE,
dimnames = list(NULL, c("lon", "lat")))

pts <- st_as_sf(
  data.frame(lon = sg_poly[,1], lat = sg_poly[,2]),
  coords = c("lon", "lat"),
  crs = 4326
)
sg <- st_convex_hull(st_union(pts))
sg <- st_transform(sg, st_crs(water))
plot(water$geometry); plot(sg, col = NA, border = "red", lwd = 2, add = TRUE)

# find the cell centroid that are inside the cockburn polygon
plot(st_centroid(water[, !sapply(water, is.list)]))
inside <- st_within(st_centroid(water), sg)
sg_cell_id <- which(lengths(inside) > 0)

intersections <- st_within(st_centroid(water$geometry), sg)
intersects_logical <- lengths(intersections) > 0

sg_cells <- water[intersects_logical, ]
plot(sg_cells[, !sapply(sg_cells, is.list)])
sg_cell_id <- sg_cells$ID # find the cell IDs

# modify cockburn cells and change their type to identify them.
water$spawning_status <- case_when(
  water$ID %in% sg_cell_id ~ TRUE,
  .default = FALSE
  )

ggplot() + # plot check AFTER adding spawning type
  geom_sf(data = water, aes(fill = spawning_status), col = NA) +
  scale_fill_manual(values = colour_palette) +
  theme_minimal()

# Calculating grid cell area and removing cells < 1m2
water <- water %>%
  mutate(cell_area = st_area(geometry),
         ID = row_number()) %>% 
  filter(as.numeric(cell_area)>1)
water <- st_make_valid(water) %>% 
  st_as_sf()


## Save files for next step ---------------------------------------------------

saveRDS(st_as_sf(water), file = "data/output_data/02_watergrid.rds")

### END ###
