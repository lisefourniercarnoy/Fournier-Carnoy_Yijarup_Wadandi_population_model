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
library(sf) # to manipulate spatial layers

rm(list = ls()) # Clean working environment

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))


## Set projection and extent --------------------------------------------------

common_crs = 4283
bbox <- st_bbox(c(xmin = 114.4, xmax = 116.0, ymin = -35.3, ymax = -32.5), crs = common_crs)


## Files used in this script --------------------------------------------------

file_NTZ    <- "data/output_data/Q_NTZ_manual_clean_01_B.shp"
file_water  <- "data/output_data/01_B_water.rds"
file_land   <- "data/input_data/Q_aus_land_high_res_no_estuary.shp"

fleets <- c("commercial", "boat_rec", "shore_rec") # these are the fleets from 01_B "MP_permissions" that describe fishing per marine park zone


## Load files -----------------------------------------------------------------

land <- st_read(file_land) %>%
  st_transform(common_crs)

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


# add a commercial closure in the metro area
water <- readRDS(file_water) %>%
  st_transform(common_crs) %>%
  st_make_valid()
plot(water)

com_metro_closure <- water %>% 
  dplyr::filter(type %in% c("shore_north", "offshore_north")) %>% 
  summarise(geometry = st_union(geometry)) %>%
  st_sf(geometry = .) %>%  # convert geometry back to sf
  mutate( # make sure names match with NTZ above
    X = NA,
    name = "Metro closure for commercial fleet",
    zone_type = NA,
    zone = NA,
    epbc = NA,
    commercial = FALSE,
    boat_rec = TRUE,
    shore_rec = TRUE,
    SC_restriction_date = "15/11/2007" # see source : https://www.abc.net.au/news/2007-11-15/commercial-fisherman-angry-about-new-bans/726090 
  )

plot(com_metro_closure$geometry)
plot(NTZ$geometry, add = T, col = "red")

com_closure_only <- st_difference(com_metro_closure, st_union(NTZ))
plot(com_closure_only$geometry, col = "gray")
plot(NTZ$geometry, add = T, col = "red")

NTZ <- rbind(com_closure_only, NTZ)
plot(NTZ$geometry, col = "red")


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
    TC_restriction_perc_fished = list(vapply(temporal_closure_details, `[`, "", 3))
  ) %>% 
  rename(
    geometry = x
  )

# 2000: 15sep - 31oct, source: https://www.wa.gov.au/government/media-statements/Court%20Coalition%20Government/Cockburn-Sound-spawning-closure-gets-go-ahead-20000830
# 2005: 01oct - 15dec source: https://www.wa.gov.au/government/media-statements/Gallop%20Labor%20Government/Spawning-closures-protect-Perth%27s-pink-snapper-20050819
# 2023: 01aug - 31jan, source: https://www.wa.gov.au/government/announcements/recreational-demersal-fishing-closure-aid-stock-recovery

# add the closures in the NTZ dataframe
closures <- bind_rows(TC, NTZ)
plot(closures[, !sapply(closures, is.list)])


## Adjusting the grid to account for temporal/spatial closures ----------------

# colnames that should be present: 
# SC_status SC_restriction_date
# TC_status TC_restriction_date, TC_restriction_months, TC_restriction_perc_fished
cols <- c("ID", "name", "zone_type", "zone", "epbc", "type", 
          "sand", "reef", "seagrass", 
          fleets, 
          "TC_status", "SC_status",
          "TC_restriction_date", "TC_restriction_months", "TC_restriction_perc_fished", "SC_restriction_date", 
          "geometry")

# obtained grid cells that are fished
fished <- st_difference(water, st_union(NTZ)) %>% # takes a while
  mutate(
    SC_status = FALSE,
    TC_status = FALSE,
    TC_restriction_date = list(NA),
    TC_restriction_months = list(NA),
    TC_restriction_perc_fished = list(NA),
    SC_restriction_date = list(NA)
  )
fished$area <- st_area(fished) 
fished <- fished[as.numeric(fished$area) > 1,] # the st_difference leaves very small features around the coast - filter them out by keeping cells > 1m2
mapview::mapview(fished[!is.list(fished)]) # check things make sense


# obtain grid cells that are temporal closure (HERE THE CLOSURE OVERLAPS COMPLETELY WITH THE NTZ. OTHERWISE ST_DIFFERENCE IS NEEDED)
TC_area <- st_intersection(water, TC, sparse = F) %>% 
  mutate(zone_type = NA, 
         zone = NA) %>% 
  st_make_valid() %>%
  st_transform(common_crs)
TC_SC_overlap <- st_intersection(TC_area, NTZ) %>% 
  mutate(
    TC_status = TRUE,
    SC_status = TRUE
    )
mapview::mapview(TC_SC_overlap[!is.list(TC_SC_overlap)])


# Make a grid for areas that are spatial closure only
NTZarea <- st_intersection(NTZ, water) %>% 
  st_make_valid() %>%
  st_transform(common_crs)
plot(NTZarea$geometry)

SC_only <- st_difference(NTZarea, st_union(TC_SC_overlap)) %>%
  mutate(
    SC_status = TRUE,
    TC_status = FALSE,
    TC_restriction_date = list(NA),
    TC_restriction_months = list(NA),
    TC_restriction_perc_fished = list(NA)
  ) %>% 
  dplyr::select(any_of(cols)) %>% 
  dplyr::filter(st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON"))
SC_only$area <- st_area(SC_only)
SC_only <- SC_only[as.numeric(SC_only$area) > 1,] # the st_difference leaves very small features around the coast - filter them out by keeping cells > 1m2
mapview::mapview(SC_only[!is.list(SC_only)])

# check that everything lines up
plot(SC_only$geometry, col = "red")
plot(fished$geometry, col = "gray", add = T)
plot(TC_SC_overlap$geometry, col = "blue", add = T)


names(SC_only)
names(TC_SC_overlap)
names(fished)



water <- bind_rows(
  SC_only,
  TC_SC_overlap,
  fished
  ) %>% 
  dplyr::filter(st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON")) %>% 
  dplyr::select(any_of(cols))
plot(water$geometry, col = "red")
mapview::mapview(water[!is.list(water)])


# Check that the NTZs are where you expect them to be
ggplot(water) +
  geom_sf(aes(fill = SC_status), col = NA) +
  theme_void() +
  scale_fill_manual(values = c(colour_palette[4], colour_palette[5], colour_palette[6], colour_palette[2]))
ggplot(water) +
  geom_sf(aes(fill = TC_status), col = NA) +
  theme_void() +
  scale_fill_manual(values = c(colour_palette[4], colour_palette[5], colour_palette[6], colour_palette[2]))

mapview::mapview(water[!is.list(water)])

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
mapview::mapview(water[!is.list(water)])

# Calculating grid cell area and removing cells < 1m2
water <- water %>%
  mutate(cell_area = st_area(geometry),
         ID = row_number()) %>% 
  filter(as.numeric(cell_area)>1)
water <- st_make_valid(water) %>% 
  st_as_sf()

# give a new cell ID to all cells, because a few were cut in two in the process
water$ID <- 1:nrow(water)

mapview::mapview(water[!is.list(water)])

## Save files for next step ---------------------------------------------------

saveRDS(st_as_sf(water), file = "data/output_data/02_watergrid.rds")

### END ###
