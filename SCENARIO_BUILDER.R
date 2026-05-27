# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    ?
# Task:    Create many scenarios from a single script
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    December 2025

# -----------------------------------------------------------------------------

# Status: Brand new start.

# -----------------------------------------------------------------------------

rm(list = ls()) # clear environment

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))

library(tidyverse) # data manipulation
library(sf) # shapefile manipulation
library(raster) # raster manipulation
library(RColorBrewer) # plotting colours
library(abind) # matrix manipulation
library(sfnetworks) # for making distance to access points
library(exactextractr) # estracting raster values at points

# The point of this script is to set up a fishing effort surface over time based
# on user-defined spatial, temporal and fleet restrictions.

## find scenario files --------------------------------------------------------

# before this script, you should set up the restrictions of each fleet in the "qgis_projects/SCENARIO_BUILDER_MAP"
# make a dedicated folder for the shapefiles (exact name format: "data/input_data/(scenario name)_restrictions)
# you should save the corresponding fleet restrictions (exact format "Q_management_scenario_(fleet)_(scenario_name)"). the script is set up to find these objects seamlessly.

scenario_name <- "SXX_test"

n_restrictions <- st_read(paste0("data/input_data/", scenario_name, "_restrictions/Q_management_scenario_", scenario_name, ".shp"))

folder_name <- paste0("data/output_data/", scenario_name)
dir.create(file.path(folder_name), showWarnings = FALSE)


# load other files

file_wa           <- "data/output_data/01_B_land.shp"
file_bathy        <- "data/input_data/SW_crop_AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
file_boat_ramps_w <- "data/input_data/wadandi_boat_ramps.shp"
file_boat_ramps_n <- "data/input_data/north_boat_ramps.shp"
file_carpark_w    <- "data/input_data/wadandi_carparks.shp"
file_carpark_n    <- "data/input_data/north_carparks.shp"
file_water        <- "data/output_data/03_B_water.rds"
file_network      <- "data/output_data/03_B_network_shapefile.shp"

year_start <- 2025 # this is different from the 04 scripts because we're setting up a future fishing surface.
year_end <- 2050 # this is different from the 04 scripts because we're setting up a future fishing surface.
n_years_tot <- year_end - year_start + 1

common_crs <- 4283
crs_raster <- "+proj=longlat +datum=WGS84 +no_defs"
bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -31), crs = common_crs) %>%
  st_as_sfc() %>%
  st_transform(crs_raster)

water <- readRDS(file_water) %>% filter(!is.na(ID))
NCELL <- nrow(water)

fishing_info_list <- list()


# COMMERCIAL ==================================================================

current_fleet <- "commercial"

# 1. travel cost = distance from access points * fuel price -------------------

## 1.a fuel price -------------------------------------------------------------

# obtain fuel prices from the literature
fuel_price <- read.csv("data/input_data/fuel_prices WA population - forecast.csv") # data from various reports from government and international agencies.

plot(fuel_price$petrol_price ~ fuel_price$year, type = "l") # plot check

fuel_price_clean <- fuel_price$petrol_price

## 1.b distance from access points --------------------------------------------

water <- readRDS(file_water)

BR_n <- st_read(file_boat_ramps_n) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  dplyr::filter(!is.na(ABS_name)) %>% # select the ramps that are mentioned in the ABS reports - these are the commercial ramps
  mutate(build_year = 1900, # all commercial ramps start in 1900 cuz some of the build years dont make sense
         build_mnth = 1,     # assume Jan if unknown
         norm_popularity = as.numeric(com_prop) / sum(as.numeric(com_prop), na.rm = TRUE)) %>%
  glimpse()


BR_w <- st_read(file_boat_ramps_w) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  mutate(build_year = 1900, # all commercial ramps start in 1900 cuz some of the build years dont make sense
         build_mnth = 1,     # assume Jan if unknown
         norm_popularity = as.numeric(com_prop) / sum(as.numeric(com_prop), na.rm = TRUE)) %>%
  glimpse()

# for commercial fishing, only a few boat ramps are used, so select only correct ramps (see Andrea Gaynor's historical fishing resources, fishing localities in ABS stats)
unique(BR_w$name)
BR_w <- BR_w %>% 
  dplyr::filter(BR_w$name %in% c("SC_Augusta_Ellis_St_Jetty", "WC_Gnarabup", "WC_Hamelin_Bay",
                                 "GB_Quindalup", "GB_Eagle_Bay", "GB_Bunbury_Stirling_St", "GB_Busselton_Georgette_Street"))
plot(water$geometry); plot(BR_w$geometry, col = colour_palette[5], pch = 16, cex = 1, add = TRUE)

glimpse(BR_n)
glimpse(BR_w)

BR <- rbind(BR_n %>% dplyr::select(name, geometry) %>% mutate(region = "north"), 
            BR_w %>% dplyr::select(name, geometry) %>% mutate(region = "wadandi")) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>% 
  glimpse()

plot(water$geometry); plot(BR$geometry, col = colour_palette[6], pch = 16, add = TRUE)

network <- st_read(file_network) %>% 
  st_transform(common_crs); plot(network$geometry)
BR <- st_as_sf(BR); BR <- BR %>% st_transform(common_crs) 

# find the centre of each grid cell
sf::sf_use_s2(FALSE)
centroids <- st_centroid(st_make_valid(water))
sf::sf_use_s2(TRUE)
points <- as.data.frame(st_coordinates(centroids))%>% # the points start at the bottom left and then work their way their way right
  mutate(ID = row_number()) 
points_sf <- st_as_sf(points, coords = c("X", "Y")) 
st_crs(points_sf) <- common_crs

network <- as_sfnetwork(network, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

net <- activate(network, "nodes")
net <- net %>% st_transform(common_crs) 

# measure the distance from access points to cell centroids
network_matrix <- st_network_cost(net, from = BR, to = points_sf)
dim(network_matrix) # number of ramps x number of cells

glimpse(network_matrix)
access_dist <- as.data.frame(t(network_matrix))
colnames(access_dist) <- BR$name
access_dist$ID <- water$ID
access_dist <- access_dist / 1000
head(access_dist) # this gives us each cell's distance to the access points

## 1.c calculate travel cost --------------------------------------------------

glimpse(access_dist) # cell x access_point
glimpse(fuel_price_clean) # vec year

# figure out, for each ramp, the travel cost (cell x access_p x years)
dist_matrix <- as.matrix(access_dist[, -ncol(access_dist)])  # drop ID column

travel_cost <- array(rep(dist_matrix, times = length(fuel_price_clean)),
                     dim = c(nrow(access_dist), 
                             ncol(access_dist)-1, 
                             length(fuel_price_clean))
)


travel_cost <- sweep(
  travel_cost,
  MARGIN = 3,
  STATS  = fuel_price_clean,
  FUN    = "*"
) # calculate the product (fun = "*") of access_dist using fuel price (stats =) across years (margin = 3), 

## sanity check
par(mfrow = c(1,1))
test1 <- travel_cost[1,1,] # cell 1 from access_p 1
test2 <- travel_cost[1,2,] # cell 1 from access_p 2

plot(test2, type = "l", col = "red")
lines(test1, type = "l")


# 2. calculate cell utility ---------------------------------------------------

## 2.a. distance to shore -----------------------------------------------------

shore <- water[water$type %in% c("shore_wadandi", "shore_north"),] %>%
  st_make_valid() %>% 
  st_union() %>% 
  st_transform(st_crs(water)) %>%  # find shore. we will use this to calculate distance of all cells to shore.
  st_as_sf()

network <- st_read(file_network); plot(network$geometry)

centroids <- st_centroid(water %>% st_make_valid())
shore_dist <- st_distance(centroids, shore) / 1000 # distance from shore in km
shore_dist <- as.data.frame(shore_dist) %>% mutate(ID = water$ID, shore_dist = as.numeric(shore_dist))
glimpse(shore_dist)

# plot check
ggplot(data = water %>% dplyr::select(!where(is.list)) %>% mutate(shore_dist = as.numeric(shore_dist$shore_dist))) + 
  geom_sf(aes(fill = shore_dist), col = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4])


## 2.b cell utility = shore_dist / (travel cost + 1) --------------------------

glimpse(shore_dist)
glimpse(travel_cost)

# we're making cell x access_p x year

utility <- array(shore_dist$shore_dist,
                 dim = c(nrow(access_dist), 
                         ncol(access_dist)-1, 
                         length(fuel_price_clean))
)

utility <- (utility / (travel_cost + 1))

## sanity check - is utility higher for areas close to access points but far from shore)
test_access = 9
ggplot() +
  geom_sf(data = water %>% dplyr::select(!where(is.list)) %>% mutate(test = utility[, test_access, 1]), 
          aes(fill = (test)), 
          col = NA) +
  geom_sf(data = BR[test_access,] %>% st_set_crs(4326) %>% st_transform(st_crs(water)), col = "red") +
  scale_fill_gradientn(colours = colour_palette[6:4])
# ayoooo im a geniuuus


# 3. fishability --------------------------------------------------------------

# fishability refers to whether this cell in this month and this year is open to fishing.
# each fleet has a different fishability because not all closures/NTZs apply to everyone.

# fishability is between 0 and 1: 
# 0 not fishable for the entire month, 
# 1 fishable for the entire month
# decimals mean either:
## a. fishable for part of the month (in the case of temporal closures), or 
## b. partially fishable in area (related to the % of the cell that's within fishable depth)

# for the final cell x month x year cube, we need: 
## depth fishability (better engines over time allow deeper fishing)
## temporal fishability (cell x month x year)
## cell area

## we will calculate: (depth fishability x temporal fishability * cell area) / total fishable area

## 3.a depth fishability ------------------------------------------------------

# this section of the fishing effort relates to boats being able to fish further over the years, with bigger boats, more powerful engines etc.
# the model takes this into account by adding a 'fishable depth' into the mix, which goes more and more offshore.
# for commercial fishing, we're assuming that on average, fishers can fish 1.3m more every year. (change as needed in parameters)
NCELL <- nrow(water)

# prepare grid cells
plot(water$geometry)
wa_mask <- st_read(file_wa); wa_mask <- st_transform(wa_mask, crs = crs_raster); wa_mask <- as(wa_mask, "Spatial"); plot(wa_mask, col = "lightgray", add = T)
plot(bbox, add = T)

# set parameters (in 2024, commercial boats could fish in 122.7. we'll continue this 1.3m/year expansion)
years <- year_start:year_end
min_depth <- 122.7 # 20m max in 1945 + ((2024-1945) * 1.3m)
increase_year <- 2025 # the year that boats start to go deeper (in the scenarios, this is straightaway)
depth_per_year <- 1.3 # fishable depth gained per year (m) - arbitrary, but consistent with the historical reconstruction.

# prepare the bathymetry layer
bathy <- terra::rast(file_bathy) %>%
  terra::project(crs_raster) %>% 
  terra::crop(as(water, "Spatial")) %>%
  abs() %>%
  terra::mask(terra::vect(wa_mask), inverse = TRUE); plot(bathy)

thresholds <- min_depth + pmax(0, depth_per_year*(years - increase_year)) # depth limits for each year

fishable_stack <- terra::rast(
  lapply(thresholds, function(th) bathy <= th)
)

names(fishable_stack) <- paste0("year_", years)

fishable_summary <- terra::extract(
  fishable_stack,
  terra::vect(water),
  fun = mean,
  na.rm = TRUE
)

fishable_summary$ID <- water$ID

fishable_depth_cell_year <- fishable_summary %>%
  arrange(ID) %>%
  dplyr::select(-ID) %>%
  as.matrix()

# strip "year_" prefix from colnames
colnames(fishable_depth_cell_year) <- gsub("year_", "", colnames(fishable_depth_cell_year))

fishable_depth_cell_month_year <- array(
  aperm(replicate(12, fishable_depth_cell_year), c(1, 3, 2)),
  dim = c(NCELL, 12, length(years)),
  dimnames = list(NULL))

fishable_depth_cell_month_year[is.nan(fishable_depth_cell_month_year)] <- 1 # replace NaN with 1 (these are shore cells that have depth 0)


## sanity check station 
test_year = 10
test <- fishable_depth_cell_month_year[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()


## 3.b Temporal fishability ---------------------------------------------------

# this is just a cell x month x year cube that tells you whether there's a closure or not in the cell.
# some cells are closed for only part of the month (ex. temporal closure that starts midway through the month) 
# so that cell would have the value 0.5 that month.

glimpse(water)

## in this section (1.) we calculate the fishable area of the cells in each time step (1.a), then divide it by the fishable area of the whole grid at each time step (1.b)
## which gives us each the portion of Catchability of each cell.
## because the area that is catchable changes (with spatial and temporal restrictions), we have to calculate the grid's catchability for each month of each year.


### 3.b.a. find the fishable area of each cell in each time step --------------

water <- readRDS(file_water) %>% filter(!is.na(ID))
water <- merge(water, st_drop_geometry(n_restrictions), by = "ID")

NCELL <- nrow(water)

# identify the important cells
temporal_cells <- water$ID[which(water$TC_status == TRUE | water$n_TC_sttus == TRUE)] # we're selecting cells that are closed currently

# setup the arrays to calculate things.
water_area <- array(0, dim = c(NCELL, 12, (year_end-year_start+1))) # cells x months x years
water_area[, 1:12, 1] <- water$cell_area # for all months of the first year, the catchable area is the cell's area.

SC_check <- function(sc_d_string, fleet, current_year, current_month) {
  
  if (fleet == TRUE) return(FALSE) # if the closure applies to a different fleet than commercial (TRUE = fleet can fish, FALSE = fleet can't fish)
  
  if (is.na(sc_d_string) || sc_d_string == "") return(FALSE) # if there is no spatial restriction, we are not in a closure, so FALSE
  
  ranges <- trimws(strsplit(sc_d_string, ",")[[1]]) # split the different closure periods
  
  for (rng in ranges) {
    
    current_ym <- current_year * 12 + current_month
    
    parts <- strsplit(rng, "to")[[1]] # obtain the parts of the closure periods (start and end)
    start_parts <- strsplit(trimws(parts[1]), "-")[[1]]
    start_month <- as.integer(start_parts[1]); start_year <- as.integer(start_parts[2])
    start_ym   <- start_year  * 12 + start_month    
    
    if (length(parts) == 1) { # if there is only one date in the SC_d, the closure is valid for the entire simulation (no end date), so TRUE (if we've entered the closure start)
      if (current_ym >= start_ym) return(TRUE)
      next
    } 
    # otherwise get the end dates.
    end_parts   <- strsplit(trimws(parts[2]), "-")[[1]]
    end_month   <- as.integer(end_parts[1]);   end_year   <- as.integer(end_parts[2])
    end_ym     <- end_year    * 12 + end_month
    
    if (current_ym >= start_ym && current_ym <= end_ym) return(TRUE) # if you're in a closure period, you're in a closure, so TRUE
  }
  return(FALSE)
} # this is a function to make sense of the start/ends of the spatial closures.


# the loop goes as follows: for each time step, 
# check whether there's a spatial restriction that year, if so the fishable area is 0
# then, if the cell has a temporal closure, restrict fishable area as needed.
current_fleet_n <- "n_com"

for (YEAR in 1:dim(water_area)[3]) {
  
  current_year <- year_start + YEAR - 1
  print(current_year)
  
  for (MONTH in 1:dim(water_area)[2]) {
    
    # spatial closures
    
    water_area[, MONTH, YEAR] <- water$cell_area
    
    fishing_allowed <- mapply(SC_check, # here, we check whether each cell is open to fishing in this MONTH and YEAR. FALSE = not open, TRUE = open
                              sc_d_string  = water$n_SC_d,
                              fleet        = water$n_com,
                              MoreArgs     = list(current_year  = current_year,
                                                  current_month = MONTH))
    
    restricted_cells <- which(
      !is.na(water[[current_fleet_n]]) & water[[current_fleet_n]] == FALSE &
        fishing_allowed
    )
    water_area[restricted_cells, MONTH, YEAR] <- 0 # ...the cell is not fishable
    
    # temporal closures
    for (CELL in temporal_cells) {
      
      # find the relevant temporal restrictions
      TC_years_list <- water$n_TC_d[[CELL]] # years the closure has changed
      TC_months_list <- water$n_TC_m[[CELL]] # months restricted
      TC_perc_list <- water$n_TC_p[[CELL]] # how much they are restricted
      
      current_TC <- which(TC_years_list <= current_year) # only select the temporal closures that exist at the current point
      
      if (length(current_TC) > 0) { # if there are temporal restrictions at this YEAR, restrict the fishable area accordingly
        current_TC <- max(current_TC)
        
        TC_months_current <- as.integer(strsplit(TC_months_list[current_TC], "-")[[1]])
        TC_perc_current   <- as.numeric(strsplit(TC_perc_list[current_TC], "-")[[1]])
        
        if (MONTH %in% TC_months_current) {
          month_idx <- match(MONTH, TC_months_current)
          water_area[CELL, MONTH, YEAR] <- water_area[CELL, MONTH, YEAR] * TC_perc_current[month_idx]
        }
      }
    }
  }
} # TAKES FOREVER - this loop calculates for every month and every year, the fishable area of each cell.

### sanity check station

test_cell <- 200
test_month <- 10
test_year <- 23

# reference cell numbers as of 13.01.2026:
## cockburn sound cell (temporal closure): 1
## SWC NTZ cell: 200
## random fished cell: 1000

# restrictions in this cell should be:
glimpse(st_drop_geometry(water[test_cell,]))

# check that it is correct
cat("In month", test_month, "of year", (year_start+test_year-1),
    ", the cell is", ifelse(water_area[test_cell, test_month, test_year]>0, "fishable", "NOT fishable"), "for commercial boats.", "Fishable area: ", water_area[test_cell, test_month, test_year]/1e06, "km2")

## sanity check station 
test_year = 20
test <- water_area[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()


### 3.b.b calculate catchability ----------------------------------------------

glimpse(water_area)
glimpse(fishable_depth_cell_month_year)

fishable_area <- water_area * fishable_depth_cell_month_year
dim(water_area)
dim(fishable_depth_cell_month_year)

## sanity check station 
test_year = 20
test <- fishable_area[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()

# now divide by the sum of fishable area.
glimpse(fishable_area)
fishable_area_sum <- colSums(fishable_area) # month x year matrix

catchability <- sweep(fishable_area, c(2, 3), fishable_area_sum, FUN = "/")
dim(fishable_area)
dim(fishable_area_sum)


## sanity check station 
test_year = 20
test <- catchability[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()


# 4. set up fishing days values -----------------------------------------------

# we need to know how much fishing occurs in the region, so that the function knows how much to distribute.
# for commercial fishing, we are using a variety of sources to reconstruct the trends (a),
# the sources are split between the Metropolitan area effort (North, 4.N) and the Southwest area effort (Wadandi, 4.W)
# we will (b.) split the yearly effort into months, then (c.) split the monthly effort into access points (boat ramps).



### 4.N. North (Metropolitan) -------------------------------------------------

#### 4.N.a. enter the overall effort values (boat days) -----------------------

# obtain boat days from forecasting (see google slides on reconstruction)
future_boat_days <- read.csv("data/input_data/YIJARUP - Fishing effort reconstruction - FUTURE_fishing_effort.csv", skip = 2) %>% 
  dplyr::select(YEAR, COM_n, COM_w)

plot(future_boat_days$YEAR, future_boat_days$COM_n, 
     pch = 19, col = colour_palette[4], xlim = c(year_start, year_end), ylim = c(0, max(future_boat_days$COM_n)), 
     xlab = "Year", ylab = "Boat Days", main = "Original + Fake Boat Days")


#### 4.N.b. split yearly fishing effort by month ------------------------------

# obtain the monthly distribution of fishing from the literature
seasonal_multipliers <- c( # see figure 21c in Ryan et al. 2022 - THIS IS FOR RECREATIONAL FISHING BUT CANT FIND COMMERCIAL EQUIVALENT
  "01" = 0.14, "02" = 0.081, "03" = 0.097, "04" = 0.081,
  "05" = 0.033, "06" = 0.033, "07" = 0.033, "08" = 0.033,
  "09" = 0.033, "10" = 0.065, "11" = 0.11, "12" = 0.26
)
seasonal_multipliers <- seasonal_multipliers / sum(seasonal_multipliers) # standardise so it adds up to 1
barplot(seasonal_multipliers, col = colour_palette[6], main = "distribution of yearly \nboat fishing effort by month in % \n(deduced from Ryan et al. 2022, fig. 21c)")

# add monthly distribution back to the timeseries
boat_effort_n <- expand.grid(
  year = min(future_boat_days$YEAR):max(future_boat_days$YEAR),
  month = sprintf("%02d", 1:12)
) %>%
  arrange(year, month) %>%
  mutate(
    annual_boat_days = rep(future_boat_days$COM_n, each = 12),
    monthly_effort = annual_boat_days * seasonal_multipliers[month]
  ) %>%
  dplyr::select(year, month, monthly_effort)

# sanity check station
ggplot(boat_effort_n, aes(x = as.Date(paste(year, month, "01", sep = "-")), y = monthly_effort)) +
  geom_line(color = colour_palette[4]) +
  labs(title = paste0("Monthly Commercial Fishing Effort ", year_start, "-", year_end),
       x = "Date", y = "Monthly Boat Days") +
  theme_minimal() +
  geom_smooth(color = colour_palette[5])


# Proportion of each month's contribution to yearly boat days
boat_month_prop <- boat_effort_n %>% 
  group_by(year) %>% 
  mutate(year_sum = sum(monthly_effort)) %>%
  mutate(month_prop = monthly_effort/year_sum) %>% 
  dplyr::select(-year_sum)

boat_month_prop <- boat_month_prop %>% 
  group_by(month) %>% 
  mutate(ave_month_prop = mean(month_prop))
prop_month_ave <- boat_month_prop[1:12, c(2, 5)]

saveRDS(prop_month_ave, paste0(folder_name, "/SB_", scenario_name, "_commercial_metro_prop_month_ave.rds")) # charlotte's 'Average_Monthly_Effort"


#### 4.N.c. distribute monthly effort into boat ramps -------------------------

BR_n <- st_read(file_boat_ramps_n) %>% 
  st_transform(4283) %>%
  st_make_valid() %>%
  dplyr::filter(!is.na(ABS_name)) %>% # select the ramps that are mentioned in the ABS reports
  mutate(build_year = 1900, # all commercial ramps start in 1900 cuz some of the build years dont make sense
         build_mnth = 1,     # assume Jan if unknown
         norm_popularity = as.numeric(com_prop) / sum(as.numeric(com_prop), na.rm = TRUE)) %>%
  glimpse()

plot(water$geometry); plot(BR_n$geometry, col = colour_palette[5], pch = 16, cex = 1, add = TRUE)


# distribute effort across all ramps, across time
boat_effort_n <- boat_effort_n %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) # make a column with year and month of ramp build
ramp_effort_n <- expand.grid(
  ramp_index = 1:nrow(BR_n),
  date = boat_effort_n$date
) %>%
  mutate(
    build_date = as.Date(paste(BR_n$build_year[ramp_index], BR_n$build_mnth[ramp_index], "01", sep = "-")),
    norm_popularity = BR_n$norm_popularity[ramp_index],
    ramp_name = BR_n$name[ramp_index]
  ) %>%
  filter(date >= build_date) %>%
  mutate(
    months_since_build = interval(build_date, date) %/% months(1),
    logistic_growth = 1 / (1 + exp(-0.1 * (months_since_build - 60))),
    ramp_weight = norm_popularity * logistic_growth
  )

# Merge in monthly effort
ramp_effort_n <- ramp_effort_n %>%
  left_join(boat_effort_n, by = "date") %>%
  group_by(date) %>%
  mutate(
    total_weight = sum(ramp_weight),
    adjusted_effort = ifelse(total_weight > 0, monthly_effort * (ramp_weight / total_weight), 0),
    year = year(date),
    month = month(date)
  ) %>%
  ungroup() %>%
  dplyr::select(year, month, boat_ramp = ramp_name, adjusted_effort)

# Check total per month equals monthly_effort
check_totals <- ramp_effort_n %>%
  group_by(year, month) %>%
  summarise(total_effort = sum(adjusted_effort), .groups = "drop") %>%
  left_join(boat_effort_n %>% mutate(year = year(date), month = month(date)), by = c("year", "month")) %>%
  mutate(diff = abs(total_effort - monthly_effort))
summary(check_totals$diff)  # should be near zero

# check each boat ramp is populated correctly
ggplot(ramp_effort_n, aes(x=year, y=adjusted_effort)) +
  geom_line(color = colour_palette[5], lwd = 1) +
  facet_wrap(~boat_ramp, ncol = 2)



### 4.W. Wadandi (Southwest) --------------------------------------------------

#### 4.W.a. enter the overall effort values (boat days) -----------------------

# obtain boat days from forecasting (see google slides on reconstruction)
future_boat_days <- read.csv("data/input_data/YIJARUP - Fishing effort reconstruction - FUTURE_fishing_effort.csv", skip = 2) %>% 
  dplyr::select(YEAR, COM_n, COM_w)

plot(future_boat_days$YEAR, future_boat_days$COM_w, 
     pch = 19, col = colour_palette[4], xlim = c(year_start, year_end), ylim = c(0, max(future_boat_days$COM_w)), 
     xlab = "Year", ylab = "Boat Days", main = "Original + Fake Boat Days")


#### 4.W.b. split yearly fishing effort by month ------------------------------

# obtain the monthly distribution of fishing from the literature
seasonal_multipliers <- c( # see figure 21c in Ryan et al. 2022 - THIS IS FOR REC FISHING BUT CANT FIND COMM EQUIVALENT
  "01" = 0.14, "02" = 0.081, "03" = 0.097, "04" = 0.081,
  "05" = 0.033, "06" = 0.033, "07" = 0.033, "08" = 0.033,
  "09" = 0.033, "10" = 0.065, "11" = 0.11, "12" = 0.26
)
seasonal_multipliers <- seasonal_multipliers / sum(seasonal_multipliers) # standardise so it adds up to 1
barplot(seasonal_multipliers, col = colour_palette[6], main = "distribution of yearly \nboat fishing effort by month in % \n(deduced from Ryan et al. 2022, fig. 21c)")

# add monthly distribution back to the timeseries
boat_effort_w <- expand.grid(
  year = min(future_boat_days$YEAR):max(future_boat_days$YEAR),
  month = sprintf("%02d", 1:12)
) %>%
  arrange(year, month) %>%
  mutate(
    annual_boat_days = rep(future_boat_days$COM_w, each = 12),
    monthly_effort = annual_boat_days * seasonal_multipliers[month]
  ) %>%
  dplyr::select(year, month, monthly_effort)

# check
ggplot(boat_effort_w, aes(x = as.Date(paste(year, month, "01", sep = "-")), y = monthly_effort)) +
  geom_line(color = colour_palette[4]) +
  labs(title = paste0("Monthly Boat Fishing Effort ", year_start, "-", year_end),
       x = "Date", y = "Monthly Boat Days") +
  theme_minimal() +
  geom_smooth(color = colour_palette[5])

# Proportion of each month's contribution to yearly boat days - as of 21.01.2026 this is the same in W as in N.
boat_month_prop <- boat_effort_w %>% 
  group_by(year) %>% 
  mutate(year_sum = sum(monthly_effort)) %>%
  mutate(month_prop = monthly_effort/year_sum) %>% 
  dplyr::select(-year_sum)

boat_month_prop <- boat_month_prop %>% 
  group_by(month) %>% 
  mutate(ave_month_prop = mean(month_prop, na.rm = TRUE))
prop_month_ave <- boat_month_prop[1:12, c(2, 5)]

plot(prop_month_ave)
saveRDS(prop_month_ave, paste0(folder_name, "/SB_", scenario_name, "_commercial_wadandi_prop_month_ave.rds")) # charlotte's 'Average_Monthly_Effort"

# sanity check station
plot(future_boat_days$YEAR, future_boat_days$COM_w, type = "l", col = colour_palette[5], lwd = 4,
     main = "Commercial Boat Days in \nWadandi Country (dashed) and Metro (solid)", 
     xlab = "Year", ylab = "Boat Days")
lines(future_boat_days$YEAR, future_boat_days$COM_n, col = colour_palette[4], lwd = 4)
legend("topleft", legend = c("metro", "wadandi"), col = c(colour_palette[4], colour_palette[5]), lwd = 4)


#### 4.W.c. distribute monthly effort into boat ramps -------------------------

# effort by boat_ramp (wadandi)
BR_w <- st_read(file_boat_ramps_w) %>% 
  st_transform(4283) %>%
  st_make_valid() %>%
  mutate(build_year = as.numeric(build_year),
         build_year = ifelse(is.na(build_year), year_start, build_year), # fill in missing dates with the start year
         build_mnth = ifelse(is.na(build_mnth), 1, build_mnth),     # assume Jan if unknown
         norm_popularity = com_prop / sum(com_prop, na.rm = TRUE)) %>%
  glimpse()

# for commercial fishing, only a few boat ramps are used (see Andrea's historical fishing resources, fishing localities in ABS stats)
unique(BR_w$name)
BR_w <- BR_w %>% 
  dplyr::filter(BR_w$name %in% c("SC_Augusta_Ellis_St_Jetty", "WC_Gnarabup", "WC_Hamelin_Bay",
                                 "GB_Quindalup", "GB_Eagle_Bay", "GB_Bunbury_Stirling_St", "GB_Busselton_Georgette_Street"))
plot(water$geometry); plot(BR_w$geometry, col = colour_palette[5], pch = 16, cex = 1, add = TRUE)


# Distribute effort across all ramps, across time
boat_effort_w <- boat_effort_w %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) # make a column with year and month of ramp build
ramp_effort_w <- expand.grid(
  ramp_index = 1:nrow(BR_w),
  date = boat_effort_w$date
) %>%
  mutate(
    build_date = as.Date(paste(BR_w$build_year[ramp_index], BR_w$build_mnth[ramp_index], "01", sep = "-")),
    norm_popularity = BR_w$norm_popularity[ramp_index],
    ramp_name = BR_w$name[ramp_index]
  ) %>%
  mutate(
    months_since_build = interval(build_date, date) %/% months(1),
    logistic_growth = 1 / (1 + exp(-0.1 * (months_since_build - 60))),
    ramp_weight = norm_popularity * logistic_growth
  )

# Merge in monthly effort
ramp_effort_w <- ramp_effort_w %>%
  left_join(boat_effort_w, by = "date") %>%
  group_by(date) %>%
  mutate(
    total_weight = sum(ramp_weight),
    adjusted_effort = ifelse(total_weight > 0, monthly_effort * (ramp_weight / total_weight), 0),
    year = year(date),
    month = month(date)
  ) %>%
  ungroup() %>%
  dplyr::select(year, month, boat_ramp = ramp_name, adjusted_effort)
#saveRDS(ramp_effort_df, "data/output_data/04A_commercial_ramp_effort.rds")

# Check total per month equals monthly_effort
check_totals <- ramp_effort_w %>%
  group_by(year, month) %>%
  summarise(total_effort = sum(adjusted_effort), .groups = "drop") %>%
  left_join(boat_effort_w %>% mutate(year = year(date), month = month(date)), by = c("year", "month")) %>%
  mutate(diff = abs(total_effort - monthly_effort))
summary(check_totals$diff)  # should be near zero

# check each boat ramp is populated correctly
ggplot(ramp_effort_w, aes(x=year, y=adjusted_effort)) +
  geom_line(color = colour_palette[5], lwd = 1) +
  facet_wrap(~boat_ramp, ncol = 2)


# okay combine into a month x access point x year cube

ramp_effort_n
ramp_effort_w
access_point_effort <- rbind(ramp_effort_n, ramp_effort_w)

# cut up this big unwieldy dataframe into a cube
access_point_effort <- access_point_effort %>% 
  pivot_wider(names_from = boat_ramp,
              values_from = adjusted_effort) %>% 
  as.data.frame()
access_point_effort <- access_point_effort %>%
  dplyr::select(-year, -month) %>%
  {abind::abind(split(., access_point_effort$year), along = 3)}

# rename dimension names so they're consistent with the rest of the objects used in the function.
dimnames(access_point_effort)[[1]] <- c(1:12)
dimnames(access_point_effort)[[3]] <- c(1:length(dimnames(access_point_effort)[[3]]))


# check that utility and access_point_effort have the same order of access point
names(access_dist %>% dplyr::select(!ID)) == dimnames(access_point_effort)[[2]]
dimnames(access_point_effort)[[2]] <- c(1:length(dimnames(access_point_effort)[[2]]))

# utility # watch out that the order of ramp effort and of utility match: otherwise you're assigning the effort to the wrong ramp in the function


# 5. create a list to use in the C++ function ---------------------------------

# little details to include in the final list
cell_area_m2 <- water$cell_area # just cell area
coef_values <- tibble(log_utility = 1,
                      expected_catch = 1,
                      expected_catch_sq = 1,
                      log_cell_area = 1
) # coefficients from Matt's 2022 paper. they weigh the relative importance of each for the distribution of effort. the paper was for rec fishing so these are all 1.


fishing_info_list <- list(catchability, utility, access_point_effort, cell_area_m2, coef_values)
names(fishing_info_list) <- c("catchability", "utility", "fishing_days", "cell_area_m2", "coef_values")
saveRDS(fishing_info_list, paste0(folder_name, "/SB_", scenario_name, "_commercial_fishing_info.rds"))


# SHORE-RECREATIONAL ==========================================================

current_fleet <- "shore_rec"

# 1. travel cost = distance from access points --------------------------------

water <- readRDS(file_water)
mapview::mapview(water[,!is.list(water)])

AP_n <- st_read(file_carpark_n) %>% 
  st_transform(common_crs) %>%
  rename(name = rough_area) %>% 
  mutate(id = 1:nrow(.),
         name = paste0(name, "_", id)) %>% # many access points have the same name, make them unique otherwise later it's hard to deal with in the access_dist df.
  st_make_valid() %>%
  glimpse()

AP_w <- st_read(file_carpark_w) %>% 
  st_transform(common_crs) %>%
  rename(name = area) %>%  
  st_make_valid() %>%
  glimpse()

plot(water$geometry); plot(AP_w$geometry, col = colour_palette[5], pch = 16, cex = 1, add = TRUE)

glimpse(AP_n)
glimpse(AP_w)

AP <- rbind(AP_n %>% dplyr::select(name, year_start, geometry) %>% mutate(region = "north"), 
            AP_w %>% dplyr::select(name, year_start, geometry) %>% mutate(region = "wadandi")) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>% 
  glimpse()

plot(water$geometry); plot(AP$geometry, col = colour_palette[6], pch = 16, add = TRUE)

network <- st_read(file_network); plot(network$geometry)
AP <- st_as_sf(AP) %>% st_transform(common_crs)

# find the centre of each grid cell
sf::sf_use_s2(FALSE)
centroids <- st_centroid(st_make_valid(water))
sf::sf_use_s2(TRUE)
points <- as.data.frame(st_coordinates(centroids))%>% # the points start at the bottom left and then work their way their way right
  mutate(ID = row_number()) 
points_sf <- st_as_sf(points, coords = c("X", "Y"))
st_crs(points_sf) <- common_crs


network <- as_sfnetwork(network, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

net <- activate(network, "nodes")
st_crs(net)
net <- net %>% st_transform(common_crs)

# measure the distance from access points to cell centroids
network_matrix <- st_network_cost(net, from = AP, to = points_sf) / 1000 #in km
dim(network_matrix) # number of ramps x number of cells

glimpse(network_matrix)
access_dist <- as.data.frame(t(network_matrix))
colnames(access_dist) <- AP$name
access_dist$ID <- water$ID
head(access_dist) # this gives us each cell's distance to the access points

access_dist <- units::drop_units(access_dist) %>% mutate_at(vars(1:(ncol(access_dist)-1)), ~replace(., .>2, NA)) # cells over 3 km from the access point aren't accessible to shore fishers.
access_dist[water$ID[!water$type %in% c("shore_wadandi", "shore_north")], 1:(ncol(access_dist)-1)] <- NA # offshore cells are not accessible to shore fishers

## sanity check station
access_dist2 <- access_dist %>% mutate(ID = centroids$ID)
water_dist_long <- water %>% left_join(access_dist2, by = "ID") %>% pivot_longer(cols = AP$name, names_to = "AP", values_to = "Distance_km")
test <- water %>% mutate(test = access_dist2[,100]) %>%   dplyr::select(test, geometry) # keep only what you need

mapview::mapview(test, zcol = "test") + mapview::mapview(AP[100,])


ggplot(water_dist_long %>% dplyr::filter(AP %in% unique(water_dist_long$AP)[1:10])) + 
  geom_sf(aes(fill = as.numeric(Distance_km)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  facet_wrap(~ AP, ncol = 5) + 
  labs(title = "Cell distance to access points", fill = "Distance (km)") +
  theme_minimal()


# 2. calculate cell utility ---------------------------------------------------

## 2.a add in build year ------------------------------------------------------

# there are so many carparks, some of which are quite recent. we'll zero the utility of access points if they haven't been built yet.

build_year <- array(0,
                    dim = c(ncol(access_dist)-1, 
                            n_years_tot)
)

for (ACCESS in 1:nrow(AP)) {
  for (YEAR in 1:n_years_tot) {
    yr <- year_start + YEAR - 1
    build_year[ACCESS, YEAR] <- ifelse(AP$year_start[ACCESS] <= yr, 1, NA)
  }
}


## 2.b cell utility = 1 / ((built * travel cost) + 1) -------------------------

glimpse(access_dist)

# we're making cell x access_p x year

utility <- array(NA,
                 dim = c(nrow(access_dist), 
                         ncol(access_dist)-1, 
                         n_years_tot)
)

for (ACCESS in 1:ncol(utility)) {
  for (YEAR in 1:n_years_tot) {
    # utility = cell x access x year
    # build_year = access x year
    # access_dist = cell x access    
    utility[, ACCESS, YEAR] <- 1 / (build_year[ACCESS, YEAR] * (as.numeric(access_dist[, ACCESS]) + 1))
  }
}
utility[,1:20,10]
utility[!is.finite(utility)] <- 0
utility[is.na(utility)] <- 0



## sanity check - is utility higher for areas close to access points but far from shore)
test_access = 10
test_year = 10
test <- water %>% dplyr::select(!where(is.list)) %>% mutate(test = utility[, test_access, test_year])
mapview::mapview(test, zcol = "test") + mapview::mapview(AP[test_access,])

ggplot() +
  geom_sf(data = water %>% dplyr::select(!where(is.list)) %>% mutate(test = utility[, test_access, test_year]), 
          aes(fill = test), 
          col = NA) +
  geom_sf(data = AP[test_access,] %>% st_set_crs(common_crs) %>% st_transform(st_crs(water)), col = "red") +
  scale_fill_gradientn(colours = colour_palette[6:4])
# ayoooo im a geniuuus


# 3. fishability --------------------------------------------------------------

# fishability refers to whether this cell in this month and this year is open to fishing.
# each fleet has a different fishability because not all closures/NTZs apply to everyone.

# fishability is between 0 and 1: 
# 0 not fishable for the entire month, 
# 1 fishable for the entire month
# decimals mean either:
## a. fishable for part of the month (in the case of temporal closures), or 
## b. partially fishable in area (related to the % of the cell that's within fishable depth)

# for the final cell x month x year cube, we need: 
## depth fishability (better engines over time allow deeper fishing)
## temporal fishability (cell x month x year)
## cell area

## we will calculate: (depth fishability x temporal fishability * cell area) / total fishable area

## 3.b Temporal fishability ---------------------------------------------------

# this is just a cell x month x year cube that tells you whether there's a closure or not in the cell.
# some cells are closed for only part of the month (ex. temporal closure that starts midway through the month) 
# so that cell would have the value 0.5 that month.

glimpse(water)

## in this section (1.) we calculate the fishable area of the cells in each time step (1.a), then divide it by the fishable area of the whole grid at each time step (1.b)
## which gives us each the portion of Catchability of each cell.
## because the area that is catchable changes (with spatial and temporal restrictions), we have to calculate the grid's catchability for each month of each year.


### 3.b.a. find the fishable area of each cell in each time step --------------

water <- readRDS(file_water) %>% filter(!is.na(ID))
water <- merge(water, st_drop_geometry(n_restrictions), by = "ID")

NCELL <- nrow(water)

# identify the important cells
temporal_cells <- water$ID[which(water$TC_status == TRUE | water$n_TC_sttus == TRUE)] # we're selecting cells that are closed currently

# setup the arrays to calculate things.
water_area <- array(0, dim = c(NCELL, 12, (year_end-year_start+1))) # cells x months x years
water_area[, 1:12, 1] <- water$cell_area # for all months of the first year, the catchable area is the cell's area.

SC_check <- function(sc_d_string, fleet, current_year, current_month) {
  
  if (fleet == TRUE) return(FALSE) # if the closure applies to a different fleet than commercial (TRUE = fleet can fish, FALSE = fleet can't fish)
  
  if (is.na(sc_d_string) || sc_d_string == "") return(FALSE) # if there is no spatial restriction, we are not in a closure, so FALSE
  
  ranges <- trimws(strsplit(sc_d_string, ",")[[1]]) # split the different closure periods
  
  for (rng in ranges) {
    
    current_ym <- current_year * 12 + current_month
    
    parts <- strsplit(rng, "to")[[1]] # obtain the parts of the closure periods (start and end)
    start_parts <- strsplit(trimws(parts[1]), "-")[[1]]
    start_month <- as.integer(start_parts[1]); start_year <- as.integer(start_parts[2])
    start_ym   <- start_year  * 12 + start_month    
    
    if (length(parts) == 1) { # if there is only one date in the SC_d, the closure is valid for the entire simulation (no end date), so TRUE (if we've entered the closure start)
      if (current_ym >= start_ym) return(TRUE)
      next
    } 
    # otherwise get the end dates.
    end_parts   <- strsplit(trimws(parts[2]), "-")[[1]]
    end_month   <- as.integer(end_parts[1]);   end_year   <- as.integer(end_parts[2])
    end_ym     <- end_year    * 12 + end_month
    
    if (current_ym >= start_ym && current_ym <= end_ym) return(TRUE) # if you're in a closure period, you're in a closure, so TRUE
  }
  return(FALSE)
} # this is a function to make sense of the start/ends of the spatial closures.

# the loop goes as follows: for each time step, 
# check whether there's a spatial restriction that year, if so the fishable area is 0
# then, if the cell has a temporal closure, restrict fishable area as needed.
current_fleet_n <- "n_sr"

for (YEAR in 1:dim(water_area)[3]) {
  
  current_year <- year_start + YEAR - 1
  print(current_year)
  
  for (MONTH in 1:dim(water_area)[2]) {
    
    # spatial closures
    
    water_area[, MONTH, YEAR] <- water$cell_area
    
    fishing_allowed <- mapply(SC_check, # here, we check whether each cell is open to fishing in this MONTH and YEAR. FALSE = not open, TRUE = open
                              sc_d_string  = water$n_SC_d,
                              fleet        = water$n_com,
                              MoreArgs     = list(current_year  = current_year,
                                                  current_month = MONTH))
    
    restricted_cells <- which(
      !is.na(water[[current_fleet_n]]) & water[[current_fleet_n]] == FALSE &
        fishing_allowed
    )
    water_area[restricted_cells, MONTH, YEAR] <- 0 # ...the cell is not fishable
    
    # temporal closures
    for (CELL in temporal_cells) {
      
      # find the relevant temporal restrictions
      TC_years_list <- water$n_TC_d[[CELL]] # years the closure has changed
      TC_months_list <- water$n_TC_m[[CELL]] # months restricted
      TC_perc_list <- water$n_TC_p[[CELL]] # how much they are restricted
      
      current_TC <- which(TC_years_list <= current_year) # only select the temporal closures that exist at the current point
      
      if (length(current_TC) > 0) { # if there are temporal restrictions at this YEAR, restrict the fishable area accordingly
        current_TC <- max(current_TC)
        
        TC_months_current <- as.integer(strsplit(TC_months_list[current_TC], "-")[[1]])
        TC_perc_current   <- as.numeric(strsplit(TC_perc_list[current_TC], "-")[[1]])
        
        if (MONTH %in% TC_months_current) {
          month_idx <- match(MONTH, TC_months_current)
          water_area[CELL, MONTH, YEAR] <- water_area[CELL, MONTH, YEAR] * TC_perc_current[month_idx]
        }
      }
    }
  }
} # TAKES FOREVER - this loop calculates for every month and every year, the fishable area of each cell.

water_area[water$ID[!water$type %in% c("shore_wadandi", "shore_north")],,] <- 0 # offshore cells are not open to shore-fishers

### sanity check station

test_cell <- 200
test_month <- 10
test_year <- 10

# reference cell numbers as of 13.01.2026:
## cockburn sound cell (temporal closure): 1
## SWC NTZ cell: 200
## random fished cell: 1000

# restrictions in this cell should be:
glimpse(st_drop_geometry(water[test_cell, c("TC_status", "SC_status", "commercial", "SC_restriction_date", "TC_restriction_date", "TC_restriction_months", "TC_restriction_perc_fished")]))

# check that it is correct
cat("In month", test_month, "of year", (year_start+test_year-1),
    ", the cell is", ifelse(water_area[test_cell, test_month, test_year]>0, "fishable", "NOT fishable"), "for shore rec fishers", "Fishable area: ", water_area[test_cell, test_month, test_year]/1e06, "km2")

## sanity check station 
test_year = 10
test <- water_area[, 1, test_year]
mapview::mapview(water %>% dplyr::select(!where(is.list)) %>% mutate(test = test), zcol = "test")


ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()


### 3.b.b calculate catchability ----------------------------------------------

glimpse(water_area)

fishable_area <- water_area

## sanity check station 
test_year = 10
test <- fishable_area[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()

# now divide by the sum of fishable area.
glimpse(fishable_area)
fishable_area_sum <- colSums(fishable_area) # month x year matrix

catchability <- sweep(fishable_area, c(2, 3), fishable_area_sum, FUN = "/")
dim(fishable_area)
dim(fishable_area_sum)


## sanity check station 
test_year = 20
test <- catchability[, 1, test_year]
mapview::mapview(water %>% dplyr::select(!where(is.list)) %>% mutate(test = test), zcol = "test")

ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()


# 4. set up fishing days values -----------------------------------------------

# we need to know how much fishing occurs in the region, so that the function knows how much to distribute.
# for commercial fishing, we are using a variety of sources to reconstruct the trends (a),
# the sources are split between the Metropolitan area effort (North, 4.N) and the Southwest area effort (Wadandi, 4.W)
# we will (b.) split the yearly effort into months, then (c.) split the monthly effort into access points (boat ramps).


## 4.a. enter the overall effort values (boat days) ---------------------------

# obtain boat days from forecasting (see google slides on reconstruction)
future_boat_days <- read.csv("data/input_data/YIJARUP - Fishing effort reconstruction - FUTURE_fishing_effort.csv", skip = 2) %>% 
  dplyr::select(YEAR, SHORE_REC_fishing_days)

plot(future_boat_days$YEAR, future_boat_days$SHORE_REC_fishing_days, 
     pch = 19, col = colour_palette[4], xlim = c(year_start, year_end), ylim = c(0, max(future_boat_days$SHORE_REC_fishing_days)), 
     xlab = "Year", ylab = "Shore fishing Days", main = "Forecasted shore fishing days")


## 4.b. split yearly fishing effort by month ----------------------------------

# obtain the monthly distribution of fishing from the literature
seasonal_multipliers <- c( # see figure 21c in Ryan et al. 2022 - THIS IS FOR RECREATIONAL FISHING BUT CANT FIND COMMERCIAL EQUIVALENT
  "01" = 0.14, "02" = 0.081, "03" = 0.097, "04" = 0.081,
  "05" = 0.033, "06" = 0.033, "07" = 0.033, "08" = 0.033,
  "09" = 0.033, "10" = 0.065, "11" = 0.11, "12" = 0.26
)
seasonal_multipliers <- seasonal_multipliers / sum(seasonal_multipliers) # standardise so it adds up to 1
barplot(seasonal_multipliers, col = colour_palette[6], main = "distribution of yearly \nboat fishing effort by month in % \n(deduced from Ryan et al. 2022, fig. 21c)")

# add monthly distribution back to the timeseries
shore_effort <- expand.grid(
  year = year_start:year_end,
  month = sprintf("%02d", 1:12)
) %>%
  arrange(year, month) %>%
  mutate(
    annual_shore_days = rep(future_boat_days$SHORE_REC_fishing_days, each = 12),
    monthly_effort = annual_shore_days * seasonal_multipliers[month]
  ) %>%
  dplyr::select(year, month, monthly_effort)

# sanity check station
ggplot(shore_effort, aes(x = as.Date(paste(year, month, "01", sep = "-")), y = monthly_effort)) +
  geom_line(color = colour_palette[4]) +
  labs(title = paste0("Monthly Commercial Fishing Effort ", year_start, "-", year_end),
       x = "Date", y = "Monthly Boat Days") +
  theme_minimal() +
  geom_smooth(color = colour_palette[5])


## 4.c. distribute monthly effort into boat ramps -----------------------------

# carparks
AP_n <- st_read(file_carpark_n) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  mutate(build_mnth = 1 # assume Jan if unknown
  ) %>%
  dplyr::select(geometry, year_start, build_mnth) %>% 
  glimpse()

AP_w <- st_read(file_carpark_w) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  mutate(build_mnth = month_strt) %>%
  dplyr::select(geometry, year_start, build_mnth) %>% 
  glimpse()

AP <- rbind(AP_n, AP_w)


plot(water$geometry); plot(AP$geometry, col = colour_palette[5], pch = 16, cex = 1, add = TRUE)

glimpse(shore_effort)
logistic_growth <- function(t, k = 0.3, t0 = 10) {
  1 / (1 + exp(-k * (t - t0)))
}
# we want a month x access point x year
distributed_effort <- array(0, 
                            dim = c(12, 
                                    nrow(AP),
                                    n_years_tot))

for (YEAR in 1:n_years_tot) {
  
  built_logic <- AP$year_start <= (year_start + YEAR - 1) # find which carparks are built in this year
  built_carparks <- which(built_logic) # and their index
  
  # for logistic growth of effort
  years_since_build <- ifelse((year_start + YEAR -1) - AP$year_start[built_carparks] < 0, 0, (year_start + YEAR -1) - AP$year_start[built_carparks])
  weights <- logistic_growth(years_since_build, k = 0.3, t0 = 10)
  
  for (MONTH in 1:12) {
    
    current_effort <- shore_effort %>% 
      dplyr::filter(year == (year_start + YEAR - 1), as.numeric(month) == MONTH) %>% 
      dplyr::pull(monthly_effort) # extract how many fishing days to distribute
    
    distributed_effort[MONTH, built_carparks, YEAR] <- (current_effort * (weights / sum(weights)))   # and distribute them
    
  }
}


plot(x = 1:n_years_tot, y = distributed_effort[1,1,])
lines(x = 1:n_years_tot, y = distributed_effort[1,50,])
points(x = 1:n_years_tot, y = distributed_effort[1,300,], col = "red", pch = 3)



# okay combine into a month x access point x year cube

access_point_effort <- distributed_effort

# rename dimension names so they're consistent with the rest of the objects used in the function.
dimnames(access_point_effort)[[1]] <- c(1:12)
dimnames(access_point_effort)[[3]] <- c(1:dim(access_point_effort)[3])


# check that utility and access_point_effort have the same order of access point
dimnames(access_point_effort)[[2]] <- names(access_dist %>% dplyr::select(!ID))

# utility # watch out that the order of ramp effort and of utility match: otherwise you're assigning the effort to the wrong ramp in the function

# 5. create a list to use in the C++ function ---------------------------------

# little details to include in the final list
cell_area_m2 <- water$cell_area # just cell area
coef_values <- tibble(log_utility = 1,
                      expected_catch = 1,
                      expected_catch_sq = 1,
                      log_cell_area = 1
) # coefficients from Matt's 2022 paper. they weigh the relative importance of each for the distribution of effort. the paper was for rec fishing so these are all 1.


fishing_info_list <- list(catchability, utility, access_point_effort, cell_area_m2, coef_values)
names(fishing_info_list) <- c("catchability", "utility", "fishing_days", "cell_area_m2", "coef_values")
saveRDS(fishing_info_list, paste0(folder_name, "/SB_", scenario_name, "_shore_rec_fishing_info.rds"))


# BOAT-RECREATIONAL ===========================================================

current_fleet <- "b_rec"

# 1. travel cost = distance from access points * fuel price -------------------

## 1.a fuel price -------------------------------------------------------------

# obtain fuel prices from the literature
fuel_price <- read.csv("data/input_data/fuel_prices WA population - forecast.csv") # data from various reports from government and international agencies.

plot(fuel_price$petrol_price ~ fuel_price$year, type = "l") # plot check

fuel_price_clean <- fuel_price$petrol_price


## 1.b distance from access points --------------------------------------------

water <- readRDS(file_water)

BR_n <- st_read(file_boat_ramps_n) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  glimpse()


BR_w <- st_read(file_boat_ramps_w) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  glimpse()

glimpse(BR_n)
glimpse(BR_w)

BR <- rbind(BR_n %>% dplyr::select(name, build_year, geometry) %>% mutate(region = "north"), 
            BR_w %>% dplyr::select(name, build_year, geometry) %>% mutate(region = "wadandi")) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>% 
  glimpse()

plot(water$geometry); plot(BR$geometry, col = colour_palette[6], pch = 16, add = TRUE)

network <- st_read(file_network) %>% 
  st_transform(common_crs); plot(network$geometry)
BR <- st_as_sf(BR); BR <- BR %>% st_transform(common_crs) %>% dplyr::filter(!st_is_empty(.))

# find the centre of each grid cell
sf::sf_use_s2(FALSE)
centroids <- st_centroid(st_make_valid(water))
sf::sf_use_s2(TRUE)
points <- as.data.frame(st_coordinates(centroids))%>% # the points start at the bottom left and then work their way their way right
  mutate(ID = row_number()) 
points_sf <- st_as_sf(points, coords = c("X", "Y")) 
st_crs(points_sf) <- common_crs

network <- as_sfnetwork(network, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

net <- activate(network, "nodes")
net <- net %>% st_transform(common_crs) 

# measure the distance from access points to cell centroids
network_matrix <- st_network_cost(net, from = BR, to = points_sf)
dim(network_matrix) # number of ramps x number of cells

glimpse(network_matrix)
access_dist <- as.data.frame(t(network_matrix))
colnames(access_dist) <- BR$name
access_dist$ID <- water$ID
access_dist <- access_dist / 1000
head(access_dist) # this gives us each cell's distance to the access points

## sanity check station
access_dist <- access_dist %>% mutate(ID = centroids$ID)
water_dist_long <- water %>% left_join(access_dist, by = "ID") %>% pivot_longer(cols = BR$name, names_to = "Ramp", values_to = "Distance_km")
ggplot(water_dist_long %>% dplyr::filter(Ramp %in% unique(water_dist_long$Ramp)[1:10])) + 
  geom_sf(aes(fill = as.numeric(Distance_km)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  facet_wrap(~ Ramp, ncol = 5) + 
  labs(title = "Cell distance to access points", fill = "Distance (km)") +
  theme_minimal()


## 1.c calculate travel cost --------------------------------------------------

glimpse(access_dist) # cell x access_point
glimpse(fuel_price_clean) # vec year

# figure out, for each ramp, the travel cost (cell x access_p x years)
dist_matrix <- as.matrix(access_dist[, -ncol(access_dist)])  # drop ID column

travel_cost <- array(rep(dist_matrix, times = length(fuel_price_clean)),
                     dim = c(nrow(access_dist), 
                             ncol(access_dist)-1, 
                             length(fuel_price_clean))
)


travel_cost <- sweep(
  travel_cost,
  MARGIN = 3,
  STATS  = fuel_price_clean,
  FUN    = "*"
) # calculate the product (fun = "*") of access_dist using fuel price (stats =) across years (margin = 3), 

## sanity check
par(mfrow = c(1,1))
test1 <- travel_cost[1,1,] # cell 1 from access_p 1
test2 <- travel_cost[1,2,] # cell 1 from access_p 2

plot(test2, type = "l", col = "red")
lines(test1, type = "l")


# 2. calculate cell utility ---------------------------------------------------

## 2.a add in build year ------------------------------------------------------

# there are so many access points, some of which are quite recent. we'll zero the utility of access points if they haven't been built yet.

build_year <- array(0,
                    dim = c(ncol(access_dist)-1, 
                            n_years_tot)
)

for (ACCESS in 1:nrow(BR)) {
  for (YEAR in 1:n_years_tot) {
    yr <- year_start + YEAR - 1
    build_year[ACCESS, YEAR] <- ifelse(BR$build_year[ACCESS] <= yr, 1, NA)
  }
}


## 2.b cell utility = 1 / ((built * travel cost) + 1) -------------------------

glimpse(access_dist)

# we're making cell x access_p x year

utility <- array(NA,
                 dim = c(nrow(access_dist), 
                         ncol(access_dist)-1, 
                         n_years_tot)
)

for (ACCESS in 1:ncol(utility)) {
  for (YEAR in 1:n_years_tot) {
    # utility = cell x access x year
    # build_year = access x year
    # access_dist = cell x access    
    utility[, ACCESS, YEAR] <- 1 / (build_year[ACCESS, YEAR] * travel_cost[, ACCESS, YEAR] + 1)
  }
}
utility[,1:20,10]
utility[!is.finite(utility)] <- 0
utility[is.na(utility)] <- 0




## sanity check - is utility higher for areas close to access points but far from shore)
test_access = 10
test_year = 10
test <- water %>% dplyr::select(!where(is.list)) %>% mutate(test = log(utility[, test_access, test_year]))
mapview::mapview(test, zcol = "test") + mapview::mapview(BR[test_access,])

ggplot() +
  geom_sf(data = water %>% dplyr::select(!where(is.list)) %>% mutate(test = utility[, test_access, 10]), 
          aes(fill = (test)), 
          col = NA) +
  geom_sf(data = BR[test_access,] %>% st_set_crs(common_crs) %>% st_transform(st_crs(water)), col = "red") +
  scale_fill_gradientn(colours = colour_palette[6:4])
# ayoooo im a geniuuus


# 3. fishability --------------------------------------------------------------

# fishability refers to whether this cell in this month and this year is open to fishing.
# each fleet has a different fishability because not all closures/NTZs apply to everyone.

# fishability is between 0 and 1: 
# 0 not fishable for the entire month, 
# 1 fishable for the entire month
# decimals mean either:
## a. fishable for part of the month (in the case of temporal closures), or 
## b. partially fishable in area (related to the % of the cell that's within fishable depth)

# for the final cell x month x year cube, we need: 
## depth fishability (better engines over time allow deeper fishing)
## temporal fishability (cell x month x year)
## cell area

## we will calculate: (depth fishability x temporal fishability * cell area) / total fishable area

## 3.a depth fishability ------------------------------------------------------

# this section of the fishing effort relates to boats being able to fish further over the years, with bigger boats, more powerful engines etc.
# the model takes this into account by adding a 'fishable depth' into the mix, which goes more and more offshore.
# for commercial fishing, we're assuming that on average, fishers can fish 1.3m more every year. (change as needed in parameters)
NCELL <- nrow(water)

# prepare grid cells
plot(water$geometry)
wa_mask <- st_read(file_wa); wa_mask <- st_transform(wa_mask, crs = crs_raster); wa_mask <- as(wa_mask, "Spatial"); plot(wa_mask, col = "lightgray", add = T)
plot(bbox, add = T)

# set parameters (in 2024, commercial boats could fish in 122.7. we'll continue this 1m/year expansion)
years <- year_start:year_end
min_depth <- 122.7 # 20m max in 1945 + ((2024-1945) * 1m)
increase_year <- 2025 # the year that boats start to go deeper (in the scenarios, this is straightaway)
depth_per_year <- 1 # fishable depth gained per year (m) - arbitrary, but consistent with the historical reconstruction.

# prepare the bathymetry layer
bathy <- terra::rast(file_bathy) %>%
  terra::project(crs_raster) %>% 
  terra::crop(as(water, "Spatial")) %>%
  abs() %>%
  terra::mask(terra::vect(wa_mask), inverse = TRUE); plot(bathy)

thresholds <- min_depth + pmax(0, depth_per_year*(years - increase_year)) # depth limits for each year

fishable_stack <- terra::rast(
  lapply(thresholds, function(th) bathy <= th)
)

names(fishable_stack) <- paste0("year_", years)

fishable_summary <- terra::extract(
  fishable_stack,
  terra::vect(water),
  fun = mean,
  na.rm = TRUE
)

fishable_summary$ID <- water$ID

fishable_depth_cell_year <- fishable_summary %>%
  arrange(ID) %>%
  dplyr::select(-ID) %>%
  as.matrix()

# strip "year_" prefix from colnames
colnames(fishable_depth_cell_year) <- gsub("year_", "", colnames(fishable_depth_cell_year))

fishable_depth_cell_month_year <- array(
  aperm(replicate(12, fishable_depth_cell_year), c(1, 3, 2)),
  dim = c(NCELL, 12, length(years)),
  dimnames = list(NULL))

fishable_depth_cell_month_year[is.nan(fishable_depth_cell_month_year)] <- 1 # replace NaN with 1 (these are shore cells that have depth 0)


## sanity check station 
test_year = 10
test <- fishable_depth_cell_month_year[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()


## 3.b Temporal fishability ---------------------------------------------------

# this is just a cell x month x year cube that tells you whether there's a closure or not in the cell.
# some cells are closed for only part of the month (ex. temporal closure that starts midway through the month) 
# so that cell would have the value 0.5 that month.

glimpse(water)

## in this section (1.) we calculate the fishable area of the cells in each time step (1.a), then divide it by the fishable area of the whole grid at each time step (1.b)
## which gives us each the portion of Catchability of each cell.
## because the area that is catchable changes (with spatial and temporal restrictions), we have to calculate the grid's catchability for each month of each year.


### 3.b.a. find the fishable area of each cell in each time step --------------

water <- readRDS(file_water) %>% filter(!is.na(ID))
water <- merge(water, st_drop_geometry(n_restrictions), by = "ID")

NCELL <- nrow(water)

# identify the important cells
temporal_cells <- water$ID[which(water$TC_status == TRUE | water$n_TC_sttus == TRUE)] # we're selecting cells that are closed currently

# setup the arrays to calculate things.
water_area <- array(0, dim = c(NCELL, 12, (year_end-year_start+1))) # cells x months x years
water_area[, 1:12, 1] <- water$cell_area # for all months of the first year, the catchable area is the cell's area.

SC_check <- function(sc_d_string, fleet, current_year, current_month) {
  
  if (fleet == TRUE) return(FALSE) # if the closure applies to a different fleet than commercial (TRUE = fleet can fish, FALSE = fleet can't fish)
  
  if (is.na(sc_d_string) || sc_d_string == "") return(FALSE) # if there is no spatial restriction, we are not in a closure, so FALSE
  
  ranges <- trimws(strsplit(sc_d_string, ",")[[1]]) # split the different closure periods
  
  for (rng in ranges) {
    
    current_ym <- current_year * 12 + current_month
    
    parts <- strsplit(rng, "to")[[1]] # obtain the parts of the closure periods (start and end)
    start_parts <- strsplit(trimws(parts[1]), "-")[[1]]
    start_month <- as.integer(start_parts[1]); start_year <- as.integer(start_parts[2])
    start_ym   <- start_year  * 12 + start_month    
    
    if (length(parts) == 1) { # if there is only one date in the SC_d, the closure is valid for the entire simulation (no end date), so TRUE (if we've entered the closure start)
      if (current_ym >= start_ym) return(TRUE)
      next
    } 
    # otherwise get the end dates.
    end_parts   <- strsplit(trimws(parts[2]), "-")[[1]]
    end_month   <- as.integer(end_parts[1]);   end_year   <- as.integer(end_parts[2])
    end_ym     <- end_year    * 12 + end_month
    
    if (current_ym >= start_ym && current_ym <= end_ym) return(TRUE) # if you're in a closure period, you're in a closure, so TRUE
  }
  return(FALSE)
} # this is a function to make sense of the start/ends of the spatial closures.


# the loop goes as follows: for each time step, 
# check whether there's a spatial restriction that year, if so the fishable area is 0
# then, if the cell has a temporal closure, restrict fishable area as needed.
current_fleet_n <- "b_rec"

for (YEAR in 1:dim(water_area)[3]) {
  
  current_year <- year_start + YEAR - 1
  print(current_year)
  
  for (MONTH in 1:dim(water_area)[2]) {
    
    # spatial closures
    
    water_area[, MONTH, YEAR] <- water$cell_area
    
    fishing_allowed <- mapply(SC_check, # here, we check whether each cell is open to fishing in this MONTH and YEAR. FALSE = not open, TRUE = open
                              sc_d_string  = water$n_SC_d,
                              fleet        = water$n_br,
                              MoreArgs     = list(current_year  = current_year,
                                                  current_month = MONTH))
    
    restricted_cells <- which(
      !is.na(water[[current_fleet_n]]) & water[[current_fleet_n]] == FALSE &
        fishing_allowed
    )
    water_area[restricted_cells, MONTH, YEAR] <- 0 # ...the cell is not fishable
    
    # temporal closures
    for (CELL in temporal_cells) {
      
      # find the relevant temporal restrictions
      TC_years_list <- water$n_TC_d[[CELL]] # years the closure has changed
      TC_months_list <- water$n_TC_m[[CELL]] # months restricted
      TC_perc_list <- water$n_TC_p[[CELL]] # how much they are restricted
      
      current_TC <- which(TC_years_list <= current_year) # only select the temporal closures that exist at the current point
      
      if (length(current_TC) > 0) { # if there are temporal restrictions at this YEAR, restrict the fishable area accordingly
        current_TC <- max(current_TC)
        
        TC_months_current <- as.integer(strsplit(TC_months_list[current_TC], "-")[[1]])
        TC_perc_current   <- as.numeric(strsplit(TC_perc_list[current_TC], "-")[[1]])
        
        if (MONTH %in% TC_months_current) {
          month_idx <- match(MONTH, TC_months_current)
          water_area[CELL, MONTH, YEAR] <- water_area[CELL, MONTH, YEAR] * TC_perc_current[month_idx]
        }
      }
    }
  }
} # TAKES FOREVER - this loop calculates for every month and every year, the fishable area of each cell.

### sanity check station

test_cell <- 200
test_month <- 10
test_year <- 23

# reference cell numbers as of 13.01.2026:
## cockburn sound cell (temporal closure): 1
## SWC NTZ cell: 200
## random fished cell: 1000

# restrictions in this cell should be:
glimpse(st_drop_geometry(water[test_cell,]))

# check that it is correct
cat("In month", test_month, "of year", (year_start+test_year-1),
    ", the cell is", ifelse(water_area[test_cell, test_month, test_year]>0, "fishable", "NOT fishable"), "for recreational boats.", "Fishable area: ", water_area[test_cell, test_month, test_year]/1e06, "km2")

## sanity check station 
test_year = 20
test <- water_area[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()


### 3.b.b calculate catchability ----------------------------------------------

glimpse(water_area)
glimpse(fishable_depth_cell_month_year)

fishable_area <- water_area * fishable_depth_cell_month_year
dim(water_area)
dim(fishable_depth_cell_month_year)

## sanity check station 
test_year = 20
test <- fishable_area[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()

# now divide by the sum of fishable area.
glimpse(fishable_area)
fishable_area_sum <- colSums(fishable_area) # month x year matrix

catchability <- sweep(fishable_area, c(2, 3), fishable_area_sum, FUN = "/")
dim(fishable_area)
dim(fishable_area_sum)


## sanity check station 
test_year = 20
test <- catchability[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()


# 4. set up fishing days values -----------------------------------------------

# we need to know how much fishing occurs in the region, so that the function knows how much to distribute.
# for commercial fishing, we are using a variety of sources to reconstruct the trends (a),
# the sources are split between the Metropolitan area effort (North, 4.N) and the Southwest area effort (Wadandi, 4.W)
# we will (b.) split the yearly effort into months, then (c.) split the monthly effort into access points (boat ramps).

## 4.a. enter the overall effort values (boat days) ---------------------------
# obtain boat days from forecasting (see google slides on reconstruction)
future_boat_days <- read.csv("data/input_data/YIJARUP - Fishing effort reconstruction - FUTURE_fishing_effort.csv", skip = 2) %>% 
  dplyr::select(YEAR, BOAT_REC)

plot(future_boat_days$YEAR, future_boat_days$BOAT_REC, 
     pch = 19, col = colour_palette[4], xlim = c(year_start, year_end), ylim = c(0, max(future_boat_days$BOAT_REC)), 
     xlab = "Year", ylab = "Boat Days", main = "Original + Fake Boat Days")

## 4.b. split yearly fishing effort by month ----------------------------------

# obtain the monthly distribution of fishing from the literature
seasonal_multipliers <- c( # see figure 21c in Ryan et al. 2022 - THIS IS FOR RECREATIONAL FISHING BUT CANT FIND COMMERCIAL EQUIVALENT
  "01" = 0.14, "02" = 0.081, "03" = 0.097, "04" = 0.081,
  "05" = 0.033, "06" = 0.033, "07" = 0.033, "08" = 0.033,
  "09" = 0.033, "10" = 0.065, "11" = 0.11, "12" = 0.26
)
seasonal_multipliers <- seasonal_multipliers / sum(seasonal_multipliers) # standardise so it adds up to 1
barplot(seasonal_multipliers, col = colour_palette[6], main = "distribution of yearly \nboat fishing effort by month in % \n(deduced from Ryan et al. 2022, fig. 21c)")

# add monthly distribution back to the timeseries
boat_effort <- expand.grid(
  year = year_start:year_end,
  month = sprintf("%02d", 1:12)
) %>%
  arrange(year, month) %>%
  mutate(
    annual_boat_days = rep(future_boat_days$BOAT_REC, each = 12),
    monthly_effort = annual_boat_days * seasonal_multipliers[month]
  ) %>%
  dplyr::select(year, month, monthly_effort)

# sanity check station
ggplot(boat_effort, aes(x = as.Date(paste(year, month, "01", sep = "-")), y = monthly_effort)) +
  geom_line(color = colour_palette[4]) +
  labs(title = paste0("Monthly Commercial Fishing Effort ", year_start, "-", year_end),
       x = "Date", y = "Monthly Boat Days") +
  theme_minimal() +
  geom_smooth(color = colour_palette[5])


## 4.c. distribute monthly effort into boat ramps -----------------------------

glimpse(BR)
plot(water$geometry); plot(BR$geometry, col = colour_palette[5], pch = 16, cex = 1, add = TRUE)

glimpse(boat_effort)
logistic_growth <- function(t, k = 0.3, t0 = 10) {
  1 / (1 + exp(-k * (t - t0)))
}
# we want a month x access point x year
distributed_effort <- array(0, 
                            dim = c(12, 
                                    nrow(BR),
                                    n_years_tot))

for (YEAR in 1:n_years_tot) {
  
  built_logic <- BR$build_year <= (year_start + YEAR - 1) # find which carparks are built in this year
  built_carparks <- which(built_logic) # and their index
  
  # for logistic growth of effort
  years_since_build <- ifelse((year_start + YEAR -1) - BR$build_year[built_carparks] < 0, 0, (year_start + YEAR -1) - BR$build_year[built_carparks])
  weights <- logistic_growth(years_since_build, k = 0.3, t0 = 10)
  
  for (MONTH in 1:12) {
    
    current_effort <- boat_effort %>% 
      dplyr::filter(year == (year_start + YEAR - 1), as.numeric(month) == MONTH) %>% 
      dplyr::pull(monthly_effort) # extract how many fishing days to distribute
    
    distributed_effort[MONTH, built_carparks, YEAR] <- (current_effort * (weights / sum(weights)))   # and distribute them
    
  }
}


plot(x = 1:n_years_tot, y = distributed_effort[1,1,])
lines(x = 1:n_years_tot, y = distributed_effort[1,15,])
lines(x = 1:n_years_tot, y = distributed_effort[1,10,], col = "red", pch = 3)
# ugly as heck but this way the effort is a bit smooth.


# okay combine into a month x access point x year cube

access_point_effort <- distributed_effort

# rename dimension names so they're consistent with the rest of the objects used in the function.
dimnames(access_point_effort)[[1]] <- c(1:12)
dimnames(access_point_effort)[[3]] <- c(1:dim(access_point_effort)[3])


# check that utility and access_point_effort have the same order of access point
dimnames(access_point_effort)[[2]] <- names(access_dist %>% dplyr::select(!ID))

# utility # watch out that the order of ramp effort and of utility match: otherwise you're assigning the effort to the wrong ramp in the function


# 5. create a list to use in the C++ function ---------------------------------

# little details to include in the final list
cell_area_m2 <- water$cell_area # just cell area
coef_values <- tibble(log_utility = -0.848,
                      expected_catch = 1.5,
                      expected_catch_sq = -1.171,
                      log_cell_area = 1.134
) # coefficients from Matt's 2022 paper. they weigh the relative importance of each for the distribution of effort. the paper was for rec fishing


fishing_info_list <- list(catchability, utility, access_point_effort, cell_area_m2, coef_values)
names(fishing_info_list) <- c("catchability", "utility", "fishing_days", "cell_area_m2", "coef_values")
saveRDS(fishing_info_list, paste0(folder_name, "/SB_", scenario_name, "_boat_rec_fishing_info.rds"))


### END ###
