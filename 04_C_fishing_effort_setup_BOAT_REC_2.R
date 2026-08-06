# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Data from project files.
# Task:    Set up info cubes for commercial fishing, to use in C++ function
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    March 2026

# -----------------------------------------------------------------------------

# intended outputs are:
# catchability coefficient (% of all individuals removed from a cell by hypothetical fishing sweeping the cell)
# cell attractiveness (which cells are more likely to be fished)
# fishing effort values (how many boat days every ramp gets in every timestep)

## ultimately to distribute fishing effort we need to have a formula of this shape: 
# (travel cost * coef1) + (log(offshore dist) * coef2) + (log(cell area) * coef3) + (expected catch * coef5) + (expected catch^2 * coef5)
# the latter two will be calculated in the C++ function but the rest we can do here.

# so the script below is in order, with:

## Load libraries -------------------------------------------------------------

rm(list = ls())

library(tidyverse) # data manipulation
library(sf) # shapefiles
library(terra) # for fast raster computations
library(RColorBrewer) # plotting colours
library(abind) # dealing with matrices
library(sfnetworks) # distance from cell to cell
library(exactextractr) # extracting raster values


## 0. files needed ----

file_water <- "data/output_data/03_B_water.rds"
file_network <- "data/output_data/03_B_network_shapefile.shp"
file_wa <- "data/output_data/01_B_land.shp"
file_bathy<- "data/input_data/SW_crop_AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
file_boat_ramps_w <- "data/input_data/access_point_shapefiles/wadandi_boat_ramps.shp"
file_boat_ramps_n <- "data/input_data/access_point_shapefiles/north_boat_ramps.shp"

year_start <- 1900
year_end <- 2024
n_years_tot <- year_end - year_start +1
current_fleet <- "boat_rec"

common_crs <- 4283
crs_raster <- "+proj=longlat +datum=WGS84 +no_defs"
bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -31), crs = common_crs) %>%
  st_as_sfc() %>%
  st_transform(crs_raster)

NCELL <- nrow(readRDS(file_water))

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))


## 1. fuel prices ----

# we want to obtain one value per year, expressed in $/km

# obtain fuel prices from the literature
fuel_price <- read.csv("data/input_data/fuel_prices.csv") # data from https://www.bitre.gov.au/sites/default/files/is_082.pdf

# adjust for inflation
base_cpi <- fuel_price$CPI[fuel_price$year == 2012]
fuel_price$price_inflation <- fuel_price$petrol_price * (base_cpi / fuel_price$CPI) # adjust the fuel prices for inflation
fuel_price$price_dollar <- fuel_price$price_inflation / 100 # convert into $

plot(fuel_price$price_dollar ~ fuel_price$year, type = "l") # plot check

fuel_price$price_per_km <- fuel_price$price_dollar / (1 / 2.54) # here we do cross-product of $, L and km. 2.54 is a value from Nicole Hamre

# fill in missing values
fuel_price_clean <- data.frame(year = year_start:year_end) %>%
  left_join(fuel_price %>% dplyr::select(year, price_per_km), by = "year") %>%
  fill(price_per_km, .direction = "downup") %>% 
  glimpse()
plot(fuel_price_clean$price_per_km ~ fuel_price_clean$year, type = "l") # plot check


## 2. distance from access point ----

water <- readRDS(file_water)

BR_n <- st_read(file_boat_ramps_n) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  glimpse()

BR_w <- st_read(file_boat_ramps_w) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  glimpse()

# for commercial fishing, only a few boat ramps are used, so select only correct ramps (see Andrea Gaynor's historical fishing resources, fishing localities in ABS stats)
unique(BR_w$name)
plot(water$geometry); plot(BR_w$geometry, col = colour_palette[5], pch = 16, cex = 1, add = TRUE)

glimpse(BR_n)
glimpse(BR_w)

BR <- rbind(BR_n %>% dplyr::select(name, geometry) %>% mutate(region = "north"), 
            BR_w %>% dplyr::select(name, geometry) %>% mutate(region = "wadandi")) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>% 
  glimpse()
BR <- BR[!sf::st_is_empty(BR), ]

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
access_dist <- access_dist / 1000 # in km
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


## 3. shore distance ----

water <- readRDS(file_water)

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


## 4. cell area ----

water <- readRDS(file_water)
cell_area <- as.numeric(water$cell_area) / 1000000


## 5. fishable depth ----

# this section of the fishing effort relates to boats being able to fish further over the years, with bigger boats, more powerful engines etc.
# the model takes this into account by adding a 'fishable depth' into the mix, which goes more and more offshore.
# for rec fishing, we're assuming that on average, fishers can fish 1m more every year. (change as needed in parameters)
NCELL <- nrow(water)

# prepare grid cells
plot(water$geometry)
wa_mask <- st_read(file_wa); wa_mask <- st_transform(wa_mask, crs = crs_raster); wa_mask <- as(wa_mask, "Spatial"); plot(wa_mask, col = "lightgray", add = T)
plot(bbox, add = T)

# set parameters
years <- year_start:year_end
min_depth <- 20 # fishable depth is about 20m in 1945 (from Gaynor 2008, p.38)
increase_year <- 1950 # the year that boats start to go deeper (post-war industrialisation)
depth_per_year <- 1 # fishable depth gained per year (m) - arbitrary to reach the continental shelf quickly.

# prepare the bathymetry layer
bathy <- rast(file_bathy) %>%
  project(crs_raster) %>% 
  crop(as(water, "Spatial")) %>%
  abs() %>%
  terra::mask(vect(wa_mask), inverse = TRUE); plot(bathy)

thresholds <- min_depth + pmax(0, depth_per_year*(years - increase_year)) # depth limits for each year

fishable_stack <- rast(
  lapply(thresholds, function(th) bathy <= th)
)

names(fishable_stack) <- paste0("year_", years)

fishable_summary <- terra::extract( # this extracts the fishability (TRUE or FALSE) at cell centroids
  fishable_stack,
  centroids,
  na.rm = TRUE
)

fishable_summary$ID <- water$ID

fishable_depth_cell_year <- fishable_summary %>%
  arrange(ID) %>%
  dplyr::select(-ID) %>%
  as.matrix()

# strip "year_" prefix from colnames
colnames(fishable_depth_cell_year) <- gsub("year_", "", colnames(fishable_depth_cell_year))

fishable_depth_cell_year[is.na(fishable_depth_cell_year)] <- TRUE # replace NA with 1, these are cells close to shore which are out of the bathymetry layer (but would be shallow enough to fish from the start)

## sanity check station 
test_year = 70
ggplot(data = water %>% mutate(test = fishable_depth_cell_year[, test_year])) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()


## XX. temporal fishability ---------------------------------------------------

# this is just a cell x month x year cube that tells you whether there's a closure or not in the cell.
# some cells are closed for only part of the month (ex. temporal closure that starts midway through the month) 
# so that cell would have the value 0.5 that month.

glimpse(water)

## in this section (1.) we calculate the fishable area of the cells in each time step (1.a), then divide it by the fishable area of the whole grid at each time step (1.b)
## which gives us each the portion of Catchability of each cell.
## because the area that is catchable changes (with spatial and temporal restrictions), we have to calculate the grid's catchability for each month of each year.

water <- readRDS(file_water) %>% filter(!is.na(ID))
water[[current_fleet]] <- ifelse(water[[current_fleet]] == "T", TRUE, 
                                 ifelse(water[[current_fleet]] == "F", FALSE, water[[current_fleet]]))
NCELL <- nrow(water)

# identify the important cells
temporal_cells <- water$ID[water$TC_status == TRUE]

# setup the arrays to calculate things.
water_area <- array(0, dim = c(NCELL, 12, (year_end-year_start+1))) # cells x months x years
water_area[, 1:12, 1] <- water$cell_area # for all months of the first year, the catchable area is the cell's area.

# the loop goes as follows: for each time step, 
# check whether there's a spatial restriction that year, if so the fishable area is 0
# then, if the cell has a temporal closure, restrict fishable area as needed.
for (YEAR in 1:dim(water_area)[3]) {
  
  current_year <- year_start + YEAR - 1
  
  for (MONTH in 1:dim(water_area)[2]) {
    
    # spatial closures
    restriction_dates <- as.numeric(substr(water$SC_restriction_date, 7, 10))
    water_area[, MONTH, YEAR] <- water$cell_area
    
    restricted_cells <- which(
      !is.na(water$boat_rec) & water[[current_fleet]] == FALSE & # where fleet is not allowed,
        !is.na(restriction_dates) & current_year >= restriction_dates # and when SC is in place...
    )
    water_area[restricted_cells, MONTH, YEAR] <- 0 # ...the cell is not fishable
    
    # temporal closures
    for (CELL in temporal_cells) {
      
      # find the relevant temporal restrictions
      TC_years_list <- water$TC_restriction_date[[CELL]] # years the closure has changed
      TC_months_list <- water$TC_restriction_months[[CELL]] # months restricted
      TC_perc_list <- water$TC_restriction_perc_fished[[CELL]] # how much they are restricted
      
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
} # this loop calculates for every month and every year, the fishable area of each cell.


### sanity check station

test_cell <- 225
test_month <- 10
test_year <- 125

# reference cell numbers as of 13.01.2026:
## cockburn sound cell (temporal closure): 1
## SWC NTZ cell: 200
## random fished cell: 1000

# restrictions in this cell should be:
glimpse(st_drop_geometry(water[test_cell, c("TC_status", "SC_status", "boat_rec", "SC_restriction_date", "TC_restriction_date", "TC_restriction_months", "TC_restriction_perc_fished")]))

# check that it is correct
cat("In month", test_month, "of year", (year_start+test_year-1),
    ", the cell is", ifelse(water_area[test_cell, test_month, test_year]>0, "fishable", "NOT fishable"), "for rec boats.", "Fishable area: ", water_area[test_cell, test_month, test_year]/1e06, "km2")

## sanity check station 
test_year = 125
test <- water_area[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  geom_sf(aes(fill = ID[ID == 200]), colour = "red") +
  
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()

ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  geom_sf(data = water %>% filter(ID == 225), fill = "red") +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()

## 6. put it all together ----

attractivity <- array(dim = c(NCELL, nrow(BR), 12, n_years_tot))

for(YEAR in 1:n_years_tot) {
  for(MONTH in 1:12) {
    for (ACCESS in 1:nrow(BR)) {
      attractivity[, ACCESS, MONTH,  YEAR] <- 
        -0.204 * fuel_price_clean$price_per_km[YEAR] * (2 * as.numeric(access_dist[, ACCESS])) + # travel cost, with round-trip distance
        -0.848 * log(shore_dist$shore_dist + 1) + # offshore dist
        1.134 * log(cell_area + 1) # cell area
      
      # and then we add the depth fishability, to turn off cells that are too deep
      attractivity[, ACCESS, MONTH, YEAR][fishable_depth_cell_year[,YEAR] == 0] <- NA # effort gets allocated to cells with high attractivity. here all cells have a negative attractivity, so we set unfishable cells to be even more negative than that
      # also turn off cells that are under closure
      attractivity[, ACCESS, MONTH, YEAR][water_area[, MONTH,YEAR] == 0] <- NA
    }
  }
}
tail(attractivity[,1,12,125])

## sanity check:
test_year = 100
test_month = 12
test_access = 1
ggplot(data = water %>% mutate(test = attractivity[, test_access, test_month, test_year])) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()

# transform to a format that C++ deals with better
attractivity <- lapply(1:12, function(m) {
  arr <- attractivity[, , m, , drop = FALSE]   # keep 4D, dim3 = 1
  dim(arr) <- dim(attractivity)[-3]             # collapse only the month axis -> (NCELL, n_access, n_years)
  arr
})

## 7. calculate catchability --------------------------------------------------

# catchability is the % of the population being harvested by 1 unit effort.
# it is unknowable, so we will need to calibrate it. 
# see script 04_D for process 

q <- readRDS("data/output_data/04_D_boat_rec_q.rds")$Q

fishable_depth_expanded <- array(NA, dim = c(NCELL, 12, n_years_tot)); for (m in 1:12) {fishable_depth_expanded[, m, ] <- fishable_depth_cell_year}
fishable_area <- water_area * fishable_depth_expanded # cell x month x year
fishable_area_sum <- colSums(fishable_area) # month x year
fishable_area_perc <- sweep(fishable_area, c(2, 3), fishable_area_sum, FUN = "/")
catchability <- array(0, dim = dim(fishable_area))
for (y in 1:dim(fishable_area)[3]) {
  for (m in 1:12) {
    catchability[, m, y] <- q[y] / fishable_area_perc[, m, y] # divide here - 1 unit of effort in a small area gets a lot of fish, but in a large area there's more 'places for fish not to be caught'
  }
}
catchability[catchability == Inf] <- 0  # if cells are not fishable, division by zero (in loop above). change to zero to avoid messing everything up later on.

## sanity check station 
test_year = 50
test <- catchability[, 1, test_year]
ggplot(data = water %>% mutate(test = (test))) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()
test_year = 125
test <- catchability[, 10, test_year]
test[test > 0.1] <- NA
ggplot(data = water %>% mutate(catchability = test)) +
  geom_sf(aes(fill = catchability), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()



# 4. set up fishing days values -----------------------------------------------

# we need to know how much fishing occurs in the region, so that the function knows how much to distribute.
# for commercial fishing, we are using a variety of sources to reconstruct the trends (a),
# the sources are split between the Metropolitan area effort (North, 4.N) and the Southwest area effort (Wadandi, 4.W)
# we will (b.) split the yearly effort into months, then (c.) split the monthly effort into access points (boat ramps).


## 4.a. enter the overall effort values (boat days) ---------------------------

g_sheets <- read.csv("data/input_data/YIJARUP - Fishing effort reconstruction - FINAL_fishing_effort (14-07-2026).csv", skip = 3) %>% 
  dplyr::select(c(YEAR, boat.days.west.coast)) %>% # select only boat rec fishing columns.
  glimpse()

plot(g_sheets$YEAR, g_sheets$boat.days.west.coast)

annual_effort_boat <- g_sheets %>% dplyr::select(c(YEAR, boat.days.west.coast)) %>% 
  rename(year = YEAR,
         boat_days = boat.days.west.coast)
annual_effort_boat$boat_days <- as.vector(approx(annual_effort_boat$boat_days, n = nrow(annual_effort_boat))$y)

# check
plot(annual_effort_boat$year, annual_effort_boat$boat_days, type = "l", col = colour_palette[5], lwd = 4,
     xlab = "Year", ylab = "Boat Days")


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
    annual_boat_days = rep(annual_effort_boat$boat_days, each = 12),
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


BR <- st_as_sf(BR); BR <- BR %>% st_transform(common_crs) %>% dplyr::filter(!st_is_empty(.))
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


plot(x = 1:125, y = distributed_effort[1,1,])
lines(x = 1:125, y = distributed_effort[1,15,])
lines(x = 1:125, y = distributed_effort[1,10,], col = "red", pch = 3)
# ugly as heck but this way the effort is a bit smooth.


# okay combine into a month x access point x year cube

access_point_effort <- distributed_effort

# rename dimension names so they're consistent with the rest of the objects used in the function.
dimnames(access_point_effort)[[1]] <- c(1:12)
dimnames(access_point_effort)[[3]] <- c(1:dim(access_point_effort)[3])


# check that utility and access_point_effort have the same order of access point
dimnames(access_point_effort)[[2]] <- names(access_dist %>% dplyr::select(!ID))

# watch out that the order of ramp effort and of utility match: otherwise you're assigning the effort to the wrong ramp in the function

# 5. create a list to use in the C++ function ---------------------------------

# little details to include in the final list
coef_values <- tibble(expected_catch = 1.5,
                      expected_catch_sq = -1.171
) # coefficients from Matt's 2022 paper. they weigh the relative importance of each for the distribution of effort. the paper was for rec fishing


fishing_info_list <- list(catchability, 
                          attractivity, 
                          access_point_effort, 
                          cell_area, 
                          coef_values)
names(fishing_info_list) <- c("catchability", 
                              "attractivity", 
                              "fishing_days", 
                              "cell_area", 
                              "coef_values")
saveRDS(fishing_info_list, "data/output_data/04_C_boat_rec_fishing_info.rds")


## 6. SETTING UP EFFORT FOR BURN IN -------------------------------------------

# the burn-in exists only to stabilise the population before the simulations start.
# to stabilise the population, we need to run the function with a small amount of fishing (smaller than what we start with)

# obtain the relative contribution of each boat ramp
BR_1900 <- BR %>%
  dplyr::filter(build_year == 1900) %>% 
  dplyr::select(name) %>% 
  glimpse()

# calculate the low level of fishing
eq.init.fish = 0.025 # this is the new level of fishing mortality, very low.
q = 0.00001 # catchability, or the risk of a fish being caught
effort = (-log(1-eq.init.fish))/q # this is the number of this fleet's fishing days in a year.


# we then split this effort in every month
seasonal_multipliers
burn_in_effort <- seasonal_multipliers*effort

burn_in_effort <- as.data.frame(burn_in_effort) %>%
  mutate(!!!setNames(rep(list(0), length(BR$name)), BR$name)) %>% 
  rename(effort = "burn_in_effort")


## Split up by boat ramp
for(ap in BR_1900$name){
  burn_in_effort[ap] = burn_in_effort$effort / length(BR_1900$name)
}

burn_in <- as.matrix(burn_in_effort[-1])

saveRDS(burn_in, file = "data/output_data/04_C_boat_rec_burn_in_effort.rds")
# for the burn in, just replace the corresponding list object from section 5 with this one.

### END ###

