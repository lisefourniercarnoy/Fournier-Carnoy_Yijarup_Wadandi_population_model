# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Various literature figures.
# Task:    Set up fishing effort (understand the split of fishing effort within a year, over space, and across years)
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    March 2026

# -----------------------------------------------------------------------------

# intended outputs are:
# catchability coefficient (% of all individuals removed from a cell by hypothetical fishing sweeping the cell)
# cell utility (which cells are more likely to be fished)
# fishing effort values (how many boat days every)

# Notes: This script goes like this:
# 1. calculate travel cost
# 2. calculate cell utility
# 3. catchability

# -----------------------------------------------------------------------------

## Load libraries -------------------------------------------------------------

rm(list = ls())

library(tidyverse) # data manipulation
library(sf) # shapefiles
library(terra) # for fast raster computations
library(RColorBrewer) # plotting colours
library(abind) # dealing with matrices
library(sfnetworks) # distance from cell to cell
library(exactextractr) # extracting raster values


## Custom plotting parameters -------------------------------------------------
colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))
source("custom_theme.R")

## 0. Files used in this script -----------------------------------------------

file_wa           <- "data/output_data/01_B_land.shp"
file_bathy        <- "data/input_data/SW_crop_AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
file_carpark_w    <- "data/input_data/wadandi_carparks.shp"
file_carpark_n    <- "data/input_data/north_carparks.shp"
file_water        <- "data/output_data/03_B_water.rds"
file_network      <- "data/output_data/03_B_network_shapefile.shp"

year_start <- 1900
year_end <- 2024
n_years_tot <- year_end - year_start +1
current_fleet <- "shore_rec"

common_crs <- 4283
crs_raster <- "+proj=longlat +datum=WGS84 +no_defs"
bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -31), crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(crs_raster)

water <- readRDS(file_water) %>% filter(!is.na(ID))
NCELL <- nrow(water)

fishing_info_list <- list()

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

ggsave("plots/checking_plots_during_setup/04_B_shore_rec_carpark_distance.png", plot = last_plot(), height = 10, width = 10)



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
utility[,1:20,100]
utility[!is.finite(utility)] <- 0
utility[is.na(utility)] <- 0



## sanity check - is utility higher for areas close to access points but far from shore)
test_access = 10
test_year = 71
test <- water %>% dplyr::select(!where(is.list)) %>% mutate(test = utility[, test_access, test_year])
mapview::mapview(test, zcol = "test") + mapview::mapview( AP[test_access,])

ggplot() +
  geom_sf(data = water %>% dplyr::select(!where(is.list)) %>% mutate(test = utility[, test_access, 125]), 
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
      !is.na(water$shore_rec) & water[[current_fleet]] == FALSE & # where fleet is not allowed,
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

water_area[water$ID[!water$type %in% c("shore_wadandi", "shore_north")],,] <- 0 # offshore cells are not open to shore-fishers

### sanity check station

test_cell <- 200
test_month <- 10
test_year <- 125

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
test_year = 124
test <- water_area[, 1, test_year]
mapview::mapview(water %>% dplyr::select(!where(is.list)) %>% mutate(test = test), zcol = "test")


ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  labs(y = "Fishable proportion", colour = "Cell ID") +
  theme_minimal()

### 3.b.b calculate catchability ----------------------------------------------

# catchability is the % of the population being harvested by 1 unit effort.
# it is unknowable, so we will need to calibrate it. 
# see script 04_D for process 

q <- readRDS("data/output_data/04_D_shore_rec_q.rds")$Q

fishable_area <- water_area
fishable_area_sum <- colSums(fishable_area)
fishable_area_perc <- sweep(fishable_area, c(2, 3), fishable_area_sum, FUN = "/")

catchability <- array(0, dim = dim(fishable_area))
for (y in 1:dim(fishable_area)[3]) {
  for (m in 1:12) {
    catchability[, m, y] <- q[y] / fishable_area_perc[, m, y] # divide here - 1 unit of effort in a small area gets a lot of fish, but in a large area there's more 'places for fish not to be caught'
  }
}
catchability[catchability == Inf] <- 0  # avoid division by zero

## sanity check station 
test_year = 100
test <- catchability[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()


# 4. set up fishing days values -----------------------------------------------

# we need to know how much fishing occurs in the region, sso that the function knows how much to distribute.
# for commercial fishing, we are using a variety of sources to reconstruct the trends (a),
# the sources are split between the Metropolitan area effort (North, 4.N) and the Southwest area effort (Wadandi, 4.W)
# we will (b.) split the yearly effort into months, then (c.) split the monthly effort into access points (boat ramps).



## 4.a. enter the overall effort values (boat days) ---------------------------

g_sheets <- read.csv("data/input_data/YIJARUP - Fishing effort reconstruction - FINAL_fishing_effort (14-07-2026).csv", skip = 3) %>% 
  dplyr::select(c(YEAR, state.pop:last_col())) %>% #select only shore fishing columns.
  glimpse()

plot(g_sheets$YEAR, g_sheets$state.pop, type = "l", col = colour_palette[5], lwd = 4,
     #main = "Shore Boat Days in \nWadandi Country (dashed) and Metro (solid) \nfrom the Literature (with some extrapolation)", 
     xlab = "Year", ylab = "WA population",
     ylim = c(0, max(g_sheets$state.pop))
)

lines(g_sheets$YEAR, g_sheets$shore.fishing.days, col = colour_palette[4], lwd = 4)

legend("topleft", legend = c("shore fishing days", "WA population"), col = c(colour_palette[4], colour_palette[5]), lwd = 4)

annual_effort_shore <- g_sheets %>% dplyr::select(c(YEAR, shore.fishing.days)) %>% 
  rename(year = YEAR,
         shore_days = shore.fishing.days)

# check
plot(annual_effort_shore$year, annual_effort_shore$shore_days, type = "l", col = colour_palette[5], lwd = 4,
     main = "Shore Days in Wadandi Country \nfrom the Literature \n(with some extrapolation)", 
     xlab = "Year", ylab = "Boat Days")
legend("bottomright", legend = c("shore days (derived from population)"), col = c(colour_palette[5]), lwd = 4)


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
    annual_shore_days = rep(annual_effort_shore$shore_days, each = 12),
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


## 4.c. distribute monthly effort into access points --------------------------

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


plot(x = 1:125, y = distributed_effort[1,1,])
lines(x = 1:125, y = distributed_effort[1,50,])
points(x = 1:125, y = distributed_effort[1,300,], col = "red", pch = 3)



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
saveRDS(fishing_info_list, "data/output_data/04_B_shore_rec_fishing_info.rds")


## 6. SETTING UP EFFORT FOR BURN IN -------------------------------------------

# the burn-in exists only to stabilise the population before the simulations start.
# to stabilise the population, we need to run the function with a small amount of fishing (smaller than what we start with)

# obtain the access points built in 1900
AP_n <- st_read(file_carpark_n) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  mutate(build_mnth = 1, # assume Jan if unknown
         id = 1:nrow(.),
         name = paste0(rough_area, "_", id)
         )

AP_w <- st_read(file_carpark_w) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  mutate(build_mnth = 1, # assume Jan if unknown
         #id = 1:nrow(.),
         name = area
  )

AP_all <- rbind(AP_n %>% dplyr::select(name, year_start, build_mnth), AP_w %>% dplyr::select(name, year_start, build_mnth))

AP_1900 <- AP_all %>% 
  dplyr::filter(year_start == 1900,
                build_mnth == 1)

# calculate the low level of fishing
eq.init.fish = 0.025 # this is the new level of fishing mortality, very low.
q = 0.0000001 # catchability, or the risk of a fish being caught
effort = (-log(1-eq.init.fish))/q # this is the number of this fleet's fishing days in a year.


# we then split this effort in every month
seasonal_multipliers
burn_in_effort <- seasonal_multipliers*effort

burn_in_effort <- as.data.frame(burn_in_effort) %>%
  mutate(!!!setNames(rep(list(0), nrow(AP_all)), colnames(access_point_effort))) %>% 
  rename(effort = "burn_in_effort")


## split up by access point (the ones that are built only)
for(ap in AP_1900$name) {
  burn_in_effort[ap] = burn_in_effort$effort * (1/nrow(access_point_effort))
}

burn_in <- as.matrix(burn_in_effort[-1])

saveRDS(burn_in, file = "data/output_data/04_B_shore_rec_burn_in_effort.rds")
# for the burn in, just replace the corresponding list object from section 5 with this one.


### END ###
