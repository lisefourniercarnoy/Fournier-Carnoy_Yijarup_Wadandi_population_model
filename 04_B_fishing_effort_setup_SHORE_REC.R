# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Various literature figures.
# Task:    Set up fishing effort (understand the split of fishing effort within a year, over space, and across years)
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    August 2026

# -----------------------------------------------------------------------------

# intended outputs are:
# catchability coefficient (% of all individuals removed from a cell by hypothetical fishing sweeping the cell)
# cell attractivity (which cells are more likely to be fished)
# fishing effort values (how many boat days every)

# Notes: This script goes like this:
# 1. calculate travel cost
# 2. calculate cell utility
# 3. catchability

## ultimately to distribute fishing effort we need to have a formula of this shape: 
# (access point built y/n) * (1 / (travel cost + 1))

# so the script below is in order, with:


# -----------------------------------------------------------------------------

# -- intended outputs are:
# -- cell attractiveness for commercial fishing (which cells are more likely to be fished)
# -- catchability of each cell (which will determine the % of all individuals removed from each cell by 1 unit of effort)
# -- fishing effort values (how many boat days every ramp gets in every timestep)

rm(list = ls())

library(tidyverse) # data manipulation
library(sf) # shapefiles
library(terra) # for fast raster computations
library(RColorBrewer) # plotting colours
library(abind) # dealing with matrices
library(sfnetworks) # distance from cell to cell
library(exactextractr) # extracting raster values
library(patchwork) # for plot arrangements

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))


## 0. Files needed ------------------------------------------------------------

file_water          <- "data/output_data/03_B_water.rds"
file_network        <- "data/output_data/03_B_network_shapefile.shp"
file_wa             <- "data/output_data/01_B_land.shp"
file_bathy          <- "data/input_data/SW_crop_AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
file_carpark_w      <- "data/input_data/wadandi_carparks.shp"
file_carpark_n      <- "data/input_data/north_carparks.shp"

year_start <- 1900
year_end <- 2024
n_years_tot <- year_end - year_start + 1
current_fleet <- "shore_rec"

common_crs = 7850
bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -31), crs = common_crs) %>%
  st_as_sfc() %>%
  st_transform(crs = common_crs)

NCELL <- nrow(readRDS(file_water))


## 1. Elements useful to multiple elements of the script ----------------------

# -- several things can be calculated once and used for multiple output calculations.

### 1.1 Cell area -------------------------------------------------------------

water <- readRDS(file_water) |> st_transform(common_crs)
cell_area <- water$cell_area / 1000000 # just cell area in km2
cell_area <- as.numeric(water$cell_area)

### 1.2 Temporal fishability --------------------------------------------------

# -- this is a cell x month x year cube that tells you whether the cell is closed (0), partially open (>0, <1) or open (1).

water <- readRDS(file_water) |> st_transform(common_crs)
glimpse(water)

water[[current_fleet]] <- ifelse(water[[current_fleet]] == "T", TRUE, 
                                 ifelse(water[[current_fleet]] == "F", FALSE, water[[current_fleet]]))
NCELL <- nrow(water)

# identify the important cells
temporal_cells <- water$ID[water$TC_status == TRUE]

# setup the arrays to calculate things.
water_area <- array(0, dim = c(NCELL, 12, (year_end-year_start+1))) # cells x months x years
water_area[, 1:12, 1] <- water$cell_area # for all months of the first year, the catchable area is the cell's area.

# -- the loop goes as follows: for each time step, 
# -- check whether there's a spatial restriction that year, if so the fishable area is 0
# -- then, if the cell has a temporal closure, restrict fishable area as needed.
for (YEAR in 1:dim(water_area)[3]) {
  
  current_year <- year_start + YEAR - 1
  
  for (MONTH in 1:dim(water_area)[2]) {
    
    # spatial closures
    restriction_dates <- as.numeric(substr(water$SC_restriction_date, 7, 10))
    water_area[, MONTH, YEAR] <- water$cell_area
    
    restricted_cells <- which(
      !is.na(water[[current_fleet]]) & water[[current_fleet]] == FALSE & # where fleet is not allowed,
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


# -- sanity check station

test_cell <- 200
test_month <- 10
test_year <- 125

# -- reference cell numbers as of 13.01.2026:
# -- cockburn sound cell (temporal closure): 1
# -- SWC NTZ cell: 200
# -- random fished cell: 1000

# restrictions in this cell should be:
glimpse(st_drop_geometry(water[test_cell, c("TC_status", "SC_status", "commercial", "SC_restriction_date", "TC_restriction_date", "TC_restriction_months", "TC_restriction_perc_fished")]))

# check that it is correct
cat("In month", test_month, "of year", (year_start+test_year-1),
    ", the cell is", ifelse(water_area[test_cell, test_month, test_year]>0, "fishable", "NOT fishable"), "for commercial boats.", "Fishable area: ", water_area[test_cell, test_month, test_year]/1e06, "km2")


# CHECK: closure of cells
test_year = c(10, 50, 100, 120)
test_month = 12
p <- ggplot(data = water %>%
              mutate(as.data.frame(water_area[, test_month, test_year]) %>%
                       setNames(paste0("year_", test_year))) %>%
              tidyr::pivot_longer(
                cols = starts_with("year_"),
                names_to = "year",
                names_prefix = "year_",
                names_transform = as.numeric,
                values_to = "test"
              )
) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[c(3, 5)]) +
  facet_wrap(~ year, ncol = length(test_year)) +
  theme_minimal()
ggsave("plots/script_plot_checks/04_B/04_B_fishable_area_shore_rec.png", plot = p, width = 10, height = 6, dpi = 1000)


## 2. Fisher cell attractivity ------------------------------------------------

# -- ultimately to distribute fishing effort we need to have a formula of this shape: 
# -- (travel cost * coef1) + (log(cell area) * coef3)
# -- ‾‾‾‾‾‾‾‾‾‾‾‾                 ‾‾‾‾‾‾‾‾‾
# --     2.2                        1.2

# -- the coefficients we do not have from published literature, and will therefore be eyeballed.
# -- in commercial fishing, we consider all ramps to be built from the start. boat rec fishing would include build date in this attractivity.


### 2.1 Ramp built ------------------------------------------------------------

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

# -- there are so many access points, some of which are quite recent. 
# -- we'll zero the attractivity of access points if they haven't been built yet.

build_year <- array(0,dim = c(nrow(AP),n_years_tot))

for (ACCESS in 1:nrow(AP)) {
  for (YEAR in 1:n_years_tot) {
    yr <- year_start + YEAR - 1
    build_year[ACCESS, YEAR] <- ifelse(AP$year_start[ACCESS] <= yr, 1, NA)
  }
}


### 2.2 Travel cost -----------------------------------------------------------

# -- unlike boat fishers, shore fishers are limited to the coast. 

#### 2.3.1 Distance from access points ----------------------------------------

# -- the intended output here is an object of size cell x access point

# obtain water and access points
water <- readRDS(file_water) |> st_transform(common_crs)

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

AP <- st_as_sf(AP) %>% st_transform(common_crs)

network <- st_read(file_network); plot(network$geometry)

# find the centre of each grid cell
sf::sf_use_s2(FALSE)
centroids <- st_centroid(st_make_valid(water))
sf::sf_use_s2(TRUE)
points <- as.data.frame(st_coordinates(centroids))%>% # the points start at the bottom left and then work their way their way right
  mutate(ID = row_number()) 
points_sf <- st_as_sf(points, coords = c("X", "Y"), crs = common_crs)


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

access_dist <- units::drop_units(access_dist) %>% mutate_at(vars(1:(ncol(access_dist)-1)), ~replace(., .>2, NA)) # cells over 3 km from the access point aren't accessible to shore fishers (extreme negative value to avoid NAs in the C++ function)
access_dist[water$ID[!water$type %in% c("shore_wadandi", "shore_north")], 1:(ncol(access_dist)-1)] <- NA # offshore cells are not accessible to shore fishers


# CHECK: access point distance
test_access <- c(100, 10, 300, 1)
access_dist2 <- access_dist %>% mutate(ID = centroids$ID)
water_dist_long <- water %>%
  left_join(access_dist2, by = "ID") %>%
  pivot_longer(cols = names(access_dist %>% dplyr::select(-ID)), names_to = "AP", values_to = "dist_km")
ap_names <- unique(water_dist_long$AP)[test_access]

plots <- lapply(ap_names, function(ap_name) {
  ap_point <- AP %>% dplyr::filter(name == ap_name)
  bb <- st_bbox(st_buffer(ap_point, dist = 4000))  # 2km buffer in metres, adjust as needed
  
  ggplot() +
    geom_sf(data = st_read(file_wa),  col = NA, fill = "gray") +
    geom_sf(data = water_dist_long %>% dplyr::filter(AP == ap_name),
            aes(fill = as.numeric(dist_km)), color = "gray") +
    geom_sf(data = ap_point, shape = 21, fill = "red") +
    scale_fill_gradientn(colours = colour_palette[4:6], na.value = "transparent") +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"])) +
    labs(title = ap_name) +
    theme_minimal()
})
p <- wrap_plots(plots, ncol = 2)
ggsave("plots/script_plot_checks/04_B/04_B_access_point_distance_shore_rec.png", plot = p, width = 10, height = 7, dpi = 500)


#### 2.3.2 Fuel prices & efficiency -------------------------------------------

# -- we want to obtain one value per year, expressed in $/km

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


### 2.3 Calculate cell attractivity -------------------------------------------

# -- okay so now that we have most of our elements, we can make a function that defines the 'attractivity' of each cell
# -- the format should be: cell x access point x year

attractivity <- array(dim = c(NCELL, nrow(AP), 12, n_years_tot))

for(YEAR in 1:n_years_tot) {
  for(MONTH in 1:12) {
    for (ACCESS in 1:nrow(AP)) {
      attractivity[, ACCESS, MONTH,  YEAR] <- 
        -0.204 * fuel_price_clean$price_per_km[YEAR] * (2 * as.numeric(access_dist[, ACCESS])) + # travel cost, with round-trip distance
        1.134 * log(cell_area + 1) # cell area
      
      # also turn off cells that are under closure
      attractivity[, ACCESS, MONTH, YEAR][water_area[, MONTH, YEAR] == 0] <- NA
    }
  }
}
tail(attractivity[,1, 12, 125])

# CHECK: attractivity of cells from the POV of few access points
test_year = 100
test_month = 12
test_access = c(1, 10, 100, 200)
ap_names <- unique(water_dist_long$AP)[test_access]
plots <- lapply(ap_names, function(ap_name) {
  ap_point <- AP %>% dplyr::filter(name == ap_name)
  bb <- st_bbox(st_buffer(ap_point, dist = 4000))  # 2km buffer in metres, adjust as needed
  
  ggplot() +
    geom_sf(data = st_read(file_wa),  col = NA, fill = "gray") +
    geom_sf(data = purrr::map_dfr(test_access, function(a) {
      water %>%
        mutate(attractivity = attractivity[, a, test_month, test_year],
               access_id = a)
    }), aes(fill = attractivity), colour = "gray") +
    geom_sf(data = ap_point, shape = 21, fill = "red") +
    scale_fill_gradientn(colours = colour_palette[4:6], na.value = "transparent") +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"])) +
    labs(title = ap_name) +
    theme_minimal()
})
p <- wrap_plots(plots, ncol = 2)
ggsave("plots/script_plot_checks/04_B/04_B_cell_attractivity_shore_rec.png", plot = p, width = 10, height = 10, dpi = 500)


# transform in a format that C++ can deal with easier
attractivity <- lapply(1:12, function(m) {
  arr <- attractivity[, , m, , drop = FALSE]   # keep 4D, dim3 = 1
  dim(arr) <- dim(attractivity)[-3]             # collapse only the month axis -> (NCELL, n_access, n_years)
  arr
})


## 3. Shore rec catchability --------------------------------------------------

# -- catchability is the susceptibility of fish to be caught by 1 unit of effort. 
# -- catchability differs for every cell, given that fish is more susceptible to being caught by 1 unit of effort in a cell of 10m2 compared to 1 unit of effort in a cell of 100km2.

# -- catchability works with fishing effort to produce fishing mortality. 
# -- it is unknowable, so we will need to calibrate it. 
# -- see script 04_D for process 

q <- readRDS("data/output_data/04_D_shore_rec_q.rds")$Q

shore_mask <- water$type %in% c("shore_wadandi", "shore_north")  # TRUE for shore cells only
fishable_area <- water_area
fishable_area[!shore_mask, , ] <- 0   # zero out offshore cells so they don't count toward the total

fishable_area_sum <- colSums(fishable_area) # month x year
fishable_area_perc <- sweep(fishable_area, c(2, 3), fishable_area_sum, FUN = "/")
catchability <- array(0, dim = dim(fishable_area)) # cell x month x year

for (y in 1:dim(fishable_area)[3]) {
  for (m in 1:12) {
    catchability[, m, y] <- q[y] / fishable_area_perc[, m, y] # divide here - 1 unit of effort in a small area gets a lot of fish, but in a large area there's more 'places for fish not to be caught'
  }
}
catchability[catchability == Inf] <- 0  # if cells are not fishable, division by zero (in loop above). change to zero to avoid messing everything up later on.

# CHECK: catchability 
test_year = c(90, 100, 110, 120)
test_access <- c(100, 10, 300, 1)
catch_mat <- catchability[,1, test_year]
ap_names <- unique(water_dist_long$AP)[test_access]

plots <- lapply(ap_names, function(ap_name) {
  ap_point <- AP %>% dplyr::filter(name == ap_name)
  bb <- st_bbox(st_buffer(ap_point, dist = 4000))
  
  plot_data <- water %>%
    bind_cols(as.data.frame(catch_mat) %>%
                setNames(paste0("year_", test_year))) %>%
    tidyr::pivot_longer(
      cols = starts_with("year_"),
      names_to = "year",
      names_prefix = "year_",
      names_transform = as.numeric,
      values_to = "catchability"
    ) %>%
    dplyr::mutate(catchability = dplyr::na_if(catchability, 0))  # zero -> NA
  
  ggplot(data = plot_data) +
    geom_sf(aes(fill = catchability), colour = NA) +
    geom_sf(data = AP[test_access, ] %>% mutate(access_id = test_access), shape = 21, fill = "red") +
    scale_fill_gradientn(colours = colour_palette[c(3, 5)], trans = "sqrt", na.value = "transparent") +
    facet_wrap(~ year) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"])) +
    labs(title = ap_name) +
    theme_minimal()
})
p <- wrap_plots(plots, ncol = 2)
ggsave("plots/script_plot_checks/04_B/04_B_access_point_catchability_shore_rec.png", plot = p, width = 10, height = 7, dpi = 500)


## 4. Fishing effort values ---------------------------------------------------

# -- we know how much fishing there is (roughly) but we now need to distribute it across access points.

### 4.1 Enter the overall effort values ---------------------------------------

g_sheets <- read.csv("data/input_data/YIJARUP - Fishing effort reconstruction - FINAL_fishing_effort (14-07-2026).csv", skip = 3) %>% 
  dplyr::select(c(YEAR, state.pop:last_col())) %>% #select only shore fishing columns.
  glimpse()

plot(g_sheets$YEAR, g_sheets$state.pop, type = "l", col = colour_palette[5], lwd = 4,
     xlab = "Year", ylab = "WA population",
     ylim = c(0, max(g_sheets$state.pop))
)

lines(g_sheets$YEAR, g_sheets$shore.fishing.days, col = colour_palette[4], lwd = 4)

legend("topleft", legend = c("shore fishing days", "WA population"), col = c(colour_palette[4], colour_palette[5]), lwd = 4)

annual_effort_shore <- g_sheets %>% dplyr::select(c(YEAR, shore.fishing.days)) %>% 
  rename(year = YEAR,
         shore_days = shore.fishing.days)

# CHECK: shore rec fishing effort
plot(annual_effort_shore$year, annual_effort_shore$shore_days, type = "l", col = colour_palette[5], lwd = 4,
     main = "Shore Days in Wadandi Country \nfrom the Literature \n(with some extrapolation)", 
     xlab = "Year", ylab = "Boat Days")
legend("bottomright", legend = c("shore days (derived from population)"), col = c(colour_palette[5]), lwd = 4)


### 4.2 Split yearly fishing effort by month ----------------------------------

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


### 4.3. Distribute monthly effort into access points -------------------------

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

## 5. Create a list to use in the C++ function --------------------------------

# little details to include in the final list
coef_values <- tibble(expected_catch = 1,
                      expected_catch_sq = 1
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

saveRDS(fishing_info_list, "data/output_data/04_B_shore_rec_fishing_info.rds")


## 6. Setting up effort for the burn-in ---------------------------------------

# -- the burn-in exists only to stabilise the population before the simulations start.
# -- to stabilise the population, we need to run the function with a small amount of fishing (smaller than what we start with)

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
# for the burn in, just replace the corresponding list object from section 6 with this one.


### END ###
