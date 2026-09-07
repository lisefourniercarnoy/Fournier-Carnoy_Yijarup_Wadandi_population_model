# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Data from project files.
# Task:    Set up info cubes for commercial fishing, to use in C++ function
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    March 2026

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

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))


## 0. files needed ------------------------------------------------------------

file_water          <- "data/output_data/03_B_water.rds"
file_network        <- "data/output_data/03_B_network_shapefile.shp"
file_wa             <- "data/output_data/01_B_land.shp"
file_bathy          <- "data/input_data/SW_crop_AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
file_boat_ramps_w   <- "data/input_data/wadandi_boat_ramps.shp"
file_boat_ramps_n   <- "data/input_data/north_boat_ramps.shp"

year_start <- 1900
year_end <- 2024
n_years_tot <- year_end - year_start + 1
current_fleet <- "boat_rec"

common_crs = 7850
bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -31), crs = common_crs) %>%
  st_as_sfc() %>%
  st_transform(crs = common_crs)

NCELL <- nrow(readRDS(file_water))


## 1. Elements useful to multiple elements of the script ----------------------

# -- several things can be calculated once and used for multiple output calculations.

### 1.1 Fishable depth --------------------------------------------------------

# -- this section of the fishing effort relates to boats being able to fish further over the years, with bigger boats, more powerful engines etc.
# -- the model takes this into account by adding a 'fishable depth' into the mix, which goes more and more offshore.
# -- for commercial fishing, we're assuming that on average, fishers can fish 1.3m more every year. (change as needed in parameters)
water <- readRDS(file_water) |> st_transform(common_crs)
centroids <- st_centroid(water %>% st_make_valid()) |> 
  st_transform(common_crs)
NCELL <- nrow(water)

# prepare grid cells
plot(water$geometry)
wa_mask <- st_read(file_wa); wa_mask <- st_as_sf(wa_mask) %>% st_transform(common_crs); plot(wa_mask, col = "lightgray", add = T)
plot(bbox, add = T)

# set parameters
years <- year_start:year_end
min_depth <- 20 # fishable depth is about 20m in 1945 (from Gaynor 2008, p.38)
increase_year <- 1950 # the year that boats start to go deeper (post-war industrialisation)
depth_per_year <- 1 # fishable depth gained per year (m) - arbitrary to reach the continental shelf quickly.

# prepare the bathymetry layer
bathy <- rast(file_bathy) %>%
  project(paste0("EPSG:", common_crs)) %>%
  crop(vect(water)) %>%
  abs() %>%
  terra::mask(vect(st_transform(wa_mask, common_crs)), inverse = TRUE)

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

# CHECK: fishable depth at different times
test_year = c(1, 50, 60, 90)
p <- ggplot(data =  water %>%
              mutate(as.data.frame(fishable_depth_cell_year[, test_year]) %>%
                       setNames(paste0("year_", test_year))) %>%
              tidyr::pivot_longer(
                cols = starts_with("year_"),
                names_to = "year",
                names_prefix = "year_",
                values_to = "test"
              )
) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[c(3, 5)]) +
  facet_wrap(~ year, ncol = length(test_year)) +
  theme_minimal()
ggsave("plots/script_plot_checks/04_C/04_C_fishable_depth_boat_rec.png", plot = p, width = 10, height = 6, dpi = 500)


### 1.2 Cell area -------------------------------------------------------------

water <- readRDS(file_water) |> st_transform(common_crs)
cell_area <- water$cell_area / 1000000 # just cell area in km2
cell_area <- as.numeric(cell_area)


### 1.3 Temporal fishability --------------------------------------------------

# -- this is a cell x month x year cube that tells you whether the cell is closed (0), partially open (>0, <1) or open (1).

water <- readRDS(file_water) |> st_transform(common_crs)
glimpse(water)

water[[current_fleet]] <- dplyr::case_when(
  water[[current_fleet]] %in% c("T", "TRUE")  ~ TRUE,
  water[[current_fleet]] %in% c("F", "FALSE") ~ FALSE,
  TRUE ~ NA
)

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
test_year = 120
test <- water_area[, 1, test_year]
ggplot(data = water %>% mutate(test = test)) +
  geom_sf(aes(fill = test), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()

test_year = c(90, 100, 110, 120)
test_month = 11
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
ggsave("plots/script_plot_checks/04_C/04_C_fishable_area_boat_rec.png", plot = p, width = 10, height = 6, dpi = 500)


## 2. Fisher cell attractivity ------------------------------------------------

# -- ultimately to distribute fishing effort we need to have a formula of this shape: 
# -- (travel cost * coef1) + (log(offshore dist) * coef2) + (log(cell area) * coef3)
# -- ‾‾‾‾‾‾‾‾‾‾‾‾                 ‾‾‾‾‾‾‾‾‾‾‾‾‾                  ‾‾‾‾‾‾‾‾‾
# --     2.3                           2.1                          1.2

# -- the coefficients we do not have from published literature, and will therefore be eyeballed.
# -- boat ramps get built over the years, so attractivity will only be non-zero after a build date. 


### 2.1 Shore distance --------------------------------------------------------

water <- readRDS(file_water) |> st_transform(common_crs)

# find shore. we will use this to calculate distance of all cells to shore.
shore <- water[water$type %in% c("shore_wadandi", "shore_north"),] %>%
  st_make_valid() %>% 
  st_union() %>% 
  st_transform(common_crs) %>% 
  st_as_sf()

network <- st_read(file_network)|> st_transform(common_crs); plot(network$geometry)

centroids <- st_centroid(water %>% st_make_valid()) |> st_transform(common_crs)

shore_dist <- st_distance(centroids, shore) / 1000 # distance from shore in km
shore_dist <- as.data.frame(shore_dist) %>% mutate(ID = water$ID, shore_dist = as.numeric(shore_dist))
glimpse(shore_dist)


# CHECK: cells' distance from shore
p <- ggplot(data = water %>% dplyr::select(!where(is.list)) %>% mutate(shore_dist = as.numeric(shore_dist$shore_dist))) + 
  geom_sf(aes(fill = shore_dist), col = NA) +
  scale_fill_gradientn(colours = colour_palette[c(3, 5)]) +
  theme_minimal()
ggsave("plots/script_plot_checks/04_C/04_C_shore_distance_commercial.png", plot = p, width = 6, height = 10, dpi = 500)


### 2.2 Ramp built ------------------------------------------------------------

# -- ramps are built over time, giving access to rec fishers. we need to account for these build dates.

water <- readRDS(file_water) |> st_transform(common_crs)

BR_n <- st_read(file_boat_ramps_n) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  glimpse()

BR_w <- st_read(file_boat_ramps_w) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  glimpse()

plot(water$geometry); plot(BR_w$geometry, col = colour_palette[5], pch = 16, cex = 1, add = TRUE)

glimpse(BR_n)
glimpse(BR_w)

BR <- rbind(BR_n %>% dplyr::select(name, build_year, build_mnth, geometry) %>% mutate(region = "north"), 
            BR_w %>% dplyr::select(name, build_year, build_mnth, geometry) %>% mutate(region = "wadandi")) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>% 
  glimpse()
BR <- BR[!st_is_empty(BR), ]
BR$build_year[is.na(BR$build_year)] <- year_start # fill unknown build years with the start year
BR$build_mnth[is.na(BR$build_mnth)] <- 1 # fill unknown build months with the start month

plot(water$geometry); plot(BR$geometry, col = colour_palette[6], pch = 16, add = TRUE)

BR <- st_as_sf(BR); BR <- BR %>% st_transform(common_crs) 

# built or not
ap_names <- sort(unique(BR$name))

built <- array(NA, dim = c(year_end - year_start + 1, 12, length(ap_names)))

for (ACCESS in 1:length(ap_names)) {
  ap_row <- BR %>% dplyr::filter(name == ap_names[ACCESS])
  build_year  <- ap_row$build_year[1]
  build_month <- ap_row$build_mnth[1]
  
  for (YEAR in 1:(year_end - year_start + 1)) {
    calendar_year <- year_start + YEAR - 1
    for (MONTH in 1:12) {
      built[YEAR, MONTH, ACCESS] <- (calendar_year >= build_year) |
        (calendar_year == build_year & MONTH >= build_month)
    }
  }
}

built[built == FALSE] = NA

### 2.3 Travel cost -----------------------------------------------------------

# -- here we account for the fact that fishers are limited in how far they can afford to go.
# -- for commercial fishermen, cells close to the ramp, but far from shore are preferred.


#### 2.3.1 Distance from access points ----------------------------------------

# -- the intended output here is an object of size cell x access point

water <- readRDS(file_water) |> st_transform(common_crs)

BR_n <- st_read(file_boat_ramps_n) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  glimpse()

BR_w <- st_read(file_boat_ramps_w) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>%
  glimpse()

plot(water$geometry); plot(BR_w$geometry, col = colour_palette[5], pch = 16, cex = 1, add = TRUE)

glimpse(BR_n)
glimpse(BR_w)

BR <- rbind(BR_n %>% dplyr::select(name, build_year, build_mnth, geometry) %>% mutate(region = "north"), 
            BR_w %>% dplyr::select(name, build_year, build_mnth, geometry) %>% mutate(region = "wadandi")) %>% 
  st_transform(common_crs) %>%
  st_make_valid() %>% 
  glimpse()
BR <- BR[!st_is_empty(BR), ]
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

net <- activate(network, "nodes")  %>% st_transform(common_crs) 

plot(centroids$geometry)

# measure the distance from access points to cell centroids
network_matrix <- st_network_cost(net, from = BR, to = points_sf)
dim(network_matrix) # number of ramps x number of cells

glimpse(network_matrix)
access_dist <- as.data.frame(t(network_matrix))
colnames(access_dist) <- BR$name
access_dist$ID <- water$ID
access_dist <- access_dist / 1000
head(access_dist) # this gives us each cell's distance to the access points


# CHECK: distance from access point
access_dist <- access_dist %>% mutate(ID = centroids$ID)
water_dist_long <- water %>% left_join(access_dist, by = "ID") %>% pivot_longer(cols = BR$name, names_to = "Ramp", values_to = "Distance_km")
BR_2 <- BR %>% rename(Ramp = name)
p <- ggplot() +
  geom_sf(data = water_dist_long %>% dplyr::filter(Ramp %in% unique(water_dist_long$Ramp)[1:10]),
          aes(fill = as.numeric(Distance_km)), color = NA) +
  geom_sf(data = BR_2 %>% dplyr::filter(Ramp %in% unique(water_dist_long$Ramp)[1:10]),
          shape = 21, fill = "red") +
  scale_fill_gradientn(colours = colour_palette[c(3, 5)]) +
  facet_wrap(~ Ramp, ncol = 5) +
  labs(fill = "Distance (km)") +
  theme_minimal()
ggsave("plots/script_plot_checks/04_C/04_C_access_point_distance_boat_rec.png", plot = p, width = 10, height = 10, dpi = 500)


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


### 2.4 Calculate cell attractivity -------------------------------------------

# -- okay so now that we have most of our elements, we can make a function that defines the 'attractivity' of each cell
# -- the format should be: cell x access point x year

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
      attractivity[, ACCESS, MONTH, YEAR][water_area[, MONTH, YEAR] == 0] <- NA
      
      # and turn off access points if they haven't been built yet
      attractivity[, ACCESS, MONTH, YEAR] <- built[YEAR, MONTH, ACCESS] * attractivity[, ACCESS, MONTH, YEAR]
    }
  }
}
tail(attractivity[,1,12,1])

# CHECK: attractivity of cells from the POV of few access points
test_month <- 12
test_access <- c(6, 10, 29)
test_years <- c(1, 50, 71, 120)

plot_data <- purrr::map_dfr(test_access, function(a) {
  purrr::map_dfr(test_years, function(y) {
    water %>%
      mutate(attractivity = attractivity[, a, test_month, y],
             access_id = a,
             year = y)
  })
})

access_labels <- setNames(
  paste0(BR_2$Ramp[test_access], " (built ", BR_2$build_year[test_access], ")"),
  test_access
)

p <- ggplot() +
  geom_sf(data = plot_data, aes(fill = attractivity), colour = NA) +
  geom_sf(data = BR_2[test_access, ] %>% mutate(access_id = test_access),
          shape = 21, fill = "red") +
  scale_fill_gradientn(colours = colour_palette[c(3, 5)]) +
  facet_grid(access_id ~ year,
             labeller = labeller(access_id = as_labeller(access_labels))) +
  theme_minimal()
ggsave("plots/script_plot_checks/04_C/04_C_cell_attractivity_boat_rec.png", plot = p, width = 7, height = 10, dpi = 500)

# transform in a format that C++ can deal with easier
attractivity <- lapply(1:12, function(m) {
  arr <- attractivity[, , m, , drop = FALSE]   # keep 4D, dim3 = 1
  dim(arr) <- dim(attractivity)[-3]             # collapse only the month axis -> (NCELL, n_access, n_years)
  arr
})


## 3. Boat rec catchability ---------------------------------------------------

# -- catchability is the susceptibility of fish to be caught by 1 unit of effort. 
# -- catchability differs for every cell, given that fish is more susceptible to being caught by 1 unit of effort in a cell of 10m2 compared to 1 unit of effort in a cell of 100km2.

# -- catchability works with fishing effort to produce fishing mortality. 
# -- it is unknowable, so we will need to calibrate it. 
# -- see script 04_D for process 

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


# CHECK: catchability 
test_year <- c(50, 100, 120)
catch_mat <- catchability[, 1, test_year]  # matrix: cells x years
catch_mat[catch_mat > 0.1] <- NA # some cells have extremely high catchability and make the plot flat - remove them to see fresh.
p <- ggplot(data = water %>%
              bind_cols(as.data.frame(catch_mat) %>%
                          setNames(paste0("year_", test_year))) %>%
              tidyr::pivot_longer(
                cols = starts_with("year_"),
                names_to = "year",
                names_prefix = "year_",
                names_transform = as.numeric,
                values_to = "catchability"
              )) +
  geom_sf(aes(fill = catchability), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[c(3, 5)]) +
  facet_wrap(~ year) +
  ggtitle("cells with catchability <0.1 are shown, for display purposes. these cells which end up with high catchability due to small size") +
  theme_minimal()
ggsave("plots/script_plot_checks/04_C/04_C_catchability_boat_rec.png", plot = p, width = 10, height = 10, dpi = 500)


## 4. Fishing effort values ---------------------------------------------------

# -- we know how much fishing there is (roughly) but we now need to distribute it across access points.


### 4.1 Enter the overall effort values (boat days) ---------------------------

# obtain literature values
years_full <- year_start:year_end
g_sheets <- read.csv("data/input_data/YIJARUP - Fishing effort reconstruction - FINAL_fishing_effort (14-07-2026).csv", skip = 3) %>% 
  dplyr::select(c(YEAR, boat.days.west.coast)) %>% # select only commercial fishing columns.
  glimpse()

# fill the unknown year values with linear interpolation
annual_effort_boat <- g_sheets %>% 
  rename(year = YEAR,
         boat_days = boat.days.west.coast)
annual_effort_boat$boat_days <- as.vector(approx(annual_effort_boat$boat_days, n = nrow(annual_effort_boat))$y)

# CHECK: commercial effort
plot(annual_effort_boat$year, annual_effort_boat$boat_days, type = "l", col = colour_palette[5], lwd = 4,
     xlab = "Year", ylab = "Boat Days")


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
boat_effort <- expand.grid(
  year = years_full,
  month = sprintf("%02d", 1:12)
) %>%
  arrange(year, month) %>%
  mutate(
    annual_boat_days = rep(annual_effort_boat$boat_days, each = 12),
    monthly_effort = annual_boat_days * seasonal_multipliers[month]
  ) %>%
  dplyr::select(year, month, monthly_effort)

# CHECK: monthly boat recreational effort
ggplot(boat_effort, aes(x = as.Date(paste(year, month, "01", sep = "-")), y = monthly_effort)) +
  geom_line(color = colour_palette[4]) +
  labs(title = paste0("Monthly Commercial Fishing Effort ", year_start, "-", year_end),
       x = "Date", y = "Monthly Boat Days") +
  theme_minimal() +
  geom_smooth(color = colour_palette[5])

# proportion of each month's contribution to yearly boat days
boat_month_prop <- boat_effort %>% 
  group_by(year) %>% 
  mutate(year_sum = sum(monthly_effort)) %>%
  mutate(month_prop = monthly_effort/year_sum) %>% 
  dplyr::select(-year_sum)

boat_month_prop <- boat_month_prop %>% 
  group_by(month) %>% 
  mutate(ave_month_prop = mean(month_prop))
prop_month_ave <- boat_month_prop[1:12, c(2, 5)]

saveRDS(prop_month_ave, "data/output_data/04_C_boat_rec_prop_month_ave.rds") # charlotte's 'Average_Monthly_Effort"


### 4.3 Distribute monthly effort into boat ramps -----------------------------

# obtain access point information
BR_n <- st_read(file_boat_ramps_n) %>% 
  st_transform(4283) %>%
  st_make_valid() %>%
  mutate(build_year = as.numeric(build_year),
         build_year = ifelse(is.na(build_year), year_start, build_year), # fill in missing dates with the start year
         build_mnth = ifelse(is.na(build_mnth), 1, build_mnth),     # assume Jan if unknown
         build_mnth = as.numeric(build_mnth)
  ) %>%
  dplyr::select(-c(com_prop, ABS_name)) |> 
  glimpse()
BR_n <- BR_n[!st_is_empty(BR_n), ]

BR_w <- st_read(file_boat_ramps_w) %>% 
  st_transform(4283) %>%
  st_make_valid() %>%
  mutate(build_year = as.numeric(build_year),
         build_year = ifelse(is.na(build_year), year_start, build_year), # fill in missing dates with the start year
         build_mnth = ifelse(is.na(build_mnth), 1, build_mnth),     # assume Jan if unknown
         build_mnth = as.numeric(build_mnth)
  ) %>%
  dplyr::select(-c(X1943, X1965, X1981, X2024, rec_visits, visit_2024, com_prop)) |> 
  glimpse()
BR_w <- BR_w[!st_is_empty(BR_w), ]

BR <- rbind(st_transform(BR_n, common_crs), st_transform(BR_w, common_crs)) |> 
  glimpse()

# distribute effort across ramps and time
glimpse(boat_effort)
boat_effort <- boat_effort %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) # make a column with year and month of ramp build
ramp_effort <- expand.grid(
  ramp_index = 1:nrow(BR),
  date = boat_effort$date
) %>%
  mutate(
    build_date = as.Date(paste(BR$build_year[ramp_index], BR$build_mnth[ramp_index], "01", sep = "-")),
    ramp_name = BR$name[ramp_index]
  ) %>%
  filter(date >= build_date) %>%
  mutate(
    months_since_build = interval(build_date, date) %/% months(1),
    logistic_growth = 1 / (1 + exp(-0.1 * (months_since_build - 60))),
    ramp_weight = logistic_growth
  )

# merge in monthly effort
ramp_effort <- ramp_effort %>%
  left_join(boat_effort, by = "date") %>%
  group_by(date) %>%
  mutate(
    total_weight = sum(ramp_weight),
    adjusted_effort = ifelse(total_weight > 0, monthly_effort * (ramp_weight / total_weight), 0),
    year = year(date),
    month = month(date)
  ) %>%
  ungroup() %>%
  dplyr::select(year, month, boat_ramp = ramp_name, adjusted_effort)

# check total per month equals monthly_effort
check_totals <- ramp_effort %>%
  group_by(year, month) %>%
  summarise(total_effort = sum(adjusted_effort), .groups = "drop") %>%
  left_join(boat_effort %>% mutate(year = year(date), month = month(date)), by = c("year", "month")) %>%
  mutate(diff = abs(total_effort - monthly_effort))
summary(check_totals$diff)  # should be near zero

# check each boat ramp is populated correctly
ggplot(ramp_effort, aes(x = year, y = adjusted_effort)) +
  geom_line(color = colour_palette[5], lwd = 1) +
  facet_wrap(~boat_ramp, ncol = 5)

# cut up this big unwieldy dataframe into a cube
access_point_effort <- ramp_effort %>% 
  pivot_wider(names_from = boat_ramp,
              values_from = adjusted_effort,
              values_fill = 0) %>% 
  as.data.frame()

access_point_effort <- access_point_effort %>%
  dplyr::select(-year, -month) %>%
  {abind::abind(split(., access_point_effort$year), along = 3)}

# rename dimension names so they're consistent with the rest of the objects used in the function.
dimnames(access_point_effort)[[1]] <- c(1:12)
dimnames(access_point_effort)[[3]] <- c(1:length(dimnames(access_point_effort)[[3]]))
access_point_effort <- access_point_effort[,order(colnames(access_point_effort)),] # arrange access point columns alphabetically to match other objects

# check that utility and access_point_effort have the same order of access point
names(access_dist %>% dplyr::select(!ID)) == dimnames(access_point_effort)[[2]]
dimnames(access_point_effort)[[2]] <- c(1:length(dimnames(access_point_effort)[[2]]))

# utility # watch out that the order of ramp effort and of utility match: otherwise you're assigning the effort to the wrong ramp in the function


## 5. Create a list to use in the C++ function --------------------------------

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


## 6. Setting up effort for the burn-in ---------------------------------------

# -- the burn-in exists only to stabilise the population before the simulations start.
# -- to stabilise the population, we need to run the function with a small amount of fishing (smaller than what we start with)

# obtain the relative contribution of each boat ramp
effort_prop <- ramp_effort %>%
  dplyr::filter(year == 1900, 
                month == 1) %>% 
  dplyr::mutate(com_prop = adjusted_effort/sum(adjusted_effort)) %>% 
  dplyr::select(boat_ramp, com_prop, adjusted_effort) %>% 
  glimpse()
barplot(effort_prop$com_prop)

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


# split up by boat ramp
for(AP in 1:(ncol(burn_in_effort)-1)){
  burn_in_effort[AP+1]= burn_in_effort$effort*effort_prop[AP, "com_prop"]
}

burn_in <- as.matrix(burn_in_effort[-1])

saveRDS(burn_in, file = "data/output_data/04_C_boat_rec_burn_in_effort.rds")
# for the burn in, just replace the corresponding list object from section 5 with this one.


### END ###
