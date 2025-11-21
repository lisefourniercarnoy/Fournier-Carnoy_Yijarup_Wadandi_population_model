# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    ?
# Task:    Create a simulated fishing effort with NO MANAGEMENT
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    September 2025

# -----------------------------------------------------------------------------

# Status: I think finished, unless fishing effort changes after 17/10/2025

# -----------------------------------------------------------------------------

rm(list = ls()) # clear environment

library(tidyverse) # data manipulation
library(sf) # shapefile manipulation
library(raster) # raster manipulatio
library(RColorBrewer) # plotting colours
library(abind) # matrix manipulation
library(sfnetworks) # for making distance to access points
library(exactextractr) # estracting raster values at points

# The point of this script is to set up a fishing effort surface over time as if 
# no Sanctuary Zones were in place.
# This is essentially a copy of the 04 fishing effort scripts, removing the NTZ
# sections.
# All 3 fishing effort sources are present here.

## Custom plotting parameters -------------------------------------------------
colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))
source("custom_theme.R")


# 04_A COMMERCIAL -------------------------------------------------------------
## 0. Files used in this script -----------------------------------------------

file_wa         <- "data/output_data/01_wadandi_land.shp"
file_bathy      <- "data/input_data/SW_crop_AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
file_boat_ramps <- "data/input_data/wadandi_boat_ramps.shp"
file_water      <- "data/output_data/03_water.rds"
file_network    <- "data/output_data/03_network_shapefile.shp"

file_com_fishable <- "data/output_data/04A_commercial_fishable_area_over_time.shp"
file_com_cell_dist <- "data/output_data/04A_commercial_cell_dist.rds"
file_com_BR_trips <- "data/output_data/04A_BR_trips.rds"
file_com_boat_days <- "data/output_data/04A_commercial_total_boat_days.rds"
file_com_ramp_effort <- "data/output_data/04A_commercial_ramp_effort.rds"

year_start <- 1945
year_end <- 2024
n_years_tot <- year_end - year_start +1
n_years_pre18 <- 2018 - year_start

crs_raster <- "+proj=longlat +datum=WGS84 +no_defs"
bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -33.2), crs = crs_raster)


## 1. BOAT Catchability -------------------------------------------------------

# original catchability has pre- and post- NTZ values. we'll remake a single unified catchability, as S00 has no management

water <- readRDS(file_water) %>% filter(!is.na(ID)); plot(water)

# We'll calculate the catchability of each cell by its area and whether it's no-take or not.
water <- water %>% 
  mutate(Area = as.vector((water$cell_area)/1000000))

# all cells have the same catchability over time. names kept from original script for clarity
water_area <- water %>% 
  dplyr::select(status, Area) %>% 
  mutate(area_2000 = Area) %>%
  dplyr::select(area_2000) %>% 
  mutate(sum_2000 = sum(area_2000)) %>% 
  st_drop_geometry() %>% 
  mutate(q_2000 = area_2000/sum_2000) %>% 
  mutate(ID = row_number())

# Create an array of catchability for each cell (rows) and each year (columns)
NCELL <- nrow(water)
spatial_q <- array(0.000006, dim = c(NCELL, n_years_tot)) # Why is the original catchability set to 0.000006 ?

# The following bit is for increasing catchability after 1990 
for (COL in ((1990-1945)+1):n_years_tot) {
  spatial_q[, COL] <- spatial_q[, COL-1] * 1.02
} # This increases cells by 2% each year, i think to account for better tech

for (COL in 1:n_years_tot) {
  for (ROW in 1:NCELL){
    spatial_q[ROW, COL] <- spatial_q[ROW, COL] / water_area[ROW, "q_2000"]
  }
} # For all cells, catchability is in proportion to the cell's area.

spatial_q[spatial_q == Inf] <- 0 # IDK what this does.
summary(spatial_q) # this is a matrix that tells us for each cell (each row), and each year (each column), how likely you'd catch a fish in that cell based on how big it is.

# Plot check
catch_df <- as.data.frame(spatial_q)
colnames(catch_df) <- paste0("Year_", year_start:(year_start + n_years_tot-1))
catch_df$ID <- water_area$ID
water_catch <- water %>% left_join(catch_df, by = "ID")
water$ID <- water_area$ID
water_long <- water_catch %>% pivot_longer(cols = starts_with("Year_"), names_to = "Year", names_prefix = "Year_", values_to = "Catchability") %>% mutate(Year = as.numeric(Year))
ggplot(water_long[water_long$Year %in% c(1975, 2020),]) +
  geom_sf(aes(fill = log(Catchability)), color = NA) +
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  facet_wrap(~ Year, ncol = 2) +
  labs(title = "S00 - Boat Catchability", fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/S00_boat_catchability.png", plot = last_plot())

# Save for future use.
saveRDS(spatial_q, file = paste0("data/output_data/S00_spatial_q_NTZ.rds"))


## 5. Create a utility function -----------------------------------------------

## Now need to create a separate fishing surface for each month of each year based on distance to boat ramp, size of each cell,
## and multiply that by the effort in the cell to spatially allocate the effort across the area. Effort is also able to go more offshore over time.
## But we need to account for the fact that there will be sanctuary zones going in and the effort that would have gone in there will get allocated somewhere else
## Will then need to put the rows/columns back in as 0s 

# What we want to do is to distribute effort across month and year based on :
# 1. distance to boat ramp, 2. size of the cell, and 3. whether the cell is 'fishable' that year and 4. fishing effort
# After the SZ comes in the effort will also be redistributed.

fishable_long_sf <- st_read(file_com_fishable); names(fishable_long_sf) <- c("ID", "year", "fishable_prop", "type", "sand", "reef", "seagrass", "status", "cell_area", "Area", "geometry")
Cell_Dist <- readRDS(file_com_cell_dist)
BR_trips <- readRDS(file_com_BR_trips)

years <- sort(unique(fishable_long_sf$year))
ramps <- unique(BR_trips$boat_ramp)
nramps <- length(ramps)

# Get full list of cell IDs from fishable_long_sf
all_cell_ids <- sort(unique(fishable_long_sf$ID))
ncells <- length(all_cell_ids)
nyears <- length(years)

BR_U_array <- array(0, dim = c(ncells, nramps, nyears),
                    dimnames = list(cell_id = all_cell_ids,
                                    ramp = ramps,
                                    year = as.character(years)))

for (y in seq_along(years)) {
  yr <- years[y]
  
  # Step 2.1: Get fishable surface for the year
  fish_yr <- fishable_long_sf %>%
    filter(year == yr) %>%
    mutate(fishable_prop = if_else(
      status %in% c("NTZ_boat", "NTZ_boat_shore") & year >= 2019, 0, fishable_prop)) %>%
    dplyr::select(ID, fishable_prop)
  
  # Step 2.2: Join distances, areas, and fishable_prop
  Vj_df <- Cell_Dist %>%
    inner_join(fish_yr, by = c("cell_id" = "ID")) %>%
    arrange(cell_id)
  
  # Step 2.3: Calculate raw utility: exp(-distance) * area * fishable proportion
  # We'll loop over ramps
  U_mat <- matrix(0, nrow = nrow(Vj_df), ncol = nramps)
  
  for (r in seq_along(ramps)) {
    ramp_name <- ramps[r]
    U_mat[, r] <- exp(-Vj_df[[ramp_name]]) * Vj_df$Area * Vj_df$fishable_prop # exp(-Vj_df...) because otherwise high utility is given to areas far from ramps
  }
  
  # Step 2.4: Normalize utilities so they sum to 1 per ramp
  U_norm <- sweep(U_mat, 2, colSums(U_mat, na.rm = TRUE), "/")
  
  # Step 2.5: Store in the 3D array
  # Align to full cell list
  row_ids <- match(Vj_df$cell_id, all_cell_ids)
  BR_U_array[row_ids, , y] <- U_norm
}
plot(BR_U_array[,,70])
head(BR_U_array[,,70])
head(water)

br_slice <- BR_U_array[,,70]

# Calculate row sums (sum of utilities across ramps for each cell)
row_sums <- rowSums(br_slice)

# Add the sums as a new column to the water sf object
water$utility_sum <- row_sums

# Plot using ggplot2, coloring polygons by the sum of utilities
ggplot(water) +
  geom_sf(aes(fill = utility_sum), color = NA) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", 
                       na.value = "grey80", name = "Sum of Utilities") +
  theme_minimal() +
  labs(title = "Sum of Utilities Across Ramps per Cell",
       subtitle = "Layer 70 of BR_U_array") +
  theme(legend.position = "right")



# plot check
ramps
ramp_check <- "GB_Bunbury_Stirling_St"
utility_check <- as.data.frame(BR_U_array[, ramp_check, ])  # dimensions: cells × years
utility_check$cell_id <- as.integer(rownames(utility_check))
head(utility_check)

water_catch <- water %>% mutate(cell_id = row_number()) %>% left_join(utility_check, by = "cell_id")
water_catch_long <- water_catch %>% pivot_longer(cols = as.character(years), names_to = "year", values_to = "Catchability")
ggplot(water_catch_long) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  facet_wrap(~ year, ncol = 16) + 
  labs(title = paste0("commercial catchability of each cell, \nby distance from ", ramp_check, " ramp, cell size and tech-fishability"), fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/S00_commercial_boat_ramp_catchability_over_time.png", plot = last_plot())

# Plot check (notice the difference between popular/non-popular ramps)
utility_check2 <- as.data.frame(BR_U_array[,,80]) %>% mutate(cell_id = as.numeric(rownames(BR_U_array[,,80])))
BR <- st_read(file_boat_ramps) %>% 
  dplyr::filter(name %in% c("SC_Augusta_Ellis_St_Jetty", "WC_Gnarabup", "WC_Hamelin_Bay",
  "GB_Quindalup", "GB_Eagle_Bay", "GB_Bunbury_Stirling_St", "GB_Busselton_Georgette_Street")); plot(BR$geometry)

water_catch <- water %>% mutate(cell_id = row_number()) %>% left_join(utility_check2, by = "cell_id")
water_catch_long <- water_catch %>% pivot_longer(cols = BR$name, names_to = "Ramp", values_to = "Catchability")
ggplot(water_catch_long) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  facet_wrap(~ Ramp, ncol = 4) + 
  labs(title = "Catchability of each cell in 2024, by boat ramp distance and cell size - \nnotice differences between popular/non popular ramps", fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/S00_commercial_catchability_surface_by_ramp.png", plot = last_plot())


## 6. Allocating effort to cells ----------------------------------------------

# We know the fishing effort over the years
# We know the fishing effort across months
# We know how 'useful' each cell is for fishing (utility function)
# Now we want to allocate fishing effort to cells over time

BR_trips <- BR_trips %>%
  dplyr::select(num_year, month, boat_ramp, adjusted_effort) %>%
  pivot_wider(
    names_from = boat_ramp,
    values_from = adjusted_effort,
    values_fill = 0  # Fill missing effort with 0s
  ) %>%
  arrange(num_year, month)


# now assign to cells

c_fishing <- array(0, dim = c(NCELL, 12, length(years))) # this array has a row for every cell, a column for every month, and a layer for every year
months <- array(0, dim = c(NCELL, 12)) # array for the number of months
ramps <- array(0, dim = c(NCELL, length(BR$name))) # array for the number of ramps
layer <- 1

head(BR_trips)
head(BR)
head(BR_U_array[,,1])

for(YEAR in 1:length(years)){ # for all years,
  
  for(MONTH in 1:12){ # run the loop for every month,
    
    for(RAMP in 1:length(BR$name)){ # and for every ramp,
      
      temp <- BR_trips %>% 
        filter(num_year == YEAR) %>% 
        dplyr::select(-c(num_year, month))
      
      temp <- as.matrix(temp)
      for(CELL in 1:NCELL){ # and assign fishing effort to every cell, based on its utility
        ramps[CELL,  RAMP] <- BR_U_array[CELL, RAMP, YEAR] * temp[MONTH, RAMP]
      }
    }
    months[, MONTH] <- rowSums(ramps)
  }
  c_fishing[ , , layer] <- months 
  layer <- layer + 1
} # this loop assigns each cell a fishing effort based on utility by month
c_fishing

# plot check
head(c_fishing[,,80])

water
year_idx <- 75  # first year
month_idx <- 1  # January

# Extract effort vector (one value per cell)
effort_vec <- c_fishing[, month_idx, year_idx]

# Add to water polygons
water$effort <- effort_vec

# Plot
ggplot(water) +
  geom_sf(aes(fill = log(effort))) +
  scale_fill_viridis_c() +
  labs(title = paste("Fishing Effort - Year", year_idx, "Month", month_idx),
       fill = "log Effort") +
  theme_minimal()

## Set up effort for burn-in --------------------------------------------------

# boat 
location_name <- unique(BR$name)

BR_trips <- readRDS(file_com_ramp_effort) %>%
  group_by(boat_ramp) %>%
  summarise(boat_days = sum(adjusted_effort, na.rm = TRUE)) %>%
  arrange(desc(boat_days))
BR_trips$effort <- 1

BR_trips <- BR_trips %>%
  mutate(trip_per_hr = as.numeric(unlist((boat_days / effort)))) %>% # Standardise the no. trips based on how much time you spent sampling
  mutate(BR_prop = trip_per_hr/sum(trip_per_hr)) #Then work out the proportion of trips each hour that leave from each boat ramp

# Fishing parameters
eq.init.fish = 0.025
q = 0.00001
effort = (-log(1 - eq.init.fish)) / q # We assume the same level of nominal effort in each year

## Split up this effort by the same proportions as before and allocate it to the different access points
# Months
burn_in_effort <- prop_month_ave[, 2] * effort

burn_in_effort <- as.data.frame(burn_in_effort) %>%
  rename(effort = "ave_month_prop")
for (location in location_name) {
  burn_in_effort[[location]] <- 0
} # Loop through each location name and create a new column with 0's

for (M in 1:12) { # for each month,
  for (i in 1:nrow(BR_trips)) { # and each boat ramp,
    location_name <- BR_trips$boat_ramp[i] # extract the location name
    effort_value <- BR_trips$effort[i]  # and effort value
    burn_in_effort[[location_name]] <- effort_value * BR_trips$BR_prop[i] # and calculate burn in effort
  }
} # this loop calculate burn-in effort for all boat ramps, based on how many visits and hours each ramp gets fished from


## Allocate to the cells using the same utilities that we set up earlier
b_burn_in_fishing <- array(0, dim = c(NCELL, 12, n_years_tot)) #This array has a row for every cell, a column for every month, and a layer for every year
months <- array(0, dim = c(NCELL, 12))
ramps <- array(0, dim = c(NCELL, length(location_name)))
yr <- 1

for(YEAR in 1:n_years_tot){
  
  print(YEAR)
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:length(location_name)){
      
      temp <- burn_in_effort %>%
        dplyr::select(-c(effort))
      
      temp <- as.matrix(temp)
      
      for(CELL in 1:NCELL){
        ramps[CELL, RAMP] <- BR_U_array[CELL, RAMP, YEAR] * temp[MONTH, RAMP]
      }
    }
    
    months[, MONTH] <-  rowSums(ramps)
  }
  b_burn_in_fishing[ , , yr] <- months
  b_burn_in_fishing[ , , yr] <- b_burn_in_fishing[ , , yr] * spatial_q[, 1]
  yr <- yr + 1
} # this loop assigns each cell an amount of fishing effort each month and across years based on how useful it is (how far from a boatramp and how frequented that boat ramp is)


saveRDS(b_burn_in_fishing, file = "data/output_data/S00_commercial_burn_in_fishing.rds")


# 04_B SHORE FISHING ----------------------------------------------------------
## 0. Files used in this script -----------------------------------------------

file_wa         <- "data/output_data/01_wadandi_land.shp"
file_ntz        <- "data/output_data/01_wadandi_NTZ.shp"
file_bathy      <- "data/input_data/wadandi_250m_bathy.tif"
file_carpark    <- "data/input_data/wadandi_carparks.shp"
file_water      <- "data/output_data/03_water.rds"
file_network    <- "data/output_data/03_network_shapefile.shp"

file_shore            <- "data/output_data/04B_shore_fishable_area_over_time.shp"
file_shore_cell_dist  <- "data/output_data/04B_shore_cell_dist.rds"
file_shore_CP_trips   <- "data/output_data/04B_CP_trips.rds"
file_shore_boat_days  <- "data/output_data/04B_shore_total_boat_days.rds"
file_carpark          <- "data/input_data/wadandi_carparks.shp"
file_shore_effort <- "data/output_data/04B_shore_effort.rds"

# load land and ntz for sanity checks throughout
land <- st_read(file_wa); plot(land$geometry)
ntz <- st_read(file_ntz); plot(ntz$geometry, add = T)

year_start <- 1945
year_end <- 2024
n_years_tot <- year_end - year_start +1
n_years_pre18 <- 2018 - year_start

crs_raster <- "+proj=longlat +datum=WGS84 +no_defs"
bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -33.2), crs = crs_raster)


## 1. Catchability ------------------------------------------------------------

water <- readRDS(file_water) %>% filter(!is.na(ID)); plot(water)
CP <- st_read(file_carpark)
shore_cells <- water$ID[water$type == "shore"]

# We'll calculate the catchability of each cell by its area and whether it's no-take or not.
water <- water %>% 
  mutate(area = as.vector((water$cell_area)/1000000))

# cells in all years will have the same catchability over time
water_area <- water %>% # this is for all cells, pelagic and shore
  dplyr::select(ID, status, area, type) %>% 
  mutate(area_2000 = ifelse(!(ID %in% shore_cells), 0, area)) %>%  # make offshore cells unfishable for all years
  dplyr::select(ID, area_2000) %>%
  mutate(sum_2000 = sum(area_2000)) %>% 
  st_drop_geometry() %>% 
  mutate(q_2000 = area_2000/sum_2000) %>% 
  glimpse()

# create an array of catchability for each cell (rows) and each year (columns)
NCELL <- length(shore_cells)
spatial_q <- array(0.000006, dim = c(NCELL, n_years_tot)) # Why is the original catchability set to 0.000006 ?

# index of the q_2000 column
ix_q2000 <- which(colnames(water_area) == "q_2000")

for (ROW in 1:NCELL) {
  id <- shore_cells[ROW]
  i_area <- match(id, water_area$ID)
  
  q_2000 <- water_area[i_area, ix_q2000]

  # Pre-NTZ (up to 2017)
  for (COL in 1:n_years_tot) {
    spatial_q[ROW, COL] <- spatial_q[ROW, COL] / q_2000
  }
} # this is a cleaner version of Charlotte's loops. there is no increase in catchability over time for now.
summary(spatial_q) # this is a matrix that tells us for each cell (each row), and each year (each column), how likely you'd catch a fish in that cell based on how big it is.


# Plot check
catch_df <- as.data.frame(spatial_q)
catch_df[catch_df == 0] <- NA # this is just for visualisation purposes, to show that there is no fishing in SZ after 2018

colnames(catch_df) <- paste0("Year_", year_start:(year_start + n_years_tot-1))
catch_df$ID <- water_area$ID[water_area$ID %in% water$ID[water$type == "shore"]]  # or just 1:NCELL if they match in order
water_catch <- water %>% left_join(catch_df, by = "ID")  # ID must be in 'water' too
water$ID <- water_area$ID  # or use `row_number()`
water_long <- water_catch %>% pivot_longer(cols = starts_with("Year_"), names_to = "Year", names_prefix = "Year_", values_to = "Catchability") %>% mutate(Year = as.numeric(Year))
ggplot(water_long[water_long$Year %in% c(2000),]) +
  geom_sf(data = land, fill = "lightgray", col = NA) +
  geom_sf(data = ntz, fill = "#E4F2FF", col = NA) +
  geom_sf(aes(fill = log(Catchability + 1e-06)), color = NA) + # adding a small value to log(catchability) because log(0) = -Inf and doesn't display well
  scale_fill_gradientn(colours = colour_palette[4:6], na.value = NA) +
  facet_wrap(~ Year, ncol = 2) +
  labs(title = "Catchability over time (pre- and post- NTZ)", 
       fill = "log(Catchability)") +
  theme_minimal()
ggsave("plots/checking_plots_during_setup/S00_shore_rec_catchability.png", plot = last_plot())

# Save for future use.
saveRDS(spatial_q, file = paste0("data/output_data/S00_shore_rec_spatial_q_NTZ.rds"))



## 4. Set up a utility function -----------------------------------------------
DistCP <- readRDS(file_shore_cell_dist)
water_shore <- water[water$type == "shore",]

Cell_Vars <- DistCP %>% 
  mutate(Area = as.vector((water_shore$cell_area)/1000000)) # Cells are now in km^2 but with no units

## Now need to create a separate fishing surface for each month of each year based on distance to access point, size of each
## cell and multiply that by the effort in the cell to spatially allocate the effort across the area
## But we need to account for the fact that there will be sanctuary zones going in and the effort that would have gone in there will get allocated somewhere else
## Will then need to put the rows/columns back in as 0s 

# What we want to do is to distribute effort across month and year based on :
# 1. distance to access point, 2. size of the cell, and 3. fishing effort.
# After the SZ comes in the effort will also be redistributed.

NCELL_pre18 <- nrow(water_shore) # number of cells you can fish in- only shore cells

Vj <- Cell_Vars %>%
  mutate(vj = rowSums(across(all_of(CP$area)), na.rm = TRUE)) %>% 
  glimpse()

Vj_pre18 <- Vj # %>% mutate(across(where(is.numeric), ~ifelse(is.na(.), 0, .))) # replace NAs with 0
# catchable cells before SZ

area_col <- grep("Area", colnames(Vj)) # extract the column index for area - important for calculations. in theory it should be the same pre- and post-NTZ.


# pre-SZ - NO NTZ BUT CANT BE BOTHERED TO CHANGE NAMES

CP_U_pre18 <- as.data.frame(matrix(0, nrow = NCELL, ncol = length(unique(CP$area)))) # Set up data frame to hold utilities of cells
colnames(CP_U_pre18) <- c(CP$area)

cellU <- matrix(NA, ncol = length(unique(CP$area)), nrow = NCELL)

for(ACCESS in 1:length(unique(CP$area))){
  for(cell in 1:NCELL){
    U <- exp(Vj_pre18[cell, ACCESS] + log(Vj_pre18[cell, area_col]))
    U <- ifelse(is.na(U), 0, U) # many cells are NA, so the calculation above will return NA. Replace them with 0 so the line below works.
    cellU[cell, ACCESS] <- U
  }
} # this loop goes over each cell for each ramp, and calculates how catchable each cell is based on how close it is to a access point and how popular that access point is.

rowU <- as.data.frame(colSums(cellU))

for(ACCESS in 1:length(unique(CP$area))){
  for(cell in 1:NCELL_pre18){
    CP_U_pre18[cell, ACCESS] <- (exp(Vj_pre18[cell, ACCESS]+log(Vj_pre18[cell, area_col])))/rowU[ACCESS, 1]
  }
} # this loop goes over each cell for each ramp, and calculates how catchable (in %) each cell is based on how close it is to a access point and how popular that boat ramp is.
colSums(CP_U_pre18, na.rm = T) # all adds up to 1, perfect.
head(CP_U_pre18)

# Plot check
CP_U_pre18 <- bind_cols(ID = water_shore$ID, CP_U_pre18)
water_catch <- water_shore %>% mutate(ID = shore_cells) %>% left_join(CP_U_pre18, by = "ID")
water_catch_long <- water_catch %>% pivot_longer(cols = CP$area, names_to = "access_point", values_to = "Catchability")
ggplot(water_catch_long %>% filter(access_point %in% unique(water_catch_long$access_point)[51:52])) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  facet_wrap(~ access_point, ncol = 16) + 
  labs(title = "(before 2018) catchability of each cell, by shore access distance and cell size", fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/S00_shore_rec_catchability_surface.png", plot = last_plot())


## 5. Allocating effort to cells ----------------------------------------------

# We know the fishing effort over the years
# We know the fishing effort across months
# We know how 'useful' each cell is for fishing (utility function)
# Now we want to allocate fishing effort to cells over time
CP_trips <- readRDS(file_shore_CP_trips)
  
CP_trips2 <- CP_trips %>%
  dplyr::select(num_year, month, carpark, adjusted_effort) %>%
  pivot_wider(
    names_from = carpark,
    values_from = adjusted_effort,
    values_fill = 0  # Fill missing effort with 0s
  ) %>%
  arrange(num_year, month)

ggplot(shore_effort_df %>% filter(CP_index %in% c(1:30)), aes(x=year, y=adjusted_effort)) +
  geom_line(color = colour_palette[5]) +
  facet_wrap(~carpark, ncol = 5)
ggsave("plots/checking_plots_during_setup/S00_shore_rec_effort_over_time_by_carpark.png", plot = last_plot())

# link the carparks to their cells
water_dist <- water_shore %>% left_join(DistCP, by = "ID")
water_dist_long <- water_dist %>% pivot_longer(cols = CP$area, names_to = "access", values_to = "Distance_km")

water_dist_long <- water_dist_long[water_dist_long$type == 'shore',]

cell_access_pairs <- water_dist_long %>%
  filter(!is.na(Distance_km)) %>%
  dplyr::select(ID, access) %>% 
  st_drop_geometry() %>% 
  glimpse()


# pre

# Unique shore cell IDs
cell_ids <- sort(unique(water$ID[water$ID %in% shore_cells]))
NCELL <- length(cell_ids)

# Map actual cell ID to array index
rownames(CP_U_pre18) <- water_shore$ID
cp_u_cell_ids <- rownames(CP_U_pre18)
cell_id_to_index <- setNames(seq_along(cp_u_cell_ids), cp_u_cell_ids)

s_fishing_pre18 <- array(0, dim = c(NCELL, 12, n_years_tot))
for (YEAR in 1:n_years_tot) {
  print(paste("Processing year", YEAR))
  
  for (MONTH in 1:12) {
    
    # for each carpark
    for (ACCESS in 1:length(CP$area)) {
      
      carpark_name <- CP$area[ACCESS]
      
      # fet the effort for this carpark in this year/month
      effort <- CP_trips %>%
        filter(num_year == YEAR, month == MONTH, carpark == carpark_name) %>%
        pull(adjusted_effort)
      
      if (length(effort) == 0 || is.na(effort)) next
      
      # find cells linked to this carpark
      linked_cells <- cell_access_pairs %>%
        filter(access == carpark_name) %>%
        pull(ID)
      
      # assign effort to each linked cell
      for (cell_id in linked_cells) {
        
        cell_index <- cell_id_to_index[as.character(cell_id)]
        utility <- CP_U_pre18[as.character(cell_id), carpark_name]
        
        if (is.na(utility)) next
        
        s_fishing_pre18[cell_index, MONTH, YEAR] <- s_fishing_pre18[cell_index, MONTH, YEAR] + (utility * effort)
      }
    }
  }
} # this loop assigns fishing effort to cells according to their utility.

summary(s_fishing_pre18[,,1])

# plot check
year_idx <- 30  # choose the year to plot
current_year <- 1944 + sort(unique(CP_trips2$num_year))[year_idx]  # 1944 + year 1 (1945) = 1945
current_month <- 1 # choose the month to plot
built_CPs <- CP %>% filter(year_strt < current_year | (year_strt == current_year & build_mnth <= current_month));st_crs(built_CPs) <- 4326 # filter CPs that are already built
effort_vec <- s_fishing_pre18[, 1, year_idx] # extract that year&month's effort as a vector (one value per cell)
water$effort[water$ID %in% shore_cells] <- effort_vec # add that effort to the grid

ggplot(water) +
  geom_sf(aes(fill = log(effort)), col = "lightgray") +
  geom_sf(data = built_CPs, color = "red", size = 0.5) +  # only plot active carparks
  scale_fill_gradientn(colours = c("blue", "blue"), na.value = NA) +
  labs(title = paste("Shore fishing effort check - Year", current_year, "Month", current_month),
       fill = "log Effort") +
  theme_minimal()

s_fishing <- s_fishing_pre18 # only 1, not pre and post
dimnames(s_fishing)[[1]] <- shore_cells # now the rownames are the cells' ID

s_fishing[,,80]



## Set up effort for burn-in --------------------------------------------------

# shore
location_name <- unique(CP$area)
NCELL_shore <- nrow(water_shore)

CP_trips <- readRDS(file_shore_effort) %>%
  group_by(carpark) %>%
  summarise(boat_days = sum(adjusted_effort, na.rm = TRUE)) %>%
  arrange(desc(boat_days))
CP_trips$effort <- 1

CP_trips <- CP_trips %>%
  mutate(trip_per_hr = as.numeric(unlist((shore_days / effort)))) %>% # Standardise the no. trips based on how much time you spent sampling
  mutate(CP_prop = trip_per_hr/sum(trip_per_hr)) #Then work out the proportion of trips each hour that leave from each boat ramp


# Fishing parameters
eq.init.fish = 0.025
q = 0.00001
effort = (-log(1 - eq.init.fish)) / q # We assume the same level of nominal effort in each year

## Split up this effort by the same proportions as before and allocate it to the different access points
# Months
burn_in_effort <- prop_month_ave[, 2] * effort

burn_in_effort <- as.data.frame(burn_in_effort) %>%
  rename(effort = "ave_month_prop")
for (location in location_name) {
  burn_in_effort[[location]] <- 0
} # Loop through each location name and create a new column with 0's

for (M in 1:12) { # for each month,
  for (i in 1:nrow(CP_trips)) { # and each boat ramp,
    location_name <- CP_trips$carpark[i] # extract the location name
    effort_value <- CP_trips$effort[i]  # and effort value
    burn_in_effort[[location_name]] <- effort_value * CP_trips$CP_prop[i] # and calculate burn in effort
  }
} # this loop calculate burn-in effort for all carparks, based on how many visits and hours each carparks gets fished from


## Allocate to the cells using the same utilities that we set up earlier
s_burn_in_fishing <- array(0, dim = c(NCELL_shore, 12, n_years_tot)) #This array has a row for every cell, a column for every month, and a layer for every year
months <- array(0, dim = c(NCELL_shore, 12))
carpark <- array(0, dim = c(NCELL_shore, length(location_name)))
yr <- 1

for(YEAR in 1:n_years_tot){
  
  print(YEAR)
  
  for(MONTH in 1:12){
    
    for(ACCESS in 1:length(location_name)){
      
      temp <- burn_in_effort %>%
        dplyr::select(-c(effort))
      
      temp <- as.matrix(temp)
      
      for(CELL in 1:NCELL_shore){
        carpark[CELL, ACCESS] <- CP_U_pre18[CELL, ACCESS] * temp[MONTH, ACCESS] # Use the same utility from before2018 as nothing should have changed
      }
    }
    
    months[, MONTH] <-  rowSums(carpark)
  }
  s_burn_in_fishing[ , , yr] <- months
  s_burn_in_fishing[ , , yr] <- s_burn_in_fishing[ , , yr] * spatial_q[ACCESS, 1]
  yr <- yr + 1
} # this loop assigns each cell an amount of fishing effort each month and across years based on how useful it is (how far from a carpark and how frequented that carpark is)

# add back the offshore cells
dim(s_burn_in_fishing)
glimpse(s_burn_in_fishing[,,1])
glimpse(water)

shore_idx <- which(water$type == "shore")  # or some other logic identifying shore cells
length(shore_idx)  # should be 255 to match s_burn_in_fishing

s_burn_in_fishing_full <- array(0, dim = c(nrow(water), dim(s_fishing)[2], dim(s_burn_in_fishing)[3]))
s_burn_in_fishing_full[shore_idx, , ] <- s_burn_in_fishing
summary(s_burn_in_fishing_full)
dim(s_burn_in_fishing_full)

saveRDS(s_burn_in_fishing_full, file = "data/output_data/S00_shore_rec_burn_in_fishing.rds")


## 04.C BOAT RECREATIONAL FISHING ---------------------------------------------
## 0. Files used in this script -----------------------------------------------

file_wa         <- "data/output_data/01_wadandi_land.shp"
file_bathy      <- "data/input_data/SW_crop_AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
file_boat_ramps <- "data/input_data/wadandi_boat_ramps.shp"
file_water      <- "data/output_data/03_water.rds"
file_network    <- "data/output_data/03_network_shapefile.shp"

file_rec_fishable <- "data/output_data/04A_commercial_fishable_area_over_time.shp"
file_rec_cell_dist <- "data/output_data/04C_boat_rec_cell_dist.rds"
file_rec_BR_trips <- "data/output_data/04C_BR_trips.rds"
file_rec_boat_days <- "data/output_data/04C_boat_rec_total_boat_days.rds"
file_rec_ramp_effort <- "data/output_data/04C_boat_rec_ramp_effort.rds"

year_start <- 1945
year_end <- 2024
n_years_tot <- year_end - year_start +1
n_years_pre18 <- 2018 - year_start

crs_raster <- "+proj=longlat +datum=WGS84 +no_defs"
bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -33.2), crs = crs_raster)


## 1. BOAT Catchability -------------------------------------------------------

# this is exactly the same as commercial boat spatial q, no need to remake it.
spatial_q <- readRDS("data/output_data/S00_spatial_q_NTZ.rds")


## 5. Create a utility function -----------------------------------------------

## Now need to create a separate fishing surface for each month of each year based on distance to boat ramp, size of each cell,
## and multiply that by the effort in the cell to spatially allocate the effort across the area. Effort is also able to go more offshore over time.
## But we need to account for the fact that there will be sanctuary zones going in and the effort that would have gone in there will get allocated somewhere else
## Will then need to put the rows/columns back in as 0s 

# What we want to do is to distribute effort across month and year based on :
# 1. distance to boat ramp, 2. size of the cell, and 3. whether the cell is 'fishable' that year and 4. fishing effort
# After the SZ comes in the effort will also be redistributed.

fishable_long_sf <- st_read(file_rec_fishable); names(fishable_long_sf) <- c("ID", "year", "fishable_prop", "type", "sand", "reef", "seagrass", "status", "cell_area", "Area", "geometry")
Cell_Dist <- readRDS(file_rec_cell_dist)
BR_trips <- readRDS(file_rec_BR_trips)

years <- sort(unique(fishable_long_sf$year))
ramps <- unique(BR_trips$boat_ramp)
nramps <- length(ramps)

# Get full list of cell IDs from fishable_long_sf
all_cell_ids <- sort(unique(fishable_long_sf$ID))
ncells <- length(all_cell_ids)
nyears <- length(years)

BR_U_array <- array(0, dim = c(ncells, nramps, nyears),
                    dimnames = list(cell_id = all_cell_ids,
                                    ramp = ramps,
                                    year = as.character(years)))

for (y in seq_along(years)) {
  yr <- years[y]
  
  # Step 2.1: Get fishable surface for the year
  fish_yr <- fishable_long_sf %>%
    filter(year == yr) %>%
    dplyr::select(ID, fishable_prop)
  
  # Step 2.2: Join distances, areas, and fishable_prop
  Vj_df <- Cell_Dist %>%
    inner_join(fish_yr, by = c("cell_id" = "ID")) %>%
    arrange(cell_id)
  
  # Step 2.3: Calculate raw utility: exp(-distance) * area * fishable proportion
  # We'll loop over ramps
  U_mat <- matrix(0, nrow = nrow(Vj_df), ncol = nramps)
  
  for (r in seq_along(ramps)) {
    ramp_name <- ramps[r]
    U_mat[, r] <- exp(-Vj_df[[ramp_name]]) * Vj_df$Area * Vj_df$fishable_prop # exp(-Vj_df...) because otherwise high utility is given to areas far from ramps
  }
  
  # Step 2.4: Normalize utilities so they sum to 1 per ramp
  U_norm <- sweep(U_mat, 2, colSums(U_mat, na.rm = TRUE), "/")
  
  # Step 2.5: Store in the 3D array
  # Align to full cell list
  row_ids <- match(Vj_df$cell_id, all_cell_ids)
  BR_U_array[row_ids, , y] <- U_norm
}
plot(BR_U_array[,,1])
head(BR_U_array[,,70])
head(water)

# quick plot check
br_slice <- BR_U_array[,,70]
row_sums <- rowSums(br_slice) # Calculate row sums (sum of utilities across ramps for each cell)
water$utility_sum <- row_sums
ggplot(water) +
  geom_sf(aes(fill = utility_sum), color = NA) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", 
                       na.value = "grey80", name = "Sum of Utilities") +
  theme_minimal() +
  labs(title = "Sum of Utilities Across Ramps per Cell",
       subtitle = "Layer 70 of BR_U_array") +
  theme(legend.position = "right")



# plot check - utility over time
ramps
ramp_check <- "GB_Bunbury_Stirling_St"
utility_check <- as.data.frame(BR_U_array[, ramp_check, ])  # dimensions: cells × years
utility_check$cell_id <- as.integer(rownames(utility_check))
head(utility_check)
water_catch <- water %>% mutate(cell_id = row_number()) %>% left_join(utility_check, by = "cell_id")
water_catch_long <- water_catch %>% pivot_longer(cols = as.character(years), names_to = "year", values_to = "Catchability")
ggplot(water_catch_long) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  facet_wrap(~ year, ncol = 16) + 
  labs(title = paste0("commercial catchability of each cell, \nby distance from ", ramp_check, " ramp, cell size and tech-fishability"), fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/S00_rec_boat_ramp_catchability_over_time.png", plot = last_plot())

# Plot check - utility across space (notice the difference between popular/non-popular ramps)
utility_check2 <- as.data.frame(BR_U_array[,,80]) %>% mutate(cell_id = as.numeric(rownames(BR_U_array[,,80])))
BR <- st_read(file_boat_ramps); plot(BR$geometry)

water_catch <- water %>% mutate(cell_id = row_number()) %>% left_join(utility_check2, by = "cell_id")
water_catch_long <- water_catch %>% pivot_longer(cols = BR$name, names_to = "Ramp", values_to = "Catchability")
ggplot(water_catch_long) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  facet_wrap(~ Ramp, ncol = 6) + 
  labs(title = "Catchability of each cell in 2024, by boat ramp distance and cell size - \nnotice differences between popular/non popular ramps", fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/S00_rec_boat_catchability_surface_by_ramp.png", plot = last_plot())


## 6. Allocating effort to cells ----------------------------------------------

# We know the fishing effort over the years
# We know the fishing effort across months
# We know how 'useful' each cell is for fishing (utility function)
# Now we want to allocate fishing effort to cells over time

BR_trips <- BR_trips %>%
  dplyr::select(num_year, month, boat_ramp, adjusted_effort) %>%
  pivot_wider(
    names_from = boat_ramp,
    values_from = adjusted_effort,
    values_fill = 0  # Fill missing effort with 0s
  ) %>%
  arrange(num_year, month)


# now assign to cells
NCELL <- nrow(water)
b_fishing <- array(0, dim = c(NCELL, 12, length(years))) # this array has a row for every cell, a column for every month, and a layer for every year
months <- array(0, dim = c(NCELL, 12)) # array for the number of months
ramps <- array(0, dim = c(NCELL, length(BR$name))) # array for the number of ramps
layer <- 1

head(BR_trips)
head(BR)
head(BR_U_array[,,1])

for(YEAR in 1:length(years)){ # for all years,
  
  print(YEAR)
  
  for(MONTH in 1:12){ # run the loop for every month,
    
    for(RAMP in 1:length(BR$name)){ # and for every ramp,
      
      temp <- BR_trips %>% 
        filter(num_year == YEAR) %>% 
        dplyr::select(-c(num_year, month))
      
      temp <- as.matrix(temp)
      for(CELL in 1:NCELL){ # and assign fishing effort to every cell, based on its utility
        ramps[CELL,  RAMP] <- BR_U_array[CELL, RAMP, YEAR] * temp[MONTH, RAMP]
      }
    }
    months[, MONTH] <- rowSums(ramps)
  }
  b_fishing[ , , layer] <- months 
  layer <- layer + 1
} # this loop assigns each cell a fishing effort based on utility by month
b_fishing

# plot check
head(b_fishing[,,80])

water
year_idx <- 75  # first year
month_idx <- 1  # January

effort_vec <- b_fishing[, month_idx, year_idx] # Extract effort vector (one value per cell)
water$effort <- effort_vec # Add to water polygons
ggplot(water) +
  geom_sf(aes(fill = log(effort))) +
  scale_fill_viridis_c() +
  labs(title = paste("Fishing Effort - Year", year_idx, "Month", month_idx),
       fill = "log Effort") +
  theme_minimal()

## Set up effort for burn-in --------------------------------------------------

# boat 
location_name <- unique(BR$name)

BR_trips <- ramp_effort_df %>%
  group_by(boat_ramp) %>%
  summarise(boat_days = sum(adjusted_effort, na.rm = TRUE)) %>%
  arrange(desc(boat_days))
BR_trips$effort <- 1

BR_trips <- BR_trips %>%
  mutate(trip_per_hr = as.numeric(unlist((boat_days / effort)))) %>% # Standardise the no. trips based on how much time you spent sampling
  mutate(BR_prop = trip_per_hr/sum(trip_per_hr)) #Then work out the proportion of trips each hour that leave from each boat ramp

# Fishing parameters
eq.init.fish = 0.025
q = 0.00001
effort = (-log(1 - eq.init.fish)) / q # We assume the same level of nominal effort in each year

## Split up this effort by the same proportions as before and allocate it to the different access points
# Months
burn_in_effort <- prop_month_ave[, 2] * effort

burn_in_effort <- as.data.frame(burn_in_effort) %>%
  rename(effort = "ave_month_prop")
for (location in location_name) {
  burn_in_effort[[location]] <- 0
} # Loop through each location name and create a new column with 0's

for (M in 1:12) { # for each month,
  for (i in 1:nrow(BR_trips)) { # and each boat ramp,
    location_name <- BR_trips$boat_ramp[i] # extract the location name
    effort_value <- BR_trips$effort[i]  # and effort value
    burn_in_effort[[location_name]] <- effort_value * BR_trips$BR_prop[i] # and calculate burn in effort
  }
} # this loop calculate burn-in effort for all boat ramps, based on how many visits and hours each ramp gets fished from


## Allocate to the cells using the same utilities that we set up earlier
b_burn_in_fishing <- array(0, dim = c(NCELL, 12, n_years_tot)) #This array has a row for every cell, a column for every month, and a layer for every year
months <- array(0, dim = c(NCELL, 12))
ramps <- array(0, dim = c(NCELL, length(location_name)))
yr <- 1

for(YEAR in 1:n_years_tot){
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:length(location_name)){
      
      temp <- burn_in_effort %>%
        dplyr::select(-c(effort))
      
      temp <- as.matrix(temp)
      
      for(CELL in 1:NCELL){
        ramps[CELL, RAMP] <- BR_U_array[CELL, RAMP, YEAR] * temp[MONTH, RAMP]
      }
    }
    
    months[, MONTH] <-  rowSums(ramps)
  }
  b_burn_in_fishing[ , , yr] <- months
  b_burn_in_fishing[ , , yr] <- b_burn_in_fishing[ , , yr] * spatial_q[, 1]
  yr <- yr + 1
} # this loop assigns each cell an amount of fishing effort each month and across years based on how useful it is (how far from a boatramp and how frequented that boat ramp is)


saveRDS(b_burn_in_fishing, file = "data/output_data/S00_boat_rec_burn_in_fishing.rds")

### END ###