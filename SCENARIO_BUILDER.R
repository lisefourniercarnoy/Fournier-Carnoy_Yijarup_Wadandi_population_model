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

library(tidyverse) # data manipulation
library(sf) # shapefile manipulation
library(raster) # raster manipulation
library(RColorBrewer) # plotting colours
library(abind) # matrix manipulation
library(sfnetworks) # for making distance to access points
library(exactextractr) # estracting raster values at points

# The point of this script is to set up a fishing effort surface over time based
# on user-defined spatial, temporal and fleet restrictions.

## SCENARIO DEFINITION --------------------------------------------------------

scenario_name <- "sXX_test"

shore_fishing <- data.frame(
  fleet_name = "shore_rec",
  
  # is this fleet allowed at all?
  fleet_allowed = "yes", 
  
  # spatial restrictions
  current_ntz = "yes", # are current NTZs protected?
  new_ntz_polygon_file = "data/input_data/new_ntz_scenario_polygons/test_NTZ.shp", # what's the spatial restriction you want?
  
  # when is the spatial restriction deviating from reality?
  onset_of_spatial_management_m = "01", 
  onset_of_spatial_management_y = "2000", 
  
  # temporal restrictions - what % of the month is open to fishing?
  jan = 1,
  feb = 1,
  mar = 1,
  apr = 1,
  may = 1,
  jun = 1,
  jul = 1,
  aug = 1,
  sep = 1,
  oct = 1,
  nov = 1,
  dec = 1,
  
  # when is the temporal restriction deviating from reality?
  onset_of_temporal_management_m = "01", 
  onset_of_temporal_management_y = "2000"
  
  # catch restrictions (idk how to do yet)
)

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
file_carpark_north    <- "data/input_data/north_carparks.shp"

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



