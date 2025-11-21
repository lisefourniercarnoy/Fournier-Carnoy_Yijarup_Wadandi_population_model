# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Various literature figures.
# Task:    Set up fishing effort (understand the split of fishing effort within a year, over space, and across years)
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    September 2025

# -----------------------------------------------------------------------------

# Notes: This script goes like this:
# 1. Calculate how 'catchable' each cell is each year (based on how big the cell is and whether it's NTZ or not)
# 2. Use literature info to hindcast fishing effort trends since 1945
# 3. Make fishability of cells (whether they are located within 2km of the access point)
# 4. Set up a utility function (how 'useful' each cell is, based on its catchability and distance to the access points)
# 5. Allocate the fishing effort (from step 2) to each cell by how 'useful' it is (step 4)

# -----------------------------------------------------------------------------

## Load libraries -------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(sf)
library(raster)
library(stringr)
library(forcats)
library(RColorBrewer)
library(geosphere)
library(abind)
library(sfnetworks)
library(purrr)
library(exactextractr)


## Custom plotting parameters -------------------------------------------------
colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))
source("custom_theme.R")

## 0. Files used in this script -----------------------------------------------

file_wa         <- "data/output_data/01_wadandi_land.shp"
file_ntz        <- "data/output_data/01_wadandi_NTZ.shp"
file_bathy      <- "data/input_data/wadandi_250m_bathy.tif"
file_carpark    <- "data/input_data/wadandi_carparks.shp"
file_water      <- "data/output_data/03_water.rds"
file_network    <- "data/output_data/03_network_shapefile.shp"

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

ntz_shore_cells <- water$ID[water$type == "shore" & water$status %in% c('NTZ_for_shore_and_boat')]
shore_cells <- water$ID[water$type == "shore"]

# We'll calculate the catchability of each cell by its area and whether it's no-take or not.
water <- water %>% 
  mutate(area = as.vector((water$cell_area)/1000000))

# Cells that are in NTZs will have catchability of 0 after designation (2018)
water_area <- water %>% # this is for all cells, pelagic and shore
  dplyr::select(ID, status, area, type) %>% 
  mutate(area_2000 = ifelse(!(ID %in% shore_cells), 0, area), # make offshore cells unfishable for pre-SZ years
         area_2018 = ifelse(!(ID %in% shore_cells) | ID %in% ntz_shore_cells, 0, area) # make offshore cells + SZ cells unfishable for post-SZ years
         ) %>% 
  dplyr::select(ID, area_2000, area_2018) %>%
  mutate(sum_2000 = sum(area_2000),
         sum_2018 = sum(area_2018)
         ) %>% 
  st_drop_geometry() %>% 
  mutate(q_2000 = area_2000/sum_2000,
         q_2018 = area_2018/sum_2018
         ) %>% 
  glimpse()

# create an array of catchability for each cell (rows) and each year (columns)
NCELL <- length(shore_cells)
spatial_q <- array(0.000006, dim = c(NCELL, n_years_tot)) # Why is the original catchability set to 0.000006 ?

# index of the q_2000 and q_2018 columns
ix_q2000 <- which(colnames(water_area) == "q_2000")
ix_q2018 <- which(colnames(water_area) == "q_2018")

for (ROW in 1:NCELL) {
  id <- shore_cells[ROW]
  i_area <- match(id, water_area$ID)
  
  q_2000 <- water_area[i_area, ix_q2000]
  q_2018 <- water_area[i_area, ix_q2018]
  
  # Pre-NTZ (up to 2017)
  for (COL in 1:n_years_pre18) {
    spatial_q[ROW, COL] <- spatial_q[ROW, COL] / q_2000
  }
  
  # Post-NTZ (from 2018 onward)
  if (id %in% ntz_shore_cells) {
    spatial_q[ROW, (n_years_pre18+1):n_years_tot] <- 0
  } else {
    for (COL in (n_years_pre18 + 1):n_years_tot) {
      spatial_q[ROW, COL] <- spatial_q[ROW, COL] / q_2018
    }
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
ggplot(water_long[water_long$Year %in% c(2000, 2020),]) +
  geom_sf(data = land, fill = "lightgray", col = NA) +
  geom_sf(data = ntz, fill = "#E4F2FF", col = NA) +
  geom_sf(aes(fill = log(Catchability + 1e-06)), color = NA) + # adding a small value to log(catchability) because log(0) = -Inf and doesn't display well
  scale_fill_gradientn(colours = colour_palette[4:6], na.value = NA) +
  facet_wrap(~ Year, ncol = 2) +
  labs(title = "Catchability over time (pre- and post- NTZ)", 
       fill = "log(Catchability)") +
  theme_minimal()
ggsave("plots/checking_plots_during_setup/04B_shore_rec_catchability_pre.post-NTZ.png", plot = last_plot())

# Save for future use.
saveRDS(spatial_q, file = paste0("data/output_data/04B_shore_rec_spatial_q_NTZ.rds"))


## 2. Fishing days ------------------------------------------------------------

# Total shore days
# Provide data from the literature (extrapolated with multipliers, see google sheets "YIJARUP - Fishing effort reconstruction")
years <- c(1987, 1989, 2000, 2020); shore_days <- c(79600, 84895, 107380, 133095)
shore_days_lit <- data.frame(years, shore_days)
plot(x = shore_days_lit$years, y = shore_days_lit$shore_days, type = "l", col = colour_palette[4], lwd = 4, main = "shore days in wadandi country from the literature")

# fit a logistic curve
start_vals <- list(a = 100000, b = 0.2, c = 2000)
logistic_model <- nls(
  shore_days ~ a / (1 + exp(-b * (years - c))),
  start = start_vals
)
predict_years <- year_start:year_end
predicted_annual_days <- predict(logistic_model, newdata = data.frame(years = predict_years))
years <- year_start:year_end

annual_effort_shore <- data.frame(year = years, predicted_shore_days = predicted_annual_days)

# check
plot(annual_effort_shore$year, annual_effort_shore$predicted_shore_days, type = "l", col = colour_palette[5], lwd = 4,
     main = "Shore Days in Wadandi Country \nfrom the Literature \n(with some extrapolation)", 
     xlab = "Year", ylab = "Boat Days")
lines(shore_days_lit$years, shore_days_lit$shore_days, col = colour_palette[4], lwd = 4)
legend("bottomright", legend = c("Observed", "Predicted (logistic)"), col = c(colour_palette[4], colour_palette[5]), lwd = 4)

# Bring in seasonal multipliers
seasonal_multipliers <- c( # see figure 21c in Ryan et al. 2022
  "01" = 0.14, "02" = 0.081, "03" = 0.097, "04" = 0.081,
  "05" = 0.033, "06" = 0.033, "07" = 0.033, "08" = 0.033,
  "09" = 0.033, "10" = 0.065, "11" = 0.11, "12" = 0.26
)
seasonal_multipliers <- seasonal_multipliers / sum(seasonal_multipliers) # standardise so it adds up to 1
barplot(seasonal_multipliers, col = colour_palette[6], main = "distribution of yearly \nboat fishing effort by month in % \n(deduced from Ryan et al. 2022, fig. 21c)")

# Add monthly distribution back to the timeseries
shore_effort <- expand.grid(
  year = predict_years,
  month = sprintf("%02d", 1:12)
) %>%
  arrange(year, month) %>%
  mutate(
    annual_shore_days = rep(predicted_annual_days, each = 12),
    monthly_effort = annual_shore_days * seasonal_multipliers[month]
  ) %>%
  dplyr::select(year, month, monthly_effort)

# check
ggplot(shore_effort, aes(x = as.Date(paste(year, month, "01", sep = "-")), y = monthly_effort)) +
  geom_line(color = colour_palette[4]) +
  labs(title = "Monthly shore Fishing Effort (1980–2024)",
       x = "Date", y = "Monthly shore fishing Days") +
  theme_minimal() +
  geom_smooth(color = colour_palette[5])

saveRDS(shore_effort, "data/output_data/04B_shore_rec_total_boat_days.rds")

# Proportion of each month's contribution to yearly boat days
shore_month_prop <- shore_effort %>% 
  group_by(year) %>% 
  mutate(year_sum = sum(monthly_effort)) %>%
  mutate(month_prop = monthly_effort/year_sum) %>% 
  dplyr::select(-year_sum)

shore_month_prop <- shore_month_prop %>% 
  group_by(month) %>% 
  mutate(ave_month_prop = mean(month_prop))
prop_month_ave <- shore_month_prop[1:12, c(2, 5)]

plot(prop_month_ave)
saveRDS(prop_month_ave, "data/output_data/04B_shore_prop_month_ave.rds") # charlotte's 'Average_Monthly_Effort"


CP <- st_read(file_carpark) %>% 
  st_transform(4283) %>%
  st_make_valid() %>%
  mutate(year_strt = as.numeric(year_start),
         year_strt = ifelse(is.na(year_strt), 1945, year_strt), # fill in missing dates with the simulation start year (here year_start is 1945, and year_strt is the build date)
         build_mnth = ifelse(is.na(month_strt), 1, month_strt), # assume Jan if unknown
         norm_popularity = 1 / length(unique(area))) %>% # popularity is set equally for now
  glimpse()
plot(water$geometry); plot(CP$geometry, col = colour_palette[5], pch = 16, cex = 1.5, add = TRUE)

# Distribute effort across all ramps, across time
shore_effort <- shore_effort %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) # make a column with year and month of carpark build

shore_effort_df <- expand.grid(
  CP_index = 1:nrow(CP),
  date = shore_effort$date
) %>%
  mutate(
    build_date = as.Date(paste(CP$year_strt[CP_index], CP$build_mnth[CP_index], "01", sep = "-")),
    norm_popularity = CP$norm_popularity[CP_index],
    carpark_name = CP$name[CP_index]
  ) %>%
  filter(date >= build_date) %>%
  mutate(
    months_since_build = interval(build_date, date) %/% months(1),
    logistic_growth = 1 / (1 + exp(-0.1 * (months_since_build - 60))),
    carpark_weight = norm_popularity * logistic_growth,
    carpark = CP$area[CP_index]
  )

# Merge in monthly effort
shore_effort_df <- shore_effort_df %>%
  left_join(shore_effort, by = "date") %>%
  group_by(date) %>%
  mutate(
    total_weight = sum(carpark_weight),
    adjusted_effort = ifelse(total_weight > 0, monthly_effort * (carpark_weight / total_weight), 0),
    year = year(date),
    month = month(date)
  ) %>%
  ungroup() %>%
  dplyr::select(carpark, CP_index, year, month, adjusted_effort)
saveRDS(shore_effort_df, "data/output_data/04B_shore_rec_effort.rds")

# check each boat ramp is populated correctly
ggplot(shore_effort_df %>% filter(CP_index %in% c(1:9)), aes(x = year, y = adjusted_effort)) +
  geom_line(color = colour_palette[5]) +
  facet_wrap(~ carpark, ncol = 3)

# check that the sum of distributed effort is the same as the whole region's effort that we predicted earlier
test_region_predicted_effort <- shore_effort %>% # the predicted boat effort, from bits of the literature
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  dplyr::select(date, monthly_effort) %>%
  rename(predicted_total = monthly_effort)
test_all_shore_effort <- shore_effort_df %>% # each ramp's effort, summed back up (SHOULD BE EXACTLY LIKE THE PREDICTED EFFORT)
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  group_by(date) %>%
  summarise(distributed_total = sum(adjusted_effort), .groups = "drop")
ggplot() +
  geom_line(data = test_region_predicted_effort, aes(x = date, y = predicted_total), color = colour_palette[4], lwd = 2) +
  geom_line(data = test_all_shore_effort, aes(x = date, y = distributed_total), color = colour_palette[5], linetype = "solid") +
  labs(
    title = "Predicted effort and effort should overlap.\nif not, the effort splitting by ramp (norm_popularity) is not working",
    y = "Monthly Boat Days", x = "Date"
  ) +
  theme_minimal()


## 3. Access point distance to cells ------------------------------------------

st_centroid_within_poly <- function (poly) { #returns true centroid if inside polygon otherwise makes a centroid inside the polygon
  
  # check if centroid is in polygon
  centroid <- poly %>% st_centroid() 
  in_poly <- st_within(centroid, poly, sparse = F)[[1]] 
  
  # if it is, return that centroid
  if (in_poly) return(centroid) 
  
  # if not, calculate a point on the surface and return that
  centroid_in_poly <- st_point_on_surface(poly) 
  return(centroid_in_poly)
} # this function finds the centre of each cell, to calculate distance to other cells

network <- st_read(file_network)

## Work out the probability of visiting a cell from each access point based on distance and size
CP <- st_as_sf(CP)
st_crs(CP) <- NA 

water_shore <- water[water$type == "shore",]

centroids <- st_centroid_within_poly(water_shore)
plot(water_shore$geometry)

points <- as.data.frame(st_coordinates(centroids))%>% # The points start at the bottom left and then work their way their way right
  mutate(ID = row_number()) 
points_sf <- st_as_sf(points, coords = c("X", "Y")) 
st_crs(points_sf) <- NA 

network <- as_sfnetwork(network, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

net <- activate(network, "nodes")
st_crs(net)
st_crs(net) <- NA 

network_matrix <- st_network_cost(net, from = CP, to = points_sf)
network_matrix <- network_matrix * 111
dim(network_matrix)

glimpse(network_matrix)

DistCP <- as.data.frame(t(network_matrix))
colnames(DistCP) <- CP$area
head(DistCP) # this gives us each cell's distance to the carparks

# make cells beyond 2km distance not fishable
head(DistCP)
DistCP <- DistCP %>%  
  mutate(across(where(is.numeric), ~ifelse(. > 2, NA, .)))


# Plot check
DistCP <- DistCP %>% mutate(ID = water_shore$ID)
saveRDS(DistCP, "data/output_data/04B_shore_cell_dist.rds")

water_dist <- water_shore %>% left_join(DistCP, by = "ID")
water_dist_long <- water_dist %>% pivot_longer(cols = CP$area, names_to = "access", values_to = "Distance_km")

water_dist_long <- water_dist_long[water_dist_long$type == 'shore',]

ggplot(water_dist_long %>% filter(access %in% unique(water_dist_long$access)[c(1, 30)])) + 
  geom_sf(aes(fill = Distance_km), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  facet_wrap(~ access, ncol = 3) + 
  labs(title = "Cell distance to each shore access point", fill = "Distance (km)")
ggsave("plots/checking_plots_during_setup/04B_shore_access_distance.png", plot = last_plot())


## 4. Set up a utility function -----------------------------------------------
Cell_Vars <- DistCP %>% 
  mutate(Area = as.vector((water_shore$cell_area)/1000000)) # Cells are now in km^2 but with no units

## Now need to create a separate fishing surface for each month of each year based on distance to access point, size of each
## cell and multiply that by the effort in the cell to spatially allocate the effort across the area
## But we need to account for the fact that there will be sanctuary zones going in and the effort that would have gone in there will get allocated somewhere else
## Will then need to put the rows/columns back in as 0s 

# What we want to do is to distribute effort across month and year based on :
# 1. distance to access point, 2. size of the cell, and 3. fishing effort.
# After the SZ comes in the effort will also be redistributed.

NCELL_pre18 <- nrow(water_shore) # number of cells you can fish in (before SZ) - only shore cells
NCELL_post18 <- NCELL_pre18 - nrow(water_shore[water_shore$status == 'NTZ_for_shore_and_boat',]) # number of cells you can fish in (after SZ)

Vj <- Cell_Vars %>%
  mutate(vj = rowSums(across(all_of(CP$area)), na.rm = TRUE)) %>% 
  glimpse()

Vj_pre18 <- Vj # %>% mutate(across(where(is.numeric), ~ifelse(is.na(.), 0, .))) # replace NAs with 0
# catchable cells before SZ

Vj_post18 <- Vj %>% filter(!ID %in% ntz_shore_cells) # %>% mutate(across(where(is.numeric), ~ifelse(is.na(.), 0, .))) # replace NAs with 0 otherwise the line below returns NA
# catchable cells after SZ

area_col <- grep("Area", colnames(Vj)) # extract the column index for area - important for calculations. in theory it should be the same pre- and post-NTZ.


# pre-SZ

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
ggsave("plots/checking_plots_during_setup/04B_shore_rec_catchability_surface_before2018.png", plot = last_plot())

# post-SZ
dim(Vj_post18)
CP_U_post18 <- as.data.frame(matrix(0, nrow = NCELL_post18, ncol = length(unique(CP$area)))) #Set up data frame to hold utilities of cells
colnames(CP_U_post18) <- c(CP$area)

cellU <- matrix(NA, ncol = length(unique(CP$area)), nrow = NCELL_post18)

for(ACCESS in 1:length(unique(CP$area))){
  for (CELL in 1:NCELL_post18){
      U <- exp(Vj_post18[CELL, ACCESS] + log(Vj_post18[CELL, area_col]))
      U <- ifelse(is.na(U), 0, U)
      cellU[CELL, ACCESS] <- U
  }
} # this loop goes over each cell for each ramp after the NTZ is in place, and calculates how catchable each cell is based on how close it is to a carpark

rowU <- as.data.frame(colSums(cellU))

for (ACCESS in 1:length(unique(CP$area))){
  for (CELL in 1:NCELL_post18){
      CP_U_post18[CELL, ACCESS] <- (exp(Vj_post18[CELL, ACCESS] + log(Vj_post18[CELL, area_col])))/rowU[ACCESS, 1]
    }
  } # this loop goes over each cell for each ramp, and calculates how catchable (in %) each cell is based on how close it is to a access point and how popular that boat ramp is.
colSums(CP_U_post18, na.rm = T) # all adds up to 1, perfect.
head(CP_U_post18)


# Plot check 
non_ntz_cells <- setdiff(shore_cells, ntz_shore_cells)
CP_U_post18 <- bind_cols(ID = non_ntz_cells, CP_U_post18)
water_catch <- water_shore[water_shore$ID %in% non_ntz_cells,] %>% left_join(CP_U_post18, by = "ID")
water_catch_long <- water_catch %>% pivot_longer(cols = CP$area, names_to = "access_point", values_to = "Catchability")
ggplot(water_catch_long %>% filter(access_point %in% unique(water_catch_long$access_point)[c(51:52)])) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  facet_wrap(~ access_point, ncol = 16) + 
  labs(title = "(after 2018) catchability of each cell, by shore access distance and cell size", fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/04B_shore_rec_catchability_surface_after2018.png", plot = last_plot())


## 5. Allocating effort to cells ----------------------------------------------

# We know the fishing effort over the years
# We know the fishing effort across months
# We know how 'useful' each cell is for fishing (utility function)
# Now we want to allocate fishing effort to cells over time
CP_trips <- shore_effort_df %>% # This is just the trips from each carpark
  arrange(year, month) %>% 
  mutate(num_year = match(year, sort(unique(year)))) %>% # This is to number the years 1 to 80 for the loop.
  glimpse()
saveRDS(CP_trips, "data/output_data/04B_CP_trips.rds")

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
ggsave("plots/checking_plots_during_setup/04B_shore_rec_effort_over_time_by_carpark.png", plot = last_plot())

# link the carparks to their cells
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

s_fishing_pre18 <- array(0, dim = c(NCELL, 12, n_years_pre18))
for (YEAR in 1:n_years_pre18) {
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


# post
# Unique shore cell IDs
ntz_shore_cells
cell_ids <- sort(unique(water$ID[water$ID %in% shore_cells & !water$ID %in% ntz_shore_cells ]))
NCELL <- length(cell_ids)

# Map actual cell ID to array index
rownames(CP_U_post18) <- water_shore$ID[water_shore$status != c("NTZ_for_shore_and_boat")]
cp_u_cell_ids <- rownames(CP_U_post18)
cell_id_to_index <- setNames(seq_along(cp_u_cell_ids), cp_u_cell_ids)
s_fishing_post18 <- array(0, dim = c(NCELL, 12, n_years_tot - n_years_pre18))
layer <- 1

for (year_index in 1:(n_years_tot - n_years_pre18)) {
  YEAR <- n_years_pre18 + year_index
  print(paste("Processing year", YEAR))
  
  for (MONTH in 1:12) {
    
    # for each carpark
    for (ACCESS in 1:length(CP$area)) {
      
      carpark_name <- CP$area[ACCESS]
      
      # get the effort for this carpark in this year/month
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
        utility <- CP_U_post18[as.character(cell_id), carpark_name]
        
        if (is.na(utility)) next
        
        s_fishing_post18[cell_index, MONTH, year_index] <- s_fishing_post18[cell_index, MONTH, year_index] + (utility * effort)
      }
    }
  }
}

glimpse(s_fishing_post18[,,1])

# add NTZ cells back in
shore_cells
ntz_shore_cells
# Initialize full array with all shore cells (including NTZ)
s_fishing_post18_full <- array(0, dim = c(length(shore_cells), 12, n_years_tot - n_years_pre18))

# Create a mapping from cell ID to its index in shore_cells
shore_index_map <- setNames(seq_along(shore_cells), shore_cells)

# Create a mapping from the cell ID to its index in the reduced (non-NTZ) array
non_ntz_cells <- setdiff(shore_cells, ntz_shore_cells)
non_ntz_index_map <- setNames(seq_along(non_ntz_cells), non_ntz_cells)

# Copy existing values into the correct positions in the full array
for (i in seq_along(non_ntz_cells)) {
  cell_id <- non_ntz_cells[i]
  full_index <- shore_index_map[as.character(cell_id)]
  reduced_index <- non_ntz_index_map[as.character(cell_id)]
  
  s_fishing_post18_full[full_index, , ] <- s_fishing_post18[reduced_index, , ]
}


# plot check
year_idx <- 1  # choose the year to plot
current_year <- 2018 + sort(unique(CP_trips2$num_year))[year_idx]  # 2018 + year 1 (2019) = 2019
current_month <- 1 # choose the month to plot
built_CPs <- CP %>% filter(year_strt < current_year | (year_strt == current_year & build_mnth <= current_month));st_crs(built_CPs) <- 4326 # filter CPs that are already built
effort_vec <- s_fishing_post18_full[, 1, year_idx] # extract that year&month's effort as a vector (one value per cell)
water$effort[water$ID %in% shore_cells] <- effort_vec # add that effort to the grid

ggplot(water) +
  geom_sf(aes(fill = log(effort)), col = "lightgray") +
  geom_sf(data = built_CPs, color = "red", size = 0.5) +  # only plot active carparks
  scale_fill_gradientn(colours = c("blue", "blue"), na.value = NA) +
  labs(title = paste("Shore fishing effort check - Year", current_year, "Month", current_month),
       fill = "log Effort") +
  theme_minimal()


# merge pre- and post- effort distributions.
glimpse(s_fishing_pre18)
glimpse(s_fishing_post18_full)

s_fishing <- abind(s_fishing_pre18, s_fishing_post18_full, along = 3)
dimnames(s_fishing)[[1]] <- shore_cells # now the rownames are the cells' ID

s_fishing[,,80]


## X. Make a GIF --------------------------------------------------------------

library(gifski)

frame_count <- 1
m <- 12 # month to plot - using december because it's the highest

# make the plot's colour limits otherwise log(0) wouldnt work in the plot
log_effort_all <- log(s_fishing[s_fishing > 0]) # to avoid -Inf as a limit in the plot which wouldn't work
global_limits <- range(log_effort_all, na.rm = TRUE)

for (y in seq_along(years)) {
  
  # Update water$effort for this year and month
  #water$effort[water$type == "shore"] <- s_fishing[, 1, y]
  effort_this_year <- s_fishing[, m, y]
  names(effort_this_year) <- dimnames(s_fishing)[[1]]  # make sure these are character IDs
  water$effort[water$ID %in% shore_cells] <- effort_this_year
  
  # Year and month labels
  year_idx <- years[y]
  month_idx <- m

  # Plot
  p <- ggplot() +
    geom_sf(data = land, fill = "lightgray", col = NA) +
    geom_sf(data = ntz, fill = "#E4F2FF", col = NA) +
    geom_sf(data = water, aes(fill = log(effort)), color = NA) +
    scale_fill_gradientn(colours = colour_palette[6:4], na.value = NA, 
                         limits = global_limits) +
    labs(
      title = paste("Fishing Effort - Year", year_idx, "Month", month_idx),
      fill = "log Effort"
    ) +
    theme_minimal()
  
  # Save frame
  ggsave(
    filename = sprintf("plots/gif_frames/shore_effort_frame_%03d.png", frame_count),
    plot = p,
    width = 6, height = 6, dpi = 150
  )
  
  frame_count <- frame_count + 1
}

# stitch the GIF frames together
png_files <- list.files("plots/gif_frames", pattern = "shore_effort_frame_\\d+\\.png", full.names = TRUE)
gifski(
  png_files,
  gif_file = "plots/gifs/shore_fishing_effort_over_time.gif",
  width = 600,
  height = 600,
  delay = 0.25  # seconds per frame (adjust as needed)
)

## Set up effort for burn-in --------------------------------------------------

# shore
location_name <- unique(CP$area)
NCELL_shore <- nrow(water_shore)

CP_trips <- shore_effort_df %>%
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

# Split up this effort by the same proportions as before and allocate it to the different access points
burn_in_effort <- prop_month_ave[, 2] * effort

burn_in_effort <- as.data.frame(burn_in_effort) %>%
  rename(effort = "ave_month_prop")
for (location in location_name) {
  burn_in_effort[[location]] <- 0
} # Loop through each location name and create a new column with 0's

for (M in 1:12) { # for each month,
  for (i in 1:nrow(CP_trips)) { # and each boat ramp,
    location_name <- CP_trips$carpark[i] # extract the location name
    burn_in_effort[[location_name]] <- effort * CP_trips$CP_prop[i] # and calculate burn in effort
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

saveRDS(s_burn_in_fishing_full, file = "data/output_data/04B_shore_rec_burn_in_fishing.rds")

# Charlotte then adds another burn-in array for high mortality, I'm not doing rn for the sake of moving along

## END ##