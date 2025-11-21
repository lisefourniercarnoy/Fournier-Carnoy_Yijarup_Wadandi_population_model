# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Various literature figures.
# Task:    Set up fishing effort (understand the split of fishing effort within a year, over space, and across years)
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    May 2025

# -----------------------------------------------------------------------------

# Notes: This script goes like this:
# 1. Calculate how 'catchable' each cell is each year (based on how big the cell is and whether it's NTZ or not)
# 2. Use literature info to hindcast fishing effort trends since 1945
# 3. Make fishability of cells (whether they can be fished based on the tech at that time) increase over time
# 4. Calculate the distance from access points to each cell (do this separately for shore and boat fishing)
# 5. Set up a utility function (how 'useful' each cell is, based on its catchability and distance to the access points)
# 6. Allocate the fishing effort (from step 2) to each cell by how 'useful' it is (step 4)

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
file_bathy      <- "data/input_data/SW_crop_AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
file_boat_ramps <- "data/input_data/wadandi_boat_ramps.shp"
file_water      <- "data/output_data/03_water.rds"
file_network    <- "data/output_data/03_network_shapefile.shp"

year_start <- 1945
year_end <- 2024
n_years_tot <- year_end - year_start +1
n_years_pre18 <- 2018 - year_start

crs_raster <- "+proj=longlat +datum=WGS84 +no_defs"
bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -33.2), crs = crs_raster)


## 1. Catchability ------------------------------------------------------------

water <- readRDS(file_water) %>% filter(!is.na(ID)); plot(water)
no_take_list <- list(water$ID[water$status %in% c('NTZ_boat_shore', 'NTZ_boat')]) # List of all the cells that are No-Take for boats

# We'll calculate the catchability of each cell by its area and whether it's no-take or not.
water <- water %>% 
  mutate(Area = as.vector((water$cell_area)/1000000))

# Cells that are in NTZs will have catchability of 0 after designation (2018)
water_area <- water %>% 
  dplyr::select(status, Area) %>% 
  mutate(area_2000 = Area,
         area_2018 = ifelse(status %in% c('NTZ_boat_shore', 'NTZ_boat'), 0, Area)) %>%
  dplyr::select(area_2000, area_2018) %>% 
  mutate(sum_2000 = sum(area_2000),
         sum_2018 = sum(area_2018)) %>% 
  st_drop_geometry() %>% 
  mutate(q_2000 = area_2000/sum_2000,
         q_2018 = area_2018/sum_2018) %>% 
  mutate(ID = row_number())

# Create an array of catchability for each cell (rows) and each year (columns)
NCELL <- nrow(water)
spatial_q <- array(0.000006, dim = c(NCELL, n_years_tot)) # Why is the original catchability set to 0.000006 ?

# The following bit is for changing catchability after NTZ implementation.
for (COL in (n_years_pre18+1):n_years_tot) { # 2018 is year number 38
  spatial_q[, COL] <- spatial_q[, COL-1] * 1.02
} # This increases post-2018 cells by 2% each year?? if the min size increases surely catchability decreases?

# Change catchability of cells pre-SZ
for (COL in 1:n_years_pre18) { # 2017 is year number 37, and the last year with fishing before NTZs
  for (ROW in 1:NCELL){
    spatial_q[ROW, COL] <- spatial_q[ROW, COL] / water_area[ROW, 5]
  }
} # For all cells pre-SZ, catchability is in proportion to the cell's area.

for (COL in (n_years_pre18+1):n_years_tot) { # +1 because n_years_pre18 relates to 2017
  for (ROW in 1:NCELL){
    spatial_q[ROW, COL] <- spatial_q[ROW, COL] / water_area[ROW, 6]
  }
} # For all cells post-SZ, catchability is in proportion to the cell's area. (Cells in NTZs are obviously zero)

spatial_q[spatial_q == Inf] <- 0 # IDK what this does.
summary(spatial_q) # this is a matrix that tells us for each cell (each row), and each year (each column), how likely you'd catch a fish in that cell based on how big it is.

# Plot check
catch_df <- as.data.frame(spatial_q)
colnames(catch_df) <- paste0("Year_", year_start:(year_start + n_years_tot-1))
catch_df$ID <- water_area$ID  # or just 1:NCELL if they match in order
water_catch <- water %>% left_join(catch_df, by = "ID")  # ID must be in 'water' too
water$ID <- water_area$ID  # or use `row_number()`
water_long <- water_catch %>% pivot_longer(cols = starts_with("Year_"), names_to = "Year", names_prefix = "Year_", values_to = "Catchability") %>% mutate(Year = as.numeric(Year))
ggplot(water_long[water_long$Year %in% c(2000, 2020),]) +
  geom_sf(aes(fill = log(Catchability)), color = NA) +
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  facet_wrap(~ Year, ncol = 2) +
  labs(title = "Catchability over time (pre- and post- NTZ)", fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/04A_commercial_catchability_pre.post-NTZ.png", plot = last_plot())

# Save for future use.
saveRDS(spatial_q, file = paste0("data/output_data/04A_commercial_spatial_q_NTZ.rds"))


## 2. Fishing days -------------------------------------------------------------

# obtain boat days from the literature
# Provide data from the literature (obtained **roughly** from the Fisheries Research Report 2007, p.32)
years <- c(1975, 1980, 1985, 1990, 1995, 2000, 2005); boat_days <- c(3100, 3000, 4500, 2750, 2000, 2000, 2250)
boat_days_lit <- data.frame(years, boat_days)

# adding fake boat days to fill out. can't make fit a curve because they're all terrible. based on no literature, just vibes.
years <- c(1945, 1960, 2020)
boat_days <- c(1000, 2000, 2000)
boat_days_fake <- data.frame(years, boat_days)

plot(boat_days_lit$years, boat_days_lit$boat_days, pch = 19, col = colour_palette[4], xlim = c(1940, 2025), ylim = c(0, 5000), xlab = "Year", ylab = "Boat Days", main = "Original + Fake Boat Days")
points(boat_days_fake$years, boat_days_fake$boat_days, pch = 17, col = colour_palette[6])
legend("topright", legend = c("Original Data", "Fake Data"), col = c(colour_palette[4], colour_palette[6]), pch = c(19, 17))

# fill in the gaps to obtain values for every year
years <- c(boat_days_lit$years, boat_days_fake$years)
boat_days <- c(boat_days_lit$boat_days, boat_days_fake$boat_days)

# Sort data (important for interpolation)
sorted_index <- order(years)
all_years_sorted <- years[sorted_index]
all_boat_days_sorted <- boat_days[sorted_index]

years_full <- 1945:2024
# linear interpolation with extrapolation
interp <- approx(x = all_years_sorted, y = all_boat_days_sorted, xout = years_full, method = "linear", rule = 2)
annual_effort_df <- data.frame(year = interp$x, boat_days = interp$y)

# check
plot(annual_effort_df$year, annual_effort_df$boat_days, type = "l", col = colour_palette[5], lwd = 4,
     main = "Boat Days in Wadandi Country \nfrom the Literature \n(with some extrapolation)", 
     xlab = "Year", ylab = "Boat Days")
lines(boat_days_lit$years, boat_days_lit$boat_days, col = colour_palette[4], lwd = 4)
legend("topleft", legend = c("Observed", "Predicted"), col = c(colour_palette[4], colour_palette[5]), lwd = 4)

# Bring in seasonal multipliers
seasonal_multipliers <- c( # see figure 21c in Ryan et al. 2022
  "01" = 0.14, "02" = 0.081, "03" = 0.097, "04" = 0.081,
  "05" = 0.033, "06" = 0.033, "07" = 0.033, "08" = 0.033,
  "09" = 0.033, "10" = 0.065, "11" = 0.11, "12" = 0.26
)
seasonal_multipliers <- seasonal_multipliers / sum(seasonal_multipliers) # standardise so it adds up to 1
barplot(seasonal_multipliers, col = colour_palette[6], main = "distribution of yearly \nboat fishing effort by month in % \n(deduced from Ryan et al. 2022, fig. 21c)")

# Add monthly distribution back to the timeseries
boat_effort <- expand.grid(
  year = years_full,
  month = sprintf("%02d", 1:12)
) %>%
  arrange(year, month) %>%
  mutate(
    annual_boat_days = rep(annual_effort_df$boat_days, each = 12),
    monthly_effort = annual_boat_days * seasonal_multipliers[month]
  ) %>%
  dplyr::select(year, month, monthly_effort)

# check
ggplot(boat_effort, aes(x = as.Date(paste(year, month, "01", sep = "-")), y = monthly_effort)) +
  geom_line(color = colour_palette[4]) +
  labs(title = "Monthly Boat Fishing Effort (1945–2024)",
       x = "Date", y = "Monthly Boat Days") +
  theme_minimal() +
  geom_smooth(color = colour_palette[5])

saveRDS(boat_effort, "data/output_data/04A_commercial_total_boat_days.rds")


# Proportion of each month's contribution to yearly boat days
boat_month_prop <- boat_effort %>% 
  group_by(year) %>% 
  mutate(year_sum = sum(monthly_effort)) %>%
  mutate(month_prop = monthly_effort/year_sum) %>% 
  dplyr::select(-year_sum)

boat_month_prop <- boat_month_prop %>% 
  group_by(month) %>% 
  mutate(ave_month_prop = mean(month_prop))
prop_month_ave <- boat_month_prop[1:12, c(2, 5)]

plot(prop_month_ave)
saveRDS(prop_month_ave, "data/output_data/04A_commercial_prop_month_ave.rds") # charlotte's 'Average_Monthly_Effort"

# effort by boat_ramp
BR <- st_read(file_boat_ramps) %>% 
  st_transform(4283) %>%
  st_make_valid() %>%
  #filter(!is.na(visit_2024)) %>% # remove ramps where there's no info on visitation (TEMPORARY, all should have data at some point)
  mutate(build_year = as.numeric(build_year),
         build_year = ifelse(is.na(build_year), year_start, build_year), # fill in missing dates with the start year
         build_mnth = ifelse(is.na(build_mnth), 1, build_mnth),     # assume Jan if unknown
         norm_popularity = com_prop / sum(com_prop, na.rm = TRUE)) %>%
  glimpse()

# for commercial fishing, only a few boat ramps are used (see Andrea's historical fishing resources, fishing localities in ABS stats)
unique(BR$name)
BR <- BR %>% 
  dplyr::filter(BR$name %in% c("SC_Augusta_Ellis_St_Jetty", "WC_Gnarabup", "WC_Hamelin_Bay",
                               "GB_Quindalup", "GB_Eagle_Bay", "GB_Bunbury_Stirling_St", "GB_Busselton_Georgette_Street"))
plot(water$geometry); plot(BR$geometry, col = colour_palette[5], pch = 16, cex = 3, add = TRUE)


# Distribute effort across all ramps, across time
boat_effort <- boat_effort %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) # make a column with year and month of ramp build
ramp_effort_df <- expand.grid(
  ramp_index = 1:nrow(BR),
  date = boat_effort$date
) %>%
  mutate(
    build_date = as.Date(paste(BR$build_year[ramp_index], BR$build_mnth[ramp_index], "01", sep = "-")),
    norm_popularity = BR$norm_popularity[ramp_index],
    ramp_name = BR$name[ramp_index]
  ) %>%
  filter(date >= build_date) %>%
  mutate(
    months_since_build = interval(build_date, date) %/% months(1),
    logistic_growth = 1 / (1 + exp(-0.1 * (months_since_build - 60))),
    ramp_weight = norm_popularity * logistic_growth
  )

# Merge in monthly effort
ramp_effort_df <- ramp_effort_df %>%
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
saveRDS(ramp_effort_df, "data/output_data/04A_commercial_ramp_effort.rds")

# Check total per month equals monthly_effort
check_totals <- ramp_effort_df %>%
  group_by(year, month) %>%
  summarise(total_effort = sum(adjusted_effort), .groups = "drop") %>%
  left_join(boat_effort %>% mutate(year = year(date), month = month(date)), by = c("year", "month")) %>%
  mutate(diff = abs(total_effort - monthly_effort))
summary(check_totals$diff)  # should be near zero

# check each boat ramp is populated correctly
ggplot(ramp_effort_df, aes(x=year, y=adjusted_effort)) +
  geom_line(color = colour_palette[5], lwd = 1) +
  facet_wrap(~boat_ramp, ncol = 1)

# check that the sum of distributed effort is the same as the whole region's effort that we predicted earlier
test_region_predicted_effort <- boat_effort %>% # the predicted boat effort, from bits of the literature
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  dplyr::select(date, monthly_effort) %>%
  rename(predicted_total = monthly_effort)
test_all_ramp_effort <- ramp_effort_df %>% # each ramp's effort, summed back up (SHOULD BE EXACTLY LIKE THE PREDICTED EFFORT)
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  group_by(date) %>%
  summarise(distributed_total = sum(adjusted_effort), .groups = "drop")
ggplot() +
  geom_line(data = test_region_predicted_effort, aes(x = date, y = predicted_total), color = colour_palette[4], lwd = 2) +
  geom_line(data = test_all_ramp_effort, aes(x = date, y = distributed_total), color = colour_palette[5], linetype = "solid") +
  labs(
    title = "Predicted effort and ramps' effort should overlap.\nif not, the effort splitting by ramp (norm_popularity) is not working",
    y = "Monthly Boat Days", x = "Date"
  ) +
  theme_minimal()



## 3. "boats can fish further over time" --------------------------------------

# prepare grid cells
plot(water$geometry)
wa_mask <- st_read(file_wa); wa_mask <- st_transform(wa_mask, crs = crs_raster); wa_mask <- as(wa_mask, "Spatial"); plot(wa_mask, col = "lightgray", add = T)

# prepare the bathymetry layer
bathy <- raster(file_bathy) %>% 
  projectRaster(crs = crs_raster) %>% 
  crop(extent(bbox)) %>% 
  abs() %>% 
  mask(wa_mask, inverse = TRUE); plot(bathy)

# parameters
years <- year_start:year_end
min_depth <- 20  # fishable depth is about 20m in 1945 (from Gaynor 2008, p.38)
increase_year <- 1950


# classify raster cells by whether they're within the fishable depth for that year or not
results <- matrix(NA, nrow = ncell(bathy), ncol = length(years)) # Empty data frame to store results
for (i in seq_along(years)) {
  year <- years[i]
  threshold <- min_depth + 1 * (ifelse(year - increase_year < 0, 0, year - increase_year))  # depth threshold per year, increases by 1m per year
  
  fishable <- bathy <= threshold # fishable mask: 1 = fishable, 0 = too deep
  
  fishable_vals <- getValues(fishable)
  results[, i] <- fishable_vals
} # this loop checks for every raster cell and every year whether the cell is fishable (by depth)
results_df <- as.data.frame(results)
colnames(results_df) <- paste0("year_", years)

# Add cell index or coordinates to join with spatial data
results_df$cell <- 1:ncell(bathy)
coords <- xyFromCell(bathy, results_df$cell)
results_df <- cbind(results_df, coords)
head(results_df)

results_df[ , grepl("^year_", names(results_df))] <- 
  lapply(results_df[ , grepl("^year_", names(results_df))], as.numeric)

#plot check
long_df <- pivot_longer(
  results_df,
  cols = starts_with("year_"),
  names_to = "year",
  names_prefix = "year_",
  values_to = "fishable"
)
ggplot(long_df %>% filter(year %in% c(1945:2024)), aes(x = x, y = y, fill = fishable)) +
  geom_raster() +
  scale_fill_gradient(low = "white", high = colour_palette[5], name = "Fishable") +
  coord_fixed() +
  facet_wrap(~year, ncol = 10) +
  theme_minimal()

# overlap with the grid
bathy[is.na(bathy)] <- 1 # these are shore cells that are not covered by bathy, as they're all <20m, they are fished from 1945.
fishable_rasters <- list()

for (year in years) {
  threshold <- min_depth + ifelse(year - increase_year < 0, 0, year - increase_year)
  fishable_rasters[[as.character(year)]] <- bathy <= threshold
} # this loop makes a raster out of fishable/non fishable cells for every year

fishable_summary <- data.frame(ID = water$ID)
water <- st_make_valid(water) # make the object consistent, otherwise the computation doesn't work
water <- st_collection_extract(water, "POLYGON")
water <- st_cast(water, "MULTIPOLYGON")

for (year in years) {
  fishable_r <- fishable_rasters[[as.character(year)]]
  
  # Compute fraction of each polygon that is fishable
  summary_vals <- exact_extract(fishable_r, water, 'mean')
  
  # Add to summary table
  fishable_summary[[paste0("year_", year)]] <- summary_vals
} # this loop calculates the area of each cell that is fishable (i.e. within the depth limit of that year)
ntz_mask <- water$status == "NTZ"

# Apply override from 2019 onward
fishable_summary[ntz_mask, paste0("year_", 2019:2024)] <- 0

fishable_long <- fishable_summary %>%
  pivot_longer(
    cols = starts_with("year_"),
    names_to = "year",
    names_prefix = "year_",
    values_to = "fishable_prop"
  ) %>%
  mutate(year = as.integer(year)) %>%
  left_join(water, by = "ID")  # add geometry back
fishable_long_sf <- st_as_sf(fishable_long)

ggplot(fishable_long_sf %>% filter(year %in% c(1960:1975))) +
  geom_sf(aes(fill = fishable_prop)) +
  scale_fill_gradient(low = colour_palette[6], high = colour_palette[4], name = "Fishable") +
  facet_wrap(~year) +
  theme_minimal()

st_write(fishable_long_sf, "data/output_data/04A_commercial_fishable_area_over_time.shp", append = FALSE)


## 4. Access point distance to cells ------------------------------------------

BR <- st_read(file_boat_ramps) %>% 
  st_transform(4283) %>%
  st_make_valid() %>%
  #filter(!is.na(visit_2024)) %>% # remove ramps where there's no info on visitation (TEMPORARY, all should have data at some point)
  mutate(build_year = as.numeric(build_year),
         build_year = ifelse(is.na(build_year), year_start, build_year), # fill in missing dates with the start year
         build_mnth = ifelse(is.na(build_mnth), 1, build_mnth),     # assume Jan if unknown
         norm_popularity = com_prop / sum(com_prop, na.rm = TRUE)) %>% # com_prop is the proportion of catch coming from each 'boat_ramp' (locality) - see effort reconstruction excel sheet
  glimpse()

# for commercial fishing, only a few boat ramps are used (see Andrea's historical fishing resources, fishing localities in ABS stats)
unique(BR$name)
BR <- BR %>% 
  dplyr::filter(BR$name %in% c("SC_Augusta_Ellis_St_Jetty", "WC_Gnarabup", "WC_Hamelin_Bay",
                               "GB_Quindalup", "GB_Eagle_Bay", "GB_Bunbury_Stirling_St", "GB_Busselton_Georgette_Street"))

plot(water$geometry); plot(BR$geometry, col = colour_palette[6], pch = 16, add = TRUE)

st_centroid_within_poly <- function (poly) { # returns true centroid if inside polygon otherwise makes a centroid inside the polygon
  
  # check if centroid is in polygon
  centroid <- poly %>% st_centroid() 
  in_poly <- st_within(centroid, poly, sparse = F)[[1]] 
  
  # if it is, return that centroid
  if (in_poly) return(centroid) 
  
  # if not, calculate a point on the surface and return that
  centroid_in_poly <- st_point_on_surface(poly) 
  return(centroid_in_poly)
}

network <- st_read(file_network)

## Work out the probability of visiting a cell from each boat ramp based on distance and size
BR <- st_as_sf(BR)
st_crs(BR) <- NA 

centroids <- st_centroid_within_poly(water)
points <- as.data.frame(st_coordinates(centroids))%>% #The points start at the bottom left and then work their way their way right
  mutate(ID = row_number()) 
points_sf <- st_as_sf(points, coords = c("X", "Y")) 
st_crs(points_sf) <- NA 

network <- as_sfnetwork(network, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

net <- activate(network, "nodes")
st_crs(net)
st_crs(net) <- NA 

network_matrix <- st_network_cost(net, from = BR, to = points_sf)
network_matrix <- network_matrix * 111
dim(network_matrix)

glimpse(network_matrix)
DistBR <- as.data.frame(t(network_matrix))
colnames(DistBR) <- BR$name
head(DistBR) # this gives us each cell's distance to the boat ramps

# Plot check
DistBR <- DistBR %>% mutate(ID = 1:nrow(.))
water_dist <- water %>% left_join(DistBR, by = "ID")
water_dist_long <- water_dist %>% pivot_longer(cols = BR$name, names_to = "Ramp", values_to = "Distance_km")
ggplot(water_dist_long) + 
  geom_sf(aes(fill = Distance_km), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  facet_wrap(~ Ramp, ncol = 4) + 
  labs(title = "Cell distance to each boat ramp", fill = "Distance (km)")
ggsave("plots/checking_plots_during_setup/04A_commercial_boat_ramp_distance.png", plot = last_plot())


## 5. Create a utility function -----------------------------------------------

## Now need to create a separate fishing surface for each month of each year based on distance to boat ramp, size of each cell,
## and multiply that by the effort in the cell to spatially allocate the effort across the area. Effort is also able to go more offshore over time.
## But we need to account for the fact that there will be sanctuary zones going in and the effort that would have gone in there will get allocated somewhere else
## Will then need to put the rows/columns back in as 0s 

# What we want to do is to distribute effort across month and year based on :
# 1. distance to boat ramp, 2. size of the cell, and 3. whether the cell is 'fishable' that year and 4. fishing effort
# After the SZ comes in the effort will also be redistributed.


years <- sort(unique(fishable_long_sf$year))
ramps <- unique(BR$name)
nramps <- length(ramps)

# Get full list of cell IDs from fishable_long_sf
all_cell_ids <- sort(unique(fishable_long_sf$ID))
ncells <- length(all_cell_ids)
nyears <- length(years)

# Distance from each cell to each ramp (matrix)
# Assuming `DistBR` has same row order as `water`
Cell_Dist <- DistBR %>%
  mutate(cell_id = water$ID,
         Area = as.vector(water$cell_area / 1e6))  # km²
saveRDS(Cell_Dist, "data/output_data/04A_commercial_cell_dist.rds")

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
BR$name
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
ggsave("plots/checking_plots_during_setup/04A_commercial_boat_ramp_catchability_over_time.png", plot = last_plot())

# Plot check (notice the difference between popular/non-popular ramps)
utility_check2 <- as.data.frame(BR_U_array[,,80]) %>% mutate(cell_id = as.numeric(rownames(BR_U_array[,,80])))
class(utility_check2$cell_id)
water_catch <- water %>% mutate(cell_id = row_number()) %>% left_join(utility_check2, by = "cell_id")
water_catch_long <- water_catch %>% pivot_longer(cols = BR$name, names_to = "Ramp", values_to = "Catchability")
ggplot(water_catch_long) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  facet_wrap(~ Ramp, ncol = 4) + 
  labs(title = "Catchability of each cell in 2024, by boat ramp distance and cell size - \nnotice differences between popular/non popular ramps", fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/04A_commercial_catchability_surface_by_ramp.png", plot = last_plot())


## 6. Allocating effort to cells ----------------------------------------------

# We know the fishing effort over the years
# We know the fishing effort across months
# We know how 'useful' each cell is for fishing (utility function)
# Now we want to allocate fishing effort to cells over time
BR_trips <- ramp_effort_df %>% # This is just the trips from each boat ramp
  arrange(year, month) %>%
  mutate(num_year = match(year, sort(unique(year)))) # This is to number the years 1 to 80 for the loop.
saveRDS(BR_trips, "data/output_data/04A_BR_trips.rds")

ggplot(BR_trips, aes(x = order(year, month), y = adjusted_effort)) +
  geom_line(color = colour_palette[5]) +
  facet_wrap(~ boat_ramp, ncol = 4) +
  labs(x = "Time (month index)", y = "Effort", title = "Commercial Boat Ramp Effort Over Time")
ggsave("plots/checking_plots_during_setup/04A_commercial_effort_over_time_by_ramp.png", plot = last_plot())

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

## X. Make a GIF for the laughs -----------------------------------------------

library(gifski)

frame_count <- 1

for (y in seq_along(years)) {
    
    # Update water$effort for this year and month
    water$effort <- c_fishing[, 1, y]
    
    # Year and month labels
    year_idx <- years[y]
    month_idx <- m
    
    # Plot
    p <- ggplot(water) +
      geom_sf(aes(fill = log(effort)), color = NA) +
      scale_fill_gradientn(colours = colour_palette[6:4], na.value = NA) +
      labs(
        title = paste("Fishing Effort - Year", year_idx, "Month", month_idx),
        fill = "log Effort"
      ) +
      theme_minimal()
    
    # Save frame
    ggsave(
      filename = sprintf("plots/gif_frames/commercial_effort_frame_%03d.png", frame_count),
      plot = p,
      width = 6, height = 6, dpi = 150
    )
    
    frame_count <- frame_count + 1
}

# stitch the GIF frames together
png_files <- list.files("plots/gif_frames", pattern = "commercial_effort_frame_\\d+\\.png", full.names = TRUE)
gifski(
  png_files,
  gif_file = "plots/gifs/commercial_fishing_effort_over_time.gif",
  width = 600,
  height = 600,
  delay = 0.5  # seconds per frame (adjust as needed)
)


## Set up effort for burn in --------------------------------------------------

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

# Split up this effort by the same proportions as before and allocate it to the different access points
burn_in_effort <- prop_month_ave[, 2] * effort

burn_in_effort <- as.data.frame(burn_in_effort) %>%
  rename(effort = "ave_month_prop")
for (location in location_name) {
  burn_in_effort[[location]] <- 0
} # Loop through each location name and create a new column with 0's

for (M in 1:12) { # for each month,
  for (i in 1:nrow(BR_trips)) { # and each boat ramp,
    location_name <- BR_trips$boat_ramp[i] # extract the location name
    burn_in_effort[[location_name]] <- effort * BR_trips$BR_prop[i] # and calculate burn in effort
  }
} # this loop calculate burn-in effort for all boat ramps, based on how many visits and hours each ramp gets fished from
# commercial boat and rec boat have different fishing days per boat ramp because commercial boats launch from a few ramps v. rec boats launch from many. the total is the same tho.

## Allocate to the cells using the same utilities that we set up earlier
c_burn_in_fishing <- array(0, dim = c(NCELL, 12, n_years_tot)) #This array has a row for every cell, a column for every month, and a layer for every year
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
  c_burn_in_fishing[ , , yr] <- months
  c_burn_in_fishing[ , , yr] <- c_burn_in_fishing[ , , yr] * spatial_q[, 1]
  yr <- yr + 1
} # this loop assigns each cell an amount of fishing effort each month and across years based on how useful it is (how far from a boatramp and how frequented that boat ramp is)


saveRDS(c_burn_in_fishing, file = "data/output_data/04A_commercial_burn_in_fishing.rds")

# Charlotte then adds another burn-in array for high mortality, I'm not doing rn for the sake of moving along

## END ##
