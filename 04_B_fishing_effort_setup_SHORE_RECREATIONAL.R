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

library(tidyverse) # data manipulation
library(sf) # shapefiles
#library(raster) # bathy rasters
library(terra) # for fast raster computations
#library(stringr)
#library(forcats)
library(RColorBrewer) # plotting colours
#library(geosphere)
library(abind) # dealing with matrices
library(sfnetworks) # distance from cell to cell
#library(purrr)
library(exactextractr) # extracting raster values


## Custom plotting parameters -------------------------------------------------
colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))
source("custom_theme.R")

## 0. Files used in this script -----------------------------------------------

file_wa         <- "data/output_data/01_B_land.shp"
#file_ntz        <- "data/output_data/01_B_wadandi_NTZ.shp"
file_bathy      <- "data/input_data/SW_crop_AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
file_carpark_w  <- "data/input_data/wadandi_carparks.shp"
file_carpark_n  <- "data/input_data/north_carparks.shp"
file_water      <- "data/output_data/03_water.rds"
file_network    <- "data/output_data/03_network_shapefile.shp"

# load land and ntz for sanity checks throughout
land <- st_read(file_wa); plot(land$geometry)
#ntz <- st_read(file_ntz); plot(ntz$geometry, add = T)

year_start <- 1900
year_end <- 2024
n_years_tot <- year_end - year_start +1

crs_raster <- "+proj=longlat +datum=WGS84 +no_defs"
bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -33.2), crs = crs_raster)

## 1. Catchability ------------------------------------------------------------


## catchability is the proportion of available fish in a population that would be captured by a unit of effort. (van Oostenbrugge et al. 2008)
## in a grid, each cell has a portion of the catchability of the whole grid, which we need to calculate.
## in this section (1.) we calculate the fishable area of the cells in each time step (1.a), then divide it by the fishable area of the whole grid at each time step (1.b)
## which gives us each the portion of Catchability of each cell.
## because the area that is catchable changes (with spatial and temporal restrictions), we have to calculate the grid's catchability for each month of each year.


### 1.a. find the fishable area of each cell in each time step ----------------

water <- readRDS(file_water)# %>% filter(!is.na(ID))
NCELL <- nrow(water)

# identify the important cells
temporal_cells <- water$ID[water$status == "TC"]
offshore_cells <- water$ID[!water$type %in% c("shore_north", "shore_wadandi")]

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
    fleet_allowed <- water$boat_rec # if FALSE, cell is not fishable for this fleet
    
    restriction_active <- !is.na(restriction_dates) & current_year >= restriction_dates
    fleet_blocked <- !is.na(fleet_allowed) & fleet_allowed == FALSE
    restricted <- restriction_active | fleet_blocked
    
    water_area[, MONTH, YEAR] <- ifelse( # if there are restrictions, fishable area = 0.
      restricted,
      0,
      water$cell_area
    )
    
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
          water_area[CELL, MONTH, YEAR] <- water$cell_area[CELL] * TC_perc_current[month_idx]
        }
      }
    }
  }
} # this loop calculates for every month and every year, the fishable area of each cell.

# override the fishability of offshore cells - here for shore fishing, none of them are fishable.
water_area[offshore_cells,,] <- 0

### sanity check station

test_cell <- 200
test_month <- 10
test_year <- 125

# reference cell numbers as of 13.01.2026:
## cockburn sound cell (temporal closure): 1
## SWC NTZ cell: 200
## random fished cell: 1000

# restrictions in this cell should be:
glimpse(st_drop_geometry(water[test_cell, c("status", "SC_restriction_date", "TC_restriction_date", "TC_restriction_months", "TC_restriction_perc_fished")]))

# check that it is correct
cat("In month", test_month, "of year", (year_start+test_year),
    ", the cell is", ifelse(water_area[test_cell, test_month, test_year]>0, "fishable.", "NOT fishable."), "Fishable area: ", water_area[test_cell, test_month, test_year]/1e06, "km2")


### 1.b. divide the fishable area by the grid sum area ------------------------

water_q <- array(0.000006, dim = c(NCELL, 12, (year_end-year_start+1)))

for (YEAR in 1:dim(water_q)[3]) {
  for (MONTH in 1:dim(water_q)[2]) {
    
    water_q[, MONTH, YEAR] <- water_area[, MONTH, YEAR] / sum(water_area[1:dim(water_q)[2]])
  }
} # this loop calculates the catchability of each cell for every year, based on the fishable area at that period.

summary(water_q)


### sanity check station 

test_cell <- 1
test_month <- 10
test_year <- 120

# reference cell numbers as of 13.01.2026:
## cockburn sound cell (temporal closure): 1
## SWC NTZ cell: 200
## random fished cell: 1000

catch_df <- as.data.frame(water_q[, test_month, test_year])
names(catch_df) <- "test_catchability"
catch_df$ID <- water$ID

water_catch <- water %>% left_join(catch_df, by = "ID")  # ID must be in 'water' too

ggplot(water_catch) +
  geom_sf(aes(fill = test_catchability), color = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  labs(title = paste0("Catchability in ", (year_start+test_year), ", month ", test_month), fill = "Catchability") +
  theme_minimal()

ggsave("plots/checking_plots_during_setup/04_B_shore_rec_catchability_test_year.png", plot = last_plot(), height = 15, width = 5)

## save for future use --------------------------------------------------------

saveRDS(water_q, file = "data/output_data/04_B_shore_rec_spatial_q_NTZ.rds")

#water_q <- readRDS("data/output_data/04_B_shore_rec_spatial_q_NTZ.rds")


## 2. Fishing days ------------------------------------------------------------

g_sheets <- read.csv("data/input_data/YIJARUP - Fishing effort reconstruction - FINAL_boat_days_total (19-01-2026).csv", skip = 1) %>% 
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
  year = year_start:year_end,
  month = sprintf("%02d", 1:12)
) %>%
  arrange(year, month) %>%
  mutate(
    annual_shore_days = rep(annual_effort_shore$shore_days, each = 12),
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

saveRDS(shore_effort, "data/output_data/04_B_shore_rec_total_boat_days.rds")

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
saveRDS(prop_month_ave, "data/output_data/04_B_shore_prop_month_ave.rds") # charlotte's 'Average_Monthly_Effort"

# unlike boat fishing i am assigning each carpark the same proportion of effort in north and wadandi, because no data anywhere.
CP_n <- st_read(file_carpark_n) %>% 
  mutate(id = 1:nrow(.),
         area = paste0(rough_area, "_", id)) %>% 
  rename(year_strt = year_start) %>% 
  dplyr::select(c(id, year_strt, area)) %>%
  glimpse()
CP_w <- st_read(file_carpark_w) %>% 
  dplyr::select(c(id, year_start, area)) %>%
  rename(year_strt = year_start) %>% 
  glimpse()


CP <- rbind(CP_n, CP_w) %>% 
  st_transform(4283) %>%
  st_make_valid() %>%
  mutate(year_strt = as.numeric(year_strt),
         year_strt = ifelse(is.na(year_strt), year_start, year_strt), # fill in missing dates with the simulation start year (here year_start is 1945, and year_strt is the build date)
         #build_mnth = ifelse(is.na(month_strt), 1, month_strt), # assume Jan if unknown
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
    #build_date = as.Date(paste(CP$year_strt[CP_index], CP$build_mnth[CP_index], "01", sep = "-")),
    build_date = as.Date(paste(CP$year_strt[CP_index], "01-01", sep = "-")), # for now build dates are all Jan, when i have more info, use line above.
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
saveRDS(shore_effort_df, "data/output_data/04_B_shore_rec_effort.rds")

# check each boat ramp is populated correctly
ggplot(shore_effort_df %>% filter(CP_index %in% c(1:9)), aes(x = year, y = adjusted_effort)) +
  geom_line(color = colour_palette[5]) +
  facet_wrap(~ CP_index, ncol = 3)

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

unique(water$type)

centroids <- st_centroid_within_poly(st_make_valid(water))
plot(water$geometry)

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

distCP <- as.data.frame(t(network_matrix))
colnames(distCP) <- CP$area
head(distCP) # this gives us each cell's distance to the carparks

# make cells beyond 2km distance not fishable
head(distCP)
distCP <- distCP %>%  
  mutate(across(where(is.numeric), ~ifelse(. > 2, NA, .)))

# make offshore cells not fishable
distCP[water$ID[!water$type %in% c("shore_north", "shore_wadandi")], ] <- NA

# Plot check
distCP <- distCP %>% mutate(ID = water$ID)
glimpse(distCP)
saveRDS(distCP, "data/output_data/04_B_shore_cell_dist.rds")

water_dist <- water %>% left_join(distCP, by = "ID")
glimpse(water_dist)
water_dist_long <- water_dist %>% pivot_longer(cols = CP$area, names_to = "access", values_to = "Distance_km")
glimpse(water_dist_long)

#water_dist_long <- water_dist_long[water_dist_long$type %in% c("shore_north", "shore_north_cockburn_warnbro", "shore_wadandi"),]
carpark_selection <- unique(water_dist_long$access)[c(55)]
st_crs(CP) <- st_crs(water_dist_long) 

ggplot() +
  geom_sf(data = water_dist_long %>% filter(access %in% carpark_selection),
          aes(fill = Distance_km), color = NA) +
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  labs(title = paste0("Cell distance to selected shore access point"), fill = "Distance (km)")

ggsave("plots/checking_plots_during_setup/04_B_shore_access_distance.png", plot = last_plot(),
       width = 10, height = 15, dpi = 2000, units = "in", device='png')


## 4. Set up a utility function -----------------------------------------------

## Now need to create a separate fishing surface for each month of each year based on distance to boat ramp, size of each cell,
## and multiply that by the effort in the cell to spatially allocate the effort across the area. Effort is also able to go more offshore over time.
## But we need to account for the fact that there will be sanctuary zones going in and the effort that would have gone in there will get allocated somewhere else
## Will then need to put the rows/columns back in as 0s 

# What we want to do is to distribute effort across month and year based on :
# 1. distance to boat ramp, 2. size of the cell, and 3. whether the cell is 'fishable' that year and 4. fishing effort
# After the SZ comes in the effort will also be redistributed.

# in boat fishing, each year a bit more area is fishable (deeper every year). this doesn't exist for shore fishing, so replace with 1.
fishable_long <- expand.grid(year = year_start:year_end,
                             ID = water$ID) %>% 
  mutate(fishable_prop = ifelse(ID %in% water$ID[water$type %in% c("shore_north", "shore_wadandi")], 1, 0)) # only shore cells are fishable.
  
years <- year_start:year_end
nyears <- length(years)

carparks <- unique(CP$area)
ncarparks <- length(carparks)

cell_ids <- sort(unique(fishable_long$ID))
ncells <- length(cell_ids)

# Get full list of cell IDs from fishable_long
all_cell_ids <- sort(unique(water$ID))
ncells <- length(all_cell_ids)
nyears <- length(years)

# Distance from each cell to each ramp (matrix)
# Assuming `distCP` has same row order as `water`
cell_dist <- distCP
saveRDS(cell_dist, "data/output_data/04_B_shore_rec_cell_dist.rds")

CP_U_array <- array(0, dim = c(ncells, ncarparks, 12, nyears),
                    dimnames = list(cell_id = all_cell_ids,
                                    ramp = carparks,
                                    month = 1:12,
                                    year = as.character(years)))

for (YEAR in seq_along(years)) {
  yr <- years[YEAR]
  
  fishable_depth <- fishable_long %>%
    filter(year == yr) %>%
    #rename(fishable_prop = fshbl_p) %>%  # fshbl_p = fishable_prop (saving a file earlier as .shp shortens colnames)
    dplyr::select(ID, fishable_prop)
  
  Vj_df <- cell_dist %>%
    inner_join(fishable_depth, by = "ID") %>%
    arrange(ID)
  
  row_ids <- match(Vj_df$ID, cell_ids)
  
  for (MONTH in 1:12) { # because each month's catchability can change (temporal closures etc.), get the correct one
    
    cell_area <- water_q[, MONTH, YEAR]
    U_mat <- matrix(0, nrow = nrow(Vj_df), ncol = ncarparks) # cells x ramps
    
    for (ACCESS in seq_along(carparks)) {
      carpark_name <- carparks[ACCESS]
      U_mat[, ACCESS] <- exp(-Vj_df[[carpark_name]]) * cell_area * Vj_df$fishable_prop # exp(-Vj_df...) because otherwise high utility is given to areas far from ramps
    }
    
    carpark_sums <- colSums(U_mat, na.rm = TRUE)
    U_norm <- sweep(U_mat, 2, carpark_sums, "/")
    
    CP_U_array[row_ids, , MONTH, YEAR] <- U_norm
    
  }
} # this loop calculates how useful a cell is to fishing, based on its area (the bigger, the more useful), its distance to ramps (the closer, the more useful), and its fishable depth status (if within fishable depth that year, useful)


plot(CP_U_array[row_ids,,,70])
head(CP_U_array[row_ids,,,70])
head(water)

# plot check
cp_slice <- CP_U_array[,,1,1]
summary(cp_slice[!is.na(cp_slice)])
row_sums <- rowSums(!is.na(cp_slice)) # calculate row sums (sum of utilities across ramps for each cell)
water$utility_sum <- row_sums # add the sums as a new column to the water sf object
ggplot(water) +
  geom_sf(aes(fill = utility_sum), color = NA) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", 
                       na.value = "grey80", name = "Sum of Utilities") +
  theme_minimal() +
  theme(legend.position = "right")


# plot check DOESNT WOTK FSR
carpark_check <- which(CP$area == "Mandurah_5")
utility_check <- data.frame(CP_U_array[, carpark_check, 1, ])  # dimensions: cells × ramp x month x years
names(utility_check) <- year_start:year_end
utility_check$ID <- water$ID#[water$type %in% c("shore_north", "shore_north_cockburn_warnbro", "shore_wadandi")]
head(utility_check)
water_catch <- water %>% left_join(utility_check, by = "ID")

water_catch_long <- water_catch %>% pivot_longer(cols = as.character(years), names_to = "year", values_to = "Catchability")
plot(water_catch_long$Catchability)
ggplot(water_catch_long %>% dplyr::filter(year %in% c(1950:1969))) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  facet_wrap(~ year, ncol = 10) + 
  labs(title = paste0("commercial catchability of each cell, \nby distance from ", carpark_check, " ramp, cell size and tech-fishability"), fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/04_A_commercial_boat_ramp_catchability_over_time.png", plot = last_plot())

# Plot check DOESNT WORK
utility_check2 <- as.data.frame(CP_U_array[,,80]) %>% mutate(ID = as.numeric(rownames(CP_U_array[,,80])))
class(utility_check2$ID)
water_catch <- water %>% left_join(utility_check2, by = "ID")
water_catch_long <- water_catch %>% pivot_longer(cols = CP$area, names_to = "carpark", values_to = "Catchability")

ggplot(water_catch_long) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  facet_wrap(~ carpark, ncol = 8) + 
  labs(title = "Catchability of each cell in year 80, by boat ramp distance and cell size", fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/04_A_commercial_catchability_surface_by_ramp.png", plot = last_plot())


## 5. Allocating effort to cells ----------------------------------------------

# We know the fishing effort over the years
# We know the fishing effort across months
# We know how 'useful' each cell is for fishing (utility function)
# Now we want to allocate fishing effort to cells over time
CP_trips <- shore_effort_df %>% # This is just the trips from each carpark
  arrange(year, month) %>% 
  mutate(num_year = match(year, sort(unique(year)))) %>% # This is to number the years 1 to 80 for the loop.
  glimpse()
saveRDS(CP_trips, "data/output_data/04_B_CP_trips.rds")
#CP_trips <- readRDS("data/output_data/04_B_CP_trips.rds")
CP_trips <- CP_trips %>%
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
ggsave("plots/checking_plots_during_setup/04_B_shore_rec_effort_over_time_by_carpark.png", plot = last_plot())

# link the carparks to their cells
cell_access_pairs <- water_dist_long %>%
  filter(!is.na(Distance_km)) %>%
  dplyr::select(ID, access) %>% 
  st_drop_geometry() %>% 
  glimpse()




s_fishing <- array(0, dim = c(NCELL, 12, length(years))) # this array has a row for every cell, a column for every month, and a layer for every year
months <- array(0, dim = c(NCELL, 12)) # array for the number of months
carparks <- array(0, dim = c(NCELL, length(CP$area))) # array for the number of ramps

head(CP_trips)
head(CP)
head(CP_U_array[,,1,1])
CP_U_array[is.na(CP_U_array)] <- 0

shore_cells <- water$type %in% c(
  "shore_north",
  "shore_wadandi"
)
for(YEAR in 1:length(years)){
  
  print(YEAR)
  
  temp <- CP_trips %>% 
    filter(num_year == YEAR) %>% 
    dplyr::select(-c(num_year, month))
  temp <- as.matrix(temp)
      
  for(MONTH in 1:12){
    
      #carparks[,] <- 0

    for(CARPARK in 1:length(CP$area)){
      
      carparks[, CARPARK] <- CP_U_array[, CARPARK, MONTH, YEAR] * temp[MONTH, CARPARK]
      
    }
    
    months[, MONTH] <- rowSums(carparks, na.rm = TRUE)
    
  }
  s_fishing[ , , YEAR] <- months
} # this loop assigns each cell a fishing effort based on utility by month
s_fishing

summary(s_fishing[,,1])


# plot check
year_idx <- 125  # choose the year to plot
current_year <- year_start-1 + sort(unique(CP_trips$num_year))[year_idx]
current_month <- 1 # choose the month to plot
built_CPs <- CP %>% filter(year_strt <= current_year);built_CPs <- st_transform(built_CPs, 4326) # filter CPs that are already built
effort_vec <- s_fishing[, 1, year_idx] # extract that year&month's effort as a vector (one value per cell)
water$effort <- effort_vec # add that effort to the grid

ggplot(water) +
  geom_sf(aes(fill = (effort)), col = NA) +
  #geom_sf(data = built_CPs, color = "red", size = 0.5) +  # only plot active carparks
  scale_fill_gradientn(colours = c("blue", "red"), na.value = NA) +
  labs(title = paste("Shore fishing effort check - Year", current_year, "Month", current_month),
       fill = "log Effort") +
  theme_minimal()


## X. Make a GIF --------------------------------------------------------------

library(gifski)


frame_count <- 1

# make the plot's colour limits otherwise log(0) wouldnt work in the plot
log_effort_all <- log(s_fishing[s_fishing > 0]) # to avoid -Inf as a limit in the plot which wouldn't work
global_limits <- range(log_effort_all, na.rm = TRUE)

for (y in seq_along(years)) {
  year_idx <- years[y]
  print(y)
  for (m in 1:12) {
    month_idx <- m
    
    # update effort for this year & month
    water$effort <- s_fishing[, m, y]
    
    p <- ggplot(water) +
      geom_sf(aes(fill = log(effort)), colour = NA) +
      scale_fill_gradientn(
        colours = colour_palette[6:4],
        limits  = global_limits,
        oob     = scales::squish,  # very important
        na.value = NA
      )
    labs(
      title = paste("shore fishing Effort – Year", year_idx, "Month", month_idx),
      fill = "log Effort"
    ) +
      theme_minimal()
    
    ggsave(
      filename = sprintf(
        "plots/gif_frames/04_B_shore_fishing_setup_gif_frames/shore_effort_frame_%03d.png",
        frame_count
      ),
      plot = p,
      width = 6,
      height = 6,
      dpi = 150
    )
    
    frame_count <- frame_count + 1
  }
}

# stitch the GIF frames together
png_files <- list.files("plots/gif_frames/04_B_shore_fishing_setup_gif_frames", pattern = "shore_effort_frame_\\d+\\.png", full.names = TRUE)
gifski(
  png_files,
  gif_file = "plots/gifs/commercial_fishing_effort_over_time.gif",
  width = 600,
  height = 600,
  delay = 0.05  # seconds per frame (adjust as needed)
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