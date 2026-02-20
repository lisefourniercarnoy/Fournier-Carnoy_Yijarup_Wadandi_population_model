# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Various literature figures.
# Task:    Set up fishing effort (understand the split of fishing effort within a year, over space, and across years)
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    May 2025

# -----------------------------------------------------------------------------

# Notes: This script goes like this:
# 1. Calculate how 'catchable' each cell is each year (based on how big the cell is and whether it's NTZ or not)
# 2. Use literature info to hindcast fishing effort trends since 1900
# 3. Make fishability of cells (whether they can be fished based on the tech at that time) increase over time
# 4. Calculate the distance from access points to each cell (do this separately for shore and boat fishing)
# 5. Set up a utility function (how 'useful' each cell is, based on its catchability and distance to the access points)
# 6. Allocate the fishing effort (from step 2) to each cell by how 'useful' it is (step 4)

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

file_wa           <- "data/output_data/01_B_land.shp"
file_bathy        <- "data/input_data/SW_crop_AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
file_boat_ramps_w <- "data/input_data/wadandi_boat_ramps.shp"
file_boat_ramps_n <- "data/input_data/north_boat_ramps.shp"
file_water        <- "data/output_data/03_water.rds"
file_network      <- "data/output_data/03_network_shapefile.shp"

year_start <- 1900
year_end <- 2024
n_years_tot <- year_end - year_start +1

crs_raster <- "+proj=longlat +datum=WGS84 +no_defs"
bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -31), crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(crs_raster)


## 1. Catchability ------------------------------------------------------------

## catchability is the proportion of available fish in a population that would be captured by a unit of effort. (van Oostenbrugge et al. 2008)
## in a grid, each cell has a portion of the catchability of the whole grid, which we need to calculate.
## in this section (1.) we calculate the fishable area of the cells in each time step (1.a), then divide it by the fishable area of the whole grid at each time step (1.b)
## which gives us each the portion of Catchability of each cell.
## because the area that is catchable changes (with spatial and temporal restrictions), we have to calculate the grid's catchability for each month of each year.


### 1.a. find the fishable area of each cell in each time step ----------------

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
      !is.na(water$boat_rec) & water$boat_rec == FALSE & # where fleet is not allowed,
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

test_cell <- 200
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
    ", the cell is", ifelse(water_area[test_cell, test_month, test_year]>0, "fishable", "NOT fishable"), "for commercial boats.", "Fishable area: ", water_area[test_cell, test_month, test_year]/1e06, "km2")


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
  labs(title = paste0("Comm. Catchability in ", (year_start+test_year), ", month ", test_month), fill = "Catchability") +
  theme_minimal()

ggsave("plots/checking_plots_during_setup/04_A_commercial_catchability_test_year.png", plot = last_plot(), height = 15, width = 5)

## save for future use --------------------------------------------------------

saveRDS(water_q, file = "data/output_data/04_A_commercial_spatial_q_NTZ.rds")
water_q <- readRDS("data/output_data/04_A_commercial_spatial_q_NTZ.rds")


## 2. Fishing days -------------------------------------------------------------

# we need to know how much fishing occurs in the region, in order to distribute it in the grid.
# for commercial fishing, we are using a variety of sources to reconstruct the trends (a),
# the sources are split between the Metropolitan area effort (North, 2.N) and the Southwest area effort (Wadandi, 2.W)
# we will (b.) split the yearly effort into months, then (c.) split the monthly effort into access points (boat ramps).



### 2.N. North (Metropolitan) -------------------------------------------------

#### 2.N.a. enter the overall effort values (boat days) -----------------------

# obtain boat days from the literature (see google slides on reconstruction)
years     <- c(1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,
               1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,
               2001,2002,2003,2004,2005)
boat_days <- c(4000,3500,3000,3500,4000,4250,4250,4250,5250,5250,5500,3500,4500,
               4500,3250,3000,2500,2750,2750,3000,2500,2500,4000,3750,3500,3000,
               3750,3750,3750,3750,3750)

boat_days_lit_n <- data.frame(years, boat_days)

plot(boat_days_lit_n$years, boat_days_lit_n$boat_days, 
     pch = 19, col = colour_palette[4], xlim = c(year_start, year_end), ylim = c(0, max(boat_days_lit_n)), 
     xlab = "Year", ylab = "Boat Days", main = "Original + Fake Boat Days")

# fill in with fake values
years     <- c(1900,1950,2024)
boat_days <- c(500, 1700,3750)
boat_days_fake <- data.frame(years, boat_days)

points(boat_days_fake$years, boat_days_fake$boat_days, 
       pch = 17, col = colour_palette[6])
legend("topright", legend = c("Original Data", "Fake Data"), col = c(colour_palette[4], colour_palette[6]), pch = c(19, 17))


# fill in the gaps to obtain values for every year
years <- c(boat_days_lit_n$years, boat_days_fake$years)
boat_days <- c(boat_days_lit_n$boat_days, boat_days_fake$boat_days)

sorted_index <- order(years) # sort data in order (important for interpolation)
all_years_sorted <- years[sorted_index]
all_boat_days_sorted <- boat_days[sorted_index]

years_full <- year_start:year_end
# linear interpolation with extrapolation
interp <- approx(x = all_years_sorted, y = all_boat_days_sorted, xout = years_full, method = "linear", rule = 2)
annual_effort_n <- data.frame(year = interp$x, boat_days = interp$y)


# sanity check station
plot(annual_effort_n$year, annual_effort_n$boat_days, type = "l", col = colour_palette[5], lwd = 4,
     main = "Boat Days in Metro area \nfrom the Literature \n(with some extrapolation)", 
     xlab = "Year", ylab = "Boat Days")
lines(boat_days_lit_n$years, boat_days_lit_n$boat_days, col = colour_palette[4], lwd = 4)
legend("topleft", legend = c("From literature", "Fake"), col = c(colour_palette[4], colour_palette[5]), lwd = 4)


#### 2.N.b. split yearly fishing effort by month ------------------------------

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
  year = years_full,
  month = sprintf("%02d", 1:12)
) %>%
  arrange(year, month) %>%
  mutate(
    annual_boat_days = rep(annual_effort_n$boat_days, each = 12),
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

plot(prop_month_ave)
saveRDS(prop_month_ave, "data/output_data/04_A_commercial_metro_prop_month_ave.rds") # charlotte's 'Average_Monthly_Effort"

#### 2.N.c. distribute monthly effort into boat ramps -------------------------

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



### 2.W. Wadandi (Southwest) --------------------------------------------------

#### 2.W.a. enter the overall effort values (boat days) -----------------------

# obtain boat days from the literature (see google slides on reconstruction)
years     <- c(1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,1985,1986,1987,
               1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,
               2001,2002,2003,2004,2005)
boat_days <- c(3000,3250,3000,3500,2750,3500,3250,3750,4000,4500,4500,3500,3250,
               2250,1500,1250,1250,1000,1000,1250,1000,1000,1500,1250,1250,1250,
               1500,1500,1500,1750,1500)
boat_days_lit_w <- data.frame(years, boat_days)

plot(boat_days_lit_w$years, boat_days_lit_w$boat_days, 
     pch = 19, col = colour_palette[4], xlim = c(year_start, year_end), ylim = c(0, max(boat_days_lit_w)), 
     xlab = "Year", ylab = "Boat Days", main = "Original + Fake Boat Days")

# adding fake boat days to fill out.
years     <- c(1900,1950,1500)
boat_days <- c(400, 1500,2024)
boat_days_fake <- data.frame(years, boat_days)

points(boat_days_fake$years, boat_days_fake$boat_days, pch = 17, col = colour_palette[6])
legend("topright", legend = c("Original Data", "Fake Data"), col = c(colour_palette[4], colour_palette[6]), pch = c(19, 17))

# fill in the gaps to obtain values for every year
years <- c(boat_days_lit_w$years, boat_days_fake$years)
boat_days <- c(boat_days_lit_w$boat_days, boat_days_fake$boat_days)

sorted_index <- order(years) # sort data in order (important for interpolation)
all_years_sorted <- years[sorted_index]
all_boat_days_sorted <- boat_days[sorted_index]

years_full <- year_start:year_end
# linear interpolation with extrapolation
interp <- approx(x = all_years_sorted, y = all_boat_days_sorted, xout = years_full, method = "linear", rule = 2)
annual_effort_w <- data.frame(year = interp$x, boat_days = interp$y)

# sanity check
plot(annual_effort_w$year, annual_effort_w$boat_days, type = "l", col = colour_palette[5], lwd = 4,
     main = "Boat Days in Wadandi Country \nfrom the Literature \n(with some extrapolation)", 
     xlab = "Year", ylab = "Boat Days")
lines(boat_days_lit_w$years, boat_days_lit_w$boat_days, col = colour_palette[4], lwd = 4)
legend("topleft", legend = c("From literature", "Fake"), col = c(colour_palette[4], colour_palette[5]), lwd = 4)


#### 2.W.b. split yearly fishing effort by month ------------------------------

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
  year = years_full,
  month = sprintf("%02d", 1:12)
) %>%
  arrange(year, month) %>%
  mutate(
    annual_boat_days = rep(annual_effort_w$boat_days, each = 12),
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
  mutate(ave_month_prop = mean(month_prop))
prop_month_ave <- boat_month_prop[1:12, c(2, 5)]

plot(prop_month_ave)
saveRDS(prop_month_ave, "data/output_data/04_A_commercial_prop_month_ave.rds") # charlotte's 'Average_Monthly_Effort"

# sanity check station
plot(annual_effort_n$year, annual_effort_n$boat_days, type = "l", col = colour_palette[5], lwd = 4,
     main = "Commercial Boat Days in \nWadandi Country (dashed) and Metro (solid)", 
     xlab = "Year", ylab = "Boat Days")
lines(boat_days_lit_n$years, boat_days_lit_n$boat_days, col = colour_palette[4], lwd = 4)
legend("topleft", legend = c("From literature", "Fake"), col = c(colour_palette[4], colour_palette[5]), lwd = 4)

lines(annual_effort_w$year[annual_effort_w$year < min(boat_days_lit_w$years)], annual_effort_w$boat_days[annual_effort_w$year < min(boat_days_lit_w$years)], lty = "dashed", col = colour_palette[5], lwd = 2)
lines(annual_effort_w$year[annual_effort_w$year > max(boat_days_lit_w$years)], annual_effort_w$boat_days[annual_effort_w$year > max(boat_days_lit_w$years)], lty = "dashed", col = colour_palette[5], lwd = 2)

lines(boat_days_lit_w$years, boat_days_lit_w$boat_days, lty = "dashed", col = colour_palette[4], lwd = 2)
legend("topleft", legend = c("From literature", "Fake"), col = c(colour_palette[4], colour_palette[5]), lwd = 4)


#### 2.W.c. distribute monthly effort into boat ramps -------------------------

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


## 3. "boats can fish further over time" --------------------------------------

# this section of the fishing effort relates to boats being able to fish further over the years, with bigger boats, more powerful engines etc.
# the model takes this into account by adding a 'fishable depth' into the mix, which goes more and more offshore.
# for commercial fishing, we're assuming that on average, fishers can fish 1.3m more every year. (change as needed in parameters)

# prepare grid cells
plot(water$geometry)
wa_mask <- st_read(file_wa); wa_mask <- st_transform(wa_mask, crs = crs_raster); wa_mask <- as(wa_mask, "Spatial"); plot(wa_mask, col = "lightgray", add = T)
plot(bbox, add = T)

# set parameters
years <- year_start:year_end
min_depth <- 20 # fishable depth is about 20m in 1945 (from Gaynor 2008, p.38)
increase_year <- 1950 # the year that boats start to go deeper (post-war industrialisation)
depth_per_year <- 1.3 # fishable depth gained per year (m)

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

fishable_summary <- terra::extract(
  fishable_stack,
  vect(water),
  fun = mean,
  na.rm = TRUE
)

fishable_summary$ID <- water$ID

fishable_long <- fishable_summary %>%
  pivot_longer(
    cols = starts_with("year_"),
    names_to = "year",
    names_prefix = "year_",
    values_to = "fishable_prop"
  ) %>%
  mutate(year = as.integer(year)) %>%
  left_join(water, by = "ID") %>%
  st_as_sf()


## sanity check station 

plot(fishable_stack[[c("year_1940", "year_1960", "year_1975", "year_1990")]])

fishable_long %>%
  filter(year == 1965) %>% 
  ggplot() +
  geom_sf(aes(fill = fishable_prop), colour = NA) +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  labs(y = "Fishable proportion", colour = "Cell ID") +
  theme_minimal()


## save for future use --------------------------------------------------------

st_write(fishable_long %>% dplyr::select(!where(is.list)), "data/output_data/04_A_commercial_fishable_area_over_time.shp", append = FALSE)

fishable_long <- st_read("data/output_data/04_A_commercial_fishable_area_over_time.shp")


## 4. Access point distance to cells ------------------------------------------

## this section (4.) calculates each cell's distance to each access point (boat ramp)

glimpse(BR_n)
glimpse(BR_w)

BR <- rbind(BR_n %>% dplyr::select(name, geometry) %>% mutate(region = "north"), 
            BR_w %>% dplyr::select(name, geometry) %>% mutate(region = "wadandi")) %>% 
  st_transform(4283) %>%
  st_make_valid() %>% 
  glimpse()

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

network <- st_read(file_network); plot(network$geometry)
BR <- st_as_sf(BR); st_crs(BR) <- NA 

# find the centre of each grid cell
sf::sf_use_s2(FALSE)
centroids <- st_centroid(st_make_valid(water))
sf::sf_use_s2(TRUE)
points <- as.data.frame(st_coordinates(centroids))%>% # the points start at the bottom left and then work their way their way right
  mutate(ID = row_number()) 
points_sf <- st_as_sf(points, coords = c("X", "Y")) 
st_crs(points_sf) <- NA 

network <- as_sfnetwork(network, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

net <- activate(network, "nodes")
st_crs(net)
st_crs(net) <- NA 

# measure the distance from access points to cell centroids
network_matrix <- st_network_cost(net, from = BR, to = points_sf)
network_matrix <- network_matrix * 111
dim(network_matrix) # number of ramps x number of cells

glimpse(network_matrix)
DistBR <- as.data.frame(t(network_matrix))
colnames(DistBR) <- BR$name
head(DistBR) # this gives us each cell's distance to the boat ramps

## sanity check station
DistBR <- DistBR %>% mutate(ID = centroids$ID)
water_dist <- water %>% left_join(DistBR, by = "ID")
water_dist_long <- water_dist %>% pivot_longer(cols = BR$name, names_to = "Ramp", values_to = "Distance_km")
ggplot(water_dist_long) + 
  geom_sf(aes(fill = Distance_km), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  facet_wrap(~ Ramp, ncol = 8) + 
  labs(title = "Cell distance to each boat ramp", fill = "Distance (km)") +
  theme_minimal()

ggsave("plots/checking_plots_during_setup/04_A_commercial_boat_ramp_distance.png", plot = last_plot(), height = 10, width = 10)


## 5. Create a utility function -----------------------------------------------

## Now need to create a separate fishing surface for each month of each year based on distance to boat ramp, size of each cell,
## and multiply that by the effort in the cell to spatially allocate the effort across the area. Effort is also able to go more offshore over time.
## But we need to account for the fact that there will be sanctuary zones going in and the effort that would have gone in there will get allocated somewhere else
## Will then need to put the rows/columns back in as 0s 

# What we want to do is to distribute effort across month and year based on :
# 1. distance to boat ramp, 2. size of the cell, and 3. whether the cell is 'fishable' that year and 4. fishing effort
# After the SZ comes in the effort will also be redistributed.

years <- year_start:year_end
nyears <- length(years)

ramps <- unique(BR$name)
nramps <- length(ramps)

cell_ids <- sort(unique(fishable_long$ID))
ncells <- length(cell_ids)

# Distance from each cell to each ramp (matrix)
summary(DistBR$ID == water$ID)# assuming `DistBR` has same row order as `water`
cell_dist <- DistBR
water_q

saveRDS(cell_dist, "data/output_data/04_A_commercial_cell_dist.rds")

BR_U_array <- array(0, dim = c(ncells, nramps, 12, nyears),
                    dimnames = list(cell_id = cell_ids,
                                    ramp = ramps,
                                    month = 1:12,
                                    year = as.character(years)))

for (YEAR in seq_along(years)) {
  yr <- years[YEAR]
  
  fishable_depth <- fishable_long %>%
    filter(year == yr) %>%
    rename(fishable_prop = fshbl_p) %>%  # fshbl_p = fishable_prop (saving a file earlier as .shp shortens colnames)
    dplyr::select(ID, fishable_prop)
  
  Vj_df <- cell_dist %>%
    inner_join(fishable_depth, by = "ID") %>%
    arrange(ID)
  
  row_ids <- match(Vj_df$ID, cell_ids)
  
  for (MONTH in 1:12) { # because each month's catchability can change (temporal closures etc.), get the correct one
    
    cell_area <- water_q[, MONTH, YEAR]
    U_mat <- matrix(0, nrow = nrow(Vj_df), ncol = nramps) # cells x ramps
    
    for (RAMP in seq_along(ramps)) {
      ramp_name <- ramps[RAMP]
      U_mat[, RAMP] <- exp(-Vj_df[[ramp_name]]) * cell_area * Vj_df$fishable_prop # exp(-Vj_df...) because otherwise high utility is given to areas far from ramps
    }
    
    ramp_sums <- colSums(U_mat, na.rm = TRUE)
    U_norm <- sweep(U_mat, 2, ramp_sums, "/")
    
    BR_U_array[row_ids, , MONTH, YEAR] <- U_norm
    
  }
} # this loop calculates how useful a cell is to fishing, based on its area (the bigger, the more useful), its distance to ramps (the closer, the more useful), and its fishable depth status (if within fishable depth that year, useful)
head(BR_U_array[,,,70])
head(water)

## sanity check plot (utility across all ramps)
br_slice <- BR_U_array[,,,125]
row_sums <- rowSums(br_slice) # calculate row sums (sum of utilities across ramps for each cell)
water$utility_sum <- row_sums # add the sums as a new column to the water sf object
ggplot(water) +
  geom_sf(aes(fill = log(utility_sum)), color = NA) +
  scale_fill_gradientn(colours = colour_palette[6:3]) +
  theme_minimal() +
  theme(legend.position = "right")


## sanity check plot (utility for a single ramp over time)
ramp_check <- "GB_Bunbury_Stirling_St"
utility_check <- as.data.frame(BR_U_array[, ramp_check, 1, ])  # dimensions: cells × ramp x month x years
utility_check$ID <- water$ID
head(utility_check)

water_catch <- water %>% left_join(utility_check, by = "ID")

water_catch_long <- water_catch %>% pivot_longer(cols = as.character(years), names_to = "year", values_to = "Catchability")
ggplot(water_catch_long %>% dplyr::filter(year %in% c(1940, 1960, 1980, 2000, 2024))) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  facet_wrap(~ year, ncol = 10) + 
  labs(title = paste0("commercial catchability, \nby cell size, depth-fishability, and \ndistance from ", ramp_check, " ramp"), fill = "log(Catchability)")+
  theme_minimal()

ggsave("plots/checking_plots_during_setup/04_A_commercial_boat_ramp_catchability_over_time.png", plot = last_plot(), height = 10, width = 10)


## 6. Allocating effort to cells ----------------------------------------------

## this section allocates the fishing effort from section 2., according to the utility of each cell.

# we know the fishing effort over the years
glimpse(ramp_effort_n)
glimpse(ramp_effort_w)
ramp_effort <- rbind(ramp_effort_n, ramp_effort_w)

# we know the fishing effort across months
# we know how 'useful' each cell is for fishing (utility function)
glimpse(BR_U_array)

# now we want to allocate fishing effort to cells over time
BR_trips <- ramp_effort %>% # This is just the trips from each boat ramp
  arrange(year, month) %>%
  mutate(num_year = match(year, sort(unique(year)))) %>% # This is to number the years 1 to 80 for the loop.
  glimpse()
saveRDS(BR_trips, "data/output_data/04_A_BR_trips.rds")

ggplot(BR_trips, aes(x = order(year, month), y = adjusted_effort)) +
  geom_line(color = colour_palette[5]) +
  facet_wrap(~ boat_ramp, ncol = 4) +
  labs(x = "Time (month index)", y = "Effort", title = "Commercial Boat Ramp Effort Over Time") +
  theme_minimal()
ggsave("plots/checking_plots_during_setup/04_A_commercial_effort_over_time_by_ramp.png", plot = last_plot())

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

head(BR_trips)
head(BR)
head(BR_U_array[,,,1])

for(YEAR in 1:length(years)){ # for all years,
  
  print(YEAR)
  
  temp <- BR_trips %>%
    filter(num_year == YEAR) %>% 
    dplyr::select(-c(num_year, month))
  temp <- as.matrix(temp)
  
  for(MONTH in 1:12){ # run the loop for every month,
    
    for(RAMP in 1:length(BR$name)){ # and for every ramp,
      
      ramps[,  RAMP] <- BR_U_array[, RAMP, MONTH, YEAR] * temp[MONTH, RAMP] # assign fishing effort to every cell, based on its utility
      
    }      
    months[, MONTH] <- rowSums(ramps)
  }
  c_fishing[ , , YEAR] <- months 
} # this loop assigns each cell a fishing effort based on utility by month
head(c_fishing)


# plot check
water
year_idx <- 125  # first year
month_idx <- 1  # January

effort_vec <- c_fishing[, month_idx, year_idx] # Extract effort vector (one value per cell)
water$effort <- effort_vec

ggplot(water) +
  geom_sf(aes(fill = log(effort)), color = NA) +
  labs(title = paste("Fishing Effort - Year", year_idx, "Month", month_idx),
       fill = "log Effort") +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()

## X. Make a GIF for the laughs -----------------------------------------------

library(gifski)

frame_count <- 1
m <- 1
global_limits <- range(log(c_fishing + 1), na.rm = TRUE)

for (y in seq_along(years)) {
  year_idx <- years[y]
  print(y)
  for (m in 1:12) {
    month_idx <- m
    
    # update effort for this year & month
    water$effort <- c_fishing[, m, y]
    
    p <- ggplot(water) +
      geom_sf(aes(fill = log(effort)), colour = NA) +
      scale_fill_gradientn(
        colours = colour_palette[6:4],
        limits  = global_limits,
        oob     = scales::squish,  # very important
        na.value = NA
      )
      labs(
        title = paste("Fishing Effort – Year", year_idx, "Month", month_idx),
        fill = "log Effort"
      ) +
      theme_minimal()
    
    ggsave(
      filename = sprintf(
        "plots/gif_frames/04_A_commercial_fishing_setup_gif_frames/commercial_effort_frame_%04d.png",
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
png_files <- list.files("plots/gif_frames/04_A_commercial_fishing_setup_gif_frames", pattern = "commercial_effort_frame_\\d+\\.png", full.names = TRUE)
gifski(
  png_files,
  gif_file = "plots/gifs/commercial_fishing_effort_over_time.gif",
  width = 600,
  height = 600,
  delay = 0.05  # seconds per frame (adjust as needed)
)


## Set up effort for burn in --------------------------------------------------

# the burn-in effort is a low level of effort, to allow the population to stabilise before the simulations start

access_points_name <- unique(BR$name)

BR_trips <- ramp_effort %>%
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
for (location in access_points_name) {
  burn_in_effort[[location]] <- 0
} # Loop through each location name and create a new column with 0's

for (M in 1:12) { # for each month,
  for (i in 1:nrow(BR_trips)) { # and each boat ramp,
    location_name <- BR_trips$boat_ramp[i] # extract the location name
    burn_in_effort[[location_name]] <- effort * BR_trips$BR_prop[i] # and calculate burn in effort
  }
} # this loop calculate burn-in effort for all boat ramps, based on how many visits and hours each ramp gets fished from

# Allocate to the cells using the same utilities that we set up earlier
n_years_burn_in <- 60 # this is the number of years the burn-in should run for, it should be >= max age of the species.
c_burn_in_fishing <- array(0, dim = c(NCELL, 12, n_years_burn_in)) #This array has a row for every cell, a column for every month, and a layer for every year
months <- array(0, dim = c(NCELL, 12))
ramps <- array(0, dim = c(NCELL, length(access_points_name)))

for(YEAR in 1:n_years_burn_in){
  
  print(YEAR)
  
  temp <- burn_in_effort %>%
    dplyr::select(-c(effort))
  temp <- as.matrix(temp)
  
  for(MONTH in 1:12){
    
    for(RAMP in 1:length(access_points_name)){
      
      for(CELL in 1:NCELL){
        
        ramps[CELL, RAMP] <- BR_U_array[CELL, RAMP, MONTH, 1] * temp[MONTH, RAMP] # we are taking the YEAR = 1 (1900) utility (so that there is no SC or TCs)
      }
    }
    
    months[, MONTH] <-  rowSums(ramps)
  }
  c_burn_in_fishing[ , , YEAR] <- months
  c_burn_in_fishing[ , , YEAR] <- c_burn_in_fishing[ , , YEAR] * water_q[, , YEAR]
  
} # this loop assigns each cell an amount of fishing effort each month and across years based on how useful it is (how far from a boatramp and how frequented that boat ramp is)

head(c_burn_in_fishing)

water$temp <- c_burn_in_fishing[,6,1]
ggplot(water) +
  geom_sf(aes(fill = log(temp)), color = NA) +
  labs(title = paste("Fishing Effort - Year", year_idx, "Month", month_idx),
       fill = "log Effort") +
  scale_fill_gradientn(colours = colour_palette[6:4]) +
  theme_minimal()

saveRDS(c_burn_in_fishing, file = "data/output_data/04_A_commercial_burn_in_fishing.rds")

# Charlotte then adds another burn-in array for high mortality, I'm not doing rn for the sake of moving along

## END ##
