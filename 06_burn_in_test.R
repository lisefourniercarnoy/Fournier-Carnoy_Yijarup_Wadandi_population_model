# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    ?
# Task:    Set up a population out of the population parameters and fishing effort
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    April 2025

# -----------------------------------------------------------------------------

# Status:

# -----------------------------------------------------------------------------

rm(list = ls())

# Load libraries
library(tidyverse) # for data manipulation
library(sf) # for spatial objects
library(MQMF)
library(Rcpp) # for reading the C++ functions in R
library(RcppArmadillo) # for reading the C++ functions in R
library(abind) # for dealing with arrays i think.

## Read in the functions ------------------------------------------------------

sourceCpp("functions/age-length_functions/run_full_model_function.cpp", verbose = TRUE)
colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))


## Read files -----------------------------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds"); plot(water$geometry)

## Get model parameters -------------------------------------------------------

n_yrs_modelled <- 60 # number of burn-in years - at least a fish's full life.

max_cell    <- nrow(water) # Number of cells in the model
max_age     <- 40 # The max age of the fish in the model, see script 05 for correct value
n_lengths   <- length(readRDS("data/output_data/05_length_bins.rds"))
max_year    <- n_yrs_modelled # Number of years the model should run for 

starting_pop  <- readRDS("data/output_data/05_starting_population.rds") %>% glimpse()
weight        <- readRDS("data/output_data/05_weight.rds") %>% glimpse()
selectivity   <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOw THEY ARE THE SAME FOR ALL FLEETS BUT SHOULD END UP BEING DIFFERENT AT SOME POINT
nat_mort      = 0.12 # from table 4.2, p12 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr

spawn_months <- c(10, 11, 12) # 1-indexed. the function deals with zero-indexing it.
BHa           = as.double(readRDS("data/output_data/05_Beverton-Holt_alpha.rds")) # see script 05
BHb           = as.double(readRDS("data/output_data/05_Beverton-Holt_beta.rds")) # see script 05
PF            = 0.5 # proportion expected to be females
hyperallo     = (1.26 + 1.14 + 1.33)/3 # average from 3 Sparids in Barneche 2018 (1.26, 1.14 and 1.33)
mature         <- readRDS("data/output_data/05_maturity.rds") %>% glimpse()
settlement     <- readRDS("data/output_data/03_B_recruitment.rds") %>% glimpse() #; settlement <- settlement[, 1] # selecting a single column because the function expects a vector
age_transition <- readRDS("data/output_data/05_age_transition_matrix.rds") %>% glimpse()
adult_movement <- readRDS("data/output_data/03_b_adult_movement_10_swim_speed.rds") %>% glimpse()

fleet_names <- c(
  "commercial", 
  "boat_rec", 
  "shore_rec")

com_info <- readRDS("data/output_data/04_A_commercial_fishing_info.rds") %>% glimpse()
com_info$fishing_days[com_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

brec_info <- readRDS("data/output_data/04_C_boat_rec_fishing_info.rds") %>% glimpse()
brec_info$fishing_days[brec_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

srec_info <- readRDS("data/output_data/04_B_shore_rec_fishing_info.rds") %>% glimpse()
srec_info$fishing_days[srec_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.


# for the burn-in, we'll use a constant low level of fishing, so replace all years with 1900 fishing effort, catchability, etc.
com_info$fishing_days[,,1:dim(com_info$fishing_days)[[3]]] <- com_info$fishing_days[,,1]
brec_info$fishing_days[,,1:dim(brec_info$fishing_days)[[3]]] <- brec_info$fishing_days[,,1]
srec_info$fishing_days[,,1:dim(srec_info$fishing_days)[[3]]] <- srec_info$fishing_days[,,1]

com_info$catchability[,,1:dim(com_info$catchability)[[3]]] <- com_info$catchability[,,1]
brec_info$catchability[,,1:dim(brec_info$catchability)[[3]]] <- brec_info$catchability[,,1]
srec_info$catchability[,,1:dim(srec_info$catchability)[[3]]] <- srec_info$catchability[,,1]

first_vals <- com_info$attractivity[[1]][, , 1]  # 1585 x 16 matrix (month 1, year 1)
com_info$attractivity <- lapply(com_info$attractivity, function(x) {
  array(first_vals, dim = dim(x))
})

first_vals <- brec_info$attractivity[[1]][, , 1]  # 1585 x 16 matrix (month 1, year 1)
brec_info$attractivity <- lapply(brec_info$attractivity, function(x) {
  array(first_vals, dim = dim(x))
})

first_vals <- srec_info$attractivity[[1]][, , 1]  # 1585 x 16 matrix (month 1, year 1)
srec_info$attractivity <- lapply(srec_info$attractivity, function(x) {
  array(first_vals, dim = dim(x))
})

fleet_info = list(com_info, 
                  brec_info,
                  srec_info
                  )


## Set up the initial population ----------------------------------------------

total <- array(0, dim = c(max_year, 1))

current_pop <- array(0, dim = c(max_cell, n_lengths, max_age)) # for every cell (row), and every length (column) across all fish ages (matrix slice), we will have a population

BURN_IN_pop <- list() # keep record of the burn-in outputs
BURN_IN_SSB <- list()
BURN_IN_catch_weight <- array(0, dim = c(n_yrs_modelled, length(fleet_names))); names(BURN_IN_catch_weight) <- fleet_names
BURN_IN_effort <- list() # cell x month x fleet, one array per year

for(AGE in 1:max_age){
  total_this_age <- starting_pop[AGE, ]
  # distribute proportionally to settlement — same habitat weighting as recruits
  settlement_prop <- settlement / sum(settlement) # length max_cell
  
  # outer product: each cell's proportion × each length-class abundance
  current_pop[, , AGE] <- outer(settlement_prop, total_this_age)
}

cat("Total fish initialised:", sum(current_pop), "\n")


## Start the burn-in ----------------------------------------------------------

Start = Sys.time()
for (YEAR in 0:(max_year-1)){ # max_year-1 because it starts at 0.
  
  # loop over all the Rcpp functions in the model
  ModelOutput <- run_full_model_function(YEAR = YEAR,
                                         max_cell = max_cell,
                                         max_age = max_age,
                                         max_year = max_year,
                                         n_lengths = n_lengths,
                                         current_pop = current_pop,
                                         weight = weight,
                                         selectivity = selectivity,
                                         age_transition = age_transition,
                                         natural_mortality = nat_mort,
                                         spawning_months = spawn_months,
                                         BHa = BHa,
                                         BHb = BHb,
                                         PF = PF,
                                         ha_scaling = hyperallo,
                                         maturity = mature,
                                         settlement = settlement,
                                         adult_movement_prob = adult_movement,
                                         fleet_names = fleet_names,
                                         fleet_info = fleet_info
  )

  # save the population at the end of the year (to give to the loop again as January population)
  current_pop <- ModelOutput$next_pop #  cell x length x age

  # fill objects to check after the burn-in
  BURN_IN_pop[[YEAR+1]]         <- current_pop # the format is (cell x length x age) x year
  BURN_IN_SSB[[YEAR+1]]         <- ModelOutput$spawning_biomass  # cell, summed across spawning months
  BURN_IN_catch_weight[YEAR+1,] <- sapply(ModelOutput$catch_weight_by_fleet, sum)
  BURN_IN_effort[[YEAR+1]] <- ModelOutput$effort_by_fleet # array: cell x month x fleet
  
  # add a plot check, to avoid wasting time on running the function
  total_by_length_age <- apply(current_pop, c(2, 3), sum)
  age_structure       <- colSums(total_by_length_age)   # summed over lengths
  size_structure      <- rowSums(total_by_length_age)   # summed over ages
  
  png(
    filename = file.path("plots/checking_plots_during_setup/06_burn_in_plot_checks/", sprintf("burnin_year_%04d.png", YEAR + 1)),
    width = 1200, height = 500, res = 120
  )
  
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  plot(age_structure,
       type = "l", lwd = 2, col = colour_palette[6],
       xlab = "Age class", ylab = "Total abundance",
       main = paste0("Age structure — Year ", YEAR + 1))
  
  plot(size_structure,
       type = "l", lwd = 2, col = colour_palette[4],
       xlab = "Length bin", ylab = "Total abundance",
       main = paste0("Size structure — Year ", YEAR + 1))
  
  dev.off()
  
}

End = Sys.time(); Runtime = End - Start; Runtime
# after the burn-in, make sure you run the sections below, where some outputs are saved for the historical reconstruction.

### CHECK: burn-in population stability over time -----------------------------

total_pop <- lapply(BURN_IN_pop, function(pop) sum(pop))
par(mfrow = c(1,1))
plot(1:max_year, total_pop, type = "l",
     xlab = "Year", ylab = "Total abundance",
     main = "total Population during burn-in")
# we are looking for a population that ends up being somewhat stable in the last few years of the burn in, with no huge dips nor spikes.

saveRDS(current_pop, file = "data/output_data/06_burn_in_population.rds")


### CHECK: age and length structure in the last time step ---------------------

total_by_length_age <- apply(current_pop, c(2,3), sum)
age_structure <- colSums(total_by_length_age)
size_structure <- rowSums(total_by_length_age)

par(mfrow = c(1, 2))
plot(age_structure, type = "l", xlab = "Age", ylab = "Total abundance",
     main = "Age structure")
plot(size_structure, type = "l", xlab = "Length bin", ylab = "Total abundance",
     main = "Size structure")
# we're looking for a nice inverse-exponential age structure.
# length structure has many small fish, many large fish (these are 20-40yo fish that have accumulated in the max length for many years)


### CHECK: spawning biomass ---------------------------------------------------

SSB_timeseries <- sapply(BURN_IN_SSB, sum)  # sum over cells, one value per year
par(mfrow = c(1, 1))
plot(SSB_timeseries, type = "l", lwd = 2, col = colour_palette[6],
     xlab = "Year", ylab = "Total SSB",
     main = "Spawning Stock Biomass over time")
# as with total population, the SSB should plateau.

# we'll also export the final year's SSB to compare the historical reconstruction period's SSB to:
SSB0 <- SSB_timeseries[length(SSB_timeseries)]
saveRDS(SSB0, file = "data/output_data/06_burn_in_SSB0.rds")


### CHECK: fish density -------------------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds")

ages_to_plot <- c(1, 5, 10, 20, 30, 40)
age_df <- lapply(ages_to_plot, function(a) {
  cell_totals <- apply(BURN_IN_pop[[60]][,, a], 1, sum)
  water$fish <- cell_totals
  water$age <- factor(paste("Age", a), levels = paste("Age", ages_to_plot))
  water
}) %>% 
  bind_rows()

ggplot(age_df) +
  geom_sf(aes(fill = fish), color = NA) +
  scale_fill_viridis_c(name = "Fish") +
  facet_wrap(~ age, nrow = 1) +
  theme_void() +
  ggtitle("Fish distribution by age (year 60)")
# here we're checking whether the spatial distribution of fish is reasonable, whether movement makes sense, and whether the number of fish per cell is also realistic. 

### CHECK: effort distribution ------------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds")

water_effort <- water %>%
  mutate(
    commercial = BURN_IN_effort[[60]][, 12, 1],
    boat_rec   = BURN_IN_effort[[60]][, 12, 2],
    shore_rec  = BURN_IN_effort[[60]][, 12, 3]
  ) %>%
  pivot_longer(
    cols = c(commercial, boat_rec, shore_rec),
    names_to = "fleet",
    values_to = "effort"
  )

ggplot(water_effort) +
  geom_sf(aes(fill = (effort)), color = NA) +
  scale_fill_viridis_c(name = "Effort") +
  facet_wrap(~ fleet, nrow = 1) +
  theme_void() +
  ggtitle("Fishing effort by fleet (Year 60, Month 12)")


plot_data <- water_effort[water_effort$fleet == "commercial", c("effort")] # keeps geometry automatically

mapview::mapview(plot_data,
        zcol = "effort",
        #col.regions = viridisLite::viridis(100),
        layer.name = "Effort")


### CHECK: relative catch of fleets -------------------------------------------

yearly_catch <- as.data.frame(BURN_IN_catch_weight)
names(yearly_catch) <- fleet_names
ggplot(yearly_catch) +
  geom_line(lwd = 1, aes(x = 1:n_yrs_modelled, y = commercial), colour = "red") +
  geom_line(lwd = 1, aes(x = 1:n_yrs_modelled, y = boat_rec), colour = "blue") +
  geom_line(lwd = 1, aes(x = 1:n_yrs_modelled, y = shore_rec), colour = "green") +
  labs(x = "Year", y = "Catch (kg)", title = "Annual catch by fleet") +
  theme_bw()

## END ##
