# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Previously simulated dataframes
# Task:    ?
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    August 2026

# -----------------------------------------------------------------------------

# Status: 

# -----------------------------------------------------------------------------

## Load libraries -------------------------------------------------------------

rm(list = ls())

library(tidyverse) # for data manipulation
library(Rcpp) # to execute C++ functions
library(RcppArmadillo) # to execute C++ functions

## Read in functions ----------------------------------------------------------

sourceCpp("functions/age-length_functions/run_full_model_function.cpp", verbose = TRUE)


## Load files  ----------------------------------------------------------------

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))
water <- readRDS("data/output_data/02_watergrid.rds"); plot(water$geometry)


## Get model parameters -------------------------------------------------------

n_yrs_modelled <- 2024-1900 # number of years to simulate

max_cell    <- nrow(water) # Number of cells in the model
max_age     <- 40 # The max age of the fish in the model, see script 05 for correct value
n_lengths   <- length(readRDS("data/output_data/05_length_bins.rds"))
max_year    <- n_yrs_modelled # Number of years the model should run for 

current_pop   <- readRDS("data/output_data/06_burn_in_population.rds") %>% glimpse()
weight        <- readRDS("data/output_data/05_weight.rds") %>% glimpse()
selectivity   <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOw THEY ARE THE SAME FOR ALL FLEETS BUT SHOULD END UP BEING DIFFERENT AT SOME POINT
nat_mort      = 0.12 # from table 4.2, p12 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr

spawn_months <- c(10, 11, 12) # 1-indexed. the function deals with zero-indexing it.
BHa           = as.double(readRDS("data/output_data/05_Beverton-Holt_alpha.rds")) # see script 05
BHb           = as.double(readRDS("data/output_data/05_Beverton-Holt_beta.rds")) # see script 05
PF            = 0.5 # proportion expected to be females
hyperallo     = (1.26 + 1.14 + 1.33)/3 # average from 3 Sparids in Barneche 2018 (1.26, 1.14 and 1.33)
mature         <- readRDS("data/output_data/05_maturity.rds") %>% glimpse()
settlement     <- readRDS("data/output_data/03_B_recruitment.rds") %>% glimpse(); settlement <- settlement[, 1] # selecting a single column because the function expects a vector
age_transition <- readRDS("data/output_data/05_age_transition_matrix.rds") %>% glimpse()
adult_movement <- readRDS("data/output_data/03_b_adult_movement_15_swim_speed.rds") %>% glimpse()
juv_movement   <- readRDS("data/output_data/03_b_juv_movement_10_swim_speed.rds") %>% glimpse()
l_spawn_movement <- readRDS("data/output_data/03_b_large_spawning_movement_15_swim_speed.rds") %>% glimpse()
s_spawn_movement <- readRDS("data/output_data/03_b_small_spawning_movement_15_swim_speed.rds") %>% glimpse()

fleet_names <- c(
  "commercial", 
  "boat_rec", 
  "shore_rec"
)
com_info <- readRDS("data/output_data/04_A_commercial_fishing_info.rds") %>% glimpse()
com_info$fishing_days[com_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

brec_info <- readRDS("data/output_data/04_C_boat_rec_fishing_info.rds") %>% glimpse()
brec_info$fishing_days[brec_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

srec_info <- readRDS("data/output_data/04_B_shore_rec_fishing_info.rds") %>% glimpse()
srec_info$fishing_days[srec_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

fleet_info = list(com_info, 
                  brec_info,
                  srec_info
)


## start the simulation -------------------------------------------------------

# fish outputs to keep
RECONS_pop <- list()
RECONS_SSB <- data.frame(
  year = 1900:2024,
  SSB = NA
)

# fishing outputs to keep
RECONS_effort <- list()
RECONS_catch_weight <- array(0, dim = c(n_yrs_modelled, length(fleet_names))); colnames(RECONS_catch_weight) <- fleet_names
RECONS_catch_length <- array(0, dim = c(n_yrs_modelled, n_lengths, length(fleet_names)))
RECONS_catch_june <- array(0, dim = c(n_yrs_modelled, n_lengths))
RECONS_welfare <- data.frame(
  year = 1900:2024,
  welfare = 0
)

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
                                         juv_movement_prob = juv_movement,
                                         small_spawn_movement_prob = s_spawn_movement,
                                         large_spawn_movement_prob = l_spawn_movement,
                                         fleet_names = fleet_names,
                                         fleet_info = fleet_info
  )
  
  # save the population at the end of the year (to give to the loop again as January population)
  current_pop <- ModelOutput$next_pop #  cell x length x age
  
  # desired outputs are: 
  
  # overall population (spatially, length-structured, over months & years)
  temp_pop_year <- array(0, dim = c(max_cell, n_lengths, 12))
  for (MONTH in 1:12) {
    temp_pop_month <- rowSums(ModelOutput$master_current_pop[[MONTH]], dims = 2) # sum across ages - we don't really care about data across years.
    temp_pop_year[,,MONTH] <- temp_pop_month
    RECONS_pop[[YEAR+1]] <- temp_pop_year
  } # -- population into a year x (cell x lengths x months) array
  
  # spawning biomass (single value per year, the average SSB during spawning months)
  RECONS_SSB[YEAR+1, "SSB"] <- ModelOutput$spawning_biomass # vector of length n_years

  # fishing effort spatial distribution
  RECONS_effort[[YEAR+1]] <- ModelOutput$effort_by_fleet # cell x month x fleet
  
  # annual catch by fleet
  RECONS_catch_weight[YEAR+1,] <- ModelOutput$yearly_catch  # year x fleet of catch sums
  
  # number of fish per length caught
  for(FLEET in 1:3) {
    RECONS_catch_length[YEAR+1,,FLEET] <- colSums(ModelOutput$catch_by_length[,,FLEET])
    } # years x length x fleets
  
  # wadandi catch (number of fish caught by shore rec fishers in June)
  RECONS_catch_june[YEAR+1,]   <- colSums(ModelOutput$catch_by_length_june) # years x lengths
  
  # welfare
  RECONS_welfare[YEAR+1, "welfare"] <- ModelOutput$welfare
}

End = Sys.time()
Runtime = End - Start
Runtime

final_output <- list(population = RECONS_pop,
                     SSB = RECONS_SSB,
                     fishing_effort = RECONS_effort,
                     catch_weight = RECONS_catch_weight,
                     catch_length = RECONS_catch_length,
                     shore_catch_june = RECONS_catch_june,
                     boat_rec_welfare = RECONS_welfare
                     )
saveRDS(final_output, file = "data/output_data/07_A_reconstruction_outputs.rds")

# save the population that will be used in the forecasting
last_population <- ModelOutput$master_current_pop[[12]]
saveRDS(last_population, file = "data/output_data/07_A_last_population.rds")

### END ###
