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
library(sf) # for dealing with shapefiles
library(terra) # for the bathy layer
library(forcats)
library(RColorBrewer)
library(Rcpp) # to execute C++ functions
library(RcppArmadillo) # to execute C++ functions
library(gmailr) # to send emails via R
library(abind) # to del with cubes 


## Set up to send emails when finished running --------------------------------

# This only needs to be done once.
# rappdirs::user_data_dir("gmailr")
# path_old <- "data/backend_files/client_secret_75883251379-aboc07p3og08ohtu8ea030l87j6fnb96.apps.googleusercontent.com.json"
# d <- fs::dir_create(rappdirs::user_data_dir("gmailr"), recurse=TRUE) # create directory for the gmailr app
# fs::file_move(path_old, d) # place JSON file where gmailr expects it
# # client secret: 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 (for use in console.cloud.google.com)
# gm_auth_configure()
# gm_oauth_client()
# gm_auth() # R is now able to compose and send emails from address: lise.fourniercarnoy@marineecology.io


## Read in functions ----------------------------------------------------------

sourceCpp("functions/age-length_functions/run_full_model_function.cpp", verbose = TRUE)

## Load files  ----------------------------------------------------------------

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))
water <- readRDS("data/output_data/02_watergrid.rds"); plot(water$geometry)

current_scenario <- "S00_business_as_usual"


## Get model parameters -------------------------------------------------------

n_yrs_modelled <- 2050-2024 # number of years to simulate

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
com_info <- readRDS(paste0("data/output_data/08_", current_scenario, "/08_commercial_fishing_info.rds")) %>% glimpse()
com_info$fishing_days[com_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

brec_info <- readRDS(paste0("data/output_data/08_", current_scenario, "/08_boat_rec_fishing_info.rds")) %>% glimpse()
brec_info$fishing_days[brec_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

srec_info <- readRDS(paste0("data/output_data/08_", current_scenario, "/08_shore_rec_fishing_info.rds")) %>% glimpse()
srec_info$fishing_days[srec_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

fleet_info = list(com_info, 
                  brec_info,
                  srec_info
)

## create a function for all your simulations ---------------------------------

run_one_simulation <- function(SIM,
                               
                               YEAR,
                               max_cell,
                               max_age,
                               max_year,
                               n_lengths,
                               initial_pop,
                               weight,
                               selectivity,
                               age_transition,
                               nat_mort,
                               spawn_months,
                               BHa,
                               BHb,
                               PF,
                               hyperallo,
                               mature,
                               settlement,
                               adult_movement,
                               juv_movement,
                               s_spawn_movement,
                               l_spawn_movement,
                               fleet_names,
                               fleet_info
                               
                               ) {
  
  initial_pop <- current_pop
  
  SCENARIO_pop <- list() # keep record of the burn-in outputs
  months_to_save <- c(Jan = 1, Feb = 2, Mar = 3, # which months to  sve for plot checks
                      Apr = 4, May = 5, Jun = 6, 
                      Jul = 7, Aug = 8, Sep = 9, 
                      Oct = 10, Nov = 11, Dec = 12)
  
  SCENARIO_effort <- list() # cell x month x fleet, one array per year
  SCENARIO_F <- list() # one 12 x n_fleets matrix per year
  SCENARIO_catch_weight <- array(0, dim = c(n_yrs_modelled, length(fleet_names))); colnames(SCENARIO_catch_weight) <- fleet_names
  SCENARIO_SSB <- list()
  
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
    
    # fill objects to check after the burn-in
    SCENARIO_pop[[YEAR+1]]         <- setNames(ModelOutput$master_current_pop[months_to_save], names(months_to_save))
    SCENARIO_effort[[YEAR+1]]      <- ModelOutput$effort_by_fleet # array: cell x month x fleet
    SCENARIO_F[[YEAR+1]]           <- ModelOutput$fishing_mortality
    SCENARIO_catch_weight[YEAR+1,] <- ModelOutput$yearly_catch
    SCENARIO_SSB[[YEAR+1]]         <- ModelOutput$spawning_biomass  # cell, summed across spawning months
    
    list(
      SIM = SIM,
      SCENARIO_pop = SCENARIO_pop,
      SCENARIO_effort = SCENARIO_effort,
      SCENARIO_F = SCENARIO_F,
      SCENARIO_catch_weight = SCENARIO_catch_weight,
      SCENARIO_SSB = SCENARIO_SSB
    )
    
    
  }
}


## start the simulation -------------------------------------------------------

n_cores <- 4
n_sim <- 4

library(parallel)
library(doParallel)

cl <- parallel::makeCluster(n_cores, outfile = file.path("data/output_data/09_test/cluster_log.txt"))
doParallel::registerDoParallel(cl)

fun_dir_cpp <- file.path("functions/age-length_functions/run_full_model_function.cpp")

# compile the functions for every core that will use it.
parallel::clusterCall(cl, function(path) {
  library(Rcpp)
  library(RcppArmadillo)
  Rcpp::sourceCpp(path)
  exists("run_full_model_function")
}, path = fun_dir_cpp)

# # Task execution: only calls the already-compiled function
all_sim_results <- foreach(
  SIM = 1:n_sim,
  .export = c("initial_pop", 
              "max_cell", 
              "max_age", 
              "max_year", 
              "n_lengths",
              "weight", 
              "selectivity", 
              "age_transition", 
              "nat_mort", 
              "spawn_months",
              "BHa", 
              "BHb", 
              "PF", 
              "hyperallo", 
              "mature",
              "settlement",
              "adult_movement",
              "juv_movement",
              "s_spawn_movement",
              "l_spawn_movement",
              "fleet_names", 
              "fleet_info", 
              "run_one_simulation"
              )
) %dopar% {
  
  run_one_simulation(
    SIM = SIM,
    initial_pop = initial_pop,
    max_cell = max_cell, 
    max_age = max_age, 
    max_year = max_year, 
    n_lengths = n_lengths,
    weight = weight, 
    selectivity = selectivity, 
    age_transition = age_transition,
    nat_mort = nat_mort, 
    spawn_months = spawn_months,
    BHa = BHa, 
    BHb = BHb, 
    PF = PF, 
    hyperallo = hyperallo,
    mature = mature, 
    settlement = settlement, 
    adult_movement_prob = adult_movement,
    juv_movement_prob = juv_movement,
    small_spawn_movement_prob = s_spawn_movement,
    large_spawn_movement_prob = l_spawn_movement,
    fleet_names = fleet_names, 
    fleet_info = fleet_info
    
  )
}

parallel::stopCluster(cl)

sapply(all_sim_results, function(x) sum(sapply(x$RECONS_pop, is.null)))

identical(all_sim_results[[1]]$RECONS_SSB, all_sim_results[[2]]$RECONS_SSB)
















SCENARIO_pop <- list() # keep record of the burn-in outputs
months_to_save <- c(Jan = 1, Feb = 2, Mar = 3, # which months to  sve for plot checks
                    Apr = 4, May = 5, Jun = 6, 
                    Jul = 7, Aug = 8, Sep = 9, 
                    Oct = 10, Nov = 11, Dec = 12)

SCENARIO_effort <- list() # cell x month x fleet, one array per year
SCENARIO_F <- list() # one 12 x n_fleets matrix per year
SCENARIO_catch_weight <- array(0, dim = c(n_yrs_modelled, length(fleet_names))); colnames(SCENARIO_catch_weight) <- fleet_names
SCENARIO_SSB <- list()

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
  
  
  # fill objects to check after the burn-in
  SCENARIO_pop[[YEAR+1]]         <- setNames(ModelOutput$master_current_pop[months_to_save], names(months_to_save))
  SCENARIO_effort[[YEAR+1]]      <- ModelOutput$effort_by_fleet # array: cell x month x fleet
  SCENARIO_F[[YEAR+1]]           <- ModelOutput$fishing_mortality
  SCENARIO_catch_weight[YEAR+1,] <- ModelOutput$yearly_catch
  SCENARIO_SSB[[YEAR+1]]         <- ModelOutput$spawning_biomass  # cell, summed across spawning months
  
  # # add a plot check, to avoid wasting time on running the function
  # total_by_length_age <- apply(ModelOutput$master_current_pop[[12]], c(2, 3), sum)
  # age_structure       <- colSums(total_by_length_age)   # summed over lengths
  # size_structure       <- rowSums(total_by_length_age)   # summed over ages
  # 
  # png(
  #   filename = file.path("plots/script_plot_checks/06/06_age_length_structures/", sprintf("SCENARIO_year_%04d.png", YEAR + 1)),
  #   width = 1200, height = 500, res = 120
  # )
  # 
  # par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  # 
  # plot(age_structure,
  #      type = "l", lwd = 2, col = colour_palette[6],
  #      xlab = "Age class", ylab = "Total abundance",
  #      main = paste0("Age structure — Year ", YEAR + 1))
  # 
  # plot(size_structure,
  #      type = "l", lwd = 2, col = colour_palette[4],
  #      xlab = "Length bin", ylab = "Total abundance",
  #      main = paste0("Size structure — Year ", YEAR + 1))
  # 
  # dev.off()
}

End = Sys.time()
Runtime = End - Start
Runtime


### CHECK: population trends over time ----------------------------------------

total_pop <- sapply(SCENARIO_pop, function(year_list) {
  sum(sapply(year_list, sum))
})

plot(1:max_year, total_pop, type = "l",
     xlab = "Year", ylab = "Total abundance",
     main = "Total Population during SCENARIOtruction")


### CHECK: age and length structure in the last time step ---------------------

total_by_length_age <- apply(current_pop, c(2,3), sum)
age_structure <- colSums(total_by_length_age)
size_structure <- rowSums(total_by_length_age)

par(mfrow = c(1, 2))
plot(age_structure, type = "l", xlab = "Age", ylab = "Total abundance",
     main = "Age structure")
plot(size_structure, type = "l", xlab = "Length bin", ylab = "Total abundance",
     main = "Size structure")


### CHECK: spawning biomass ---------------------------------------------------

SSB0 <- readRDS("data/output_data/06_burn_in_SSB0.rds")
SSB_timeseries <- sapply(SCENARIO_SSB, sum)  # sum over cells, one value per year
SSB_ss <- read.csv("data/input_data/digitised_plots_for_checking/SSB_digitised_from_stock_assessment.csv") |> 
  arrange(year)

B_rel <- SSB_timeseries / SSB0
par(mfrow = c(1, 1))

# below, the red line should align on the blue as best as possible.
plot(x = SSB_ss$year, SSB_ss$SSB, lwd = 2, col = "steelblue", type = "l", ylim = c(0, 1),
     xlab = "Year", ylab = "SSB relative to 1900",
     main = "Effective reproductive output")
lines(x = 1975:2025, y = B_rel[75:125], col = "firebrick", lwd = 2)
# here we are checking whether the relative SSB (this year's SSB/SSB0) matches the stock assessment.

### fish density --------------------------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds")
ages_to_plot <- c(1, 2)
months_to_plot <- 1:12
year_to_plot <- 120

age_month_df <- lapply(months_to_plot, function(m) {
  lapply(ages_to_plot, function(a) {
    cell_totals <- apply(SCENARIO_pop[[year_to_plot]][[m]][,, a], 1, sum)
    water$fish <- cell_totals
    water$age <- factor(paste("Age", a), levels = paste("Age", ages_to_plot))
    water$month <- factor(month.abb[m], levels = month.abb)
    water
  }) %>% bind_rows()
}) %>% bind_rows()

ggplot(age_month_df) +
  geom_sf(aes(fill = fish), color = NA) +
  scale_fill_viridis_c(name = "Fish") +
  facet_grid(age ~ month) +
  theme_void() +
  ggtitle(paste0("Fish distribution by age and month (year ", year_to_plot, ")"))


### relative catch of fleets --------------------------------------------------

yearly_catch <- as.data.frame(SCENARIO_catch_weight)
names(yearly_catch) <- fleet_names
ggplot(yearly_catch[75:124,]) +
  geom_line(lwd = 1, aes(x = 75:124, y = commercial), colour = "red") +
  geom_line(lwd = 1, aes(x = 75:124, y = boat_rec), colour = "blue") +
  geom_line(lwd = 1, aes(x = 75:124, y = shore_rec), colour = "green") +
  labs(x = "Year", y = "Catch (kg)", title = "Annual catch by fleet") +
  theme_bw()


### CHECK: effort distribution ------------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds")
year_to_plot = 120
water_effort <- water %>%
  mutate(
    commercial = SCENARIO_effort[[year_to_plot]][, 12, 1],
    boat_rec   = SCENARIO_effort[[year_to_plot]][, 12, 2],
    shore_rec  = SCENARIO_effort[[year_to_plot]][, 12, 3]
  ) %>%
  pivot_longer(
    cols = c(commercial, boat_rec, shore_rec),
    names_to = "fleet",
    values_to = "effort"
  )

ggplot(water_effort) +
  geom_sf(aes(fill = effort), color = NA) +
  scale_fill_viridis_c(name = "Effort") +
  facet_wrap(~ fleet, nrow = 1) +
  theme_void() +
  ggtitle("Fishing effort by fleet (ultimate year, Month 12)")


plot_data <- water_effort[water_effort$fleet == "commercial", c("effort")] # keeps geometry automatically

mapview::mapview(plot_data,
                 zcol = "effort",
                 #col.regions = viridisLite::viridis(100),
                 layer.name = "Effort")


### END ###
