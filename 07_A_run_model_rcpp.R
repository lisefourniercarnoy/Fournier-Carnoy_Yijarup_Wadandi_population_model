# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Previously simulated dataframes
# Task:    ?
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    June 2026

# -----------------------------------------------------------------------------

# Status: THE BURN-IN NEEDS SOME WORK

# -----------------------------------------------------------------------------

## Load libraries -------------------------------------------------------------

rm(list = ls())

library(tidyverse) # for data manipulation
library(sf) # for dealing with shapefiles
library(terra) # for the bathy layer
library(forcats)
library(RColorBrewer)
# library(MQMF)
library(Rcpp) # to execute C++ functions
library(RcppArmadillo) # to execute C++ functions
library(gmailr) # to send emails via R
library(abind)
# library(beepr)
# library(FishPopPackage)


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
settlement     <- readRDS("data/output_data/03_B_recruitment.rds") %>% glimpse() #; settlement <- settlement[, 1] # selecting a single column because the function expects a vector
age_transition <- readRDS("data/output_data/05_age_transition_matrix.rds") %>% glimpse()
adult_movement <- readRDS("data/output_data/03_b_adult_movement_10_swim_speed.rds") %>% glimpse()

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

# add scaling for each fleet, to account for different gears being less efficient, compared to commercial. see g-sheets Fishing effort reconstruction
brec_info$fishing_days <- brec_info$fishing_days * 0.1
srec_info$fishing_days <- srec_info$fishing_days * 0.017

fleet_info = list(com_info, 
                  brec_info,
                  srec_info
)


## start the simulation -------------------------------------------------------

RECONS_pop <- list() # keep record of the burn-in outputs
RECONS_catch_weight <- list()
RECONS_catch_number <- list()
RECONS_F   <- list()
RECONS_SSB <- list()
RECONS_catch_by_fleet <- list()


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
  fishing_mortality <- ModelOutput$fishing_mortalities
  
  RECONS_pop[[YEAR+1]] <- current_pop # the format is (cell x length x age) x year
  fishing_mortality  <- ModelOutput$fishing_mortalities
  
  RECONS_pop[[YEAR+1]] <- current_pop
  RECONS_F[[YEAR+1]]   <- ModelOutput$annual_F          # cell x length
  RECONS_SSB[[YEAR+1]] <- ModelOutput$spawning_biomass  # cell, summed across spawning months
  fleet_catch_this_year <- sapply(1:length(fleet_names), function(f) {
    sum(ModelOutput$catch_weight_by_fleet[[f]])  # sum over all cells and ages
  })
  RECONS_catch_by_fleet[[YEAR+1]] <- fleet_catch_this_year
  
  # add a plot check, to avoid wasting time on running the function
  total_by_length_age <- apply(current_pop, c(2, 3), sum)
  age_structure       <- colSums(total_by_length_age)   # summed over lengths
  size_structure      <- rowSums(total_by_length_age)   # summed over ages
  
  png(
    filename = file.path("plots/checking_plots_during_setup/07_A_recons_plot_checks/", sprintf("recons_year_%04d.png", YEAR + 1)),
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

End = Sys.time()
Runtime = End - Start
Runtime

### simulated population stability over time ----------------------------------


## a few sanity checks --------------------------------------------------------

### burn-in population stability over time ------------------------------------

total_pop <- lapply(RECONS_pop, function(pop) sum(pop))
par(mfrow = c(1,1))
plot(1:max_year, total_pop, type = "l",
     xlab = "Year", ylab = "Total abundance",
     main = "total Population during reconstruction")
# we are looking for a population that ends up being somewhat stable in the last few years of the burn in, with no huge dips of spikes.


### finite f ------------------------------------------------------------------

# finite fishing mortality (f) is the % of fish removed from the population over any time period.
F_by_length_year <- sapply(RECONS_F, function(f) colMeans(f))  # n_lengths x n_years
matplot(t(F_by_length_year), type = "l", lty = 1,
        xlab = "Year", ylab = "Mean F (across cells)",
        main = "Fishing mortality by length bin over time",
        col = colorRampPalette(c(colour_palette[3], colour_palette[6]))(nrow(F_by_length_year)))

legend("topright", legend = c("small fish", "large fish"), 
       col = c(colour_palette[3], colour_palette[6]), lty = 1, title = "Length bin")
# fishing mortality should be lower for small fish, and vice versa, and higher where there's loads of fishing and vice versa


### age and length structure in the last time step ----------------------------

total_by_length_age <- apply(current_pop, c(2,3), sum)

age_structure <- colSums(total_by_length_age)
size_structure <- rowSums(total_by_length_age)

par(mfrow = c(1, 2))
plot(age_structure, type = "l", xlab = "Age", ylab = "Total abundance",
     main = "Age structure")
plot(size_structure, type = "l", xlab = "Length bin", ylab = "Total abundance",
     main = "Size structure")


### spawning biomass ----------------------------------------------------------

SSB_timeseries <- sapply(RECONS_SSB, sum)  # sum over cells, one value per year
par(mfrow = c(1, 1))
plot(SSB_timeseries, type = "l", lwd = 2, col = colour_palette[6],
     xlab = "Year", ylab = "Total SSB",
     main = "Effective reproductive output")


### fish density --------------------------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds")

ages_to_plot <- c(1, 5, 10, 20, 30, 40)
age_df <- lapply(ages_to_plot, function(a) {
  cell_totals <- apply(RECONS_pop[[60]][,, a], 1, sum)
  water$fish <- cell_totals
  water$age <- factor(paste("Age", a), levels = paste("Age", ages_to_plot))
  water
}) %>% bind_rows()

ggplot(age_df) +
  geom_sf(aes(fill = fish), color = NA) +
  scale_fill_viridis_c(name = "Fish") +
  facet_wrap(~ age, nrow = 1) +
  theme_void() +
  ggtitle("Fish distribution by age (year 60)")


### relative catch of fleets --------------------------------------------------

catch_df <- do.call(rbind, lapply(seq_along(RECONS_catch_by_fleet), function(yr) {
  data.frame(
    year  = 1900 + yr,
    fleet = fleet_names,
    catch_kg = RECONS_catch_by_fleet[[yr]]
  )
}))

ggplot(catch_df, aes(x = year, y = catch_kg, colour = fleet)) +
  geom_line(lwd = 1) +
  labs(x = "Year", y = "Catch (kg)", title = "Annual catch by fleet") +
  theme_bw()

## ### BELOW IS ARCHIVE BUT KEEPING IN CASE THE STRUCTURE IS DIFFERNET FROM BURN IN ####
## Set up initial population --------------------------------------------------

# We need to create loads of objects to save different elements of the model output:
sim_n <- 1 # number of simulations to run - for now it's just 1 but eventually it'll be more

# population information
pop_total       <- array(0, dim = c(max_cell, 12, max_year)) # number of fish of all ages in our population, in each cell (row), each month (column) and each year simulated (matrix slice)
total           <- array(NA, dim = c(max_year, 1)) # for plotting purposes, population summed for each year
pop_total_dist  <- array(0, dim = c(max_cell, max_year)) #???

# catch information
age_catch       <- array(0, dim = c(12, max_age, max_year))
catch_by_cell   <- array(0, dim = c(max_cell, max_year))
catch_by_age    <- array(0, dim = (c(max_age, max_year)))
catch_by_weight <- array(0, dim = (c(max_cell, max_year)))

## Save all information by simulation
sim_pop   <- array(0, dim = c(max_year, sim_n))
sim_ages  <- array(0, dim = c(max_age, max_year, sim_n))

sp_pop_f <- array(0, dim = c(length(shallow_fished_id), max_age, max_year))
sp_pop_ntz <- array(0, dim = c(length(shallow_ntz_id), max_age, max_year))


## Run model simulations ------------------------------------------------------

start = Sys.time()
for (SIM in 1:sim_n){ # Simulation loop - CHARLOTTE HAD 100, I'M STARTING WITH 10

  ## Set up initial population ------------------------------------------------
  pop_total       <- array(0, dim = c(max_cell, 12, max_year)) # number of fish of all ages in our population, in each cell (row), each month (column) and each year simulated (matrix slice)
  total           <- array(NA, dim = c(max_year, 1)) # for plotting purposes, population summed for each year
  
  print(paste0("Simulation number ", SIM)) # progress update
  
  yearly_total <- readRDS("data/output_data/06_RECONS__population.rds")
  
  for (YEAR in 14:(max_year-1)){ # Start of model year loop.
    
    print(paste0("Year ", YEAR))
    
    ### Loop over all the Rcpp functions in the model -------------------------
    
    ModelOutput <- RunModelfunc_cpp(YEAR = YEAR,                                   
                                    MaxCell = max_cell,
                                    MaxYear = max_year, 
                                    
                                    MaxAge = max_age, 
                                    NatMort = nat_mort, 
                                    BHa = BHa, 
                                    BHb = BHb, 
                                    PF = PF, 
                                    AdultMove = adult_movement, 
                                    Mature = maturity,
                                    ha_scaling = hyperallo,
                                    Weight = weight, 
                                    Settlement = settlement, 
                                    
                                    YearlyTotal = yearly_total, 
                                    
                                    Selectivity_com = selectivity_com, 
                                    Selectivity_b_rec = selectivity_b_rec, 
                                    Selectivity_s_rec = selectivity_s_rec, 
                                    
                                    Effort_com = effort_com, # commercial effort of the simulation
                                    Effort_b_rec = effort_b_rec, # boat rec effort of the simulation
                                    Effort_s_rec = effort_s_rec # shore rec effort of the simulation
    )
    print("model success")

    ### Get outputs from the model --------------------------------------------
    
    # Have to add 1 to all YEAR because the loop is now starting at 0
    
    # Abundance in different areas
    pop_total[ , , YEAR+1] <- rowSums(ModelOutput$YearlyTotal[, , 1:max_age], dim = 2) # This flattens the matrix to give you the number of fish present in the population each month in each cell, with layers representing the year

    # Whole area
    water$pop <- pop_total[ , 12, YEAR+1] # We just want the population at the end of the year
    total[YEAR+1, 1] <- sum(water$pop) # Add this to a dataframe we can then use later

    # By zone
    sp_pop_f[, , YEAR+1] <- ModelOutput$YearlyTotal[c(shallow_fished_id), 12, ] # Saving the population at the end of the year in cells <30m depth for plots
    sp_pop_ntz[, , YEAR+1] <- ModelOutput$YearlyTotal[c(shallow_ntz_id), 12, ] # Saving the population at the end of the year in cells <30m depth for plots

    # By cell so we can get distances to boat ramps
    pop_total_dist[ , YEAR+1] <- pop_total[, 12, YEAR+1]

    # Catch data
    monthly_catch               <- ModelOutput$month_catch
    age_catch[,,YEAR+1]         <- colSums(ModelOutput$month_catch) #This is the number of fish in each age class caught in each month

    catch_by_cell[, YEAR+1]     <- rowSums(monthly_catch[, , 3:max_age], dims = 1) # Number of legal size fish caught in each cell (age 3+)
    catch_by_age[, YEAR+1]      <- colSums(age_catch[, , YEAR+1]) # number of fish caught at by the end of the year in each age class

    monthly_catch_weight        <- ModelOutput$month_catch_weight
    catch_by_weight[ , YEAR+1]  <- rowSums(monthly_catch_weight[,,3:max_age], dims = 1)

    sim_pop[YEAR+1, SIM]        <- total[YEAR+1, 1]
    sim_ages[ , YEAR+1, SIM]    <- colSums(ModelOutput$YearlyTotal[, 12, 1:max_age]) # number of fish present in age age group at the end of the year
  } # this loop runs the population model for every year, obtaining a simulated population over 20 or 30 years 

  ## Population in different zones
  SIM_sp_f[[SIM]] <- sp_pop_f
  SIM_sp_ntz[[SIM]] <- sp_pop_ntz
  SIM_n_dist[[SIM]] <- pop_total_dist

  ## Catches
  SIM_n_catches[[SIM]] <- catch_by_cell # Catches in each cell
  SIM_age_catches[[SIM]] <- catch_by_age # Catches by age in each month of the year in each year
  SIM_weight_catches[[SIM]] <- catch_by_weight

  ## Save the last simulation? ------------------------------------------------

  if(SIM == sim_n){ # Saving if statement
    #print(Movement_Speed)
    print(scenario)

    ### Population ------------------------------------------------------------
    
    # Total population
    saveRDS(sim_pop, file = paste0("simulations/dummy_run/", scenario, "_total_population_try1.rds"))

    # Numbers of each age that make it to the end of each year
    saveRDS(sim_ages, file = paste0("simulations/dummy_run/", scenario, "_age_distribution_try1.rds"))
    saveRDS(SIM_sp_ntz, file = paste0("simulations/dummy_run/", scenario, "_sp_population_ntz_try1.rds")) # Numbers of fish of each age, inside sanctuary zones
    saveRDS(SIM_sp_f, file = paste0("simulations/dummy_run/", scenario, "_sp_population_fished_try1.rds")) # Numbers of fish of each age, outside sanctuary zones
    saveRDS(SIM_n_dist, file = paste0("simulations/dummy_run/", scenario, "_cell_population_try1.rds")) # Number of fish in each cell at the end of each year

    ### Catches ---------------------------------------------------------------

    saveRDS(SIM_age_catches, file = paste0("simulations/dummy_run/", scenario, "_catch_by_age_baranov_try1.rds")) # Catch in each year by age
    saveRDS(SIM_n_catches, file = paste0("simulations/dummy_run/", scenario, "_catch_by_cell_baranov_try1.rds")) # catch in each cell across the year
    saveRDS(SIM_weight_catches, file = paste0("simulations/dummy_run/", scenario, "_catch_by_weight_try1.rds")) # Catch in each cell by weight across the year

  } else { } # End saving if statement

} # this loop makes new simulations using the selected scenario

end = Sys.time() 
runtime = end - start
runtime

finished_email <- gm_mime() %>%
  gm_to("lise.fourniercarnoy@research.uwa.edu.au") %>%
  gm_from("lise.fourniercarnoy@marineecology.io") %>%
  gm_subject("Model code is done running") %>%
  gm_text_body(paste("Feckin finally! This run took ", runtime, "minutes."))

d <- gm_create_draft(finished_email)
gm_send_draft(d)

### END ###
