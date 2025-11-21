# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Previously simulated dataframes
# Task:    ?
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    May 2025

# -----------------------------------------------------------------------------

# Status: Starting out... NOT ADAPTED VERY MUCH AT ALL

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

sourceCpp("functions/C_Model_RcppArm_test.cpp")
source("functions/X_Functions.R")

## Load files  ----------------------------------------------------------------

movement_Speed = "fast_movement"
colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))

#### LOAD FILES ####

adult_movement    <- readRDS("data/output_data/03_adult_movement_10_swim_speed.rds") %>% glimpse()
settlement        <- readRDS("data/output_data/03_recruitment.rds") %>% glimpse()
selectivity_com   <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOw THEY ARE THE SAME BUT SHOULD END UP BEING DIFFERENT AT SOME POINT
selectivity_b_rec <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOW THEY ARE THE SAME BUT SHOULD END UP BEING DIFFERENT AT SOME POINT
selectivity_s_rec <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOW THEY ARE THE SAME BUT SHOULD END UP BEING DIFFERENT AT SOME POINT
maturity          <- readRDS("data/output_data/05_maturity.rds") %>% glimpse()
weight            <- readRDS("data/output_data/05_weight.rds") %>% glimpse()

water             <- readRDS("data/output_data/02_watergrid.rds") %>% st_make_valid() %>% glimpse(); plot(water)

## make selectivity an array, as for burn-in
selectivity_com <- selectivity_com[, , 44] # selecting the most recent selectivity-retention
for(i in 1:6){
  selectivity_com <- abind(selectivity_com, selectivity_com, along=3)
}
selectivity_b_rec <- selectivity_b_rec[, , 44] # selecting the most recent selectivity-retention
for(i in 1:6){
  selectivity_b_rec <- abind(selectivity_b_rec, selectivity_b_rec, along=3)
}
selectivity_s_rec <- selectivity_s_rec[, , 44] # selecting the most recent selectivity-retention
for(i in 1:6){
  selectivity_s_rec <- abind(selectivity_s_rec, selectivity_s_rec, along=3)
}


# crop and mask the bathymetry raster to the water extent - takes a while, re-run if extent changes after 18/09/2025
# bathy <- rast("data/input_data/AusBathyTopo__Australia__2024_250m_MSL_cog.tif") %>% glimpse()
# bathy <- project(bathy, st_crs(water)$wkt)
# water_vect <- vect(water); bathy_cropped <- crop(bathy, water_vect)
# bathy_wadandi  <- mask(bathy_cropped, water_vect)
# saveRDS(bathy_wadandi, "data/output_data/07_wadandi_bathymetry.rds")
bathy <- readRDS("data/output_data/07_wadandi_bathymetry.rds")
ntz_list <- readRDS("data/output_data/02_no_take_list.rds")

## Setup scenario names -------------------------------------------------------

# for each scenario, we need to get the right burn-in files
burn_in_pop       <- readRDS(paste0("data/output_data/06_burn_in_population.rds")) %>% glimpse() # this is the same for all scenarios.

#scenario <- "07_A_SC" # 07_A for the script we're in, and SC for the scenario we're testing
#scenario <- "07_A_S00" # 07_A for the script we're in, and S00 for the scenario we're testing
scenario <- "07_A_S01" # 07_A for the script we're in, and S01 for the scenario we're testing

scenario_root_file <- case_when(scenario == "07_A_SC" ~ "04",
                                scenario == "07_A_S00" ~ "S00",
                                scenario == "07_A_S01" ~ "S01"
)

if (scenario == "07_A_S01") {
  effort_com <- array(0, dim = c(1063, 12, 80))
} else {
  effort_com <- readRDS(paste0(
    "data/output_data/",
    scenario_root_file,
    if (scenario == "07_A_SC") "A",
    "_commercial_burn_in_fishing.rds"
  ))
}
effort_s_rec      <- readRDS(paste0("data/output_data/", scenario_root_file, if(scenario == "07_A_SC") {"B"}, "_shore_rec_burn_in_fishing.rds")) %>% glimpse()
if (scenario == "07_A_S01") {
  effort_b_rec <- array(0, dim = c(1063, 12, 80))
} else {
  effort_b_rec <- readRDS(paste0(
    "data/output_data/",
    scenario_root_file,
    if (scenario == "07_A_SC") "C",
    "_boat_rec_burn_in_fishing.rds"
  ))
}



## Crop water ??? -------------------------------------------------------------
#  FOR SOME REASON SOME NEARSHORE CELLS ARE NOT INCLUDED. NEEDS FIXING SOMEHOW.
# We'll select the cells that are in <30m water (this will allow us to plot catch in different depths??)

# only done once, valid for all simulations

# plot(bathy)
# plot(water$geometry, add = T)
# 
# water_points <- st_centroid_within_poly(water); plot(water_points) # make centroids
# water_bathy <- raster::extract(bathy, water_points, fun = mean, df = TRUE); plot(water_bathy) # get the depth at the centroid location
# water_bathy <- water_bathy %>% mutate(ID = as.factor(ID))
# model_WHA <- water %>% st_intersects(., water) %>% as.data.frame()
# water_WHA <- water[c(as.numeric(model_WHA$row.id)), ]
# 
# water_shallow <- water_WHA %>%
#   mutate(ID = as.factor(ID)) %>%
#   left_join(., water_bathy, by="ID") %>%
#   rename(bathy = "AusBathyTopo__Australia__2024_250m_MSL_cog") %>%
#   filter(bathy >= c(-30)) %>% 
#   filter(!is.na(bathy))
# summary(is.na(water_shallow$bathy))
# plot(water_shallow$geometry)
# 
# shallow_cells_ntz <- water_shallow %>%
#   filter(ID %in% ntz_list[[1]]) %>%
#   distinct(ID, .keep_all = TRUE)
# 
# shallow_cells_fished <- water_shallow %>%
#   filter(!ID %in% ntz_list[[1]]) %>%
#   distinct(ID, .keep_all = TRUE)
# 
# plot(shallow_cells_fished$geometry, col = colour_palette[4]) # check that they make sense
# plot(shallow_cells_ntz$geometry, col = colour_palette[6], add = T)
# 
# shallow_ntz_id <- as.numeric(levels(shallow_cells_ntz$ID))[as.integer(shallow_cells_ntz$ID)]
# shallow_fished_id <- as.numeric(levels(shallow_cells_fished$ID))[as.integer(shallow_cells_fished$ID)]
# 
# saveRDS(shallow_ntz_id, "data/output_data/07A_shallow_ntz_cell_id.rds")
# saveRDS(shallow_fished_id, "data/output_data/07A_shallow_fished_cell_id.rds")

shallow_ntz_id <- readRDS("data/output_data/07A_shallow_ntz_cell_id.rds")
shallow_fished_id <- readRDS("data/output_data/07A_shallow_fished_cell_id.rds")

## Set model parameters -------------------------------------------------------

# Natural Mortality
nat_mort <- 0.12 # table 9, p.95 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1206&context=fr_rr
hyperallo <- 1.24 # average from 3 Sparids in Barneche 2018 (1.26, 1.14 and 1.33)


# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
BHa = 0.4344209 # NOT CHANGED FROM CHARLOTTE, NEED TO FIND SOMEWHERE
BHb = 0.0002349538 # NOT CHANGED FROM CHARLOTTE, NEED TO FIND SOMEWHERE
PF = 0.5 # proportion expected to be females

# Model settings
max_cell    <- nrow(water) # Number of cells in the model
max_age     <- 30 # the max age of the fish in the model (-1 to account for the fact that Rcpp functions start from 0)
max_year    <- 49-1 # number of years the model should run for (1945 + number of burn-in years = 1975, 2024-1975 = 49 years of non-burn-in simulation, -1 to account for the fact that Rcpp functions start from 0)
#plot_total  <- T # T if you want a line plot of the total or F for the map,

pop_groups  <- seq(1, 12)

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

# lists to hold data for plots
SIM_sp_f <- list()
SIM_sp_ntz <- list()
SIM_n_dist <- list()
SIM_n_catches <- list()
SIM_age_catches <- list()
SIM_weight_catches <- list()


## Run model simulations ------------------------------------------------------

start = Sys.time()
for (SIM in 1:sim_n){ # Simulation loop - CHARLOTTE HAD 100, I'M STARTING WITH 10

  ## Set up initial population ------------------------------------------------
  pop_total       <- array(0, dim = c(max_cell, 12, max_year)) # number of fish of all ages in our population, in each cell (row), each month (column) and each year simulated (matrix slice)
  total           <- array(NA, dim = c(max_year, 1)) # for plotting purposes, population summed for each year
  
  print(paste0("Simulation number ", SIM)) # progress update
  
  yearly_total <- readRDS("data/output_data/06_burn_in_population.rds")
  
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
