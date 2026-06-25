# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    ?
# Task:    Set up a population out of the population parameters and fishing effort
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    April 2025

# -----------------------------------------------------------------------------

# Status: Doing checks here and there but it's functioning. might need to play with the initial recruitment a bit

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
  "shore_rec"
                 )
com_info <- readRDS("data/output_data/04_A_commercial_fishing_info.rds") %>% glimpse()
com_info$fishing_days[com_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

brec_info <- readRDS("data/output_data/04_C_boat_rec_fishing_info.rds") %>% glimpse()
brec_info$fishing_days[brec_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

srec_info <- readRDS("data/output_data/04_B_shore_rec_fishing_info.rds") %>% glimpse()
srec_info$fishing_days[srec_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.


# for the burn-in, we'll use a constant low level of fishing, so replace all years with 1900 fishing effort
com_info$fishing_days[,,1:dim(com_info$fishing_days)[[3]]] <- com_info$fishing_days[,,1]
brec_info$fishing_days[,,1:dim(brec_info$fishing_days)[[3]]] <- brec_info$fishing_days[,,1]
srec_info$fishing_days[,,1:dim(srec_info$fishing_days)[[3]]] <- srec_info$fishing_days[,,1]

# add scaling for each fleet, to account for different gears being less efficient, compared to commercial. see g-sheets Fishing effort reconstruction
brec_info$fishing_days <- brec_info$fishing_days * 0.1
srec_info$fishing_days <- srec_info$fishing_days * 0.017

fleet_info = list(com_info, 
                  brec_info,
                  srec_info
                  )


## Set up the initial population ----------------------------------------------

total <- array(0, dim = c(max_year, 1))

current_pop <- array(0, dim = c(max_cell, n_lengths, max_age)) # for every cell (row), and every month (column) across all fish ages (matrix slice), we will have a population

BURN_IN_pop <- list() # keep record of the burn-in outputs
BURN_IN_catch_weight <- list()
BURN_IN_catch_number <- list()
BURN_IN_F   <- list()
BURN_IN_SSB <- list()

for(AGE in 1:max_age){
  total_this_age <- starting_pop[AGE, ]
  # distribute proportionally to settlement — same habitat weighting as recruits
  settlement_prop <- settlement / sum(settlement) # length max_cell
  
  # outer product: each cell's proportion × each length-class abundance
  current_pop[, , AGE] <- outer(settlement_prop, total_this_age)
}

cat("Total fish initialised:", sum(current_pop), "\n")
cat("Age structure check:\n")
print(round(colSums(current_pop[,,1])))


## start the burn-in ----------------------------------------------------------

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
  
  BURN_IN_pop[[YEAR+1]] <- current_pop # the format is (cell x length x age) x year
  fishing_mortality  <- ModelOutput$fishing_mortalities
  
  BURN_IN_pop[[YEAR+1]] <- current_pop
  BURN_IN_F[[YEAR+1]]   <- ModelOutput$annual_F          # cell x length
  BURN_IN_SSB[[YEAR+1]] <- ModelOutput$spawning_biomass  # cell, summed across spawning months
  
  # then for a quick annual check you can print
  cat("Year", YEAR+1, 
      "| Total SSB:", sum(BURN_IN_SSB[[YEAR+1]]),
      "| Mean F (fished lengths):", mean(BURN_IN_F[[YEAR+1]]),
      "\n")
  
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

End = Sys.time()
Runtime = End - Start
Runtime

## save the burn-in population ------------------------------------------------

saveRDS(current_pop, file = "data/output_data/06_burn_in_population.rds")
SSB0 <- mean(SSB_timeseries[(max_year-10):max_year])  # mean of last 10 burn-in years


## a few sanity checks --------------------------------------------------------

### burn-in population stability over time ------------------------------------

total_pop <- lapply(BURN_IN_pop, function(pop) sum(pop))
par(mfrow = c(1,1))
plot(1:max_year, total_pop, type = "l",
     xlab = "Year", ylab = "Total abundance",
     main = "total Population during burn-in")
# we are looking for a population that ends up being somewhat stable in the last few years of the burn in, with no huge dips of spikes.


### finite f ------------------------------------------------------------------

# finite fishing mortality (f) is the % of fish removed from the population over any time period.
F_by_length_year <- sapply(BURN_IN_F, function(f) colMeans(f))  # n_lengths x n_years
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

SSB_timeseries <- sapply(BURN_IN_SSB, sum)  # sum over cells, one value per year
par(mfrow = c(1, 1))
plot(SSB_timeseries, type = "l", lwd = 2, col = colour_palette[6],
     xlab = "Year", ylab = "Total SSB",
     main = "Spawning Stock Biomass over time")


### fish density --------------------------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds")

ages_to_plot <- c(1, 5, 10, 20, 30, 40)
age_df <- lapply(ages_to_plot, function(a) {
  cell_totals <- apply(BURN_IN_pop[[60]][,, a], 1, sum)
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

avg_weight <- sum(ModelOutput$catch_weight_by_fleet[1,][[1]]) / 
  sum(ModelOutput$catch_number_by_fleet[1,][[1]])

n_fleets <- length(fleet_names)
fleet_totals <- sapply(1:n_fleets, function(f) {
  arr <- ModelOutput$catch_weight_by_fleet[[f, 1]]
  sum(arr)  # sum everything: cells, months, ages/years
})

plot_df <- data.frame(
  fleet = fleet_names[1:n_fleets],
  catch_kg = fleet_totals
)

ggplot(plot_df, aes(x = fleet, y = catch_kg, fill = fleet)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(legend.position = "none")



## END ##
