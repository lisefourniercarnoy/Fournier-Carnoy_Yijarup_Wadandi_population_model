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

## Read files -----------------------------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds"); plot(water$geometry)

## Get model parameters -------------------------------------------------------

n_yrs_modelled <- 5 # number of burn-in years - at least a fish's full life.

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

  # save the population at the end of the month (to give to the loop again as January population)
  current_pop <- ModelOutput$next_pop #  cell x length x age
  fishing_mortality <- ModelOutput$fishing_mortalities
  
  BURN_IN_pop[[YEAR+1]] <- current_pop # the format is (cell x length x age) x year
  
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
       type = "l", lwd = 2, col = "steelblue",
       xlab = "Age class", ylab = "Total abundance",
       main = paste0("Age structure — Year ", YEAR + 1))
  
  plot(size_structure,
       type = "l", lwd = 2, col = "darkred",
       xlab = "Length bin", ylab = "Total abundance",
       main = paste0("Size structure — Year ", YEAR + 1))
  
  dev.off()
  
}
End = Sys.time()
Runtime = End - Start
Runtime

## a few sanity checks --------------------------------------------------------

### finite f ------------------------------------------------------------------

# finite fishing mortality (f) is the % of fish removed from the population over any time period.
# (fishing mortality should be lower for small fish, and vice versa, and higher where there's loads of fishing and vice versa)
f <- 1 - exp(-(fishing_mortality)) 

age_to_plot <- 3
age_f <- f[,, age_to_plot]
water2 <- readRDS("data/output_data/02_watergrid.rds")
water2$f <- age_f

mapview::mapview(water2[,c("geometry", "f")], zcol = "f")


### burn-in population stability over time ------------------------------------

total_pop <- lapply(BURN_IN_pop, function(pop) sum(pop))

plot(1:max_year, total_pop, type = "l",
     xlab = "Year", ylab = "Total abundance",
     main = "Population during burn-in")
BURN_IN_pop[[60]][1,1,1]

### age structure in the last time step ---------------------------------------

total_by_length_age <- apply(current_pop, c(2,3), sum)

age_structure <- colSums(total_by_length_age)
size_structure <- rowSums(total_by_length_age)

par(mfrow = c(1,2))

plot(age_structure, type = "l", xlab = "Age", ylab = "Total abundance",
     main = "Age structure")

plot(size_structure, type = "l", xlab = "Length bin", ylab = "Total abundance",
     main = "Size structure")

### fish density --------------------------------------------------------------

year_to_plot <- 59
age_to_plot <- 10
month_index <- 12
monthly_fish_df <- as.data.frame(BURN_IN_pop[[year_to_plot]][,,age_to_plot])

colnames(monthly_fish_df) <- paste0("month_", 1:12)
water <- readRDS("data/output_data/02_watergrid.rds")
water2 <- cbind(water, monthly_fish_df)
month_col <- paste0("month_", month_index)

water2$test <- water2[[month_col]]/water2$cell_area
mapview::mapview(water2[,c("geometry", "test")], zcol = "test")


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


## save the burn-in population ------------------------------------------------

saveRDS(yearly_pop, file = "data/output_data/06_burn_in_population.rds")


## BELOW IS ARCHIVE -----------------------------------------------------------

## Colours --------------------------------------------------------------------

## Create colours for the plot
colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))
my.colours <- "PuBu"

## Load files -----------------------------------------------------------------

adult_movement <- readRDS("data/output_data/03_b_adult_movement_10_swim_speed.rds") %>% glimpse()
effort_com <- readRDS("data/output_data/04_A_commercial_burn_in_fishing.rds") %>% glimpse() ## MUST NOT HAVE ZEROES OR NAs
effort_s_rec <- readRDS("data/output_data/04_B_shore_rec_burn_in_fishing.rds") %>% glimpse() ## MUST NOT HAVE ZEROES OR NAs
effort_b_rec <- readRDS("data/output_data/04_C_boat_rec_burn_in_fishing.rds") %>% glimpse() ## MUST NOT HAVE ZEROES OR NAs

# replacing effort zero to a small positive value, otherwise the function can't compute things well.
effort_com[effort_com == 0] <- 1e-10 
effort_b_rec[effort_b_rec == 0] <- 1e-10
effort_s_rec[effort_s_rec == 0] <- 1e-10

effort_com[is.na(effort_com)] <- 1e-10 
effort_b_rec[is.na(effort_b_rec)] <- 1e-10
effort_s_rec[is.na(effort_s_rec)] <- 1e-10

# no_take <- readRDS("data/output_data/02_no_take_list.rds") %>% glimpse()
water             <- readRDS("data/output_data/02_watergrid.rds"); plot(water)
starting_pop      <- readRDS("data/output_data/05_starting_population.rds") %>% glimpse()
selectivity_com   <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOw THEY ARE THE SAME BUT SHOULD END UP BEING DIFFERENT AT SOME POINT
selectivity_b_rec <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOW THEY ARE THE SAME BUT SHOULD END UP BEING DIFFERENT AT SOME POINT
selectivity_s_rec <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOW THEY ARE THE SAME BUT SHOULD END UP BEING DIFFERENT AT SOME POINT

mature          <- readRDS("data/output_data/05_maturity.rds") %>% glimpse()
weight          <- readRDS("data/output_data/05_weight.rds") %>% glimpse()
settlement      <- readRDS("data/output_data/03_B_recruitment.rds") %>% glimpse()

n_yrs_modelled <- 10 # number of burn-in years - at least a fish's full life. doing 50 here like Charlotte

# selectivity changes with length limits- charlotte had a couple of length limit changes, so she chose the most recent length limit to run the model
# she then makes enough matrix layers for the model to run (in my case, for 45 years- so dim() should be 30 age classes x 12 months x at least 45 years)
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
dim(selectivity_s_rec)[3] > n_yrs_modelled


## Model settings -------------------------------------------------------------

# Natural Mortality
nat_mort = 0.12 # from table 4.2, p12 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
BHa = 1.201632 # see the second alpha in script 5
BHb = 0.0001889668 # see the second beta in script 6
PF = 0.5 # proportion expected to be females

# Model settings
max_cell    <- nrow(water) # Number of cells in the model
max_age     <- 40-1 # The max age of the fish in the model, see script 05 for correct value (-1 to account for the fact that Rcpp functions start from 0)
max_year    <- n_yrs_modelled # Number of years the model should run for 
plot_total  <- T # T if you want a line plot of the total or F for the map,

pop_groups  <- seq(1, 12)


## Set up the initial population ----------------------------------------------

total <- array(0, dim = c(max_year+1, 1))

yearly_total <- array(0, dim = c(max_cell, 12, max_age)) # for every cell (row), and every month (column) across all fish ages (matrix slice), we will have a population

start_pop_year <- starting_pop %>% 
  slice(which(row_number() %% 12 == 1)) # This sets the population in January to be the same population as the December just before

for(d in 1:dim(yearly_total)[3]){ # for every age of the fish...
  for(N in 1:starting_pop[d, 1]){ # and every number of fish of that age...
    cellID <- ceiling((runif(n=1, min = 0, max = 1))*max_cell) # select a random cell...
    yearly_total[cellID, 1, d] <- yearly_total[cellID, 1, d] + 1 # put a fish in the cell
  }
} # this loop randomly puts fish from the starting population (fish of different ages) in the cells


## Run the model for the burn-in ----------------------------------------------

pop_total <- array(0, dim = c(max_cell, 12, max_year)) # This is our total population, all ages are summed and each column is a month (each layer is a year)
dim(selectivity_b_rec)[3] > max_year # check there are enough layers to run. otherwise change in section 'load files'

## TEST

# test dynamic spawning months - master function is modified for this
spawn_months <- c(9, 10, 11) # october-november
spawn_months-1
## END TEST

Start = Sys.time()
for (YEAR in 0:(max_year-1)){ # max_year-1 because it starts at 0.
  
  print(YEAR)
  
  # loop over all the Rcpp functions in the model
  ModelOutput <- RunModelfunc_cpp(YEAR = YEAR,                                   
                                  MaxCell = max_cell,
                                  MaxYear = max_year, 
                                  
                                  MaxAge = max_age, 
                                  NatMort = nat_mort, 
                                  BHa = BHa, 
                                  BHb = BHb, 
                                  PF = PF, 
                                  AdultMove = adult_movement, 
                                  Mature = mature, 
                                  Weight = weight, 
                                  spawn_months = spawn_months,
                                  Settlement = settlement, 
                                  ha_scaling = hyperallo,
                                  
                                  YearlyTotal = yearly_total, 
                                  
                                  Selectivity_com = selectivity_com, 
                                  Selectivity_b_rec = selectivity_b_rec, 
                                  Selectivity_s_rec = selectivity_s_rec, 
                                  
                                  Effort_com = effort_com,
                                  Effort_b_rec = effort_b_rec,
                                  Effort_s_rec = effort_s_rec
  )

  yearly_total <- ModelOutput$YearlyTotal # this is an array of how many fish of each age exist in the first year of the model.
  
  # Save some outputs from the model 
  pop_total[ , , YEAR+1] <- rowSums(ModelOutput$YearlyTotal[, , 1:max_age], dim = 2)
  
  water$pop <- pop_total[ , 12, YEAR + 1] # just keep the population at the end of the year
  
  total[YEAR+1, 1] <- sum(water$pop) # store the total population at the end of the year
  
}

End = Sys.time()
Runtime = End - Start
Runtime

## plot checks
# the burn-in total population
total <- as.data.frame(total)
plot(x = seq(1, max_year+1, 1), y = total$V1)

# number of fish spatially
year_to_plot <- 10
monthly_fish_df <- as.data.frame(pop_total[,,year_to_plot])
colnames(monthly_fish_df) <- paste0("month_", 1:12)
water <- readRDS("data/output_data/02_watergrid.rds"); plot(water)
water2 <- cbind(water, monthly_fish_df)
month_index <- 1
month_col <- paste0("month_", month_index)
ggplot(water2) +
  geom_sf(aes(fill = .data[[month_col]])) +
  scale_fill_viridis_c() +
  labs(
    title = paste0("Fish population - Year ", year_to_plot, ", Month ", month_index),
    fill = "Population"
  ) +
theme_minimal()

## Save burn in population for use in the actual model
saveRDS(yearly_total, file = "data/output_data/06_burn_in_population.rds")


## Make a GIF for shits n gigs ------------------------------------------------

library(gifski)

# Parameters
month_to_plot <- 1
frame_count <- 1
years <- seq_len(dim(pop_total)[3])  # 1 to number of years

# Loop over years and create plots
for (y in seq_along(years)) {
  
  print(y)
  
  # Get the fish population for this year and selected month
  colnames(pop_total) <- paste0("month_", 1:12)
  water$fish_pop <- pop_total[, month_to_plot, y]
  
  overall_min <- min(pop_total, na.rm = TRUE)
  overall_max <- max(pop_total, na.rm = TRUE)
  
  # Plot
  p <- ggplot(water) +
    geom_sf(aes(fill = fish_pop), color = NA) +
    scale_fill_gradientn(
      colours = colour_palette[6:4],
      na.value = NA,
      limits = c(overall_min, overall_max)  # <- This fixes the scale across frames
    ) +
    labs(
      title = paste("Fish Population - Year", years[y], "Month", month_to_plot),
      fill = "Population"
    ) +
    theme_minimal()
  
  # Save the plot as a frame
  ggsave(
    filename = sprintf("plots/gif_frames/06_burn_in_population_gif_frames/06_burn_in_population_frame_%03d.png", frame_count),
    plot = p,
    width = 6, height = 6, dpi = 150
  )
  
  frame_count <- frame_count + 1
}

# Create the animated GIF
png_files <- list.files("plots/gif_frames/06_burn_in_population_gif_frames/", pattern = "06_burn_in_population_frame_\\d+\\.png", full.names = TRUE)

gifski(
  png_files,
  gif_file = "plots/gifs/06_burn_in_population.gif",
  width = 600,
  height = 600,
  delay = 0.25  # Adjust for speed (0.4 = slower)
)

## CHARLOTTE DOES ANOTHER RUN WITH HIGH LEVELS OF FISHING MORTALITY, I WONT FOR NOW - this is to check how sensitive the population is to diff levels of fishing. she found it wasn't sensitive so didn't bother too much.

### END ###


## ARCHIVE TRYING WITH NEW FUNCTION ? -------------------------------------------------

## Load files -----------------------------------------------------------------

adult_movement <- readRDS("data/output_data/03_b_adult_movement_10_swim_speed.rds") %>% glimpse()
selectivity <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()
com_fishing_info <- readRDS("data/output_data/04_A_commercial_fishing_info.rds")

# replacing effort zero to a small positive value, otherwise the function can't compute things well.
com_fishing_info$fishing_days[com_fishing_info$fishing_days == 0] <- 1e-10 

com_fishing_info$fishing_days[is.na(com_fishing_info$fishing_days)] <- 1e-10 


# no_take <- readRDS("data/output_data/02_no_take_list.rds") %>% glimpse()
water             <- readRDS("data/output_data/02_watergrid.rds"); plot(water)
starting_pop      <- readRDS("data/output_data/05_starting_population.rds") %>% glimpse()

mature          <- readRDS("data/output_data/05_maturity.rds") %>% glimpse()
weight          <- readRDS("data/output_data/05_weight.rds") %>% glimpse()
settlement      <- readRDS("data/output_data/03_B_recruitment.rds") %>% glimpse()

n_yrs_modelled <- 10 # number of burn-in years - at least a fish's full life. doing 50 here like Charlotte

# selectivity changes with length limits- charlotte had a couple of length limit changes, so she chose the most recent length limit to run the model for the burn-in
# she then makes enough matrix layers for the model to run (in my case, for n_yrs_modelled years- so dim() should be 30 age classes x 12 months x at least n_yrs_modelled years)
selectivity <- selectivity[, , 44] # selecting the most recent selectivity-retention
for(i in 1:6){
  selectivity <- abind(selectivity, selectivity, along=3)
}

dim(selectivity)[3] > n_yrs_modelled


## Model settings -------------------------------------------------------------

# Natural Mortality
nat_mort = 0.12 # from table 4.2, p12 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
BHa = 1.201632 # see the second alpha in script 5
BHb = 0.0001889668 # see the second beta in script 6
PF = 0.5 # proportion expected to be females

# Model settings
max_cell    <- nrow(water) # Number of cells in the model
max_age     <- 40-1 # The max age of the fish in the model, see script 05 for correct value (-1 to account for the fact that Rcpp functions start from 0)
max_year    <- n_yrs_modelled # Number of years the model should run for 
plot_total  <- T # T if you want a line plot of the total or F for the map,

pop_groups  <- seq(1, 12)


## Set up the initial population ----------------------------------------------

total <- array(0, dim = c(max_year+1, 1))

yearly_total <- array(0, dim = c(max_cell, 12, max_age)) # for every cell (row), and every month (column) across all fish ages (matrix slice), we will have a population

start_pop_year <- starting_pop %>% 
  slice(which(row_number() %% 12 == 1)) # This sets the population in January to be the same population as the December just before

for(d in 1:dim(yearly_total)[3]){ # for every age of the fish...
  for(N in 1:starting_pop[d, 1]){ # and every number of fish of that age...
    cellID <- ceiling((runif(n=1, min = 0, max = 1))*max_cell) # select a random cell...
    yearly_total[cellID, 1, d] <- yearly_total[cellID, 1, d] + 1 # put a fish in the cell
  }
} # this loop randomly puts fish from the starting population (fish of different ages) in the cells



## Run the model for the burn-in ----------------------------------------------

pop_total <- array(0, dim = c(max_cell, 12, max_year)) # This is our total population, all ages are summed and each column is a month (each layer is a year)
dim(selectivity)[3] > max_year # check there are enough layers to run. otherwise change in section 'load files'

## TEST

# test dynamic spawning months - master function is modified for this
spawn_months <- c(9, 10, 11) # october-november
spawn_months-1
## END TEST

Start = Sys.time()
for (YEAR in 0:(max_year-1)){ # max_year-1 because it starts at 0.
  
  print(YEAR)
  
  # loop over all the Rcpp functions in the model
  ModelOutput <- RunModelfunc_cpp(YEAR = YEAR,                                   
                                  MaxCell = max_cell,
                                  MaxYear = max_year, 
                                  
                                  MaxAge = max_age, 
                                  NatMort = nat_mort, 
                                  BHa = BHa, 
                                  BHb = BHb, 
                                  PF = PF, 
                                  AdultMove = adult_movement, 
                                  Mature = mature, 
                                  Weight = weight, 
                                  spawn_months = spawn_months,
                                  Settlement = settlement, 
                                  ha_scaling = hyperallo,
                                  
                                  YearlyTotal = yearly_total, 
                                  
                                  Selectivity_com = selectivity_com,
                                  
                                  Effort_com = effort_com,
                                  Effort_b_rec = effort_b_rec,
                                  Effort_s_rec = effort_s_rec
  )
  
  yearly_total <- ModelOutput$YearlyTotal # this is an array of how many fish of each age exist in the first year of the model.
  
  # Save some outputs from the model 
  pop_total[ , , YEAR+1] <- rowSums(ModelOutput$YearlyTotal[, , 1:max_age], dim = 2)
  
  water$pop <- pop_total[ , 12, YEAR + 1] # just keep the population at the end of the year
  
  total[YEAR+1, 1] <- sum(water$pop) # store the total population at the end of the year
  
}

End = Sys.time()
Runtime = End - Start
Runtime

## plot checks
# the burn-in total population
total <- as.data.frame(total)
plot(x = seq(1, max_year+1, 1), y = total$V1)

# number of fish spatially
year_to_plot <- 10
monthly_fish_df <- as.data.frame(pop_total[,,year_to_plot])
colnames(monthly_fish_df) <- paste0("month_", 1:12)
water <- readRDS("data/output_data/02_watergrid.rds"); plot(water)
water2 <- cbind(water, monthly_fish_df)
month_index <- 1
month_col <- paste0("month_", month_index)
ggplot(water2) +
  geom_sf(aes(fill = .data[[month_col]])) +
  scale_fill_viridis_c() +
  labs(
    title = paste0("Fish population - Year ", year_to_plot, ", Month ", month_index),
    fill = "Population"
  ) +
  theme_minimal()

## Save burn in population for use in the actual model
saveRDS(yearly_total, file = "data/output_data/06_burn_in_population.rds")

