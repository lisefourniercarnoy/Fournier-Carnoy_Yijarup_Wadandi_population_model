# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    ?
# Task:    Set up a population out of the population parameters and fishing effort
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    April 2025

# -----------------------------------------------------------------------------

# Status: THIS IS  A TEST MODIFICATION OF THE SINGLE-FLEET BURN IN CODE. SO FAR IT WORKS FINE

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

sourceCpp("functions/C_model_RcppArm_test.cpp") # currently using a test version of the model - which has 3 sources of effort. 
source("functions/X_Functions.R")


## Colours --------------------------------------------------------------------

## Create colours for the plot
colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))
my.colours <- "PuBu"

## Load files -----------------------------------------------------------------

adult_movement <- readRDS("data/output_data/03_adult_movement_10_swim_speed.rds") %>% glimpse()
effort_com <- readRDS("data/output_data/04A_commercial_burn_in_fishing.rds") %>% glimpse()
effort_s_rec <- readRDS("data/output_data/04B_shore_rec_burn_in_fishing.rds") %>% glimpse()
effort_b_rec <- readRDS("data/output_data/04C_boat_rec_burn_in_fishing.rds") %>% glimpse()

# replacing effort zero to a small positive value, otherwise the function can't compute things well.
effort_com[effort_com == 0] <- 1e-10 
effort_b_rec[effort_b_rec == 0] <- 1e-10
effort_s_rec[effort_s_rec == 0] <- 1e-10

# no_take <- readRDS("data/output_data/02_no_take_list.rds") %>% glimpse()
water             <- readRDS("data/output_data/02_watergrid.rds"); plot(water)
starting_pop      <- readRDS("data/output_data/05_starting_population.rds") %>% glimpse()
selectivity_com   <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOw THEY ARE THE SAME BUT SHOULD END UP BEING DIFFERENT AT SOME POINT
selectivity_b_rec <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOW THEY ARE THE SAME BUT SHOULD END UP BEING DIFFERENT AT SOME POINT
selectivity_s_rec <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOW THEY ARE THE SAME BUT SHOULD END UP BEING DIFFERENT AT SOME POINT

mature          <- readRDS("data/output_data/05_maturity.rds") %>% glimpse()
weight          <- readRDS("data/output_data/05_weight.rds") %>% glimpse()
settlement      <- readRDS("data/output_data/03_recruitment.rds") %>% glimpse()

n_yrs_modelled <- 50 # number of burn-in years - at least a fish's full life. doing 50 here like Charlotte

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
nat_mort = 0.146 # NOT CHANGED FROM CHARLOTTE, NEED TO FIND SOMEWHERE

# Beverton-Holt Recruitment Values - Have sourced the script but need to check that alpha and beta are there
BHa = 0.4344209 # NOT CHANGED FROM CHARLOTTE, NEED TO FIND SOMEWHERE
BHb = 0.0002349538 # NOT CHANGED FROM CHARLOTTE, NEED TO FIND SOMEWHERE
PF = 0.5 # proportion expected to be females
hyperallo <- 1.24 # average from 3 Sparids in Barneche 2018 (1.26, 1.14 and 1.33)

# Model settings
max_cell    <- nrow(water) # Number of cells in the model
max_age     <- 30 # The max age of the fish in the model (-1 to account for the fact that Rcpp functions start from 0)
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
year_to_plot <- 30
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