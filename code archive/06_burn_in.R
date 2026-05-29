# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    ?
# Task:    Set up a population out of the population parameters and fishing effort
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    April 2025

# -----------------------------------------------------------------------------

# Status: THIS IS AN ARCHIVE, USING CHARLOTTE'S FUNCTIONS, NOT MY NEW ONES

# -----------------------------------------------------------------------------

rm(list = ls())

# Load libraries
library(tidyverse) # for data manipulation
library(sf)
# library(forcats)
# library(RColorBrewer)
library(MQMF)
library(Rcpp)
library(RcppArmadillo)
library(abind)

## Read in the functions ------------------------------------------------------
sourceCpp("functions/X_Model_RccpArm.cpp")
source("functions/X_Functions.R")

## Colours --------------------------------------------------------------------

## Create colours for the plot
pop.groups <- c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150)
my.colours <- "PuBu"

## Load files -----------------------------------------------------------------

adult_movement <- readRDS("data/output_data/03_adult_movement_10_swim_speed.rds") %>% glimpse()
effort <- readRDS("data/output_data/04A_commercial_burn_in_fishing.rds") %>% glimpse()

effort[effort == 0] <- 1e-10 # replacing effort zero to a small positive value, otherwise the function can't compute things well.

# no_take <- readRDS("data/output_data/02_no_take_list.rds") %>% glimpse()
water           <- readRDS("data/output_data/02_watergrid.rds"); plot(water)
starting_pop    <- readRDS("data/output_data/05_starting_population.rds") %>% glimpse()
selectivity     <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()
mature          <- readRDS("data/output_data/05_maturity.rds") %>% glimpse()
weight          <- readRDS("data/output_data/05_weight.rds") %>% glimpse()
settlement      <- readRDS("data/output_data/03_recruitment.rds") %>% glimpse()

n_yrs_modelled <- 80 # number of burn-in years

# selectivity changes with length limits- charlotte had a couple of length limit changes, so she chose the most recent length limit to run the model
# she then makes enough matrix layers for the model to run (in my case, for n_yrs_modelled- so dim() should be 30 age classes x 12 months x at least n_yrs_modelled)
selectivity <- selectivity[, , 44] # selecting the most recent selectivity-retention
for(i in 1:7){
  selectivity <- abind(selectivity, selectivity, along=3)
}
dim(selectivity)[3] > n_yrs_modelled


## Plot checks before running the burn-in -------------------------------------

# age structure of the starting population
ggplot(data = starting_pop) +
  geom_line(aes(x = starting_pop$age, y = starting_pop$N)) +
  ggtitle("age structure of the starting population") + xlab("age") + ylab("number of fish") +
  theme_minimal()

# Selectivity of the fish size over time
sel_df <- as.data.frame(selectivity[, 1, ])  # dimensions: age x month x year, keep only jan to flatten
colnames(sel_df) <- paste0("year_", 1:ncol(sel_df))
sel_df$age <- 1:nrow(sel_df)
sel_long <- sel_df %>% pivot_longer(cols = starts_with("year_"), names_to = "year", values_to = "selectivity") %>% mutate(year = as.numeric(gsub("year_", "", year)))

ggplot(sel_long, aes(x = year, y = selectivity, group = age, color = factor(age))) +
  geom_line() +
  scale_color_viridis_d(option = "plasma") +
  theme_minimal() +
  labs(
    y = "Selectivity",
    color = "Age"
  )

# maturity % with age
ggplot(data = as.data.frame(mature) %>% mutate(age = 1:nrow(mature))) +
  geom_line(aes(x = age, y = V1)) +
  ggtitle("age structure of the starting population") + xlab("age") + ylab("% mature fish") +
  theme_minimal()

# settlement probability
water2 <- water %>% mutate(settlement = settlement)
plot(water2["settlement"], 
     main = "settlement probability")

# weight with age
ggplot(data = as.data.frame(weight) %>% mutate(age = 1:nrow(mature))) +
  geom_line(aes(x = age, y = V1)) +
  ggtitle("weight and fish age") + xlab("age") + ylab("weight kg") +
  theme_minimal()


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

total <- array(0, dim = c(max_year, 1))

yearly_total <- array(0, dim = c(max_cell, 12, max_age)) # for every cell (row), and every month (column) across all fish ages (matrix slice), we will have a population

start_pop_year <- starting_pop %>% 
  slice(which(row_number() %% 12 == 1)) # This sets the population in January to be the same population as the December just before

for(d in 1:dim(yearly_total)[3]){ # for every age of the fish...
  for(N in 1:starting_pop[d, 1]){ # and every number of fish of that age...
    cellID <- ceiling((runif(n=1, min = 0, max = 1))*max_cell) # select a random cell...
    yearly_total[cellID, 1, d] <- yearly_total[cellID, 1, d] + 1 # put a fish in the cell
  }
} # this loop randomly puts fish from the starting population (fish of different ages) in the cells


## Run the burn-in period -----------------------------------------------------

pop_total <- array(0, dim = c(max_cell, 12, max_year)) # This is our total population, all ages are summed and each column is a month (each layer is a year)
dim(selectivity)[3] > max_year # check there are enough layers to run. otherwise change line 53

Start=Sys.time()
for (YEAR in 0:(max_year-1)){ #max_year-1 because it starts at 0.
  
  print(YEAR)

  # Loop over all the Rcpp functions in the model
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
                                    
                                    YearlyTotal = yearly_total, 
                                    Select = selectivity, 
                                    Effort = effort
    )

  yearly_total <- ModelOutput$YearlyTotal
  
  # Save some outputs from the model 
  pop_total[ , , YEAR+1] <- rowSums(ModelOutput$YearlyTotal[, , 1:max_age], dim = 2)
  
  water$pop <- pop_total[ , 12, YEAR + 1] # We just want the population at the end of the year
  
  total[YEAR+1, 1] <- sum(water$pop)
}

End=Sys.time()
Runtime = End - Start
Runtime

# Plot check
# the population should be stable by the end of the burn-in
total <- as.data.frame(total)
ggplot(total, aes(x = seq(1, max_year, 1), y = V1)) +
  geom_line() +
  theme_minimal()

# Save burn in population for use in the actual model
saveRDS(yearly_total, file = "data/output_data/06_burn_in_population.rds")

## CHARLOTTE DOES ANOTHER RUN WITH HIGH LEVELS OF FISHING MORTALITY, I WONT FOR NOW

### END ###