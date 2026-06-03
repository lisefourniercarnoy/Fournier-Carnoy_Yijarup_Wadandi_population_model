# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Values from the literature
# Task:    Set up populations
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    May 2026

# -----------------------------------------------------------------------------

# Status: May 2026 update: trying to switch from age-based to length-based system

# -----------------------------------------------------------------------------

rm(list = ls())

# Load libraries

library(tidyverse)
library(sf)
library(raster)
library(stringr)
library(forcats)
library(RColorBrewer)
library(geosphere)
library(abind)

# we need a population to start the C++ model.
# the start population must be stable and realistic, given the species, the size of the simulation area etc.
# for this, we must know (1.) how the species grows/matures, then (2.) what population size/structure the area *can* support (carrying capacity),
# which will allow us to calculate (3.) what the population size/structure would be under a stable amount of fishing.

# the first two steps are calculated in a 'per recruit' unit, which takes a hypothetical R0=1 recruit through all ages, calculating the hypothetical survival probability, maturity probability and mature biomass
# we then rescale the values from these two first steps to real numbers of fish.

# from this, we can calculate (4.) our start population.
# we need step 2. because the Beverton-Holt equation (which the C++ code runs) calculates density-dependent recruitment (limited by the carrying capacity)


## STEP 0: life history values ================================================

year_end <- 2024
year_start <- 1900
n_yrs <- 2024-1900 # number of years modelled

step = 1/12 # monthly timestep
max_age <- 40 # table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15


### length parameters ----

# Von Bertalanffy Parameters (from Table 2, lower west coast, both sexes combined, https://academic.oup.com/icesjms/article/74/1/180/2669555?login=false#186901916)
# we're not using this anymore, see schnute below
# Linf  <- 1136 # hypothetical asymptotic length at infinite age (mm)
# k     <- 0.12 # growth coefficient, (/year)
# t0    <- -0.42 # hypothetical age at length = 0 (y)

# Schnute equation parameters - found to fit Yijarup better, table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
l1 <- 0
l2 <- (927.2 + 915.9)/2 # average of females and males
t1 <- 0
t2 <- 20
a <- (0.160 + 0.168)/2 # average of females and males
b <- (1.239 + 1.206)/2

max_length <- 1095 # in mm, table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15

# age-length standard deviation parameters, for equation 3 in Francis et al. 2016, https://www.sciencedirect.com/science/article/pii/S0165783615000909
al_a <- 50 # manually selected to kinda fit figure 6 in Wakefield et al. 2017  https://academic.oup.com/icesjms/article/74/1/180/2669555?login=true&guestAccessKey=
al_b <- 0.2 # manually selected to kinda fit figure 6 in Wakefield et al. 2017  https://academic.oup.com/icesjms/article/74/1/180/2669555?login=true&guestAccessKey=


### weight parameters ----

# weight-length relationship (W = WLa*(forklength^WLb)) , table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
WLa <- 3.67*10^-5
WLb <- 2.83


### maturity parameters ----

# age at maturity
M50 <- (6.7 + 5.5)/2 # female male average, A50 in table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
M95 <- (12.9 + 12.3)/2 # female male average, A95 table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15

# length at maturity
L50 <- 660 # females only available, L50 in table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
L95 <- 900 # females only available, L95 in table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15


### spawning/recruitment parameters ----

prop_f <- 0.5 # proportion expected to be females

# Beverton-Holt Parameters
h   <- 0.75 # table 7.4 p.68 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
R0  <- 1 # initial recruitment, hypothetical number. do not change!!

# batch fecundity = a * length^b
fecundity_a <- 9.436*10^-5 # table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
fecundity_b <- 3.359 # table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15

hyperallo <- (1.26 + 1.14 + 1.33)/3 # THIS VALUE IS CORRECT FOR WEIGHT, NOT AGE OR LENGTH !! average from 3 Sparids in Barneche 2018 (1.26, 1.14 and 1.33)


### mortality/fishing parameters ----

M <- 0.12 # yearly natural mortality, table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15

# selectivity - for now considered as the same for rec and commercial fishing, as even Wise et al. 2007 group them together, p.99. https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1206&context=fr_rr#page=94.59 
# parameter of selectivity ogive
V50 <- 3.701 # Table 6.3.8, v50 from https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
V95 <-  5.290 # Table 6.3.8, v95 from https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr

PRM <- mean(c(0.15, 0.14, 0.16, 0.13, 0.17, 0.22, 0.20, 0.23, 0.19, 0.25, 0.15, 0.13, 0.16, 0.11, 0.19)) # post release mortality for our retention function, recreational south, p.109 in Fairclough 2021 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1108&context=fr_rr#page=118.15

# McLennan et al. 2014 finds post-release survival after barotrauma treatment to be 88%, same for all depths, all sizes of fish.
# Maggs et al. 2024 finds the hook site (lip, body or gut) and depth both sig. affects PRM.
# Grixti et al. 2010 finds that PRM in deep waters (48%) is lower than shallow (97%).

eq.init.fish <-  0.025 # this is equilibrium instantaneous fishing mortality?? a low level of fishing effort to apply to the hypothetical fish population so that it stabilises for the burn-in.


## STEP 1: set up life history value ==========================================

# below we use the published life history values to obtain mean values of length-at-age and the weight-at-mean-length-at-age.
mean_life_hist <- data.frame(
  age = seq(1, 1 + step * (max_age * 12 - 1), by = step)) %>% # each age (1, 2, 3 year old) is split into 12
  mutate(
    length = ((l1^b+(l2^b-l1^b))*((1-exp(-a*(age-t1)))/(1-exp(-a*(t2-t1)))))^(1/b), # this is the Schnute equation, replaces the VB equation.
    # length_vb = Linf * (1 - exp(-k * (age - t0))), # this is the old VB equation
    weight = WLa * (length^WLb) / 1000 # divide by 1000 to get in kg 
  ) %>% 
  glimpse()

# from this, let's add uncertainty (because not every 1yo fish will be the same length)
# we'll make a age x length matrix that tells us how likely it is that a fish of age x is of length a, b, c, etc.
length_age_matrix <- matrix(nrow = length(mean_life_hist$age),
                            ncol = length(mean_life_hist$length),
                            dimnames = list(as.numeric(mean_life_hist$age), 
                                            as.numeric(mean_life_hist$length)))

for (i in 1:nrow(length_age_matrix)) {
  
  mean_length = mean_life_hist$length
  sd_length = al_a + al_b * mean_length # eq. 3 from Francis et al. 2016, https://www.sciencedirect.com/science/article/pii/S0165783615000909
  
  length_age_matrix[i, ] <- dnorm(mean_life_hist$length, mean = mean_length, sd = sd_length)
} # this loops calculates how likely it is that a fish of age x is of each length y

length_age_matrix <- length_age_matrix / rowSums(length_age_matrix) # normalise so it all adds up to 1

dimnames(length_age_matrix) <- list(mean_life_hist$age, mean_life_hist$length)
summary(rowSums(length_age_matrix))  # should all be 1

# okay so this matrix is going to be the base on which most things get calculated downstream.


## STEP 2: set up a hypothetical unfished population ==========================

# the goal here is to calculate what the distribution of mature biomass looks like for a hypothetical population.
# what this means is that for a hypothetical recruitment of R0 = 1, where is the bulk of biomass across lengths and ages?
# this whole setup is made to obtain a single value: unfished_female_spawning_biomass.
# we use this downstream.

# calculate maturity and mature biomass across all lengths
length_mature_biomass <- data.frame(
  length = mean_life_hist$length,
  weight = mean_life_hist$weight
  ) %>% 
  mutate(
    maturity = 1/(1+(exp(-log(19)*((length-L50)/(L95-L50))))), # you can replace length, l50, and l95 with age, M50, and M95 if you wish, but length determines maturity better than age
    mature_biomass = maturity * weight^hyperallo
  ) %>% 
  glimpse()

plot(x = length_mature_biomass$length, 
     y = length_mature_biomass$mature_biomass, 
     main = "mature biomass by length", 
     type = "l")


# calculate natural survival in each age group
age_survival <- data.frame(
  age = mean_life_hist$age,
  survival = NA
) %>% 
  glimpse()
age_survival[1, "survival"] <- R0*prop_f # survival probability on age 1
for (r in 2:nrow(age_survival)){
  age_survival[r, "survival"] <- age_survival[r-1, "survival"]*exp(-M/12) 
} # this calculates the survival in the next age based on the previous age

plot(x = age_survival$age, 
     y = age_survival$survival,
     main = "probability of surviving to this age (unfished population", 
     type = "l")


# put it all together
unfished_pop_ssb <- sweep(length_age_matrix, MARGIN = 1, length_mature_biomass$mature_biomass, FUN = "*")
test <- sweep(unfished_pop_ssb, MARGIN = 2, age_survival$survival, FUN = "*")

# below we are taking a hypothetical female, which we weigh every year (12th month), to figure out how much total spawning biomass it contributes over its lifetime.
# with the matrix, we're kind of collapsing the uncertainty (summing), so that we get a single value.
unfished_female_spawning_biomass <- sum(test[seq(12, length(test), by = 12)])
unfished_female_spawning_biomass # about 64kg

### Calculate Beverton-Holt parameters for the unfished population ----

alpha <- (unfished_female_spawning_biomass/R0)*((1-h)/(4*h)) 
beta <- ((h-0.2)/(0.8*h*R0))


## STEP 3: set up a hypothetical fished population ============================

# we need to calculate things for a hypothetical fished population.
# survival is a little different when you add fishing mortality, so we need to recalculate.

# we have biomass across lengths from the previous step
glimpse(length_mature_biomass)

fished_pop_setup <- data.frame(age = mean_life_hist$age) %>% 
  mutate(selectivity = 1 / (1 + exp(-log(19) * ((age - V50) / (V95 - V50)))), # selectivity at each age (this is probably a length-based process but the stock assessment reports age-based values for V50 and V95)
         fishing_mortality = selectivity * eq.init.fish, # fishing mortality
         total_mortality = fishing_mortality + M
         )
head(fished_pop_setup)

# Calculate survival for the fished population
fished_pop_setup <- fished_pop_setup %>% 
  mutate(survival = NA) %>% 
  glimpse()

fished_pop_setup[1, "survival"] <- R0*prop_f # survival probability on age 1

for (r in 2:nrow(fished_pop_setup)){
  fished_pop_setup[r, "survival"] <- fished_pop_setup[r-1, "survival"]*exp(-fished_pop_setup[r-1, "total_mortality"]*step) # Divide by the time step here
} # this calculates the survival in the next age based on the previous age

# put it all together
fished_pop_ssb <- sweep(length_age_matrix, MARGIN = 1, length_mature_biomass$mature_biomass, FUN = "*")
fished_pop_ssb <- sweep(fished_pop_ssb, MARGIN = 2, fished_pop_setup$survival, FUN = "*")

# below we are taking a hypothetical female, which we weigh every year (12th month), to figure out how much total spawning biomass it contributes over its lifetime.
# with the matrix, we're kind of collapsing the uncertainty (summing), so that we get a single value.
fished_female_spawning_biomass <- sum(fished_pop_ssb[seq(12, length(fished_pop_ssb), by = 12)])
fished_female_spawning_biomass # about 58kg, lower than unfished, which makes sense because a fished population will have fewer big fish (spawners)


### Spawner per recruit and equilibrium recruitment ----

spr <- fished_female_spawning_biomass/unfished_female_spawning_biomass # this says 'under fished equilibrium, 1 recruit loses ~10% (1-spr) of its reproductive output over her lifetime compared to unfished equilibrium'
equil_recr <- (fished_female_spawning_biomass-alpha)/(beta*fished_female_spawning_biomass) # this says 'however, under fished equilibrium, 1 recruit produces ~99% of its unfished recruitment (because there is lower density of fish)
# the density-dependence is visible here: fewer spawners = population below carrying capacity = more recruits


## STEP 4: set up the starting population =====================================

## so far the fished and unfished populations and spawning biomasses were hypothetical, for a single recruit.
## here we scale the spawning biomasses to the level of recruitment, to make a starting population.
## first it's age-based, then we transform it into a age-length-based structure.

init_recr <- 5000 # in thousands - 4000 is too small

# calculate initial fished recruitment (how many new fish from the fished population)
init_fished_recr <- (fished_female_spawning_biomass-alpha) / (fished_female_spawning_biomass*beta) * init_recr 

# calculate initial unfished spawning biomass - this is the carrying capacity for the size of our model, defined by init_recr
init_unfished_f_sb <- init_fished_recr * unfished_female_spawning_biomass

# calculate initial fished spawning biomass - this is the spawning biomass of our starting population
init_fished_f_sb <- init_fished_recr * fished_female_spawning_biomass

# calculate new Beverton-Holt parameters - because we have a different initial recruitment value, alpha and beta will change.
alpha <- (init_unfished_f_sb / init_fished_recr) * ((1-h) / (4*h))
beta <- ((h-0.2) / (0.8*h*init_fished_recr))

# calculate number of female recruits in the next time step
n_f_recr <- (init_fished_f_sb/(alpha+beta*init_fished_f_sb))*prop_f # originally Charlotte had init_unfished_sb/(alpha...) but i think that's incorrect

age_starting_pop <- data.frame(
  age = mean_life_hist$age,
  n = NA
  )

age_starting_pop[1, "n"] <- n_f_recr

for (r in 2:(max_age*12)){ 
  age_starting_pop[r, "n"] <- age_starting_pop[r-1, "n"] * exp(-fished_pop_setup[r-1, "total_mortality"] * step) # divide by the time step here
} # this calculates the survival in the next age based on the previous age using total mortality from our fished population


# okay so now we know that our starting population has ~2456 fish of age 1, but not all will be the same length. we'll use our age-length matrix to figure out how many there are of each length.
age_length_starting_pop <- sweep(length_age_matrix, MARGIN = 2, age_starting_pop[, "n"], FUN = "*")

## These fish form our starting population for the model
## At the end of the year the spawning stock biomass of all females will be calculated to generate recruitment for the next year
## Alpha and Beta need to be recorded for use in the next step of the model

## Selectivity-retention for the fish in the model ----------------------------

# we want to make a matrix that tells us for every month of every year, how likely is it that a fish of length x gets caught and kept by a fisher.

ret <- array(0,
             dim = c(ncol(length_age_matrix)), 
             dimnames = list(colnames(length_age_matrix)))

sel_ret <- array(0,
                 dim = c(ncol(length_age_matrix), n_yrs+1), 
                 dimnames = list(colnames(length_age_matrix), year_start:year_end))

mll_changes <- data.frame( # minimum legal length (MLL) changes
  year = c(1913, 1977, 1988, 2008, 2009),
  mll =  c(279,  380,  410,  450,  500)
)

for (YEAR in year_start:year_end) {
  
  current_mll <- if (YEAR < min(mll_changes$year)) 0 else (
    tail(mll_changes$mll[mll_changes$year <= YEAR], 1)
  )
  
  # retention: the % of fish of each size that are kept by fishers. changes with MLL legislation
  lengths <- as.numeric(rownames(sel_ret))
  ret <- ifelse(lengths < current_mll, 0, # fish that are too small will not be kept
                ifelse(lengths < current_mll + 5, 0.5, # fish around the MLL will not all be kept
                       0.95)) # and the larger fish will mostly be kept, but not all
  
  # landings: the % of fish of each length that can be caught and will be retained
  landings = fished_pop_setup$selectivity * ret
  # discards: the % of fish of each length that can be caught but will be chucked back
  discards = fished_pop_setup$selectivity * (1 - landings)
  
  # selectivity-retention: the % of fish of each age (in each year) that will die from fishing activities (caught, and either kept, or chucked back and die)
  sel_ret[, as.character(YEAR)] = landings + (PRM * discards)
  
}
head(sel_ret)


## Saving files ---------------------------------------------------------------

# initial population
saveRDS(age_length_starting_pop, file = paste0("data/output_data/05_starting_population.rds"))

# selectivity
selectivity <- fished_pop_setup$selectivity
saveRDS(selectivity, file = "data/output_data/05_selectivity.rds")
saveRDS(sel_ret, file = "data/output_data/05_selectivity_retention.rds")

# maturity
maturity <- length_mature_biomass$maturity 
maturity <- array(maturity, dim = c(12, max_age))
maturity <- t(maturity)
saveRDS(maturity, file = "data/output_data/05_maturity.rds")

# weight of each age group
weight <- mean_life_hist$weight
weight <- array(weight, dim = c(12, max_age))
weight <- t(weight)
saveRDS(weight, file = "data/output_data/05_weight.rds")

# Beverton-Holt parameters
saveRDS(alpha, file = "data/output_data/05_Beverton-Holt_alpha.rds")
saveRDS(beta, file = "data/output_data/05_Beverton-Holt_beta.rds")

### END ###
