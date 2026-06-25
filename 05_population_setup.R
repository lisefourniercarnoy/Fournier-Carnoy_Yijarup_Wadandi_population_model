# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    ?
# Task:    Set up populations
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    April 2025

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

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))

# we need a population to start the C++ model.
# the start population must be stable and realistic, given the species, the size of the simulation area etc.
# for this, we must know (1.) what population size/structure the area *can* support (carrying capacity),
# which will allow us to calculate (2.) what the population size/structure would be under a stable amount of fishing.

# the first two steps are calculated in a 'per recruit' unit, which takes a hypothetical R0=1 recruit through all ages, calculating the hypothetical survival probability, maturity probability and mature biomass
# we then rescale the values from these two first steps to real numbers of fish.

# from this, we can calculate (3.) our start population.
# we need step 1. because the Beverton-Holt equation (which the C++ code runs) calculates density-dependent recruitment (limited by the carrying capacity)


## STEP 0: set up life history ------------------------------------------------

year_end <- 2024
year_start <- 1900

### life history parameters ----

# set timestep 
step = 1/12 # monthly timestep

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

# weight-length relationship (W = WLa*(forklength^WLb)) , table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
WLa <- 3.67*10^-5
WLb <- 2.83

prop_f <- 0.5 # proportion expected to be females
M <- 0.12 # yearly natural mortality, table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15

# Beverton-Holt Parameters
h   <- 0.75 # table 7.4 p.68 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
R0  <- 1 # Initial recruitment, hypothetical number. do not change.

max_age <- 40 # table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
max_length <- 1095 # in mm, table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15

# age at Maturity
M50 <- (6.7 + 5.5)/2 # female male average, A50 in table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
M95 <- (12.9 + 12.3)/2 # female male average, A95 table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
hyperallo <- (1.26 + 1.14 + 1.33)/3 # average from 3 Sparids in Barneche 2018 (1.26, 1.14 and 1.33)

L50 <- 660 # females only available, L50 in table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
L95 <- 900 # females only available, L95 in table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15

# fecundity parameters 
fec_a <- 9.436*10^-5 # table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
fec_b <- 3.359 # table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15
# batch fecundity = a * length^b


### fishing parameters ----

# selectivity - for now considered as the same for rec and commercial fishing, as even Wise et al. 2007 group them together, p.99. https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1206&context=fr_rr#page=94.59 
# parameter of selectivity ogive
A50 <- 3.701 # Table 6.3.8, v50 from https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr
A95 <-  5.290 # Table 6.3.8, v95 from https://researchlibrary.agric.wa.gov.au/cgi/viewcontent.cgi?article=1029&context=fr_rr

n_yrs <- 2024-1900 # number of years modelled

PRM <- 0.25 # post release mortality for our retention function - NOT CHANGED FROM CHARLOTTE
# McLennan et al. 2014 finds post-release surival after barotrauma treatment to be 88%, same for all depths, all sizes of fish.
# Maggs et al. 2024 finds the hook site (lip, body or gut) and depth both sig. affects PRM.
# Grixti et al. 2010 finds that PRM in deep waters (48%) is lower than shallow (97%).


# fishing parameters
eq.init.fish <-  0.025 # NOT CHANGED FROM CHARLOTTE, dont know what that is


### set up life history values ------------------------------------------------

life_hist <- as.data.frame(array(1, dim = c((max_age*12), 3))) %>%
  rename(age = "V1") %>% 
  rename(length = "V2") %>% 
  
  rename(weight = "V3") %>% 
  glimpse()

ages <- seq(1, 1 + step * (max_age * 12 - 1), by = step) # ages monthly from year 1 to max year (40) in decimals

life_hist <- data.frame(age = ages) %>%
  mutate(length = ((l1^b+(l2^b-l1^b))*((1-exp(-a*(age-t1)))/(1-exp(-a*(t2-t1)))))^(1/b), # this is the Schnute equation, replaces the VB equation.
         #length_vb = Linf * (1 - exp(-k * (age - t0))), #this is the old VB equation
         weight = WLa * (length^WLb) / 1000) # Divide by 1000 to get in kg 
glimpse(life_hist)

plot(x = life_hist$age, y = life_hist$length)


## STEP 1: set up a hypothetical unfished population --------------------------

unfished_pop_setup <- as.data.frame(array(1, dim = c((max_age*12), 3))) %>% 
  rename(unfished_surv = "V1") %>% # survival
  rename(unfished_mat = "V2") %>% # maturity
  rename(unfished_bio = "V3") %>% # mature biomass
  glimpse()

# Survival in each age group
unfished_pop_setup[1, 1] <- R0*prop_f # survival probability on age 1
for (r in 2:(max_age*12)){ # This calculates the survival in the next age based on the previous age
  unfished_pop_setup[r, 1] <- unfished_pop_setup[r-1, 1]*exp(-M/12) 
}
head(unfished_pop_setup) # survival of the unfished population decreases a bit at every age, with an initial value of 0.5 on day 1

# Proportion of mature individuals per recruit in each age group
unfished_pop_setup <- unfished_pop_setup %>% 
  mutate(age = life_hist$age) %>% 
  mutate(unfished_mat = 1/(1+(exp(-log(19)*((age-M50)/(M95-M50))))))
head(unfished_pop_setup)

# Mature biomass per recruit in each age group
unfished_pop_setup <- unfished_pop_setup %>% 
  mutate(unfished_bio = unfished_mat * unfished_surv * life_hist$weight^hyperallo)
head(unfished_pop_setup)


# quick explanation: survival decreases with age, expected with natural mortality
# maturity increases then plateaus
# biomass increases (fish grow) then decreases again (when natural mortality wipes out a good number of them)
matplot(x = unfished_pop_setup[4], y = unfished_pop_setup[1:3], type = "l", lty = 1, lwd = 2,
        col = 1:ncol(unfished_pop_setup), xlab = "age", ylab = "probability")
legend("right", legend = colnames(unfished_pop_setup)[1:3],
       col = 1:ncol(unfished_pop_setup[1:3]), lty = 1, lwd = 2)


# total mature biomass per recruit
unfished_f_sb <- unfished_pop_setup %>%  # This is the total mature female spawning biomass per recruit for the unfished population
  slice(which(row_number() %% 12 == 1)) %>% # ages 1-40 are in monthly increments, so we obtain the spawning biomass in December for each age
  summarise(., sum = sum(unfished_bio))
unfished_f_sb # if you recruited 1 fish (weighed her every year), over her lifetime (summed for every year), she would contribute 23kg of mature biomass (weighted by survival probability)


## Calculate Beverton-Holt parameters for the unfished population -------------

alpha <- (unfished_f_sb/R0)*((1-h)/(4*h)) 
beta <- ((h-0.2)/(0.8*h*R0))


## STEP 2: set up a hypothetical fished population ----------------------------

fished_pop_setup <- as.data.frame(array(1, dim  = c((max_age*12), 5))) %>% 
  rename(selectivity = "V1",
         fishing_mort = "V2",
         tot_mort = "V3",
         fished_surv = "V4",
         fished_mat = "V5") %>% 
  
  mutate(age = life_hist$age,
         length = life_hist$length,
         fished_mat = unfished_pop_setup$unfished_mat,
         
         selectivity = 1/(1+(exp(-log(19)*((age-A50)/(A95-A50))))), # selectivity for each age group - This is for the equilibrium population
         fishing_mort = selectivity * eq.init.fish, # fishing mortality 
         tot_mort = fishing_mort + M) # and total mortality
head(fished_pop_setup)


# Calculate Fished Survival
fished_pop_setup[1, "fished_surv"] <- R0*prop_f # survival probability on day 1
for (AGE in 2:(max_age*12)){ # survival of day + 1 based on survival of the previous day.
  fished_pop_setup[AGE, "fished_surv"] <- fished_pop_setup[AGE-1, "fished_surv"]*exp(-fished_pop_setup[AGE-1, "tot_mort"]*step) # Divide by the time step here
}
head(fished_pop_setup)


# Calculate Mature Fished Biomass per Recruit
fished_pop_setup <- fished_pop_setup %>% 
  mutate(fished_bio = fished_mat * fished_surv * life_hist$weight^hyperallo)
head(fished_pop_setup)

# below shows the difference between the fished and unfished populations: 
# fished pop has lower biomass (makes sense if you're fishing)
# maturity (makes sense, it's not affected by fishing)
# survival is a steeper decline in the fished population (again makes sense if you're fishing) - slight but there is a diff.
matplot(
  x = fished_pop_setup$age, 
  y = fished_pop_setup[,c("selectivity", "fishing_mort", "tot_mort", "fished_surv", "fished_mat", "fished_bio")],
  type = "l", lty = 1, lwd = 2,
  col = 1:6,  # 6 columns being plotted
  xlab = "Age", ylab = "Value"
)
legend("right", legend = colnames(fished_pop_setup[c("selectivity", "fishing_mort", "tot_mort", "fished_surv", "fished_mat", "fished_bio")]),
       col = 1:ncol(fished_pop_setup), lty = 1, lwd = 2)


# plot check comparisons
plot(x = unfished_pop_setup$age, y = unfished_pop_setup$unfished_surv); points(x = fished_pop_setup$age, y = fished_pop_setup$fished_surv, col = "red", add = T)
legend("topright", legend = c("survival from unfished population", "survival from fished population"), col = c("black", "red"), lty = 1, lwd = 2)

plot(x = unfished_pop_setup$age, y = unfished_pop_setup$unfished_bio); points(x = fished_pop_setup$age, y = fished_pop_setup$fished_bio, col = "red", add = T)
legend("topright", legend = c("biomass from unfished population", "biomass from fished population"), col = c("black", "red"), lty = 1, lwd = 2)


# Calculate the Total Mature Fished Biomass per Recruit
fished_f_sb <- fished_pop_setup %>%  # This is the total mature female spawning biomass per recruit for the fished population
  slice(which(row_number() %% 12 == 1)) %>% # select december of each age (ages 1-40 in monthly steps)
  summarise(., sum = sum(fished_bio))



## Spawner per recruit and equilibrium recruitment ----------------------------

spr <- fished_f_sb/unfished_f_sb # this says 'under fished equilibrium, 1 recruit loses 26% (1-spr) of its reproductive output over her lifetime compared to unfished equilibrium'
equil_recr <- (fished_f_sb-alpha)/(beta*fished_f_sb) # this says 'however, under fished equilibrium, 1 recruit produces 97% of its unfished recruitment (because there is lower density)
# the density-dependence is visible here: fewer spawners = population under carrying capacity = more recruits


## STEP 3: set up the starting population -------------------------------------

## so far the fished and unfished populations and spawning biomasses were hypothetical, for a single recruit.
## here we scale the spawning biomasses to the level of recruitment, to make a starting population.

init_recr <- 5000 # in thousands - 4000 is too small

# calculate initial fished recruitment (how many new fish from the fished population)
init_fished_recr <- (fished_f_sb-alpha) / (fished_f_sb*beta) * init_recr 

# calculate initial unfished spawning biomass - this is the carrying capacity for the size of our model, defined by init_recr
init_unfished_f_sb <- init_fished_recr * unfished_f_sb

# calculate initial fished spawning biomass - this is the spawning biomass of our starting population
init_fished_f_sb <- init_fished_recr * fished_f_sb

# calculate new Beverton-Holt parameters - because we have a different initial recruitment value, alpha and beta will change.
alpha <- (init_unfished_f_sb / init_fished_recr) * ((1-h) / (4*h))
beta <- ((h-0.2) / (0.8*h*init_fished_recr))

# calculate number of female recruits in the next time step
n_f_recr <- (init_fished_f_sb/(alpha+beta*init_fished_f_sb))*prop_f # originally Charlotte had init_unfished_sb/(alpha...) but i think that's incorrect


# calculate survival into the next time step to create a full age structured population
starting_pop <- as.data.frame(array(1, dim=c((max_age*12), 1)))

starting_pop <- starting_pop %>% 
  mutate(age = life_hist$age) 

starting_pop[1,1] <- n_f_recr

for (r in 2:(max_age*12)){ # This calculates the survival in the next age based on the previous age using total mortality from our fished population
  starting_pop[r, 1] <- starting_pop[r-1, 1]*exp(-fished_pop_setup[r-1, 3]*step) # Divide by the time step here
}

starting_pop <- starting_pop %>% 
  rename(N = "V1")

## These fish form our starting population for the model
## At the end of the year the spawning stock biomass of all females will be calculated to generate recruitment for the next year
## Alpha and Beta need to be recorded for use in the next step of the model




## Selectivity-retention for the fish in the model ----------------------------

# Retention for each age group
# This is for when you actually run the model you don't use this for setting up the initial population
# Trying to account for the fact that fish that are below the legal size limit are likely to be thrown back and so won't necessarily die
fished_pop_setup <- fished_pop_setup %>% 
  mutate(
    retention_pre1913 = 1,
    retention1913 = ifelse(length <= 279, 0, ifelse(length > 279 & length < 279+5, 0.5, ifelse(length >= 279+5, 0.95, 0))),
    retention1977 = ifelse(length <= 380, 0, ifelse(length > 380 & length < 380+5, 0.5, ifelse(length >= 380+5, 0.95, 0))),
    retention1988 = ifelse(length <= 410, 0, ifelse(length > 410 & length < 410+5, 0.5, ifelse(length >= 410+5, 0.95, 0))),
    retention2008 = ifelse(length <= 450, 0, ifelse(length > 450 & length < 450+5, 0.5, ifelse(length >= 450+5, 0.95, 0))),
    retention2009 = ifelse(length <= 500, 0, ifelse(length > 500 & length < 500+5, 0.5, ifelse(length >= 500+5, 0.95, 0)))
    ) # account for minimum legal length changes

## Landings and discards
# This gives us the proportion of fish that are kept and the proportion that are thrown back
fished_pop_setup <- fished_pop_setup %>% 
  mutate(landings_pre1913 = selectivity * retention_pre1913,
         discards_pre1913 = selectivity*(1 - landings_pre1913)) %>% 
  
  mutate(landings_1913_1977 = selectivity * retention1913,
         discards_1913_1977 = selectivity*(1 - landings_1913_1977)) %>% 
  
  mutate(landings_1977_1988 = selectivity * retention1977,
         discards_1977_1988 = selectivity*(1 - landings_1977_1988)) %>% 
  
  mutate(landings_1988_2008 = selectivity * retention1988,
         discards_1988_2008 = selectivity*(1 - landings_1988_2008)) %>% 
  
  mutate(landings_2008_2009 = selectivity * retention2008,
         discards_2008_2009 = selectivity*(1 - landings_2008_2009)) %>%
  
  mutate(landings_post2009 = selectivity * retention2009,
         discards_post2009 = selectivity*(1 - landings_post2009))
  

## Selectivity-Retention Values Including post-release mortality
fished_pop_setup <- fished_pop_setup %>% 
  mutate(sel_ret_pre1913 = landings_pre1913 + (PRM * discards_pre1913),
         sel_ret_1913_1977 = landings_1913_1977 + (PRM * discards_1913_1977),
         sel_ret_1977_1988 = landings_1977_1988 + (PRM * discards_1977_1988),
         sel_ret_1988_2008 = landings_1988_2008 + (PRM * discards_1988_2008),
         sel_ret_2008_2009 = landings_2008_2009 + (PRM * discards_2008_2009),
         sel_ret_post2009 = landings_post2009 + (PRM * discards_post2009)
         )

sel_ret_pre1913 <- fished_pop_setup$sel_ret_pre1913
sel_ret_pre1913 <- array(sel_ret_pre1913, dim = c(12 , max_age)) # 12 months, 30 years old max
sel_ret_pre1913 <- t(sel_ret_pre1913)

sel_ret_1913_1977 <- fished_pop_setup$sel_ret_1913_1977
sel_ret_1913_1977 <- array(sel_ret_1913_1977, dim = c(12 , max_age)) # 12 months, 30 years old max
sel_ret_1913_1977 <- t(sel_ret_1913_1977)

sel_ret_1977_1988 <- fished_pop_setup$sel_ret_1977_1988
sel_ret_1977_1988 <- array(sel_ret_1977_1988, dim = c(12 , max_age)) # 12 months, 30 years old max
sel_ret_1977_1988 <- t(sel_ret_1977_1988)

sel_ret_1988_2008 <- fished_pop_setup$sel_ret_1988_2008
sel_ret_1988_2008 <- array(sel_ret_1988_2008, dim = c(12 , max_age)) # 12 months, 30 years old max
sel_ret_1988_2008 <- t(sel_ret_1988_2008)

sel_ret_2008_2009 <- fished_pop_setup$sel_ret_2008_2009
sel_ret_2008_2009 <- array(sel_ret_2008_2009, dim = c(12 , max_age)) # 12 months, 30 years old max
sel_ret_2008_2009 <- t(sel_ret_2008_2009)

sel_ret_post2009 <- fished_pop_setup$sel_ret_post2009
sel_ret_post2009 <- array(sel_ret_post2009, dim = c(12 , max_age)) # 12 months, 30 years old max
sel_ret_post2009 <- t(sel_ret_post2009)
# i interpret sel_ret as: how likely a fish is to be caught + kept from recruitment (1, 1) to the last month of its first year (1, 12) to the last month of its 30th year (30, 12)

# Then create a matrix of selectivity-retention, one slice per year.
sel_ret <- NULL

mll_change_1988_year <- 1988 - 1945
mll_change_2008_year <- 2008 - 1945
mll_change_2009_year <- 2009 - 1945
max_year <- year_end - year_start

for(i in 1:(mll_change_1988_year-1)){
  sel_ret <- abind(sel_ret, sel_ret_1977_1988, along = 3)
}

for(i in mll_change_1988_year:(mll_change_2008_year-1)){
  sel_ret <- abind(sel_ret, sel_ret_1988_2008, along = 3)
}

for(i in mll_change_2008_year:(mll_change_2009_year-1)){
  sel_ret <- abind(sel_ret, sel_ret_2008_2009, along = 3)
}

for(i in mll_change_2009_year:max_year){
  sel_ret <- abind(sel_ret, sel_ret_post2009, along = 3)
}

head(sel_ret)


## Saving files ---------------------------------------------------------------

# initial population
saveRDS(starting_pop, file = paste0("data/output_data/05_starting_population.rds"))

# selectivity
selectivity <- fished_pop_setup$selectivity
saveRDS(selectivity, file = "data/output_data/05_selectivity.rds")
saveRDS(sel_ret, file = "data/output_data/05_selectivity_retention.rds")

# maturity
maturity <- fished_pop_setup$fished_mat 
maturity <- array(maturity, dim = c(12, max_age))
maturity <- t(maturity)
saveRDS(maturity, file = "data/output_data/05_maturity.rds")

# weight of each age group
weight <- life_hist$weight
weight <- array(weight, dim = c(12, max_age))
weight <- t(weight)
saveRDS(weight, file="data/output_data/05_weight.rds")

# Beverton-Holt parameters
saveRDS(alpha, file = "data/output_data/05_Beverton-Holt_alpha.rds")
saveRDS(beta, file = "data/output_data/05_Beverton-Holt_beta.rds")


### END ###
