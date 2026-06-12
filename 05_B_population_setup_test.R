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
a <- (0.160 + 0.168)/2 # "constant relative rate of relative growth rate", Schnute 1981, average of females and males
b <- (1.239 + 1.206)/2 # "incremental relative rate of relative growth rate", Schnute 1981, average of females and males

max_length <- 1095 # in mm, table 4.2 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr#page=12.15

# age-length standard deviation parameters, for equation 3 in Francis et al. 2016, https://www.sciencedirect.com/science/article/pii/S0165783615000909
al_a <- 20 # manually selected to kinda fit figure 6 in Wakefield et al. 2017  https://academic.oup.com/icesjms/article/74/1/180/2669555?login=true&guestAccessKey=
al_b <- 0.1 # manually selected to kinda fit figure 6 in Wakefield et al. 2017  https://academic.oup.com/icesjms/article/74/1/180/2669555?login=true&guestAccessKey=


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
lengths <- (10:(max_length + 100)) # all lengths starting at 50mm (see here for post settlement juvenile size ref https://rsnz.onlinelibrary.wiley.com/doi/10.1080/00288330.2014.892013)(+ a little bit to encompass extra big fish just in case)
lengths_1cm <- lengths[seq(1, length(lengths), 10)] # 1cm bins
length_age_matrix <- matrix(0,
                            nrow = length(mean_life_hist$age),
                            ncol = length(lengths_1cm), # 1cm length bins
                            dimnames = list(as.numeric(mean_life_hist$age), 
                                            as.numeric(lengths_1cm)))

for (i in 1:nrow(length_age_matrix)) {
  
  mean_length = round(mean_life_hist$length[i])
  sd_length = al_a + al_b * mean_length # eq. 3 from Francis et al. 2016, https://www.sciencedirect.com/science/article/pii/S0165783615000909
  
  length_age_matrix[i, ] <- dnorm(lengths_1cm, mean = mean_length, sd = sd_length)
} # this loops calculates how likely it is that a fish of age x is of each length y

length_age_matrix <- length_age_matrix / rowSums(length_age_matrix) # normalise so it all adds up to 1

dimnames(length_age_matrix) <- list(mean_life_hist$age, lengths_1cm)
summary(rowSums(length_age_matrix))  # should all be 1

# okay so this matrix is going to be the base on which most things get calculated downstream.
library(ggridges)

# sanity check : the distribution of lengths at various ages
as.data.frame(length_age_matrix) %>% 
  mutate(age = as.numeric(rownames(length_age_matrix))) %>% 
  pivot_longer(cols = -age, names_to = "length_bin", values_to = "density") %>% 
  mutate(length_bin = as.numeric(length_bin)) %>% 
  filter(age %in% as.numeric(rownames(length_age_matrix))[seq(1, nrow(length_age_matrix), by = 12)]) %>% 
  group_by(age) %>% 
  mutate(density_scaled = density / max(density)) %>% # normalise each age to peak = 1
  ungroup() %>% 
  
  ggplot(aes(x = length_bin, y = factor(age), height = density_scaled, group = age)) +
  geom_ridgeline(scale = 2, alpha = 0.8, linewidth = 0.3) +
  labs(x = "Length", y = "Age") +
  theme_ridges(grid = TRUE)


## STEP 1.5: make a length transition matrix ==================================

# the length age matrix tells us the probability of a fish of age x to be of length y, but that doesnt tell us what length a fish of age x and length y will be at age x+1
# for this we need an age-transition matrix, of format length x length.
# age transition matrices are usually calculated with tag-recapture data. 
# i don't have that data, so i'll simulate it.

### A. release the fish! ----

set.seed(123)

# get some functions sorted
schnute_length <- function(AGE, l1, l2, t1, t2, a, b) {
    L <- (l1^b + (l2^b - l1^b) * (1 - exp(-a * (AGE - t1))) / (1 - exp(-a * (t2 - t1))))^(1 / b)
  pmin(L, max_length)
} # this function is to obtain the expected length at AGE years old 

sd_at_length <- function(mu_length, al_a, al_b) {
  al_a + al_b * mu_length
} # this is to add some variability, eq.3 in Francis 2016

# some mean life history
ages = as.numeric(rownames(length_age_matrix))
lengths = as.numeric(colnames(length_age_matrix))

# catch, tag, and release some fish!
n_fish  <- 100000

tag_recap <- data.frame(
  release_length  = sample(lengths, size = n_fish, replace = TRUE)
  ) %>% 
  glimpse()

### B. pretend you know the age of these released fish ----

mu_by_age <- schnute_length(ages, l1, l2, t1, t2, a, b) # this is the mean length for all ages that a fish can be (monthly)
sd_by_age <- sd_at_length(mu_by_age, al_a, al_b) # this is the variability around the mean length, for all ages that a fish can be

age_likelihoods <- function(LENGTH) {
  dnorm(LENGTH, mean = mu_by_age, sd = sd_by_age)
} # this gives us the probability that a fish of a length x is of each age y.

tag_recap$release_age <- sapply(tag_recap$release_length, function(LENGTH) {
  len_capped <- min(LENGTH, max(lengths_1cm)) # prevent fish from being larger than the limit
  liks <- age_likelihoods(len_capped) # finds the likely ages for that length
  liks <- liks / sum(liks) # normalise to get weights
  sum(ages * liks) # get the weighted average age
}) # this obtains the age (as an average from the likely ages) for every fish we released

head(tag_recap)


### C. one month later, catch the fish again! ----

# one year later, we go back and capture the fish again
tag_recap$recapture_age <- tag_recap$release_age + (1/12) # they have all aged 1 month

# their lengths follow Schnute again, based on their recapture age (with some variation)
mu_recap <- schnute_length(tag_recap$recapture_age, l1, l2, t1, t2, a, b)
sd_recap <- sd_at_length(mu_recap, al_a, al_b)

# we have to randomly decide what length the fish will be one year after release, according to mu & sd,
# keeping in mind that fish can't shrink.
library(truncnorm)
# so rtruncnorm will do "okay this fish was release at __mm, one year later its length will be randomly selected from a normal distribution that makes sense for its age (mu & sd), keeping in mind the min length (release length) and max length."
tag_recap$recapture_length <- rtruncnorm(
  n     = n_fish,
  a     = tag_recap$release_length, # lower bound (no shrinkage)
  b     = max(lengths_1cm),# upper bound
  mean  = mu_recap,
  sd    = sd_recap
) # boom shakalaka

# for very big fish, the rtruncnorm can't calculate a growth (because it's very close to the max length) so just assign them the release length.
tag_recap$recapture_length <- ifelse(is.na(tag_recap$recapture_length), tag_recap$release_length, tag_recap$recapture_length) 

tag_recap$growth <- tag_recap$recapture_length - tag_recap$release_length

head(tag_recap)
summary(tag_recap)

# cool beans, we now have a fully structured tag-recapture dummy dataset to base our age-transition matrix on.

# check it over
ggplot(tag_recap, aes(x = release_length, y = recapture_length)) +
  geom_bin2d(bins = 60) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  labs(x = "Release Length (mm)", y = "Recapture Length (mm)", fill = "Count") +
  theme_minimal()

# here we have quick and homogenous growth for little fish (hurry up otherwise you're prey), then the medium fish can grow fast or not grow fast, and the big fish can't grow bigger than the max
tag_recap %>%
  filter(!is.na(growth)) %>%
  mutate(length_bin = factor(
    paste0(floor(release_length / 50) * 50, "–", floor(release_length / 50) * 50 + 49),
    levels = paste0(
      seq(floor(min(release_length, na.rm=TRUE)/50)*50,
          floor(max(release_length, na.rm=TRUE)/50)*50, by=50),
      "–",
      seq(floor(min(release_length, na.rm=TRUE)/50)*50,
          floor(max(release_length, na.rm=TRUE)/50)*50, by=50) + 49
    )
  )) %>%
  ggplot(aes(x = growth, y = length_bin)) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  labs(x = "growth (mm)", y = "release length (mm)") +
  theme_ridges()


### XX okay now the length transition matrix ----

glimpse(tag_recap)

# define your length bins (e.g. 50mm bins)
bin_breaks <- seq(5, max_length+100, by = 10) # +100 for extra big fish
bin_labels <- seq_along(bin_breaks[-1])  # indices 1 to n_bins

bin_mids <- bin_breaks[-length(bin_breaks)] + 5  # midpoints
rownames(ATM) <- bin_mids
colnames(ATM) <- bin_mids

n_bins <- length(bin_labels)
ATM <- matrix(0, nrow = n_bins, ncol = n_bins)

# assign each fish to release and recapture bins
tag_recap$release_bin <- cut(tag_recap$release_length, 
                             breaks = bin_breaks, labels = FALSE)
tag_recap$recapture_bin <- cut(tag_recap$recapture_length, 
                               breaks = bin_breaks, labels = FALSE)


for (i in 1:n_bins) {
  current_length_bin <- bin_breaks[i] + 5 # bin midpoint
  fish <- tag_recap[tag_recap$release_bin == i, ] # find all the fish that were released at this length bin.
  
  mu_recap <- mean(fish$recapture_length) # find their average recapture length
  sigma_recap <- sd(fish$recapture_length) # and the SD of the recapture length
  
  # for each bin i, what is the probability of being that length, given the recapture data (p(length bin) is p(upper bound) - p(lower bound))
  probs <- pnorm(bin_breaks[-1], mu_recap, sigma_recap) - pnorm(bin_breaks[-length(bin_breaks)], mu_recap, sigma_recap)
  
  probs[bin_breaks[-length(bin_breaks)] < current_length_bin] <- 0  # no shrinkage, zero out smaller bins.
  
  ATM[i,] <- probs / sum(probs) # scale so it all adds up to 1
}

head(ATM) # and that's the age transition matrix!

# rename for clarity
colnames(ATM) <- lengths
rownames(ATM) <- lengths


saveRDS(ATM, "data/output_data/05_age_transition_matrix.rds")


## STEP 2: set up a hypothetical unfished population ==========================

# the goal here is to calculate what the distribution of mature biomass looks like for a hypothetical population.
# what this means is that for a hypothetical recruitment of R0 = 1, where is the bulk of biomass across lengths and ages?
# this whole setup is made to obtain a single value: unfished_female_spawning_biomass.
# we use this downstream.

# calculate maturity and mature biomass across all lengths
length_mature_biomass <- data.frame(
  length = lengths_1cm
  ) %>% 
  mutate(
  weight = WLa * (length^WLb) / 1000
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
age_length_starting_pop <- sweep(length_age_matrix, MARGIN = 1, age_starting_pop[, "n"], FUN = "*")

## These fish form our starting population for the model
## At the end of the year the spawning stock biomass of all females will be calculated to generate recruitment for the next year
## Alpha and Beta need to be recorded for use in the next step of the model

## Selectivity-retention for the fish in the model ----------------------------

# we want to make a matrix that tells us for every month of every year, how likely is it that a fish of length x gets caught and kept by a fisher.

# first we translate selectivity to a length-based vector, so that the selectivity-retention has the same number of values to calculate things on.
fished_pop_setup$selectivity

# we'll use the selectivity-at-age to approximate a selectivity-at-length
sel_at_age <- data.frame(age = mean_life_hist$age, selectivity = fished_pop_setup$selectivity)
sel_at_length <- data.frame(
  length     = lengths_1cm,
  selectivity = NA
)

sel_at_age <- approx(
  x = mean_life_hist$age,
  y = fished_pop_setup$selectivity,
  xout = mean_life_hist$age,
  rule = 2
)$y

length_to_sel <- approxfun( # linear interpolation between selectivity-at-age and -at-length
  x = mean_life_hist$length,
  y = sel_at_age,
  rule = 2
)

sel_at_length$selectivity <- length_to_sel(as.numeric(sel_at_length$length) + 5) # apply the interpolation to the length bins

# manually fill bins that are too-big and too-small (they're not observed in the sel-at-age so need to be manually filled)
sel_at_length$selectivity[as.numeric(sel_at_length$length) + 5 < min(mean_life_hist$length)] <- 0
sel_at_length$selectivity[as.numeric(sel_at_length$length) + 5 > max(mean_life_hist$length)] <- 1



## okay now to calculate things
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
  landings = sel_at_length$selectivity * ret
  # discards: the % of fish of each length that can be caught but will be chucked back
  discards = sel_at_length$selectivity * (1 - ret)
  
  # selectivity-retention: the % of fish of each age (in each year) that will die from fishing activities (caught, and either kept, or chucked back and die)
  sel_ret[, as.character(YEAR)] = landings + (PRM * discards)
  
}
head(sel_ret)


## Saving files ---------------------------------------------------------------

# initial population
saveRDS(age_length_starting_pop, file = paste0("data/output_data/05_starting_population.rds"))

# number of length bins
lengths
saveRDS(lengths, file = "data/output_data/05_length_bins.rds")

# selectivity
selectivity <- fished_pop_setup$selectivity
saveRDS(selectivity, file = "data/output_data/05_selectivity.rds")
saveRDS(sel_ret, file = "data/output_data/05_selectivity_retention.rds")

# maturity
maturity <- length_mature_biomass$maturity
saveRDS(maturity, file = "data/output_data/05_maturity.rds")

# weight of each age group
weight <- mean_life_hist$weight
saveRDS(weight, file = "data/output_data/05_weight.rds")

# Beverton-Holt parameters
saveRDS(alpha, file = "data/output_data/05_Beverton-Holt_alpha.rds")
saveRDS(beta, file = "data/output_data/05_Beverton-Holt_beta.rds")

### END ###
