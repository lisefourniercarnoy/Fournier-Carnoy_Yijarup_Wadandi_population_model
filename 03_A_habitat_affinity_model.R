# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    BRUV data from GB and SW
# Task:    Figure out habitat affinities of mature and immature snapper
# Author:  Lise Fournier-Carnoy + Harry Carmody
# Date:    February 2026

# -----------------------------------------------------------------------------

# Status:  Did a good first try, need to see what soops says

# -----------------------------------------------------------------------------

## Load libraries -------------------------------------------------------------

rm(list = ls())

library(mgcv)
library(lme4)
library(MASS)
library(CheckEM) # obtain data from EM
library(tidyverse) # data handling
library(sf) # to handle spatial objects
library(glmmTMB) # to fit GLMM

## Files needed ---------------------------------------------------------------

file_gb_obs_habitat <- "data/input_data/habitat_affinity/D01_waatern_tidy_habitat.rds"
file_sw_obs_habitat <- "data/input_data/habitat_affinity/D01_waatu_tidy_habitat.rds"

file_gb_lengths <- "data/input_data/habitat_affinity/2024_waatern_commonwealth/2024-04_Geographe_stereo-BRUVs_Lengths.txt"
file_sw_lengths <- "data/input_data/habitat_affinity/2024_waatu_commonwealth/2024-10_SwC_stereo-BRUVs_Lengths.txt"

file_gb_metadata <- "data/input_data/habitat_affinity/2024_waatern_commonwealth/2024-04_Geographe_stereo-BRUVs_Metadata.csv"
file_sw_metadata <- "data/input_data/habitat_affinity/2024_waatu_commonwealth/2024-10_SwC_stereo-BRUVs_Metadata.csv"


## Load data ------------------------------------------------------------------

# Habitat
hab_gb <- readRDS(file_gb_obs_habitat) %>% 
  mutate(sand = sand/total_pts, # standardise habitat
         reef = reef/total_pts,
         seagrass = seagrass/total_pts
         ) %>% 
  glimpse()

hab_sw <- readRDS(file_sw_obs_habitat) %>% 
  mutate(sand = sand/total_pts, # standardise habitat
         reef = reef/total_pts,
         seagrass = seagrass/total_pts
  ) %>% 
  glimpse()

hab <- rbind(hab_gb, hab_sw)


# Lengths
length_gb <- read.table(
  file_gb_lengths,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  stringsAsFactors = FALSE
) %>% 
  dplyr::filter(Species == "auratus") %>% 
  group_by(OpCode) %>% 
  summarise(
    n_mature   = sum(Length > 660, na.rm = TRUE),
    n_immature = sum(Length <= 660, na.rm = TRUE),
    total   = n(),
    .groups = "drop"
  ) %>% 
  rename(opcode = OpCode) %>% 
  glimpse()

length_sw <- read.table(
  file_sw_lengths,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  stringsAsFactors = FALSE
) %>% 
  dplyr::filter(Species == "auratus") %>% 
  group_by(OpCode) %>% 
  summarise(
    n_mature   = sum(Length > 660, na.rm = TRUE),
    n_immature = sum(Length <= 660, na.rm = TRUE),
    total   = n(),
    .groups = "drop"
  ) %>% 
  rename(opcode = OpCode) %>% 
  glimpse()

length <- rbind(length_gb, length_sw)


# Add metadata back into it
meta_gb <- read.csv(file_gb_metadata, sep = ",", header = T, fill = T) %>% 
  dplyr::select(
    opcode,
    depth = depth_m
    ) %>% 
  glimpse()

meta_sw <- read.csv(file_sw_metadata, sep = ",", header = T, fill = T) %>% 
  dplyr::select(
    opcode,
    depth = depth_m
  ) %>% 
  glimpse()

metadata <- rbind(meta_gb, meta_sw)


# put it all together
glimpse(metadata)
glimpse(length)
dat <- left_join(metadata, hab, by = "opcode") %>% 
  dplyr::filter(!is.na(reef)) %>% # remove rows with missing habitat
  glimpse()
dat <- left_join(dat, length, by = "opcode") %>% 
  glimpse()

dat[is.na(dat)] <- 0 # fill NAs with true zeroes 

dat <- dat %>%
  pivot_longer(c(n_mature, n_immature)) %>%
  dplyr::rename(size_class = name, count = value) %>% 
  dplyr::select(!c("total")) %>% 
  mutate(location = ifelse(grepl("GB", opcode), "GB", "SW")) %>% 
  glimpse()

saveRDS(dat, file = "data/output_data/03_A_habitat_affinity_outputs/03_A_tidy_data.rds")

## Explore the data -----------------------------------------------------------

glimpse(dat)

## test
# following Harrison et al. 2018 - https://peerj.com/articles/4794/#p-36

glimpse(dat)

par(mfrow = c(2, 2))
plot(dat$count ~ dat$reef)
plot(dat$count ~ dat$sand)
plot(dat$count ~ dat$seagrass)

ggplot(dat, aes(x = reef, y = count)) + 
  geom_point(alpha = 0.2) +
  facet_grid(location~size_class) +
  geom_smooth(method = lm, se = T) +
  ggtitle("reef")
ggplot(dat, aes(x = sand, y = count)) + 
  geom_point(alpha = 0.2) +
  facet_grid(location~size_class) +
  geom_smooth(method = lm, se = T) +
  ggtitle("sand")
ggplot(dat, aes(x = seagrass, y = count)) + 
  geom_point(alpha = 0.2) +
  facet_grid(location~size_class) +
  geom_smooth(method = lm, se = T) +
  ggtitle("seagrass")

ggplot(dat, aes(x = depth, y = count)) + 
  geom_point(alpha = 0.2) +
  facet_grid(location~size_class) +
  geom_smooth(method = lm, se = T) +
  ggtitle("depth")


## Actually try some models ---------------------------------------------------

## the aim of the analysis is to understand the relative affinity of mature and immature snapper to different habitat
# dat$depth <- scale(dat$depth) # center + scale depth # i think this may be breaking the prediction later.
dat$depth2 <- dat$depth^2 # square depth

# check for correlation between predictors. nothing should be above 0.7.
cor <- dat %>%
  dplyr::select(depth, depth2, sand, reef, seagrass) %>%
  cor(use = "complete.obs")
corrplot::corrplot(cor, addCoef.col = T)
# we're all good here

### 1. try the Poisson distribution -------------------------------------------
glm_p <- glm(data = dat, # 'count is dependent on depth, size class, and habitat, and the effect of habitat depends on size class.'
             count ~ depth + depth2 +
               size_class +
               reef + sand +
               size_class:depth +
               size_class:reef + 
               size_class:sand, # not adding seagrass because otherwise habitat is super predictable and model doesn't run.
             family = poisson(link = "log")
)
summary(glm_p)

# overdispersion: conditional variance = conditional mean ?
mean(glm_p$fitted.values)
var(glm_p$residuals) # lol very much not. var>mean, overdispersion detected.
# this means that Poisson is not suitable. Ignoring this leads to inflated SEs, and therefore wrong p-values


### 2. try the Negative Binomial distribution ---------------------------------

glm_nb <- glm.nb(data = dat, # 'count is dependent on depth, size class, and habitat, and the effect of habitat depends on size class.'
                 count ~ depth + depth2 +
                   size_class +
                   reef + sand +
                   size_class:depth +
                   size_class:depth2 +
                   size_class:reef + 
                   size_class:sand # not adding seagrass because otherwise habitat is super predictable and model doesn't run.
)
summary(glm_nb) # notice Theta: Theta = 1/dispersion parameter. It's meant to capture the dispersion
1/glm_nb$theta # this is the dispersion parameter

# likelihood ratio test: is the nb model sig better than the poisson model
lmtest::lrtest(glm_p, glm_nb) # very much significant. this would justify the use of nb model.


# however the diagnostics look criminal so we'll try a zero-inflated model
par(mfrow = c(2, 2))
plot(glm_nb) #top left has 2 clouds (instead of a single shapeless cloud), qq-plot is not sitting on the line...

summary(glm_nb)
saveRDS(glm_nb, "data/output_data/03_A_habitat_affinity_outputs/03_A_model.rds")

# AS OF 02/03/2026 I AM USING THE ABOVE MODEL. 

### 3. test for zero-inflation ------------------------------------------------

# BRUV data may be zero-inflated (more zeroes than expected)
# we will test for zero-inflation by calculating how many zeroes our model expects
pred <- predict(glm_nb, type = "response")
theta <- sigma(glm_nb)
expected_zero_prob <- (theta / (theta + pred))^theta
sum(expected_zero_prob) # about 570 zeroes are expected
sum(dat$count == 0) # we have 600 zeroes ... it's a bit over what the model expects, let's test whether this is significant

sim_nb <- DHARMa::simulateResiduals(fittedModel = glm_nb)
plot(sim_nb)
DHARMa::testZeroInflation(sim_nb)  # p-value is non-significant, meaning that our data is not zero-inflated. 
# we can use a normal negative binomial model to model this data.


## below is archive, when i thought my data was overinflated

# # zero-inflated models allow you to specify the count model and the zero model. if unspecified, they have the exact same predictors.
# # for the count model, you need to think about what may influence abundance.
# # for the zero model, you need to think about what may cause structural zeroes (= in what conditions is it impossible to observe fish)

# # therefore let's try a zero-inflated model (try poisson first)
# glm_zip <- pscl::zeroinfl(data = dat, # 'count is dependent on depth, size class, and habitat, and the effect of habitat depends on size class.'
#                           count ~ depth + depth2 +
#                             size_class +
#                             reef + sand +
#                             size_class:depth +
#                             size_class:reef + 
#                             size_class:sand, # not adding seagrass because otherwise habitat is super predictable and model doesn't run.
#                           dist = "poisson", link = "logit"
# )
# summary(glm_zip) # it looks like 2 model outputs
# # count model coefficients (poisson with log link) is basically a normal glm (poisson) that tells you how much more counts are expected along each predictor.
# # zero-inflation model coefficients (binomial with logit link) the odds that the response takes the value of 0. these are odds ratios. not useful for us here, but they account for the zero-inflation.
# 
# # let's try a zero-inflated model (negative binomial) to see if it's better.
# glm_zinb <- pscl::zeroinfl(data = dat, # 'count is dependent on depth, size class, and habitat, and the effect of habitat depends on size class.'
#                            count ~ depth + depth2 +
#                              size_class +
#                              reef + sand + 
#                              size_class:depth +
#                              size_class:depth2 +
#                              size_class:reef + 
#                              size_class:sand, # not adding seagrass because otherwise habitat is super predictable and model doesn't run.
#                            dist = "negbin", link = "logit"
# )
# hist(glm_zinb$residuals) # residuals are not perfectly normally-distributed, but good enough for our purposes
# plot(glm_zinb$residuals~glm_zinb$fitted.values) # homogeneity of variance - looks a bit wonky but good enough for our purposes.
# 
# 
# summary(glm_zinb) # it looks like 2 model outputs
# # count model coefficients (poisson with log link) is basically a normal glm (poisson)
# # zero-inflation model coefficients (binomial with logit link) the odds that the response takes the value of 0. these are odds ratios.
# 
# # let's test whether the zero-inflated and non-zero-inflated are statistically the same
# pscl::vuong(glm_p, glm_zip) # here the 'Raw' p-value is highly significant, meaning the zi model is better than the non-zi
# 
# pscl::vuong(glm_nb, glm_zinb) # here the 'Raw' p-value is highly significant, meaning the zi model is better than the non-zi
# AIC(glm_nb, glm_zinb) # AIC of the zinb model is lower than non-zi model, so zinb model is better.
# 
# lmtest::lrtest(glm_zip, glm_zinb) # highly significant, the zinb is better. 
# 
# summary(glm_zinb)
# 
# # the exponentiated coefficients give us the IRR (incidence rate ratio) - i.e. how much more/less likely the expected count is compared to the reference
# exp(coef(glm_zinb))


### 4. Try a GLMM to account for location -------------------------------------

# there is a nested structure to the data (2 sampling sites: GB and SW), which violates the assumption of independence of errors.
# a GLMM accounts for this with random effects for location (by adding + (1|location)).

glmm_zinb <- glmmTMB(data = dat, # 'count is dependent on depth, size class, and habitat, and the effect of habitat depends on size class.'
                     count ~ depth + depth2 +
                       size_class +
                       reef + sand + 
                       size_class:depth +
                       size_class:depth2 +
                       size_class:reef + 
                       size_class:sand +
                       (1 | location), # the random effect accounting for nesting.
                     family = nbinom2,
                     ziformula = ~1) # zero inflation model
summary(glmm_zinb)
# the model is full of NaN - because 1. there isnt much difference in count between the locations - 2. there aren't enough sites (>5 recommended for GLMM)
# it could be an over-parameteristion problem but removing effects is not a good move here because we need all those covariates for the purposes of our model.

# we'll go back to our glm_nb, with the assumption of independence of errors kind of violated-but-not-too-much
# the large difference in habitat between GB and SW doesn't need to be accounted for in a random effect, because habitat (fixed effects) already capture this difference.



### END ###
