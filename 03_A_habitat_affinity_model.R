# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    BRUV data from GB and SW
# Task:    Figure out habitat affinities of mature and immature snapper
# Author:  Lise Fournier-Carnoy + Harry Carmody
# Date:    August 2026

# -----------------------------------------------------------------------------

# Status:  Finished, clean code, all good.

# -----------------------------------------------------------------------------

## Load libraries -------------------------------------------------------------

rm(list = ls())

library(tidyverse) # data handling
library(CheckEM) # obtain data from EM
library(sf) # to handle spatial objects
library(terra) # to handle spatial objectslibrary(mgcv)
library(lme4) # to model
library(MASS)
library(glmmTMB) # to fit GLMM (section 2.4)

## 0. Files needed ------------------------------------------------------------

file_gb_obs_habitat   <- "data/input_data/habitat_affinity/D01_waatern_tidy_habitat.rds"
file_sw_obs_habitat   <- "data/input_data/habitat_affinity/D01_waatu_tidy_habitat.rds"

file_gb_lengths       <- "data/input_data/habitat_affinity/2024-04_Geographe_stereo-BRUVs_Lengths.txt"
file_sw_lengths       <- "data/input_data/habitat_affinity/2024-10_SwC_stereo-BRUVs_Lengths.txt"

file_gb_maxn          <- "data/input_data/habitat_affinity/2024-04_Geographe_stereo-BRUVs_Points.txt"
file_sw_maxn          <- "data/input_data/habitat_affinity/2024-10_SwC_stereo-BRUVs_Points.txt"

file_gb_metadata      <- "data/input_data/habitat_affinity/2024-04_Geographe_stereo-BRUVs_Metadata.csv"
file_sw_metadata      <- "data/input_data/habitat_affinity/2024-10_SwC_stereo-BRUVs_Metadata.csv"

file_pred_hab         <- "data/output_data/01_A_bathymetry_habitat_rasters.rds"


## 1. Load data ---------------------------------------------------------------

# -- load habitat, lengths, metadata, and maxn, then remove false zeroes from the dataset

### 1.1 Habitat ---------------------------------------------------------------

# observed
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

hab <- rbind(hab_gb, hab_sw) %>% glimpse()

# predicted
pred_hab <- readRDS(file_pred_hab) %>% 
  glimpse()

# combine predicted and observed habitats
hab_sf <- st_as_sf(hab)
hab_vect <- vect(hab_sf)
test <- terra::extract(terra::unwrap(pred_hab), hab_vect)

test2 <- cbind(hab, test) %>% dplyr::select(-ID)
glimpse(test2)
hab <- test2


### 1.2 Metadata --------------------------------------------------------------

# add metadata back into it
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


### 1.3 MaxN counts -----------------------------------------------------------

maxn_gb <- read.table(
  file_gb_maxn,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  stringsAsFactors = FALSE
) %>% 
  dplyr::filter(Species == "auratus") %>% 
  mutate(opcode = OpCode) %>% 
  group_by(opcode, Frame) %>% 
  summarise(maxn = sum(Number)) %>% 
  glimpse() 

maxn_sw <- read.table(
  file_sw_maxn,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  stringsAsFactors = FALSE,
  quote = ""
) %>% 
  dplyr::filter(Species == "auratus") %>% 
  mutate(opcode = OpCode) %>% 
  group_by(opcode, Frame) %>% 
  summarise(maxn = sum(Number)) %>% 
  #dplyr::filter(!opcode %in% c('SWC-BV-016', "SWC-BV-084", "SWC-BV-131", 'SWC-BV-153')) %>% # these drops somehow don't have any fish in 
  glimpse()

maxn <- rbind(maxn_gb, maxn_sw) %>% 
  full_join(metadata, by = "opcode") %>% 
  mutate(maxn = replace_na(maxn, 0)) %>%  # fill true zeroes
  dplyr::select(-Frame) %>% 
  unique() %>% 
  group_by(opcode) %>% 
  slice_max(maxn, n = 1, with_ties = TRUE) %>%
  glimpse()


### 1.4 Lengths ---------------------------------------------------------------

length_split <- 375 
length_gb <- read.table(
  file_gb_lengths,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  stringsAsFactors = FALSE
) %>% 
  dplyr::filter(Species == "auratus",
                Length >= 171, 
                Length <= 1142) %>% # remove fish that are too big or too small to be correctly lengthed
  group_by(OpCode) %>% 
  summarise(
    n_mature   = sum(Length >= length_split, na.rm = TRUE),
    n_immature = sum(Length < length_split, na.rm = TRUE),
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
  dplyr::filter(Species == "auratus",
                Length >= 171, 
                Length <= 1142) %>% # remove fish that are too big or too small to be correctly lengthed
  group_by(OpCode) %>% 
  summarise(
    n_mature   = sum(Length > length_split, na.rm = TRUE),
    n_immature = sum(Length <= length_split, na.rm = TRUE),
    total   = n(),
    .groups = "drop"
  ) %>% 
  rename(opcode = OpCode) %>% 
  glimpse()

length_dat <- rbind(length_gb, length_sw)
length_dat <- full_join(metadata, length_dat, by = "opcode") %>% 
  glimpse()


### 1.5 Find false zeroes -----------------------------------------------------

length_zeroes <- unique(length_dat$opcode[is.na(length_dat$total)]) # opcodes where no fish were lengthed. calculate here because this is zeroes not split by mature/immature which is what we want
length_zeroes

true_zeroes <- unique(maxn$opcode[maxn$maxn == 0]) # opcodes with no fish counted
true_zeroes # opcodes where no fish were counted.

false_zeroes <- setdiff(length_zeroes, true_zeroes) # opcodes where fish were counted, but not lengthed
false_zeroes

# put it all together
glimpse(metadata)
glimpse(length_dat)
dat <- left_join(metadata, hab, by = "opcode") %>% 
  glimpse()
dat <- left_join(dat, length_dat %>% dplyr::select(-depth), by = "opcode") %>% 
  glimpse()

# remove false zeroes
dat <- dat %>% 
  dplyr::filter(!opcode %in% false_zeroes) %>% 
  glimpse()

# remove the couple drops that don't have habitat for some reason
dat <- dat %>% 
  dplyr::filter(!dplyr::if_any(c(sand, seagrass, reef, # remove drops where there is no habitat. otherwise the model removes them and outputs don't match with Harry's
                                 preef.fit, preef.se.fit,
                                 psand.fit, psand.se.fit,
                                 pseagrass.fit, pseagrass.se.fit), is.na))
dat[is.na(dat)] <- 0 # for drops where no snapper was seen, fill in with true zeroes.

# pivot longer so we can use the data in the model
dat_final <- dat %>%
  pivot_longer(c(n_mature, n_immature)) %>%
  dplyr::rename(size_class = name, count = value) %>% 
  dplyr::select(!c("total")) %>% 
  dplyr::filter(!opcode %in% c("GB-BV-182")) %>% # remove manually because doesn't have count data.
  glimpse()

saveRDS(dat, file = "data/output_data/03_A_habitat_affinity_outputs/03_A_tidy_data.rds")


## 2. Actually try some models ------------------------------------------------

## the aim of the analysis is to understand the relative affinity of mature and immature snapper to different habitat
dat_final$depth2 <- dat_final$depth^2 # square depth

# check for correlation between predictors. nothing should be above 0.7.
cor <- dat_final %>%
  dplyr::select(depth = depth, depth2, sand, reef, seagrass) %>%
  cor(use = "complete.obs")
corrplot::corrplot(cor, addCoef.col = T)
# we're all good here (nothing over 0.7 except for depth and depth2 which is okay)


### 2.1 try the Poisson distribution ------------------------------------------

glm_p <- glm(data = dat_final, # 'count is dependent on depth, size class, and habitat, and the effect of habitat depends on size class.'
             count ~ depth + depth2 +
               size_class +
               reef + sand +
               size_class:depth +
               size_class:reef + 
               size_class:sand, # not adding seagrass because otherwise habitat is super predictable and model doesn't run.
             family = poisson(link = "log")
)
summary(glm_p)

# overdispersion: does the conditional variance = conditional mean ?
mean(glm_p$fitted.values)
var(glm_p$residuals) # lol very much not. var>mean, overdispersion detected.
# this means that Poisson is not suitable. Ignoring this leads to inflated SEs, and therefore wrong p-values


### 2.2 Try the Negative Binomial distribution --------------------------------

# -- to model habitat affinity, we'll use the predicted habitat (from script 01_A) layers instead of the observed habitat.
# -- this is because the predicted habitat represents a broader area (250m resolution), whereas the observed habitat is only within 20m.
# -- by using predicted habitat we're trying to understand the macro-scale affinity, not the micro-scale affinity.

# try predicted habitat
glm_nb_p <- glm.nb(data = dat_final, # "count is dependent on depth, size class, and habitat, and the effect of habitat depends on size class."
                 count ~ depth + depth2 +
                   size_class +
                   preef.fit + psand.fit +
                   size_class:depth +
                   size_class:depth2 +
                   size_class:preef.fit + 
                   size_class:psand.fit, # not adding seagrass because otherwise habitat is super predictable and model doesn't run.
                 link = log
                 )

# have a look at the diagnostics
par(mfrow = c(2, 2))
plot(glm_nb_p) # that looks decent

print(summary(glm_nb_p))
saveRDS(glm_nb_p, "data/output_data/03_A_habitat_affinity_outputs/03_A_model.rds") # AS OF 06/03/2026 I AM USING THIS MODEL. 


### 2.3 Test for zero-inflation -----------------------------------------------

# -- BRUV data may be zero-inflated (more zeroes than expected)
# -- we will test for zero-inflation by calculating how many zeroes our model expects
pred <- predict(glm_nb_p, type = "response")
theta <- sigma(glm_nb_p)
expected_zero_prob <- (theta / (theta + pred))^theta
sum(expected_zero_prob) # about 500 zeroes are expected
sum(dat_final$count == 0) # we have 550 zeroes ... it's a bit over what the model expects, let's test whether this is significant

sim_nb <- DHARMa::simulateResiduals(fittedModel = glm_nb_p)
plot(sim_nb)
DHARMa::testZeroInflation(sim_nb)  # p-value is non-significant, meaning that our data is not zero-inflated. 
# -- we can use a normal negative binomial model to model this data.


## below is not needed if the model is not overinflated. 

# # -- zero-inflated models allow you to specify the count model and the zero model. if unspecified, they have the exact same predictors.
# # -- for the count model, you need to think about what may influence abundance.
# # -- for the zero model, you need to think about what may cause structural zeroes (= in what conditions is it impossible to observe fish)
# 
# # -- therefore let's try a zero-inflated model (try poisson first)
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
# # -- count model coefficients (poisson with log link) is basically a normal glm (poisson) that tells you how much more counts are expected along each predictor.
# # -- zero-inflation model coefficients (binomial with logit link) the odds that the response takes the value of 0. these are odds ratios. not useful for us here, but they account for the zero-inflation.
# 
# # -- let's try a zero-inflated model (negative binomial) to see if it's better.
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
# # -- count model coefficients (poisson with log link) is basically a normal glm (poisson)
# # -- zero-inflation model coefficients (binomial with logit link) the odds that the response takes the value of 0. these are odds ratios.
# 
# # -- let's test whether the zero-inflated and non-zero-inflated are statistically the same
# pscl::vuong(glm_p, glm_zip) # here the 'Raw' p-value is highly significant, meaning the zi model is better than the non-zi
# 
# pscl::vuong(glm_nb, glm_zinb) # here the 'Raw' p-value is highly significant, meaning the zi model is better than the non-zi
# AIC(glm_nb, glm_zinb) # AIC of the zinb model is lower than non-zi model, so zinb model is better.
# 
# lmtest::lrtest(glm_zip, glm_zinb) # highly significant, the zinb is better.
# 
# summary(glm_zinb)
# 
# # -- the exponentiated coefficients give us the IRR (incidence rate ratio) - i.e. how much more/less likely the expected count is compared to the reference
# exp(coef(glm_zinb))


### 2.4 Try a GLMM to account for location ------------------------------------

# -- there is a nested structure to the data (2 sampling sites: GB and SW), which violates the assumption of independence of errors.
# -- a GLMM accounts for this with random effects for location (by adding + (1 | location)).

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
# -- the model is full of NaN - because 1. there isnt much difference in count between the locations - 2. there aren't enough sites (>5 recommended for GLMM)
# -- it could be an over-parameteristion problem but removing effects is not a good move here because we need all those covariates for the purposes of our model.

# -- we'll go back to our glm_nb, with the assumption of independence of errors kind of violated-but-not-too-much
# -- the large difference in habitat between GB and SW doesn't need to be accounted for in a random effect, because habitat (fixed effects) already capture this difference.


### END ###
