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
dat$depth <- scale(dat$depth) # center + scale depth
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
# this means that Poisson is not suitable. Ignoring thi leds to inflated SEs, and therefore wrong p-values


### 2. try the Negative Binomial distribution ---------------------------------

glm_nb <- glm.nb(data = dat, # 'count is dependent on depth, size class, and habitat, and the effect of habitat depends on size class.'
                 count ~ depth + depth2 +
                   size_class +
                   reef + sand +
                   size_class:depth +
                   size_class:reef + 
                   size_class:sand # not adding seagrass because otherwise habitat is super predictable and model doesn't run.
)
summary(glm_nb) # notice Theta: Theta = 1/dispersion parameter. It's meant to capture the dispersion
1/glm_nb$theta # this is the dispersion parameter

# likelihood ratio test: is the nb model sig better than the poisson model
lmtest::lrtest(glm_p, glm_nb) # very much significant. this would justify the use of nb model.


# however the diagnostics look criminal so we'll try a zero-inflated model
plot(glm_nb) #top left has 2 clouds (instead of a single shapeless cloud), qq-plot is not sitting on the line...


### 3. try a zero-inflated model ----------------------------------------------

# we have loads of zeroes
(sum(dat$count == 0)/nrow(dat))*100 # 72% of zeroes... that's wayyy too much to ignore.
pr <- predict(glm_nb, type = "response")
exp(-1*mean(pr))*500 # we would expect around 230 zeroes for our dataset but there's...
sum(dat$count == 0) # close to 650...

# therefore let's try a zero-inflated model (try poisson first)
glm_zip <- pscl::zeroinfl(data = dat, # 'count is dependent on depth, size class, and habitat, and the effect of habitat depends on size class.'
                          count ~ depth + depth2 +
                            size_class +
                            reef + sand +
                            size_class:depth +
                            size_class:reef + 
                            size_class:sand, # not adding seagrass because otherwise habitat is super predictable and model doesn't run.
                          dist = "poisson", link = "logit"
)
summary(glm_zip) # it looks like 2 model outputs
# count model coefficients (poisson with log link) is basically a normal glm (poisson) that tells you how much more counts are expected along each predictor.
# zero-inflation model coefficients (binomial with logit link) the odds that the response takes the value of 0. these are odds ratios. not useful for us here, but they account for the zero-inflation.

# let's try a zero-inflated model (negative binomial) to see if it's better.
glm_zinb <- pscl::zeroinfl(data = dat, # 'count is dependent on depth, size class, and habitat, and the effect of habitat depends on size class.'
                           count ~ depth + depth2 +
                             size_class +
                             reef + sand + 
                             size_class:depth +
                             size_class:reef + 
                             size_class:sand, # not adding seagrass because otherwise habitat is super predictable and model doesn't run.
                           dist = "negbin", link = "logit"
)
hist(glm_zinb$residuals) # residuals are not perfectly normally-distributed, but good enough for our purposes
plot(glm_zinb$residuals~glm_zinb$fitted.values) # homogeneity of variance - looks a bit wonky but good enough for our purposes.


summary(glm_zinb) # it looks like 2 model outputs
# count model coefficients (poisson with log link) is basically a normal glm (poisson)
# zero-inflation model coefficients (binomial with logit link) the odds that the response takes the value of 0. these are odds ratios.

# let's test whether the zero-inflated and non-zero-inflated are statistically the same
pscl::vuong(glm_p, glm_zip) # here the 'Raw' p-value is highly significant, meaning the zi model is better than the non-zi

pscl::vuong(glm_nb, glm_zinb) # here the 'Raw' p-value is highly significant, meaning the zi model is better than the non-zi

lmtest::lrtest(glm_zip, glm_zinb) # highly significant, the zinb is better. 

summary(glm_zinb)

# the exponentiated coefficients give us the IRR (incidence rate ratio) - i.e. how much more/less likely the expected count is compared to the reference
exp(coef(glm_zinb))



### END ###


## below is archive ----

## Remove false zeroes in lengths ---------------------------------------------

# Not all fish that are MaxN'ed will be lengthed. Sometimes MaxN is >0, but no fish are lengthed. We need to filter those out.

test <- dat %>% 
  mutate( # where MaxN is 0 (true 0), length counts are also 0, not NA
    immature_snapper = ifelse(Total_count == 0, 0, immature_snapper),
    mature_snapper = ifelse(Total_count == 0, 0, mature_snapper)
  ) %>% 
  dplyr::filter(!is.na(immature_snapper)) %>% # all other NAs need to be filtered out.
  dplyr::select(-c(sample_url))

dat_op1 <- test %>%
  dplyr::rename(response = Total_count) %>%
  dplyr::mutate(depth_c = as.numeric(scale(depth_m, center = TRUE, scale = FALSE)),
                depth2  = depth_c^2, 
                Location = case_when(
                  grepl("2024-04_Geographe_stereo-BRUVs", campaignid) ~ "Geographe Bay", 
                  grepl("2020-06_south-west_stereo-BRUVs", campaignid) ~ "South-West", 
                  grepl("2020-10_south-west_stereo-BRUVs", campaignid) ~"South-West",
                  grepl("2023-03_SwC_stereo-BRUVs", campaignid) ~"South-West",
                  TRUE ~ NA_character_
                ),
                Location = factor(Location))



tidy_all_length <- test %>%
  pivot_longer(c(immature_snapper, mature_snapper)) %>%
  dplyr::rename(size_class = name, count = value) %>% 
  glimpse()


dat_op2 <- tidy_all_length %>%
  dplyr::rename(response = count) %>%
  dplyr::mutate(depth_c = as.numeric(scale(depth_m, center = TRUE, scale = FALSE)),
                depth2  = depth_c^2, 
                Location = case_when(
                  grepl("2024-04_Geographe_stereo-BRUVs", campaignid) ~ "Geographe Bay", 
                  grepl("2020-06_south-west_stereo-BRUVs", campaignid) ~ "South-West", 
                  grepl("2020-10_south-west_stereo-BRUVs", campaignid) ~"South-West",
                  grepl("2023-03_SwC_stereo-BRUVs", campaignid) ~"South-West",
                  grepl("2024-10_SwC_stereo-BRUVs", campaignid) ~"South-West",
                  
                  TRUE ~ NA_character_
                ),
                Location = factor(Location),
                size_class = factor(size_class))


## Option 1: does MaxN abundance depend on habitat/depth? ---------------------

glimpse(dat_op1)

# Charlotte says: comparing option1 and option2 is not reasonable, (predictors dont predict for differences in measurability)
# Charlotte also says: it's okay for model assumptions to be a bit shaky given the sample size and the fact that it's ecological data.

## Option 3: does the habitat affinity depend on fish size ? ------------------

glimpse(dat_op2)


# first try: lm()

test1 <- lm(data = dat_op2,
            response ~ depth_c + depth2 + 
              size_class * (preef.fit + psand.fit + pseagrass.fit) + 
              Location)

# assumption 1: independence of samples (yes)

# assumption 2: normality of residuals
hist(test1$residuals) # looks a bit skewed, but not horrible

# assumption 3: homoscedasticity
plot(test1$residuals ~ test1$fitted.values) # looks a bit skewed, but not horrible
par(mfrow=c(2,2))
plot(test1)

# conclusion: assumptions are meh-validated, maybe a glm() would be better

library(statmod)

# second try: glm()
test2 <- glm(data = dat_op2,
             response ~ depth_c + depth2 + 
               size_class * (preef.fit + psand.fit + pseagrass.fit) + 
               Location,
             #family = tweedie(link.power=0, var.power = 1)
             family = quasipoisson(link = "log") # potential? 
             #family = poisson(link = "log") # not suitable because data is overdispersed (variance > mean)
             #family = gaussian(link = "identity) # not suitable because our data is not continuous
             #family = Gamma(link = "inverse) # not suitable because our data is not continuous
             #family = binomial(link = "logit") # not suitable because our data is not presence/absence
             )

test2 <- glm.nb(data = dat_op2,
                response ~ depth_c + depth2 + 
                  size_class * (preef.fit + psand.fit + pseagrass.fit) + 
                  Location,
                )

# assumption 1: independence of samples (yes)

# assumption 2: normality of residuals
hist(test2$residuals) # looks a bit skewed, but not horrible

# assumption 3: homoscedasticity
plot(test2$residuals ~ test2$fitted.values) # looks a bit skewed, but not horrible

par(mfrow=c(2,2))
plot(test2)

summary(test2)



# test3: simpler glm
test3 <- glm(data = dat_op2,
             response ~ depth_c + depth2 + pseagrass.fit,
             #family = tweedie(link.power=0, var.power = 1)
             family = quasipoisson(link = "log") # potential? 
             #family = poisson(link = "log") # not suitable because data is overdispersed (variance > mean)
             #family = gaussian(link = "identity) # not suitable because our data is not continuous
             #family = Gamma(link = "inverse) # not suitable because our data is not continuous
             #family = binomial(link = "logit") # not suitable because our data is not presence/absence
)

hist(test3$residuals)
summary(test3)


# glmm is not for us, because no random effect (maybe location) ??




### END ###
