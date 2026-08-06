# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Values from the literature
# Task:    Set up populations
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    June 2026

# -----------------------------------------------------------------------------

# Status: First checks

# -----------------------------------------------------------------------------

rm(list = ls())

# Load libraries

library(tidyverse)
library(RColorBrewer)

# fishing mortality F = fishing effort E * catchability Q. 
# F is the proportion of the population which is harvested. 

# E is findable (data from papers, see historical reconstructions in 04_A to 04_C, section 4)
# F is estimated in the stock assessment (see Fig. 3.14 in https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1001&context=fish_rar)
# Q is the only one that is unknowable, and therefore calibrate-able.

## read in files --------------------------------------------------------------

water   <- readRDS("data/output_data/02_watergrid.rds")
dat_com <- readRDS("data/output_data/04_A_commercial_fishing_info.rds")
dat_brec <- readRDS("data/output_data/04_C_boat_rec_fishing_info.rds")
dat_srec <- readRDS("data/output_data/04_B_shore_rec_fishing_info.rds")
names(dat_com)

## target F from stock assessment ---------------------------------------------

full_years <- 1975:2024
F_ss <- read.csv("data/input_data/digitised_plots_for_checking/F_digitised_from_stock_assessment.csv") %>%
  glimpse()
F_ss <- data.frame(
  x = full_years,
  y = approx(x = F_ss$x, y = F_ss$y, xout = full_years)$y
)
plot(F_ss, type = "l", lwd = 2, col = "steelblue")

## COMMERCIAL fishing mortality -----------------------------------------------

# obtain 1 fishing effort value per year
c_effort <- apply(dat_com$fishing_days, 3, sum)  # vector, length 125
plot(c_effort[75:125], type = "l")

# catchability likely changed over time as better gear and better tech allowed fishers to fish more efficiently
# see Marriott et al. 2010 for details https://academic.oup.com/icesjms/article/68/1/76/631662?login=false 
# we are going to try to recreate the commercials catch curve from Figure 3.8.d in Fisher et al. 2025, by adjusting the catchability.

# in Marriott et al. 2010, table 3, we take the mean of the GPS and colour sounder values for a logistic curve.
# we then fit an 'tech adoption' curve as they did in the paper
alpha = mean(c(1, 0.918)) # % of all fishers who adopt the tech, mean of values in paper
beta = mean(c(-0.531, -0.493)) # shape of curve, mean of values in paper
delta = mean(c(3.821, 2.251)) + 88 # x midpoint of curve, mean of values in paper + 86 (since the midpoint of fig 2a and 2b is 1984 and 1990)

logistic_fn <- function(y) {
  alpha / (1 + exp(beta * (y - delta)))
}
y <- 1:125
Py <- logistic_fn(y)
plot(1900:2024, Py, type = "l", ylab = "Proportion adopted", ylim = c(0, 1))


# in table 2 we are given the mean efficiency increase of fishers when they use GPS and colour sounders.
# we will use these to get a pattern of Q over time for commercial fishers. i've increased them a bit to  make the catchability curve fit.
# we'll also use a base Q (catchability) value - to be calibrated to obtain a sensical fishing mortality
base_Q <- 5e-06
Q_inc = (mean(c(133.2, 39.2))/100) + 5
Q_commercial <- data.frame(
  year = 1900:2024, 
  Q = base_Q * (1 + (Py / alpha) * Q_inc)
)

plot(Q_commercial$year, Q_commercial$Q, type = "l") # catchability here increases logistically from the 80s. 

F_commercial <- data.frame(
  year = 1900:2024, 
  F = Q_commercial$Q * c_effort
)

plot(F_commercial$year[75:125], F_commercial$F[75:125], type = "l") # this should somewhat match the commercial catch. see Fisher et al. 2025, figure 3.8d


## BOAT REC fishing mortality -------------------------------------------------

# obtain 1 fishing effort value per year
b_effort <- apply(dat_brec$fishing_days, 3, sum)  # vector, length 125


# catchability likely changed over time as better gear and better tech allowed fishers to fish more efficiently
# see Marriott et al. 2010 for details https://academic.oup.com/icesjms/article/68/1/76/631662?login=false 
# in table 3, we take the mean of the GPS and colour sounder values for a logistic curve.
# we then fit an 'tech adoption' curve as they did in the paper
alpha = mean(c(1, 0.918)) # mean of values in paper
beta = mean(c(-0.531, -0.493)) # mean of values in paper
delta = mean(c(3.821, 2.251)) + 95 # mean of values in paper + 86 (since the midpoint of fig 2a and 2b is 1984 and 1990)

logistic_fn <- function(y) {
  alpha / (1 + exp(beta * (y - delta)))
}
y <- 1:125
Py <- logistic_fn(y)
plot(1900:2024, Py, type = "l", ylab = "Proportion adopted", ylim = c(0, 1))

# in table 2 we are given the mean efficiency increase of fishers when they use GPS and colour sounders.
# we will use these to get a pattern of Q over time for commercial fishers
# we'll also use a base Q (catchability) value - to be calibrated to obtain a sensical fishing mortality
base_Q <- 5e-08
Q_inc = mean(c(133.2, 39.2))/100
Q_boat_rec <- data.frame(
  year = 1900:2024, 
  Q = base_Q * (1 + Py * Q_inc)
)

plot(x = Q_boat_rec$year, Q_boat_rec$Q, type = "l") # catchability here increases logistically from the 80s. 

F_boat_rec <- data.frame(
  year = 1900:2024, 
  F = Q_boat_rec$Q * b_effort
)

plot(F_boat_rec$year[75:125], F_boat_rec$F[75:125], type = "l") # this should somewhat match the commercial catch. see Fisher et al. 2025, figure 3.8d


## SHORE REC effort -----------------------------------------------------------

s_effort <- apply(dat_srec$fishing_days, 3, sum)  # vector, length 125

start_Q <- 4e-08
Q_shore_rec <- data.frame(year = 1900:2024, 
                          Q = start_Q) %>%
  glimpse()

for (Y in 1950:2024) {
  Q_shore_rec$Q[Q_shore_rec$year == Y] <- 1.02 * Q_shore_rec$Q[Q_shore_rec$year == Y-1]
}


F_shore_rec <- data.frame(
  year = 1900:2024, 
  F = Q_shore_rec$Q * s_effort
)
plot(F_shore_rec[75:125,], type = "l", col = "steelblue", lwd = 2)


## see overall ----------------------------------------------------------------

ggplot() + # this should somewhat match Fisher et al. 2025, figure 3.8d, shore rec excluded
  geom_line(data = F_commercial[75:125,], aes(x = year, y = F), linewidth = 1.5, linetype = "dashed") +
  geom_line(data = F_boat_rec[75:125,], aes(x = year, y = F), linewidth = 1.5, linetype = "dotted") +
  geom_line(data = F_shore_rec[75:125,], aes(x = year, y = F), linewidth = 1.5, linetype = "solid") +
  
  theme_minimal()

# below, the red line should align on the blue as best as possible.
F_all <- (F_shore_rec$F[75:125] + F_boat_rec$F[75:125] + F_commercial$F[75:125])
plot(x = F_ss$x, F_ss$y, lwd = 2, col = "steelblue", type = "l")
lines(x = 1975:2025, y = F_all, col = "firebrick", lwd = 2)

# make sure the whole time period makes sense
ggplot() + # this should somewhat match Fisher et al. 2025, figure 3.8d, shore rec excluded
  geom_line(data = F_commercial, aes(x = year, y = F), linewidth = 1.5, linetype = "dashed") +
  geom_line(data = F_boat_rec, aes(x = year, y = F), linewidth = 1.5, linetype = "dotted") +
  geom_line(data = F_shore_rec, aes(x = year, y = F), linewidth = 1.5, linetype = "solid") +
  theme_minimal()



## save outputs to use in 04_A-C ----------------------------------------------

saveRDS(Q_commercial, "data/output_data/04_D_commercial_q.rds")
saveRDS(Q_boat_rec, "data/output_data/04_D_boat_rec_q.rds")
saveRDS(Q_shore_rec, "data/output_data/04_D_shore_rec_q.rds")
