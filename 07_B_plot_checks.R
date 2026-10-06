# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Previously simulated dataframes
# Task:    Plot-check the outputs of the 07_A reconstruction script
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    October 2026

# -----------------------------------------------------------------------------

# Status: 

# -----------------------------------------------------------------------------

## Load libraries -------------------------------------------------------------

library(tidyverse) # for data manipulation


### CHECK: population trends over time ----------------------------------------

total_pop <- sapply(RECONS_pop, function(year_list) {
  sum(sapply(year_list, sum))
})

p_pop <- ggplot() +
  geom_line(aes(x = 1900:(1900+max_year-1), y = total_pop), linewidth = 1) +
  theme_minimal()


### CHECK: age and length structure in the last time step ---------------------

total_by_length_age <- apply(current_pop, c(2,3), sum)
age_structure <- colSums(total_by_length_age)
size_structure <- rowSums(total_by_length_age)

ggplot() +
  geom_line(aes(x = 1:length(age_structure), 
                y = age_structure), 
            linewidth = 1) +
  labs(x = "Age", y = "Total abundance", title = "age structure") +
  theme_minimal()

ggplot() +
  geom_line(aes(x = 1:length(size_structure), 
                y = size_structure), 
            linewidth = 1) +
  labs(x = "Size", y = "Total abundance", title = "size structure") +
  theme_minimal()


### CHECK: spawning biomass ---------------------------------------------------

SSB0 <- readRDS("data/output_data/06_burn_in_SSB0.rds")
SSB_timeseries <- RECONS_SSB
SSB_ss <- read.csv("data/input_data/digitised_plots_for_checking/SSB_digitised_from_stock_assessment.csv") |> 
  arrange(year)

B_rel <- SSB_timeseries$SSB / SSB0
par(mfrow = c(1, 1))

# below, the red line should align on the blue as best as possible.
p_ssb <- ggplot() +
  geom_line(aes(x = SSB_ss$year, y = SSB_ss$SSB), 
            linewidth = 1, 
            col = "steelblue") +
  ylim(c(0, 1)) +
  theme_minimal() +
  geom_line(aes(x = 1975:2025, y = B_rel[74:124]), 
            linewidth = 1, 
            col = "firebrick") +
  labs(x = "", y = "SSB relative to 1900", title = "SSB relative to 1900")
# here we are checking whether the relative SSB (this year's SSB/SSB0) matches the stock assessment.


### fish density --------------------------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds")
ages_to_plot <- c(1, 2)
months_to_plot <- 1:12
year_to_plot <- 120

age_month_df <- lapply(months_to_plot, function(m) {
  lapply(ages_to_plot, function(a) {
    cell_totals <- apply(RECONS_pop[[year_to_plot]][[m]][,, a], 1, sum)
    water$fish <- cell_totals
    water$age <- factor(paste("Age", a), levels = paste("Age", ages_to_plot))
    water$month <- factor(month.abb[m], levels = month.abb)
    water
  }) %>% bind_rows()
}) %>% bind_rows()

ggplot(age_month_df) +
  geom_sf(aes(fill = fish), color = NA) +
  scale_fill_viridis_c(name = "Fish") +
  facet_grid(age ~ month) +
  theme_void() +
  ggtitle(paste0("Fish distribution by age and month (year ", year_to_plot, ")"))


### relative catch of fleets --------------------------------------------------

yearly_catch <- as.data.frame(RECONS_catch_weight)
names(yearly_catch) <- fleet_names
ggplot(yearly_catch[75:124,]) +
  geom_line(lwd = 1, aes(x = 75:124, y = commercial), colour = "red") +
  geom_line(lwd = 1, aes(x = 75:124, y = boat_rec), colour = "blue") +
  geom_line(lwd = 1, aes(x = 75:124, y = shore_rec), colour = "green") +
  labs(x = "Year", y = "Catch (kg)", title = "Annual catch by fleet") +
  theme_minimal()


### CHECK: wadandi metric -----------------------------------------------------

test <- t(RECONS_catch_june) # shore fishing catches

for (YEAR in 1:n_yrs_modelled) {
  
  plot_dat <- data.frame(
    length = (1:length(test[, YEAR]))*50,
    catch_june = test[, YEAR]) |> 
    mutate(perc = catch_june/sum(catch_june))
  
  p <- ggplot(data = plot_dat) +
    geom_line(aes(x = length, y = catch_june)) +
    ylim(c(0, max(test))) +
    #ylim(c(0, 1)) +
    labs(x = "fish length (mm)", y = "catch number", title = 1900+YEAR-1) + 
    theme_minimal()
  ggsave(p, filename = paste0("plots/script_plot_checks/07/07_catch_lengths_boat_rec/08_length_freq_", YEAR, ".jpg"))
}


### CHECK: effort distribution ------------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds")
year_to_plot = 110
water_effort <- water %>%
  mutate(
    #commercial = RECONS_effort[[year_to_plot]][, 12, 1],
    boat_rec   = RECONS_effort[[year_to_plot]][, 12, 2]
    #shore_rec  = RECONS_effort[[year_to_plot]][, 12, 3]
  ) %>%
  pivot_longer(
    cols = c(
      #commercial, 
      boat_rec 
      #shore_rec
    ),
    names_to = "fleet",
    values_to = "effort"
  )

ggplot(water_effort) +
  geom_sf(aes(fill = effort), color = NA) +
  scale_fill_viridis_c(name = "Effort") +
  facet_wrap(~ fleet, nrow = 1) +
  theme_void() +
  ggtitle(paste0("Fishing effort by fleet (year ", year_to_plot, " Month 12)"))


plot_data <- water_effort[water_effort$fleet == "commercial", c("effort")] # keeps geometry automatically

mapview::mapview(plot_data,
                 zcol = "effort",
                 #col.regions = viridisLite::viridis(100),
                 layer.name = "Effort")


### END ###

