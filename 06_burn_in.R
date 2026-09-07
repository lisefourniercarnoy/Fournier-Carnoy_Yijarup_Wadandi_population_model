# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    ?
# Task:    Set up a population out of the population parameters and fishing effort
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    August 2025

# -----------------------------------------------------------------------------

# Status:

# -----------------------------------------------------------------------------

rm(list = ls())

# Load libraries
library(tidyverse) # for data manipulation
library(sf) # for spatial objects
#library(MQMF)
library(Rcpp) # for reading the C++ functions in R
library(RcppArmadillo) # for reading the C++ functions in R
library(abind) # for dealing with arrays i think.

## 1. Set up ------------------------------------------------------------------

### 1.1 Read in the functions -------------------------------------------------

sourceCpp("functions/age-length_functions/run_full_model_function.cpp", verbose = TRUE)
colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))


### 1.2 Get model parameters --------------------------------------------------

n_yrs_modelled <- 60 # number of burn-in years - at least a fish's full life.

max_cell    <- nrow(readRDS("data/output_data/02_watergrid.rds")) # Number of cells in the model
max_age     <- 40 # The max age of the fish in the model, see script 05 for correct value
n_lengths   <- length(readRDS("data/output_data/05_length_bins.rds"))
max_year    <- n_yrs_modelled # Number of years the model should run for 

starting_pop  <- readRDS("data/output_data/05_starting_population.rds") %>% glimpse()
weight        <- readRDS("data/output_data/05_weight.rds") %>% glimpse()
selectivity   <- readRDS("data/output_data/05_selectivity_retention.rds") %>% glimpse()# FOR NOw THEY ARE THE SAME FOR ALL FLEETS BUT SHOULD END UP BEING DIFFERENT AT SOME POINT
nat_mort      = 0.12 # from table 4.2, p12 https://library.dpird.wa.gov.au/cgi/viewcontent.cgi?article=1240&context=fr_rr

spawn_months <- c(10, 11, 12) # 1-indexed. the function deals with zero-indexing it.
BHa           = as.double(readRDS("data/output_data/05_Beverton-Holt_alpha.rds")) # see script 05
BHb           = as.double(readRDS("data/output_data/05_Beverton-Holt_beta.rds")) # see script 05
PF            = 0.5 # proportion expected to be females
hyperallo     = (1.26 + 1.14 + 1.33)/3 # average from 3 Sparids in Barneche 2018 (1.26, 1.14 and 1.33)
mature         <- readRDS("data/output_data/05_maturity.rds") %>% glimpse()
settlement     <- readRDS("data/output_data/03_B_recruitment.rds") %>% glimpse(); settlement <- settlement[, 1] # selecting a single column because the function expects a vector
age_transition <- readRDS("data/output_data/05_age_transition_matrix.rds") %>% glimpse()
adult_movement <- readRDS("data/output_data/03_b_adult_movement_15_swim_speed.rds") %>% glimpse()
juv_movement   <- readRDS("data/output_data/03_b_juv_movement_10_swim_speed.rds") %>% glimpse()
l_spawn_movement <- readRDS("data/output_data/03_b_large_spawning_movement_15_swim_speed.rds") %>% glimpse()
s_spawn_movement <- readRDS("data/output_data/03_b_small_spawning_movement_15_swim_speed.rds") %>% glimpse()


fleet_names <- c(
  "commercial", 
  "boat_rec", 
  "shore_rec")

com_info <- readRDS("data/output_data/04_A_commercial_fishing_info.rds") %>% glimpse()
com_info$fishing_days[com_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

brec_info <- readRDS("data/output_data/04_C_boat_rec_fishing_info.rds") %>% glimpse()
brec_info$fishing_days[brec_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.

srec_info <- readRDS("data/output_data/04_B_shore_rec_fishing_info.rds") %>% glimpse()
srec_info$fishing_days[srec_info$fishing_days == 0] <- 1e-10 # replace zero with small number to avoid calculations freaking out.


# for the burn-in, we'll use a constant low level of fishing, so replace all years with 1900 fishing effort, catchability, etc.
com_info$fishing_days[,,1:dim(com_info$fishing_days)[[3]]] <- com_info$fishing_days[,,1]
brec_info$fishing_days[,,1:dim(brec_info$fishing_days)[[3]]] <- brec_info$fishing_days[,,1]
srec_info$fishing_days[,,1:dim(srec_info$fishing_days)[[3]]] <- srec_info$fishing_days[,,1]

com_info$catchability[,,1:dim(com_info$catchability)[[3]]] <- com_info$catchability[,,1]
brec_info$catchability[,,1:dim(brec_info$catchability)[[3]]] <- brec_info$catchability[,,1]
srec_info$catchability[,,1:dim(srec_info$catchability)[[3]]] <- srec_info$catchability[,,1]

first_vals <- com_info$attractivity[[1]][, , 1]  # 1585 x 16 matrix (month 1, year 1)
com_info$attractivity <- lapply(com_info$attractivity, function(x) {
  array(first_vals, dim = dim(x))
})

first_vals <- brec_info$attractivity[[1]][, , 1]  # 1585 x 16 matrix (month 1, year 1)
brec_info$attractivity <- lapply(brec_info$attractivity, function(x) {
  array(first_vals, dim = dim(x))
})

first_vals <- srec_info$attractivity[[1]][, , 1]  # 1585 x 16 matrix (month 1, year 1)
srec_info$attractivity <- lapply(srec_info$attractivity, function(x) {
  array(first_vals, dim = dim(x))
})

fleet_info = list(com_info, 
                  brec_info,
                  srec_info
                  )


### 1.3 Set up the initial population -----------------------------------------

BURN_IN_pop <- list() # keep record of the burn-in outputs
months_to_save <- c(Jan = 1, Feb = 2, Mar = 3, # which months to  sve for plot checks
                    Apr = 4, May = 5, Jun = 6, 
                    Jul = 7, Aug = 8, Sep = 9, 
                    Oct = 10, Nov = 11, Dec = 12)

BURN_IN_effort <- list() # cell x month x fleet, one array per year
BURN_IN_F <- list() # one 12 x n_fleets matrix per year
BURN_IN_catch_weight <- array(0, dim = c(n_yrs_modelled, length(fleet_names))); colnames(BURN_IN_catch_weight) <- fleet_names
BURN_IN_SSB <- list()

current_pop <- array(0, dim = c(max_cell, n_lengths, max_age)) # for every cell (row), and every length (column) across all fish ages (matrix slice), we will have a population

# distribute initial population
for(AGE in 1:max_age){
  total_this_age <- starting_pop[AGE, ]
  # distribute proportionally to settlement — same habitat weighting as recruits
  settlement_prop <- settlement / sum(settlement) # length max_cell
  
  # outer product: each cell's proportion × each length-class abundance
  current_pop[, , AGE] <- outer(settlement_prop, total_this_age)
}

cat("Total fish initialised:", sum(current_pop), "\n")


## 2. Start the burn-in -------------------------------------------------------

Start = Sys.time()
for (YEAR in 0:(max_year-1)){ # max_year-1 because it starts at 0.
  
  # loop over all the Rcpp functions in the model
  ModelOutput <- run_full_model_function(YEAR = YEAR,
                                         max_cell = max_cell,
                                         max_age = max_age,
                                         max_year = max_year,
                                         n_lengths = n_lengths,
                                         current_pop = current_pop,
                                         weight = weight,
                                         selectivity = selectivity,
                                         age_transition = age_transition,
                                         natural_mortality = nat_mort,
                                         spawning_months = spawn_months,
                                         BHa = BHa,
                                         BHb = BHb,
                                         PF = PF,
                                         ha_scaling = hyperallo,
                                         maturity = mature,
                                         settlement = settlement,
                                         adult_movement_prob = adult_movement,
                                         juv_movement_prob = juv_movement,
                                         small_spawn_movement_prob = s_spawn_movement,
                                         large_spawn_movement_prob = l_spawn_movement,
                                         fleet_names = fleet_names,
                                         fleet_info = fleet_info
  )

  # save the population at the end of the year (to give to the loop again as January population)
  current_pop <- ModelOutput$next_pop #  cell x length x age
  
  
  # fill objects to check after the burn-in
  BURN_IN_pop[[YEAR+1]]         <- setNames(ModelOutput$master_current_pop[months_to_save], names(months_to_save))
  BURN_IN_effort[[YEAR+1]]      <- ModelOutput$effort_by_fleet # array: cell x month x fleet
  BURN_IN_F[[YEAR+1]]           <- ModelOutput$fishing_mortality
  BURN_IN_catch_weight[YEAR+1,] <- ModelOutput$yearly_catch
  BURN_IN_SSB[[YEAR+1]]         <- ModelOutput$spawning_biomass  # cell, summed across spawning months
  
  
  
  
  
  # # add a plot check, to avoid wasting time on running the function
  # total_by_length_age <- apply(ModelOutput$master_current_pop[[12]], c(2, 3), sum)
  # age_structure       <- colSums(total_by_length_age)   # summed over lengths
  # size_structure       <- rowSums(total_by_length_age)   # summed over ages
  # 
  # png(
  #   filename = file.path("plots/script_plot_checks/06/06_age_length_structures/", sprintf("BURN_IN_year_%04d.png", YEAR + 1)),
  #   width = 1200, height = 500, res = 120
  # )
  # 
  # par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  # 
  # plot(age_structure,
  #      type = "l", lwd = 2, col = colour_palette[6],
  #      xlab = "Age class", ylab = "Total abundance",
  #      main = paste0("Age structure — Year ", YEAR + 1))
  # 
  # plot(size_structure,
  #      type = "l", lwd = 2, col = colour_palette[4],
  #      xlab = "Length bin", ylab = "Total abundance",
  #      main = paste0("Size structure — Year ", YEAR + 1))
  # 
  # dev.off()
  
}

End = Sys.time(); Runtime = End - Start; Runtime
# after the burn-in, make sure you run the sections below, where some outputs are saved for the historical reconstruction.


## 3. Plot checks -------------------------------------------------------------
### CHECK: burn-in population stability over time -----------------------------

total_pop <- sapply(BURN_IN_pop, function(year_list) {
  sum(sapply(year_list, sum))
})

plot(1:max_year, total_pop, type = "l",
     xlab = "Year", ylab = "Total abundance",
     main = "Total Population during burn-in")

saveRDS(BURN_IN_pop[[60]], "data/output_data/06_burn_in_population.rds")


### CHECK: spawning biomass ---------------------------------------------------

SSB_timeseries <- sapply(BURN_IN_SSB, sum)  # sum over cells, one value per year

par(mfrow = c(1, 1))

# below, the red line should align on the blue as best as possible.
plot(x = 1:length(SSB_timeseries), y = SSB_timeseries, lwd = 2, col = "steelblue", type = "l",
     xlab = "Year", ylab = "SSB",
     main = "Effective reproductive output")

# we'll also export the final year's SSB to compare the historical reconstruction period's SSB to:
SSB0 <- SSB_timeseries[length(SSB_timeseries)]
saveRDS(SSB0, file = "data/output_data/06_burn_in_SSB0.rds")


### CHECK: fish density -------------------------------------------------------

library(purrr)
library(sf)
library(patchwork)

water <- readRDS("data/output_data/02_watergrid.rds")

lengths_to_plot <- c(5, 10, 15)
months_to_plot  <- 1:12
years_to_plot   <- 54:59

output_dir <- "plots/gif_frames/06_burn_in_movement"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Step 1: compute per-length global min/max (fixed across all frames) ---
length_limits <- lapply(lengths_to_plot, function(l) {
  totals <- unlist(lapply(years_to_plot, function(y) {
    lapply(months_to_plot, function(m) {
      apply(BURN_IN_pop[[y]][[m]][, l, ], 1, sum)
    })
  }))
  c(min = min(totals), max = max(totals))
})
names(length_limits) <- paste("Length", lengths_to_plot * 50)

# --- Step 2: build one frame as a set of independently-scaled panels ---
make_frame <- function(y, m) {
  
  panels <- lapply(lengths_to_plot, function(l) {
    cell_totals <- apply(BURN_IN_pop[[y]][[m]][, l, ], 1, sum)
    water$fish  <- cell_totals
    label <- paste("Length", l * 50)
    lims  <- length_limits[[label]]
    
    ggplot(water) +
      geom_sf(aes(fill = fish), color = NA) +
      scale_fill_gradientn(
        colors = colour_palette,
        limits = lims,
        name = "Fish",
        guide = guide_colorbar(direction = "horizontal", title.position = "top")
      ) +
      theme_void() +
      theme(
        legend.position = "none"
      ) +
      ggtitle(label)
  })

  combined <- wrap_plots(panels, nrow = 1) +
    plot_annotation(title = paste0("Year ", y, " - ", month.abb[m]))
  
  fname <- file.path(output_dir,
                     sprintf("frame_y%02d_m%02d.png", y, m))
  ggsave(fname, combined,
         width = 5, height = 6, dpi = 300)
}

# --- Step 3: loop over every year x month combination ---
walk(years_to_plot, function(y) {
  walk(months_to_plot, function(m) make_frame(y, m))
})


### CHECK: effort distribution ------------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds")

year_to_plot <- n_yrs_modelled

water_effort <- water %>%
  mutate(
    commercial = BURN_IN_effort[[year_to_plot]][, 12, 1],
    boat_rec   = BURN_IN_effort[[year_to_plot]][, 12, 2],
    shore_rec  = BURN_IN_effort[[year_to_plot]][, 12, 3]
  ) %>%
  pivot_longer(
    cols = c(commercial, boat_rec, shore_rec),
    names_to = "fleet",
    values_to = "effort"
  )

ggplot(water_effort) +
  geom_sf(aes(fill = (effort)), color = NA) +
  scale_fill_viridis_c(name = "Effort") +
  facet_wrap(~ fleet, nrow = 1) +
  theme_void() +
  ggtitle(paste0("Fishing effort by fleet (Year ", year_to_plot, ", Month 12)"))


### CHECK: relative catch of fleets -------------------------------------------

yearly_catch <- as.data.frame(BURN_IN_catch_weight)
names(yearly_catch) <- fleet_names
ggplot(yearly_catch) +
  geom_line(lwd = 1, aes(x = 1:n_yrs_modelled, y = commercial), colour = "red") +
  geom_line(lwd = 1, aes(x = 1:n_yrs_modelled, y = boat_rec), colour = "blue") +
  geom_line(lwd = 1, aes(x = 1:n_yrs_modelled, y = shore_rec), colour = "green") +
  labs(x = "Year", y = "Catch (kg)", title = "Annual catch by fleet") +
  theme_bw()

## END ##


### CHECK: mega-plot ----------------------------------------------------------

#### SECTION 1: effort distribution -------------------------------------------

library(patchwork) # for combining plots with independent legends — verify installed

water <- readRDS("data/output_data/02_watergrid.rds")
bbox <- st_bbox(water)  # compute once outside the loop, since geometry doesn't change

wa_map <- st_read("data/input_data/Q_aus_land_high_res_no_estuary.shp") |> st_as_sf()
out_dir <- "plots/script_plot_checks/06/06_effort_distribution_plots/"

# compute global max effort PER FLEET, across all years, months, and cells
fleet_idx <- 1:3 # matches order in BURN_IN_effort's 3rd dimension: commercial, boat_rec, shore_rec
global_max_effort <- sapply(fleet_idx, function(f) {
  max(sapply(1:n_yrs_modelled, function(y) max(BURN_IN_effort[[y]][, , f])))
})
names(global_max_effort) <- fleet_names

for (YEAR in 1:n_yrs_modelled) {
  for (MONTH in 1:12) {
    
    water_year_month <- water %>%
      mutate(
        commercial = BURN_IN_effort[[YEAR]][, MONTH, 1],
        boat_rec   = BURN_IN_effort[[YEAR]][, MONTH, 2],
        shore_rec  = BURN_IN_effort[[YEAR]][, MONTH, 3]
      )
    
    # build one plot per fleet, each with its own fixed scale
    plot_list <- lapply(fleet_names, function(f) {
      ggplot() +
        geom_sf(data = wa_map, fill = "grey90", color = "grey40") +
        geom_sf(data = water_year_month, aes(fill = .data[[f]]), color = NA) +
        coord_sf(
          xlim = c(bbox["xmin"], bbox["xmax"]),
          ylim = c(bbox["ymin"], bbox["ymax"]),
          expand = FALSE
        ) +
        scale_fill_gradient(
          name = "Effort", low = "white", high = colour_palette[6],
          limits = c(0, global_max_effort[[f]])
        ) +
        theme_minimal() +
        ggtitle(f)
    })

    p <- wrap_plots(plot_list, nrow = 1) +
      plot_annotation(title = paste0("Fishing effort by fleet — Year ", YEAR, ", Month ", MONTH))
    
    ggsave(
      filename = file.path(
        out_dir,
        sprintf("effort_year%02d_month%02d.png", YEAR, MONTH)
      ),
      plot = p, width = 10, height = 7, dpi = 300
    )
    
  }
}


#### SECTION 2: length distributions ------------------------------------------

n_years_done <- length(BURN_IN_pop_monthly)

# precompute December age/size structure for every year first
age_structure_all  <- matrix(NA_real_, nrow = max_age, ncol = n_years_done)
size_structure_all <- matrix(NA_real_, nrow = n_lengths, ncol = n_years_done)

for (YEAR in 1:n_years_done) {
  total_by_length_age <- apply(BURN_IN_pop_monthly[[YEAR]][[12]], c(2, 3), sum)
  age_structure_all[, YEAR]  <- colSums(total_by_length_age)  # summed over lengths
  size_structure_all[, YEAR] <- rowSums(total_by_length_age)  # summed over ages
}

# fixed axis limits across all years, so scale doesn't jump between plots
age_ylim  <- c(0, max(age_structure_all, na.rm = TRUE))
size_ylim <- c(0, max(size_structure_all, na.rm = TRUE))

for (YEAR in 1:n_years_done) {
  
  png(
    filename = file.path(
      "plots/script_plot_checks/06/06_age_length_structures/",
      sprintf("BURN_IN_year_%04d.png", YEAR)
    ),
    width = 1200, height = 500, res = 120
  )
  
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  # --- age structure panel ---
  plot(age_structure_all[, YEAR],
       type = "n", # empty frame first, so we can layer lines under it
       ylim = age_ylim,
       xlab = "Age class", ylab = "Total abundance",
       main = paste0("Age structure — Year ", YEAR))
  
  if (YEAR > 1) {
    for (past_yr in 1:(YEAR - 1)) {
      # older years get progressively lighter/more transparent
      fade <- 0.15 + 0.55 * (past_yr / (YEAR - 1))  # ranges ~0.15 (oldest) to ~0.7 (most recent past)
      lines(age_structure_all[, past_yr],
            lwd = 1.5,
            col = adjustcolor(colour_palette[6], alpha.f = fade))
    }
  }
  lines(age_structure_all[, YEAR], lwd = 2, col = colour_palette[6])
  
  # --- size structure panel ---
  plot(size_structure_all[, YEAR],
       type = "n",
       ylim = size_ylim,
       xlab = "Length bin", ylab = "Total abundance",
       main = paste0("Size structure — Year ", YEAR))
  
  if (YEAR > 1) {
    for (past_yr in 1:(YEAR - 1)) {
      fade <- 0.15 + 0.55 * (past_yr / (YEAR - 1))
      lines(size_structure_all[, past_yr],
            lwd = 1.5,
            col = adjustcolor(colour_palette[4], alpha.f = fade))
    }
  }
  lines(size_structure_all[, YEAR], lwd = 2, col = colour_palette[4])
  
  dev.off()
}


#### SECTION 3: Annual catch --------------------------------------------------

out_dir_catch <- "plots/script_plot_checks/06/06_catch/"
dir.create(out_dir_catch, recursive = TRUE, showWarnings = FALSE)

yearly_catch <- as.data.frame(BURN_IN_catch_weight)
names(yearly_catch) <- fleet_names
yearly_catch$year <- 1:n_yrs_modelled

catch_df <- yearly_catch %>%
  pivot_longer(cols = all_of(fleet_names), names_to = "fleet", values_to = "catch")

global_max_catch <- max(catch_df$catch, na.rm = TRUE)
global_max_year  <- max(catch_df$year, na.rm = TRUE)

frame_i <- 0

for (YEAR in 1:n_yrs_modelled) {

    frame_i <- frame_i + 1
    frame_df <- catch_df %>% filter(year <= YEAR) # only full years are "known"
    
    p <- ggplot(frame_df, aes(x = year, y = catch, colour = fleet)) +
      geom_line(lwd = 1) +
      scale_colour_manual(values = c(commercial = "red", boat_rec = "blue", shore_rec = "green")) +
      coord_cartesian(xlim = c(1, global_max_year), ylim = c(0, global_max_catch)) +
      labs(x = "Year", y = "Catch (kg)",
           title = paste0("Catch by fleet — Year ", YEAR)) +
      theme_bw()
    
    ggsave(
      filename = file.path(out_dir_catch, sprintf("catch_frame_%04d.png", frame_i)),
      plot = p, width = 10, height = 5, dpi = 120
    )
}


#### SECTION 4: fish distribution ---------------------------------------------

water <- readRDS("data/output_data/02_watergrid.rds")
bbox  <- st_bbox(water)

wa_map <- st_read("data/input_data/Q_aus_land_high_res_no_estuary.shp") |> st_as_sf()
out_dir_length <- "plots/script_plot_checks/06/06_movement/"
dir.create(out_dir_length, recursive = TRUE, showWarnings = FALSE)

length_bins <- readRDS("data/output_data/05_length_bins.rds") # numeric vector, length n_lengths

# group the 50mm bins into 200mm groups (4 bins per group), based on actual bin values
group_breaks <- seq(floor(min(length_bins) / 200) * 200,
                    ceiling(max(length_bins) / 200) * 200,
                    by = 200)

length_group <- cut(length_bins,
                    breaks = group_breaks,
                    include.lowest = TRUE, right = FALSE)

group_levels <- levels(length_group)

# compute global max abundance per 200mm group, across ALL years AND months, for a fixed colour scale
cell_totals_all <- lapply(1:n_yrs_modelled, function(YEAR) {
  lapply(1:12, function(MONTH) {
    pop_this_month <- BURN_IN_pop_monthly[[YEAR]][[MONTH]] # cell x length x age
    sapply(group_levels, function(g) {
      length_idx <- which(length_group == g)
      apply(pop_this_month[, length_idx, , drop = FALSE], 1, sum)
    })
  })
}) %>% unlist(recursive = FALSE)

global_max_by_group <- apply(
  do.call(rbind, cell_totals_all), 2, max, na.rm = TRUE
)
names(global_max_by_group) <- group_levels
global_max_all <- max(global_max_by_group) # single scale across ALL groups, years, months

for (YEAR in 1:n_yrs_modelled) {
  for (MONTH in 1:12) {
    
    pop_this_month <- BURN_IN_pop_monthly[[YEAR]][[MONTH]] # cell x length x age
    
    length_df <- lapply(group_levels, function(g) {
      length_idx <- which(length_group == g)
      cell_totals <- apply(pop_this_month[, length_idx, , drop = FALSE], 1, sum)
      water_g <- water
      water_g$fish <- cell_totals
      water_g$length_group <- factor(g, levels = group_levels)
      water_g
    }) %>%
      bind_rows()
    
    p <- ggplot() +
      geom_sf(data = wa_map, fill = "grey90", color = "grey40") +
      geom_sf(data = length_df, aes(fill = fish), color = NA) +
      coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]),
        expand = FALSE
      ) +
      scale_fill_gradient(
        name = "Fish", low = "white", high = colour_palette[5],
        limits = c(0, global_max_all)
      ) +
      facet_wrap(~ length_group, nrow = 1) +
      theme_minimal() +
      ggtitle(paste0("Fish distribution by length group (200mm) — Year ", YEAR, ", Month ", MONTH))
    
    ggsave(
      filename = file.path(
        out_dir_length,
        sprintf("length_dist_year%02d_month%02d.png", YEAR, MONTH)
      ),
      plot = p, width = 14, height = 8, dpi = 500
    )
  }
}
