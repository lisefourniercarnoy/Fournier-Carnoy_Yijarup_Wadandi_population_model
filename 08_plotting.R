# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    ?
# Task:    Set up a population out of the population parameters and fishing effort
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    May 2025

# -----------------------------------------------------------------------------

# Status: Starting out...

# -----------------------------------------------------------------------------

rm(list = ls()) # clear environment

## Load libraries -------------------------------------------------------------

library(tidyverse)
library(sf)
# library(sfnetworks)
# library(raster)
# library(stringr)
# library(forcats)
# library(RColorBrewer)
# library(geosphere)
# library(forcats)
# library(ggridges)
# library(grid)
# library(gridExtra)
# library(gtable)
# library(purrr)
# library(matrixStats)
# library(rcartocolor)
# library(ggtext)
# library(abind)
# library(scales)

## Load functions -------------------------------------------------------------
source("XX_functions.R")

## Load files -----------------------------------------------------------------

bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -33.2), crs = 4326) %>% st_as_sfc()

maturity          <- readRDS("data/output_data/05_maturity.rds") %>% glimpse()
weight            <- readRDS("data/output_data/05_weight.rds") %>% glimpse()

water             <- readRDS("data/output_data/02_watergrid.rds") %>% st_make_valid() %>% glimpse(); plot(water)
no_take_list      <- readRDS("data/output_data/02_no_take_list.rds") %>% glimpse()
ntz               <- st_read("data/output_data/01_wadandi_NTZ.shp") %>% glimpse(); plot(ntz)
BR                <- st_read("data/input_data/wadandi_boat_ramps.shp") %>% glimpse(); plot(BR$geometry)
network           <- st_read("data/output_data/03_network_shapefile.shp") %>% glimpse(); plot(network$geometry)
wa_map            <- st_read("data/output_data/01_wadandi_land.shp") %>% glimpse(); plot(wa_map)

tot_pop_list <- list(readRDS("simulations/dummy_run/age_distribution_try1.rds"), readRDS("simulations/dummy_run/age_distribution_try1.rds"))
total_pop <- total.pop.format(pop.file.list = tot_pop_list, 
                              scenario.names = c("first try", "dummy"), 
                              nsim = 2, 
                              nyears = 44, 
                              startyear = 15, 
                              maxage = 30, 
                              mat = maturity, 
                              kg = Weight)

