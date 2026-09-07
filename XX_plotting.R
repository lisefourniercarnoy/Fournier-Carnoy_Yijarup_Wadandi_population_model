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


## Load files -----------------------------------------------------------------

library(sf)
library(ggplot2)
library(dplyr)
colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))

common_crs <- 4326

bbox <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -31),
                crs = common_crs) %>%
  st_as_sfc()

# read + fix + reproject water grid
water <- readRDS("data/output_data/02_watergrid.rds") %>%
  st_make_valid() %>%
  st_transform(common_crs)

# read + fix + reproject land layer
wa_map <- st_read("data/output_data/01_B_land.shp") %>%
  st_make_valid() %>%
  st_transform(common_crs)

# confirm they now match
st_crs(water) == st_crs(wa_map)  # should be TRUE

# crop both to your bbox (st_crop is much faster than st_intersection
# for a simple rectangular clip — it doesn't compute new boundary
# geometry along cut edges, it just discards what's outside)
water_crop  <- st_crop(water, bbox)
wa_map_crop <- st_crop(wa_map, bbox)

ggplot() +
  geom_sf(data = wa_map_crop, 
          aes(fill = "Land"), color = NA) +

  # zones
  geom_sf(data = water_crop[!is.na(water_crop$zone), ],
          aes(fill = "No-Take Zone"), color = NA, alpha = 0.75) +
  
  # temporal closure
  geom_sf(data = water_crop[water_crop$TC_status == TRUE, ],
          aes(fill = "Temporal closure"), col = NA, alpha = 0.75) +

  # base grid
  geom_sf(data = water_crop,
          aes(color = "Grid cells"), fill = NA) +
  
  # spawning
  geom_sf(data = water_crop[water_crop$spawning_status == TRUE, ],
          aes(col = "Spawning cells"), fill = NA, alpha = 0.75) +
  
  # legend
  scale_fill_manual(
    name = NULL,
    values = c("No-Take Zone" = colour_palette[5], 
               "Land" = colour_palette[3],
               "Temporal closure" = colour_palette[6]
               )
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Grid cells" = colour_palette[5],
               "Spawning cells" = colour_palette[4]
               )
  ) +
  
  # themes
  theme_minimal() +
  theme(axis.text = element_blank()) +
  labs(x = NULL, y = NULL)


"#08415C" "#3E6990" "#F9EBE0" "#F18805" "#A3320B" "#6B0504"



