# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Outputs from multiple scripts
# Task:    Produce tables for the manuscript and appendices
# Author:  Lise Fournier-Carnoy
# Date:    June 2026

# -----------------------------------------------------------------------------

# Status:  

# -----------------------------------------------------------------------------

library(tidyverse) # for data manipulation
library(ggplot2) # for plotting
library(gridExtra) # for plot arranging
library(sf) # for dealing with shapefiles
library(sfnetworks) # to create geospatial networks

rm(list = ls())
par(mfrow = c(1, 1))

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))
common_crs = 7850 

## Table A2_1: habitat affinity model output ----------------------------------

size_split <- 375

hab_aff <- readRDS( "data/output_data/03_A_habitat_affinity_outputs/03_A_model.rds")
hab_aff

library(flextable)
library(officer)  # <-- needed for fp_border



# Pull full summary table from the model
coef_table <- summary(hab_aff)$coefficients

# Build data frame with all relevant columns
model_table <- data.frame(
  Parameter = c(
    "Intercept",
    "Depth",
    "Depth²",
    paste0("Size class (<", length_split, "mm)"),    "P(reef)",
    "P(sand)",
    paste0("Depth × Size class (<", length_split, "mm)"),
    paste0("Depth² × Size class (<", length_split, "mm)"),
    paste0("P(reef) × Size class (<", length_split, "mm)"),
    paste0("P(sand) × Size class (<", length_split, "mm)")
  ),
  Estimate  = round(coef_table[, "Estimate"],2),
  SE        = round(coef_table[, "Std. Error"],2),
  Z         = round(coef_table[, "z value"],2),
  P         = round(coef_table[, "Pr(>|z|)"],2)
)

# Helper to format p-values
fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0.001", formatC(p, digits = 3, format = "f"))
}

model_table$P_fmt <- fmt_p(model_table$P)

# Build flextable
ft <- flextable(model_table[, c("Parameter", "Estimate", "SE", "Z", "P_fmt")]) |>
  set_header_labels(
    Parameter = "Parameter",
    Estimate  = "Estimate",
    SE        = "SE",
    Z         = "z",
    P_fmt     = "p-value"
  ) |>
  add_header_lines("Table 1. Negative binomial GLM coefficient estimates") |>
  hline(i = 6, border = fp_border(width = 0.5)) |>   # split main effects / interactions
  bold(part = "header") |>
  italic(i = 1, part = "header") |>                  # italicise sub-header row (optional)
  align(j = c("Estimate", "SE", "Z", "P_fmt"), align = "right",  part = "all") |>
  align(j = "Parameter",                             align = "left",  part = "all") |>
  colformat_num(j = c("Estimate", "SE", "Z"), digits = 4) |>
  font(fontname = "Times New Roman", part = "all") |>
  fontsize(size = 11, part = "all") |>
  autofit() |>
  theme_booktabs()

ft

save_as_docx(ft, path = "plots/XX_table_A2_1.docx")



