# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    MEGlab BOSS and BRUV habitat data (workflow adapted from Claude)
# Task:    Create a habitat prediction
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    December 2025

# -----------------------------------------------------------------------------

# Status:  First try

# -----------------------------------------------------------------------------

rm(list = ls()) # Clear working environment

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))

# Load libraries
library(tidyverse) # for data manipulation
library(FSSgam) # for selecting the models
library(mgcv) # for making the models
library(predicts) # for making the prediction
library(terra) # for extracting bathy values at opcode locations
library(raster) # to deal with rasters
library(sf) # for dealing with polygons
library(stars) # for detrended bathymetry
library(starsExtra) # for detrended bathymetry

# extent of the area we're predicting to
bbox_whole <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -31), crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(4326)
bbox_whole <- extent(st_bbox(bbox_whole))


## Files used in this script --------------------------------------------------

file_obs_hab            <- "data/external_data/sw-network_broad-habitat.csv" # Observed habitat values, cleaned in 01_prep
file_bathy_predictors   <- "data/external_data/sw-network_bathymetry.rds" # made in 01_prep


# 250m habitat predictions ----------------------------------------------------

# -- Let's start with the 250m predictions. 

## Load data ------------------------------------------------------------------

habi <- read.csv(file_obs_hab) %>% 
  mutate(reef = inverts + macroalgae + rock) %>% 
  dplyr::select(-c(inverts, macroalgae, rock, campaignid)) %>% 
  glimpse()
habi$diff <- (habi$reef + habi$seagrass + habi$sand) - habi$total_pts
habi$diff # There should only be zeroes.


# Extract bathymetry values to the location of the opcodes
bathy <- readRDS(file_bathy_predictors) %>% 
  rasterFromXYZ()
crs(bathy) <- "+proj=longlat +datum=WGS84 +no_defs"
bathy <- crop(bathy, bbox_whole)
plot(bathy)

# make derivatives

aspect <- terra::terrain(bathy, "aspect", unit = "degrees", neighbors = 8); plot(aspect)

roughness <- terrain(bathy, "roughness", neighbors = 8); plot(roughness)

zstar <- st_as_stars(bathy)
detrended <- detrend(zstar, parallel = 8) %>%
  terra::rast()
names(detrended) <- c("detrended", "lineartrend")
detrended <- detrended[["detrended"]]
detrended <- raster(detrended); plot(detrended)

# gather with all predictors
predictors <- list(aspect = aspect, depth = bathy, roughness = roughness, detrended = detrended)

for (i in 1:length(predictors)) { # Unifying extents and resampling so that they can be stacked
  reference_raster <- predictors[[1]] # Use the raster with the highest resolution as a reference. Otherwise the re-sampling will apply lower resolutions to everything.
  predictors[[i]] <- terra::resample(predictors[[i]], reference_raster, method = "bilinear")
  cat(paste("New extent of", i, ":", ext(predictors[[i]]), "\n"))
}

list2env(predictors, envir = .GlobalEnv) # Put the elements of the list in the environment

predictors <- raster::stack(depth, aspect, roughness, detrended) # 250m data

names(predictors) <- c("bathy", "aspect", "roughness", "detrended")
plot(predictors)


# extract
coords <- habi[, c("longitude", "latitude")]
bathy_values <- terra::extract(predictors, coords)
bathy_values <- cbind(coords, bathy_values)
bathy_values <- bathy_values[complete.cases(bathy_values), ] # NAs exist for points outside the extent. remove them.

hab <- habi %>%
  inner_join(bathy_values, by = c("longitude", "latitude")) %>% 
  pivot_longer(cols = c("sand", "reef", "seagrass"), # specify the columns to pivot
               values_to = "number",
               names_to = "taxa") %>% 
  glimpse()
hab <- as.data.frame(hab) %>% glimpse()

names(hab)
pred.vars <- c("bathy", "aspect", 
               "roughness", "detrended")

colSums(is.na(hab)) # There should be no NAs. We've removed them just above, so that only LiDAR extent data exists.

# Check for correlation of predictor variables- remove anything highly correlated (>0.95)
round(cor(hab[, pred.vars]), 2) # All good, no high correlations

# Check to make sure Response vector has not more than 80% zeros
unique.vars = unique(as.character(hab$taxa))
head(hab)
unique.vars.use = character()
for(i in 1:length(unique.vars)){
  temp.dat = hab[which(hab$taxa == unique.vars[i]),]
  if(length(which(temp.dat$taxa == 0))/nrow(temp.dat)<0.9){
    unique.vars.use = c(unique.vars.use,unique.vars[i])}
}

unique.vars.use # All remain


## Run the full subset model selection ----------------------------------------
outdir    <- ("data/output_data/01_A_habitat_modelling_outputs/")
resp.vars <- unique.vars.use
out.all   <- list()
var.imp   <- list()


## Loop through the FSS function for each Abiotic taxa ------------------------
hab <- hab %>% na.omit() # removing lidar columns otherwise it wont run
colSums(is.na(hab))
for(i in 1:length(resp.vars)){
  print(resp.vars[i])
  use.dat <- hab[hab$taxa == resp.vars[i],]
  use.dat   <- as.data.frame(use.dat)
  
  # Basic model to compare all the other combinations with
  Model1  <- gam(cbind(number, (total_pts - number)) ~ # Success and failure counts for each habitat
                   s(bathy, bs = 'cr', k = 3),
                 family = binomial(link = "logit"),  data = use.dat)
  
  # Generate variable combinations to test in models
  model.set <- generate.model.set(use.dat = use.dat,
                                  test.fit = Model1,
                                  pred.vars.cont = pred.vars,
                                  cyclic.vars = c("aspect"),
                                  k = 5,
                                  cov.cutoff = 0.7,
                                  max.predictors = 3
  )
  
  # Fit models to all variable combinations
  out.list <- fit.model.set(model.set,
                            max.models = 600,
                            parallel = T)
  names(out.list)
  
  # Examine outcomes
  out.list$failed.models # examine the list of failed models
  mod.table <- out.list$mod.data.out  # look at the model selection table
  mod.table <- mod.table[order(mod.table$AICc), ]
  mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
  out.i     <- mod.table[which(mod.table$delta.AICc <= 2), ] # Top models
  out.all   <- c(out.all, list(out.i))
  var.imp   <- c(var.imp, list(out.list$variable.importance$aic$variable.weights.raw))
  
  # Plot the top models
  for(m in 1:nrow(out.i)){
    best.model.name <- as.character(out.i$modname[m])
    
    png(file = paste0(outdir, resp.vars[i], "_mod_fits.png"))
    if(best.model.name != "null"){
      par(mfrow = c(3, 1), mar = c(9, 4, 3, 1))
      best.model = out.list$success.models[[best.model.name]]
      plot(best.model, all.terms = T, pages = 1, residuals = T, pch = 16)
      mtext(side = 2, text = resp.vars[i], outer = F)}
    dev.off()
  }
}

## Model fits and importance --------------------------------------------------

names(out.all) <- resp.vars
names(var.imp) <- resp.vars
all.mod.fits <- list_rbind(out.all, names_to = "response")
all.var.imp  <- do.call("rbind", var.imp)
out.all


## reload & format the data ---------------------------------------------------

hab_test <- hab %>% 
  pivot_wider(names_from = taxa, values_from = number, values_fill = list(number = 0)) %>% 
  glimpse()

hab_model <- hab_test %>%
  # group_by(sample) %>%
  # summarize(
  #   across(c(#bathy_lidar, roughness_lidar, aspect_lidar, 
  #     bathy, aspect, roughness, detrended), 
  #     ~ first(.x, order_by = .x),  # Keep the first non-NA value
  #     .names = "{.col}"),  # This makes sure column names are kept
  #   reef = sum(reef, na.rm = F),
  #   sand = sum(sand, na.rm = F),
  #   seagrass = sum(seagrass, na.rm = F),
  #   total_pts = first(total_pts, order_by = total_pts),  # Keeps the first non-NA value
  #   .groups = "drop"  # Remove grouping after summarization
  # ) %>% 
  glimpse()

hab_model$diff <- (hab_model$reef + hab_model$seagrass + hab_model$sand) - hab_model$total_pts
hab_model$diff # There should only be zeroes.

hab_model <- hab_model[hab_model$diff == 0, ] # a few points don't add up right, so im just removing them because icbf and also not crucial

preds <- predictors; plot(preds)
preddf <- as.data.frame(preds, xy = TRUE, na.rm = TRUE) %>% glimpse()


## Make models for each habitat and predict extent ----------------------------

fam <- binomial("logit")

make_vars <- function(vars) {
  sapply(vars, function(var) {
    bs_type <- ifelse(grepl("aspect", var), "cc", "cr")
    paste0("s(", var, ", k = 5, bs = '", bs_type, "')")
  }) |> paste(collapse = " + ")
}

# Function to fit the optimal model
make_optimal_model <- function(response_var, predictor_terms, data) {
  formula_str <- paste0(
    "cbind(", response_var, ", total_pts - ", response_var, ") ~ ", predictor_terms
  )
  
  gam(
    as.formula(formula_str),
    data = data,
    method = "REML",
    family = binomial
  )
}

# Container for storing the fitted models
fitted_models <- list()

for (habitat in names(out.all)) {
  
  # Step 1: Extract best model variable names
  best_model_name <- out.all[[habitat]]$modname[1]
  best_model_vars <- strsplit(best_model_name, "\\+")[[1]]
  
  # Step 2: Convert variables into GAM smooth terms
  best_vars <- make_vars(best_model_vars)
  
  # Step 3: Fit the model using cbind(successes, failures)
  fitted_model <- make_optimal_model(
    response_var = habitat,
    predictor_terms = best_vars,
    data = hab_model
  )
  
  # Step 4: Store model
  fitted_models[[habitat]] <- fitted_model
  
  # Optional: print quick summary
  message("Fitted binomial model for habitat: ", habitat)
  print(summary(fitted_model))
}

plot(fitted_models$reef, page = 1, residuals = T, cex = 5)
plot(fitted_models$sand, page = 1, residuals = T, cex = 5)
plot(fitted_models$seagrass, page = 1, residuals = T, cex = 5)




# 
# # Save each model (based on top models above)
# 
# # Reef
# out.all$reef[[1]]
# m_reef_250m <- gam(cbind(reef, total_pts - reef) ~
#                      s(aspect_250m, k = 5, bs = "cc")  +
#                      s(detrended_250m, k = 5, bs = "cr")  +
#                      s(bathy_250m, k = 5, bs = "cr"),
#                      s(roughness_250m, k = 5, bs = "cr"),
#                    data = hab_model, method = "REML", family = binomial("logit"))
# summary(m_reef_250m)
# plot(m_reef_250m, pages = 1, residuals = T, cex = 5)
# 
# # Seagrass
# out.all$seagrass[[1]]
# m_seagrass_250m <- gam(cbind(seagrass, total_pts - seagrass) ~
#                          s(aspect_250m, k = 5, bs = "cc")  +
#                          #s(detrended_250m, k = 5, bs = "cr")  +
#                          s(bathy_250m, k = 5, bs = "cr") +
#                          s(roughness_250m, k = 5, bs = "cr"),
#                        data = hab_model, method = "REML", family = binomial("logit"))
# summary(m_seagrass_250m)
# plot(m_seagrass_250m, pages = 1, residuals = T, cex = 5)
# 
# # Sand
# out.all$sand[[1]]
# m_sand_250m <- gam(cbind(sand, total_pts - sand) ~
#                      s(aspect_250m, k = 5, bs = "cc")  +
#                      #s(detrended_250m, k = 5, bs = "cr")  +
#                      s(bathy_250m, k = 5, bs = "cr") +
#                      s(roughness_250m, k = 5, bs = "cr"),
#                    data = hab_model, method = "REML", family = binomial("logit"))
# summary(m_sand_250m)
# plot(m_sand_250m, pages = 1, residuals = T, cex = 5)


## predict, rasterise and plot ------------------------------------------------

preddf<- cbind(preddf,
               "preef"     = predict(fitted_models$reef, preddf, type = "response", se.fit = T),
               "psand"     = predict(fitted_models$sand, preddf, type = "response", se.fit = T),
               "pseagrass" = predict(fitted_models$seagrass, preddf, type = "response", se.fit = T)) %>%
  glimpse()

prasts <- rast(preddf %>% dplyr::select(x, y, preef.fit, psand.fit, pseagrass.fit))

prasts <- terra::crop(prasts, st_bbox(bbox_whole)); plot(prasts, range = c(0, 1))


## NOT DOING Removing parts of the prediction that are not observed ---------------------

# using MESS (Elith et al. 2010), we'll remove parts of the prediction that have
# combinations of predictors that were not observed. This makes sure we're not
# predicting overly confidently.

habi_sf <- st_as_sf(habi, coords = c("longitude", "latitude"), crs = 4326)
xy <- as.data.frame(habi_sf)
xy$geometry <- st_coordinates(habi_sf)

xy <- xy %>% 
  mutate(
    x = xy$geometry[, 1],
    y = xy$geometry[, 2]
  ) %>% 
  dplyr::select(x, y) %>% 
  glimpse()


### NOT DOING 250m prediction predicts::MESS --------------------------------------------

resp.vars <- c("preef", "psand", "pseagrass")
models <- list(fitted_models$reef, fitted_models$sand, fitted_models$seagrass); summary(models)

predhab <- preddf

# Extract predictor variables from the model's terms
model_terms <- terms(models[[1]])
model_vars <- attr(model_terms, "term.labels")

# Check the model variables
print(model_vars)
predictors_terra <- rast(predictors) # mess wants a terra object

# loop for each 250m habitat model
for(i in 1:length(resp.vars)) {
  # select a model from the list of models above
  print(resp.vars[i])
  mod <- models[[i]]
  
  # select the prediction raster
  temppred <- predhab %>% 
    dplyr::select(x, y, paste0(resp.vars[i], '.fit'),
                  paste0(resp.vars[i], '.se.fit')) %>% 
    rast(crs = "epsg:4326")
  
  # select the predictors of the model within the whole set of predictors
  model_vars <- attr(model_terms, "term.labels")
  dat <- terra::extract(subset(predictors_terra, model_vars), xy, ID = F)
  
  # run the MESS function
  messrast <- predicts::mess(subset(predictors_terra, model_vars), dat) %>% 
    terra::clamp(lower = -0.01, values = F)
  
  # remove unobserved areas
  messrast <- terra::crop(messrast, temppred)
  temppred_m <- terra::mask(temppred, messrast)
  
  # then add the final rasters together
  if (i == 1) {
    preddf_m <- temppred_m
  }
  else {
    preddf_m <- rast(list(preddf_m, temppred_m))
  }
}
plot(preddf_m) # check where things were removed


## Make a big file stack with the cleaned predictions and the bathy layers ----

#if not using MESS, run the chunk below: 
reef <- predhab %>% 
  dplyr::select(x, y, paste0("preef", '.fit'),
                paste0("preef", '.se.fit')) %>% 
  rast(crs = "epsg:4326")

sand <- predhab %>% 
  dplyr::select(x, y, paste0("psand", '.fit'),
                paste0("psand", '.se.fit')) %>% 
  rast(crs = "epsg:4326")

seagrass <- predhab %>% 
  dplyr::select(x, y, paste0("pseagrass", '.fit'),
                paste0("pseagrass", '.se.fit')) %>% 
  rast(crs = "epsg:4326")

predictors_stack <- list(aspect = predictors[["aspect"]], depth = predictors[["bathy"]], roughness = predictors[["roughness"]], detrended = predictors[["detrended"]],
                         reef = reef, sand = sand, seagrass = seagrass
)

#if using MESS, run below
# predictors_stack <- list(aspect = predictors[["aspect"]], depth = predictors[["bathy"]], roughness = predictors[["roughness"]], detrended = predictors[["detrended"]],
#                          reef = preddf_m[["preef.fit"]], sand = preddf_m[["psand.fit"]], seagrass = preddf_m[["pseagrass.fit"]]
# )

# Convert all raster::RasterLayer to SpatRaster
for (i in seq_along(predictors_stack)) {
  if (inherits(predictors_stack[[i]], "RasterLayer")) {
    predictors_stack[[i]] <- rast(predictors_stack[[i]])
  }
}


for (i in 1:length(predictors_stack)) { # Unifying extents and resampling so that they can be stacked
  reference_raster <- rast(predictors_stack[[1]]) # Use the raster with the highest resolution as a reference. Otherwise the re-sampling will apply lower resolutions to everything.
  predictors_stack[[i]] <- terra::resample(predictors_stack[[i]], reference_raster, method = "bilinear")
  cat(paste("New extent of", i, ":", ext(predictors_stack[[i]]), "\n"))
}

list2env(predictors_stack, envir = .GlobalEnv) # Put the elements of the list in the environment

predictors <- c(
  #depth, aspect, roughness, detrended, # 250m data
  reef, sand, seagrass                 # Habitat predictions using 250m bathymetry
)

plot(predictors[[1]])

#predictors <- raster::stack(predictors) # convert to the right format before saving
#predictors <-  rast(predictors) # convert to rast
saveRDS(predictors, file = "data/output_data/01_A_bathymetry_habitat_rasters.rds") # Takes a while if running with LiDAR


### END ###
