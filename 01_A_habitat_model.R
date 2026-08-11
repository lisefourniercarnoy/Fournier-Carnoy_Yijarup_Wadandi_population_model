# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    MEGlab BOSS and BRUV habitat data
# Task:    Create a habitat prediction
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    August 2026

# -----------------------------------------------------------------------------

# Status:  Complete and final. 

# -----------------------------------------------------------------------------

rm(list = ls()) # clear working environment

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))

# load libraries
library(tidyverse) # for data manipulation
library(FSSgam) # for selecting the models
library(mgcv) # for making the models
library(predicts) # for making the prediction
library(terra) # for extracting bathy values at opcode locations
library(tidyterra) # for plotting rasters
library(raster) # to deal with rasters
library(sf) # for dealing with polygons
library(stars) # for detrended bathymetry
library(starsExtra) # for detrended bathymetry


## 0. Files used in this script -----------------------------------------------

file_obs_hab            <- "data/external_data/sw-network_broad-habitat.csv" # Observed habitat values, cleaned in 01_prep
file_bathy_predictors   <- "data/external_data/sw-network_bathymetry.rds" # made in 01_prep

# -- make a common extent of the area we're predicting to, to use throughout the script
bbox_whole <- st_bbox(c(xmin = 114.4, ymin = -34.75, xmax = 116.0, ymax = -31), crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(4326)
bbox_whole <- extent(st_bbox(bbox_whole))


## 1. Load data + make derivatives --------------------------------------------

habi <- read.csv(file_obs_hab) %>% 
  mutate(reef = inverts + macroalgae + rock) %>% 
  dplyr::select(-c(inverts, macroalgae, rock, campaignid)) %>% 
  glimpse()
habi$diff <- (habi$reef + habi$seagrass + habi$sand) - habi$total_pts
habi$diff # There should only be zeroes.


# extract bathymetry values to the location of the opcodes
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

# check for correlation of predictor variables- remove anything highly correlated (>0.95)
round(cor(hab[, pred.vars]), 2) # All good, no high correlations

# check to make sure Response vector has not more than 80% zeros
unique.vars = unique(as.character(hab$taxa))
head(hab)
unique.vars.use = character()
for(i in 1:length(unique.vars)){
  temp.dat = hab[which(hab$taxa == unique.vars[i]),]
  if(length(which(temp.dat$taxa == 0))/nrow(temp.dat)<0.9){
    unique.vars.use = c(unique.vars.use,unique.vars[i])}
}

unique.vars.use # all remain


## 2. Run the full subset model selection -------------------------------------

outdir    <- ("data/output_data/01_A_habitat_modelling_outputs/")
resp.vars <- unique.vars.use
out.all   <- list()
var.imp   <- list()


### -- 2.1 Loop through the FSS function for each Abiotic taxa ----------------

hab <- hab %>% na.omit() # removing lidar columns otherwise it wont run
colSums(is.na(hab)) # should be all zeroes

for(i in 1:length(resp.vars)){
  print(resp.vars[i])
  use.dat <- hab[hab$taxa == resp.vars[i],]
  use.dat   <- as.data.frame(use.dat)
  
  # basic model to compare all the other combinations with
  Model1  <- gam(cbind(number, (total_pts - number)) ~ # Success and failure counts for each habitat
                   s(bathy, bs = 'cr', k = 3),
                 family = binomial(link = "logit"),  data = use.dat)
  
  # generate variable combinations to test in models
  model.set <- generate.model.set(use.dat = use.dat,
                                  test.fit = Model1,
                                  pred.vars.cont = pred.vars,
                                  cyclic.vars = c("aspect"),
                                  k = 5,
                                  cov.cutoff = 0.7,
                                  max.predictors = 3
  )
  
  # fit models to all variable combinations
  out.list <- fit.model.set(model.set,
                            max.models = 600,
                            parallel = T)
  names(out.list)
  
  # examine outcomes
  out.list$failed.models # examine the list of failed models
  mod.table <- out.list$mod.data.out  # look at the model selection table
  mod.table <- mod.table[order(mod.table$AICc), ]
  mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
  out.i     <- mod.table[which(mod.table$delta.AICc <= 2), ] # Top models
  out.all   <- c(out.all, list(out.i))
  var.imp   <- c(var.imp, list(out.list$variable.importance$aic$variable.weights.raw))
  
  # plot the top models
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


### -- 2.2 Model fits and importance ------------------------------------------

names(out.all) <- resp.vars
names(var.imp) <- resp.vars
all.mod.fits <- list_rbind(out.all, names_to = "response")
all.var.imp  <- do.call("rbind", var.imp)
out.all


## 3. Reload & format the data ------------------------------------------------

hab_model <- hab %>% 
  pivot_wider(names_from = taxa, values_from = number, values_fill = list(number = 0)) %>% 
  glimpse()

# check that the habitat points add up correctly, otherwise the model mis-counts.
hab_model$diff <- (hab_model$reef + hab_model$seagrass + hab_model$sand) - hab_model$total_pts
hab_model$diff[hab_model$diff != 0] # there should only be zeroes.
hab_model <- hab_model[hab_model$diff == 0, ] # a few points don't add up right, so im just removing them because i can't figure out what happened, and they're not crucial

preds <- predictors; plot(preds)
predictions <- as.data.frame(preds, xy = TRUE, na.rm = TRUE) %>% glimpse()


## 4. Make models for each habitat and predict extent -------------------------

fam <- binomial("logit")

make_vars <- function(vars) {
  sapply(vars, function(var) {
    bs_type <- ifelse(grepl("aspect", var), "cc", "cr")
    paste0("s(", var, ", k = 5, bs = '", bs_type, "')")
  }) |> paste(collapse = " + ")
}

# make a function to fit the optimal model
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

fitted_models <- list() # object for storing the fitted models

for (habitat in names(out.all)) {
  
  # extract best model variable names
  best_model_name <- out.all[[habitat]]$modname[1]
  best_model_vars <- strsplit(best_model_name, "\\+")[[1]]
  
  # convert variables into GAM smooth terms
  best_vars <- make_vars(best_model_vars)
  
  # fit the model using cbind(successes, failures)
  fitted_model <- make_optimal_model(
    response_var = habitat,
    predictor_terms = best_vars,
    data = hab_model
  )
  
  # store model
  fitted_models[[habitat]] <- fitted_model
  
  # print quick summary
  message("Fitted binomial model for habitat: ", habitat)
  print(summary(fitted_model))
} # -- this loop fits the best model (identified in step 1.2.2) for each habitat type (reef, sand, seagrass)

plot(fitted_models$reef, page = 1, residuals = T, cex = 5)
plot(fitted_models$sand, page = 1, residuals = T, cex = 5)
plot(fitted_models$seagrass, page = 1, residuals = T, cex = 5)


## 5. Predict, rasterise and plot ---------------------------------------------

# -- for each raster cell in the 250m resolution raster, we'll predict each habitat, based on the bathy and derivatives

predictions <- cbind(predictions,
                "preef"     = predict(fitted_models$reef, predictions, type = "response", se.fit = T),
                "psand"     = predict(fitted_models$sand, predictions, type = "response", se.fit = T),
                "pseagrass" = predict(fitted_models$seagrass, predictions, type = "response", se.fit = T)) %>%
  glimpse()

pred_rast <- rast(predictions %>% dplyr::select(x, y, preef.fit, psand.fit, pseagrass.fit))

pred_rast <- terra::crop(pred_rast, st_bbox(bbox_whole)); plot(pred_rast, range = c(0, 1))


## 6. Make a file stack with predictions + bathy layers -----------------------

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

# convert all raster::RasterLayer to SpatRaster
for (i in seq_along(predictors_stack)) {
  if (inherits(predictors_stack[[i]], "RasterLayer")) {
    predictors_stack[[i]] <- rast(predictors_stack[[i]])
  }
}

# unify extents and resample so that they can be stacked
for (i in 1:length(predictors_stack)) { 
  reference_raster <- rast(predictors_stack[[1]])
  predictors_stack[[i]] <- terra::resample(predictors_stack[[i]], reference_raster, method = "bilinear")
  cat(paste("New extent of", i, ":", ext(predictors_stack[[i]]), "\n"))
}

list2env(predictors_stack, envir = .GlobalEnv) # put the elements of the list in the environment

predictors <- c(
  reef, sand, seagrass # habitat predictions using 250m bathymetry
)

plot(predictors[[1]])


## 7. Save outputs ------------------------------------------------------------

# CHECK
p <- ggplot() +
  geom_spatraster(data = pred_rast) +
  facet_wrap(~lyr, ncol = 3) +
  scale_fill_gradientn(colours = colour_palette[c(3:6)], na.value = "transparent") +
  coord_sf() +
  theme_minimal() +
  labs(fill = "Probability")
ggsave("plots/script_plot_checks/01_A/01_A_predicted_habitat.png", plot = p, width = 10, height = 6, dpi = 500)

# save
saveRDS(predictors, file = "data/output_data/01_A_bathymetry_habitat_rasters.rds")


### END ###
