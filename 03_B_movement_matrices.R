# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Spatial layers from 01_Cleaning-files.R
# Task:    Set fish movement in grid cells
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    August 2026

# -----------------------------------------------------------------------------

# Status:  Semi-finalised

# -----------------------------------------------------------------------------

rm(list = ls())

library(tidyverse) # for data manipulation
library(ggplot2) # for plotting
library(gridExtra) # for plot arranging
library(sf) # for dealing with shapefiles
library(sfnetworks) # to create geospatial networks

par(mfrow = c(1, 1))
colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))
common_crs = 7850 


## 0. Files used in this script -----------------------------------------------

file_water      <- "data/output_data/02_watergrid.rds"
file_land       <- "data/input_data/Q_aus_land_high_res_no_estuary.shp"
file_bathy      <- "data/input_data/SW_crop_AusBathyTopo__Australia__2024_250m_MSL_cog.tif"
file_hab_model  <- "data/output_data/03_A_habitat_affinity_outputs/03_A_model.rds"


## 1. Load original files -----------------------------------------------------

water <- readRDS(file_water) %>% 
  st_transform(common_crs) |> 
  st_make_valid()

land <- st_read(file_land) |> 
  st_transform(common_crs)

bathy <- raster::raster(file_bathy) |> 
  raster::projectRaster(crs = common_crs)


water <- water %>% 
  mutate(
    depth = exactextractr::exact_extract(bathy, water, 'mean') # extract depth info onto the grid
  ) |> 
  mutate(
    depth = ifelse(depth > 0, 0, -depth)
  ) # remove above-water depths, and make depth a positive (to match the habitat affinity data)

depth <- as.vector(water$depth)

ggplot(water) +
  geom_sf(aes(fill = depth), colour = NA) +
  theme_void()

reef <- water$reef
sand <- water$sand
seagrass <- water$seagrass


## 2. Setup the network for movement ------------------------------------------

# -- this consists in: figuring out how far cells are from each other,
# -- making a network to connect them

### 2.1 Calculate cell centroids and distance-to-neighbour-cells --------------

# -- create a function to calculate cell centroids

st_centroid_within_poly <- function (poly) { # This returns the centre of the ploygon, but if it's on land it will create a new centroid
  
  # check if centroid is in polygon
  centroid <- poly %>% st_centroid() 
  in_poly <- st_within(centroid, poly, sparse = F)[[1]] 
  
  # if it is, return that centroid
  if (in_poly) return(centroid) 
  
  # if not, calculate a point on the surface and return that
  centroid_in_poly <- st_point_on_surface(poly) 
  return(centroid_in_poly)
} # this returns the centre of the polygon, but if it's on land it will create a new centroid

# get centroids for the grid cells
centroids <- st_centroid_within_poly(water)
plot(centroids[, !sapply(centroids, is.list)], cex=0.3) # plotting all but list-columns

# get the number of cells in the model, this will allow to to calculate distances and habitat cover differences.
points <- as.data.frame(st_coordinates(centroids))%>%
  mutate(ID = row_number())
NCELL <- nrow(points)

# convert the points in the centroids of the polygon to a spatial points file
points$ID <- as.integer(points$ID)
points_sf <- st_as_sf(points, coords = c("X", "Y"), crs = st_crs(water)) %>%
  st_transform(common_crs)
points_sp <- st_cast(st_geometry(points_sf), "POINT")

# calculate the distance from each point to other points
dist_matrix_full <- st_distance(points_sp)
plot(points_sp)

# get the IDS for neighbour cells
n.closest <- 8 # number of neighbours, 6 if cell shape is hexagon, 8 if square
neighbours <- as.data.frame(array(0, dim = c(NCELL, n.closest)))

for (i in 1:n.closest){
  neighbours[,i] <- apply(dist_matrix_full, 1, function(x) {
    order(x, decreasing=F)[i+1] })
} # so you end up with a distance matrix of each grid cell centroid's distance to its 8 nearest neighbours


# give the neighbouring points geometry based on the original set of points
point.list <- list()
for (i in 1:n.closest){
  
  temp1 <- as.data.frame(neighbours[,i])
  
  temp2 <- temp1 %>%
    rename(ID = "neighbours[, i]") %>% 
    inner_join(., points, by="ID") %>% 
    st_as_sf(., coords = c("X", "Y"), crs = st_crs(water)) %>%
    st_transform(common_crs)
  
  temp3 <- st_cast(st_geometry(temp2), "POINT")
  
  point.list[[i]] <- temp3
  
}


### 2.2 Connect points to neighbours in a network -----------------------------

n <- nrow(points)

multilinestrings <- list()
for(i in 1:n.closest){
  
  linestring <- point.list[[i]]
  
  temp <- lapply(X = 1:n, FUN = function(x) {
    pair <- st_combine(c(points_sp[x], linestring[x]))
    line <- st_cast(pair, "LINESTRING")
    return(line)
  })
  
  temp2 <- st_multilinestring(do.call("rbind", temp))
  
  multilinestrings[[i]] <- temp2
  
} # this loop connects each cell to its neighbour cells

connected <- st_combine(c(multilinestrings[[1]], multilinestrings[[2]], multilinestrings[[3]], 
                          multilinestrings[[4]], multilinestrings[[5]], multilinestrings[[6]], 
                          multilinestrings[[7]], multilinestrings[[8]])) # should have as many objects to st_combine as n.closest
connected <- st_cast(connected, "LINESTRING") # needs to be a line string rather than multiline for the next step
plot(connected)

connected <- st_set_crs(connected, common_crs)

edges_sf <- st_as_sf(st_cast(connected, "LINESTRING")) %>%
  mutate(id = row_number())

st_write(connected, "data/output_data/03_B_network_shapefile.shp", delete_layer = T)


### 2.3 Set up an sf network and a distance matrix ----------------------------

network <- as_sfnetwork(connected, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

# calculate the distances from each point to every other point on the network
net <- activate(network, "nodes")
network_matrix <- st_network_cost(net, from=points_sf, to=points_sf)
network_matrix <- network_matrix / 1000  # convert to km if needed

dim(network_matrix); nrow(points_sf) # check that the dimensions match up to how many points you think you should have in the network

# checking that short distances between points are the same e.g. 1 -> 2
test_point_1 <- 1
test_point_2 <- 2
plot(points_sf, col = ifelse(points_sf$ID == test_point_1 | points_sf$ID == test_point_2, "red", "gray"), 
     main = paste0("Distance between red points is ", round((st_distance(points_sf[test_point_1,2], points_sf[test_point_2,2]))/1000), "km, which is realistic."))

# distance should not be measured in a straight line, it should go around land (if there is land between two points).
test_point_1 <- 500
test_point_2 <- max(points_sf$ID)
plot(points_sf, col = ifelse(points_sf$ID == test_point_2 | points_sf$ID == test_point_1, "red", "gray"))
mtext(paste0("Distance around land is ", round((st_network_cost(net, from=points_sf[test_point_2,2], to=points_sf[test_point_1,2])) / 1000), " km",
             ", Distance in a straight line is ", round((st_distance(points_sf[test_point_2,2], points_sf[test_point_1,2])) / 1000), " km"), 
      side = 3, line = -1.5, at = par("usr")[1] + 0.4*diff(par("usr")[1:2]),
      cex = 1.2, col = colour_palette[5])


### 2.4 Adding habitat to the cells -------------------------------------------

# the habitat predictions are probabilities, and we want to normalise them so they add up to 1.
habitat_cols <- c("reef", "seagrass", "sand")
summary(rowSums(st_drop_geometry(water[habitat_cols]))) # here each cell's habitat doesn't add up to 1

# normalize each row so that the sum of habitat types equals 1
water[habitat_cols] <- as.data.frame(st_drop_geometry(water[habitat_cols])) / rowSums(st_drop_geometry(water[habitat_cols]))

summary(rowSums(st_drop_geometry(water[habitat_cols]))) # and now they do!


### 2.5 Save files to use in the next step ------------------------------------

saveRDS(network_matrix, file = "data/output_data/03_B_network_matrix.rds")
saveRDS(water, file = "data/output_data/03_B_water.rds")


### 2.6 Create connectivity matrix for fish movement --------------------------

# -- fish will tend to move from cells with habitat they have low affinity with, to cells with habitat they have high affinity with.

network_matrix <- readRDS("data/output_data/03_B_network_matrix.rds")

# # below is hashed out because it takes a while to run (5-10min). rerun if needed.
# 
# pDist <- matrix(NA, ncol=NCELL, nrow=NCELL)
# for(r in 1:NCELL){
#   for(c in 1:NCELL){
#     p <- network_matrix[r,c]
#     pDist[r,c] <- p
#   }
# } # this loop compares the distance of cell 1 with every other cell, cell 2 with every other cell, cell 3....


### 2.7 Save all the files you need to remake these matrices ------------------

# saveRDS(pDist, "data/output_data/03_B_pDist.rds")


## 3. Fish movement -----------------------------------------------------------

# -- assume the fish will try and swim the shortest path between locations.
# -- calculate the probability a fish moves to this site in a given time step using a swimming speed.
# -- this creates a dispersal kernel based on the negative exponential distribution.

### 3.1 Non-migrating adult movement probability ------------------------------

water <- water %>% mutate(cell_index = seq_len(nrow(water))) # this gives each row an ID that matches with the matrices' ID. Otherwise the movements make no sense

# -- movement probability depends on three things:

# 1. habitat affinity
hab_aff_mod <- readRDS(file_hab_model)
summary(hab_aff_mod)
hab_aff <- predict(hab_aff_mod, # this will use colnames that are in the model predictors to make a prediction about how many fish we can expect to see in this cell. 
                   newdata = water %>% mutate(depth2 = depth^2,
                                              depth_m = depth,
                                              size_class = "n_mature",
                                              preef.fit = reef,
                                              psand.fit = sand,
                                              pseagrass.fit = seagrass
                                              ),
                   type = "response")
#list2env(habitat_perc, envir = .GlobalEnv)# this brings all the list objects into the environment (seagrass, sand, reef)


# 2. distance from other cells
pDist <- readRDS("data/output_data/03_B_pDist.rds")


# 3. swimming speed 
# -- Small movement Swim Speed = 2.5 95% within approx 10km 
# -- Medium movement Swim Speed = 5 95% within approx 25km
# -- Big movement Swim Speed = 10 95% within approx 45km 
swim_speed_adult <- 10
a = -(1 / swim_speed_adult)

# difference in habitat affinity from one cell to all others
hab_aff_diff <- matrix(0,nrow = NCELL, ncol = NCELL)
for (i in 1:NCELL) {
  for (j in 1:NCELL) {
    
    from_aff <- hab_aff[i]
    to_aff <- hab_aff[j]
    
    hab_aff_diff[i,j] <- to_aff - from_aff
  }
} # this loop calculates the difference in habitat affinity between cells.

# from this we can determine the utility of each of the cells.
# -- this is very sensitive to changes in the habitat values.
adult_hab_attractivity <- hab_aff_diff + (a*pDist) # here we need to make cells with suitable habitat MORE attractive, to square it. i think if you dont, the habitat makes very little difference as to where the population is.

# calculate the summed utility across the rows 
rowU <- matrix(NA, ncol = 1, nrow = NCELL)
cell_utility <- matrix(NA, ncol = NCELL, nrow = NCELL)

water$cell_area <- as.numeric(water$cell_area) * 1e-6 + 1 # convert m² to km², add a constant

# -- difference in attractivity between a cell of 1km2 and 2km2 (+1) is huge (doubling) but a +1 in area from 18km2 to 19km2 is much less (proportionately)
# -- logging accounts for this different relationship. otherwise large cells are wayyyy too attractive.
cell_area_km2 <- log(water$cell_area) 
cell_utility <- exp(adult_hab_attractivity) * matrix(cell_area_km2, 
                                                     nrow = NCELL, 
                                                     ncol = NCELL, 
                                                     byrow = TRUE) # this calculates the likelihood of moving from cell x to any other cell based on its attractivity and distance to it.
glimpse(cell_utility)

rowU <- as.data.frame(rowSums(cell_utility))
summary(rowU)

# CHECK: habitat affinity overall
water_test <- water
water_test$test <- hab_aff
summary(water_test$test)
p <- ggplot() +
  geom_sf(data = water_test, aes(fill = test), colour = NA) +
  scale_fill_gradient(low = colour_palette[3], high = colour_palette[5]) +
  theme_minimal()
ggsave("plots/script_plot_checks/03_B/03_B_adult_habitat_affinity.png", plot = p, width = 6, height = 10, dpi = 500)


# calculate the probability that the fish will move to this site
adult_cell_movement_probability <- matrix(NA, ncol = NCELL, nrow = NCELL)
adult_cell_movement_probability <- cell_utility / rowU[, 1] # this calculates the probability of moving to a certain cell based on all other possible moves.
rowSums(adult_cell_movement_probability) # should be full of 1, because cell 1's probability of moving to any other cell (all the row) is 1.
sum(is.na(adult_cell_movement_probability)) # there should be no NAs, otherwise the model can't calculate things correctly.


# CHECK: movement from test cells
num_samples <- 3  # Number of cells to visualize
plot_list <- list()  # Store all plots here

for (i in 1:num_samples) {
  random_point <- sample(1:NCELL, 1)  # Pick a random cell
  movement <- adult_cell_movement_probability[random_point, ]
  
  # create a dataframe with movement probabilities
  water_2 <- water %>%
    mutate(movement_prob = movement,
           test_point = (cell_index == random_point))  # use cell_index from earlier

  # movement probability map
  movement_plot <- ggplot() +
    geom_sf(data = water_2, aes(fill = movement_prob), color = NA, lwd = 0) +
    geom_sf(data = water_2 %>% filter(test_point), fill = "red", color = NA, lwd = 0) +
    scale_fill_gradient(low = colour_palette[3], high = colour_palette[5]) +
    ggtitle(paste("Movement Prob. - Cell", random_point)) +
    theme_minimal()

  # store both plot in list
  plot_list[[i]] <- movement_plot
} # this loop creates, for a random cell, a map of which cells are most likely to be travelled to, and a plot of how cumulative probability of travel by distance

# plot all in a 2-row, 3-column layout (fits 6 plots, 3 pairs)
p <- do.call(grid.arrange, c(plot_list, ncol = 3))
ggsave("plots/script_plot_checks/03_B/03_B_adult_movement_probability_test_cell.png", plot = p, width = 10, height = 7, dpi = 500)


### 3.2 Migrating adult movement probability ----------------------------------

# -- okay so for large adults (say 600mm+) they need to move to embayments to spawn...
# -- movement probability depends on three things:

# 1. habitat affinity
sg_list <- water$ID[water$spawning_status == TRUE]
spawning <- data.frame(
  perc_habitat = seagrass,
  ID = water$ID,
  area_km2 = as.vector(water$cell_area * 1e-6)
) %>% 
  mutate(
    spawning = ifelse(ID %in% sg_list, 1, 0)
  )
glimpse(spawning)
sum(is.na(spawning$perc_habitat)) # no NAs, all good.

ggplot() +
  geom_sf(data = water, aes(fill = spawning$spawning), col = NA)
ggplot() +
  geom_sf(data = water, aes(fill = spawning$perc_habitat), col = NA)


# okay now we make a map of affinity for adults during spawning
spawn_aff <- array(0, dim = c(nrow(spawning), 2))
spawn_aff[ ,2] <- as.numeric(spawning$ID) 

cell_utility <- matrix(0, ncol = 2, nrow=(nrow(spawning)))
cell_utility[,2] <- as.numeric(spawning$ID) 

for(cell in 1:nrow(spawn_aff)){
  U <- (exp(spawning[cell, "perc_habitat"]) + 0.5 * exp(spawning[cell, "spawning"]))# * spawning[cell, "area_km2"] # each cell is exponentially more attractive the more seagrass it has, AND if it's a spawning ground cell. also weighted by cell size.
  cell_utility[cell, 1] <- as.numeric(U)
} # this loop makes spawning ground cells and seagrass cells exponentially more attractive.

rowU <- as.data.frame(sum(cell_utility[,1]))

for (cell in 1:nrow(spawn_aff)){
  spawn_aff[cell, 1] <- cell_utility[cell, 1] / rowU[1, 1]
} # this loop calculates attractivity for each cell for juveniles. (proportion of that cell's attractivity to the total attractivity of the area)

# 2. distance from other cells
pDist <- readRDS("data/output_data/03_B_pDist.rds")


# 3. swimming speed 
# -- small movement swim speed = 2.5 95% within approx 10km 
# -- medium movement swim speed = 5 95% within approx 25km
# -- big movement swim speed = 10 95% within approx 45km 
swim_speed_spawner <- 10
a = -(1 / swim_speed_spawner)

# difference in habitat affinity from one cell to all others
hab_aff_diff <- matrix(0,nrow = NCELL, ncol = NCELL)
for (i in 1:NCELL) {
  for (j in 1:NCELL) {
    
    from_aff <- spawn_aff[i]
    to_aff <- spawn_aff[j]
    
    hab_aff_diff[i,j] <- to_aff - from_aff
  }
} # -- this loop calculates the difference in habitat affinity between cells.

# -- from this we can determine the utility of each of the cells.
# -- this is very sensitive to changes in the habitat values.

# -- here we need to make habitat affinity much more attractive. when checking the ranges of the two elements we're using to make the movement prob, hab_aff definitely needs to be bumped up.
summary(as.vector(hab_aff_diff))
summary(as.vector(a * pDist)) 

spawning_hab_attractivity <- 10000*(hab_aff_diff) + (a*pDist) # here we have to make spawning cells overpower the distance between cells

# calculate the summed utility across the rows 
rowU <- matrix(NA, ncol = 1, nrow = NCELL)
cell_utility <- matrix(NA, ncol = NCELL, nrow = NCELL)

water$cell_area <- as.numeric(water$cell_area) * 1e-6 + 1 # convert m² to km², add a constant

# difference in attractivity between a cell of 1km2 and 2km2 (+1) is huge (doubling) but a +1 in area from 18km2 to 19km2 is much less (proportionately)
# logging accounts for this different relationship. otherwise large cells are wayyyy too attractive.
cell_area_km2 <- log(water$cell_area) 
cell_utility <- exp(spawning_hab_attractivity) * matrix(cell_area_km2, 
                                                     nrow = NCELL, 
                                                     ncol = NCELL, 
                                                     byrow = TRUE) # this calculates the likelihood of moving from cell x to any other cell based on its attractivity and distance to it.
glimpse(cell_utility)

rowU <- as.data.frame(rowSums(cell_utility))
summary(rowU)

# CHECK: spawner habitat affinity overall
water_test <- water
water_test$test <- spawn_aff[,1]
p <- ggplot() +
  geom_sf(data = water_test, aes(fill = test), colour = NA) +
  scale_fill_gradient(low = colour_palette[3], high = colour_palette[5]) +
  theme_minimal()
ggsave("plots/script_plot_checks/03_B/03_B_spawner_habitat_affinity.png", plot = p, width = 6, height = 10, dpi = 500)


# calculate the probability that the fish will move to this site
spawning_cell_movement_probability <- matrix(NA, ncol = NCELL, nrow = NCELL)
spawning_cell_movement_probability <- cell_utility / rowU[, 1] # this calculates the probability of moving to a certain cell based on all other possible moves.
rowSums(spawning_cell_movement_probability) # should be full of 1, because cell 1's probability of moving to any other cell (all the row) is 1.
sum(is.na(spawning_cell_movement_probability)) # There should be no NAs, otherwise the model can't calculate things correctly.


# CHECK: spawner movement probability from test cells
num_samples <- 3  # Number of cells to visualize
plot_list <- list()  # Store all plots here

for (i in 1:num_samples) {
  random_point <- sample(1:NCELL, 1)  # Pick a random cell
  movement <- spawning_cell_movement_probability[random_point, ]
  
  # create a dataframe with movement probabilities
  water_2 <- water %>%
    mutate(movement_prob = movement,
           test_point = (cell_index == random_point))  # use cell_index from earlier
  
  # movement probability map
  movement_plot <- ggplot() +
    geom_sf(data = water_2, aes(fill = movement_prob), color = NA, lwd = 0) +
    geom_sf(data = water_2 %>% filter(test_point), fill = "red", color = NA, lwd = 0) +
    scale_fill_gradient(low = colour_palette[3], high = colour_palette[5]) +
    ggtitle(paste("Movement Prob. - Cell", random_point)) +
    theme_minimal()
  
  # store both plot in list
  plot_list[[i]] <- movement_plot
  
} # this loop creates, for a random cell, a map of which cells are most likely to be travelled to, and a plot of how cumulative probability of travel by distance

# plot all in a 2-row, 3-column layout (fits 6 plots, 3 pairs)
p <- do.call(grid.arrange, c(plot_list, ncol = 3))
ggsave("plots/script_plot_checks/03_B/03_B_spawner_movement_probability_test_cell.png", plot = p, width = 10, height = 7, dpi = 500)


### 3.3 Recruit probability ---------------------------------------------------

# -- we want the recruits to be in seagrass and spawning ground, and then move out from there 

# setup the correct elements
sg_list <- water$ID[water$spawning_status == TRUE]
dispersal <- data.frame(
  perc_habitat = seagrass,
  ID = water$ID,
  area_km2 = as.vector(water$cell_area*0.000001)
  ) %>% 
  mutate(
  spawning = ifelse(ID %in% sg_list, 1, 0)
)
glimpse(dispersal)
sum(is.na(dispersal$perc_habitat)) # no NAs, all good.

ggplot() +
  geom_sf(data = water, aes(fill = dispersal$spawning), col = NA)
ggplot() +
  geom_sf(data = water, aes(fill = dispersal$perc_habitat), col = NA)


# okay now we choose how much seagrass and spawning cells contribute to the map
recruitment <- array(0, dim = c(nrow(dispersal), 2))
recruitment[ ,2] <- as.numeric(dispersal$ID) 

cell_utility <- matrix(0, ncol = 2, nrow=(nrow(dispersal)))
cell_utility[,2] <- as.numeric(dispersal$ID) 

for(cell in 1:nrow(recruitment)){
  U <- (exp(dispersal[cell, "perc_habitat"]) + 0.5 * exp(dispersal[cell, "spawning"]))# * dispersal[cell, "area_km2"] # each cell is exponentially more attractive the more seagrass it has, AND if it's a spawning ground cell. also weighted by cell size.
  cell_utility[cell, 1] <- as.numeric(U)
} # this loop makes spawning ground cells and seagrass cells exponentially more attractive.

rowU <- as.data.frame(sum(cell_utility[,1]))

for (cell in 1:nrow(recruitment)){
  recruitment[cell, 1] <- cell_utility[cell, 1] / rowU[1, 1]
} # this loop calculates attractivity for each cell for juveniles. (proportion of that cell's attractivity to the total attractivity of the area)

# check that recruitment probability makes sense with where seagrass is predicted to be
water_recruitment <- water %>%
  left_join(as.data.frame(recruitment) %>%
              rename(recruitment_prob = V1,
                     ID = V2), by = "ID")

# CHECK: recruitment probability map
p <- ggplot() +
  geom_sf(data = water_recruitment, aes(fill = recruitment_prob), color = NA, lwd = 0) +
  scale_fill_gradient(low = colour_palette[3], high = colour_palette[5]) +
  theme_minimal() # looks okay.
ggsave("plots/script_plot_checks/03_B/03_B_recruitment_habitat_affinity.png", plot = p, width = 6, height = 10, dpi = 500)


# CHECK:: "Based on the percentage frequency of females in each length class, 
# -- the relative contribution of Cockburn Sound females, in terms of batch 
# -- fecundity, was ~1.6 times that of snapper in oceanic waters in the metro and 
# -- south-west" from https://wamsi.org.au/app/uploads/2025/02/WWMSP-4.1-Snapper-connectivity-and-juvenile-stocking.pdf#page=4.35

cs_spawn <- sum(water_recruitment$recruitment_prob[water_recruitment$spawning_status==TRUE])
non_cs_spawn <- sum(water_recruitment$recruitment_prob[water_recruitment$spawning_status==FALSE])

cs_spawn / non_cs_spawn
# okay so - problem here is that no matter how important i make CS cells, they don't reach the 1.6x thing we're lookin for
# putting this on the shelf for nowwwwwww...... that'll bite me in the arse later but it's a calculated decision..?


### 3.4 Juvenile movement -----------------------------------------------------

water <- water %>% mutate(cell_index = seq_len(nrow(water))) # this gives each row an ID that matches with the matrices' ID. Otherwise the movements make no sense

# -- movement probability depends on three things:

# 1. habitat affinity
hab_aff_mod <- readRDS(file_hab_model)
summary(hab_aff_mod)
hab_aff <- predict(hab_aff_mod, # this will use colnames that are in the model predictors to make a prediction about how many fish we can expect to see in this cell. 
                   newdata = water %>% mutate(depth2 = depth^2,
                                              depth_m = depth,
                                              size_class = "n_immature",
                                              preef.fit = reef,
                                              psand.fit = sand,
                                              pseagrass.fit = seagrass
                                              ),
                   type = "response")
#list2env(habitat_perc, envir = .GlobalEnv)# this brings all the list objects into the environment (seagrass, sand, reef)


# 2. distance from other cells
pDist <- readRDS("data/output_data/03_B_pDist.rds")


# 3. swimming speed 
# -- small movement Swim Speed = 2.5 95% within approx 10km 
# -- medium movement Swim Speed = 5 95% within approx 25km
# -- big movement Swim Speed = 10 95% within approx 45km 
swim_speed_juv <- 5
a = -(1/swim_speed_juv)


# difference in habitat affinity from one cell to all others
hab_aff_diff <- matrix(0,nrow = NCELL, ncol = NCELL)
for (i in 1:NCELL) {
  for (j in 1:NCELL) {
    
    from_aff <- hab_aff[i]
    to_aff <- hab_aff[j]
    
    hab_aff_diff[i,j] <- to_aff - from_aff
  }
} # this loop calculates the difference in habitat affinity between cells.

# from this we can determine the utility of each of the cells.
# -- this is very sensitive to changes in the habitat values.
juv_hab_attractivity <- hab_aff_diff + (a*pDist) # here we need to make cells with suitable habitat MORE attractive, to square it. i think if you dont, the habitat makes very little difference as to where the population is.

# calculate the summed utility across the rows 
rowU <- matrix(NA, ncol = 1, nrow = NCELL)
cell_utility <- matrix(NA, ncol = NCELL, nrow = NCELL)

water$cell_area <- as.numeric(water$cell_area) * 1e-6 + 1 # convert m² to km², add a constant

# -- difference in attractivity between a cell of 1km2 and 2km2 (+1) is huge (doubling) but a +1 in area from 18km2 to 19km2 is much less (proportionately)
# -- logging accounts for this different relationship. otherwise large cells are wayyyy too attractive.
cell_area_km2 <- log(water$cell_area) 
cell_utility <- exp(juv_hab_attractivity) * matrix(cell_area_km2, 
                                                   nrow = NCELL, 
                                                   ncol = NCELL, 
                                                   byrow = TRUE) # this calculates the likelihood of moving from cell x to any other cell based on its attractivity and distance to it.
glimpse(cell_utility)

rowU <- as.data.frame(rowSums(cell_utility))
summary(rowU)

# CHECK: habitat affinity overall
water_test <- water
water_test$test <- hab_aff
summary(water_test$test)
p <- ggplot() +
  geom_sf(data = water_test, aes(fill = test), colour = NA) +
  scale_fill_gradient(low = colour_palette[3], high = colour_palette[5]) +
  theme_minimal()
ggsave("plots/script_plot_checks/03_B/03_B_juvenile_habitat_affinity.png", plot = p, width = 6, height = 10, dpi = 500)


# calculate the probability that the fish will move to this site
juv_cell_movement_probability <- matrix(NA, ncol = NCELL, nrow = NCELL)
juv_cell_movement_probability <- cell_utility / rowU[, 1] # this calculates the probability of moving to a certain cell based on all other possible moves.
rowSums(juv_cell_movement_probability) # should be full of 1, because cell 1's probability of moving to any other cell (all the row) is 1.
sum(is.na(juv_cell_movement_probability)) # There should be no NAs, otherwise the model can't calculate things correctly.


# CHECK: movement from test cells
num_samples <- 3  # Number of cells to visualize
plot_list <- list()  # Store all plots here

for (i in 1:num_samples) {
  random_point <- sample(1:NCELL, 1)  # Pick a random cell
  movement <- juv_cell_movement_probability[random_point, ]
  
  # create a dataframe with movement probabilities
  water_2 <- water %>%
    mutate(movement_prob = movement,
           test_point = (cell_index == random_point))  # use cell_index from earlier
  
  # movement probability map
  movement_plot <- ggplot() +
    geom_sf(data = water_2, aes(fill = movement_prob), color = NA, lwd = 0) +
    geom_sf(data = water_2 %>% filter(test_point), fill = "red", color = NA, lwd = 0) +
    scale_fill_gradient(low = colour_palette[3], high = colour_palette[5]) +
    ggtitle(paste("Movement Prob. - Cell", random_point)) +
    theme_minimal()
  
  # store both plot in list
  plot_list[[i]] <- movement_plot
} # this loop creates, for a random cell, a map of which cells are most likely to be travelled to, and a plot of how cumulative probability of travel by distance

# plot all in a 2-row, 3-column layout (fits 6 plots, 3 pairs)
p <- do.call(grid.arrange, c(plot_list, ncol = 3))
ggsave("plots/script_plot_checks/03_B/03_B_juvenile_movement_probability_test_cell.png", plot = p, width = 10, height = 7, dpi = 500)


### 3.5 Save files for next step ----------------------------------------------

saveRDS(adult_cell_movement_probability, paste0("data/output_data/03_B_adult_movement_", swim_speed_adult, "_swim_speed.rds"))
saveRDS(spawning_cell_movement_probability, paste0("data/output_data/03_B_spawning_movement_", swim_speed_adult, "_swim_speed.rds"))
saveRDS(juv_cell_movement_probability, paste0("data/output_data/03_B_juv_movement_", swim_speed_juv, "_swim_speed.rds"))
saveRDS(recruitment, "data/output_data/03_B_recruitment.rds")

### END ###