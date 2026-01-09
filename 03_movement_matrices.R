# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    Spatial layers from 01_Cleaning-files.R
# Task:    Set fish movement in grid cells
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    August 2024

# -----------------------------------------------------------------------------

# Status:  Mostly done, not sure whether movement is working ??

# -----------------------------------------------------------------------------

library(tidyverse) # for data manipulation
library(ggplot2) # for plotting
library(gridExtra) # for plot arranging
library(sf) # for dealing with shapefiles
library(sfnetworks) # to create geospatial networks

rm(list = ls())
par(mfrow = c(1, 1))

colour_palette <- eval(parse(text = readLines("yijarup_chapter_colours.txt")))

## Files used in this script --------------------------------------------------

file_water <- "data/output_data/02_watergrid.rds"

## Centroid Function ----------------------------------------------------------

# This returns the centre of the polygon, but if it's on land it will create a new centroid
st_centroid_within_poly <- function (poly) { # This returns the centre of the ploygon, but if it's on land it will create a new centroid
  
  # check if centroid is in polygon
  centroid <- poly %>% st_centroid() 
  in_poly <- st_within(centroid, poly, sparse = F)[[1]] 
  
  # if it is, return that centroid
  if (in_poly) return(centroid) 
  
  # if not, calculate a point on the surface and return that
  centroid_in_poly <- st_point_on_surface(poly) 
  return(centroid_in_poly)
}


## Load original files --------------------------------------------------------

water <- readRDS(file_water) %>% 
  st_make_valid()

ggplot(water) +
  geom_sf(aes(fill = status), colour = NA) +
  theme_void() +
  scale_fill_manual(values=c(colour_palette[4], colour_palette[6], colour_palette[5], colour_palette[1]))


## Calculate cell centroids and distance-to-neighbour-cells -------------------

# Get centroids for the grid cells - CHARLOTTE HAS A LOT MORE COLUMNS TO HER CENTROIDS DATASET ???
centroids <- st_centroid_within_poly(water)
plot(centroids[, !sapply(centroids, is.list)], cex=0.3) #plotting all but list-columns

# Get the number of cells in the model, this will allow to to calculate distances and habitat cover differences.
points <- as.data.frame(st_coordinates(centroids))%>%
  mutate(ID=row_number())
NCELL <- nrow(points)

# Convert the points in the centroids of the polygon to a spatial points file
points$ID <- as.integer(points$ID)
points_sf <- st_as_sf(points, coords = c("X", "Y"))
points_sp <- st_cast(st_geometry(points_sf), "POINT")

# Calculate the distance from each point to other points
dist.mat <- st_distance(points_sp)

# Get the IDS for neighbour cells
n.closest <- 8 # Number of neighbours, 6 if cell shape is hexagon, 8 if square
neighbours <- as.data.frame(array(0, dim=c(NCELL, n.closest)))

for (i in 1:n.closest){
  neighbours[,i] <- apply(dist.mat, 1, function(x) {
    order(x, decreasing=F)[i+1] })
} # So you end up with a distance matrix of each grid cell centroid's distance to its 8 nearest neighbours


# Give the neighbouring points geometry based on the original set of points
point.list <- list()
for (i in 1:n.closest){
  
  temp1 <- as.data.frame(neighbours[,i])
  
  temp2 <- temp1 %>%
    rename(ID = "neighbours[, i]") %>% 
    inner_join(., points, by="ID") %>% 
    st_as_sf(., coords = c("X", "Y")) 
  
  temp3 <- st_cast(st_geometry(temp2), "POINT")
  
  point.list[[i]] <- temp3
  
}


## Connect points to neighbours in a network ----------------------------------
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

connected <- st_combine(c(multilinestrings[[1]], multilinestrings[[2]], multilinestrings[[3]], multilinestrings[[4]], multilinestrings[[5]], multilinestrings[[6]]))
connected <- st_cast(connected, "LINESTRING") # Needs to be a line string rather than multiline for the next step
plot(connected)

st_write(connected, "data/output_data/03_network_shapefile.shp", delete_layer = T)


## Set up an sf network and a distance matrix ---------------------------------

network <- as_sfnetwork(connected, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

# Calculate the distances from each point to every other point on the network
net <- activate(network, "nodes")
network_matrix <- st_network_cost(net, from=points_sf, to=points_sf)
network_matrix <- network_matrix * 111 # Multiple by 111 to get from degrees to kms
dim(network_matrix); nrow(points_sf) # Check that the dimensions match up to how many points you think you should have in the network

# Checking that short distances between points are the same e.g. 1 -> 2
test_point_1 <- 1
test_point_2 <- 2
plot(points_sf, col = ifelse(points_sf$ID == test_point_1 | points_sf$ID == test_point_2, "red", "gray"), 
     main = paste0("Distance between red points is ", round((st_distance(points_sf[test_point_1,2], points_sf[test_point_2,2])) * 111), "km, which is realistic."))

# Distance should not be measured in a straight line, it should go around land (if there is land between two points).
test_point_1 <- 10
test_point_2 <- max(points_sf$ID)
plot(points_sf, col = ifelse(points_sf$ID == test_point_2 | points_sf$ID == test_point_1, "red", "gray"), 
     main = paste0("Distance between red points is ", round((st_distance(points_sf[test_point_2,2], points_sf[test_point_1,2])) * 111), "km, which is *not* realistic."))
mtext(paste0("Distance around land is ", round((st_network_cost(net, from=points_sf[test_point_2,2], to=points_sf[test_point_1,2])) * 111), " km",
             ", Distance in a straight line is ", round((st_distance(points_sf[test_point_2,2], points_sf[test_point_1,2])) * 111), " km"), 
      side = 3, line = -1.5, at = par("usr")[1] + 0.4*diff(par("usr")[1:2]),
      cex = 1.2, col = colour_palette[5])


## Adding habitat to the cells ------------------------------------------------

# the habitat predictions are probabilities, and we want to normalise them so they add up to 1.
habitat_cols <- c("reef", "seagrass", "sand")
summary(rowSums(st_drop_geometry(water[habitat_cols]))) # here each cell's habitat doesn't add up to 1

# Normalize each row so that the sum of habitat types equals 1
water[habitat_cols] <- as.data.frame(st_drop_geometry(water[habitat_cols])) / rowSums(st_drop_geometry(water[habitat_cols]))

summary(rowSums(st_drop_geometry(water[habitat_cols]))) # and now they do!
plot(water[, !sapply(water, is.list)]) # plot all but list-columns

## Save files to use in the next step -----------------------------------------

saveRDS(network_matrix, file = "data/output_data/03_network_matrix.rds")
saveRDS(water, file="data/output_data/03_water.rds")

## Create connectivity matrix for fish movement -------------------------------

# Assume the fish will try and swim the shortest path between locations.
# Calculate the probability a fish moves to this site in a given time step using a swimming speed.
# This creates a dispersal kernel based on the negative exponential distribution.

network_matrix <- readRDS("data/output_data/03_network_matrix.rds")

pDist <- matrix(NA, ncol=NCELL, nrow=NCELL)
for(r in 1:NCELL){
  for(c in 1:NCELL){
    p <- network_matrix[r,c]*1
    pDist[r,c] <- p
  } 
} # this loop compares the distance of cell 1 with every other cell, cell 2 with every other cell, cell 3....

# Calculate the difference in habitat types between each of the cells i.e. will there be an increase in reef % if you go from cell 1 to cell 2
habitat_types <- c("reef", "seagrass", "sand")
habitat_perc <- list(
  reef = water$reef,
  seagrass = water$seagrass,
  sand = water$sand
)
glimpse(habitat_perc)

# We want to calculate the difference in habitat cover between each grid cell. This will determine how likely a fish is to move from one grid cell to another.
p_habitat <- list() # Initialize a list to store the matrices

# loop over each habitat type in habitat_perc
for (habitat in names(habitat_perc)) {
  
  print(habitat) # see where you are
  p_matrix <- matrix(NA, ncol = NCELL, nrow = NCELL) # Create an empty matrix
  
  # compute the difference matrix
  for (r in 1:NCELL) { 
    for (c in 1:NCELL) {
      p_matrix[r, c] <- as.numeric(habitat_perc[[habitat]][c] - habitat_perc[[habitat]][r])
    } 
  } 
  
  # store the matrix in the list
  p_habitat[[habitat]] <- p_matrix
} # the loop compare the habitat cover of cell 1 to every other cell, then cell 2 to every other cell, then cell 3, etc. for every habitat type
glimpse(p_habitat)

## Save all the files you need to remake these matrices -----------------------
saveRDS(pDist, "data/output_data/03_pDist.rds")
saveRDS(p_habitat, "data/output_data/03_p_habitat.rds")


## Create adult movement probability using utility function -------------------

# Load files if need be
pDist <- readRDS("data/output_data/03_pDist.rds")

list2env(habitat_perc, envir = .GlobalEnv)# this brings all the list objects into the environment (seagrass, sand, reef)
water <- water %>% mutate(cell_index = seq_len(nrow(water))) # this gives each row an ID that matches with the matrices' ID. Otherwise the movements make no sense

# First determine the utility of each of the cells 
# This is very sensitive to changes in the habitat values particularly for reef 

## Small movement Swim Speed = 2.5 95% within approx 10km 
## Medium movement Swim Speed = 5 95% within approx 25km
## Big movement Swim Speed = 10 95% within approx 45km 

swim_speed_adult <- 10

a = -(1 / swim_speed_adult)
b =  0.150      # Attractiveness of reef habitat
c =  0.010      # Attractiveness of seagrass habitat
d =  0.005      # Attractiveness of sand habitat
e =  0.2        # attractiveness of cockburn sound cells

# From the attractivity of each habitat type, and the cell's distance to other cells, each cell's attractiveness is calculated (based on the % of each habitat in that cell)
adult_hab_attractivity <- (a * pDist) + (b * reef) + (c * seagrass) + (d * sand)
glimpse(adult_hab_attractivity)

# Calculate the summed utility across the rows 
rowU <- matrix(NA, ncol = 1, nrow = NCELL)
cell_utility <- matrix(NA, ncol = NCELL, nrow = NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    U <- exp(adult_hab_attractivity[r, c])
    cell_utility[r, c] <- U
  }
} # this loop calculates the likelihood of moving from cell x to any other cell based on its attractivity and distance to it.
glimpse(cell_utility)
rowU <- as.data.frame(rowSums(cell_utility))

# Calculate the probability that the fish will move to this site
adult_cell_movement_probability <- matrix(NA, ncol = NCELL, nrow = NCELL)
for (r in 1:NCELL){
  for (c in 1:NCELL){
    adult_cell_movement_probability[r, c] <- (exp(adult_hab_attractivity[r, c]))/rowU[r, 1]
  }
} # this loop calculates the probability of moving to a certain cell based on all other possible moves.
rowSums(adult_cell_movement_probability) # should be full of 1, because cell 1's probability of moving to any other cell (all the row) is 1.
sum(is.na(adult_cell_movement_probability)) # There should be no NAs, otherwise the model can't calculate things correctly.
adult_cell_movement_probability[1:10, 1:10]


# We'll look at whether the movement makes sense.
num_samples <- 3  # Number of cells to visualize
plot_list <- list()  # Store all plots here

for (i in 1:num_samples) {
  random_point <- sample(1:NCELL, 1)  # Pick a random cell
  movement <- adult_cell_movement_probability[random_point, ]
  
  # Create a dataframe with movement probabilities
  water_2 <- water %>%
    mutate(movement_prob = movement,
           test_point = (cell_index == random_point))  # use cell_index from earlier

  # Movement probability map
  movement_plot <- ggplot() +
    geom_sf(data = water_2, aes(fill = movement_prob), color = NA, lwd = 0) +
    geom_sf(data = water_2 %>% filter(test_point), fill = "red", color = NA, lwd = 0) +
    scale_fill_gradient(low = "#EAD1DC", high = "#B95F89") +
    ggtitle(paste("Movement Prob. - Cell", random_point)) +
    theme_void()
  
  # Cumulative probability plot
  data <- as.data.frame(movement)
  colnames(data) <- "Prob"
  
  data <- data %>%
    mutate(Distance = pDist[random_point, ]) %>%
    arrange(Distance) %>%
    mutate(Cumulative_Prob = cumsum(Prob))
  
  cumulative_plot <- ggplot(data, aes(x = Distance, y = Cumulative_Prob)) +
    geom_line(color = "#B95F89") +
    ggtitle(paste("Cumulative Prob. - Cell", random_point)) +
    xlab("Distance") + ylab("Cumulative Probability") +
    theme_minimal()
  
  # Store both plots in list
  plot_list[[2 * i - 1]] <- movement_plot
  plot_list[[2 * i]] <- cumulative_plot
} # this loop creates, for a random cell, a map of which cells are most likely to be travelled to, and a plot of how cumulative probability of travel by distance

# Plot all in a 2-row, 3-column layout (fits 6 plots, 3 pairs)
do.call(grid.arrange, c(plot_list, ncol = 2, nrow = 3))


## Recruitment matrix ---------------------------------------------------------

# we want the recruits to be in seagrass, and then move out from there 
cockburn_list <- water$ID[water$type == "shore_north_cockburn_warnbro" | 
                            water$type =="north_cockburn_warnbro"]

dispersal <- as.data.frame(seagrass) %>%
  rename(perc_habitat = seagrass) %>%
  mutate(ID = water$ID, # so that it lines up with the seagrass ID and also cockburn sound cells
         perc_habitat = ifelse(is.na(perc_habitat), 0, perc_habitat),
         cockburn = ifelse(ID %in% cockburn_list, 1, 0)) # replace NAs with zeroes for now.
dispersal$area_km2 <- as.vector(water$cell_area*0.000001) # convert cell_area to km2
glimpse(dispersal)
sum(is.na(dispersal$perc_habitat)) # no NAs, all good.

dispersal <- dispersal %>% 
  #filter(perc_habitat != 0) %>% 
  mutate(area_km2 = area_km2 * 0.1, # why??
         perc_habitat = perc_habitat * 0.5) %>% 
  glimpse()

# check that water type is correct
ggplot() +
  geom_sf(data = water, aes(fill = dispersal$cockburn))


recruitment <- array(0, dim = c(nrow(dispersal), 2))
recruitment[ ,2] <- as.numeric(dispersal$ID) 

cell_utility <- matrix(0, ncol = 2, nrow=(nrow(dispersal)))
cell_utility[,2] <- as.numeric(dispersal$ID) 

for(cell in 1:nrow(recruitment)){
  U <- exp(dispersal[cell, "perc_habitat"]) + 0.25*exp(dispersal[cell, "cockburn"]) # each cell is exponentially more attractive the more seagrass it has, AND if it's a cockburn sound cell.
  cell_utility[cell, 1] <- as.numeric(U)
} # this loop makes Cockburn sound cells and seagrass cells exponentially more attractive.

rowU <- as.data.frame(sum(cell_utility[,1]))

for (cell in 1:nrow(recruitment)){
  recruitment[cell, 1] <- cell_utility[cell, 1] / rowU[1, 1]
} # this loop calculates attractivity for each cell for juveniles. (proportion of that cell's attractivity to the total attractivity of the area)

# check that recruitment probability makes sense with where seagrass is predicted to be
water_recruitment <- water %>%
  left_join(as.data.frame(recruitment) %>%
              rename(recruitment_prob = V1,
                     ID = V2), by = "ID")

ggplot() +
  geom_sf(data = water_recruitment, aes(fill = recruitment_prob), color = NA, lwd = 0) +
  scale_fill_gradient(low = "#EAD1DC", high = "#B95F89") +
  ggtitle("Recruitment Probability Map") +
  theme_minimal() # looks okay.

recruitment <- as.vector(recruitment[,1])

# technically cockburn cells should contribute to 1.6x the batch fecundity of metro offshore and southwest cells (https://researchportal.murdoch.edu.au/esploro/outputs/report/Snapper-connectivity-and-evaluation-of-juvenile/991005792873207891)
# but recruitment probability is different from batch fecundity i think so skipping that for now.


## Recruit movement -----------------------------------------------------------

# Want the recruits to stay in seagrass until they mature and move to the reef

swim_speed_juv <- 5

a = -(1/swim_speed_juv)
b =  0.150     # Attractiveness of reef habitat
c =  0.010     # Attractiveness of seagrass habitat
d =  0.005     # Attractiveness of sand habitat

juv_hab_attractivity <- (a * pDist) + (b * reef) + (c * seagrass) + (d * sand)

# Calculate the summed utility across the rows 
rowU <- matrix(NA, ncol=1, nrow=NCELL)
cell_utility <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    U <- exp(juv_hab_attractivity[r,c])
    cell_utility[r,c] <- U
  }
} # this loop makes --??

rowU <- as.data.frame(rowSums(cell_utility))


# Calculate the probability that the fish will move to this site

juv_cell_movement_probability <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    juv_cell_movement_probability[r,c] <- (exp(juv_hab_attractivity[r,c]))/rowU[r,1]
  }
}
rowSums(juv_cell_movement_probability) # should be all 1, because the probability of moving to any other cell (all rows) is 1.

num_samples <- 3  # Number of cells to visualize
plot_list <- list()  # Store all plots here

for (i in 1:num_samples) {
  random_point <- sample(1:NCELL, 1)  # Pick a random cell
  movement <- juv_cell_movement_probability[random_point, ]
  
  # Create a dataframe with movement probabilities
  water_2 <- water %>%
    mutate(movement_prob = movement,
           test_point = (cell_index == random_point))  # use cell_index from earlier
  
  # Movement probability map
  movement_plot <- ggplot() +
    geom_sf(data = water_2, aes(fill = movement_prob), color = NA, lwd = 0) +
    geom_sf(data = water_2 %>% filter(test_point), fill = "red", color = NA, lwd = 0) +
    scale_fill_gradient(low = "#EAD1DC", high = "#B95F89") +
    ggtitle(paste("Movement Prob. - Cell", random_point)) +
    theme_void()
  
  # Cumulative probability plot
  data <- as.data.frame(movement)
  colnames(data) <- "Prob"
  
  data <- data %>%
    mutate(Distance = pDist[random_point, ]) %>%
    arrange(Distance) %>%
    mutate(Cumulative_Prob = cumsum(Prob))
  
  cumulative_plot <- ggplot(data, aes(x = Distance, y = Cumulative_Prob)) +
    geom_line(color = "#B95F89") +
    ggtitle(paste("Cumulative Prob. - Cell", random_point)) +
    xlab("Distance") + ylab("Cumulative Probability") +
    theme_minimal()
  
  # Store both plots in list
  plot_list[[2 * i - 1]] <- movement_plot
  plot_list[[2 * i]] <- cumulative_plot
} # this loop creates, for a random cell, a map of which cells are most likely to be travelled to, and a plot of how cumulative probability of travel by distance

# Plot all in a 2-row, 3-column layout (fits 6 plots, 3 pairs)
do.call(grid.arrange, c(plot_list, ncol = 2, nrow = 3))


## Save files for next step ---------------------------------------------------

saveRDS(adult_cell_movement_probability, paste0("data/output_data/03_adult_movement_", swim_speed_adult, "_swim_speed.rds"))
saveRDS(juv_cell_movement_probability, paste0("data/output_data/03_juv_movement_", swim_speed_juv, "_swim_speed.rds"))
saveRDS(recruitment, "data/output_data/03_recruitment.rds")

### END ###