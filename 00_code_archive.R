## code archive




## 01_cleaning_files making habitat polygons based on a threshold value

# Making polygons out of the rasters
all_polygons <- list()
par(mfrow = c(2, 3))
for (i in 1:nlyr(habitat_ras)) {

  class <- c(0, 0.5, 1) # Choose your inclusion threshold for the raster
  reclass <- classify(habitat_ras[[i]], rcl = class) # Classify the habitat layers
  polygons_filtered <- droplevels(reclass, level = class[1:2]) # Remove the lower level of habitat probability

  polygons <- as.polygons(polygons_filtered) # Make polygons with your classification
  plot(polygons, main = paste(names(habitat_ras[[i]]))) # Check to see if they look right

  # Save outputs
  output_filename <- paste0("Spatial_Data/polygon_", names(habitat_ras[[i]]), ".shp")
  st_write(st_as_sf(polygons), output_filename, delete_layer = TRUE)
  cat("Layer", i, "saved to", output_filename, "\n")
  all_polygons[[i]] <- st_as_sf(polygons)
}
graphics.off()

# Combine all polygons into a single sf object
habitat <- bind_rows(all_polygons)
habitat <- habitat %>% # Formatting so it's usable for 03_movement-matrices
  mutate(type = case_when(
    !is.na(reef) ~ "reef",
    !is.na(invert) ~ "invert",
    !is.na(seagrass) ~ "seagrass",
    !is.na(macro) ~ "macro",
    !is.na(rock) ~ "rock",
    !is.na(sand) ~ "sand",
    TRUE ~ NA_character_ )) %>%
  dplyr::select(type, geometry) %>%
  st_transform(crs)

habitat
st_write(habitat, "Spatial_Data/full_habitat.shp", delete_layer = T)


# Merge polygons into a single feature
col_palette <- brewer.pal(n = nrow(habitat), name = "RdBu")
plot(st_geometry(habitat), main = "Combined Habitats", col = col_palette[row(habitat)])
plot(wa_map$geometry, col = "#eeeeeeff", add = T)

legend("topright",
       legend = paste(rownames(habitat)),
       fill = col_palette,
       title = "Polygons",
       cex = 0.8)
# There is some overlap between habitat polygons, may have to adjust things. Also some blank spaces to adjust for



# Code below is for adjusting the whole SW habitat predictions. File is too big to work with using Git
# I'm cropping it and saving so that it's workable. RUN AGAIN IF EXTENT OR HABITATS CHANGE.

#habitat_ras <- readRDS("Spatial_Data/sw-network_predicted-habitat.RDS") %>%
#  dplyr::select(x, y, p_reef.fit, p_inverts.fit, p_seagrass.fit, p_macro.fit, p_rock.fit, p_sand.fit) %>%
#  rast(crs = "EPSG:4326", extent = bbox)
#names(habitat_ras) <- c("reef", "invert", "seagrass", "macro", "rock", "sand") # Rename layer names to make it easier

#saveRDS(habitat_ras, file = "Spatial_Data/wadandi_predicted-habitat.RDS")
#plot(readRDS(file = "wadandi_predicted-habitat.RDS"))



## 01_cleaning_files whole grid creation


distance <- 500 # Buffer distance in meters (adjust as needed)

land_buffer <- st_buffer(wa_map, dist = distance) %>%
  st_make_valid() %>%
  st_crop(bbox) %>%
  st_transform(crs)
plot(land_buffer$geometry)

crop_buf <- st_difference(land_buffer$geometry, wa_map$geometry) %>%
  st_make_valid() %>%
  st_transform(crs)
plot(crop_buf)

SmallGrd <- st_make_grid(crop_buf, cellsize=0.1, square=T, crs = 4283) %>%
  st_crop(bbox) %>% 
  st_make_valid() %>%
  st_transform(crs)

plot(SmallGrd)

SmallGrd <- st_intersection(SmallGrd, crop_buf) %>% 
  st_make_valid() %>%
  st_transform(crs)
plot(SmallGrd)


# Merge the small grid with the large grid
BigGrd <- st_difference(water, SmallGrd) %>% # takes a while to run
  st_make_valid() %>%
  st_transform(crs)
plot(water)
plot(BigGrd)
plot(SmallGrd, add=T)

str(SmallGrd)

water <- st_union(BigGrd, SmallGrd) %>% 
  st_make_valid()

plot(water) # check that it looks right

water_combined <- st_union(BigGrd, SmallGrd) %>% 
  st_make_valid()
plot(water_combined)

water_tbl <- water %>% 
  st_as_sf() %>% 
  mutate(cell_area=st_area(.)) %>% 
  filter(as.numeric(cell_area)>1) %>%  # get rid of any tiny cells that are less than 1m^2
  glimpse()

water <- water_tbl %>% 
  dplyr::select(x)

water_cleaned <- water_sf %>%
  filter(st_geometry_type(geometry) %in% c("POLYGON", "MULTIPOLYGON"))

plot(water_cleaned$geometry)

st_write(st_as_sf(water), "Spatial_Data/water.shp", delete_layer = T)

## making habitat types


# Calculate the percentage of each habitat type in each of the different cells
habitat.types <- c("reef", "invert", "seagrass", "macro", "rock", "sand")

habitat_inter <- list()

for(i in 1:length(habitat.types)){
  
  hab <- water %>% 
    filter(type %in% habitat.types[i])
  
  inter <- st_intersection(water, hab) %>% 
    st_make_valid() 
  #mutate(ID = as.numeric(water$ID))
  inter$hab_area <- as.numeric(inter$cell_area)*0.000001
  
  habitat_inter[[i]] <- inter
  
} # We end up with data that tells us for each cell (ID), what habitat is in it, and how much area it makes up in the cell.

# Let's check whether the cell IDs match between that new habitat area data and the water we used previously
check_hab <- rbind(habitat_inter[[1]], habitat_inter[[2]], habitat_inter[[3]], 
                   habitat_inter[[4]], habitat_inter[[5]], habitat_inter[[6]]) %>% 
  mutate(cell_area = as.numeric(cell_area*0.000001))

water1 <- water %>% # Making a new one with cell area in km2 
  mutate(cell_area = cell_area*0.000001)

cell1_water <- water1 %>% 
  filter(ID == 9) # Selecting ID for a visible cell
cell1_hab <- check_hab %>% 
  filter(ID == 9) # The same cell ID

plot(water1$geometry)
plot(cell1_water$geometry, add = T, col = "orange")
plot(cell1_hab$geometry, add = T, col = "red") # Orange and red cells should be the same.

# Let's calculate percentage of each habitat in each cell.
water_hab <- as.data.frame(check_hab) %>% 
  pivot_wider(id_cols=ID, names_from=type, values_from=hab_area) %>% 
  left_join(check_hab, water, by="ID") %>% 
  mutate(cell_area = (as.numeric(cell_area)*0.000001)) %>% 
  rowwise() %>% 
  mutate(hab_sum = (sum(reef, invert, seagrass, macro, rock, sand, na.rm = T)),
         correct_hab = ifelse(hab_sum>as.numeric(cell_area), "N", "Y"),
         difference = hab_sum-cell_area) %>% 
  ungroup()

check <- water_hab %>% 
  filter(reef>0)
plot(check$geometry) # Cells that have at least a little bit of reef.

hab_perc <- water_hab %>% 
  dplyr::select(-geometry) %>% 
  pivot_longer(cols = c(reef, invert, seagrass, macro, rock, sand), names_to = "type1", values_to = "hab_area1") %>% 
  mutate(perc_habitat = ((hab_area/cell_area)*100)) %>%  # This tells you how much of each water grid cell is made up of the different habitat types 
  mutate(perc_habitat = ifelse(perc_habitat > 100, 100, perc_habitat)) # The intersection has given a couple of places where the % is just over 100 so just round these back to 100

# Create separate data frames for each habitat type and fill in for where cells don't have a certain habitat type
perc_by_hab <- list()

for(i in 1:length(habitat.types)){
  
  temp <- hab_perc%>%
    filter(type==paste0(habitat.types[i]))%>%
    dplyr::select(ID, type, perc_habitat)%>%
    mutate(ID = as.numeric(ID))%>%
    complete(ID = 1:NCELL, fill = list(perc_habitat=0))%>%
    mutate(type = replace_na(habitat.types[i]))
  
  perc_by_hab[[i]] <- temp
}


#### Save files to use in the next step ####

saveRDS(network_matrix, file="Spatial_Data/network_matrix")
saveRDS(perc_by_hab[[1]], file="Spatial_Data/reef_perc")
saveRDS(perc_by_hab[[2]], file="Spatial_Data/invert_perc")
saveRDS(perc_by_hab[[3]], file="Spatial_Data/seagrass_perc")
saveRDS(perc_by_hab[[4]], file="Spatial_Data/macro_perc")
saveRDS(perc_by_hab[[5]], file="Spatial_Data/rock_perc")
saveRDS(perc_by_hab[[6]], file="Spatial_Data/sand_perc")

saveRDS(water, file="Spatial_Data/water")

#### Create connectivity matrix for fish movement ####

# Assume the fish will try and swim the shortest path between locations.
# Calculate the probability a fish moves to this site in a given time step using a swimming speed.
# This creates a dispersal kernel based on the negative exponential distribution.

network_matrix <- readRDS("Spatial_Data/network_matrix")
reef_perc <- readRDS("Spatial_Data/reef_perc")
invert_perc <- readRDS("Spatial_Data/invert_perc")
seagrass_perc <- readRDS("Spatial_Data/seagrass_perc")
macro_perc <- readRDS("Spatial_Data/macro_perc")
rock_perc <- readRDS("Spatial_Data/rock_perc")
sand_perc <- readRDS("Spatial_Data/sand_perc")

pDist <- matrix(NA, ncol=NCELL, nrow=NCELL)
for(r in 1:NCELL){
  for(c in 1:NCELL){
    p <- network_matrix[r,c]*1
    pDist[r,c] <- p
  }
}

## Calculate the difference in habitat types between each of the cells i.e. will there be an increase in reef % if you go from site 1 to site 2
habitat.types <- c("reef", "invert", "seagrass", "macro", "rock", "sand")
habitat.perc <- list(reef_perc, invert_perc, seagrass_perc, macro_perc, rock_perc, sand_perc)

# We want to calculate the difference in habitat cover between each grid cell. This will determine how likely a fish is to move from one grid cell to another.
# Each loop takes forever to run.
preef <- matrix(NA, ncol=NCELL, nrow=NCELL)
for (r in 1:NCELL){
  for (c in 1:NCELL){
    p <- as.numeric((reef_perc[c,3]) - (reef_perc[r,3]))
    preef[r,c] <- p
  }
}

pinvert <- matrix(NA, ncol=NCELL, nrow=NCELL)
for (r in 1:NCELL){
  for (c in 1:NCELL){
    p <- as.numeric((invert_perc[c,3]) - (invert_perc[r,3]))
    pinvert[r,c] <- p
  }
}

pseagrass <- matrix(NA, ncol=NCELL, nrow=NCELL)
for (r in 1:NCELL){
  for (c in 1:NCELL){
    p <- as.numeric((seagrass_perc[c,3]) - (seagrass_perc[r,3]))
    pseagrass[r,c] <- p
  }
}

pmacro<- matrix(NA, ncol=NCELL, nrow=NCELL)
for (r in 1:NCELL){
  for (c in 1:NCELL){
    p <- as.numeric((macro_perc[c,3]) - (macro_perc[r,3]))
    pmacro[r,c] <- p
  }
}

psand<- matrix(NA, ncol=NCELL, nrow=NCELL)
for (r in 1:NCELL){
  for (c in 1:NCELL){
    p <- as.numeric((sand_perc[c,3]) - (sand_perc[r,3]))
    psand[r,c] <- p
  }
}


### Save all the files you need to remake these matrices ###

saveRDS(pDist, "Spatial_Data/pDist.RDS")
saveRDS(preef, "Spatial_Data/preef.RDS")
saveRDS(pinvert, "Spatial_Data/pinvert.RDS")
saveRDS(pseagrass, "Spatial_Data/pseagrass.RDS")
saveRDS(pmacro, "Spatial_Data/pmacro.RDS")
saveRDS(psand, "Spatial_Data/psand.RDS")


#### Create movement probability using utility function ####

# Load files if need be
pDist <- readRDS("Spatial_Data/pDist.RDS")
preef <- readRDS("Spatial_Data/preef.RDS")
pinvert <- readRDS("Spatial_Data/pinvert.RDS")
pseagrass <- readRDS("Spatial_Data/pseagrass.RDS")
pmacro <- readRDS("Spatial_Data/pmacro.RDS")
psand <- readRDS("Spatial_Data/psand.RDS")

# First determine the utility of each of the sites 
# This is very sensitive to changes in the values particularly for reef 
# PROBABLY ALSO NEED TO PUT DEPTH IN HERE

## Small movement Swim Speed = 2.5 95% within approx 10km 
## Medium movement Swim Speed = 5 95% within approx 25km
## Big movement Swim Speed = 10 95% within approx 45km 

SwimSpeed <-  1

a = -(1/SwimSpeed)
b =  0.0150 #0.15     # Attractiveness of reef habitat
c =  0.0001 #0.001    # Attractiveness of invertebrate habitat
d =  0.0100 #0.1      # Attractiveness of seagrass habitat
e =  0.0010 #0.010    # Attractiveness of macroalgae habitat
f =  0.0001 #0.001    # Attractiveness of sand habitat

# From the attractivity of each habitat type, and the cell's distance to other cells, each cell's attractiveness is calculated (based on the % of each habitat in that cell)
Vj <- (a * pDist) + (b * preef) + (c * pinvert) + (d * pseagrass) + (e * pmacro) + (f * psand)

# Calculate the summed utility across the rows 
rowU <- matrix(NA, ncol=1, nrow=NCELL)
cellU <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    U <- exp(Vj[r,c])
    cellU[r,c] <- U
  }
}

rowU <- as.data.frame(rowSums(cellU))

# Calculate the probability that the fish will move to this site
Pj <- matrix(NA, ncol=NCELL, nrow=NCELL)
for (r in 1:NCELL){
  for (c in 1:NCELL){
    Pj[r,c] <- (exp(Vj[r,c]))/rowU[r,1]
  }
}
rowSums(Pj)


## Visualising the movement to double check it 
movement <- Pj[190, ] # Why are we looking at this row only ???

water_2 <- water %>%
  mutate(Move.Prob = movement)

map <- ggplot() +
  geom_sf(data=water_2, aes(fill=Move.Prob), color = NA, lwd=0)+
  #scale_fill_carto_c(palette="BluYl", direction=-1)+
  #annotate("text", x = 113.45, y = -21.5, colour = "black", size = 6, label=Years[YEAR])+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
map

data <- as.data.frame(Pj[190, ])
data <- data %>%
  mutate(Distance = pDist[190, ]) %>%
  arrange(Distance) %>%
  mutate(Cumulative_Prob = cumsum(`Pj[190, ]`))

plot(y=data$Cumulative_Prob, x=data$Distance, type="l")

#### Recruitment matrix ####
## Want the recruits to be in the lagoons and then move out from there 
dispersal <- seagrass_perc %>% 
  dplyr::select(perc_habitat, ID) %>%
  distinct(ID, .keep_all = TRUE)
dispersal$area <- as.vector(water$cell_area*0.000001)


dispersal <- dispersal %>% 
  filter(perc_habitat!=0) %>% 
  mutate(area = area*0.1,
         perc_habitat = perc_habitat*0.5)

recruitment <- array(0, dim=c(nrow(dispersal), 2))
recruitment[ ,2] <- as.numeric(dispersal$ID) 

cellU <- matrix(0, ncol=2, nrow=(nrow(dispersal)))
cellU[,2] <- as.numeric(dispersal$ID) 


for(cell in 1:nrow(recruitment)){
  U <- exp(dispersal[cell,1])
  cellU[cell,1] <- as.numeric(U)
}

rowU <- as.data.frame(sum(cellU[,1]))

for (cell in 1:nrow(recruitment)){
  recruitment[cell,1] <- cellU[cell, 1]/rowU[1,1]
}
colSums(recruitment)

recruitment <- as.data.frame(recruitment) %>% 
  rename(ID = "V2")

recruitment <- merge(recruitment, seagrass_perc, by="ID", all=T) %>% #check that cells with no lagoon habitat have 0 probability of recruitment
  mutate_all(~replace(., is.na(.), 0)) #For cells where there was no lagoon habitat put probability of recruitment as 0

recruitment <- as.vector(recruitment[,2])

#### Recruit movement ####
## Want the recruits to stay in the lagoon until they mature and move to the reef

SwimSpeed <-  1

a = -(1/SwimSpeed)
b =  0.0150 #0.15     # Attractiveness of reef habitat
c =  0.0001 #0.001    # Attractiveness of invertebrate habitat
d =  0.0100 #0.1      # Attractiveness of seagrass habitat
e =  0.0010 #0.010    # Attractiveness of macroalgae habitat
f =  0.0001 #0.001    # Attractiveness of sand habitat

Recj <- (a * pDist) + (b * preef) + (c * pinvert) + (d * pseagrass) + (e * pmacro) + (f * psand)
# Calculate the summed utility across the rows 
rowU <- matrix(NA, ncol=1, nrow=NCELL)
cellU <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    U <- exp(Recj[r,c])
    cellU[r,c] <- U
  }
}

rowU <- as.data.frame(rowSums(cellU))

# Calculate the probability that the fish will move to this site

ProbRec <- matrix(NA, ncol=NCELL, nrow=NCELL)

for (r in 1:NCELL){
  for (c in 1:NCELL){
    ProbRec[r,c] <- (exp(Recj[r,c]))/rowU[r,1]
  }
}
rowSums(ProbRec)

#### Save files for next step ####
saveRDS(Pj, "Staging/movement_really_slow")
saveRDS(ProbRec, "Staging/juvmove")
saveRDS(recruitment, "Staging/recruitment")
saveRDS(water, "Staging/water")



### 04_A commercial fishing effort, utility function --------------------------


## 4. Set up a utility function -----------------------------------------------

# (this is for having 2 different utilities, before and after the SZ. replaced by 80 utilities, one per year)

Cell_Vars <- DistBR %>% 
  mutate(Area = as.vector((water$cell_area)/1000000)) # Cells are now in km^2 but with no units

NCELL_pre18 <- NCELL # number of cells you can fish in (before SZ)
NCELL_post18 <- NCELL - nrow(water[water$status %in% c('NTZ_boat', 'NTZ_boat_shore'),]) # number of cells you can fish in (after SZ)

Vj <- Cell_Vars %>%
  mutate(vj = rowSums(across(all_of(BR$name)), na.rm = TRUE)) %>% 
  glimpse()

Vj_pre18 <- Vj # catchable cells before SZ
Vj_post18 <- Vj[-c(water$ID[water$status %in% c('NTZ_boat', 'NTZ_boat_shore')]), ] # catchable cells after SZ

area_col <- grep("Area", colnames(Vj)) # extract the column index for area - important for calculations. in theory it should be the same pre- and post-NTZ.


# pre-SZ

BR_U_pre18 <- as.data.frame(matrix(0, nrow = NCELL_pre18, ncol = length(unique(BR$name)))) # Set up data frame to hold utilities of cells
colnames(BR_U_pre18) <- c(BR$name)

cellU <- matrix(NA, ncol = length(unique(BR$name)), nrow = NCELL_pre18)

for(RAMP in 1:length(unique(BR$name))){
  for(cell in 1:NCELL_pre18){
    U <- exp(Vj_pre18[cell, RAMP] + log(Vj_pre18[cell, area_col]))
    cellU[cell, RAMP] <- U
  }
} # this loop goes over each cell for each ramp, and calculates how catchable each cell is based on how close it is to a boat ramp and how popular that boat ramp is.

rowU <- as.data.frame(colSums(cellU))

for(RAMP in 1:length(unique(BR$name))){
  for(cell in 1:NCELL_pre18){
    BR_U_pre18[cell,RAMP] <- (exp(Vj_pre18[cell,RAMP]+log(Vj_pre18[cell, area_col])))/rowU[RAMP, 1]
  }
} # this loop goes over each cell for each ramp, and calculates how catchable (in %) each cell is based on how close it is to a boat ramp and how popular that boat ramp is.
colSums(BR_U_pre18) # all adds up to 1, perfect.
head(BR_U_pre18)

# Plot check (notice the difference between popular/non-popular ramps)
BR_U_pre18 <- BR_U_pre18 %>% mutate(ID = row_number())
water_catch <- water %>% mutate(ID = row_number()) %>% left_join(BR_U_pre18, by = "ID")
water_catch_long <- water_catch %>% pivot_longer(cols = BR$name, names_to = "Ramp", values_to = "Catchability")
ggplot(water_catch_long) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  facet_wrap(~ Ramp, ncol = 4) + 
  labs(title = "(before 2018) catchability of each cell, by boat ramp distance and cell size", fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/04A_commercial_catchability_surface_before2018.png", plot = last_plot())


# post-SZ

BR_U_post18 <- as.data.frame(matrix(0, nrow = NCELL_post18, ncol = length(unique(BR$name)))) #Set up data frame to hold utilities of cells
colnames(BR_U_post18) <- c(BR$name)

cellU <- matrix(NA, ncol = length(unique(BR$name)), nrow = NCELL_post18)

for(RAMP in 1:length(unique(BR$name))){
  for(cell in 1:NCELL_post18){
    U <- exp(Vj_post18[cell, RAMP] + log(Vj_post18[cell, area_col]))
    cellU[cell, RAMP] <- U
  }
} # this loop goes over each cell for each ramp after the NTZ is in place, and calculates how catchable each cell is based on how close it is to a boat ramp and how popular that boat ramp is.

rowU <- as.data.frame(colSums(cellU))

for (RAMP in 1:length(unique(BR$name))){
  for (cell in 1:NCELL_post18){
    BR_U_post18[cell, RAMP] <- (exp(Vj_post18[cell, RAMP] + log(Vj_post18[cell, area_col])))/rowU[RAMP, 1]
  }
} # this loop goes over each cell for each ramp after the NTZ is in place, and calculates how catchable each cell is based on how close it is to a boat ramp and how popular that boat ramp is.
colSums(BR_U_post18) # adds up to 1, perfect.
head(BR_U_post18)

# Plot check (notice the difference between popular/non-popular ramps)
BR_U_post18 <- BR_U_post18 %>% mutate(ID = row_number())
water_catch <- water[water$status == "Fished",] %>% mutate(ID = row_number()) %>% left_join(BR_U_post18, by = "ID")
water_catch_long <- water_catch %>% pivot_longer(cols = BR$name, names_to = "Ramp", values_to = "Catchability")
ggplot(water_catch_long) + 
  geom_sf(aes(fill = log(Catchability)), color = NA) + 
  scale_fill_gradientn(colours = colour_palette[4:6]) +
  facet_wrap(~ Ramp, ncol = 4) + 
  labs(title = "(after 2018) catchability of each cell, by boat ramp distance and cell size", fill = "log(Catchability)")
ggsave("plots/checking_plots_during_setup/04A_catchability_surface_after2018.png", plot = last_plot())


