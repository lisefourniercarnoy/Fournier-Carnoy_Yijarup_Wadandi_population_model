# -----------------------------------------------------------------------------

# Project: Wadandi Pink Snapper Population Model
# Data:    ?
# Task:    Set up a population out of the population parameters and fishing effort
# Author:  Lise Fournier-Carnoy / adapted from Charlotte Aston
# Date:    September 2025

# -----------------------------------------------------------------------------

# Status: Starting out...

# -----------------------------------------------------------------------------

rm(list = ls()) # clear environment

## Load libraries -------------------------------------------------------------

library(tidyverse) # for data manipulation + ggplot
library(sf) # for dealing with shapefiles
#library(sfnetworks) # for dealing with networks
#library(raster) # for dealing with rasters
#library(forcats)
library(RColorBrewer)

library(MQMF)
library(Rcpp)
library(RcppArmadillo)

library(gmailr)
library(beepr)
library(pwr)



## Read in functions ----------------------------------------------------------

sourceCpp("functions/C_Model_RcppArm_test.cpp")
source("functions/X_Functions.R")


#### PRE-SETS ####

## Create colours for the plot
pop.groups <- c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150)
my.colours <- "PuBu"

model.name <- "wadandi"

names <- c("Current NTZs", "No management", "Wadandi managemenent")

n_yrs <- 2024-1945 +1

## LOAD FILES -----------------------------------------------------------------

# read files from the real-life simulation
adult_movement  <- readRDS("data/output_data/03_adult_movement_10_swim_speed.rds")
settlement      <- readRDS("data/output_data/03_recruitment.rds") 
effort          <- readRDS("data/output_data/04A_commercial_burn_in_fishing.rds")
no_take_list    <- readRDS("data/output_data/02_no_take_list.rds")
water           <- readRDS("data/output_data/03_water.rds") %>% st_make_valid()
burn_in_pop     <- readRDS("data/output_data/06_burn_in_population.rds")
selectivity     <- readRDS("data/output_data/05_selectivity_retention.rds")
mature          <- readRDS("data/output_data/05_maturity.rds")
weight          <- readRDS("data/output_data/05_weight.rds")

NCELL <- nrow(water)


## WHOLE POPULATION LMs -------------------------------------------------------

# get the simulations we want to test
total_pop_list <- list()
total_pop_list[[1]] <-  readRDS("simulations/dummy_run/07_A_SC_age_distribution_try1.rds")
total_pop_list[[2]] <-  readRDS("simulations/dummy_run/07_A_S00_age_distribution_try1.rds")
total_pop_list[[3]] <-  readRDS("simulations/dummy_run/07_A_S01_age_distribution_try1.rds")

# NEED TO TURN OFF THE BIT OF THE FUNCITON THAT CREATES THE MEDIANS AS YOU WANT ALL OF THE SEPARATE RUNS

# this function does ___
total_pop <- total.pop.format.full(pop.file.list = total_pop_list, 
                                   scenario.names = names, 
                                   nsim = dim(total_pop_list[[1]])[3], # total_pop_list is dimensioned as (max age of the fish) x (years modelled) x (simulations)
                                   nyears = n_yrs, 
                                   startyear = n_yrs,
                                   maxage = 30, # see script 05 for the correct age.
                                   mat = mature, 
                                   kg = weight)
ggplot(total_pop, aes(x = Year)) +
  geom_point(aes(y = V1, colour = Scenario), alpha = 0.5) +
  geom_point(aes(y = V2, colour = Scenario), alpha = 0.5) +
  geom_point(aes(y = V3, colour = Scenario), alpha = 0.5) +
  geom_point(aes(y = V4, colour = Scenario), alpha = 0.5) +
  geom_point(aes(y = V5, colour = Scenario), alpha = 0.5)
  


total_pop <- total_pop %>% 
  filter(Year %in% c(2018)) %>% 
  #filter(Scenario %in% c("Historical and Current NTZs", "Temporal Management Only")) %>% 
  mutate(Year = as.factor(Year),
         Scenario = as.factor(Scenario)) %>% 
  mutate(Movement = "Medium")

n_sim <- 5 # number of simulations in  script 07_A
n_scenario <- 3 # number of scenarios (here i have current, no management, and wadandi management)

differences <- total_pop %>% 
  mutate(ID = rep(1:n_sim, times = n_scenario)) %>%
  #mutate(MatBio = log(MatBio)) %>% 
  pivot_wider(id_cols=c(ID),names_from = Scenario, values_from = MatBio) %>% 
  mutate(Comp.1.2 = `Current NTZs` - `No management`,
         Comp.1.3 = `Current NTZs`- `Wadandi managemenent`,
         Comp.2.3 = `No management` - `Wadandi managemenent`
         ) %>% 
  mutate(Comp.1.2.Perc = ((`Current NTZs` - `No management`)/((`Current NTZs`+`No management`)/2))*100,
         Comp.1.3.Perc = ((`Current NTZs`- `Wadandi managemenent`)/((`Current NTZs`+`Wadandi managemenent`)/2))*100,
         Comp.2.3.Perc = ((`No management` - `Wadandi managemenent`)/((`No management`+`Wadandi managemenent`)/2))*100
         )

Summary.Perc.Change <- differences %>% 
  summarise_at(vars(Comp.1.2.Perc:Comp.2.3.Perc), median)

test_dat <- differences 


wilcox.test(test_dat$Comp.1.2, conf.int = TRUE, conf.level = 0.95)


# mod <- aov(log(MatBio) ~ Scenario, data = total_pop)
# summary(mod)
# 
# total_pop %>% group_by(Scenario) %>% 
#   summarise(mean = mean(MatBio),
#             sd = sd(MatBio))
# 
# ggplot()+
#   geom_boxplot(data=total_pop, aes(x=Scenario, y=log(MatBio)), notch=T)
# 
# 
# summary(mod)
# plot(mod$residuals)
# 
# mod.tukey <- TukeyHSD(mod)
# mod.tukey



#* ZONE MODELS ####
SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_medium_movement"))
SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_medium_movement"))

SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_medium_movement"))
SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_medium_movement"))

SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_medium_movement"))
SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_medium_movement"))

SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_medium_movement"))
SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_medium_movement"))


#* Format zone data ####

## S00 - Normal

NTZ.F.Ages.S00 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S00, f.list = SP_Pop_NTZ_S00, scenario.name = Names[1], nsim = 200)

NTZ.S00 <- NTZ.F.Ages.S00[[1]] %>% 
  mutate(Movement = "Medium")
F.S00 <- NTZ.F.Ages.S00[[2]]


## S01 - No NTZs (and no temporal closure)
NTZ.F.Ages.S01 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S01, f.list = SP_Pop_NTZ_S01, scenario.name = Names[2], nsim = 200)

NTZ.S01 <- NTZ.F.Ages.S01[[1]] %>% 
  mutate(Movement = "Medium")
F.S01 <- NTZ.F.Ages.S01[[2]]


## S02 - Temporal Closure, no NTZs
NTZ.F.Ages.S02 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S02, f.list = SP_Pop_NTZ_S02, scenario.name = Names[3], nsim = 200)

NTZ.S02 <- NTZ.F.Ages.S02[[1]] %>% 
  mutate(Movement = "Medium")
F.S02 <- NTZ.F.Ages.S02[[2]]


## S03 - Temporal Closure and NTZs
NTZ.F.Ages.S03 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S03, f.list = SP_Pop_NTZ_S03, scenario.name = Names[4], nsim = 200)

NTZ.S03 <- NTZ.F.Ages.S03[[1]] %>% 
  mutate(Movement = "Medium")
F.S03 <- NTZ.F.Ages.S03[[2]]


## Put everything together
Whole_Pop_Ages_NTZ <- rbind(NTZ.S00, NTZ.S01, NTZ.S02, NTZ.S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F.S00, F.S01, F.S02, F.S03) %>% 
  mutate(Zone = "F")

NTZ_data <- Whole_Pop_Ages_NTZ %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")

F_data <- Whole_Pop_Ages_F %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")


## NTZ models 
differences <- NTZ_data %>% 
  filter(Mod_Year %in% 2018) %>% 
  mutate(ID = rep(1:200, times=4)) %>%
  #mutate(MatBio = log(MatBio)) %>% 
  pivot_wider(id_cols=c(ID,Stage),names_from=Scenario, values_from=Number) %>% 
  mutate(Comp.1.2 = `Current NTZs` - `No temporal management or NTZs`,
         Comp.1.3 = `Current NTZs`- `Temporal management and NTZs`,
         Comp.1.4 = `Current NTZs` - `Temporal management`,
         Comp.2.3 = `No temporal management or NTZs` - `Temporal management and NTZs`,
         Comp.2.4 = `No temporal management or NTZs` - `Temporal management`,
         Comp.3.4 = `Temporal management and NTZs` - `Temporal management`) %>% 
  mutate(Comp.1.2.Perc = ((`Current NTZs` - `No temporal management or NTZs`)/((`Current NTZs`+`No temporal management or NTZs`)/2))*100,
         Comp.1.3.Perc = ((`Current NTZs`- `Temporal management and NTZs`)/((`Current NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.1.4.Perc = ((`Current NTZs` - `Temporal management`)/((`Current NTZs`+`Temporal management`)/2))*100,
         Comp.2.3.Perc = ((`No temporal management or NTZs` - `Temporal management and NTZs`)/((`No temporal management or NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.2.4.Perc = ((`No temporal management or NTZs` - `Temporal management`)/((`No temporal management or NTZs`+`Temporal management`)/2))*100,
         Comp.3.4.Perc = ((`Temporal management and NTZs` - `Temporal management`)/((`Temporal management and NTZs`+`Temporal management`)/2)*100))

Summary.Perc.Change <- differences %>% 
  group_by(Stage) %>% 
  summarise_at(vars(Comp.1.2.Perc:Comp.3.4.Perc), median)

test_dat <- differences %>% 
  filter(Stage %in% "Large Legal")


wilcox.test(test_dat$Comp.3.4, conf.int = TRUE, conf.level = 0.95)


# Recruits 
NTZ_Recruits <- NTZ_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.NTZ <- aov(log(Number) ~ Scenario, dat=NTZ_Recruits)
summary(mod1.NTZ)
plot(mod1.NTZ$residuals)

mod1.NTZ.tukey <- TukeyHSD(mod1.NTZ)
mod1.NTZ.tukey

F_Recruits <- F_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.F <- aov(log(Number) ~ Scenario, dat=F_Recruits)
summary(mod1.F)
plot(mod1.F$residuals)

mod1.F.tukey <- TukeyHSD(mod1.F)
mod1.F.tukey

# Sublegal
NTZ_Sublegal <- NTZ_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.NTZ <- aov(log(Number) ~ Scenario, dat=NTZ_Sublegal)
summary(mod2.NTZ)
plot(mod2.NTZ$residuals)

mod2.NTZ.tukey <- TukeyHSD(mod2.NTZ)
mod2.NTZ.tukey


F_Sublegal <- F_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.F <- lm(log(Number) ~ Scenario, dat=F_Sublegal)
summary(mod2.F)
plot(mod2.F$residuals)

# Legal
NTZ_Legal <- NTZ_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.NTZ <- aov(log(Number) ~ Scenario, dat=NTZ_Legal)
summary(mod3.NTZ)
plot(mod3.NTZ$residuals)

mod3.NTZ.tukey <- TukeyHSD(mod3.NTZ)
mod3.NTZ.tukey


F_Legal <- F_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.F <- lm(log(Number) ~ Scenario, dat=F_Legal)
summary(mod3.F)
plot(mod3.F$residuals)

mod3.F.aov <- aov(mod3.F)
mod3.F.tukey <- TukeyHSD(mod3.F.aov)
mod3.F.tukey



# Large Legal
NTZ_Large <- NTZ_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario)) 

mod4.NTZ <- aov(log(Number) ~ Scenario, dat=NTZ_Large)
summary(mod4.NTZ)
plot(mod4.NTZ$residuals)

mod4.NTZ.tukey <- TukeyHSD(mod4.NTZ)
mod4.NTZ.tukey


F_Large <- F_data %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod4.F <- lm(log(Number) ~ Scenario, dat=F_Large)
summary(mod4.F)
plot(mod4.F$residuals)

mod4.F.aov <- aov(mod4.F)
mod4.F.tukey <- TukeyHSD(mod4.F.aov)
mod4.F.tukey


#### DISTANCE ABUNDANCE AND CATCH GLMS ####
#* Work out which cells are within 100km of each of the boat ramps ####
setwd(sp_dir)
BR <- st_read("Boat_Ramps.shp") %>% 
  st_transform(4283)%>%
  st_make_valid() 
BR <- BR[1:4,]

network <- st_read(paste0(model.name, sep="_","network.shapefile.shp"))

WHA <- st_read("2013_02_WorldHeritageMarineProgramme.shp") %>% 
  st_transform(4283)%>%
  st_make_valid %>% 
  st_crop(xmin=112.5, xmax=114.7, ymin=-24, ymax=-20.5) 

setwd(sg_dir)

water <- readRDS(paste0(model.name, sep="_","water"))
BR <- st_as_sf(BR)
st_crs(BR) <- NA 

BR <- BR %>% 
  mutate(name = c("Bundegi","Exmouth","Tantabiddi","CoralBay"))

centroids <- st_centroid_within_poly(water)
points <- as.data.frame(st_coordinates(centroids))%>% #The points start at the bottom left and then work their way their way right
  mutate(ID=row_number()) 
points_sf <- st_as_sf(points, coords = c("X", "Y")) 
st_crs(points_sf) <- NA


model_WHA <- water %>% 
  st_intersects(., WHA) %>% 
  as.data.frame()

network <- as_sfnetwork(network, directed = FALSE) %>%
  activate("edges") %>%
  mutate(weight = edge_length())

net <- activate(network, "nodes")
network_matrix <- st_network_cost(net, from=BR, to=points_sf)
network_matrix <- network_matrix*111
dim(network_matrix)

DistBR <- as.data.frame(t(network_matrix)) %>% 
  rename("Bd_BR"=V1) %>% 
  rename("ExM_BR" = V2) %>% 
  rename("Tb_BR" = V3) %>% 
  rename("CrB_BR"=V4) %>% 
  mutate(CellID = row_number()) 

DistBR.WHA <- DistBR[c(as.numeric(model_WHA$row.id)), ]

NCELL.WHA <- nrow(DistBR.WHA)
BR.10km <- NULL
temp <- array(0, dim=c(NCELL.WHA,2)) %>% 
  as.data.frame() %>% 
  rename(.,"Distance" = V1,
         "CellID" = V2)

for(RAMP in 1:4){
  temp[,1:2] <- DistBR.WHA[,c(RAMP, 5)]
  
  Dist10 <- temp %>% 
    filter(Distance <=10)
  
  BR.10km <- rbind(BR.10km, Dist10)
}

BR.10km <- unique(BR.10km$CellID)

BR.50km <- NULL
temp <- array(0, dim=c(NCELL.WHA,2)) %>% 
  as.data.frame() %>% 
  rename(.,"Distance" = V1,
         "CellID" = V2)

for(RAMP in 1:4){
  temp[,1:2] <- DistBR.WHA[,c(RAMP, 5)]
  
  Dist50 <- temp %>% 
    filter(Distance>10 & Distance <=50)
  
  BR.50km <- rbind(BR.50km, Dist50)
}

BR.50km <- unique(BR.50km$CellID)

BR.100km <- NULL
temp <- array(0, dim=c(NCELL.WHA,2)) %>% 
  as.data.frame() %>% 
  rename(.,"Distance" = V1,
         "CellID" = V2)

for(RAMP in 1:4){
  temp[,1:2] <- DistBR.WHA[,c(RAMP, 5)]
  
  Dist100 <- temp %>% 
    filter(Distance >50 & Distance <=100)
  
  BR.100km <- rbind(BR.100km, Dist100)
}

BR.100km <- unique(BR.100km$CellID)

Distances <- list()

Distances[[1]] <- BR.10km
Distances[[2]] <- BR.50km
Distances[[3]] <- BR.100km

Dist.Names <- as.character(c("0-10 km", "10-50 km", "50-100 km"))

#* Abundance ####
setwd(pop_dir)

# Each layer is a simulation, rows are cells and columns are years
Pop.Dist <- list()

Pop.Dist[[1]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S00", sep="_", "medium_movement"))
Pop.Dist[[2]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S01", sep="_", "medium_movement"))
Pop.Dist[[3]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S02", sep="_", "medium_movement"))
Pop.Dist[[4]] <- readRDS(paste0(model.name, sep="_", "Cell_Population", sep="_", "S03", sep="_", "medium_movement"))

Names <- c("Historical and Current NTZs", "Neither NTZs nor Temporal Management", 
           "Temporal and Spatial Management","Temporal Management Only")

Dists.res <- dist.abundance.full(Pops = Pop.Dist, n.year=59, n.sim=200, n.cell=NCELL,
                                 n.scenario=4, fished.cells=water, distances=Distances,
                                 dist.names=Dist.Names, scen.names=Names, mod.years = seq(1960,2018,1), target.year = 59)


S00.Abundance <- Dists.res[[1]]
S01.Abundance <- Dists.res[[2]]
S02.Abundance <- Dists.res[[3]]
S03.Abundance <- Dists.res[[4]]

# Function is doing something really weird so I've done it manually by running the inside of the function
# 
# S00.Abundance <- Pop.Abundance[[1]]
# S01.Abundance <- Pop.Abundance[[2]]
# S02.Abundance <- Pop.Abundance[[3]]
# S03.Abundance <- Pop.Abundance[[4]]

Dist.Abundance <- rbind(S00.Abundance, S01.Abundance, S02.Abundance, S03.Abundance)

## Abundance model 10km
Dist.Abundance.10km <- Dist.Abundance %>% 
  filter(Distance %in% "0-10 km") %>% 
  mutate(Scenario = as.factor(Scenario)) 

mod1.10km <- lm(log(Total) ~ Scenario, dat=Dist.Abundance.10km)
summary(mod1.10km)
plot(mod1.10km$residuals)

## Abundance model 50km
Dist.Abundance.50km <- Dist.Abundance %>% 
  filter(Distance %in% "10-50 km") %>% 
  mutate(Scenario = as.factor(Scenario)) 

mod1.50km <- lm(log(Total) ~ Scenario, dat=Dist.Abundance.50km)
summary(mod1.50km)
plot(mod1.10km$residuals)

## Abundance model 100km
Dist.Abundance.100km <- Dist.Abundance %>% 
  filter(Distance %in% "50-100 km") %>% 
  mutate(Scenario = as.factor(Scenario)) 

mod1.100km <- lm(log(Total) ~ Scenario, dat=Dist.Abundance.100km)
summary(mod1.100km)
plot(mod1.10km$residuals)

#* Catch ####
Effort_Scen <- list()
Spatial_Qs <- list()

setwd(sg_dir)
Effort_Scen[[1]] <- readRDS(paste0(model.name, sep="_", "fishing"))
Spatial_Qs[[1]] <- readRDS(paste0(model.name, sep="_", "Spatial_q_NTZ"))
Spatial_Qs[[2]] <- readRDS(paste0(model.name, sep="_", "Spatial_q_No_NTZ"))

setwd(sim_dir)
Effort_Scen[[2]] <- readRDS(paste0(model.name, sep="_", "S01_fishing"))
Effort_Scen[[3]] <- readRDS(paste0(model.name, sep="_", "S02_fishing"))
Effort_Scen[[4]] <- readRDS(paste0(model.name, sep="_", "S03_fishing"))

## Convert back to effort in boat days 

Boat_Days <- list()

Boat_Days_Scen <- array(0, dim=c(NCELL, 12, 59))

for(S in 1:4){
  
  if(S==1|S==3){
    for (YEAR in 1:59){
      Boat_Days_Scen[,,YEAR] <- Effort_Scen[[S]][,,YEAR] / Spatial_Qs[[1]][,YEAR]
    }
  } else {
    for (YEAR in 1:59){
      Boat_Days_Scen[,,YEAR] <- Effort_Scen[[S]][,,YEAR] / Spatial_Qs[[2]][,YEAR]
    }
  }
  Boat_Days_Scen[is.nan(Boat_Days_Scen)] <- 0
  
  Boat_Days[[S]] <- Boat_Days_Scen
}

#* Sum up effort for every year in each of the scenarios
Effort_Dist <-list()

for(D in 1:3){
  
  Boat_Days_sum <- NULL
  
  for(S in 1:4){
    temp <- Boat_Days[[S]]
    temp2 <- temp[as.numeric(Distances[[D]]), ,] %>% 
      colSums(., dim=2) %>% 
      as.data.frame() %>% 
      mutate(Scenario = Names[S]) %>% 
      mutate(Distance = Dist.Names[D]) %>% 
      mutate(Year = seq(1960,2018,1))
    
    Boat_Days_sum <- rbind(Boat_Days_sum, temp2)
  }
  
  Boat_Days_sum <- Boat_Days_sum %>% 
    rename(Effort = ".")
  Effort_Dist[[D]] <- Boat_Days_sum
} 



setwd(pop_dir)

Pop.Catch <- list()

# Each layer is a simulation, rows are cells and columns are years
Pop.Catch[[1]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S00", sep="_", "medium_movement"))
Pop.Catch[[2]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S01", sep="_", "medium_movement"))
Pop.Catch[[3]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S02", sep="_", "medium_movement"))
Pop.Catch[[4]] <- readRDS(paste0(model.name, sep="_", "Catch_by_Cell", sep="_", "S03", sep="_", "medium_movement"))

## NEED TO TURN INTO CPUE
Catch.res <- dist.catch.full(Pops = Pop.Catch, n.year=59, n.sim=200, n.cell=NCELL,
                             n.scenario=4, fished.cells=water, distances=Distances,
                             dist.names=Dist.Names, scen.names=Names, mod.years = seq(1960,2018,1), target.year=59)


Pop.Catch.S00 <- Catch.res[[1]]
Pop.Catch.S01 <- Catch.res[[2]]
Pop.Catch.S02 <- Catch.res[[3]]
Pop.Catch.S03 <- Catch.res[[4]]

## Function is running weirdly so I've run the inside of the function manually as I can't be bothered to fix it right now

# Pop.Catch.S00 <- Pop.Catch[[1]]
# Pop.Catch.S01 <- Pop.Catch[[2]]
# Pop.Catch.S02 <- Pop.Catch[[3]]
# Pop.Catch.S03 <- Pop.Catch[[4]]

Dist.Catch <- rbind(Pop.Catch.S00, Pop.Catch.S01, Pop.Catch.S02, Pop.Catch.S03)

differences <- Dist.Catch %>% 
  filter(Year %in% 2018) %>% 
  mutate(ID = rep(1:200, times=4*3)) %>%
  #mutate(MatBio = log(MatBio)) %>% 
  pivot_wider(id_cols=c(ID, Distance), names_from=Scenario, values_from=Total) %>% 
  mutate(Comp.1.2 = `Current NTZs` - `No temporal management or NTZs`,
         Comp.1.3 = `Current NTZs`- `Temporal management and NTZs`,
         Comp.1.4 = `Current NTZs` - `Temporal management`,
         Comp.2.3 = `No temporal management or NTZs` - `Temporal management and NTZs`,
         Comp.2.4 = `No temporal management or NTZs` - `Temporal management`,
         Comp.3.4 = `Temporal management and NTZs` - `Temporal management`) %>% 
  mutate(Comp.1.2.Perc = ((`Current NTZs` - `No temporal management or NTZs`)/((`Current NTZs`+`No temporal management or NTZs`)/2))*100,
         Comp.1.3.Perc = ((`Current NTZs`- `Temporal management and NTZs`)/((`Current NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.1.4.Perc = ((`Current NTZs` - `Temporal management`)/((`Current NTZs`+`Temporal management`)/2))*100,
         Comp.2.3.Perc = ((`No temporal management or NTZs` - `Temporal management and NTZs`)/((`No temporal management or NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.2.4.Perc = ((`No temporal management or NTZs` - `Temporal management`)/((`No temporal management or NTZs`+`Temporal management`)/2))*100,
         Comp.3.4.Perc = ((`Temporal management and NTZs` - `Temporal management`)/((`Temporal management and NTZs`+`Temporal management`)/2)*100))

# Could work out the proportion of one out of the other and see if that is different from one 

Summary.Perc.Change <- differences %>% 
  group_by(Distance) %>% 
  summarise_at(vars(Comp.1.2.Perc:Comp.3.4.Perc), median)

test_dat <- differences %>% 
  filter(Distance %in% "50-100 km")

wilcox.test(test_dat$Comp.3.4, conf.int = TRUE, conf.level = 0.95)

## Catch model 10km
Effort.10km <- Effort_Dist[[1]]

Dist.Catch.10km <- Dist.Catch %>% 
  filter(Year %in% 2018) %>% 
  filter(Distance %in% "0-10 km") %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  left_join(., Effort.10km, by=c("Scenario", "Distance", "Year")) %>% 
  mutate(CPUE = Total/Effort)

mod1.10km <- aov(log(CPUE) ~ Scenario, dat=Dist.Catch.10km)
summary(mod1.10km)
plot(mod1.10km$residuals)

mod1.10km.Tukey <- TukeyHSD(mod1.10km)
mod1.10km.Tukey 

## Catch model 50km
Effort.50km <- Effort_Dist[[2]]

Dist.Catch.50km <- Dist.Catch %>% 
  filter(Year %in% 2018) %>% 
  filter(Distance %in% "10-50 km") %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  left_join(., Effort.50km, by=c("Scenario", "Distance", "Year")) %>% 
  mutate(CPUE = Total/Effort)

mod2.50km <- aov(log(CPUE) ~ Scenario, dat=Dist.Catch.50km)
summary(mod2.50km)
plot(mod2.50km$residuals)

mod2.50km.Tukey <- TukeyHSD(mod2.50km)
mod2.50km.Tukey 


## Catch model 100km
Effort.100km <- Effort_Dist[[3]]

Dist.Catch.100km <- Dist.Catch %>% 
  filter(Year %in% 2018) %>% 
  filter(Distance %in% "50-100 km") %>% 
  mutate(Scenario = as.factor(Scenario)) %>% 
  left_join(., Effort.100km, by=c("Scenario", "Distance", "Year")) %>% 
  mutate(CPUE = Total/Effort)

mod3.100km <- aov(log(CPUE) ~ Scenario, dat=Dist.Catch.100km)
summary(mod3.100km)
plot(mod3.100km$residuals)

mod3.100km.Tukey <- TukeyHSD(mod3.100km)
mod3.100km.Tukey 


#### MOVEMENT SCENARIO GLMS ####
setwd(pop_dir)

#* Slow ####
total_pop_list_slow <- list()

total_pop_list_slow[[1]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S00_slow_movement"))
total_pop_list_slow[[2]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S01_slow_movement"))
total_pop_list_slow[[3]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S02_slow_movement"))
total_pop_list_slow[[4]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S03_slow_movement"))

# NEED TO TURN OFF THE BIT OF THE FUNCITON THAT CREATES THE MEDIANS AS YOU WANT ALL OF THE SEPARATE RUNS
total_pop_slow <- total.pop.format.full(pop.file.list = total_pop_list, scenario.names = Names, nsim=100, nyears=59, startyear=26, maxage=30, mat = Mature, kg=Weight)

differences <- total_pop_slow %>% 
  filter(Year %in% 2018) %>% 
  mutate(ID = rep(1:100, times=4)) %>%
  #mutate(MatBio = log(MatBio)) %>% 
  pivot_wider(id_cols=c(ID), names_from=Scenario, values_from=MatBio) %>% 
  mutate(Comp.1.2 = `Current NTZs` - `No temporal management or NTZs`,
         Comp.1.3 = `Current NTZs`- `Temporal management and NTZs`,
         Comp.1.4 = `Current NTZs` - `Temporal management`,
         Comp.2.3 = `No temporal management or NTZs` - `Temporal management and NTZs`,
         Comp.2.4 = `No temporal management or NTZs` - `Temporal management`,
         Comp.3.4 = `Temporal management and NTZs` - `Temporal management`) %>% 
  mutate(Comp.1.2.Perc = ((`Current NTZs` - `No temporal management or NTZs`)/((`Current NTZs`+`No temporal management or NTZs`)/2))*100,
         Comp.1.3.Perc = ((`Current NTZs`- `Temporal management and NTZs`)/((`Current NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.1.4.Perc = ((`Current NTZs` - `Temporal management`)/((`Current NTZs`+`Temporal management`)/2))*100,
         Comp.2.3.Perc = ((`No temporal management or NTZs` - `Temporal management and NTZs`)/((`No temporal management or NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.2.4.Perc = ((`No temporal management or NTZs` - `Temporal management`)/((`No temporal management or NTZs`+`Temporal management`)/2))*100,
         Comp.3.4.Perc = ((`Temporal management and NTZs` - `Temporal management`)/((`Temporal management and NTZs`+`Temporal management`)/2)*100))

Summary.Perc.Change <- differences %>% 
  #group_by(Distance) %>% 
  summarise_at(vars(Comp.1.2.Perc:Comp.3.4.Perc), median)

test_dat <- differences 


wilcox.test(test_dat$Comp.1.4, conf.int = TRUE, conf.level = 0.95)


total_pop_slow <- total_pop_slow %>% 
  filter(Year %in% c(2018)) %>% 
  mutate(Year = as.factor(Year),
         Scenario = as.factor(Scenario)) %>% 
  mutate(Movement = "Slow")

mod.slow <- lm(log(MatBio) ~ Scenario, data=total_pop_slow)
summary(mod.slow)
plot(mod.slow$residuals)

mod.slow.aov <- aov(mod.slow)
tukey.slow <- TukeyHSD(mod.slow.aov)
tukey.slow 

#* Fast ####
total_pop_list_fast <- list()

total_pop_list_fast[[1]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S00_fast_movement"))
total_pop_list_fast[[2]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S01_fast_movement"))
total_pop_list_fast[[3]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S02_fast_movement"))
total_pop_list_fast[[4]] <-  readRDS(paste0(model.name, sep="_","Age_Distribution_S03_fast_movement"))

total_pop_fast <- total.pop.format.full(pop.file.list = total_pop_list_fast, scenario.names = Names, nsim=100, nyears=59, startyear=26, maxage=30, mat = Mature, kg=Weight)

differences <- total_pop_fast %>% 
  filter(Year %in% 2018) %>% 
  mutate(ID = rep(1:100, times=4)) %>%
  #mutate(MatBio = log(MatBio)) %>% 
  pivot_wider(id_cols=c(ID), names_from=Scenario, values_from=MatBio) %>% 
  mutate(Comp.1.2 = `Current NTZs` - `No temporal management or NTZs`,
         Comp.1.3 = `Current NTZs`- `Temporal management and NTZs`,
         Comp.1.4 = `Current NTZs` - `Temporal management`,
         Comp.2.3 = `No temporal management or NTZs` - `Temporal management and NTZs`,
         Comp.2.4 = `No temporal management or NTZs` - `Temporal management`,
         Comp.3.4 = `Temporal management and NTZs` - `Temporal management`) %>% 
  mutate(Comp.1.2.Perc = ((`Current NTZs` - `No temporal management or NTZs`)/((`Current NTZs`+`No temporal management or NTZs`)/2))*100,
         Comp.1.3.Perc = ((`Current NTZs`- `Temporal management and NTZs`)/((`Current NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.1.4.Perc = ((`Current NTZs` - `Temporal management`)/((`Current NTZs`+`Temporal management`)/2))*100,
         Comp.2.3.Perc = ((`No temporal management or NTZs` - `Temporal management and NTZs`)/((`No temporal management or NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.2.4.Perc = ((`No temporal management or NTZs` - `Temporal management`)/((`No temporal management or NTZs`+`Temporal management`)/2))*100,
         Comp.3.4.Perc = ((`Temporal management and NTZs` - `Temporal management`)/((`Temporal management and NTZs`+`Temporal management`)/2)*100))

Summary.Perc.Change <- differences %>% 
  #group_by(Distance) %>% 
  summarise_at(vars(Comp.1.2.Perc:Comp.3.4.Perc), median)

test_dat <- differences 


wilcox.test(test_dat$Comp.1.4, conf.int = TRUE, conf.level = 0.95)

total_pop_fast <- total_pop_fast %>% 
  filter(Year %in% c(2018)) %>% 
  mutate(Year = as.factor(Year),
         Scenario = as.factor(Scenario)) %>% 
  mutate(Movement = "Fast")

mod.fast <- lm(log(MatBio) ~ Scenario, data=total_pop_fast)
summary(mod.fast)
plot(mod.fast$residuals)

## Comparing same scenario but different movement
total_pop_movement <- rbind(total_pop, total_pop_slow, total_pop_fast)

S00_movement <- total_pop_movement %>% 
  filter(Scenario %in% "Historical and Current NTZs")

mod.S00 <- lm(log(MatBio) ~ Movement, data=S00_movement)

mod.S00.aov <- aov(mod.S00)
tukey.S00 <- TukeyHSD(mod.S00.aov)
tukey.S00 

S01_movement <- total_pop_movement %>% 
  filter(Scenario %in% "Neither NTZs nor Temporal Management")

mod.S01 <- lm(log(MatBio) ~ Movement, data=S01_movement)

mod.S01.aov <- aov(mod.S01)
tukey.S01 <- TukeyHSD(mod.S01.aov)
tukey.S01 

S02_movement <- total_pop_movement %>% 
  filter(Scenario %in% "Temporal Management Only")

mod.S02 <- lm(log(MatBio) ~ Movement, data=S02_movement)

mod.S02.aov <- aov(mod.S02)
tukey.S02 <- TukeyHSD(mod.S02.aov)
tukey.S02 

S03_movement <- total_pop_movement %>% 
  filter(Scenario %in% "Temporal and Spatial Management")

mod.S03 <- lm(log(MatBio) ~ Movement, data=S03_movement)

mod.S03.aov <- aov(mod.S03)
tukey.S03 <- TukeyHSD(mod.S03.aov)
tukey.S03 


# #* ZONE LINEAR MODELS SLOW 
# setwd(pop_dir)
# SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_slow_movement"))
# SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_fast_movement"))
# 
# SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_slow_movement"))
# SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_fast_movement"))
# 
# SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_slow_movement"))
# SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_slow_movement"))
# 
# SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_slow_movement"))
# SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_slow_movement"))
# 


#* ZONE LINEAR MODELS - SLOW ####
setwd(pop_dir)
SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_slow_movement"))
SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_slow_movement"))

SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_slow_movement"))
SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_slow_movement"))

SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_slow_movement"))
SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_slow_movement"))

SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_slow_movement"))
SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_slow_movement"))


#* Format zone data ####

## HAVE TURNED OFF THE PART OF THE FUNCTION THAT CREATES THE MEDIANS ##
## S00 - Normal

NTZ.F.Ages.S00 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S00, f.list = SP_Pop_NTZ_S00, scenario.name = Names[1], nsim = 100)

NTZ.S00 <- NTZ.F.Ages.S00[[1]] %>% 
  mutate(Movement = "Slow")
F.S00 <- NTZ.F.Ages.S00[[2]]


## S01 - No NTZs (and no temporal closure)
NTZ.F.Ages.S01 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S01, f.list = SP_Pop_NTZ_S01, scenario.name = Names[2], nsim = 100)

NTZ.S01 <- NTZ.F.Ages.S01[[1]] %>% 
  mutate(Movement = "Slow")
F.S01 <- NTZ.F.Ages.S01[[2]]


## S02 - Temporal Closure, no NTZs
NTZ.F.Ages.S02 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S02, f.list = SP_Pop_NTZ_S02, scenario.name = Names[3], nsim = 100)

NTZ.S02 <- NTZ.F.Ages.S02[[1]] %>% 
  mutate(Movement = "Slow")
F.S02 <- NTZ.F.Ages.S02[[2]]


## S03 - Temporal Closure and NTZs
NTZ.F.Ages.S03 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S03, f.list = SP_Pop_NTZ_S03, scenario.name = Names[4], nsim = 100)

NTZ.S03 <- NTZ.F.Ages.S03[[1]] %>% 
  mutate(Movement = "Slow")
F.S03 <- NTZ.F.Ages.S03[[2]]


## Put everything together
Whole_Pop_Ages_NTZ <- rbind(NTZ.S00, NTZ.S01, NTZ.S02, NTZ.S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F.S00, F.S01, F.S02, F.S03) %>% 
  mutate(Zone = "F")

NTZ_data_slow <- Whole_Pop_Ages_NTZ %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")

F_data_slow <- Whole_Pop_Ages_F %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")


differences <- NTZ_data_slow %>% 
  filter(Mod_Year %in% 2018) %>% 
  mutate(ID = rep(1:100, times=4)) %>%
  #mutate(MatBio = log(MatBio)) %>% 
  pivot_wider(id_cols=c(ID, Stage), names_from=Scenario, values_from=Number) %>% 
  mutate(Comp.1.2 = `Current NTZs` - `No temporal management or NTZs`,
         Comp.1.3 = `Current NTZs`- `Temporal management and NTZs`,
         Comp.1.4 = `Current NTZs` - `Temporal management`,
         Comp.2.3 = `No temporal management or NTZs` - `Temporal management and NTZs`,
         Comp.2.4 = `No temporal management or NTZs` - `Temporal management`,
         Comp.3.4 = `Temporal management and NTZs` - `Temporal management`) %>% 
  mutate(Comp.1.2.Perc = ((`Current NTZs` - `No temporal management or NTZs`)/((`Current NTZs`+`No temporal management or NTZs`)/2))*100,
         Comp.1.3.Perc = ((`Current NTZs`- `Temporal management and NTZs`)/((`Current NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.1.4.Perc = ((`Current NTZs` - `Temporal management`)/((`Current NTZs`+`Temporal management`)/2))*100,
         Comp.2.3.Perc = ((`No temporal management or NTZs` - `Temporal management and NTZs`)/((`No temporal management or NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.2.4.Perc = ((`No temporal management or NTZs` - `Temporal management`)/((`No temporal management or NTZs`+`Temporal management`)/2))*100,
         Comp.3.4.Perc = ((`Temporal management and NTZs` - `Temporal management`)/((`Temporal management and NTZs`+`Temporal management`)/2)*100))

Summary.Perc.Change <- differences %>% 
  ungroup() %>% 
  #group_by(Stage) %>% 
  summarise_at(vars(Comp.1.2.Perc:Comp.3.4.Perc), median)

## NTZ models 
# Recruits 
NTZ_Recruits_slow <- NTZ_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.NTZ.slow <- lm(log(Number) ~ Scenario, dat=NTZ_Recruits_slow)
summary(mod1.NTZ.slow)
plot(mod1.NTZ.slow$residuals)

mod1.NTZ.slow.aov <- aov(mod1.NTZ.slow)
mod1.NTZ.slow.tukey <- TukeyHSD(mod1.NTZ.slow.aov)
mod1.NTZ.slow.tukey

F_Recruits_slow <- F_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.F.slow <- lm(log(Number) ~ Scenario, dat=F_Recruits_slow)
summary(mod1.F.slow)
plot(mod1.F.slow$residuals)

mod1.F.slow.aov <- aov(mod1.F.slow)
mod1.F.slow.tukey <- TukeyHSD(mod1.F.slow.aov)
mod1.F.slow.tukey

# Sublegal
NTZ_Sublegal_slow <- NTZ_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.NTZ.slow <- lm(log(Number) ~ Scenario, dat=NTZ_Sublegal_slow)
summary(mod2.NTZ.slow)
plot(mod2.NTZ.slow$residuals)

mod2.NTZ.slow.aov <- aov(mod2.NTZ.slow)
mod2.NTZ.slow.tukey <- TukeyHSD(mod2.NTZ.slow.aov)
mod2.NTZ.slow.tukey


F_Sublegal_slow <- F_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.F.slow <- lm(log(Number) ~ Scenario, dat=F_Sublegal_slow)
summary(mod2.F.slow)
plot(mod2.F.slow$residuals)

mod2.F.slow.aov <- aov(mod2.F.slow)
mod2.F.slow.tukey <- TukeyHSD(mod2.F.slow.aov)
mod2.F.slow.tukey

# Legal
NTZ_Legal_slow <- NTZ_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.NTZ.slow <- lm(log(Number) ~ Scenario, dat=NTZ_Legal_slow)
summary(mod3.NTZ.slow)
plot(mod3.NTZ.slow$residuals)

mod3.NTZ.slow.aov <- aov(mod3.NTZ.slow)
mod3.NTZ.slow.tukey <- TukeyHSD(mod3.NTZ.slow.aov)
mod3.NTZ.slow.tukey


F_Legal_slow <- F_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.F.slow <- lm(log(Number) ~ Scenario, dat=F_Legal_slow)
summary(mod3.F.slow)
plot(mod3.F.slow$residuals)

mod3.F.slow.aov <- aov(mod3.F.slow)
mod3.F.slow.tukey <- TukeyHSD(mod3.F.slow.aov)
mod3.F.slow.tukey



# Large Legal
NTZ_Large_slow <- NTZ_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario)) 

mod4.NTZ.slow <- lm(log(Number) ~ Scenario, dat=NTZ_Large_slow)
summary(mod4.NTZ.slow)
plot(mod4.NTZ.slow$residuals)

mod4.NTZ.slow.aov <- aov(mod4.NTZ.slow)
mod4.NTZ.slow.tukey <- TukeyHSD(mod4.NTZ.slow.aov)
mod4.NTZ.slow.tukey


F_Large_slow <- F_data_slow %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod4.F.slow <- lm(log(Number) ~ Scenario, dat=F_Large_slow)
summary(mod4.F.slow)
plot(mod4.F.slow$residuals)

mod4.F.slow.aov <- aov(mod4.F.slow)
mod4.F.slow.tukey <- TukeyHSD(mod4.F.slow.aov)
mod4.F.slow.tukey

#* ZONE LINEAR MODELS - FAST ####
setwd(pop_dir)
SP_Pop_NTZ_S00 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S00_fast_movement"))
SP_Pop_F_S00 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S00_fast_movement"))

SP_Pop_NTZ_S01 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S01_fast_movement"))
SP_Pop_F_S01 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S01_fast_movement"))

SP_Pop_NTZ_S02 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S02_fast_movement"))
SP_Pop_F_S02 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S02_fast_movement"))

SP_Pop_NTZ_S03 <- readRDS(paste0(model.name, sep="_", "Sp_Population_NTZ_S03_fast_movement"))
SP_Pop_F_S03 <- readRDS(paste0(model.name, sep="_","Sp_Population_F_S03_fast_movement"))


#* Format zone data ####

## S00 - Normal

NTZ.F.Ages.S00 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S00, f.list = SP_Pop_NTZ_S00, scenario.name = Names[1], nsim = 100)

NTZ.S00 <- NTZ.F.Ages.S00[[1]] %>% 
  mutate(Movement = "Fast")
F.S00 <- NTZ.F.Ages.S00[[2]]


## S01 - No NTZs (and no temporal closure)
NTZ.F.Ages.S01 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S01, f.list = SP_Pop_NTZ_S01, scenario.name = Names[2], nsim = 100)

NTZ.S01 <- NTZ.F.Ages.S01[[1]] %>% 
  mutate(Movement = "Fast")
F.S01 <- NTZ.F.Ages.S01[[2]]


## S02 - Temporal Closure, no NTZs
NTZ.F.Ages.S02 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S02, f.list = SP_Pop_NTZ_S02, scenario.name = Names[3], nsim = 100)

NTZ.S02 <- NTZ.F.Ages.S02[[1]] %>% 
  mutate(Movement = "Fast")
F.S02 <- NTZ.F.Ages.S02[[2]]


## S03 - Temporal Closure and NTZs
NTZ.F.Ages.S03 <- zone.pop.format.full(ntz.list = SP_Pop_NTZ_S03, f.list = SP_Pop_NTZ_S03, scenario.name = Names[4], nsim = 100)

NTZ.S03 <- NTZ.F.Ages.S03[[1]] %>% 
  mutate(Movement = "Fast")
F.S03 <- NTZ.F.Ages.S03[[2]]


## Put everything together
Whole_Pop_Ages_NTZ <- rbind(NTZ.S00, NTZ.S01, NTZ.S02, NTZ.S03) %>% 
  mutate(Zone = "NTZ")

Whole_Pop_Ages_F <- rbind(F.S00, F.S01, F.S02, F.S03) %>% 
  mutate(Zone = "F")

NTZ_data_fast <- Whole_Pop_Ages_NTZ %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")

F_data_fast <- Whole_Pop_Ages_F %>% 
  pivot_longer(cols=starts_with("V"), names_to = "Simulation",values_to="Number")


differences <- NTZ_data_fast %>% 
  filter(Mod_Year %in% 2018) %>% 
  mutate(ID = rep(1:100, times=4)) %>%
  #mutate(MatBio = log(MatBio)) %>% 
  pivot_wider(id_cols=c(ID, Stage), names_from=Scenario, values_from=Number) %>% 
  mutate(Comp.1.2 = `Current NTZs` - `No temporal management or NTZs`,
         Comp.1.3 = `Current NTZs`- `Temporal management and NTZs`,
         Comp.1.4 = `Current NTZs` - `Temporal management`,
         Comp.2.3 = `No temporal management or NTZs` - `Temporal management and NTZs`,
         Comp.2.4 = `No temporal management or NTZs` - `Temporal management`,
         Comp.3.4 = `Temporal management and NTZs` - `Temporal management`) %>% 
  mutate(Comp.1.2.Perc = ((`Current NTZs` - `No temporal management or NTZs`)/((`Current NTZs`+`No temporal management or NTZs`)/2))*100,
         Comp.1.3.Perc = ((`Current NTZs`- `Temporal management and NTZs`)/((`Current NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.1.4.Perc = ((`Current NTZs` - `Temporal management`)/((`Current NTZs`+`Temporal management`)/2))*100,
         Comp.2.3.Perc = ((`No temporal management or NTZs` - `Temporal management and NTZs`)/((`No temporal management or NTZs`+`Temporal management and NTZs`)/2))*100,
         Comp.2.4.Perc = ((`No temporal management or NTZs` - `Temporal management`)/((`No temporal management or NTZs`+`Temporal management`)/2))*100,
         Comp.3.4.Perc = ((`Temporal management and NTZs` - `Temporal management`)/((`Temporal management and NTZs`+`Temporal management`)/2)*100))

Summary.Perc.Change <- differences %>% 
  ungroup() %>% 
  #group_by(Stage) %>% 
  summarise_at(vars(Comp.1.2.Perc:Comp.3.4.Perc), median)

## NTZ models 
# Recruits 
NTZ_Recruits_fast <- NTZ_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.NTZ.fast <- lm(log(Number) ~ Scenario, dat=NTZ_Recruits_fast)
summary(mod1.NTZ.fast)
plot(mod1.NTZ.fast$residuals)

mod1.NTZ.fast.aov <- aov(mod1.NTZ.fast)
mod1.NTZ.fast.tukey <- TukeyHSD(mod1.NTZ.fast.aov)
mod1.NTZ.fast.tukey

F_Recruits_fast <- F_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Recruit")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod1.F.fast <- lm(log(Number) ~ Scenario, dat=F_Recruits_fast)
summary(mod1.F.fast)
plot(mod1.F.fast$residuals)

mod1.F.fast.aov <- aov(mod1.F.fast)
mod1.F.fast.tukey <- TukeyHSD(mod1.F.fast.aov)
mod1.F.fast.tukey

# Sublegal
NTZ_Sublegal_fast <- NTZ_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.NTZ.fast <- lm(log(Number) ~ Scenario, dat=NTZ_Sublegal_fast)
summary(mod2.NTZ.fast)
plot(mod2.NTZ.fast$residuals)

mod2.NTZ.fast.aov <- aov(mod2.NTZ.fast)
mod2.NTZ.fast.tukey <- TukeyHSD(mod2.NTZ.fast.aov)
mod2.NTZ.fast.tukey


F_Sublegal_fast <- F_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Sublegal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod2.F.fast <- lm(log(Number) ~ Scenario, dat=F_Sublegal_fast)
summary(mod2.F.fast)
plot(mod2.F.fast$residuals)

mod2.F.fast.aov <- aov(mod2.F.fast)
mod2.F.fast.tukey <- TukeyHSD(mod2.F.fast.aov)
mod2.F.fast.tukey

# Legal
NTZ_Legal_fast <- NTZ_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.NTZ.fast <- lm(log(Number) ~ Scenario, dat=NTZ_Legal_fast)
summary(mod3.NTZ.fast)
plot(mod3.NTZ.fast$residuals)

mod3.NTZ.fast.aov <- aov(mod3.NTZ.fast)
mod3.NTZ.fast.tukey <- TukeyHSD(mod3.NTZ.fast.aov)
mod3.NTZ.fast.tukey


F_Legal_fast <- F_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod3.F.fast <- lm(log(Number) ~ Scenario, dat=F_Legal_fast)
summary(mod3.F.fast)
plot(mod3.F.fast$residuals)

mod3.F.fast.aov <- aov(mod3.F.fast)
mod3.F.fast.tukey <- TukeyHSD(mod3.F.fast.aov)
mod3.F.fast.tukey



# Large Legal
NTZ_Large_fast <- NTZ_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario)) 

mod4.NTZ.fast <- lm(log(Number) ~ Scenario, dat=NTZ_Large_fast)
summary(mod4.NTZ.fast)
plot(mod4.NTZ.fast$residuals)

mod4.NTZ.fast.aov <- aov(mod4.NTZ.fast)
mod4.NTZ.fast.tukey <- TukeyHSD(mod4.NTZ.fast.aov)
mod4.NTZ.fast.tukey


F_Large_fast <- F_data_fast %>% 
  filter(Mod_Year %in% c(2018)) %>% 
  filter(Stage %in% c("Large Legal")) %>% 
  mutate(Mod_Year = as.factor(Mod_Year),
         Scenario = as.factor(Scenario))

mod4.F.fast <- lm(log(Number) ~ Scenario, dat=F_Large_fast)
summary(mod4.F.fast)
plot(mod4.F.fast$residuals)

mod4.F.fast.aov <- aov(mod4.F.fast)
mod4.F.fast.tukey <- TukeyHSD(mod4.F.fast.aov)
mod4.F.fast.tukey

#* Comparing movement across scenarios ####
All_movement_NTZ <- rbind(NTZ_data_Medium, NTZ_data_slow, NTZ_data_fast) %>% 
  filter(Mod_Year %in% 2018)

S00_movement <- All_movement_NTZ %>% 
  filter(Scenario %in% "Historical and Current NTZs") %>% 
  filter(Stage %in% "Large Legal")

mod.S00 <- lm(log(Number) ~ Movement, data=S00_movement)

mod.S00.aov <- aov(mod.S00)
tukey.S00 <- TukeyHSD(mod.S00.aov)
tukey.S00 

S01_movement <- All_movement_NTZ %>% 
  filter(Scenario %in% "Neither NTZs nor Temporal Management") %>% 
  filter(Stage %in% "Large Legal")

mod.S01 <- lm(log(Number) ~ Movement, data=S01_movement)

mod.S01.aov <- aov(mod.S01)
tukey.S01 <- TukeyHSD(mod.S01.aov)
tukey.S01 

S02_movement <- All_movement_NTZ %>% 
  filter(Scenario %in% "Temporal Management Only") %>% 
  filter(Stage %in% "Large Legal")

mod.S02 <- lm(log(Number) ~ Movement, data=S02_movement)

mod.S02.aov <- aov(mod.S02)
tukey.S02 <- TukeyHSD(mod.S02.aov)
tukey.S02 

S03_movement <- All_movement_NTZ %>% 
  filter(Scenario %in% "Temporal and Spatial Management")%>% 
  filter(Stage %in% "Large Legal")

mod.S03 <- lm(log(Number) ~ Movement, data=S03_movement)

mod.S03.aov <- aov(mod.S03)
tukey.S03 <- TukeyHSD(mod.S03.aov)
tukey.S03 