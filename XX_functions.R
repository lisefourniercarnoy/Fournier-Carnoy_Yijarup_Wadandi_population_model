###################################################

# Script that has all the functions
# Functions are only needed when running the model
# Have functions for movement, mortality, recruitment
# and for plotting

###################################################

st_centroid_within_poly <- function (poly) { #returns true centroid if inside polygon otherwise makes a centroid inside the polygon
  
  # check if centroid is in polygon
  centroid <- poly %>% st_centroid() 
  in_poly <- st_within(centroid, poly, sparse = F)[[1]] 
  
  # if it is, return that centroid
  if (in_poly) return(centroid) 
  
  # if not, calculate a point on the surface and return that
  centroid_in_poly <- st_point_on_surface(poly) 
  return(centroid_in_poly)
}


#### Plotting ####

total.plot.func <- function (pop) {
  TimeSeries <- ggplot(pop)+
    geom_line(aes(x=Year, y=Tot.Pop)) +
    scale_x_continuous("Year", breaks = c(1960, 1970, 1980, 1990, 2000, 2010, 2020))+
    xlab("Year")+
    ylab("Total Population")+
    theme_classic()+
    geom_vline(xintercept=1987, linetype="dashed", color="grey20")+
    geom_vline(xintercept=2005, colour="grey20")+
    geom_vline(xintercept=2017, linetype="dashed", colour="grey20")
  
  return(TimeSeries)
}

spatial.plot.func <- function (area, pop, pop.breaks, pop.labels, colours){
  
  water <- area %>%
    mutate(pop = round(pop, digits=0)) %>% 
    mutate(pop_level = cut(pop, pop.breaks, include.lowest=T)) 
  
  nb.cols <- length(pop.breaks)
  mycols <- colorRampPalette(rev(brewer.pal(8, colours)))(nb.cols)
  
  map <- ggplot(water)+
    geom_sf(aes(fill=pop_level), lwd=0)+
    scale_fill_manual(name="Population", values= mycols, drop=FALSE)+
    #scale_color_manual(name="Fished", values=c("black", "black"))+
    theme_void()
  
  return(map)
  
}

age.plot.func <- function (pop, NTZs){
  
  if(YEAR >=27 & YEAR <= 45){
    NoTakeAges <- pop[c(NTZs[[1]]),12, ]
    FishedAges <- pop[-c(NTZs[[1]]),12, ]
    
  } else if(YEAR>45 & YEAR<=57){
    NoTakeAges <- pop[c(NTZs[[2]]),12, ]
    FishedAges <- pop[-c(NTZs[[2]]),12, ]
    
  } else if (YEAR >57){
    NoTakeAges <- pop[c(NTZs[[3]]),12, ]
    FishedAges <- pop[-c(NTZs[[3]]),12, ]
    
  } else if (YEAR<27){
    NoTakeAges <- as.data.frame(array(0, dim=c(100, 30)))
    FishedAges <- pop[ ,12, ]
  }
  
  NoTakeAges <- as.data.frame(colSums(NoTakeAges)) %>% 
    rename(Number = "colSums(NoTakeAges)") %>% 
    mutate(Status = "NTZ") %>% 
    mutate(Age = seq(1:30))
  
  
  FishedAges <- as.data.frame(colSums(FishedAges)) %>% 
    rename(Number = "colSums(FishedAges)") %>% 
    mutate(Status = "Fished") %>% 
    mutate(Age = seq(1:30))
  
  AllAges <- rbind(NoTakeAges, FishedAges) %>% 
    mutate(Age = as.factor(Age)) %>%
    group_by(Status) %>% 
    mutate(TotalZone=sum(Number)) %>% 
    ungroup() %>% 
    mutate(Proportion = (Number/TotalZone))
  
  barplot.ages <- ggplot(AllAges)+
    geom_bar(aes(y=Proportion, x=Age, fill=Status), position="dodge", stat="identity")+
    theme_classic()
  
  return(barplot.ages)
  
}

msy.plot.func <- function(yield, biomass,spawning){
  
  yield.plot <- ggplot()+
    geom_line(aes(x=yield[,1], y=yield[,2]))+
    theme_classic()+
    xlab("Fishing Mortality")+
    ylab("Yeild per Recruit")
  
  biomass.plot <- ggplot()+
    geom_line(aes(x=biomass[,1], y=biomass[,2]))+
    theme_classic()+
    xlab("Fishing Mortality")+
    ylab("Total Biomass")
  
  spawning.plot <- ggplot()+
    geom_line(aes(x=spawning[,1], y=spawning[,2]))+
    theme_classic()+
    xlab("Fishing Mortality")+
    ylab("Total Spawning Biomass")
  
  msy.plots <- list(yield.plot, biomass.plot, spawning.plot)
  return(msy.plots)
  
}

#### FORMATTING FUNCTIONS ####
#* Fish by age and by zone ####


zone.fish.age <- function(NTZ.Ages, F.Ages, nsim, max.age, n.years, n.scenario, start.year, end.year){
  
  NTZ_Ages <- array(0, dim=c(max.age,n.years,nsim))
  F_Ages <- array(0, dim=c(max.age,n.years,nsim))
  
  for(S in 1:n.scenario){
    
    SP_Pop_NTZ <- NTZ.Ages[[S]]
    SP_Pop_F <- F.Ages[[S]]
    
    for(SIM in 1:length(SP_Pop_NTZ)){
      
      temp <- as.data.frame(colSums(SP_Pop_NTZ[[SIM]])) %>% 
        mutate(Age = seq(1:max.age)) %>% 
        pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
        mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", ""))) %>% 
        mutate(Mod_Year = rep(start.year:end.year, length.out=nrow(.))) %>% 
        group_by(Age, Mod_Year) %>% 
        summarise(across(where(is.numeric), sum)) %>% 
        ungroup() %>%
        dplyr::select(!Num_Year) %>%
        pivot_wider(names_from="Mod_Year",id_col="Age", values_from="Number",values_fn = list(count=list)) %>% 
        dplyr::select(!Age) %>% 
        unlist()
      
      temp2 <- array(temp, dim=c(max.age,n.years))
      
      NTZ_Ages[,,SIM] <- temp2
      
      temp <- as.data.frame(colSums(SP_Pop_F[[SIM]])) %>% 
        mutate(Age = seq(1:max.age)) %>% 
        pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
        mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", ""))) %>% 
        mutate(Mod_Year = rep(start.year:end.year, length.out=nrow(.))) %>% 
        group_by(Age, Mod_Year) %>% 
        summarise(across(where(is.numeric), sum)) %>% 
        ungroup() %>% 
        pivot_wider(names_from="Mod_Year",id_col="Age", values_from="Number") %>% 
        dplyr::select(!Age) %>% 
        unlist()
      
      temp2 <- array(temp, dim=c(max.age,n.years))
      
      F_Ages[,,SIM] <- temp2
    }
    
    Age.Dist.NTZ[[S]] <- NTZ_Ages
    Age.Dist.F[[S]] <- F_Ages
    
  }
  
  Age.Dist[[1]] <- Age.Dist.NTZ
  Age.Dist[[2]] <- Age.Dist.F
  
  return(Age.Dist)
} 

#* Catch by age ####

age.catch <- function(Catches, NTZ.Cells, F.Cells, nsim){
  Age.Catch <- list()
  temp4 <- list()
  temp5 <- list()
  catch <- array(0, dim=c(NCELL,59,nsim))
  
  for (S in 1:4){
    
    temp <- Catches[[S]] # NTZ and F are the same files 
    
    for(SIM in 1:nsim){
      catch[,,SIM] <- temp[[SIM]]
    }
    
    temp2 <- catch[as.numeric(NTZ.Cells),,]
    temp3 <- catch[as.numeric(F.Cells),,]
    
    temp4[[S]] <- temp2
    
    temp5[[S]] <- temp3
    
  }
  
  Age.Catch[[1]] <- temp4
  Age.Catch[[2]] <- temp5
  
  return(Age.Catch)
  
}


#* Means  ####

median.biomass.func <- function(Age, Catch, Zones, Age.NTZ, Age.F, Catch.NTZ, Catch.F, nsim, nyears, yearstart, maxage, scenario.names, mat, kg){
  Ages.All <- list()
  Ages.F <- list()
  Ages.NTZ <- list()
  
  Ages.All <- NULL
  Ages.F <- NULL
  Ages.NTZ <- NULL
  
  # Catches.All <- list()
  # Catches.NTZ <- list()
  # Catches.F <- list()
  
  # catch <- array(0, dim=c(maxage,nyear,nsim))
  # catch.NTZ <- array(0, dim=c(maxage,nyear,nsim))
  # catch.F <- array(0, dim=c(maxage,nyear,nsim))
  
  for(S in 1:4){
    
    #*************AGES***************
    total_pop <- NULL
    mod_year <- seq(1960,2018,1)
    
    temp <- Age[[S]]
    
    median.ages <- NULL
    
    for(year in yearstart:nyears){
      
      age <- temp[,year,]
      
      age.median <- age %>% 
        as.data.frame()
      
      age.median <- age.median[ , ]*mat[,12]
      age.median <- age.median[ , ]*kg[,12]
      
      age.median <- colSums(age.median[ ,1:nsim]) %>%   
        as.data.frame() %>% 
        mutate(Year = mod_year[year]) %>% 
        rename(MatBio = ".")
      
      median.ages <- rbind(median.ages, age.median)
    }
    
    median.ages <- median.ages %>% 
      group_by(Year) %>% 
      summarise(Median_MatBio = median(MatBio), 
                P_0.025 = quantile(MatBio, probs=c(0.025)),
                P_0.975 = quantile(MatBio, probs=c(0.975))) %>% 
      mutate(Scenario = scenario.names[S])
    
    total_pop <- rbind(total_pop, median.ages)
    
    #Ages.All[[S]] <- total_pop
    
    
    if(Zones == TRUE){
      # Ages NTZ
      NTZ_pop <- NULL
      mod_year <- seq(1960,2018,1)
      
      temp <- Age.NTZ[[S]]
      
      median.ages <- NULL
      
      for(year in yearstart:nyears){
        
        age <- temp[,year,]
        
        age.median <- age %>% 
          as.data.frame()
        
        age.median <- age.median[ , ]*mat[,12]
        age.median <- age.median[ , ]*kg[,12]
        
        age.median <- colSums(age.median[ ,1:nsim]) %>%   
          as.data.frame() %>% 
          mutate(Year = mod_year[year]) %>% 
          rename(MatBio = ".")
        
        median.ages <- rbind(median.ages, age.median)
      }
      
      median.ages <- median.ages %>% 
        group_by(Year) %>% 
        summarise(Median_MatBio = median(MatBio), 
                  P_0.025 = quantile(MatBio, probs=c(0.025)),
                  P_0.975 = quantile(MatBio, probs=c(0.975))) %>% 
        mutate(Scenario = scenario.names[S])
      
      NTZ_pop <- rbind(NTZ_pop, median.ages)
      
      #Ages.NTZ[[S]] <- NTZ_pop
      
      # Ages F
      F_pop <- NULL
      mod_year <- seq(1960,2018,1)
      
      temp <- Age.F[[S]]
      
      median.ages <- NULL
      
      for(year in yearstart:nyears){
        
        age <- temp[,year,]
        
        age.median <- age %>% 
          as.data.frame()
        
        age.median <- age.median[ , ]*mat[,12]
        age.median <- age.median[ , ]*kg[,12]
        
        age.median <- colSums(age.median[ ,1:nsim]) %>%   
          as.data.frame() %>% 
          mutate(Year = mod_year[year]) %>% 
          rename(MatBio = ".")
        
        median.ages <- rbind(median.ages, age.median)
      }
      
      median.ages <- median.ages %>% 
        group_by(Year) %>% 
        summarise(Median_MatBio = median(MatBio), 
                  P_0.025 = quantile(MatBio, probs=c(0.025)),
                  P_0.975 = quantile(MatBio, probs=c(0.975))) %>% 
        mutate(Scenario = scenario.names[S])
      
      F_pop <- rbind(F_pop, median.ages)
      
      #Ages.F[[S]] <- F_pop
      
    } else {}
    
    
    # #*********************************
    # #*************CATCHES*************
    # ## All
    # temp.c <- Catch[[S]]
    # for(SIM in 1:nsim){
    #   
    #   catch[,,SIM] <- temp.c[[SIM]]
    #   
    # }
    # 
    # temp2.c <- catch[,yearstart:nyear, ] %>% 
    #   rowMedians(., dim=2) %>% 
    #   as.data.frame(.) 
    # 
    # Catches.All[[S]] <- temp2.c
    # 
    # if(Zones == TRUE){
    #   ## NTZ
    #   temp.c.NTZ <- Catch.NTZ[[S]]
    #   for(SIM in 1:nsim){
    #     
    #     catch.NTZ[,,SIM] <- temp.c.NTZ[[SIM]]
    #     
    #   }
    #   
    #   temp2.c.NTZ <- catch.NTZ[,yearstart:nyear, ] %>% 
    #     rowMedians(., dim=2) %>% 
    #     as.data.frame(.) 
    #   
    #   if(S==1|S==3){
    #     temp2.c.NTZ[,] <- 0
    #   } else { }
    #   
    #   Catches.NTZ[[S]] <- temp2.c.NTZ
    #   
    #   ## F
    #   temp.c.F <- Catch.F[[S]]
    #   for(SIM in 1:nsim){
    #     
    #     catch.F[,,SIM] <- temp.c.F[[SIM]]
    #     
    #   }
    #   
    #   temp2.c.F <- catch.F[,yearstart:nyear, ] %>% 
    #     rowMedians(., dim=2) %>% 
    #     as.data.frame(.) 
    #   
    #   Catches.F[[S]] <- temp2.c.F
    # } else { }
    # 
    # 
    # #*********************************  
    
    Ages.All <- rbind(Ages.All, total_pop)
    Ages.F <- rbind(Ages.F, F_pop)
    Ages.NTZ <- rbind(Ages.NTZ, NTZ_pop)
    
  }
  
  
  Age.and.Catch <- list()
  
  Age.and.Catch[[1]] <- Ages.All
  Age.and.Catch[[2]] <- Ages.F
  Age.and.Catch[[3]] <- Ages.NTZ
  
  # Age.and.Catch[[4]] <- Catches.All
  # Age.and.Catch[[5]] <- Catches.NTZ
  # Age.and.Catch[[6]] <- Catches.F
  
  return(Age.and.Catch)
  
}

#* Biomass function #####

biomass.func <- function(Age, Zones, Age.NTZ, Age.F, Weights){
  Biomass.All <- list()
  Biomass.NTZ <- list()
  Biomass.F <- list()
  
  for(S in 1:4){
    ## All
    temp.All <- Age[[S]]
    
    temp3 <- temp.All %>% 
      rowMeans(., dim=2) %>% 
      as.data.frame(.) %>% 
      mutate_all(.,function(col){Weights$V12*col}) 
    
    Biomass.All[[S]] <- temp3
    
    if(Zones == TRUE){
      ## NTZ
      temp.NTZ <- Age.NTZ[[S]]
      
      temp3.NTZ <- temp.NTZ %>% 
        rowMeans(., dim=2) %>% 
        as.data.frame(.) %>% 
        mutate_all(.,function(col){Weights$V12*col}) 
      
      Biomass.NTZ[[S]] <- temp3.NTZ
      
      ## F
      temp.F <- Age.F[[S]]
      
      temp3.F<- temp.F %>% 
        rowMeans(., dim=2) %>% 
        as.data.frame(.) %>% 
        mutate_all(.,function(col){Weights$V12*col}) 
      
      Biomass.F[[S]] <- temp3.F
    } else {  }
  }
  
  
  Biomass <- list()
  Biomass[[1]] <- Biomass.All
  Biomass[[2]] <- Biomass.NTZ
  Biomass[[3]] <- Biomass.F
  
  return(Biomass)
  
}

#### Data formatting ####

## Whole population

total.pop.format <- function(pop.file.list, scenario.names, nsim, nyears, startyear, maxage, mat, kg){
  
  total_pop <- NULL
  mod_year <- seq(1980,2024,1)
  
  for(i in 1:length(pop.file.list)){
    temp <- pop.file.list[[i]]
    
    median.ages <- NULL
    
    for(year in startyear:nyears){
      
      age <- temp[,year,]
      
      age.median <- age %>% 
        as.data.frame()
      
      age.median <- age.median[ , ]*mat[,12]
      age.median <- age.median[ , ]*kg[,12]
      
      age.median <- colSums(age.median[ ,1:nsim]) %>%   
        as.data.frame() %>% 
        mutate(Year = mod_year[year]) %>% 
        rename(MatBio = ".")
      
      median.ages <- rbind(median.ages, age.median)
    }
    median.ages <- median.ages %>%
      group_by(Year) %>%
      summarise(Median_MatBio = median(MatBio),
                P_0.025 = quantile(MatBio, probs=c(0.025)),
                P_0.975 = quantile(MatBio, probs=c(0.975))) %>%
      mutate(Scenario = scenario.names[i])
    
    total_pop <- rbind(total_pop, median.ages)
    
  }
  
  return(total_pop)
  
}

## Whole population

total.pop.format.full <- function(pop.file.list, scenario.names, nsim, nyears, startyear, maxage, mat, kg){
  
  total_pop <- NULL
  mod_year <- seq(1960,2018,1)
  
  for(i in 1:length(pop.file.list)){
    temp <- pop.file.list[[i]]
    
    median.ages <- NULL
    
    for(year in startyear:nyears){
      
      age <- temp[,year,]
      
      age.median <- age %>% 
        as.data.frame()
      
      age.median <- age.median[ , ]*mat[,12]
      age.median <- age.median[ , ]*kg[,12]
      
      age.median <- colSums(age.median[ ,1:nsim]) %>%   
        as.data.frame() %>% 
        mutate(Year = mod_year[year]) %>% 
        rename(MatBio = ".")
      
      median.ages <- rbind(median.ages, age.median)
    }
    median.ages <- median.ages %>%
      group_by(Year) %>%
      # summarise(Median_MatBio = median(MatBio),
      #           P_0.025 = quantile(MatBio, probs=c(0.025)),
      #           P_0.975 = quantile(MatBio, probs=c(0.975))) %>%
      mutate(Scenario = scenario.names[i])
    
    total_pop <- rbind(total_pop, median.ages)
    
  }
  
  return(total_pop)
  
}

# Separate NTZs and Fished areas 

zone.pop.format <- function(ntz.list, f.list, scenario.name, nsim){
  
  NTZ_Ages_S00 <- NULL
  F_Ages_S00 <- NULL
  
  for(SIM in 1:nsim){
    
    temp <- as.data.frame(colSums(ntz.list[[SIM]])) %>% 
      mutate(Age = seq(1:30)) %>% 
      pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
      mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
    
    NTZ_Ages_S00 <- cbind(NTZ_Ages_S00, temp$Number)
    
    temp <- as.data.frame(colSums(f.list[[SIM]])) %>% 
      mutate(Age = seq(1:30)) %>% 
      pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
      mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
    
    F_Ages_S00 <- cbind(F_Ages_S00,temp$Number)
  }
  
  NTZ_Ages_S00 <- as.data.frame(NTZ_Ages_S00) %>% 
    mutate(Age = rep(1:30, each=59)) %>% 
    mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
    mutate(Stage = ifelse(Age==1, "Recruit",
                          ifelse(Age>1 & Age<3, "Sublegal",
                                 ifelse(Age>=3 & Age<=10, "Legal",
                                        ifelse(Age>10, "Large Legal",NA))))) %>% 
    group_by(Stage, Mod_Year) %>% 
    summarise(across(where(is.numeric) & !Age, sum)) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric) & !Mod_Year, ~./AreaNT)) %>%
    mutate(Median_Pop = rowMedians(as.matrix(.[,3:nsim]))) %>%
    mutate(P_0.025 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.025))) %>%
    mutate(P_0.975 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.975))) %>%
    mutate(Scenario = scenario.name)  %>% 
    dplyr::select(Mod_Year, Stage, Median_Pop, P_0.025, P_0.975, Scenario)
  
  
  F_Ages_S00 <- as.data.frame(F_Ages_S00) %>% 
    mutate(Age = rep(1:30, each=59)) %>% 
    mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
    mutate(Stage = ifelse(Age==1, "Recruit",
                          ifelse(Age>1 & Age<3, "Sublegal",
                                 ifelse(Age>=3 & Age<=10, "Legal",
                                        ifelse(Age>10, "Large Legal",NA))))) %>% 
    group_by(Stage, Mod_Year) %>% 
    summarise(across(where(is.numeric) & !Age, sum)) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric) & !Mod_Year, ~./AreaFished)) %>%
    mutate(Median_Pop = rowMedians(as.matrix(.[,3:nsim]))) %>%
    mutate(P_0.025 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.025))) %>%
    mutate(P_0.975 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.975))) %>%
    mutate(Scenario = scenario.name)  %>% 
    dplyr::select(Mod_Year, Stage, Median_Pop, P_0.025, P_0.975, Scenario)
  
  NTZ.F.Ages <- list()
  NTZ.F.Ages[[1]] <- NTZ_Ages_S00
  NTZ.F.Ages[[2]] <- F_Ages_S00
  
  return(NTZ.F.Ages)
  
}

zone.pop.format.full <- function(ntz.list, f.list, scenario.name, nsim){
  
  NTZ_Ages_S00 <- NULL
  F_Ages_S00 <- NULL
  
  for(SIM in 1:nsim){
    
    temp <- as.data.frame(colSums(ntz.list[[SIM]])) %>% 
      mutate(Age = seq(1:30)) %>% 
      pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
      mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
    
    NTZ_Ages_S00 <- cbind(NTZ_Ages_S00, temp$Number)
    
    temp <- as.data.frame(colSums(f.list[[SIM]])) %>% 
      mutate(Age = seq(1:30)) %>% 
      pivot_longer(cols=V1:V59, names_to="Num_Year", values_to="Number") %>% 
      mutate(Num_Year = as.numeric(str_replace(Num_Year, "V", "")))
    
    F_Ages_S00 <- cbind(F_Ages_S00,temp$Number)
  }
  
  NTZ_Ages_S00 <- as.data.frame(NTZ_Ages_S00) %>% 
    mutate(Age = rep(1:30, each=59)) %>% 
    mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
    mutate(Stage = ifelse(Age==1, "Recruit",
                          ifelse(Age>1 & Age<3, "Sublegal",
                                 ifelse(Age>=3 & Age<=10, "Legal",
                                        ifelse(Age>10, "Large Legal",NA))))) %>% 
    group_by(Stage, Mod_Year) %>% 
    summarise(across(where(is.numeric) & !Age, sum)) %>% 
    # ungroup() %>% 
    # mutate(across(where(is.numeric) & !Mod_Year, ~./AreaNT)) %>%
    # mutate(Median_Pop = rowMedians(as.matrix(.[,3:nsim]))) %>%
    # mutate(P_0.025 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.025))) %>%
    # mutate(P_0.975 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.975))) %>%
    mutate(Scenario = scenario.name) # %>% 
  # dplyr::select(Mod_Year, Stage, Median_Pop, P_0.025, P_0.975, Scenario)
  
  
  F_Ages_S00 <- as.data.frame(F_Ages_S00) %>% 
    mutate(Age = rep(1:30, each=59)) %>% 
    mutate(Mod_Year = rep(1960:2018, length.out=nrow(.))) %>% 
    mutate(Stage = ifelse(Age==1, "Recruit",
                          ifelse(Age>1 & Age<3, "Sublegal",
                                 ifelse(Age>=3 & Age<=10, "Legal",
                                        ifelse(Age>10, "Large Legal",NA))))) %>% 
    group_by(Stage, Mod_Year) %>% 
    summarise(across(where(is.numeric) & !Age, sum)) %>% 
    ungroup() %>% 
    # mutate(across(where(is.numeric) & !Mod_Year, ~./AreaFished)) %>%
    # mutate(Median_Pop = rowMedians(as.matrix(.[,3:nsim]))) %>%
    # mutate(P_0.025 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.025))) %>%
    # mutate(P_0.975 = rowQuantiles(as.matrix(.[,3:nsim]), probs=c(0.975))) %>%
    mutate(Scenario = scenario.name)  #%>% 
  # dplyr::select(Mod_Year, Stage, Median_Pop, P_0.025, P_0.975, Scenario)
  
  NTZ.F.Ages <- list()
  NTZ.F.Ages[[1]] <- NTZ_Ages_S00
  NTZ.F.Ages[[2]] <- F_Ages_S00
  
  return(NTZ.F.Ages)
  
}


distance.abundance.format <- function(Pops, n.year, n.sim, n.row, n.dist, n.cell, n.scenario, fished.cells, 
                                      LQ, HQ, distances, dist.names, scen.names, start.year, max.year,
                                      start.year.mod, end.year.mod){
  
  temp <- array(0, dim=c(n.cell, max.year, n.sim))
  
  Pop.Dist.Mean <- list()
  Pop.Dist.SD <- list()
  Pop.Dist.Median <- list()
  Pop.Dist.Quant <- list()
  
  Quantiles.10 <- array(0, dim=c(n.row, 2))
  Quantiles.50 <- array(0, dim=c(n.row, 2))
  Quantiles.100 <- array(0, dim=c(n.row, 2))
  
  
  for(S in 1:n.scenario){
    
    Scenario <- Pops[[S]]
    
    Means <- array(0, dim=c(n.dist, max.year)) %>% 
      as.data.frame(.)
    SDs <- array(0, dim=c(n.dist, max.year))%>% 
      as.data.frame(.)
    Medians <- array(0, dim=c(n.dist, max.year)) %>% 
      as.data.frame(.)
    Quantiles <- array(0, dim=c((n.dist*2), max.year)) %>% 
      as.data.frame(.)
    
    for(SIM in 1:n.sim){
      temp[,,SIM] <- Scenario[[SIM]]
    }
    for(YEAR in start.year:max.year){
      
      temp1 <- temp[,YEAR,]
      
      temp2 <- as.data.frame(temp1) %>% 
        mutate(CellID = row_number()) %>% 
        mutate(Fished_17 = fished.cells$Fished_2017) 
      
      Dist.10km <- temp2 %>% 
        filter(CellID %in% c(distances[[1]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(across(where(is.numeric), sum)) %>% 
        mutate(Mean_Pop = rowMeans(.[1:n.sim])) %>% 
        mutate(SD_Pop = rowSds(as.matrix(.[,1:n.sim]))) %>% 
        mutate(Median_Pop = rowMedians(as.matrix(.[1:n.sim]))) %>% 
        mutate(Distance = dist.names[1]) 
      
      Dist.50km <- temp2 %>% 
        filter(CellID %in% c(distances[[2]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(across(where(is.numeric), sum)) %>% 
        mutate(Mean_Pop = rowMeans(.[1:n.sim])) %>% 
        mutate(SD_Pop = rowSds(as.matrix(.[,1:n.sim]))) %>% 
        mutate(Median_Pop = rowMedians(as.matrix(.[1:n.sim]))) %>% 
        mutate(Distance = dist.names[2])
      
      Dist.100km <- temp2 %>% 
        filter(CellID %in% c(distances[[3]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(across(where(is.numeric), sum)) %>% 
        mutate(Mean_Pop = rowMeans(.[1:n.sim])) %>% 
        mutate(SD_Pop = rowSds(as.matrix(.[,1:n.sim]))) %>% 
        mutate(Median_Pop = rowMedians(as.matrix(.[1:n.sim]))) %>% 
        mutate(Distance = dist.names[3]) 
      
      temp2 <- as.matrix(Dist.10km[ , 1:n.sim])
      
      Quantiles.10 <- quantile(temp2, probs=c(LQ, HQ)) %>% 
        as.data.frame() %>% 
        mutate(Distance = dist.names[1]) 
      
      temp2 <- as.matrix(Dist.50km[ , 1:n.sim])
      
      Quantiles.50 <-quantile(temp2, probs=c(LQ, HQ)) %>% 
        as.data.frame() %>% 
        mutate(Distance = dist.names[2]) 
      
      temp2 <- as.matrix(Dist.100km[ , 1:n.sim])
      
      Quantiles.100 <-quantile(temp2, probs=c(LQ, HQ)) %>% 
        as.data.frame() %>% 
        mutate(Distance = dist.names[3]) 
      
      Means.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
        dplyr::select(Distance, Mean_Pop)
      SDs.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
        dplyr::select(Distance, SD_Pop)
      Medians.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>% 
        dplyr::select(Distance, Median_Pop)
      Quantiles.Full <- rbind(Quantiles.10, Quantiles.50, Quantiles.100) %>%
        as.data.frame() 
      Quantiles.Full$Quantile <- rownames(Quantiles.Full)
      
      Means[,YEAR] <- Means.Full$Mean_Pop
      
      SDs[ ,YEAR] <- SDs.Full$SD_Pop
      
      Medians[ ,YEAR] <- Medians.Full$Median_Pop
      
      Quantiles[, YEAR] <- Quantiles.Full$. 
      Quantiles$Quantile <- Quantiles.Full$Quantile
      
      
    }
    Means <- as.data.frame(Means) %>% 
      mutate(Scenario = scen.names[S]) %>% 
      mutate(Distances = dist.names) %>% 
      pivot_longer(cols=-c("Scenario", "Distances"), values_to = "Mean_Abundance", names_to = "Year") %>% 
      mutate(Year = rep(seq(start.year.mod,end.year.mod,1), times = 3))
    
    SDs <- as.data.frame(SDs) %>% 
      mutate(Scenario = scen.names[S]) %>% 
      mutate(Distances = dist.names)%>% 
      pivot_longer(cols=-c("Scenario", "Distances"), values_to = "SD_Abundance", names_to = "Year") %>% 
      mutate(Year = rep(seq(start.year.mod,end.year.mod,1), times = 3))
    
    Medians <- as.data.frame(Medians) %>% 
      mutate(Scenario = scen.names[S]) %>% 
      mutate(Distances = dist.names)%>% 
      pivot_longer(cols=-c("Scenario", "Distances"), values_to = "Median_Abundance", names_to = "Year") %>% 
      mutate(Year = rep(seq(start.year.mod,end.year.mod,1), times = 3))
    
    Quantiles <- as.data.frame(Quantiles) %>%
      mutate(Scenario = rep(scen.names[S], 6)) %>%
      mutate(Distances = rep(Dist.Names, each=2)) %>%
      pivot_longer(cols=-c("Scenario", "Distances", "Quantile"), values_to = "Quantile_Abundance", names_to = "Year") %>%
      mutate(Year = rep(seq(start.year.mod,end.year.mod,1), times = 6)) %>% 
      mutate(Quantile = str_replace(Quantile, "%[[:digit:]]$", "%")) %>% 
      pivot_wider(names_from = "Quantile", values_from="Quantile_Abundance")
    
    Pop.Dist.Mean[[S]] <- Means
    Pop.Dist.SD[[S]] <- SDs
    Pop.Dist.Median[[S]] <- Medians
    Pop.Dist.Quant[[S]] <- Quantiles
    
  }
  results <- list()
  results[[1]] <- Pop.Dist.Median
  results[[2]] <- Pop.Dist.Quant
  results[[3]] <- Pop.Dist.Mean
  results[[4]] <- Pop.Dist.SD
  
  return(results)
}

distance.catch.format <- function(Pops, n.year, n.sim, n.row, n.dist, n.cell, n.scenario, fished.cells, 
                                  LQ, HQ, distances, dist.names, scen.names, start.year, max.year,
                                  start.year.mod, end.year.mod){
  
  temp <- array(0, dim=c(n.cell, max.year, n.sim))
  
  Pop.Catch.Mean <- list()
  Pop.Catch.SD <- list()
  Pop.Catch.Median <- list()
  Pop.Catch.Quant<- list()
  
  for(S in 1:4){
    
    Scenario <- Pops[[S]]
    
    Means <- array(0, dim=c(n.dist, max.year)) %>%
      as.data.frame(.)
    SDs <- array(0, dim=c(n.dist, max.year))%>%
      as.data.frame(.)
    Medians <- array(0, dim=c(n.dist, n.row))%>% 
      as.data.frame(.)
    Quantiles <- array(0, dim=c((n.dist*2), n.row))%>% 
      as.data.frame(.)
    
    for(SIM in 1:n.sim){
      temp[,,SIM] <- Scenario[[SIM]]
    }
    
    for(YEAR in 1:max.year){
      
      temp1 <- temp[,YEAR,]
      
      temp2 <- as.data.frame(temp1) %>% 
        mutate(CellID = row_number()) %>% 
        mutate(Fished_17 = fished.cells$Fished_2017) 
      
      Dist.10km <- temp2 %>% 
        filter(CellID %in% c(distances[[1]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(across(where(is.numeric), sum)) %>%
        mutate(Mean_Pop = rowMeans(.[1:n.sim])) %>%
        mutate(SD_Pop = rowSds(as.matrix(.[,1:n.sim]))) %>%
        mutate(Median_Pop = rowMedians(as.matrix(.[1:n.sim]))) %>%
        mutate(Distance = dist.names[1])
      
      Dist.50km <- temp2 %>% 
        filter(CellID %in% c(distances[[2]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(across(where(is.numeric), sum)) %>%
        mutate(Mean_Pop = rowMeans(.[1:n.sim])) %>%
        mutate(SD_Pop = rowSds(as.matrix(.[,1:n.sim]))) %>%
        mutate(Median_Pop = rowMedians(as.matrix(.[1:n.sim]))) %>%
        mutate(Distance = dist.names[2])
      # 
      Dist.100km <- temp2 %>% 
        filter(CellID %in% c(distances[[3]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(across(where(is.numeric), sum)) %>%
        mutate(Mean_Pop = rowMeans(.[1:n.sim])) %>%
        mutate(SD_Pop = rowSds(as.matrix(.[,1:n.sim]))) %>%
        mutate(Median_Pop = rowMedians(as.matrix(.[1:n.sim]))) %>%
        mutate(Distance = dist.names[3])
      
      temp2 <- as.matrix(Dist.10km[ , 1:n.sim])
      
      Catch.Quantiles.10 <- quantile(temp2, probs=c(LQ, HQ)) %>% 
        as.data.frame() %>% 
        mutate(Distance = "10 km") 
      
      temp2 <- as.matrix(Dist.50km[ , 1:n.sim])
      
      Catch.Quantiles.50 <-quantile(temp2, probs=c(LQ, HQ)) %>% 
        as.data.frame() %>% 
        mutate(Distance = "50 km") 
      
      temp2 <- as.matrix(Dist.100km[ , 1:n.sim])
      
      Catch.Quantiles.100 <-quantile(temp2, probs=c(LQ, HQ)) %>% 
        as.data.frame() %>% 
        mutate(Distance = "100 km") 
      
      Means.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>%
        dplyr::select(Distance, Mean_Pop)
      SDs.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>%
        dplyr::select(Distance, SD_Pop)
      Medians.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) %>%
        dplyr::select(Distance, Median_Pop)
      Quantiles.Full <- rbind(Catch.Quantiles.10, Catch.Quantiles.50, Catch.Quantiles.100) %>%
        as.data.frame()
      Quantiles.Full$Quantile <- rownames(Quantiles.Full)
      
      Means[,YEAR] <- Means.Full$Mean_Pop
      
      SDs[ ,YEAR] <- SDs.Full$SD_Pop
      
      Medians[ ,YEAR] <- Medians.Full$Median_Pop
      
      Quantiles[, YEAR] <- Quantiles.Full$.
      Quantiles$Quantile <- Quantiles.Full$Quantile
    }
    
    Means <- as.data.frame(Means) %>%
      mutate(Scenario = scen.names[S]) %>%
      mutate(Distances = dist.names) %>%
      pivot_longer(cols=-c("Scenario", "Distances"), values_to = "Mean_Catch", names_to = "Year") %>%
      mutate(Year = rep(seq(start.year.mod,end.year.mod,1), times = 3))
    
    SDs <- as.data.frame(SDs) %>%
      mutate(Scenario = scen.names[S]) %>%
      mutate(Distances = dist.names) %>%
      pivot_longer(cols=-c("Scenario", "Distances"), values_to = "SD_Catch", names_to = "Year") %>%
      mutate(Year = rep(seq(start.year.mod,end.year.mod,1), times = 3))
    
    Medians <- as.data.frame(Medians) %>% 
      mutate(Scenario = scen.names[S]) %>%
      mutate(Distances = dist.names) %>%
      pivot_longer(cols=-c("Scenario", "Distances"), values_to = "Median_Catch", names_to = "Year") %>%
      mutate(Year = rep(seq(start.year.mod,end.year.mod,1), times = 3))
    
    Quantiles <- as.data.frame(Quantiles) %>%
      mutate(Scenario = rep(scen.names[S], 6)) %>%
      mutate(Distances = rep(dist.names, each=2, times=1)) %>%
      pivot_longer(cols=-c("Scenario", "Distances", "Quantile"), values_to = "Quantile_Catch", names_to = "Year") %>%
      mutate(Year = rep(seq(start.year.mod,end.year.mod,1), times = 6)) %>%
      mutate(Quantile = str_replace(Quantile, "%[[:digit:]]$", "%")) %>%
      pivot_wider(names_from = "Quantile", values_from="Quantile_Catch")
    
    Pop.Catch.Mean[[S]] <- Means
    Pop.Catch.SD[[S]] <- SDs
    Pop.Catch.Median[[S]] <- Medians
    Pop.Catch.Quant[[S]] <- Quantiles
    
  }
  results <- list()
  results[[1]] <- Pop.Catch.Median
  results[[2]] <- Pop.Catch.Quant
  results[[3]] <- Pop.Catch.Mean
  results[[4]] <- Pop.Catch.SD
  
  return(results)
}


dist.abundance.full <- function(Pops, n.year, n.cell, n.sim, n.scenario, fished.cells, 
                                distances, dist.names, scen.names, mod.years, target.year){
  Pop.Abundance <- list()
  
  for(S in 1:n.scenario){
    
    Scenario <- Pops[[S]]
    
    Abundance <- NULL
    
    temp <- array(0, dim=c(n.cell, n.year, n.sim))
    
    for(SIM in 1:n.sim){
      temp[,,SIM] <- Scenario[[SIM]]
    }
    for(YEAR in target.year:target.year){
      
      temp1 <- temp[,59,]
      
      temp2 <- as.data.frame(temp1) %>% 
        mutate(CellID = row_number()) %>% 
        mutate(Fished_17 = fished.cells$Fished_2017) 
      
      Dist.10km <- temp2 %>% 
        filter(CellID %in% c(distances[[1]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(Total=colSums(.[1:n.sim])) %>% 
        mutate(Distance = dist.names[1]) %>% 
        mutate(Scenario = scen.names[S]) %>%
        mutate(Year = mod.years[target.year])
      
      Dist.50km <- temp2 %>% 
        filter(CellID %in% c(distances[[2]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(Total=colSums(.[1:n.sim])) %>% 
        mutate(Scenario = scen.names[S]) %>%
        mutate(Distance = dist.names[2]) %>% 
        mutate(Year = mod.years[YEAR])
      
      Dist.100km <- temp2 %>% 
        filter(CellID %in% c(distances[[3]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(Total=colSums(.[1:n.sim])) %>% 
        mutate(Scenario = scen.names[S]) %>%
        mutate(Distance = dist.names[3])%>% 
        mutate(Year = mod.years[YEAR])
      
      Abundance.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) 
      
      Abundance<- rbind(Abundance, Abundance.Full)
      
    }
    Pop.Abundance[[S]] <- Abundance
  }
  return(Pop.Abundance) 
}

dist.catch.full <- function(Pops, n.year, n.cell, n.sim, n.scenario, fished.cells, 
                            distances, dist.names, scen.names, mod.years, target.year){
  Pop.Catch <- list()
  
  for(S in 1:n.scenario){
    
    Scenario <- Pops[[S]]
    
    Catch <- NULL
    temp <- array(0, dim=c(n.cell, n.year, n.sim))
    
    for(SIM in 1:n.sim){
      temp[,,SIM] <- Scenario[[SIM]]
    }
    for(YEAR in target.year:target.year){
      
      temp1 <- temp[,59,]
      
      temp2 <- as.data.frame(temp1) %>% 
        mutate(CellID = row_number()) %>% 
        mutate(Fished_17 = fished.cells$Fished_2017) 
      
      Dist.10km <- temp2 %>% 
        filter(CellID %in% c(distances[[1]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(Total=colSums(.[1:n.sim])) %>% 
        mutate(Distance = dist.names[1]) %>% 
        mutate(Scenario = scen.names[S]) %>%
        mutate(Year = mod.years[target.year])
      
      Dist.50km <- temp2 %>% 
        filter(CellID %in% c(distances[[2]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(Total=colSums(.[1:n.sim])) %>% 
        mutate(Scenario = scen.names[S]) %>%
        mutate(Distance = dist.names[2]) %>% 
        mutate(Year = mod.years[YEAR])
      
      Dist.100km <- temp2 %>% 
        filter(CellID %in% c(distances[[3]])) %>% 
        {if (S>0) filter(., Fished_17 =="Y") else .} %>% 
        summarise(Total=colSums(.[1:n.sim])) %>% 
        mutate(Scenario = scen.names[S]) %>%
        mutate(Distance = dist.names[3])%>% 
        mutate(Year = mod.years[YEAR])
      Catch.Full <- rbind(Dist.10km, Dist.50km, Dist.100km) 
      
      Catch <- rbind(Catch, Catch.Full)
      
    }
    Pop.Catch[[S]] <- Catch
  }
  return(Pop.Catch) 
}


#### PLOTTING FUNCTIONS ####

age.group.plots <- function(age.group, data.to.plot, plot.label.1, plot.label.2, label.pos.y, label.pos.x){
  
  plot.NTZ <- data.to.plot %>% 
    filter(Stage %in% c(age.group)) %>% 
    filter(Zone %in% c("NTZ")) %>% 
    filter(Mod_Year >=1986) %>% 
    mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Current NTZs") & Mod_Year>1985, "Current NTZs", 
                                                                   ifelse(Scenario %in% c("Temporal management") & Mod_Year>1985, "Temporal management", 
                                                                          ifelse(Scenario %in% c("No temporal management or NTZs"), "No temporal management\nor NTZs", "Temporal management\nand NTZs")))))%>% 
    mutate(ColourGroup = as.factor(ColourGroup)) %>% 
    ggplot(.)+
    geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
    geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
    scale_fill_manual(values= c("Current NTZs"="#36753B", "No temporal management\nor NTZs"="#302383" ,"Temporal management\nand NTZs"="#66CCEE",
                                "Temporal management"="#BBCC33"),
                      guide="none")+
    scale_colour_manual(values = c("Temporal management\nand NTZs"="#66CCEE","Current NTZs"="#36753B", "Temporal management"="#BBCC33", 
                                   "No temporal management\nor NTZs"="#302383"), breaks= c("Current NTZs", "No temporal management\nor NTZs", "Temporal management", "Temporal management\nand NTZs"),name= "Spatial and temporal\nmanagement scenario")+ 
    theme_classic()+
    xlab(NULL)+
    ylab(NULL)+
    xlim(1986,2020)+
    ylim(0,NA)+
    scale_y_continuous(breaks = pretty_breaks(), limits=c(0,NA))+
    scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),guide = "none")+
    theme(legend.title = element_text(size=9), #change legend title font size
          legend.text = element_text(size=8), #change legend text font size
          legend.spacing.y = unit(0.1, "cm"),
          legend.key.size = unit(2,"line")) +
    guides(color = guide_legend(byrow = TRUE))+
    theme(axis.text=element_text(size=8))+
    geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
    geom_vline(xintercept=2005, colour="grey20")+
    geom_vline(xintercept=2017, linetype="dotted", colour="grey20") +
    ggplot2::annotate("text", x=label.pos.x, y=label.pos.y, label=plot.label.1, size = 2.5, fontface=1)
  plot.NTZ
  
  plot.F <- data.to.plot %>% 
    filter(Stage %in% c(age.group)) %>% 
    filter(Zone %in% c("F")) %>% 
    filter(Mod_Year >=1986) %>% 
    mutate(ColourGroup = ifelse(Mod_Year<=1985, "Pre-1987", ifelse(Scenario %in% c("Current NTZs") & Mod_Year>1985, "Current NTZs", 
                                                                   ifelse(Scenario %in% c("Temporal management") & Mod_Year>1985, "Temporal management", 
                                                                          ifelse(Scenario %in% c("No temporal management or NTZs"), "No temporal management\nor NTZs", "Temporal management\nand NTZs")))))%>% 
    mutate(ColourGroup = as.factor(ColourGroup)) %>% 
    ggplot(.)+
    geom_line(aes(x=Mod_Year, y=Median_Pop, group=interaction(Zone,Scenario), colour=ColourGroup, linetype=Zone), size=0.7)+
    geom_ribbon(aes(x=Mod_Year, y=Median_Pop, ymin=P_0.025, ymax=P_0.975, fill=ColourGroup, group=interaction(Zone,Scenario)), alpha=0.2)+
    scale_fill_manual(values= c("Current NTZs"="#36753B", "No temporal management\nor NTZs"="#302383" ,"Temporal management\nand NTZs"="#66CCEE",
                                "Temporal management"="#BBCC33"),
                      guide="none")+
    scale_colour_manual(values = c("Temporal management\nand NTZs"="#66CCEE","Current NTZs"="#36753B", "Temporal management"="#BBCC33", 
                                   "No temporal management\nor NTZs"="#302383"), breaks= c("Temporal management\nand NTZs", "Current NTZs", "Temporal management", "No temporal management\nor NTZs"),name= "Spatial and temporal\nmanagement scenario")+ 
    theme_classic()+
    xlab(NULL)+
    ylab(NULL)+
    xlim(1986,2020)+
    ylim(0,NA)+
    scale_y_continuous(breaks = pretty_breaks(), limits=c(0,NA))+
    scale_linetype_manual(values = c("solid", "longdash" ), breaks=c("NTZ", "F") ,labels=c("NTZ area", "Always fished"),guide = "none")+
    theme(legend.title = element_text(size=9), #change legend title font size
          legend.text = element_text(size=8), #change legend text font size
          legend.spacing.y = unit(0.1, "cm"),
          legend.key.size = unit(2,"line")) +
    guides(color = guide_legend(byrow = TRUE))+
    theme(axis.text=element_text(size=8))+
    geom_vline(xintercept=1986, linetype="dashed", color="grey20")+
    geom_vline(xintercept=2005, colour="grey20")+
    geom_vline(xintercept=2017, linetype="dotted", colour="grey20") +
    ggplot2::annotate("text", x=label.pos.x, y=label.pos.y, label=plot.label.2, size = 2.5, fontface=1)
  plot.F
  
  plots <- list()
  plots[[1]] <- plot.NTZ
  plots[[2]] <- plot.F
  
  return(plots)
  
}

## Work out abundance of each age group 

age.to.length <- function(NTZ.ages, F.ages, max.age, n.scenarios, Linf, k, t0, LM, BigLM){
  
  lengths.out.NTZ <- list()
  lengths.out.F <- list()
  lengths.out <- list()
  
  for(S in 1:n.scenarios){
    ages.NTZ <- NTZ.ages[[S]][,,50]
    ages.F <- F.ages[[S]][,,50]
    
    temp.NTZ <- as.data.frame(colSums(ages.NTZ))
    temp.F <- as.data.frame(colSums(ages.F))
    
    lengths.NTZ <- temp.NTZ %>% 
      mutate(Age = seq(1,max.age,1)) %>% 
      mutate(Length = Linf*(1-exp(-k*(Age-t0)))) %>% 
      mutate(Length.Group = ifelse(Length<LM, "Too Small",
                                   ifelse(Length>BigLM, "> Big LM", "> LM"))) %>% 
      rename(Count = "colSums(ages.NTZ)")
    
    lengths.F <- temp.F %>% 
      mutate(Age = seq(1,max.age,1)) %>% 
      mutate(Length = Linf*(1-exp(-k*(Age-t0)))) %>% 
      mutate(Length.Group = ifelse(Length<LM, "Too Small",
                                   ifelse(Length>BigLM, "> Big LM", "> LM"))) %>% 
      rename(Count = "colSums(ages.F)")
    
    lengths.out.NTZ[[S]] <- lengths.NTZ
    lengths.out.F[[S]] <- lengths.F
    
  }
  
  lengths.out[[1]] <- lengths.out.NTZ
  lengths.out[[2]] <- lengths.out.F
  
  return(lengths.out)
}





