
###########################################################################
#### Create matrices of aggregation indices
###########################################################################

# Import data -------------------------------------------------------------
cpue.obs <- read.table(here::here("data","processed","cpue_obs.txt"), h = T, sep = ";")
traits <- read.csv(here::here("data","raw","species_traits.csv"), h = T)

# Load functions ----------------------------------------------------------
source(here::here("R","functions","function_species_selection.R"))
source(here::here("R","functions","function_aggregation_descriptive_table.R"))

## Select species
species <- species_sel(abundance = cpue.obs)
species <- species[-which(species == "Hopli.spp")]

# Function aggregation ----------------------------------------------------

###########################################################################
#### Function aggregation index
###########################################################################

##
# abundance         : Species abundance data matrix
# species           : Vector with species name
# season            : Season of the year
# index             : Aggregation index (IM = morisita, Cx = Green index, k = inverse k of the negative binomial, VM = variance/mean and Ip = Lloyd’s index of patchiness)
##

# Function aggregation ----------------------------------------------------
id.aggregation <- function(abundance, species, season, index = "IM"){
  
  ## Load package
  library(stringr)
  library(vegan)
  
  u.season <- unique(season)                        # removes unique season values
  n.lakes <- length(season) / length(u.season)      # calculate the number of lakes
  n.spp <- length(species)                          # species numbers
  n.season <- length(u.season)                      # seasons numbers
  
  ## Create matrix results
  tab.resu <- as.data.frame(matrix(nrow = n.season, ncol = (n.spp + 1)))    
  colnames(tab.resu)[n.spp + 1] <- "ano.mes"                                
  
  ## Define specie column 
  for(c in 1:n.spp){
    
    p.specie <- species[c]                  # specie name
    colnames(tab.resu)[c] <- p.specie       # rename column
    
    # Calculate aggregation for floodplain
    for(l in 1:n.season){
      
      p <- which(abundance[,"ano.mes"] == u.season[l])          # find cpue position in matrix
      fish.abund <- abundance[p, p.specie]                      # abundance values
      fish.mean <- mean(fish.abund)                             # mean of the abundance in lakes
      fish.var <- var(fish.abund)                               # variance of the abundance in lakes
      tab.resu[l, "ano.mes"] <- u.season[l]                     # rename the "ano.mes" column
      
      
      ## Choose the Morisita index
      if(index == "IM"){
        test.im <- sum(is.na(fish.abund))

        if(test.im > 0){
          id.morisita <- NA
          tab.resu[l, p.specie] <- id.morisita
        }
        if(test.im == 0){
          id.morisita <- vegan::dispindmorisita(fish.abund, unique.rm = F, na.rm = F)
          im <- id.morisita[1,1]
          ifelse(is.na(im),
                 tab.resu[l, p.specie] <- NA,
                 tab.resu[l, p.specie] <- im)
                 
          na <- sum(is.na(im))
          if(na == 0){
            max.value <- length(fish.abund)
            min.value <- 1-((length(fish.abund) - 1) / (sum(fish.abund, na.rm = T) - 1))
            if(id.morisita[1,1] > max.value | id.morisita[1,1] < min.value){
              print(paste("The", p.specie, "species presents values outside the range allowed by the",u.season[l],"index"))
            }
          }
        }
      }
      
      ## Choose the k of the negative binomial
      if(index == "k"){
        k <- (fish.mean ^ 2) / (fish.var - fish.mean)
        ifelse(is.na(k),
               tab.resu[l, p.specie] <- NA,
               tab.resu[l, p.specie] <- 1/k)
        
        na <- sum(is.na(k))
        if(na == 0){
          max.value <- length(fish.abund) - (1/fish.mean) + 0.001
          min.value <- - (1/fish.mean)
          if((1/k) > max.value | (1/k) < min.value){
            print(paste("The", p.specie, "species presents values outside the range allowed by the",u.season[l],"index"))
          }
        }
      }
      
      ## Choose the Green index
      if(index == "Cx"){                
        sup <- (fish.var / fish.mean) - 1
        inf <- sum(fish.abund) - 1
        cx <- sup / inf
        
        ifelse(is.na(cx),
               tab.resu[l, p.specie] <- NA,
               tab.resu[l, p.specie] <- cx)
        
        na <- sum(is.na(cx))
        if(na == 0){
          max.value <- 1 + 0.001
          min.value <- -(1/(sum(fish.abund) - 1))
          if(cx > max.value | cx < min.value){
            print(paste("The", p.specie, "species presents values outside the range allowed by the",u.season[l],"index"))
          }
        }
      }
      
      ## Choose the variance/mean index
      if(index == "VM"){
        vm <- fish.var / fish.mean
        
        ifelse(is.na(vm),
               tab.resu[l, p.specie] <- NA,
               tab.resu[l, p.specie] <- vm)
        
        na <- sum(is.na(vm))
        if(na == 0){
          max.value <- sum(fish.abund) + 0.001
          min.value <- 0
          if(vm > max.value | vm < min.value){
            print(paste("The", p.specie, "species presents values outside the range allowed by the",u.season[l],"index"))
          }
        }
      }
      
      ## Choose Lloyd’s index of patchiness
      if(index == "Ip"){
        lloyd <- 1 + ((fish.var  - fish.mean)/(fish.mean^2))
        ifelse(is.na(lloyd),
               tab.resu[l, p.specie] <- NA,
               tab.resu[l, p.specie] <- lloyd)
      }
    }
  }
  
  ## Import water level parana river
  lv.parana <- read.table(here::here("data","raw","mean_hydro.txt"), h = T)
  
  ## Function to generate water level delay matrix
  source(here::here("R","functions","function_hydrological_averages.R"))
  water.lag <- matriz_lag(lv.parana)
  
  ## Merge aggregation matrix and level water matrix
  aggre <- dplyr::left_join(tab.resu, water.lag, by = "ano.mes")
  
  return(aggre)
}



# Calculate aggregation ---------------------------------------------------
## Aggregation for cpue observed
var_med.obs <- id.aggregation(cpue.obs, species, cpue.obs[,"ano.mes"], index = "VM")
morisita.obs <- id.aggregation(cpue.obs, species, cpue.obs[,"ano.mes"], index = "IM")

indexes <- list(im = morisita.obs,
                vm = var_med.obs)

# Aggregation descriptive table -------------------------------------------
tab.aggregation <- tab_descriptive(indexes, traits)
tab.aggregation[,"scientific_name"] <- rownames(tab.aggregation)
tab.aggregation <- dplyr::left_join(traits[,c("scientific_name", "ancient_scientific_name", 
                                              "mean_regional_size", "spawning_type", "migration", 
                                              "parental_care", "diet", "order", "family")], tab.aggregation,  by = "scientific_name")


# Procrustes analysis -----------------------------------------------------
im.pro <- morisita.obs[,species]
vm.pro <- var_med.obs[,species]

## Data imputation
for(i in 1:ncol(im.pro)){
  mean.im <- mean(im.pro[,i], na.rm = T)
  mean.vm <- mean(vm.pro[,i], na.rm = T)
  
  im.pro[is.na(im.pro[,i]),i] <- mean.im
  vm.pro[is.na(vm.pro[,i]),i] <- mean.vm
}

## Standardization
im.pro <- vegan::decostand(im.pro, method = "standardize")
vm.pro <- vegan::decostand(vm.pro, method = "standardize")

## Procrustes
concord <- procrustes(im.pro, vm.pro)
protest(im.pro, vm.pro)


# Export aggregation index ------------------------------------------------
for(e in 1:length(indexes)){
  name <- paste("aggregation_",names(indexes)[e],".txt", sep = "")
  write.table(indexes[[e]], file = here::here("data","processed", name), sep = ";")
}


# Export descriptive table ------------------------------------------------
write.table(tab.aggregation, file = here::here("output","tables","aggregation_descriptive_table.txt"), sep = ";")

# End
