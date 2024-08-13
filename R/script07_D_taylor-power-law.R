
###########################################################################
### Script Taylor power law
###########################################################################

# Import data -------------------------------------------------------------
cpue.obs <- read.table(here::here("data","processed","cpue_obs.txt"), h = T, sep = ";")
traits <- read.csv(here::here("data","raw","species_traits.csv"), h = T)

# Load functions ----------------------------------------------------------
source(here::here("R","functions","function_species_selection.R"))

# Define selected species -------------------------------------------------
species <- species_sel(abundance = cpue.obs)
species <- species[species != "Hopli.spp"]

# Function to estimate the Taylor exponent --------------------------------

###########################################################################
#### Function TPL
###########################################################################

##
# data                : Species abundance data matrix
# species             : Vector with species name
# season              : Season of the year
##

b.taylor <- function(data, species, season){
  
  u_season <- unique(season)                      # unique values
  n_lagos <- length(season) / length(u_season)    # calculates the number of locations
  
  n_col <- length(species)     # number of columns in the matrix
  n_row <- length(u_season)      # number of rows in the matrix
  
  tab_resu <- as.data.frame(matrix(nrow = length(species), ncol = 5))        # create matrix to save results
  colnames(tab_resu) <- c("specie","slope", "p.value", "r2", "n")             # rename columns in matrix
  
  logs <- as.data.frame(matrix(nrow = n_row, ncol = 2))                       # create matrix to save the mean log
  colnames(logs) <- c("log.med", "log.var")
  rownames(logs) <- u_season
  
  ## fill columns
  for(c in 1:n_col){
    
    p.especie <- species[c]
    tab_resu[c,"specie"] <- p.especie       
    
    ## fill rows
    for(l in 1:n_row){
      
      p <- which(season == u_season[l])       # cpue position per month
      
      log.mean <- NA
      log.variance <- NA
      
      mean <- mean(data[p, p.especie], na.rm = TRUE)      # mean
      variance <- var(data[p, p.especie], na.rm = TRUE)   # variance
      
      # Check if mean and variance are not NA before logging
      if(!is.na(mean) && mean > 0) log.mean <- log10(mean)
      if(!is.na(variance) && variance > 0) log.variance <- log10(variance)
      
      logs[l,"log.med"] <- log.mean        
      logs[l,"log.var"] <- log.variance    
    }
    
    ols <- lm(data = logs, log.var ~ log.med, na.action = na.exclude)   # regression between mean and variance
    s.model <- summary(ols)                                             # model summary
    
    b <- s.model$coefficients["log.med","Estimate"]      # slope
    p <- s.model$coefficients["log.med","Pr(>|t|)"]      # p-value
    r2 <- s.model$r.squared                              # r-square
    n <- nobs(ols)[1]                                    # number of observations
    
    tab_resu[c,"slope"] <- b                              
    tab_resu[c,"p.value"] <- p
    tab_resu[c,"r2"] <- r2
    tab_resu[c,"n"] <- n
  }
  
  return(tab_resu)
}

tpl <- b.taylor(data = cpue.obs, species = species, season = cpue.obs$ano.mes)

# Join matrices -----------------------------------------------------------
tpl <- dplyr::left_join(traits, tpl, by = "specie")

# Export tpl value --------------------------------------------------------
write.table(tpl, file = here::here("data","processed","slope_taylor.txt"), sep = ";", row.names = FALSE)
