
###########################################################################
### Organization of cpue data
###########################################################################

# Import data -------------------------------------------------------------
cpue.peld <- read.table(here::here("data","raw","cpue_peld.txt"), header = TRUE)
peld.complete <- read.table(here::here("data","raw","peld_complete.txt"), header = TRUE)
traits <- read.csv(here::here("data","raw","species_traits.csv"), header = TRUE)

# Load functions ----------------------------------------------------------
source(here::here("R","functions","function_species_selection.R"))
source(here::here("R","functions","function_add_lines.R"))

# Abundance (CPUE) --------------------------------------------------------
## Select only the lakes
sites <- c("lfec","lgar","lgua","lpat","lpve","lven")
cpue.obs <- cpue.peld[cpue.peld$locais %in% sites,]

## Select species for n > 30
species <- species_sel(cpue.obs)

## Add missing lines
years <- c(3, 17)
months <- c(3, 9)
cpue.obs <- add_lines(cpue.obs, species, sites, years, months)

## Export abundance
write.table(cpue.obs, file = here::here("data","processed","cpue_obs.txt"), sep = ";", row.names = FALSE)


###########################################################################
### Calculate the proportion of immature individuals
###########################################################################

##
# data              : Fish tabular matrix
# traits            : Matrix with species traits
# species           : Vector with species name
# sites             : Vector with sites name
##

# Function ----------------------------------------------------------------
prop_immature <- function(data, traits, species, sites){
  
  ## Filter sites
  data <- data[data$local %in% sites,]
  
  ## Filter species
  data <- data[data$sp %in% species,]
  
  ## Standardize column "ano.mes"
  data[,"ano.mes"] <- stringr::str_pad(data[,"ano.mes"], width = 4, pad = "0")
  data[,"ano"] <- substr(data[,"ano.mes"], 1, 2)
  data[,"ano"] <- as.numeric(data[,"ano"]) + 2000
  data[,"mes"] <- substr(data[,"ano.mes"], 3, 4)
  data[,"ano.mes"] <- paste(data[,"ano"], data[,"mes"], sep = "/")
  
  ## Create matrix results
  season <- unique(data[,"ano.mes"])
  results <- as.data.frame(matrix(ncol = (length(species) + 1), nrow = length(season)))
  colnames(results)[length(species) + 1] <- "ano.mes"
  
  ## Create for
  for (c in 1:length(species)) {
    
    sp <- species[c]
    colnames(results)[c] <- sp
    
    for (i in 1:length(season)) {
      p.sp <- which(traits[,"specie"] == sp)     # locate species position
      l50 <- traits[p.sp,"L50_mean"]             # fish maturation length
      u.season <- season[i]                               # season
      results[i,"ano.mes"] <- u.season
      
      ## Counting the total number of fish
      posi.total <- which(data[,"ano.mes"] == u.season & data[,"sp"] == sp & (data[,"LT"] > 0 | data[,"LS"] > 0))
      qtd.total <- length(posi.total)
      
      ## Counting the number of fish below the LM
      posi.lm <- which(data[,"ano.mes"] == u.season & data[,"sp"] == sp & (data[,"LS"] < l50 | data[,"LT"] < l50) & (data[,"LT"] > 0 | data[,"LS"] > 0)) 
      qtd.lm <- length(posi.lm)
      
      results[i,c] <- qtd.lm / qtd.total
    }
  }
  
  return(results)
}

## Calculate proportion of immature fish
immature <- prop_immature(peld.complete, traits, species, sites)

## Export proportion of immature
write.table(immature, file = here::here("data","processed","prop_immature.txt"), sep = ";", row.names = FALSE)

# End
