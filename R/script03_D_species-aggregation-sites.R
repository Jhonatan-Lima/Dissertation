
###########################################################################
### Species aggregation sites
###########################################################################

# Import data -------------------------------------------------------------
cpue.obs <- read.table(here::here("data","processed","cpue_obs.txt"), h = T, sep = ";")
traits <- read.csv(here::here("data","raw","species_traits.csv"), h = T)

# Load functions ----------------------------------------------------------
source(here::here("R","functions","function_species_selection.R"))

## Select species for n > 30
species <- species_sel(cpue.obs)
species <- species[-which(species == "Hopli.spp")]

## Cpue
cpue.obs <- cpue.obs[,c(species,"locais","ano.mes")]

# Function ----------------------------------------------------------------
ocorr_um <- function(abundance) {
  
  m_binary <- matrix(NA, nrow = nrow(abundance), ncol = ncol(abundance))
  
  for (i in 1:nrow(abundance)) {
    
    value <- abundance[i,]
    teste <- sum(value > 0, na.rm = T)
    
    if ( teste > 0 ) {
      
      h_value <- max(value, na.rm = T)
      m_binary[i,] <- ifelse(abundance[i,] == h_value, 1, 0)
    }
    else{ next }
  }
  return(m_binary)
}


resu <- as.data.frame(matrix(ncol=6,nrow = length(species)))
colnames(resu) <- c("lfec", "lgar", "lgua", "lpat", "lpve", "lven")

for(u in 1:length(species)){
  
  sp <- species[u]
    data <- reshape(cpue.obs[,c(sp,"locais","ano.mes")],idvar = "ano.mes",timevar= "locais", direction = "wide")
  row.names(data) <- data[,"ano.mes"]
  data <- data[,-1]
  colnames(data) <- c("lfec", "lgar", "lgua", "lpat", "lpve", "lven")
  data <- as.matrix(data)
  
  m_binary <- ocorr_um(data)

  line<- colSums(m_binary, na.rm = T)
  
  rownames(resu)[u] <- sp
  resu[u,] <- line
  
}

## By region
region <- as.data.frame(matrix(ncol=3,nrow = length(species)))
colnames(region) <- c("Baia", "Parana", "Ivinhema")

for(a in 1:length(species)){
  
  sp <- species[a]
  
  rownames(region)[a] <- sp
  region[a,"Baia"] <- sum(resu[a,c("lfec","lgua")], na.rm = T)
  region[a,"Parana"] <- sum(resu[a,c("lgar","lpve")], na.rm = T)
  region[a,"Ivinhema"] <- sum(resu[a,c("lpat","lven")], na.rm = T)
  
}

region$specie <- rownames(region)
oc.region <- dplyr::left_join(traits[,c("specie","scientific_name")], region, by = "specie")

## By specie
por_sp <- ocorr_um(region[,1:3])
colnames(por_sp) <- c("Baia", "Parana", "Ivinhema")
colSums(por_sp)

# Export descriptive table ------------------------------------------------
write.table(oc.region, file = here::here("output","tables","ocorr_region.txt"), sep = ";")


## End