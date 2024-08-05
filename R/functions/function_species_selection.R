
###########################################################################
#### Function to Select Species with n > 30 in the Floodplain
###########################################################################

##
# abundance     : Species CPUE matrix
## 

# Function ------------------------------------------------------------------
species_sel <- function(abundance) {
  
  ## Add "ano.mes" column if it doesn't exist
  if (!"ano.mes" %in% colnames(abundance)) {
    abundance$anos <- as.numeric(abundance$anos) + 2000
    abundance$meses <- stringr::str_pad(abundance$meses, width = 2, pad = "0")
    abundance$ano.mes <- paste(abundance$anos, abundance$meses, sep = "/")
  }
  
  ## Define species names
  cols_to_remove <- c("ano.mes", "anos", "meses", "locais", 
                      paste0("nv.df", 0:5))
  species_names <- setdiff(colnames(abundance), cols_to_remove)
  seasons <- unique(abundance$ano.mes)
  
  ## Vector to save the names of selected species
  selected_species <- character()
  
  ## Loop through species
  for (species in species_names) {
    count <- 0
    
    for (season in seasons) {
      observations <- sum(abundance[abundance$ano.mes == season, species] > 0)
      observations <- ifelse(is.na(observations), 0, observations)
      if (observations > 0) {
        count <- count + 1
      }
    }
    
    if (count >= 30) {
      selected_species <- c(selected_species, species)
    }
  }
  
  return(selected_species)
}

# End
