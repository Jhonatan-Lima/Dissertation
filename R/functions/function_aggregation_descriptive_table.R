
###########################################################################
#### Aggregation Descriptive Table
###########################################################################

##
# index                 : List containing all indexes (data frames or matrices with species data)
# characteristics       : Data frame with species traits, where the first column contains species names
##

# Function ----------------------------------------------------------------
tab_descriptive <- function(index, characteristics) {
  
  ## Extract unique species names
  species <- unique(characteristics[, 1])
  
  ## Initialize result matrix
  num_indexes <- length(index)
  num_species <- length(species)
  aggregation <- as.data.frame(matrix(ncol = (4 * num_indexes), nrow = num_species))
  
  ## Generate column names for the result matrix based on the list names
  col_names <- vector("character", length = 4 * num_indexes)
  for (c in seq_len(num_indexes)) {
    list_name <- names(index)[c]
    col_names[(4 * (c - 1) + 1)] <- paste("min", list_name, sep = "_")
    col_names[(4 * (c - 1) + 2)] <- paste("max", list_name, sep = "_")
    col_names[(4 * (c - 1) + 3)] <- paste("med", list_name, sep = "_")
    col_names[(4 * (c - 1) + 4)] <- paste("sd", list_name, sep = "_")
  }
  colnames(aggregation) <- col_names
  
  ## Populate the result matrix with descriptive statistics
  for (i in seq_len(num_species)) {
    sp <- species[i]  # Define species
    sp_matrix <- characteristics[, 1] == sp  # Identify rows corresponding to the species
    sp_complete <- characteristics[sp_matrix, "scientific_name"]  # Get full species name
    rownames(aggregation)[i] <- sp_complete  # Set row name in the result matrix
    
    for (p in seq_len(num_indexes)) {
      presence <- colnames(index[[p]]) == sp
      if (any(presence)) {
        aggregation[i, paste("min", names(index)[p], sep = "_")] <- min(index[[p]][, sp], na.rm = TRUE)
        aggregation[i, paste("max", names(index)[p], sep = "_")] <- max(index[[p]][, sp], na.rm = TRUE)
        aggregation[i, paste("med", names(index)[p], sep = "_")] <- mean(index[[p]][, sp], na.rm = TRUE)
        aggregation[i, paste("sd", names(index)[p], sep = "_")] <- sd(index[[p]][, sp], na.rm = TRUE)
      }
    }
  }
  
  return(aggregation)
}

# End