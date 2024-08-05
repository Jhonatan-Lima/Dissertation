
###########################################################################
### Function to Add Missing Lines to a Data Matrix
###########################################################################

##
# data      : Matrix where new rows must be inserted
# columns   : Vector with column names to preserve
# sites     : Vector with names of lakes
# years     : Vector with years where rows are missing
# months    : Vector with missing months in the matrix
##

# Function ----------------------------------------------------------------
add_lines <- function(data, columns, sites, years, months) {
  
  ## Filter the matrix to keep only the relevant columns
  data <- data[, c(columns, "locais", "anos", "meses")]
  
  ## Loop through each site, year, and month to add missing rows
  for (i in seq_along(sites)) {
    lake <- sites[i]  # Define the name of the lake
    
    for (a in seq_along(years)) {
      year <- years[a]  # Define the year
      
      for (m in seq_along(months)) {
        month <- months[m]  # Define the month
        
        ### Find the index for the specified lake, year, and month
        p <- which(data[,"locais"] == lake & data[,"anos"] == year & data[,"meses"] == month)
        
        line <- rep(NA, ncol(data))  # Create a row with NA values
        line[length(columns) + 1] <- lake  # Add lake name
        line[length(columns) + 2] <- year  # Add year
        line[length(columns) + 3] <- month + 3  # Add a placeholder for month
        
        # Split the matrix into upper and lower parts
        sup <- data[1:p,]  # Upper part
        inf <- data[(p+1):nrow(data),]  # Lower part
        
        ### Insert the new line at the appropriate position
        if (p == nrow(data)) {
          data <- rbind(sup, line)  # Append at the end
        } else {
          data <- rbind(sup, line, inf)  # Insert in the middle
        }
      }
    }
  }
  
  ## Adjust "anos" and "meses" columns to the desired format
  data[,"anos"] <- as.numeric(data[,"anos"]) + 2000
  data[,"meses"] <- stringr::str_pad(data[,"meses"], width = 2, pad = "0")
  data[,"ano.mes"] <- paste(data[,"anos"], data[,"meses"], sep = "/")
  data[, columns] <- apply(data[, columns], 2, as.numeric)
  
  return(data)
}

# End