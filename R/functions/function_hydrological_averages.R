
############################################################################
#### Create the Lag Matrix for Water Level
############################################################################

##
# lv.water  : Data frame with monthly water level data. 
#             Columns should include "ano.mes" for year-month and "media" for the water level.
##

# Function ------------------------------------------------------------------
matriz_lag <- function(lv.water) {
  
  # Load required packages
  library(tidyr)
  library(DataCombine)
  library(stringr)
  
  # Standardize "ano.mes" column
  agua <- separate(lv.water, ano.mes, into = c("ano", "mes"), sep = "/")
  agua$mes <- str_pad(agua$mes, width = 2, pad = "0")
  agua$ano.mes <- paste(agua$ano, agua$mes, sep = "/")
  
  # Select columns ("ano.mes" and "media")
  nv_mensal <- agua[, c("ano.mes", "media")]
  
  # Add new rows with missing months
  missing_rows <- data.frame(
    ano.mes = c("1999/09", "1999/10", "1999/11", "1999/12"),
    media = NA
  )
  nv_mensal <- rbind(missing_rows, nv_mensal)
  
  # Define collection months
  coletas <- c("2000/02", "2000/05", "2000/08", "2000/11", "2001/02", "2001/05",
               "2001/08", "2001/10", "2002/03", "2002/05", "2002/08", "2002/11",
               "2003/03", "2003/06", "2003/09", "2003/12", "2004/03", "2004/06",
               "2004/09", "2004/12", "2005/03", "2005/06", "2005/09", "2005/12",
               "2006/03", "2006/06", "2006/09", "2006/12", "2007/03", "2007/06",
               "2007/09", "2007/11", "2008/02", "2008/06", "2008/09", "2008/11",
               "2009/03", "2009/06", "2009/09", "2009/12", "2010/03", "2010/06",
               "2010/09", "2010/12", "2011/03", "2011/06", "2011/09", "2011/12",
               "2012/02", "2012/06", "2012/09", "2012/12", "2013/03", "2013/06",
               "2013/09", "2013/12", "2014/03", "2014/06", "2014/09", "2014/12",
               "2015/03", "2015/06", "2015/09", "2015/11", "2016/02", "2016/06",
               "2016/09", "2016/11", "2017/03", "2017/06", "2017/09", "2017/12")
  
  # Create the result matrix
  medias_agua <- as.data.frame(matrix(ncol = 7, nrow = length(coletas)))
  colnames(medias_agua) <- c("ano.mes", "nv.df0", "nv.df1", "nv.df2", "nv.df3", "nv.df4", "nv.df5")
  
  # Fill the result matrix with lagged values
  for (i in seq_along(coletas)) {
    
    mes_coleta <- coletas[i]                           # Define the collection month
    posicao <- which(nv_mensal$ano.mes == mes_coleta)  # Find the position of the month
    
    # Fill results
    medias_agua[i, "ano.mes"] <- mes_coleta
    medias_agua[i, "nv.df0"] <- nv_mensal[posicao, "media"]
    medias_agua[i, "nv.df1"] <- ifelse(posicao > 1, nv_mensal[posicao - 1, "media"], NA)
    medias_agua[i, "nv.df2"] <- ifelse(posicao > 2, nv_mensal[posicao - 2, "media"], NA)
    medias_agua[i, "nv.df3"] <- ifelse(posicao > 3, nv_mensal[posicao - 3, "media"], NA)
    medias_agua[i, "nv.df4"] <- ifelse(posicao > 4, nv_mensal[posicao - 4, "media"], NA)
    medias_agua[i, "nv.df5"] <- ifelse(posicao > 5, nv_mensal[posicao - 5, "media"], NA)
  }
  
  # Convert columns to numeric
  medias_agua[, 2:7] <- lapply(medias_agua[, 2:7], as.numeric)
  
  return(medias_agua)
}

# End

# Demonstration of how to run the function ------------------------------

# lv.water = data frame with monthly water level data including "ano.mes" and "media"

