
################################################################################
## Compare AIC of Models
################################################################################

##
# aggregation       : Aggregation matrix
# proportion        : Matrix of immature proportion
# species           : Vector of species names
# water_level       : Water level of Paraná River
##

# Function ----------------------------------------------------------------
compare_aic <- function(aggregation, proportion, species, water_level) {
  
  tab_res <- as.data.frame(matrix(ncol = 4, nrow = length(species)))  # Create results matrix
  
  ## Loop through species
  for (s in seq_along(species)) {
    sp <- species[s]  # Species name
    
    ## Create data matrix for regression
    data <- dplyr::left_join(aggregation[, c("ano.mes", sp, water_level)], 
                             proportion[, c(sp, "ano.mes")], by = "ano.mes")
    
    ## Rename columns
    data$tempo <- 1:nrow(data)
    colnames(data)[2] <- "aggre"
    colnames(data)[3] <- "level"
    colnames(data)[4] <- "prop"
    
    ## Models
    formula <- as.formula("log(aggre) ~ log(level) + prop")
    
    model_0 <- nlme::gls(formula, data = data, na.action = na.omit)
    model_1 <- nlme::gls(formula, data = data, na.action = na.omit,
                         correlation = nlme::corARMA(p = 1, q = 0, form = ~ tempo))
    model_2 <- nlme::gls(formula, data = data, na.action = na.omit,
                         correlation = nlme::corARMA(p = 2, q = 0, form = ~ tempo))
    
    ## Number of observations
    n_obs <- nobs(model_0)[1]
    
    ## AIC values
    aic_0 <- AIC(model_0)
    aic_1 <- AIC(model_1)
    aic_2 <- AIC(model_2)
    
    ## Add results to the matrix
    tab_res[s, ] <- c(aic_0, aic_1, aic_2, n_obs)
    rownames(tab_res)[s] <- sp
  }
  
  ## Rename columns
  colnames(tab_res) <- c("ar0", "ar1", "ar2", "n")
  
  ## Compare AIC values
  min_aic <- apply(tab_res[, 1:3], 1, min)  # Minimum AIC for each species
  comp_aic <- tab_res[, 1:3]
  comp_aic[] <- NA  # Initialize comparison matrix
  
  ## Determine best model considering AIC differences
  for (i in seq_len(nrow(tab_res))) {
    for (j in 1:3) {
      if ((min_aic[i] + 2) > tab_res[i, j]) {
        comp_aic[i, j] <- 1
      } else {
        comp_aic[i, j] <- 0
      }
    }
  }
  
  ## Delta AIC
  tab_aic <- tab_res[, 1:3] - min_aic
  tab_aic <- as.data.frame(tab_aic)
  tab_aic$species <- rownames(tab_aic)
  tab_aic <- tab_aic[, c("species", "ar0", "ar1", "ar2")]
  
  ## Print the number of times each model was selected as best
  print(colSums(comp_aic, na.rm = TRUE))
  
  return(tab_aic)
}

# End
