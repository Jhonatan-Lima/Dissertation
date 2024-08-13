
############################################################################
#### Create the GLS Matrix of Models
############################################################################

##
# aggregation       : Aggregation matrix with species and temporal data
# proportion        : Matrix of immature proportion with species and temporal data
# species           : Vector of species names
# water_level       : Column name representing the water level of the Paraná River
##

# Function ----------------------------------------------------------------
corr_gls <- function(aggregation, proportion, species, water_level) {
  
  # Load species traits data
  traits <- read.csv(here::here("data", "raw", "species_traits.csv"), header = TRUE)
  
  library(nlme)
  
  # Initialize a data frame to store results
  tab_resu <- as.data.frame(matrix(ncol = 12, nrow = 0))
  
  # Loop through each species to fit GLS models
  for (s in seq_along(species)) {
    
    sp <- species[s]
    
    # Create a combined data frame for regression analysis
    data <- dplyr::left_join(
      aggregation[, c("ano.mes", sp, water_level)],
      proportion[, c(sp, "ano.mes")],
      by = "ano.mes"
    )
    
    # Create the 'months' column based on the row index
    data$months <- seq(nrow(data))
    
    # Rename columns for clarity
    colnames(data)[2:4] <- c("aggre", "level", "prop")
    
    # Log-transform the water level
    data$level <- log(data$level)
    
    # Remove na
    data <- na.omit(data)
    
    # Define the model formula
    fmla <- as.formula("log(aggre) ~ level + prop")
    
    # Fit the GLS model
    model <- nlme::gls(
      fmla,
      data = data,
      na.action = na.omit,
      correlation = nlme::corARMA(p = 1, q = 0, form = ~ months)
    )
    
    # Extract confidence intervals
    intervals <- nlme::intervals(model)
    
    # Extract model diagnostics
    slope_level <- intervals$coef["level", "est."]
    ci_lower_level <- intervals$coef["level", "lower"]
    ci_upper_level <- intervals$coef["level", "upper"]
    p_level <- coef(summary(model))["level", "p-value"]
    slope_prop <- intervals$coef["prop", "est."]
    ci_lower_prop <- intervals$coef["prop", "lower"]
    ci_upper_prop <- intervals$coef["prop", "upper"]
    p_prop <- coef(summary(model))["prop", "p-value"]
    
    # Perform Shapiro-Wilk test for normality of residuals
    shp <- shapiro.test(resid(model))
    p_shapiro <- shp$p.value
    
    # Calculate AIC for model selection
    aic <- AIC(model)
    
    # Compile results into a single row
    line <- c(sp, slope_level, ci_lower_level, ci_upper_level, p_level,
              slope_prop, ci_lower_prop, ci_upper_prop, p_prop, p_shapiro, aic, nobs(model))
    
    # Append the results to the results data frame
    tab_resu <- rbind(tab_resu, line)
  }
  
  # Set column names for the results data frame
  colnames(tab_resu) <- c("specie", "slope.level", "ci.lower.level", "ci.upper.level", 
                          "p.level", "slope.prop", "ci.lower.prop", "ci.upper.prop", 
                          "p.prop", "p.shapiro", "aic", "n")
  
  # Merge the results with species traits
  results <- dplyr::left_join(
    traits[, c("specie", "scientific_name", "short_name", "spawning_type", "migration", "diet", "parental_care")],
    tab_resu,
    by = "specie"
  )
  
  return(results)
}

# End