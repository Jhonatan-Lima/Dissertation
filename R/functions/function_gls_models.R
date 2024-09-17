
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
  
  # Initialize a data frame to store results with explicit column types
  tab_resu <- data.frame(
    specie = rep(NA_character_, length(species)),
    slope.level = rep(NA_real_, length(species)),
    ci.lower.level = rep(NA_real_, length(species)),
    ci.upper.level = rep(NA_real_, length(species)),
    p.level = rep(NA_real_, length(species)),
    slope.prop = rep(NA_real_, length(species)),
    ci.lower.prop = rep(NA_real_, length(species)),
    ci.upper.prop = rep(NA_real_, length(species)),
    p.prop = rep(NA_real_, length(species)),
    p.shapiro = rep(NA_real_, length(species)),
    aic = rep(NA_real_, length(species)),
    n = rep(NA_real_, length(species)),
    stringsAsFactors = FALSE
  )
  
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
    
    # Perform Shapiro-Wilk test for normality of residuals
    shp <- shapiro.test(resid(model))
    p_shapiro <- shp$p.value
    
    # Calculate AIC for model selection
    aic <- AIC(model)
    
    # Append the results to the results data frame
    tab_resu[s,"specie"] <- sp
    tab_resu[s,"slope.level"] <- intervals$coef["level", "est."]
    tab_resu[s,"ci.lower.level"] <- intervals$coef["level", "lower"]
    tab_resu[s,"ci.upper.level"] <- intervals$coef["level", "upper"]
    tab_resu[s,"p.level"] <- coef(summary(model))["level", "p-value"]
    tab_resu[s,"slope.prop"] <- intervals$coef["prop", "est."]
    tab_resu[s,"ci.lower.prop"] <- intervals$coef["prop", "lower"]
    tab_resu[s,"ci.upper.prop"] <- intervals$coef["prop", "upper"]
    tab_resu[s,"p.prop"] <- coef(summary(model))["prop", "p-value"]
    tab_resu[s,"p.shapiro"] <- p_shapiro
    tab_resu[s,"aic"] <- aic
    tab_resu[s,"n"] <- nobs(model)
  }
  
  # Merge the results with species traits
  results <- dplyr::left_join(
    traits[, c("specie", "scientific_name", "short_name", "spawning_type", "migration", "diet", "parental_care")],
    tab_resu,
    by = "specie"
  )
  
  return(results)
}

# End