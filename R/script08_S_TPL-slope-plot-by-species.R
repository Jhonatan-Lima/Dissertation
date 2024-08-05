
###########################################################################
#### TPL Slope Plot by Species
###########################################################################

# Import data -------------------------------------------------------------
cpue.obs <- read.table(here::here("data", "processed", "cpue_obs.txt"), header = TRUE, sep = ";")

# Load required packages and functions -------------------------------------
library(ggplotify)
library(cowplot)
source(here::here("R", "functions", "function_species_selection.R"))

# Define selected species -------------------------------------------------
species <- species_sel(abundance = cpue.obs)
species <- species[-which(species == "Hopli.spp")]

# Function to generate TPL slope plots -----------------------------------
# Parameters:
# data                : Matrix of species abundance data
# species             : Vector of species names
# season              : Season of the year
# labels.axis         : Boolean to determine if axis labels should be included
#
# Returns:
# A list of ggplot objects, each representing the TPL slope plot for a species
graph.tpl <- function(data, species, season, labels.axis = FALSE){
  
  # Load species traits
  traits <- read.csv(here::here("data", "raw", "species_traits.csv"), header = TRUE)
  
  library(ggplot2)
  library(dplyr)
  
  # Unique seasons and dimensions
  u_season <- unique(season)
  n_lagos <- length(season) / length(u_season)
  n_col <- length(species)
  n_row <- length(u_season)
  
  # Prepare data frame to store log-transformed mean and variance
  logs <- as.data.frame(matrix(nrow = n_row, ncol = 2))
  colnames(logs) <- c("log.med", "log.var")
  rownames(logs) <- u_season
  
  # Initialize list to store plots
  graphic_list <- vector("list", length = 49)
  names(graphic_list) <- species
  
  # Generate plots for each species
  for(c in 1:n_col){
    
    p.especie <- species[c]
    
    # Calculate log-transformed mean and variance for each season
    for(l in 1:n_row){
      
      p <- which(season == u_season[l])
      
      mean <- mean(data[p, p.especie], na.rm = TRUE)
      variance <- var(data[p, p.especie], na.rm = TRUE)
      
      log.mean <- ifelse(mean > 0, log10(mean), NA)
      log.variance <- ifelse(variance > 0, log10(variance), NA)
      
      logs[l, "log.med"] <- log.mean        
      logs[l, "log.var"] <- log.variance    
    }
    
    # Fit linear model to log-transformed data
    ols <- lm(data = logs, log.var ~ log.med, na.action = na.exclude)
    s.model <- summary(ols)
    
    # Extract model coefficients
    sp.name <- traits[which(traits[,"specie"] == p.especie), "scientific_name"]
    b.tpl <- round(s.model$coefficients["log.med", "Estimate"], 2)
    a.tpl <- round(s.model$coefficients["(Intercept)", "Estimate"], 2)
    logs.na <- na.omit(logs)
    
    # Create ggplot for the species
    tpl <- ggplot(logs.na, aes(x = log.med, y = log.var)) +
      geom_point(shape = 16, size = 4, color = "black", alpha = 0.5) +
      theme_classic() +
      geom_smooth(method = "lm", formula = 'y ~ x', color = "red", se = FALSE, fullrange = TRUE) +
      theme(text = element_text(family = "serif", size = 24),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_line(color = "black", linewidth = 1),
            axis.ticks.x = element_line(color = "black", linewidth = 1), 
            axis.text.x = element_text(color = "black", size = 24),
            axis.text.y = element_text(color = "black", size = 24),
            axis.line = element_line(color = "black", linewidth = 1)
      ) +
      labs(y = NULL, x = NULL) +
      annotate("text", x = -Inf, y = Inf, label = sp.name, hjust = -0.1, vjust = 2, size = 8, family = "serif", fontface = "italic") +
      annotate("text", x = -Inf, y = Inf - 0.1, label = paste("a = ", a.tpl, sep = ""), hjust = -0.3, vjust = 5, size = 8, family = "serif", fontface = "italic") +
      annotate("text", x = -Inf, y = Inf - 0.2, label = paste("b = ", b.tpl, sep = ""), hjust = -0.3, vjust = 6.5, size = 8, family = "serif", fontface = "italic") +
      coord_cartesian(xlim = c(-1, 2.5), ylim = c(-1, 5))
    
    # Add axis labels if required
    if(labels.axis) {
      tpl <- tpl + 
        labs(y = "Log(variância)", x = "Log(média)")
    } else {
      tpl <- tpl + 
        labs(x = NULL, y = NULL)
    }
    
    # Store the plot in the list
    graphic_list[[c]] <- tpl
  }
  
  return(graphic_list)
}

# Generate TPL slope plots for selected species
g <- graph.tpl(data = cpue.obs, species = species, season = cpue.obs$ano.mes, labels.axis = FALSE)

## Plot all graphs
for(i in 1:49){
  print(g[[i]])
}

## Combine plots into grid layouts
sp18 <- plot_grid(
  plot_grid(g[[1]], g[[2]], g[[3]], ncol = 3),
  plot_grid(g[[4]], g[[5]], g[[6]], ncol = 3),
  plot_grid(g[[7]], g[[8]], g[[9]], ncol = 3),
  plot_grid(g[[10]], g[[11]], g[[12]], ncol = 3),
  plot_grid(g[[13]], g[[14]], g[[15]], ncol = 3),
  plot_grid(g[[16]], g[[17]], g[[18]], ncol = 3),
  nrow = 6
)

sp36 <- plot_grid(
  plot_grid(g[[19]], g[[20]], g[[21]], ncol = 3),
  plot_grid(g[[22]], g[[23]], g[[24]], ncol = 3),
  plot_grid(g[[25]], g[[26]], g[[27]], ncol = 3),
  plot_grid(g[[28]], g[[29]], g[[30]], ncol = 3),
  plot_grid(g[[31]], g[[32]], g[[33]], ncol = 3),
  plot_grid(g[[34]], g[[35]], g[[36]], ncol = 3),
  nrow = 6
)

sp49 <- plot_grid(
  plot_grid(g[[37]], g[[38]], g[[39]], ncol = 3),
  plot_grid(g[[40]], g[[41]], g[[42]], ncol = 3),
  plot_grid(g[[43]], g[[44]], g[[45]], ncol = 3),
  plot_grid(g[[46]], g[[47]], g[[48]], ncol = 3),
  plot_grid(g[[49]], ncol = 3),
  nrow = 5
)

# Export combined plots to files -----------------------------------------
## Save plots as PDF
ggsave(here::here("output", "figures", "slope_sp18.pdf"), plot = sp18, width = 19, height = 27)
ggsave(here::here("output", "figures", "slope_sp36.pdf"), plot = sp36, width = 19, height = 27)
ggsave(here::here("output", "figures", "slope_sp49.pdf"), plot = sp49, width = 19, height = 22.5)

# End
