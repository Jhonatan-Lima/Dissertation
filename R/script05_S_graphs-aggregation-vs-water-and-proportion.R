
###########################################################################
#### Create species aggregation graphs
###########################################################################

# Import data -------------------------------------------------------------
gls.vm <- list(
  df0 = read.table(here::here("data","gls","vm_df0.txt"), h = T, sep = ";"),
  df1 = read.table(here::here("data","gls","vm_df1.txt"), h = T, sep = ";"),
  df2 = read.table(here::here("data","gls","vm_df2.txt"), h = T, sep = ";"),
  df3 = read.table(here::here("data","gls","vm_df3.txt"), h = T, sep = ";")
)
gls.im <- list(
  df0 = read.table(here::here("data","gls","im_df0.txt"), h = T, sep = ";"),
  df1 = read.table(here::here("data","gls","im_df1.txt"), h = T, sep = ";"),
  df2 = read.table(here::here("data","gls","im_df2.txt"), h = T, sep = ";"),
  df3 = read.table(here::here("data","gls","im_df3.txt"), h = T, sep = ";")
)
vm.obs <- read.table(here::here("data","processed","aggregation_vm.txt"), h = T, sep = ";")
im.obs <- read.table(here::here("data","processed","aggregation_im.txt"), h = T, sep = ";")
traits <- read.csv(here::here("data","raw","species_traits.csv"), h = T)

## Species
source(here::here("R","functions","function_species_selection.R"))
species <- species_sel(vm.obs)

## Load packages
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggplotify)
library(ggplot2)
library(gridExtra)

# Create function to generate the graph -----------------------------------

###########################################################################
#### Function generate the graph
###########################################################################

##
# index               : List with matrix with gls slope coefficient data
# type                : Defines whether the generated graph is of the hydrological level (level) or the proportion of immature species (prop)
# title.axis.x        : Defines whether you want the x-axis label to appear
# labels.axis.y       : Defines whether you want the species name to appear on the y-axis
# scale.x             : Defines the scale for the x-axis
##

## Obs: Migration (Yes = triangle, No = circle)
## black point represent significant species

graph <- function(index, type = "prop", title.axis.x = T, labels.axis.y = T, scale.x){
  
  library(ggplot2)
  
  ## Lag
  lag <- c("Lag0", "Lag1", "Lag2", "Lag3")
  
  ## Results
  graphic_list <- list()
  
  for(i in 1:length(index)){
    
    ## Which variable
    if(type == "level"){
      matriz <- index[[i]]
      matriz <- matriz[,c("scientific_name","short_name","slope.level","ci.lower.level","ci.upper.level","p.level","spawning_type","migration","diet")]
    }
    if(type == "prop"){
      matriz <- index[[i]]
      matriz <- matriz[,c("scientific_name","short_name","slope.prop","ci.lower.prop","ci.upper.prop","p.prop","spawning_type","migration","diet")]
    }
    
    colnames(matriz) <- c("scientific_name","short_name","slope","ci.lower","ci.upper","p","spawning_type","migration","diet")
    matriz[,"sig"] <- ifelse(matriz[,"p"] < 0.05, "Yes", "No")
    matriz$short_name <- factor(matriz$short_name, levels = matriz$short_name)
    
    
    ## Graphic
    graphic <- ggplot() +
      theme_bw() + 
      geom_vline(xintercept = 0, color = "red", linetype = "longdash", linewidth = 0.5) +
      geom_errorbarh(data = matriz, aes(y = short_name,, xmax = ci.upper, xmin = ci.lower), linetype = "solid", color = "black") +
      geom_point(
        data = matriz,
        size = 3,
        aes(
          x = slope,
          y = short_name,,
          shape = as.factor(migration),
          fill = as.factor(sig)
        )) +
      scale_x_continuous(
        limits = c(scale.x[1],scale.x[2]),
        breaks = seq(scale.x[1],scale.x[2], by = scale.x[3])
      ) +
      theme(
        text = element_text(size = 22, color = "black", family = "serif"),
        plot.title = element_text(hjust = 0.05, vjust = -7, size = 12),
        legend.position = "none",
        plot.margin = margin(5, 20, 5, 20)
      ) +
      scale_shape_manual(values = c("Yes" = 24, "No" = 21)) +
      scale_fill_manual(values = c("Yes" = "black", "No" = "white"))
    
    
    ## Control y-axis labels
    if(labels.axis.y == T) {
      graphic <- graphic + 
        theme(
          axis.text.y = element_text(color = "black", face = "italic")
        ) +
        labs(x = "", y = "")
    }
    if(labels.axis.y == F) {
      graphic <- graphic + 
        theme(
          axis.text.y = element_blank()
        ) +
        labs(x = "", y = "")
    }
    
    ## Lag
    l <- lag[i]
    
    ## Control x-axis title
    if(title.axis.x == T) {
      graphic <- graphic + 
        theme(
          axis.text.x = element_text(family = "serif", size = 22, color = "black"),
          axis.title.x = element_text(family = "serif", size = 22),
          plot.title = element_text(family = "serif", size = 22, face = "bold")
        ) +
        labs(x = "Slope", y = "", title = l)
    }
    if(title.axis.x == F) {
      graphic <- graphic +
        theme(plot.title = element_text(family = "serif", size = 22, face = "bold")) +
        labs(x = "", y = "", title = l)
    }
    
    graphic_list[[i]] <- graphic
  }
  
  # Return graphic
  names(graphic_list) <- lag
  return(graphic_list)
  
}

g.im.lv <- graph(index = gls.im, type = "level", title.axis.x = T, labels.axis.y = T, scale.x = c(-2,2,0.5))
g.vm.lv <- graph(index = gls.vm, type = "level", title.axis.x = T, labels.axis.y = T, scale.x = c(-4,5,1))


# And


# Scatter plot between aggregation and water level -------------------------
## Reorganize data
im.long <- im.obs[,c(species,"ano.mes")] %>%
  pivot_longer(
    cols = -ano.mes,               
    names_to = "specie",        
    values_to = "im"          
  )

vm.long <- vm.obs[,c(species,"ano.mes")] %>%
  pivot_longer(
    cols = -ano.mes,               
    names_to = "specie",        
    values_to = "vm"          
  )

aggregation <- im.long %>%
  left_join(vm.long, by = c("ano.mes", "specie"))

aggregation <- aggregation %>%
  left_join(im.obs[,c("ano.mes","nv.df0","nv.df1","nv.df2","nv.df3")], by = c("ano.mes"))

aggregation <- aggregation %>%
  left_join(traits[,c("specie","family")], by = c("specie"))

aggregation <- na.omit(aggregation)

## Scatter plot
sc.plot <- function(data, index, lag) {
  
  if(lag == "nv.df0") { df <- "Lag0" }
  if(lag == "nv.df1") { df <- "Lag1" }
  if(lag == "nv.df2") { df <- "Lag2" }
  if(lag == "nv.df3") { df <- "Lag3" }
  
  if(index == "im") { id <- "Índice de Morisita" }
  if(index == "vm") { id <- "Variância/Média" }
  
  cor <- c(brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 8, name = "Dark2"))
  lines <- rep("solid", length(unique(data$specie)))
  
  plot <- ggplot(data = data, aes(x = log(get(lag)), y = log(get(index)))) +
    geom_point(size = 2, color = "gray") +
    geom_smooth(aes(colour = family, linetype = specie), method = "lm", se = FALSE) +
    scale_color_manual(values = cor) +
    scale_linetype_manual(values = lines) +
    labs(
      x = paste("Log(Nível hidrológico - ", df, ")", sep = ""), 
      y = paste("Log(", id, ")", sep = ""), 
      color = "Família"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 22, color = "black", family = "serif"),
      plot.title = element_text(hjust = 0.5, size = 22),
      axis.text.y = element_text(color = "black", size = 22),
      axis.text.x = element_text(color = "black", size = 22),
      legend.position = "bottom",
      legend.title = element_text(
        hjust = 0.5,
        vjust = 0.5,  
        size = 20
      ),
      plot.margin = margin(30, 20, 20, 20),
      legend.box.margin = margin(t = 20)
    ) +
    guides(
      linetype = guide_none(),
      color = guide_legend(
        title.position = "top",  
        title.hjust = 0.5,
        keyheight = unit(1, "cm"), 
        keywidth = unit(1.2, "cm"),
        ncol = 3,
        byrow = TRUE
      )
    )
  
  return(plot)
}

sc.im.df0 <- sc.plot(data = aggregation, index = "im", lag = "nv.df0")
sc.vm.df0 <- sc.plot(data = aggregation, index = "vm", lag = "nv.df0")



# Join graphs -------------------------------------------------------------
im.df0 <- grid.arrange(g.im.lv[[1]], sc.im.df0, ncol = 2, widths = c(8,7))
vm.df0 <- grid.arrange(g.vm.lv[[1]], sc.vm.df0, ncol = 2, widths = c(8,7))


# Export graph ------------------------------------------------------------
ggsave(here::here("output","figures","df0_morisita.pdf"), plot = im.df0, width = 17, height = 13)
ggsave(here::here("output","figures","df0_var_med.pdf"), plot = vm.df0, width = 17, height = 13)
