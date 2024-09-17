
###########################################################################
### Taylor's power law graph
###########################################################################

# Import data -------------------------------------------------------------
tpl <- read.table(here::here("data", "processed", "slope_taylor.txt"), h = TRUE, sep = ";")
tree <- ape::read.nexus(here::here("data", "raw", "fish_phylo.nex"))

# Load functions ----------------------------------------------------------
library(ggtree)
library(ggplot2)
library(ggplotify)
library(aplot)
library(RColorBrewer)
library(ape)
library(dplyr)


# Graphics for tpl slope --------------------------------------------------
## Reorder positions based on phylogeny
tpl.g <- tpl
tpl.g[,"tip.label"] <- gsub(" ", "", tpl.g[, "ancient_scientific_name"])
tpl.g <- tpl.g[match(tree$tip.label, tpl.g[, "tip.label"]),]
tree$tip.label <- tpl.g[,"scientific_name"]
row.names(tpl.g) <- tpl.g[,"scientific_name"]

## Match the tree and the data
d.traits <- tpl.g[, c("scientific_name", "mean_regional_size", "spawning_type", "reproductive_guild", "origin")]
d.traits$spawning_type <- ifelse(d.traits$spawning_type == "Partial", "Parcial", "Total")
d.traits$origin <- ifelse(d.traits$origin == "Native", "Nativa", "Invasora")

tree_data <- fortify(tree) %>%
  left_join(d.traits, by = c("label" = "scientific_name"))

## Plotting tree
g <- ggtree(tree) %<+% d.traits +
  
  geom_tippoint(aes(size = mean_regional_size), shape = 1, x = 186) +
  scale_size_continuous(range = c(2, 8)) +
  
  geom_tippoint(aes(color = spawning_type), shape = 15, size = 5, x = 196) +
  scale_color_manual(values = c("Parcial" = "#ffb74d", "Total" = "#7fafff")) +
  
  geom_tippoint(aes(fill = origin), shape = 24, size = 4, x = 206) +
  scale_fill_manual(values = c("Invasora" = "#F44336", "Nativa" = "#4CAF50")) +
  
  geom_tippoint(aes(shape = reproductive_guild), size = 4, x = 216) +
  scale_shape_manual(values = c("LMEF" = 19, "NEFP" = 8, "NEFW" = 0, "NIF" = 2)) +
  
  coord_cartesian(xlim = c(0, 214), ylim = c(0, 51)) +
  
  annotate("text", x = 186, y = length(tree$tip.label) + 1, label = "TC", color = "black", 
           size = 5, angle = 90, hjust = 0, family = "serif") +
  annotate("text", x = 196, y = length(tree$tip.label) + 1, label = "TD", color = "black", 
           size = 5, angle = 90, hjust = 0, family = "serif") +
  annotate("text", x = 206, y = length(tree$tip.label) + 1, label = "OR", color = "black", 
           size = 5, angle = 90, hjust = 0, family = "serif") +
  annotate("text", x = 216, y = length(tree$tip.label) + 1, label = "GR", color = "black", 
           size = 5, angle = 90, hjust = 0, family = "serif") +
  
  guides(
    size = guide_legend(title = "Tamanho do corpo", title.theme = element_text(family = "serif", size = 16),
                        label.theme = element_text(family = "serif", size = 16)),
    shape = guide_legend(title = "Guilda reprodutiva", title.theme = element_text(family = "serif", size = 16),
                         label.theme = element_text(family = "serif", size = 16)),
    fill = guide_legend(title = "Origem", title.theme = element_text(family = "serif", size = 16),
                        label.theme = element_text(family = "serif", size = 16)),
    color = guide_legend(title = "Tipo de desova", title.theme = element_text(family = "serif", size = 16),
                         label.theme = element_text(family = "serif", size = 16))
  ) +
  
  theme(
    legend.text = element_text(family = "serif", size = 16),
    legend.title = element_text(family = "serif", size = 16)
  )

## Plotting tree with tpl slope
colors.p <- alpha(brewer.pal(5, "Set2"), 0.9)

p <- ggplot(tpl.g, aes(x = slope, y = scientific_name, fill = order)) +
  theme_tree2() +
  geom_bar(stat = "identity", alpha = 0.5) + 
  geom_text(aes(label = scientific_name, x = slope), 
            hjust = -0.1, vjust = 0.3,
            family = "serif", fontface = "italic", size = 6) +
  labs(x = "Coeficiente angular da TPL", y = NULL, fill = "Ordem") +
  scale_fill_manual(values = colors.p) +
  coord_cartesian(xlim = c(1, 2.7), ylim = c(0, 51)) +
  theme(
    text = element_text(family = "serif", size = 20),
    legend.position = "inside",
    legend.position.inside = c(0.12, 0.75),
    legend.text = element_text(size = 16),
    legend.title = element_text(family = "serif", size = 16),
    axis.ticks.y = element_line(color = "black"),
    axis.line.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.5),
    axis.ticks.length.y = unit(-0.1, "cm"),
    axis.text.x = element_text(color = "black"),
    plot.margin = unit(c(1, 4, 1, 1), "lines")
  ) +
  geom_segment(aes(x = 0.915, xend = 0.915, y = -0.6, yend = 50), color = "black", size = 0.5)

gt <- p %>% insert_left(g, width = 1)

print(gt)

# Export graphics taylor --------------------------------------------------
## pdf
ggsave(here::here("output", "figures", "slope_taylor.pdf"), plot = gt, width = 17, height = 12)

# End