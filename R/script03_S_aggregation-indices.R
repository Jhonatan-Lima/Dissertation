
###########################################################################
#### Graphs of aggregation indices
###########################################################################

# Import data -------------------------------------------------------------
tab.aggregation <- read.table(here::here("output","tables","aggregation_descriptive_table.txt"), h = TRUE, sep = ";")
im.obs <- read.table(here::here("data","processed","aggregation_im.txt"), h = TRUE, sep = ";")
vm.obs <- read.table(here::here("data","processed","aggregation_vm.txt"), h = TRUE, sep = ";")

# Load packages -----------------------------------------------------------
library(ggplot2)
library(ggplotify)
library(gridExtra)
library(RColorBrewer)
library(tidyr)


# Bar chart ---------------------------------------------------------------
index.tab <- tab.aggregation
index.tab$scientific_name <- factor(index.tab$scientific_name, levels = index.tab$scientific_name)
colors.g <- alpha(brewer.pal(5, "Set2"), 0.9)

im <- ggplot(index.tab, aes(x = med_im, y = scientific_name, fill = order)) +
  theme_bw() +
  geom_bar(stat = "identity", alpha = 0.5, size = 0.5, width = 0.8, linewidth = 0.3, color = "black") +
  labs(x = "Índice de Morisita", y = NULL) +
  theme(
    text = element_text(family = "serif", size = 22, color = "black"),
    legend.position = "none",
    axis.text.x = element_text(color = "black"),  
    axis.text.y = element_text(color = "black", face = "italic", vjust = 0.5)
  ) +
  scale_fill_manual(values = colors.g) +
  coord_cartesian(xlim = c(0.2, 6))

vm <- ggplot(index.tab, aes(x = med_vm, y = scientific_name, fill = order)) +
  theme_bw() +
  geom_bar(stat = "identity", alpha = 0.5, size = 0.5, width = 0.8, linewidth = 0.3, color = "black") +
  labs(x = "Razão variância/média", y = NULL) +
  theme(
    text = element_text(family = "serif", size = 22, color = "black"),
    legend.position = "right",
    legend.spacing.y = unit(2, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1, "cm"), 
    axis.text = element_text(color = "black"),
    axis.text.y = element_blank()
  ) +
  scale_fill_manual(values = colors.g) +
  coord_cartesian(xlim = c(5, 160)) +
  guides(color = "none",
         fill = guide_legend(
          title = "Ordem"))

bar_chart <- grid.arrange(im, vm, ncol = 2, widths = c(8,7))


# Aggregation index time series graph -------------------------------------
## Function
g.temp <- function(index, specie, scientific_name, index_name){
  
  index <- tidyr::separate(index, ano.mes, into = c("ano", "mes"), sep = "/")
  index[,"ano"] <- as.numeric(index[,"ano"])
  index[,"ano.mes"] <- index[,"ano"] + c(0,0.25,0.5,0.75)
  
  plot <- ggplot(data = index) +
    theme_bw() +
    geom_line(aes(y = get(specie), x = ano.mes), color = "black", linetype = "dotted", na.rm = TRUE) +
    geom_point(aes(y = get(specie), x = ano.mes), size = 2, color = "black", alpha = 0.8, na.rm = TRUE) +
    scale_x_continuous(limits = c(2000, 2018), breaks = seq(2000, 2018, by = 3)) +
    scale_y_continuous(limits = c(1, 4.5), breaks = seq(1, 5, by = 1)) +
    geom_hline(yintercept = 1, color = "red", linetype = "longdash", linewidth = 0.5) +
    theme(
      text = element_text(family = "serif", size = 22, color = "black"),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black")
    ) +
    labs(y = index_name, x = "Ano de coleta") +
    annotate("text", x = -Inf, y = Inf, label = scientific_name, hjust = -0.1, vjust = 2, size = 7, family = "serif", fontface = "italic")
  
  return(plot)
}

l.platyme <- g.temp(index = im.obs, specie = "L.platyme", scientific_name = "Loricariichthys platymetopon", index_name = "Índice de Morisita")

# Export graphics aggregation ---------------------------------------------
## pdf
ggsave(here::here("output","figures","mean_aggregation_bar_chart.pdf"), plot = bar_chart, width = 17, height = 14)
ggsave(here::here("output","figures","time_series_L.platyme.pdf"), plot = l.platyme, width = 6, height = 4)
