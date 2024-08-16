
###########################################################################
### Script for PGLS Analysis Between TPL and Species Traits
###########################################################################

# Import data -------------------------------------------------------------
tpl <- read.table(here::here("data", "processed", "slope_taylor.txt"), header = TRUE, sep = ";")
tree <- ape::read.nexus(here::here("data", "raw", "fish_phylo.nex"))

# Load required packages -------------------------------------------------
library(ape)
library(car)
library(visreg)
library(sjPlot)
library(ggplot2)
library(ggplotify)

# Clean species names by removing spaces ---------------------------------
tpl[, "tip_label"] <- gsub(" ", "", tpl[, "ancient_scientific_name"])

# Check and reorder species based on phylogeny ----------------------------
# Compare and reorder the species in tpl to match the phylogenetic tree
tree$tip.label
tpl$tip_label

# Reorder tpl to match tree
reorder <- match(tree$tip.label, tpl[, "tip_label"])
tpl.reordered <- tpl[reorder, ]
tree$tip.label == tpl.reordered[,"tip_label"]

# Update tree tip labels with scientific names from tpl
tree$tip.label <- tpl.reordered[,"scientific_name"]

# Subset relevant columns and factorize categorical variables -------------
sub.tpl <- subset(tpl.reordered, select = c("slope", "mean_regional_size", "spawning_type", "reproductive_guild", "origin", "order"))
sub.tpl[,"spawning_type"] <- as.factor(sub.tpl[,"spawning_type"])
sub.tpl[,"reproductive_guild"] <- as.factor(sub.tpl[,"reproductive_guild"])
sub.tpl[,"origin"] <- as.factor(sub.tpl[,"origin"])
sub.tpl[,"order"] <- as.factor(sub.tpl[,"order"])

# Log-transform mean_regional_size ---------------------------------------
sub.tpl[,"mean_regional_size"] <- log(sub.tpl[,"mean_regional_size"])

# Fit PGLS models ---------------------------------------------------------
# Define formula
fmla <- as.formula("slope ~ mean_regional_size + spawning_type + reproductive_guild + origin")

# Fit models with different evolutionary correlations
pgls.mod <- nlme::gls(fmla, data = sub.tpl)
pgls.brownian <- nlme::gls(fmla, data = sub.tpl, correlation = corBrownian(1, phy = tree))
pgls.grafen <- nlme::gls(fmla, data = sub.tpl, correlation = corGrafen(1, phy = tree))
pgls.martins <- nlme::gls(fmla, data = sub.tpl, correlation = corMartins(1, phy = tree))
pgls.pagel <- nlme::gls(fmla, data = sub.tpl, correlation = corPagel(1, phy = tree))
pgls.blomberg <- nlme::gls(fmla, data = sub.tpl, correlation = corBlomberg(1, phy = tree, fixed = TRUE))

# Plot diagnostic plots for each model -----------------------------------
plot(pgls.mod)
plot(pgls.brownian)
plot(pgls.grafen)
plot(pgls.martins)
plot(pgls.pagel)
plot(pgls.blomberg)

# Compare models using AIC ------------------------------------------------
aic <- AIC(pgls.mod, pgls.brownian, pgls.grafen, pgls.martins, pgls.pagel, pgls.blomberg)
aic[,"dAIC"] <- aic[,"AIC"] - aic[1,"AIC"]

# Check the best models ---------------------------------------------------
summary(pgls.mod)
summary(pgls.pagel)

# Perform ANOVA on models
car::Anova(pgls.mod)
car::Anova(pgls.pagel)

# Create table for Pagel model results -------------------------------------
tab.pagel <- tab_model(pgls.mod, pgls.pagel)

# Partial regression plot for mean_regional_size ---------------------------
partial <- visreg(pgls.pagel, "mean_regional_size", partial = TRUE, band = TRUE, gg = TRUE)
visreg_data <- partial$data

size_plot <- ggplot(visreg_data, aes(x = x, y = y)) +
  geom_point(shape = 16, size = 2, color = "black", alpha = 0.5) +
  theme_bw() +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  theme(text = element_text(family = "serif", size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black", size = 13)
  ) +
  labs(y = "Coeficiente angular da TPL", x = "Log(Tamanho do corpo)")

print(size_plot)

# Export plots and tables -------------------------------------------------
## Save plot as PDF
ggsave(here::here("output", "figures", "pgls_taylor.pdf"), plot = size_plot, width = 8, height = 5)

## Export AIC results to text file
write.table(aic, file = here::here("output", "tables", "aic_pgls.txt"), sep = ";")

## Export Pagel model results to Word document
tab_model(pgls.pagel, file = here::here("output", "tables", "tab_pagel.doc"))

# End
