

###########################################################################
#### Create matrices of regression models
###########################################################################

# Import data -------------------------------------------------------------
vm.obs <- read.table(here::here("data","processed","aggregation_vm.txt"), h = T, sep = ";")
im.obs <- read.table(here::here("data","processed","aggregation_im.txt"), h = T, sep = ";")
immature <- read.table(here::here("data","processed","prop_immature.txt"), h = T, sep = ";")


# Load functions ----------------------------------------------------------
library(dplyr)
source(here::here("R","functions","function_selection_models_gls.R"))
source(here::here("R","functions","function_gls_models.R"))
source(here::here("R","functions","function_species_selection.R"))

## Select species
species <- species_sel(vm.obs)

## VIF between water level and proportion
vif <- as.data.frame(matrix(ncol = 2, nrow = length(species)))
colnames(vif) <- c("species", "vif")

for(i in 1:length(species)){
  
  ## specie
  sp <- species[i]
  
  ## data
  data <- dplyr::left_join(
    vm.obs[, c("ano.mes", sp, "nv.df0")],
    immature[, c(sp, "ano.mes")],
    by = "ano.mes"
  )
  colnames(data)[2:4] <- c("aggre", "level", "prop")
  data$months <- seq(nrow(data))
  
  ## model
  model <- nlme::gls(
    log(aggre) ~ level + prop,
    data = data,
    na.action = na.omit,
    correlation = nlme::corARMA(p = 1, q = 0, form = ~ months)
  )
  
  ## vif
  vif[i,"species"] <- sp
  vif[i,"vif"] <- car::vif(model)[1]
}


# Aggregation for variance/mean -------------------------------------------
# AIC models --------------------------------------------------------------
vm.df0 <- compare_aic(vm.obs, immature, species, water_level = "nv.df0")
vm.df1 <- compare_aic(vm.obs, immature, species, water_level = "nv.df1")
vm.df2 <- compare_aic(vm.obs, immature, species, water_level = "nv.df2")
vm.df3 <- compare_aic(vm.obs, immature, species, water_level = "nv.df3")

tab_aic_vm <- vm.df0 %>%
  left_join(vm.df1, by = "species", suffix = c("", ".df1")) %>%
  left_join(vm.df2, by = "species", suffix = c("", ".df2")) %>%
  left_join(vm.df3, by = "species", suffix = c("", ".df3"))


# GLS matrix --------------------------------------------------------------
gls.vm <- list(
  vm_df0 = corr_gls(vm.obs, immature, species, water_level = "nv.df0"),
  vm_df1 = corr_gls(vm.obs, immature, species, water_level = "nv.df1"),
  vm_df2 = corr_gls(vm.obs, immature, species, water_level = "nv.df2"),
  vm_df3 = corr_gls(vm.obs, immature, species, water_level = "nv.df3"))


# Export gls (vm x water) -------------------------------------------------
for(e in 1:length(gls.vm)){
  name <- paste(names(gls.vm)[e],".txt", sep = "")
  write.table(gls.vm[[e]], file = here::here("data","gls", name), sep = ";")
}


## And



# Aggregation for Morisita index ------------------------------------------
# AIC models --------------------------------------------------------------
im.df0 <- compare_aic(im.obs, immature, species, water_level = "nv.df0")
im.df1 <- compare_aic(im.obs, immature, species, water_level = "nv.df1")
im.df2 <- compare_aic(im.obs, immature, species, water_level = "nv.df2")
im.df3 <- compare_aic(im.obs, immature, species, water_level = "nv.df3")

tab_aic_im <- im.df0 %>%
  left_join(im.df1, by = "species", suffix = c("", ".df1")) %>%
  left_join(im.df2, by = "species", suffix = c("", ".df2")) %>%
  left_join(im.df3, by = "species", suffix = c("", ".df3"))

# GLS matrix --------------------------------------------------------------
gls.im <- list(
  im_df0 = corr_gls(im.obs, immature, species, water_level = "nv.df0"),
  im_df1 = corr_gls(im.obs, immature, species, water_level = "nv.df1"),
  im_df2 = corr_gls(im.obs, immature, species, water_level = "nv.df2"),
  im_df3 = corr_gls(im.obs, immature, species, water_level = "nv.df3"))


# Export gls (im x water) -------------------------------------------------
for(e in 1:length(gls.im)){
  name <- paste(names(gls.im)[e],".txt", sep = "")
  write.table(gls.im[[e]], file = here::here("data","gls", name), sep = ";")
}


# Export aic gls ----------------------------------------------------------
write.table(tab_aic_im, file = here::here("output","tables","aic_model_im.txt"), sep = ";")
write.table(tab_aic_vm, file = here::here("output","tables","aic_model_vm.txt"), sep = ";")


## End