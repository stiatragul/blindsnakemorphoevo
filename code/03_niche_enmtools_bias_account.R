# filename: 05_niche_enmtools_bias_account.R
# Sarin Tiatragul
# March 2023; mod Oct 2023

# Packages ----------------------------------------------------------------
options(java.parameters = "-Xmx8g") # for MaxEnt java call
# library(devtools)
# install_github("danlwarren/ENMTools") ; install_github('luismurao/ntbox')
library(ENMTools); library(ntbox); library(rJava)
library(terra); library(raster); 
library(dplyr); library(tidyr)
library(ggplot2); library(ape); library(phytools)
source('code/utility/func_bias.sample.R')
'%notin%' <- Negate('%in%')

# Data --------------------------------------------------------------------

# Read in data locality data
load('data/script_generated_data/soil_bulk_density_data.rda')

#####
##### niche quantification 
#####

# Prepared this raster in code/utility/01_soil_combine_rasters.R already changed the extent to match study area
# Biolayers2 uses merged resampled bio1 annual mean temp, SLGA bulk density and global aridity index (Aridity)
# Aridity index needs to be re-scaled to make sense
Biolayers2 <- terra::rast("./data/worldclim2_30s/merged_bio1_bdod_gai.tif")
names(Biolayers2) <- c("MeanTemp", "BulkDensity", "Aridity")
Biolayers2$Aridity <- Biolayers2$Aridity * 0.0001; Biolayers2 <- setMinMax(Biolayers2)

# Check that climate layers are all there and named correctly
# plot(Biolayers2)

# Check environmental data to make sure NA in an environmental layer propagates to the other layers
Biolayers2 <- check.env(Biolayers2)

# Make a dataframe of occurrence for Anilios
snakes <- soil_df %>% 
  dplyr::filter(-43 < latdec & latdec < -6.3) %>% # 
  dplyr::filter(species %notin% c("sp.", "braminus", "", "polygrammicus", "schlegelii", "kunuaensis", "infralabialis", 
                                  "curtus", "cumingii", "angusticeps", "bibroni", "exocoeti", "fornasinii",
                                  "willeyi", "solomonis", "depressus")) %>% 
  dplyr::filter(!is.na(BIO1)) %>% dplyr::filter(!is.na(BD)) %>% dplyr::filter(!is.na(ARID))

# How many species in total in our study ?
n_sample <- snakes %>% dplyr::group_by(species) %>% summarise(n=n())
n_sample %>% arrange(desc(n)); sum(n_sample$n)

# Write full and filtered occurrence to supplement table
supp_snake_occ <- snakes %>% dplyr::select(catalogue, genus, species, latdec, longdec, institution, extra_info)

# write.csv(snakes, file = "manuscript_word/supplements/supp_table_filtered_occ.csv", row.names = F, na = '')

supp_snake_occ %>% dplyr::group_by(species) %>% 
  count() %>% arrange(n) %>%  head(20)
nrow(supp_snake_occ)
### Input phylogeny

# Subset for only species that have >= 5 occurrences to use all available species
species_count <- table(snakes$species)
subset_snakes <- snakes[snakes$species %in% names(species_count[species_count >= 5]),]
subset_snakes_sp <- unique(subset_snakes$species)

# Load phylogeny
load('data/script_generated_data/subset_mcmctree_shape.rda')

# Subset phylogeny
sub_phy$tip.label <- gsub(x = sub_phy$tip.label, pattern = "Anilios_", replacement = "")
keepers <- sub_phy$tip.label[sub_phy$tip.label %in% subset_snakes_sp]
tr_anilios <- ape::keep.tip(phy = sub_phy, tip  = keepers)
write.tree(tr_anilios, file = "tr_anilios.tree")
#####
##### Check correlation matrix to see if we can reduce collinearity
#####

raster.cor.matrix(Biolayers2)

#####
##### Create enmtools.species object with loop
#####

# Make occurrence density grid for each species
# Grid for the niches for each species that have 
# Species need at least 5 points to estimate niche

# Get vector of species that have at least 5 points
n_sample %>% dplyr::filter(n >= 5) %>% unique() %>% arrange(n) 
sp_niche_calc <- n_sample %>% dplyr::filter(n >= 5) %>% unique() %>% arrange(species) 

### 
### Create sample bias layer
###
# using biaslayer() from package ntbox

## Occurrences
locations <- snakes[, c("species", "longdec", "latdec"),]
locations <- locations[complete.cases(locations),]

snake_loc_df <- snakes[, c("catalogue", "genus", "species", "longdec", "latdec"),]
# write.csv(snake_loc_df, file = "data/script_generated_data/blindsnake_locality_for_AS.csv", row.names = FALSE)

biasBlindsnakes <- ntbox::biaslayer(occs_df = locations, 
                             longitude = "longdec",
                             latitude =  "latdec",
                             raster_mold = raster(Biolayers2))

save(biasBlindsnakes, file = "data/script_generated_data/enmtools/bias_layer.rda")
load("data/script_generated_data/enmtools/bias_layer.rda")
biasBlindsnakes

# List of species object --------------------------------------------------
# create a list of species names
snake_species <- sp_niche_calc$species

# create an empty list to store the output
sp_list <- list()

# loop through the species names to build emntools.species object
for (i in seq_along(snake_species)) {
  
  # Get locality data
  xy <- snakes[which(snakes$species == snake_species[i]), c("longdec", "latdec"),]
  xy <- xy[complete.cases(xy),]
  xy <- vect(xy, geom = c("longdec", "latdec"), crs = terra::crs(Biolayers2))
  
  sp.obj <- enmtools.species()
  sp.obj$species.name <- snake_species[i]
  sp.obj$presence.points <- xy
  
  # Include buffer for range 
  sp.obj$range <- background.buffer(points = sp.obj$presence.points, 
                                    buffer.type = "circles",
                                    buffer.width = 50000,  # 50km radius          
                                    n = 1000,
                                    mask = Biolayers2[[1]],
                                    return.type = "raster")
  
  # and background points manually account for sampling bias
  sp.obj$background.points <- terra::vect(bias.sample(biasBlindsnakes, npoints = 1000) , geom= c("x","y"), crs = "EPSG:4326")
  
  # Check the species object are fine
  sp.obj <- check.species(sp.obj, env = Biolayers2)
  
  # append the sp.obj object to the sp_list
  sp_list[[i]] <- sp.obj
}

# Name items in the list of enmtools.species object
names(sp_list) <- snake_species



### 
### Build enmtools.clade object
###

### Build species.clade object
anilios.clade.bias <- enmtools.clade(species = sp_list[c(tr_anilios$tip.label)], tree = tr_anilios)
anilios.clade.bias <- check.clade(anilios.clade.bias)

#### 
#### Age-overlap correlation test (AOC)
####

## GEOGRAPHICAL RANGE OVERLAP

### Range overlap (Geographical range overlap)
range.aoc.bias <- enmtools.aoc(clade = anilios.clade.bias, nreps = 100, overlap.source = "range")
range.aoc.bias$coefficients
save(range.aoc.bias, file = 'data/script_generated_data/range_aoc_bias.Rdata')

# Point overlap (can go in supplementary)
point.aoc.bias <- enmtools.aoc(clade = anilios.clade.bias, nreps = 100, overlap.source = "points")
save(point.aoc.bias, file = 'data/script_generated_data/point_aoc_bias.Rdata')

# ENM overlap with GLM (non linear)
glm.aoc.d.bias <- enmtools.aoc(clade = anilios.clade.bias,  env = Biolayers2, nreps = 100, overlap.source = "glm", f = pres ~ poly(MeanTemp, 2) + poly(BulkDensity, 2) + poly(Aridity, 2), metric = "D")
save(glm.aoc.d.bias, file = 'data/script_generated_data/glm_aoc_d_bias.Rdata')

glm.aoc.i.bias <- enmtools.aoc(clade = anilios.clade.bias,  env = Biolayers2, nreps = 100, overlap.source = "glm", f = pres ~ poly(MeanTemp, 2) + poly(BulkDensity, 2) + poly(Aridity, 2), metric = "I")
save(glm.aoc.i.bias, file = 'data/script_generated_data/glm_aoc_i_bias.Rdata')

# ENM Overlap with MaxEnt and calculating 
mx.aoc.i.bias <- enmtools.aoc(clade = anilios.clade.bias,  env = Biolayers2, nreps = 100, overlap.source = "mx", metric = "I")
save(mx.aoc.i.bias, file = 'data/script_generated_data/mx_aoc_i_bias.Rdata')
mx.aoc.i.bias

mx.aoc.d.bias <- enmtools.aoc(clade = anilios.clade.bias,  env = Biolayers2, nreps = 100, overlap.source = "mx", metric = "D")
save(mx.aoc.d.bias, file = 'data/script_generated_data/mx_aoc_d_bias.Rdata')



# Building Environmental Niche Models for each species --------------------

### MAXENT fitter function
maxent.enm.maker <- function(.x){
  sp.mx <- enmtools.maxent(.x, env = Biolayers2, test.prop = 0.2,
                           bg.source = "point", env.nback = 10000)
  return(sp.mx)
}

## lapply MaxEnt fitter
sp.list.enms <- lapply(sp_list, FUN = maxent.enm.maker)
# test <- lapply(sp_list, function(x) enmtools.maxent(x, env = Biolayers2, test.prop = 0.2, bg.source = "point", env.nback = 10000, verbose = TRUE))

# Use lapply to save the rest of maxent object without the suitability raster
lapply(names(sp.list.enms), function(species_name) {
  # Extract the suitability data for the current species
  species.object.data <- sp.list.enms[[which(names(sp.list.enms) == species_name)]][c("analysis.df", "species.name", "training.evaluation", "test.evaluation", "env.training.evaluation", "env.test.evaluation")]
  # Create a file name based on the species name
  file_name <- paste0("data/script_generated_data/enmtools_sp_object_", species_name, ".rds")
  # Save the suitability as an RDS file
  saveRDS(species.object.data, file = file_name)
})

# Use lapply to save the suitability for each species
lapply(names(sp.list.enms), function(species_name) {
  # Extract the suitability data for the current species
  suitability_data <- sp.list.enms[[which(names(sp.list.enms) == species_name)]]$suitability
  # Create a file name based on the species name
  file_name <- paste0("data/script_generated_data/max_suitability_", species_name, ".rds")
  # Save the suitability as an RDS file
  terra::saveRDS(suitability_data, file = file_name)
})

pair_suit_plotter <- function(sp1, sp2=NULL, sp3=NULL){
  # sp1 sp2 sp3 has to be enmtools species object at least with occurrence data
  # will need to change this if sp.list.enms need to be loaded in a separate R session 
  plot(sp1$suitability); points(sp1$presence.points); text(x = 120, y = -40, label = sp1$species.name)
  if (!is.null(sp2)) {
    points(sp2$presence.points, col = "red"); 
  }
  
  if (!is.null(sp2)) {
    plot(sp2$suitability); points(sp2$presence.points); 
    points(sp1$presence.points, pch = 6, col = "red"); text(x = 120, y = -40, label = sp2$species.name)
  }

    if (!is.null(sp3)) {
    plot(sp3$suitability); points(sp3$presence.points); text(x = 120, y = -40, label = sp3$species.name)
    # points(sp1$presence.points, pch = 6, col = "red"); text(x = 120, y = -40, label = sp2$species.name)
  }
  
}

  
par(mfrow=c(2,2))

pair_suit_plotter(sp1 = sp.list.enms$affinis, sp2 = sp.list.enms$grypusET)

# save(sp.list.enms, file = 'data/script_generated_data/enmtools/sp_list_enms.Rdata')
# load('data/script_generated_data/enmtools/sp_list_enms.Rdata')
names(sp.list.enms)

### Can check fit using
# plot(sp.list.enms$silvia$test.evaluation, "ROC")
# sp.list.enms$silvia$training.evaluation@auc

### Caluclate AUC diff score using for loop and make data frame
names(sp.list.enms)
sp_df <- data.frame(species = snake_species)
sp_df$AUC_diff <- NULL
for (sp in 1:length(sp_df$species)){
  sp_df$AUC_test[[sp]] <- sp.list.enms[[sp]]$test.evaluation@auc
  sp_df$AUC_diff[[sp]] <- sp.list.enms[[sp]]$training.evaluation@auc - sp.list.enms[[sp]]$test.evaluation@auc
}

### Function to calculate Variable Importance Permutation
enmtools.vipper <- function(.model, nsim = 5){
  
  model <- .model
  thismodel <- model$model
  train <- rbind(attr(thismodel, "presence"), attr(thismodel, "absence"))
  feature_names <- colnames(train)
  train$presence <- as.factor(c(rep(1, nrow(attr(thismodel, "presence"))), 
                                rep(0, nrow(attr(thismodel, "absence")))))
  target <- "presence"
  pred_wrapper <- function(object, newdata) predict(object, newdata)
  
  vip.table <- vip::vi_permute(object = thismodel, 
                  feature_names = feature_names, 
                  train = train, 
                  target = target, 
                  metric = "roc_auc", 
                  pred_wrapper = pred_wrapper, 
                  event_level = "second",
                  nsim = nsim, keep = TRUE,
                  smaller_is_better = FALSE)

  # Calculate % of permutation 
  vip.table$Percent <- vip.table$Importance / sum(vip.table$Importance) * 100
  return(vip.table)
}

# Calculate Variable Importance using permutations
sp.list.enms.vip <- lapply(sp.list.enms, FUN = enmtools.vipper, nsim = 100)
sp.list.enms.vip

for (sp in 1:length(sp_df$species)){
  sp_df$MeanTemp[[sp]] <- unlist(sp.list.enms.vip[[sp]]$Percent[1], use.names = FALSE)
  sp_df$BulkDensity[[sp]] <-unlist(sp.list.enms.vip[[sp]]$Percent[2], use.names = FALSE)
  sp_df$Aridity[[sp]] <- unlist(sp.list.enms.vip[[sp]]$Percent[3], use.names = FALSE)
}

sp_df

# sp_df$AUC_diff <- unlist(sp_df$AUC_diff, use.names = FALSE)
# sp_df$AUC_test <- unlist(sp_df$AUC_test, use.names = FALSE)
# sp_df$MeanTemp <- unlist(sp_df$MeanTemp, use.names = FALSE)
# sp_df$MeanTemp <- unlist(sp_df$MeanTemp, use.names = FALSE)
# sp_df$BulkDensity <- unlist(sp_df$BulkDensity, use.names = FALSE)
# sp_df$Aridity <- unlist(sp_df$Aridity, use.names = FALSE)
# 
# sp_df$species <- paste("Anilios", sp_df$species)

### Write csv
# write.csv(sp_df, file = 'manuscript_word/supplements/species_enm_auc_vip.csv', row.names = FALSE)





