# 04_morpho_pairwise.R
# Sarin Tiatragul
# Oct 2023

# Do pariwise comparisons

# ENMTOOLS ----------------------------------------------------------------
library(devtools)
# install_github("danlwarren/ENMTools")
library(ENMTools); library(beepr)
library(dplyr); library("gridExtra")
library(ggpubr); library(ggplot2)
'%notin%' <- Negate('%in%')
library(ape); library(phangorn)

# Read in data locality data
load('data/script_generated_data/soil_bulk_density_data.rda')

# morpho data
load('data/script_generated_data/three_d_mshape.rda')

# Make a dataframe of occurrence for Anilios
snakes <- soil_df %>% 
  dplyr::filter(-43 < latdec & latdec < -6.3) %>% # 
  dplyr::filter(species %notin% c("sp.", "braminus", "", "polygrammicus", "schlegelii", "kunuaensis", "infralabialis", 
                                  "curtus", "cumingii", "angusticeps", "bibroni", "exocoeti", "fornasinii",
                                  "willeyi", "solomonis", "depressus")) %>% 
  dplyr::filter(!is.na(BIO1)) %>% dplyr::filter(!is.na(BD)) %>% dplyr::filter(!is.na(ARID))

###
### Input phylogeny
###

# Subset for only species that have >= 5 occurrences to use all available species
subset_snakes <- snakes
subset_snakes_sp <- unique(subset_snakes$species)

# Load phylogeny
load('data/script_generated_data/subset_mcmctree_shape.rda')

# Subset phylogeny
sub_phy$tip.label <- gsub(x = sub_phy$tip.label, pattern = "Anilios_", replacement = "")
keepers <- sub_phy$tip.label[sub_phy$tip.label %in% subset_snakes_sp]
tr_anilios <- ape::keep.tip(phy = sub_phy, tip  = keepers)

# create a list of species names
snake_species <- tr_anilios$tip.label

# create an empty list to store the output
sp_list <- list()

# loop through the species names to build emntools.species object
for (i in seq_along(snake_species)) {
  
  sp.obj <- enmtools.species()
  sp.obj$species.name <- snake_species[i]
  
  # Check the species object are fine
  # sp.obj <- check.species(sp.obj, env = Biolayers2)
  
  # append the sp.obj object to the sp_list
  sp_list[[i]] <- sp.obj
}

# Name items in the list of enmtools.species object
names(sp_list) <- snake_species


# Identify species pair and closest related species -----------------------

dd <- lapply(1:tr_anilios$Nnode + Ntip(tr_anilios), function(n, t)
    phangorn::Descendants(t, n)[[1]], t = tr_anilios)
nodes <- c(1:tr_anilios$Nnode + Ntip(tr_anilios))[which(sapply(dd, length) == 2)]
sisters <- t(sapply(nodes, function(n, t)
  t$tip.label[Descendants(t, n)[[1]]], t = tr_anilios))
rownames(sisters) <- nodes
sisters