# 03_morpho_overlap.R
# Sarin Tiatragul
# Oct 2023

# Combine results from AOC into figure

# ENMTOOLS ----------------------------------------------------------------
library(devtools)
# install_github("danlwarren/ENMTools")
library(ENMTools); library(beepr)
library(dplyr); library("gridExtra")
library(ggpubr); library(ggplot2)
'%notin%' <- Negate('%in%')

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

### 
### Build enmtools.clade object
###

### Build species.clade object
anilios.clade <- enmtools.clade(species = sp_list[c(tr_anilios$tip.label)], tree = tr_anilios)
anilios.clade <- check.clade(anilios.clade)

## HEAD SHAPE PAIRWISE MATRIX

# Filter only species we have
three_d_mshape # this object from 00_tps_landmark_to_geomorph_data.R

sub_shape <- three_d_mshape[c(tr_anilios$tip.label)]
names(sub_shape)
two_d_sp_coords <- do.call(rbind, lapply(sub_shape, c))
pairwise_headshape <- dist(two_d_sp_coords, diag = T, upper = T)  %>% as.matrix()
inv_head <- 1 / (pairwise_headshape + 1)

# AOC for head shape
head_aoc <- enmtools.aoc(clade = anilios.clade, nrep = 100, overlap.source = "matrix", overlap.matrix = inv_head)
head_aoc$regressions.plot 
save(head_aoc, file = 'data/script_generated_data/head_aoc.Rdata')

### for body shape
body_traits <- read.csv(file = 'data/script_generated_data/shape_ratios/blindsnake_sp_pc_logbodyrat.csv')

# Make pairwise matrix with midbody width
rownames(body_traits) <- sub("^[^ ]+ ", "", body_traits$species)
sub_body <- body_traits[tr_anilios$tip.label,]
pairwise_body <- dist(sub_body["sh_mbd"]  , diag = T, upper = T) %>% as.matrix()

inv <- 1 / (pairwise_body + 1)

mbd_aoc <- enmtools.aoc(clade = anilios.clade, nrep = 100, overlap.source = "matrix", overlap.matrix = inv)
save(mbd_aoc, file = 'data/script_generated_data/mbd_aoc.Rdata')

## Load the niche and point overlap
load(file = 'data/script_generated_data/range_aoc_bias.Rdata')
load(file = 'data/script_generated_data/mx_aoc_i_bias.Rdata')

## ARRANGE GGPLOT

r <- range.aoc.bias$regressions.plot + theme_classic() + xlab("Age") + ylab("Geographical overlap (%)")
mx <- mx.aoc.i.bias$regressions.plot + theme_classic() + xlab("Age") + ylab("Niche overlap (%)")
m <- mbd_aoc$regressions.plot + theme_classic() + xlab("Age") + ylab("Body shape overlap (%)")
h <- head_aoc$regressions.plot + theme_classic() + xlab("Age") + ylab("Head shape overlap (%)") 

one_page <- ggarrange(mx, r,
                      m, h, ncol = 2, nrow = 2)
ggexport(one_page, filename = "output/AOC_plots.pdf")


summary(mx.aoc.i.bias)
mx.aoc.i.bias$coefficients

summary(range.aoc.bias)
range.aoc.bias$coefficients

mbd_aoc$coefficients
head_aoc$coefficients



