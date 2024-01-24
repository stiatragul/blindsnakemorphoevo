# 04_pairwise_distance_test.R
# Sarin Tiatragul
# Oct 2023

# PAIRWISE DISTANCE ANALYSES BETWEEN SISTER SPECIES PAIRS AND NON-SISTER
# For Geographical Overlap and Niche Overlap


# Packages ----------------------------------------------------------------
# library(devtools)
# install_github("danlwarren/ENMTools") ; install_github('luismurao/ntbox')
library(terra); library(raster);
library(dplyr); library(tidyr)
library(ggplot2); library(ape); 
library(phytools); library(phangorn)
'%notin%' <- Negate('%in%')

# Data --------------------------------------------------------------------

### Data for this analysis has been generated elsewhere

# Read in data locality data
load('data/script_generated_data/soil_bulk_density_data.rda')

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

### Input phylogeny

# Subset for only species that have >= 5 occurrences to use all available species
species_count <- table(snakes$species)
subset_snakes <- snakes[snakes$species %in% names(species_count[species_count >= 5]),]
subset_snakes_sp <- unique(subset_snakes$species)

# Load phylogeny
# Non-time calibrated for genetic distance
astral_tree <- read.tree(file = "data/script_generated_data/astral_tree.tre")
load('data/script_generated_data/subset_mcmctree_shape.rda')

# Subset phylogeny
sub_phy$tip.label <- gsub(x = sub_phy$tip.label, pattern = "Anilios_", replacement = "")
astral_tree$tip.label <- gsub(x = astral_tree$tip.label, pattern = "Anilios_", replacement = "")
keepers <- sub_phy$tip.label[sub_phy$tip.label %in% subset_snakes_sp]
tr_anilios <- ape::keep.tip(phy = sub_phy, tip  = keepers)
keepers.morph <- sub_phy$tip.label[sub_phy$tip.label %in% unique(snakes$species)]
tr_anilios.morph <- ape::keep.tip(phy = sub_phy, tip  = keepers.morph)
ast_tree <- ape::keep.tip(phy = astral_tree, tip  = keepers)

plot(tr_anilios)
plot(tr_anilios.morph)
plot(ast_tree)

### Make sure we use the correct species pairs for niche/range and morpho because we have different amount of data 
# Get vector of species that have at least 5 points for niche
n_sample %>% dplyr::filter(n >= 5) %>% unique() %>% arrange(n)
sp_niche_calc <- n_sample %>% dplyr::filter(n >= 5) %>% unique() %>% arrange(species)

# List of species object --------------------------------------------------
# create a list of species names
snake_species <- sp_niche_calc$species

# Create all pairwise combination
pair_sp_pre <- combn(snake_species, 2) %>% t() %>% as.data.frame()

### Load in MaxEnt species object
load('data/script_generated_data/enmtools/sp_list_enms.Rdata')

pair_sp <- pair_sp_pre

# Identify species pair and closest related species -----------------------

dd <- lapply(1:tr_anilios$Nnode + Ntip(tr_anilios), function(n, t)
  phangorn::Descendants(t, n)[[1]], t = tr_anilios)
nodes <- c(1:tr_anilios$Nnode + Ntip(tr_anilios))[which(sapply(dd, length) == 2)]
sisters <- t(sapply(nodes, function(n, t)
  t$tip.label[Descendants(t, n)[[1]]], t = tr_anilios))
rownames(sisters) <- nodes
sister_pairs <- as.data.frame(sisters)

## For morphology data
dd <- lapply(1:tr_anilios.morph$Nnode + Ntip(tr_anilios.morph), function(n, t)
  phangorn::Descendants(t, n)[[1]], t = tr_anilios.morph)
nodes <- c(1:tr_anilios.morph$Nnode + Ntip(tr_anilios.morph))[which(sapply(dd, length) == 2)]
sisters.m <- t(sapply(nodes, function(n, t)
  t$tip.label[Descendants(t, n)[[1]]], t = tr_anilios.morph))
rownames(sisters.m) <- nodes
sister_pairs.morph <- as.data.frame(sisters.m)

# Make data frame where non-sister species will find their second closest relatives
# non_sister <- data.frame(sp = c(sister_pairs$V1, sister_pairs$V2), pair_id = c(rownames(sister_pairs), rownames(sister_pairs)))
#### Get figure out which genetic distance is most similar
# get_pairwise_value <- function(row) {
#   value <- genetic_matrix[row["V1"], row["V2"]]
#   return(value)
# }
# 
# for (i in 1:length(non_sister$sp)) {
#   
# non_sister$sp2[i] <- names(which(genetic_matrix[, non_sister$sp[i]] == sort(unique(genetic_matrix[, non_sister$sp[i]]))[2]))
#   
# }

### MANUAL Method of finding species pair based on topology
# I selected immediate outgroups as the congener. If the outgroup is a species pair then the second species is also included 
# as non-sister2. This allow us to do some sort of sensitivity check.
# We will take the average distance between both sister to the congener

### Manually make species pair file based on topology of the tree
manual_pair <- read.csv(file = 'data/manual_sister_nonsister_pairwise_matrix.csv')
manual_pair.morph <- read.csv(file = 'data/manual_sister_nonsister_pairwise_matrix_morph.csv')

### Genetic distance
genetic_matrix <- cophenetic.phylo(ast_tree)
manual_pair$genetic_dist <- genetic_matrix[cbind(manual_pair$V1, manual_pair$V2)]

### Morphological overlap
load(file = 'data/script_generated_data/mbd_aoc.Rdata')
load(file = 'data/script_generated_data/head_aoc.Rdata')

manual_pair.morph$head_dist <- head_aoc$empirical.overlap[cbind(manual_pair.morph$V1, manual_pair.morph$V2)]
manual_pair.morph$mbd_dist <- mbd_aoc$empirical.overlap[cbind(manual_pair.morph$V1, manual_pair.morph$V2)]

### Geographical overlap
load('data/script_generated_data/range_aoc_bias.Rdata')
manual_pair$range_overlap <- range.aoc.bias$empirical.overlap[cbind(manual_pair$V1, manual_pair$V2)]

### Ecological Niche overlap
load('data/script_generated_data/mx_aoc_i_bias.Rdata')
load('data/script_generated_data/mx_aoc_d_bias.Rdata')

manual_pair$niche_i_overlap <- unlist(mx.aoc.i.bias$empirical.overlap[cbind(manual_pair$V1, manual_pair$V2)])
manual_pair$niche_d_overlap <- unlist(mx.aoc.d.bias$empirical.overlap[cbind(manual_pair$V1, manual_pair$V2)])

# Calculate averaged non-sister values
averaged_congener <- manual_pair %>% 
  dplyr::filter(other %in% c("non-sister1", "non-sister2")) %>% 
  dplyr::group_by(pairID, other) %>% 
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  dplyr::arrange(other)

averaged_congener.morph <- manual_pair.morph %>% 
  dplyr::filter(other %in% c("non-sister1", "non-sister2")) %>% 
  dplyr::group_by(pairID, other) %>% 
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  dplyr::arrange(other)

cong.1 <- averaged_congener %>% dplyr::filter(other == "non-sister1") 
cong.2 <- averaged_congener %>% dplyr::filter(other == "non-sister2") 
cong.morph.1 <- averaged_congener.morph %>% dplyr::filter(other == "non-sister1") 
cong.morph.2 <- averaged_congener.morph %>% dplyr::filter(other == "non-sister2") 

# Trim sister pair as well
sister_overlaps <- manual_pair %>% dplyr::filter(other == "sister") %>% dplyr::select(-V1, -V2)
sister_overlaps.morph <- manual_pair.morph %>% dplyr::filter(other == "sister") %>% dplyr::select(-V1, -V2)

### Write raw metric comparison
sp_pair_overlap_1 <- rbind(sister_overlaps, cong.1)
sp_pair_overlap_1$sp1 <- c(manual_pair$V1[1:11], paste(manual_pair$V1[1:11], "+", paste(manual_pair$V2[1:11])))
sp_pair_overlap_1$sp2 <- c(manual_pair$V2[1:11], manual_pair$V2[12:22])
sp_pair_overlap_2 <- rbind(sister_overlaps, cong.2)
sp_pair_overlap_2$sp1 <- c(manual_pair$V1[1:11], paste(manual_pair$V1[1:11], "+", paste(manual_pair$V2[1:11])))
sp_pair_overlap_2$sp2 <- c(manual_pair$V2[1:11], manual_pair$V2[34:44])

sp_pair_overlap_morph_1 <- rbind(sister_overlaps.morph, cong.morph.1)
sp_pair_overlap_morph_1$sp1 <- c(manual_pair.morph$V1[1:12], paste(manual_pair.morph$V1[1:12], "+", paste(manual_pair.morph$V2[1:12])))
sp_pair_overlap_morph_1$sp2 <- c(manual_pair.morph$V2[1:12], manual_pair.morph$V2[13:24])

sp_pair_overlap_morph_2 <- rbind(sister_overlaps.morph, cong.morph.2)
sp_pair_overlap_morph_2$sp1 <- c(manual_pair.morph$V1[1:12], paste(manual_pair.morph$V1[1:12], "+", paste(manual_pair.morph$V2[1:12])))
sp_pair_overlap_morph_2$sp2 <- c(manual_pair.morph$V2[1:12], manual_pair.morph$V2[37:48])

write.csv(sp_pair_overlap_1, 'manuscript_word/supplements/sp_overlap_1_env.csv',row.names = F)
write.csv(sp_pair_overlap_2, 'manuscript_word/supplements/sp_overlap_2_env.csv',row.names = F)
write.csv(sp_pair_overlap_morph_1, 'manuscript_word/supplements/sp_overlap_1_morph.csv',row.names = F)
write.csv(sp_pair_overlap_morph_2, 'manuscript_word/supplements/sp_overlap_2_morph.csv',row.names = F)

# Exact binomial test -----------------------------------------------------

### MORPHOLOGY BINOMIAL TEST
diff.ov.morph.1 <- data.frame(pairID = sister_overlaps.morph$pairID)
diff.ov.morph.2 <- data.frame(pairID = sister_overlaps.morph$pairID)

# Calculate metric difference between sister and non-sister
diff.ov.morph.1$head <- sister_overlaps.morph$head_dist - cong.morph.1$head_dist
diff.ov.morph.1$mbd <- sister_overlaps.morph$mbd_dist - cong.morph.1$mbd_dist

diff.ov.morph.2$head <- sister_overlaps.morph$head_dist - cong.morph.2$head_dist
diff.ov.morph.2$mbd <- sister_overlaps.morph$mbd_dist - cong.morph.2$mbd_dist

# Perform exact binomial test for morph differences
binom.head.test.1 <- binom.test(sum(diff.ov.morph.1$head > 0), length(diff.ov.morph.1$head))
binom.mbd.test.1 <- binom.test(sum(diff.ov.morph.1$mbd > 0), length(diff.ov.morph.1$mbd))
binom.head.test.2 <- binom.test(sum(diff.ov.morph.2$head > 0), length(diff.ov.morph.2$head))
binom.mbd.test.2 <- binom.test(sum(diff.ov.morph.2$mbd > 0), length(diff.ov.morph.2$mbd))

# Check if the results are about the same to see if test is sensitive to which is considered congener species
binom.head.test.1; binom.head.test.2
binom.mbd.test.1; binom.mbd.test.2

### Interpretation
# Sister species resemble head shape more than non-sister species p = 0.006 (11/12)
# Sister species do not significantly differ to non-sister species in log relative mbd / svl p = 0.38 (4/12)

### NICHE AND RANGE BINOMIAL TEST

diff.ov.1 <- data.frame(pairID = sister_overlaps$pairID)
diff.ov.2 <- data.frame(pairID = sister_overlaps$pairID)

diff.ov.1$range <- sister_overlaps$range_overlap - cong.1$range_overlap
diff.ov.1$niche_i <- sister_overlaps$niche_i_overlap - cong.1$niche_i_overlap
diff.ov.1$niche_d <- sister_overlaps$niche_d_overlap - cong.1$niche_d_overlap

diff.ov.2$range <- sister_overlaps$range_overlap - cong.2$range_overlap
diff.ov.2$niche_i <- sister_overlaps$niche_i_overlap - cong.2$niche_i_overlap
diff.ov.2$niche_d <- sister_overlaps$niche_d_overlap - cong.2$niche_d_overlap

binom.range.test.1 <- binom.test(sum(diff.ov.1$range > 0), length(diff.ov.1$range))
binom.niche_i.test.1 <- binom.test(sum(diff.ov.1$niche_i > 0), length(diff.ov.1$range))
binom.niche_d.test.1 <- binom.test(sum(diff.ov.1$niche_d > 0), length(diff.ov.1$range))

binom.range.test.2 <- binom.test(sum(diff.ov.2$range > 0), length(diff.ov.2$range))
binom.niche_i.test.2 <- binom.test(sum(diff.ov.2$niche_i > 0), length(diff.ov.2$range))
binom.niche_d.test.2 <- binom.test(sum(diff.ov.2$niche_d > 0), length(diff.ov.2$range))

# Check for sensitivity with congener selection...not different
binom.range.test.1; binom.range.test.2
binom.niche_i.test.1; binom.niche_i.test.2
binom.niche_d.test.1; binom.niche_d.test.2

### Interpretation 
# Sister species do not have greater range overlap than non-sister species; p = 1, 6/11
# Sister species do not have significantly greater niche overlap (I) than non-sister species; p = 0.23; 8/11
# Sister species do not have significantly greater niche overlap (D) than non-sister species; p = 0.23/ 8/11

### Make results table
sister_mean <- sister_overlaps %>% dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, mean, na.rm=TRUE) 
sister_sd <- sister_overlaps %>% dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, sd, na.rm=TRUE)
sister_mean.morph <- sister_overlaps.morph %>% dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, mean, na.rm=TRUE)
sister_sd.morph <- sister_overlaps.morph %>% dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, sd, na.rm=TRUE)

congener_mean <- cong.1 %>%dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, mean, na.rm=TRUE)
congener_sd <- cong.1 %>% dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, sd, na.rm=TRUE)
congener_mean.morph <- cong.morph.1 %>%dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, mean, na.rm=TRUE)
congener_sd.morph <- cong.morph.1 %>% dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, sd, na.rm=TRUE)

congener_mean2 <- cong.1 %>%dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, mean, na.rm=TRUE)
congener_sd2 <- cong.1 %>% dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, sd, na.rm=TRUE)
congener_mean.morph2 <- cong.morph.2 %>%dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, mean, na.rm=TRUE)
congener_sd.morph2 <- cong.morph.2 %>% dplyr::group_by(other) %>% dplyr::summarise_if(is.numeric, sd, na.rm=TRUE)


### Result data frame to include in supplement

pairwise_dist_results <- data.frame(metric = c("range", "niche_i", "niche_d", "head", "body"),
                                    number_pairs = c(rep(length(diff.ov.1$range), 3), 12, 12),
                                    sister_means = c(unlist(sister_mean[1, 4:6]), unlist(sister_mean.morph[1, 3:4])),
                                    sister_sds = c(unlist(sister_sd[1, 4:6]), unlist(sister_sd.morph[1, 3:4])),
                                    congener_means = c(unlist(congener_mean[1, 4:6]), unlist(congener_mean.morph[1, 3:4])),
                                    congener_sds = c(unlist(congener_sd[1, 4:6]), unlist(congener_sd.morph[1, 3:4])),
                                    freq = c(binom.range.test.1$statistic/11,
                                             binom.niche_i.test.1$statistic/11,
                                             binom.niche_d.test.1$statistic/11,
                                             binom.head.test.1$statistic/12,
                                             binom.mbd.test.1$statistic/12),
                                    p_value = c(binom.range.test.1$p.value,
                                                binom.niche_i.test.1$p.value,
                                                binom.niche_d.test.1$p.value,
                                                binom.head.test.1$p.value,
                                                binom.mbd.test.1$p.value))

write.csv(pairwise_dist_results, file = 'manuscript_word/supplements/sister_v_congener_binom_test.csv', row.names = FALSE)


congener_mean

pairwise_dist_results2 <- data.frame(metric = c("range", "niche_i", "niche_d", "head", "body"),
                                    number_pairs = c(rep(length(diff.ov.1$range), 3), 12, 12),
                                    sister_means = c(unlist(sister_mean[1, 4:6]), unlist(sister_mean.morph[1, 3:4])),
                                    sister_sds = c(unlist(sister_sd[1, 4:6]), unlist(sister_sd.morph[1, 3:4])),
                                    congener_means = c(unlist(congener_mean2[1, 4:6]), unlist(congener_mean.morph2[1, 3:4])),
                                    congener_sds = c(unlist(congener_sd2[1, 4:6]), unlist(congener_sd.morph2[1, 3:4])),
                                    freq = c(binom.range.test.2$statistic/11,
                                             binom.niche_i.test.2$statistic/11,
                                             binom.niche_d.test.2$statistic/11,
                                             binom.head.test.2$statistic/12,
                                             binom.mbd.test.2$statistic/12),
                                    p_value = c(binom.range.test.2$p.value,
                                                binom.niche_i.test.2$p.value,
                                                binom.niche_d.test.2$p.value,
                                                binom.head.test.2$p.value,
                                                binom.mbd.test.2$p.value))

write.csv(pairwise_dist_results2, file = 'manuscript_word/supplements/sister_v_congener_binom_test_2.csv', row.names = FALSE)
