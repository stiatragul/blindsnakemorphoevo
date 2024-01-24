# 02_convergence.R
# Sarin Tiatragul 
# March 2023; modified October 2023

# libraries ---------------------------------------------------------------
# Load packages with two lines
pkgs = c("phytools", "ape", "geomorph", "ggplot2", "dplyr", "stringr", "corHMM", "geiger", "nlme", "convevol", "geiger", "scico") # package names
inst = lapply(pkgs, library, character.only = TRUE) # load them
'%notin%' <- Negate('%in%')
source("code/utility/func_pwCheck_2.1.R") # pre-release code from William Brightly for convevol 2.1

# data --------------------------------------------------------------------
load('data/script_generated_data/subset_mcmctree_shape.rda')
# Load two D array of y symmetric component of dorsal head shape
load('data/script_generated_data/dorsal_head_shape.rda')

blindsnake_data <- read.csv(file = 'data/script_generated_data/mean_blindsnake_data.csv', row.names = 1)

# aridty type data
arid_type_df <- read.csv(file = 'data/blindsnake_arid_strict.csv', row.names = 1)

# Drop outgroups from analyses 
anilios_tree <- ape::drop.tip(phy = sub_phy,
                              tip = c("Ramphotyphlops_multilineatus",
                                      "Acutotyphlops_subocularis"))
anilios_data <- blindsnake_data[rownames(blindsnake_data) %in% anilios_tree$tip.label,]

# Change level names
anilios_data$mean_bulk_cat <- dplyr::recode_factor(.x = anilios_data$mean_bulk_cat, "1" = "soft", "2" = "med", "3" = "hard")

# Match order of tips and data
anilios_data <- anilios_data[match(anilios_tree$tip.label, anilios_data$species), ]

# Add arid type
anilios_data$arid_type <- arid_type_df[anilios_tree$tip.label,]$arid_type

# CONVERGENCE C(t)-MEASURES -----------------------------------------------
## Stayton 2015. and follow up Grossnickle et al. 2023

# C1 = magnitude of of morph convergence in focal extant taxa relative to the max divergence in their ancestral values.
# Large values represent greater degree of phenotypic convergence.
# C1 = 1-(current distance / maximum ancestral distance); 
# C2 = maximum ancestral distance - current distance; 
# C3=C2/(total phenotypic evolution in the clade defined by the two taxa); 
# C4 = C2/(total amount of phenotypic evolution in the entire tree). 
# If more than two convergent taxa are input, then C1-C4 are calculated for all 
# possible pairs of taxa, and averaged.
# C5 = freq of convergence to particular region of morphospace, estimated by number of focal lineages whose evolution transects the boundary region defined by the focal taxa

# Prepare bodyshape traits data (same with 04_model_fitting.R for mvMORPH)
BodyShape <- as.matrix(anilios_data[, c("sh_mbd", "sh_mtw", "sh_hwe", "sh_hda", "sh_tl")])
BodyShapePC <- as.matrix(anilios_data[, c("fil_PC1", "fil_PC2")])

name.check(phy = anilios_tree, BodyShape); name.check(phy = anilios_tree, BodyShapePC)

# Prepare head shape trait data

rownames(sp_means_dorsal) <- gsub(x = rownames(sp_means_dorsal), pattern = "^", replacement = "Anilios_")
means_2d_dorsal <- sp_means_dorsal[rownames(sp_means_dorsal) %in% anilios_tree$tip.label,]
means_2d_dorsal <- means_2d_dorsal[match(anilios_tree$tip.label, rownames(means_2d_dorsal)), ]

HeadLMsM <- as.matrix(anilios_data[, c("mean_PC1", "mean_PC2")])
head(HeadLMsM)

### Test for convergence in aridity category
arid_sp <- rownames(anilios_data)[which(anilios_data$arid_category == "Arid")]

# stricter category of arid based on me looking at satellite images of where most the locality came from
arid_strict <- read.csv('data/blindsnake_arid_strict.csv', row.names = 2)
arid_sp_strict <- rownames(arid_strict)[which(arid_strict$aridity_strict == "arid")]

## We are also comparing between species that are in different desert types (with dune) according to Mabbutt 1988
arid_uplan <- rownames(arid_strict)[which(arid_strict$upland == "upland")]
arid_dunes <- rownames(arid_strict)[which(arid_strict$sand == "sand")]
arid_stony <- rownames(arid_strict)[which(arid_strict$stony == "stony")]
arid_shied <- rownames(arid_strict)[which(arid_strict$shield == "shield")]

# semi_sp <- rownames(anilios_data)[which(anilios_data$arid_category == "Semi-Arid")]
humd_sp <- rownames(anilios_data)[which(anilios_data$arid_category == "Humid")]

### Check which group is uninformative for comparison and group them. 
# Using script for pwCheck() provided by Will Brightly on 2023-10-24 before release of new convevol version

arid_check <- pwCheck(phy = anilios_tree, focaltaxa = arid_sp)
arid_check

arid_strict_check <- pwCheck(phy = anilios_tree, focaltaxa = arid_sp_strict)
trim_arid_strict <- arid_sp_strict[arid_sp_strict %notin% c("Anilios_waitii", "Anilios_endoterus", "Anilios_diversus")]

# make new groups for arid species
arid_group <- rep("G1", length(arid_sp))
names(arid_group) <- arid_sp

# check humid group
pwCheck(phy = anilios_tree, focaltaxa = humd_sp)

# assign grouping variables according to the output from pwCheck()
arid_group[which(names(arid_group) %in%  c("Anilios_longissimus", "Anilios_grypusW", "Anilios_obtusifrons", "Anilios_leptosoma", "Anilios_systenos"))] <- "G1"
arid_group[which(names(arid_group) ==  "Anilios_bicolor")] <- "G2"
arid_group[which(names(arid_group) %in%  c("Anilios_waitii", "Anilios_centralis"))] <- "G3"
arid_group[which(names(arid_group) %in%  c("Anilios_pilbarensis", "Anilios_hamatus", "Anilios_endoterus"))] <- "G4"
arid_group[which(names(arid_group) == "Anilios_ganei")] <- "G5"
arid_group[which(names(arid_group) %in% c("Anilios_margaretae", "Anilios_bituberculatus"))] <- "G6"
arid_group[which(names(arid_group) == "Anilios_aspina")] <- "G7"
arid_group[which(names(arid_group) == "Anilios_grypusET")] <- "G8"
arid_group[which(names(arid_group) == "Anilios_ammodytes")] <- "G9"
arid_group

pwCheck(phy = anilios_tree, focaltaxa = arid_sp, groups = arid_group)

## Arid type based on map by Mabbutt
# check upland
pwCheck(phy = anilios_tree, focaltaxa = arid_uplan)
uplan_group <- rep("G1", length(arid_uplan))
names(uplan_group) <- arid_uplan
uplan_group[c(1,2,5:8)] <- c("gryw", "lep", "pil", "end", "gan", "amm")
uplan_group[which(names(uplan_group) %in% c("Anilios_waitii", "Anilios_centralis"))] <- "waicen"

# check dune
pwCheck(phy = anilios_tree, focaltaxa = arid_dunes)
dune_group <- rep("G1", length(arid_dunes))
names(dune_group) <- arid_dunes
dune_group[c(1:5, 8)] <- c("grpw", "bic", "wai", "pil", "end", "div")
dune_group[which(names(dune_group) %in%  c("Anilios_margaretae", "Anilios_bituberculatus"))] <- "marbit"
pwCheck(phy = anilios_tree, focaltaxa = arid_dunes, groups = dune_group)

# check stony
pwCheck(phy = anilios_tree, focaltaxa = arid_stony)

# check shield
pwCheck(phy = anilios_tree, focaltaxa = arid_shied)

# check humid taxa too, doesn't seem to be any uninformative 
pwCheck(phy = anilios_tree, focaltaxa = humd_sp)
humd_group <- gsub(pattern = "Anilios_", replacement = "", humd_sp)
names(humd_group) <- humd_sp

dune_sp <- arid_dunes
shie_sp <- arid_shied
ston_sp <- arid_stony
land_sp <- arid_uplan

# Calc Convt --------------------------------------------------------------

# Test for convergence by aridity category
body.arid.arid <- convevol::calcConvCt(phy = anilios_tree, traits =BodyShape, focaltaxa = arid_sp, groups = arid_group)
body.arid.dune <- convevol::calcConvCt(phy = anilios_tree, traits =BodyShape, focaltaxa = dune_sp, groups = dune_group)
body.arid.land <- convevol::calcConvCt(phy = anilios_tree, traits =BodyShape, focaltaxa = land_sp, groups = uplan_group)
body.arid.ston <- convevol::calcConvCt(phy = anilios_tree, traits =BodyShape, focaltaxa = ston_sp)
body.arid.shie <- convevol::calcConvCt(phy = anilios_tree, traits =BodyShape, focaltaxa = shie_sp)
body.arid.humd <- convevol::calcConvCt(phy = anilios_tree, traits =BodyShape, focaltaxa = humd_sp)

body.arid.arid.sig <- convevol::convSigCt(phy = anilios_tree, traits =BodyShape, focaltaxa = arid_sp, groups = arid_group, nsim = 500)
body.arid.dune.sig <- convevol::convSigCt(phy = anilios_tree, traits =BodyShape, focaltaxa = dune_sp, groups = dune_group, nsim = 500)
body.arid.land.sig <- convevol::convSigCt(phy = anilios_tree, traits =BodyShape, focaltaxa = land_sp, groups = uplan_group, nsim = 500)
body.arid.ston.sig <- convevol::convSigCt(phy = anilios_tree, traits =BodyShape, focaltaxa = ston_sp, nsim = 500)
body.arid.shie.sig <- convevol::convSigCt(phy = anilios_tree, traits =BodyShape, focaltaxa = shie_sp, nsim = 500)
body.arid.humd.sig <- convevol::convSigCt(phy = anilios_tree, traits =BodyShape, focaltaxa = humd_sp, groups = humd_group, nsim = 500)

## arid category for head shape
head.arid.arid <- convevol::calcConvCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = arid_sp, groups = arid_group)
head.arid.dune <- convevol::calcConvCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = dune_sp, groups = dune_group)
head.arid.land <- convevol::calcConvCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = land_sp, groups = uplan_group)
head.arid.ston <- convevol::calcConvCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = ston_sp)
head.arid.shie <- convevol::calcConvCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = shie_sp)
head.arid.humd <- convevol::calcConvCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = humd_sp)

head.arid.arid.sig <- convevol::convSigCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = arid_sp, groups = arid_group, nsim = 500)
head.arid.dune.sig <- convevol::convSigCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = dune_sp, groups = dune_group, nsim = 500)
head.arid.land.sig <- convevol::convSigCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = land_sp, groups = uplan_group, nsim = 500)
head.arid.ston.sig <- convevol::convSigCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = ston_sp, nsim = 500)
head.arid.shie.sig <- convevol::convSigCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = shie_sp, nsim = 500)
head.arid.humd.sig <- convevol::convSigCt(phy = anilios_tree, traits =HeadLMsM, focaltaxa = humd_sp, groups = humd_group, nsim = 1000)

# save(head.arid.arid.sig, head.arid.humd.sig, body.arid.arid.sig, body.arid.humd.sig, head.arid.dune.sig, head.arid.land.sig, body.arid.dune.sig, body.arid.land.sig, file = 'data/script_generated_data/convergence_ct_score.Rdata')
# load(file = 'data/script_generated_data/convergence_ct_score.Rdata')

#### p-values for Ct1 - Ct4

### ARID

## ARID BODY
body.arid.arid.sig$pvals
body.arid.arid.sig$grp.mean; body.arid.arid.sig$grp.pvals
c_measures_result_body_arid <- rbind(body.arid.arid.sig$grp.mean, body.arid.arid.sig$grp.pvals)
write.csv(c_measures_result_body_arid, file = 'output/c-measures-groupmeans-arid.csv')
plotCt(body.arid.arid.sig, phy = anilios_tree, focaltaxa = arid_sp, groups = arid_group, nsim = 500, col = scico(n = 30, palette = "roma"))

## DUNE BODY
body.arid.dune.sig$pvals
body.arid.dune.sig$grp.mean; body.arid.dune.sig$grp.pvals
c_measures_result_body_dune <- rbind(body.arid.dune.sig$grp.mean, body.arid.dune.sig$grp.pvals)
# write.csv(c_measures_result_body_dune, file = 'output/c-measures-groupmeans-dune.csv')
plotCt(body.arid.dune.sig, phy = anilios_tree, focaltaxa = dune_sp, groups = dune_group, nsim = 500, col = scico(n = 30, palette = "roma"))

## UPLAND BODY
body.arid.land.sig$pvals
body.arid.land.sig$grp.mean; body.arid.land.sig$grp.pvals
c_measures_result_body_land <- rbind(body.arid.land.sig$grp.mean, body.arid.land.sig$grp.pvals)
# write.csv(c_measures_result_body_land, file = 'output/c-measures-groupmeans-land.csv')
plotCt(body.arid.land.sig, phy = anilios_tree, focaltaxa = land_sp, groups = uplan_group, nsim = 500, col = scico(n = 6, palette = "roma"))

# STONY BODY
body.arid.ston.sig$pvals
c_measures_result_body_ston <- rbind(body.arid.ston.sig$grp.mean, body.arid.ston.sig$grp.pvals)
# write.csv(c_measures_result_body_land, file = 'output/c-measures-groupmeans-land.csv')
plotCt(body.arid.ston.sig, phy = anilios_tree, focaltaxa = ston_sp, nsim = 500, col = scico(n = 6, palette = "roma"))

# SHIELD BODY
body.arid.shie.sig$pvals
c_measures_result_body_shie <- rbind(body.arid.ston.sig$grp.mean, body.arid.ston.sig$grp.pvals)
# write.csv(c_measures_result_body_shie, file = 'output/c-measures-groupmeans-land.csv')
plotCt(body.arid.shie.sig, phy = anilios_tree, focaltaxa = land_sp, groups = land_group, nsim = 500, col = scico(n = 6, palette = "roma"))


### HEAD

### ARID HEAD
head.arid.arid.sig$pvals
head.arid.arid.sig$grp.mean; head.arid.arid.sig$pvals
c_measures_result_head_arid <- rbind(head.arid.arid.sig$grp.mean, head.arid.arid.sig$grp.pvals)
# write.csv(c_measures_result_head_arid, file = 'output/c-measures-groupmeans-arid-head.csv')
plotCt(head.arid.arid.sig, phy = anilios_tree, focaltaxa = arid_sp, groups = arid_group, nsim = 500, col = scico(n = 30, palette = "roma"))

### DUNE HEAD
head.arid.dune.sig$pvals
head.arid.dune.sig$grp.mean; head.arid.dune.sig$pvals
c_measures_result_head_dune <- rbind(head.arid.dune.sig$grp.mean, head.arid.dune.sig$grp.pvals)
# write.csv(c_measures_result_head_dune, file = 'output/c-measures-groupmeans-dune-head.csv')
plotCt(head.arid.dune.sig, phy = anilios_tree, focaltaxa = dune_sp, groups = dune_group, nsim = 500, col = scico(n = 30, palette = "roma"))

### UPLAND HEAD
head.arid.land.sig$pvals
head.arid.land.sig$grp.mean
head.arid.land.sig$grp.mean; head.arid.land.sig$pvals
c_measures_result_head_land <- rbind(head.arid.land.sig$grp.mean, head.arid.land.sig$grp.pvals)
# write.csv(c_measures_result_head_land, file = 'output/c-measures-groupmeans-land-head.csv')
plotCt(head.arid.land.sig, phy = anilios_tree, focaltaxa = land_sp, groups = land_group, nsim = 500, col = scico(n = 30, palette = "roma"))

# STONY head
head.arid.ston.sig$pvals
# c_measures_result_body_ston <- rbind(body.arid.ston.sig$grp.mean, body.arid.ston.sig$grp.pvals)
# write.csv(head.arid.ston.sig$pvals, file = 'output/c-measures-groupmeans-stony-head.csv')
# plotCt(head.arid.ston.sig, phy = anilios_tree, focaltaxa = ston_sp, nsim = 500, col = scico(n = 6, palette = "roma"))

# SHIELD head
head.arid.shie.sig$pvals
# c_measures_result_body_shie <- rbind(body.arid.ston.sig$grp.mean, body.arid.ston.sig$grp.pvals)
# write.csv(head.arid.shie.sig$pvals, file = 'output/c-measures-groupmeans-shield-head.csv')
plotCt(head.arid.shie.sig, phy = anilios_tree, focaltaxa = shie_sp, nsim = 500, col = scico(n = 6, palette = "roma"))

### ARID GROUP INTERPRETATION
# Ct1 scores for comparison within arid species were negative, indicating no convergence.
# For body shape ratio, this negative value is statistically significant (P=0.01 overall).
# while the remaining  7/28 have positive values indicating convergence but these patterns were not significantly different than the null expectation (Table S c-measures-grbbody).
# None of the Ct1 values for head shape within arid species were statistically significant.

### HUMID SPECIES
body.arid.humd.sig$pvals
body.arid.humd.sig$grp.mean; body.arid.humd.sig$grp.pvals
c_measures_result_body_humd <- rbind(body.arid.humd.sig$grp.mean, body.arid.humd.sig$grp.pvals)
# write.csv(c_measures_result_body_humd, file = 'output/c-measures-groupmeans-humd-body.csv')

head.arid.humd.sig$pvals
head.arid.humd.sig$grp.mean; head.arid.humd.sig$grp.pvals
c_measures_result_head_humd <- rbind(head.arid.humd.sig$grp.mean, head.arid.humd.sig$grp.pvals)
# write.csv(c_measures_result_head_humd, file = 'output/c-measures-groupmeans-humd-head.csv')

plotCt(body.arid.humd.sig, phy = anilios_tree, focaltaxa = humd_sp, nsim = 500, groups = humd_group, col = scico(n = 10, palette = "roma"))
plotCt(head.arid.humd.sig, phy = anilios_tree, focaltaxa = humd_sp, nsim = 500, groups = humd_group, col = scico(n = 10, palette = "roma"))

### HUMID INTERPRETATION
# Within humid lineages Ct1 for body shape ratio was -5.4 (P = 0.084) and head shape was -2.5 (P = 0.002),
# also indicating divergence rather than convergence in these traits.

# Negative Ct1 values suggest parallel evolutionary changes rather than convergence for all focal groups. 

# Frequency-based C5 ------------------------------------------------------

### BodyShapePC because have to use less variable than tips
# c5_for_arid.arid <- convevol::convnumsig(phy = anilios_tree, phendata = BodyShapePC, convtips = arid_sp, nsim = 500, plot = F)
c5_for_arid.dune <- convevol::convnumsig(phy = anilios_tree, phendata = BodyShapePC, convtips = dune_sp, nsim = 500, plot = F)
c5_for_arid.land <- convevol::convnumsig(phy = anilios_tree, phendata = BodyShapePC, convtips = land_sp, nsim = 500, plot = F)
c5_for_arid.humd <- convevol::convnumsig(phy = anilios_tree, phendata = BodyShapePC, convtips = humd_sp, nsim = 500, plot = F)

c5_for_arid.ston <- convevol::convnumsig(phy = anilios_tree, phendata = BodyShapePC, convtips = ston_sp, nsim = 500, plot = F)
c5_for_arid.shie <- convevol::convnumsig(phy = anilios_tree, phendata = BodyShapePC, convtips = shie_sp, nsim = 500, plot = F)

# C5.head.arid.arid <- convevol::convnumsig(phy = anilios_tree, phendata = HeadLMsM, convtips = arid_sp, nsim = 500)
C5.head.arid.dune <- convevol::convnumsig(phy = anilios_tree, phendata = HeadLMsM, convtips = dune_sp, nsim = 500, plot = F)
C5.head.arid.land <- convevol::convnumsig(phy = anilios_tree, phendata = HeadLMsM, convtips = land_sp, nsim = 500, plot = F)
C5.head.arid.humd <- convevol::convnumsig(phy = anilios_tree, phendata = HeadLMsM, convtips = humd_sp, nsim = 500)

C5.head.arid.ston <- convevol::convnumsig(phy = anilios_tree, phendata = HeadLMsM, convtips = ston_sp, nsim = 500, plot = F)
C5.head.arid.shie <- convevol::convnumsig(phy = anilios_tree, phendata = HeadLMsM, convtips = shie_sp, nsim = 500, plot = F)

# Interpretation for C5
# No convergence based on frequency-based measure of convergence. 

c5_for_arid.dune[[1]]; mean(c5_for_arid.dune[[2]])
c5_for_arid.land[[1]]; mean(c5_for_arid.land[[2]])
c5_for_arid.ston[[1]]; mean(c5_for_arid.ston[[2]])
c5_for_arid.shie[[1]]; mean(c5_for_arid.shie[[2]])


C5.head.arid.dune[[1]]; mean(C5.head.arid.dune[[2]])
C5.head.arid.land[[1]]; mean(C5.head.arid.land[[2]])
C5.head.arid.ston[[1]]; mean(C5.head.arid.ston[[2]])
C5.head.arid.shie[[1]]; mean(C5.head.arid.shie[[2]])


# Humid species some convergence? 
c5_for_arid.humd[[1]]; mean(c5_for_arid.humd[[2]])
C5.head.arid.humd[[1]]
mean(C5.head.arid.humd[[2]])


save(head.arid.arid.sig, head.arid.humd.sig, body.arid.arid.sig, 
     body.arid.humd.sig, head.arid.dune.sig, head.arid.land.sig, head.arid.ston.sig, head.arid.shie.sig, 
     body.arid.dune.sig, body.arid.land.sig, body.arid.ston.sig, body.arid.shie.sig,
     c5_for_arid.arid, c5_for_arid.dune, c5_for_arid.land, c5_for_arid.humd, 
     C5.head.arid.arid, C5.head.arid.dune, C5.head.arid.land, C5.head.arid.humd,
     file = 'data/script_generated_data/convergence_ct_score.Rdata')


# RESULTS TABLE -----------------------------------------------------------

results_matrix <- rbind(body.arid.dune.sig$pvals[,1], 
                        body.arid.dune.sig$pvals[,2], 
                        body.arid.land.sig$pvals[,1], 
                        body.arid.land.sig$pvals[,2],
                        body.arid.ston.sig$pvals[,1], 
                        body.arid.ston.sig$pvals[,2],
                        body.arid.shie.sig$pvals[,1], 
                        body.arid.shie.sig$pvals[,2],
                        body.arid.humd.sig$pvals[,1], 
                        body.arid.humd.sig$pvals[,2], 
                        head.arid.dune.sig$pvals[,1], 
                        head.arid.dune.sig$pvals[,2], 
                        head.arid.land.sig$pvals[,1], 
                        head.arid.land.sig$pvals[,2], 
                        head.arid.ston.sig$pvals[,1], 
                        head.arid.ston.sig$pvals[,2], 
                        head.arid.shie.sig$pvals[,1], 
                        head.arid.shie.sig$pvals[,2], 
                        head.arid.humd.sig$pvals[,1],
                        head.arid.humd.sig$pvals[,2]
                        )


results_table <- as.data.frame(round(results_matrix, 3))

rownames(results_table) <- c("Body_dune", "Body_duneP", 
                             "Body_land", "Body_landP", 
                             "Body_stony", "Body_stonyP", 
                             "Body_shield", "Body_shieldP", 
                             "Body_humid", "Body_humidP",
                             "head_dune", "head_duneP", 
                             "head_land", "head_landP", 
                             "head_stony", "head_stonyP", 
                             "head_shield", "head_shieldP", 
                             "head_humid", "head_humidP")

### C5 results
results_table_c5 <- data.frame(C5 = c(mean(c5_for_arid.dune[[2]]),
                                      mean(c5_for_arid.land[[2]]),
                                      mean(c5_for_arid.ston[[2]]),
                                      mean(c5_for_arid.shie[[2]]),
                                      mean(c5_for_arid.humd[[2]]),
                                      mean(C5.head.arid.dune[[2]]), 
                                      mean(C5.head.arid.land[[2]]), 
                                      mean(C5.head.arid.ston[[2]]), 
                                      mean(C5.head.arid.shie[[2]]), 
                                      mean(C5.head.arid.humd[[2]])),
                               P_value = c(c5_for_arid.dune[[1]],
                                           c5_for_arid.land[[1]], 
                                           c5_for_arid.ston[[1]], 
                                           c5_for_arid.shie[[1]], 
                                           c5_for_arid.humd[[1]],
                                           C5.head.arid.dune[[1]], 
                                           C5.head.arid.land[[1]], 
                                           C5.head.arid.ston[[1]], 
                                           C5.head.arid.shie[[1]], 
                                           C5.head.arid.humd[[1]]))
                                 
rownames(results_table_c5) <- c("Body_dune", "Body_land", "Body_stone", "Body_shield", "Body_humid",
                                "head_dune", "head_land", "head_stone", "head_shield", "head_humid")

# write.csv(results_table, file = "manuscript_word/supplements/supp_c_measures.csv", row.names = T)
# write.csv(results_table_c5, file = "manuscript_word/supplements/supp_c_measures_c5.csv", row.names = T)



# PLOTTING ----------------------------------------------------------------

# Supplementary figures

library(viridis)

source('code/utility/func_plotctnew.R')
dev.off()
pdf(file = "output/body_c_measures.pdf", height = 11, width = 8)
layout.matrix <- matrix(c(1:15), nrow = 5, ncol = 3)

layout(mat = layout.matrix, 
       heights = c(1,1,1,1,1),
       widths = c(1,2,2))
plotCt.phylo(body.arid.dune.sig, phy = anilios_tree, focaltaxa = dune_sp)
plotCt.phylo(body.arid.land.sig, phy = anilios_tree, focaltaxa = land_sp)
plotCt.phylo(body.arid.ston.sig, phy = anilios_tree, focaltaxa = ston_sp)
plotCt.phylo(body.arid.shie.sig, phy = anilios_tree, focaltaxa = shie_sp)
plotCt.phylo(body.arid.humd.sig, phy = anilios_tree, focaltaxa = humd_sp)

plotCt.new(body.arid.dune.sig, phy = anilios_tree, focaltaxa = dune_sp, nsim = 500, col = "black", .ylabel = "Body shape")
plotCt.new(body.arid.land.sig, phy = anilios_tree, focaltaxa = land_sp, nsim = 500, col = "black", .ylabel = "Body shape")
plotCt.new(body.arid.ston.sig, phy = anilios_tree, focaltaxa = ston_sp, nsim = 500, .ylabel = "Body shape")
plotCt.new(body.arid.shie.sig, phy = anilios_tree, focaltaxa = shie_sp, nsim = 500, .ylabel = "Body shape")
plotCt.new(body.arid.humd.sig, phy = anilios_tree, focaltaxa = humd_sp, nsim = 500, .ylabel = "Body shape")

plotCt.new(head.arid.dune.sig, phy = anilios_tree, focaltaxa = dune_sp, groups = NULL, nsim = 500, .ylabel = "Snout shape")
plotCt.new(head.arid.land.sig, phy = anilios_tree, focaltaxa = land_sp, groups = NULL, nsim = 500, .ylabel = "Snout shape")
plotCt.new(head.arid.ston.sig, phy = anilios_tree, focaltaxa = ston_sp, nsim = 500, .ylabel = "Snout shape")
plotCt.new(head.arid.shie.sig, phy = anilios_tree, focaltaxa = shie_sp, nsim = 500, .ylabel = "Snout shape")
plotCt.new(head.arid.humd.sig, phy = anilios_tree, focaltaxa = humd_sp, nsim = 500, .ylabel = "Snout shape")

dev.off()

# pdf(file = "output/head_c_measures.pdf", height = 6, width = 11)
plotCt(head.arid.dune.sig, phy = anilios_tree, focaltaxa = dune_sp, groups = dune_group, nsim = 500, col = viridis(n = 36, option = "A"))
plotCt(head.arid.land.sig, phy = anilios_tree, focaltaxa = land_sp, groups = uplan_group, nsim = 500, col = viridis(n = 36, option = "A"))
plotCt(head.arid.ston.sig, phy = anilios_tree, focaltaxa = ston_sp, nsim = 500)
plotCt(head.arid.shie.sig, phy = anilios_tree, focaltaxa = shie_sp, nsim = 500)
plotCt(head.arid.humd.sig, phy = anilios_tree, focaltaxa = humd_sp, nsim = 500)
# plotCt(head.arid.arid.sig, phy = anilios_tree, focaltaxa = arid_sp, groups = arid_group, nsim = 500, col = viridis(n = 36, option = "A"))
# plotCt(head.arid.humd.sig, phy = anilios_tree, focaltaxa = humd_sp, nsim = 500, groups = humd_group, col = viridis(n = 10, option = "D"))
dev.off()
