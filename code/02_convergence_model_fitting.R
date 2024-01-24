# 02_convergence_model_fitting.R
# Sarin Tiatragul 
# October 2023

# Results can be loaded from 
load(file = 'data/script_generated_data/06_convergence_model_fitting.rda')

# libraries ---------------------------------------------------------------
# Load packages with two lines
pkgs = c("phytools", "ape", "geomorph", "ggplot2", "dplyr", "stringr", "corHMM", "geiger", "nlme", "convevol", "geiger", "scico") # package names
inst = lapply(pkgs, library, character.only = TRUE) # load them
library(nlme)
'%notin%' <- Negate('%in%')

# data --------------------------------------------------------------------
load('data/script_generated_data/subset_mcmctree_shape.rda')
# Load two D array of y symmetric component of dorsal head shape
load('data/script_generated_data/dorsal_head_shape.rda')

blindsnake_data <- read.csv(file = 'data/script_generated_data/mean_blindsnake_data.csv', row.names = 1)

# Drop outgroups from analyses 
anilios_tree <- ape::drop.tip(phy = sub_phy, tip = c("Ramphotyphlops_multilineatus", "Acutotyphlops_subocularis"))
anilios_data <- blindsnake_data[rownames(blindsnake_data) %in% anilios_tree$tip.label,]

# Change level names
anilios_data$mean_bulk_cat <- dplyr::recode_factor(.x = anilios_data$mean_bulk_cat, "1" = "soft", "2" = "med", "3" = "hard")

# Match order of tips and data
anilios_data <- anilios_data[match(anilios_tree$tip.label, anilios_data$species), ]

# Prepare bodyshape traits data (same with 04_model_fitting.R for mvMORPH)
BodyShape <- as.matrix(anilios_data[, c("sh_mbd", "sh_mtw", "sh_hwe", "sh_hda", "sh_tl")])
BodyShapePC <- as.matrix(anilios_data[, c("fil_PC1", "fil_PC2")])

rownames(sp_means_dorsal) <- gsub(x = rownames(sp_means_dorsal), pattern = "^", replacement = "Anilios_")
means_2d_dorsal <- sp_means_dorsal[rownames(sp_means_dorsal) %in% anilios_tree$tip.label,]
means_2d_dorsal <- means_2d_dorsal[match(anilios_tree$tip.label, rownames(means_2d_dorsal)), ]

# Turn head landmark traits 
Y_head <- as.matrix(means_2d_dorsal[-c(1,2)])
colnames(Y_head) <- c(paste0("comp", 1:52))

Y_head[, 1:10]

HeadLMsM <- as.matrix(anilios_data[, c("mean_PC1", "mean_PC2")])

# stricter category of arid based on me looking at satellite images of where most the locality came from and map from Mabbutt 1988
arid_strict <- read.csv('data/blindsnake_arid_strict.csv', row.names = 2)

### Make simmap 
state <- anilios_data[anilios_tree$tip.label,]$arid_category
names(state) <- anilios_tree$tip.label
state[state == "Arid"]

plot(anilios_tree)
nodelabels()
tiplabels()
which(anilios_tree$tip.label %in% names(state[state == "Arid"]))

names(state[state == "Arid"])
painted_tree <- paintSubTree(anilios_tree, node = c(70), state = "arid", anc.state = "unspecified")
painted_tree <- paintSubTree(painted_tree, node = c(66), state = "arid", anc.state = "unspecified")
painted_tree <- paintSubTree(painted_tree, node = 45, state = "arid", anc.state = "unspecified")
plot(painted_tree)

painted_tree <- paintBranches(painted_tree, edge = c(8, 10, 11), state = "arid", anc.state = "unspecified") # bicolor, waitii, centralis
painted_tree <- paintBranches(painted_tree, edge = c(13, 14, 15, 24, 27, 28, 36), state = "arid", anc.state = "unspecified")
painted_tree <- paintBranches(painted_tree, edge = 33, state = "unspecified", anc.state = "unspecified")
names(state[state == "Humid"])
painted_tree <- paintBranches(painted_tree, edge = c(1, 2, 19, 20, 23), state = "humid", anc.state = "unspecified") # "humid species"

# Colour the tree
col <- scico(n = 3, palette = "hawaii", categorical = TRUE); names(col)<-c("arid","humid", "unspecified")
scicol_paintedtree <- plotSimmap(painted_tree, col, fsize=0.6, node.numbers=FALSE,lwd=3, pts=FALSE)


# Another tree ------------------------------------------------------------
### Make simmap 2
state <- arid_strict[anilios_tree$tip.label,]$conv_type
names(state) <- anilios_tree$tip.label
state[state == "sand"]; state[state == "humid"]; state[state == "upland"]; state[state == "other"]

plot(anilios_tree)
nodelabels()
tiplabels()
which(anilios_tree$tip.label %in% names(state[state == "upland"]))


painted_tree2 <- paintBranches(anilios_tree, edge = which(anilios_tree$tip.label %in% names(state[state == "sand"])), 
                               state = "sand", anc.state = "unspecified") # endoterus, margaretae, bituberculatus
painted_tree2 <- paintBranches(painted_tree2, edge = which(anilios_tree$tip.label %in% names(state[state == "humid"])), 
                               state = "humid", anc.state = "unspecified")
painted_tree2 <- paintBranches(painted_tree2, edge = which(anilios_tree$tip.label %in% names(state[state == "upland"])),
                               state = "upland", anc.state = "unspecified")
painted_tree2 <- paintBranches(painted_tree2, edge = which(anilios_tree$tip.label %in% names(state[state == "other"])),
                               state = "unspecified", anc.state = "unspecified")
# painted_tree2 <- paintBranches(painted_tree2, edge = c(1, 2, 19, 20, 23), state = "humid", anc.state = "unspecified") # "humid species"

plot(painted_tree2)

# Colour the tree
col <- viridis::viridis(4, option = "D"); names(col)<-c("sand","humid", "upland", "unspecified")
scicol_paintedtree2 <- plotSimmap(painted_tree2, col, fsize=0.6, node.numbers=FALSE,lwd=3, pts=FALSE)

# Single rate models ------------------------------------------------------

# BM1 analysis with unique rate matrix
mvBM1 <- mvMORPH::mvBM(tree = anilios_tree, data = BodyShape, model = "BM1", method = "rpf")
mvBM1.head <- mvMORPH::mvBM(tree = anilios_tree, data = HeadLMsM, model = "BM1", method = "rpf")

# EB model analysis 
mvEB <- mvMORPH::mvEB(tree = anilios_tree, data = BodyShape, method = "rpf", param=list(up=0, low=-10))
mvEB.head <- mvMORPH::mvEB(tree = anilios_tree, data = HeadLMsM, method = "rpf", param=list(up=0, low=-10))

# OU model analysis with a unique optimum
mvOU1 <- mvMORPH::mvOU(tree = painted_tree, data = BodyShape, model = "OU1", method = "rpf")
mvOU1.conv <- mvMORPH::mvOU(tree = painted_tree2, data = BodyShape, model = "OU1", method = "rpf")
mvOU1.head <- mvMORPH::mvOU(tree = anilios_tree, data = HeadLMsM, model = "OU1", method = "rpf")

# Multi-regime ------------------------------------------------------------
### Allow arid and humid to exhibit different trait optima
### Allow distinct trait optimum
### Support for this would provide evidence of convergence by indicating selective forces are driving 'arid' lineages to shared adaptive peak. 

# mvOU2 <- mvMORPH::mvOU(tree = painted_tree, data = BodyShape, model = "OUM", method = "pseudoinverse") ## takes longer...
mvOU2 <- mvMORPH::mvOU(tree = painted_tree, data = BodyShape, model = "OUM", method = "rpf") ## using less variables result in more "reliable' solution
mvOU2.head <- mvMORPH::mvOU(tree = painted_tree, data = HeadLMsM, model = "OUM", method = "rpf")

# BMM multi regime 
mvBMM <-  mvMORPH::mvBM(tree = painted_tree, data = BodyShape, model = "BMM")
mvBMM.head <-  mvMORPH::mvBM(tree = painted_tree, data = HeadLMsM, model = "BMM")


# mvOU2 with multiple categories of aridity 
mvOU2.conv <- mvMORPH::mvOU(tree = painted_tree2, data = BodyShape, model = "OUM", method = "rpf") ## using less variables result in more "reliable' solution
mvOU2.conv.head <- mvMORPH::mvOU(tree = painted_tree2, data = HeadLMsM, model = "OUM", method = "rpf")

# BMM multi regime 
mvBMM.conv <-  mvMORPH::mvBM(tree = painted_tree2, data = BodyShape, model = "BMM")
mvBMM.conv.head <-  mvMORPH::mvBM(tree = painted_tree2, data = HeadLMsM, model = "BMM")

# AICc compare ------------------------------------------------------------

results_table1 <- mvMORPH::aicw(list(mvBM1, mvEB, mvOU1, mvOU2, mvBMM), aicc = TRUE)
results_table2 <- mvMORPH::aicw(list(mvBM1.head, mvEB.head, mvOU1.head, mvOU2.head, mvBMM.head), aicc = TRUE)

results_table1.conv <- mvMORPH::aicw(list(mvBM1, mvEB, mvOU1, mvOU2.conv, mvBMM.conv), aicc = TRUE)
results_table2.conv <- mvMORPH::aicw(list(mvBM1.head, mvEB.head, mvOU1.head, mvOU2.conv.head, mvBMM.conv.head), aicc = TRUE)

### Interpretation of the result:
# OUM does not provide better fit to the data than other evolutionary models. 

# Write results table -----------------------------------------------------
results_body <- data.frame(models=results_table1.conv$models,
                           AICc=results_table1.conv$AIC,
                           deltaaic=results_table1.conv$diff,
                           weights=results_table1.conv$wi,
                           aicw=results_table1.conv$aicweights,
                           trait = "body")

results_head <- data.frame(models=results_table2.conv$models,
                           AICc=results_table2.conv$AIC,
                           deltaaic=results_table2.conv$diff,
                           weights=results_table2.conv$wi,
                           aicw=results_table2.conv$aicweights,
                           trait = "head")

write.csv(results_body, file = "manuscript_word/supplements/mvMORPH-modelfit-body.csv", row.names = FALSE)
write.csv(results_head, file = "manuscript_word/supplements/mvMORPH-modelfit-head.csv", row.names = FALSE)