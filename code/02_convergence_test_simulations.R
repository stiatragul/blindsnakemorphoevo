# 02_convergence_test_simulations.R
# Sarin Tiatragul 
# April 2024
# Simulations to test the power of the Ct-metrics to detect convergence under different scenarios
# For more information, check out the online supplementary material for Tiatragul et al. 2024 (Evolution)

# libraries ---------------------------------------------------------------
# Packages for analyses
library(phytools); library(ape); 
library(mvMORPH); library(geiger)
library(R.utils); 
library(convevol)
library(dplyr)

source("code/utility/func_rescaleTree_phytools_2012.R")
'%notin%' <- Negate('%in%')

# pre-release code from William Brightly for convevol v2.1+ 
source("code/utility/func_pwCheck_2.1.R") # 
source("code/utility/convSigCt_2.2.R")
source("code/utility/calcConvCt_2.2.R")

# data --------------------------------------------------------------------

# Load in our observed tree full tree, subsetted tree, and subset 
load('data/script_generated_data/subset_mcmctree_shape.rda')
anilios_tree <- ape::drop.tip(sub_phy, tip = c("Ramphotyphlops_multilineatus","Acutotyphlops_subocularis"))

# Load two D array of y symmetric component of dorsal head shape
load('data/script_generated_data/dorsal_head_shape.rda')

# Load data frame of summarised traits
anilios_data <- read.csv(file = 'data/script_generated_data/anilios_summary_data.csv', row.names = 1)

# Load in dataframe with discrete habitat characters 
arid_strict <- read.csv('data/blindsnake_arid_strict.csv', row.names = 2)
arid_sp_strict <- rownames(arid_strict)[which(arid_strict$aridity_strict == "arid")]

## We are also comparing between species that are in different desert types (with dune) according to Mabbutt 1988
arid_uplan <- rownames(arid_strict)[which(arid_strict$upland == "upland")]
arid_dunes <- rownames(arid_strict)[which(arid_strict$sand == "sand")]
arid_stony <- rownames(arid_strict)[which(arid_strict$stony == "stony")]
arid_shied <- rownames(arid_strict)[which(arid_strict$shield == "shield")]
humd_sp <- rownames(anilios_data)[which(anilios_data$arid_category == "Humid")]

# Paint observed tree -----------------------------------------------------

## To fit mvMORPH::mvOU using OUM (multi-peak OU model), we need a painted simmap
## We paint the tree using phytools::paintSubTree manually at nodes leading to "arid" species
## We specify ancestral state is unspecified (other)

painted_tree <- phytools::paintSubTree(anilios_tree, node = c(70), state = "arid", anc.state = "unspecified")
painted_tree <- phytools::paintSubTree(painted_tree, node = c(66), state = "arid", anc.state = "unspecified")
painted_tree <- phytools::paintSubTree(painted_tree, node = 45, state = "arid", anc.state = "unspecified")

# Paint additional branches
painted_tree <- phytools::paintBranches(painted_tree, edge = c(8, 10, 11), state = "arid", anc.state = "unspecified") # bicolor, waitii, centralis
painted_tree <- phytools::paintBranches(painted_tree, edge = c(13, 14, 15, 24, 27, 28, 36), state = "arid", anc.state = "unspecified")
painted_tree <- phytools::paintBranches(painted_tree, edge = 33, state = "unspecified", anc.state = "unspecified")

plot(painted_tree)

# Body traits -------------------------------------------------------------

## prep body traits
Y_body <- as.matrix(anilios_data[, c("sh_mbd", "sh_mtw", "sh_hwe", "sh_hda", "sh_tl")])
body_sh_traits <- list(Y_body=Y_body) # prepared as list for mvMORPH


### Visualise the empirical traits
# Using PCA
PC_body <- prcomp(Y_body, scale = TRUE, center = TRUE)

phytools::phylomorphospace(tree = painted_tree, X = PC_body$x[,1:2], bty = "n")
ss<-getStates(painted_tree,"tips")
colors<-setNames(c("black","red"),c("arid","unspecified"))
barcols<-setNames(sapply(ss,function(x,y) y[which(names(y)==x)],y=colors),names(ss))

# Your list of colors
barcols <- rev(c("red", "red", "black", "black", "black", "black", "black", "red", "black", 
                 "black", "black", "red", "black", "black", "black", "red", "red", "red", 
                 "red", "red", "red", "red", "red", "red", "black", "red", "black", "black", "red", 
                 "black", "red", "red", "red", "red", "black", "red", "black"))

plotTree.barplot(painted_tree, PC_body$x[,1][painted_tree$tip.label],
                 args.plotTree=list(),
                 args.barplot=list(col=barcols,xlab="Body PC1"))
legend(x="bottomright",legend=names(colors),pch=22,pt.cex=2,pt.bg=colors, box.col="transparent")


# Step 1 ------------------------------------------------------------------

## Fit model to observed data

fit_phy_BM1 <- mvMORPH::mvBM(Y_body, tree = painted_tree, model = "BM1", method = "rpf", scale.height = TRUE)
fit_phy_BMM <- mvMORPH::mvBM(Y_body, tree = painted_tree, model = "BMM", method = "rpf", scale.height = TRUE)
fit_phy_OUM <- mvMORPH::mvOU(Y_body, tree = painted_tree, model = "OUM", method = "rpf", scale.height = TRUE)
fit_phy_OU1 <- mvMORPH::mvOU(Y_body, tree = painted_tree, model = "OU1", method = "rpf", scale.height = TRUE)
fit_phy_EB <- mvMORPH::mvEB(Y_body, tree = painted_tree, method = "rpf", scale.height = TRUE, param=list(up=0, low=-10))

fit_phy_AIC_table <- mvMORPH::aicw(list(fit_phy_BM1, fit_phy_EB, fit_phy_OU1, fit_phy_OUM, fit_phy_BMM), aicc = TRUE)
fit_phy_AIC_table

### In our observed data, BM clearly and other models outperform OUM

# Step 2 ------------------------------------------------------------------

## Simulate an ideal situation where arid taxa are convergent
## By making sure that the two peaks are very different optima

# theta (θ) determines the optimum trait values 
# Very large optima suggest species tend to evolve the same traits given their discrete character (e.g. 'habitat').
theta <-rep(c(1,2), times = 5) # number of trait(s) x discrete character state(s), in this case 2 states and 5 so length of vector is 10

# sigma (σ^2) determines drift or rate of phenotypic evolution under a Brownian motion model of evolution
# In this case we will use the sigma estimated on the empirical data advised by Julien Clavel 2024-05-03. 
# Or simulate can make the rate constant at 0.1 for all traits,
# sigma <- matrix(c(rep(0.1,25)), nrow = 5)
# diag(sigma) <- 1
sigma <- fit_phy_OUM$sigma

# alpha (α) is strength of constraints or rate of adaptation. 
# Larger values of alpha indicate that the character is pulled more strongly toward θ ("rubber band")
# alpha<-diag(1,nrow=5)
alpha <- fit_phy_OUM$alpha

# Simulate multivariate continuous traits on simulated phylogeny according to OUM model
sim_Y_OUM <- mvMORPH::mvSIM(tree = painted_tree, nsim=1, model = "OUM", 
                            param = list(theta = theta, sigma=sigma, ntraits=5, alpha=alpha, root=FALSE))

### Visualise the simulated traits under OU

# Using PCA
PC_sim <- prcomp(sim_Y_OUM, scale = TRUE, center = TRUE)
phytools::phylomorphospace(tree = painted_tree, X = PC_sim$x[,1:2], bty = "n")
dev.off()
plotTree.barplot(painted_tree, PC_sim$x[,1], args.barplot=list(col=barcols,xlab="Body PC1"))
legend(x="topright",legend=names(colors),pch=22,pt.cex=2,pt.bg=colors, box.col="transparent")

## Fit model to simulated data
fit_sim_BM1 <- mvMORPH::mvBM(sim_Y_OUM, tree = painted_tree, model = "BM1", method = "rpf", scale.height = TRUE)
fit_sim_BMM <- mvMORPH::mvBM(sim_Y_OUM, tree = painted_tree, model = "BMM", method = "rpf", scale.height = TRUE)
fit_sim_OUM <- mvMORPH::mvOU(sim_Y_OUM, tree = painted_tree, model = "OUM", method = "rpf", scale.height = TRUE, param=list(decomp="diagonal"))
fit_sim_OU1 <- mvMORPH::mvOU(sim_Y_OUM, tree = painted_tree, model = "OU1", method = "rpf", scale.height = TRUE, param=list(decomp="diagonal"))
fit_sim_EB <- mvMORPH::mvEB(sim_Y_OUM, tree = painted_tree, method = "rpf", scale.height = TRUE, param=list(up=0, low=-10))

fit_sim_AIC_table <- mvMORPH::aicw(list(fit_sim_BM1, fit_sim_EB, fit_sim_OU1, fit_sim_OUM, fit_sim_BMM), aicc = TRUE)
fit_sim_AIC_table
## As expected, OUM model outperforms BM and other models because the dataset is simulated under OUM


## Functiont to fit different theta to show that strength is important

mv_fitter <- function(theta_1, theta_2){
  
  theta <- rep(c(theta_1, theta_2), times = 5)
  sigma <- fit_phy_OUM$sigma
  alpha <- fit_phy_OUM$alpha
  sim_Y_OUM <- mvMORPH::mvSIM(tree = painted_tree, nsim=1, model = "OUM", param = list(theta = theta, sigma=sigma, alpha=alpha, ntraits=5, root=FALSE))
  
  fit_sim_BM1 <- mvMORPH::mvBM(sim_Y_OUM, tree=painted_tree, model="BM1", method="rpf", scale.height=TRUE)
  fit_sim_BMM <- mvMORPH::mvBM(sim_Y_OUM, tree=painted_tree, model="BMM", method="rpf", scale.height=TRUE)
  fit_sim_OUM <- mvMORPH::mvOU(sim_Y_OUM, tree=painted_tree, model="OUM", method="rpf", scale.height=TRUE, param=list(decomp="diagonal"))
  fit_sim_OU1 <- mvMORPH::mvOU(sim_Y_OUM, tree=painted_tree, model="OU1", method="rpf", scale.height=TRUE, param=list(decomp="diagonal"))
  fit_sim_EB <- mvMORPH::mvEB(sim_Y_OUM, tree=painted_tree, method="rpf", scale.height = TRUE, param=list(up=0, low=-10))
  
  fit_sim_AIC_table <- mvMORPH::aicw(list(fit_sim_BM1, fit_sim_EB, fit_sim_OU1, fit_sim_OUM, fit_sim_BMM), aicc = TRUE)
  fit_sim_AIC_table
  
  PC_sim <- prcomp(sim_Y_OUM, scale = TRUE, center = TRUE)
  phytools::plotTree.barplot(painted_tree, PC_sim$x[,1], args.barplot=list(col=barcols,xlab="Body PC1"))
  
  return(list(fit_sim_AIC_table, sim_Y_OUM, PC_sim))
}

# Simulate under different theta (optima values)

sim_AIC_0 <- mv_fitter(NULL,NULL)
sim_AIC_1 <- mv_fitter(1,1)
sim_AIC_2 <- mv_fitter(1,2)
# sim_AIC_3 <- mv_fitter(1,3)
# sim_AIC_4 <- mv_fitter(1,4)
sim_AIC_5 <- mv_fitter(1,5)
sim_AIC_7 <- mv_fitter(1,7)
sim_AIC_9 <- mv_fitter(1,9)

sim_AIC_0[[1]]
sim_AIC_1[[1]]
sim_AIC_2[[1]]
sim_AIC_5[[1]]
sim_AIC_9[[1]]

### Make csv tables for supp

AIC_tablr <- function(.table){
  .table <- .table
  tmp <- data.frame(models=.table$models, AICc=.table$AIC, deltaaic=.table$diff,weights=.table$wi, aicw=.table$aicweights)
  tmp <-  tmp %>% dplyr::arrange(desc(weights))
  return(tmp)
}

write.csv(AIC_tablr(sim_AIC_1[[1]]), file = "output/supp_sim_conv_tab1.csv", row.names = F)
write.csv(AIC_tablr(sim_AIC_2[[1]]), file = "output/supp_sim_conv_tab2.csv", row.names = F)
write.csv(AIC_tablr(sim_AIC_5[[1]]), file = "output/supp_sim_conv_tab5.csv", row.names = F)
write.csv(AIC_tablr(sim_AIC_9[[1]]), file = "output/supp_sim_conv_tab9.csv", row.names = F)

## Make figure showing different simulated traits for different theta

stree <- painted_tree
stree$tip.label <- gsub(pattern = "Anilios_", replacement = "", x=stree$tip.label)
rownames(PC_body$x) <- gsub(pattern = "Anilios_", replacement = "", x=rownames(PC_body$x))

pdf(file = "output/phylo_barplot_sim_traits.pdf", height = 8, width = 11.33)
par(mfrow=c(1,6))
plotTree.barplot(stree, PC_body$x[,1], add=TRUE, args.barplot=list(col=barcols,xlab="Body PC1"),mar=c(5.1,0,2.1,2.1))
plotTree.barplot(painted_tree, sim_AIC_1[[3]]$x[,1], args.barplot=list(col=barcols,xlab="Sim.PC1 theta(1,1)",xlim=c(-4,4)),args.plotTree=list(plot=F),add=T)
plotTree.barplot(painted_tree, sim_AIC_2[[3]]$x[,1], args.barplot=list(col=barcols,xlab="Sim.PC1 theta(1,2)",xlim=c(-4,4)),args.plotTree=list(plot=F),add=T)
plotTree.barplot(painted_tree, sim_AIC_5[[3]]$x[,1], args.barplot=list(col=barcols,xlab="Sim.PC1 theta(1,5)",xlim=c(-4,4)),args.plotTree=list(plot=F),add=T)
plotTree.barplot(painted_tree, sim_AIC_9[[3]]$x[,1], args.barplot=list(col=barcols,xlab="Sim.PC1 theta(1,9)",xlim=c(-4,4)),args.plotTree=list(plot=F),add=T)
dev.off()


# Using convevol ----------------------------------------------------------

plot(painted_tree)

# Choose tips that are distantly related and suspect to converge
focal_tx9 <- c("Anilios_grypusW", "Anilios_systenos", "Anilios_endoterus", "Anilios_aspina",
               "Anilios_centralis", "Anilios_bicolor", "Anilios_ammodytes", "Anilios_grypusET", "Anilios_bituberculatus")

focal_tx7 <- c("Anilios_grypusW", "Anilios_endoterus", "Anilios_centralis", "Anilios_bicolor", "Anilios_ammodytes", 
               "Anilios_grypusET", "Anilios_bituberculatus")

focal_tx <- c("Anilios_grypusW", "Anilios_grypusET", "Anilios_ammodytes", "Anilios_centralis")

## Convergence as determined by Ct-metrics are based primarily on the focal taxa and their relationships with one another. 
## Ct metrics work well when the tree is large because it allows for the 'true' D.maxt to be calculated. 

pwCheck(phy = anilios_tree, focaltaxa = focal_tx)

## Fit observed data a few times
obs.sig.ct <- convevol::convSigCt(phy=painted_tree, traits = Y_body, focaltaxa = focal_tx, nsim = 100)
obs.sig.ct$pvals
plotCt(sim.sig.ct.OU, phy = anilios_tree, focaltaxa = focal_tx, nsim = 100, size = 2)

# Do this multiple times to check
obs_ct_fittr <- function(focal_tx, grps=NULL){
  obs.sig.ct <- convevol::convSigCt(phy=painted_tree, traits = Y_body, focaltaxa = focal_tx, nsim = 100, groups=grps)
  tmp <- as.data.frame(obs.sig.ct$pvals)
  result <- data.frame(Ct1 = tmp$V1[1], Ct2 = tmp$V1[2], Ct3 = tmp$V1[3], Ct4 = tmp$V1[4], p_value = tmp$pvals[1])
  return(list(result, obs.sig.ct))
}

# Replicate fit for these
arid_uplan
arid_dunes
arid_stony
arid_shied
humd_sp

## Grouping
uplan_group <- rep("G1", length(arid_uplan))
names(uplan_group) <- arid_uplan
uplan_group[c(1,2,5:8)] <- c("gryw", "lep", "pil", "end", "gan", "amm")
uplan_group[which(names(uplan_group) %in% c("Anilios_waitii", "Anilios_centralis"))] <- "waicen"

dune_group <- rep("G1", length(arid_dunes))
names(dune_group) <- arid_dunes
dune_group[c(1:5, 8)] <- c("grpw", "bic", "wai", "pil", "end", "div")
dune_group[which(names(dune_group) %in%  c("Anilios_margaretae", "Anilios_bituberculatus"))] <- "marbit"


obs_df_4tx <- replicate(10, obs_ct_fittr(focal_tx), simplify = FALSE)
obs_df_sand <- replicate(10, obs_ct_fittr(arid_dunes, grps=dune_group), simplify = FALSE)
obs_df_upland <- replicate(10, obs_ct_fittr(arid_uplan, grps=uplan_group), simplify = FALSE)
obs_df_stony <- replicate(10, obs_ct_fittr(focal_tx), simplify = FALSE)
obs_df_humid <- replicate(10, obs_ct_fittr(humd_sp), simplify = FALSE)
obs_df_shield <- replicate(10, obs_ct_fittr(arid_shied), simplify = FALSE)

## Ideal case with OU evolution
sim.sig.ct.OU <- convevol::convSigCt(phy=painted_tree, traits = sim_AIC_5[[2]], focaltaxa = focal_tx, nsim = 100)
sim.sig.ct.OU$pvals
plotCt(sim.sig.ct.OU, phy = anilios_tree, focaltaxa = focal_tx, nsim = 100, size = 2)


## When ancestral states are known
# specify user.ace when ancestral states are known (such as the case for simulations)
sim.sig.ct.OU <- convevol::convSigCt(phy=painted_tree, traits = sim_AIC_5[[2]], focaltaxa = focal_tx, nsim = 100)
sim.sig.ct.OU$pvals

## Loop through this multiple times
# 
# result_df <- data.frame(Ct1 = numeric(), Ct2  = numeric(), Ct3  = numeric(), Ct4  = numeric(), p_value = numeric())
# 
# sim_Y_OUM <- mvMORPH::mvSIM(tree = painted_tree, nsim = 100, model = "OUM", param = list(theta = theta, sigma = sigma, ntraits = 5, alpha = alpha, root = FALSE))
# 
# # Loop through the code
# for (i in 1:1) {  
#   # Perform simulation
#   sim.sig.ct.OU <- convevol::convSigCt(phy = painted_tree, traits = sim_Y_OUM[[i]], focaltaxa = focal_tx, nsim = 50)
#   tmp <- as.data.frame(sim.sig.ct.OU$pvals)
#   
#   # Extract Ct value and p-value from the result
#   Ct1 <- tmp$V1[1]
#   Ct2 <- tmp$V1[2]
#   Ct3 <- tmp$V1[3]
#   Ct4 <- tmp$V1[4]
#   
#   p_value <- tmp$pvals[1]
#   
#   # Append the results to the data frame
#   result_df <- rbind(result_df, data.frame(Ct = Ct, p_value = p_value))
# }

# Define a function to perform a single simulation and extract results
perform_single_simulation <- function(painted_tree, focal_tx, theta, sigma, alpha) {
  # Perform simulation
  sim_Y_OUM <- mvMORPH::mvSIM(tree = painted_tree, nsim=1, model="OUM", param=list(theta=theta, sigma=sigma, ntraits=5, alpha=alpha, root = FALSE))
  
  # Fit to convevol::convSigCt
  sim.sig.ct.OU <- convevol::convSigCt(phy = painted_tree, traits = sim_Y_OUM, focaltaxa = focal_tx, nsim = 100)
  
  # Extract Ct values and p-value from the result
  tmp <- as.data.frame(sim.sig.ct.OU$pvals)
  result <- data.frame(Ct1 = tmp$V1[1], Ct2 = tmp$V1[2], Ct3 = tmp$V1[3], Ct4 = tmp$V1[4], p_value = tmp$pvals[1])
  
  # Return the result data frame
  return(list(result, sim.sig.ct.OU))
}

# Perform the simulations replicate
num_simulations <- 100  # Number of simulations

# sigma (σ^2) determines drift or rate of phenotypic evolution under a Brownian motion model of evolution
# sigma <- matrix(c(rep(0.1,25)), nrow = 5)
# alpha (α) is strength of constraints or rate of adaptation. 
# alpha<-diag(1,nrow=5)

# Specify the sigma and alpha values
sigma <- fit_phy_OUM$sigma
alpha <- fit_phy_OUM$alpha
focal_tx
# theta (θ) determines the optimum trait values 
theta<-rep(c(1,1), times = 5)
simulation_results_1 <- replicate(num_simulations, perform_single_simulation(painted_tree, focal_tx, theta, sigma, alpha), simplify = FALSE)
simulation_results_1_9 <- replicate(1, perform_single_simulation(painted_tree, focal_tx9, theta, sigma, alpha), simplify = FALSE)

theta<-rep(c(1,2), times = 5)
simulation_results_2_4 <- replicate(num_simulations, perform_single_simulation(painted_tree, focal_tx, theta, sigma, alpha), simplify = FALSE)
simulation_results_2_7 <- replicate(1, perform_single_simulation(painted_tree, focal_tx7, theta, sigma, alpha), simplify = FALSE)
simulation_results_2_9 <- replicate(1, perform_single_simulation(painted_tree, focal_tx9, theta, sigma, alpha), simplify = FALSE)

theta<-rep(c(1,5), times = 5)
simulation_results_5 <- replicate(num_simulations, perform_single_simulation(painted_tree, focal_tx, theta, sigma, alpha), simplify = FALSE)
simulation_results_5_7 <- replicate(1, perform_single_simulation(painted_tree, focal_tx7, theta, sigma, alpha), simplify = FALSE)
simulation_results_5_9 <- replicate(1, perform_single_simulation(painted_tree, focal_tx9, theta, sigma, alpha), simplify = FALSE)

theta<-rep(c(1,9), times = 5)
simulation_results_9_4 <- replicate(num_simulations, perform_single_simulation(painted_tree, focal_tx, theta, sigma, alpha), simplify = FALSE)
simulation_results_9_7 <- replicate(num_simulations, perform_single_simulation(painted_tree, focal_tx7, theta, sigma, alpha), simplify = FALSE)
simulation_results_9_9 <- replicate(num_simulations, perform_single_simulation(painted_tree, focal_tx9, theta, sigma, alpha), simplify = FALSE)

# Save the results
sim_df1 <- data.frame(Ct1 = numeric(), Ct2  = numeric(), Ct3  = numeric(), Ct4  = numeric(), p_value = numeric())
for (i in 1:100){ sim_df1 <- rbind(sim_df1, simulation_results_1[[i]][[1]])}
sim_df2 <- data.frame(Ct1 = numeric(), Ct2  = numeric(), Ct3  = numeric(), Ct4  = numeric(), p_value = numeric())
for (i in 1:100){ sim_df2 <- rbind(sim_df2, simulation_results_2[[i]][[1]])}
sim_df5 <- data.frame(Ct1 = numeric(), Ct2  = numeric(), Ct3  = numeric(), Ct4  = numeric(), p_value = numeric())
for (i in 1:100){ sim_df5 <- rbind(sim_df5, simulation_results_5[[i]][[1]])}
# sim_df7 <- data.frame(Ct1 = numeric(), Ct2  = numeric(), Ct3  = numeric(), Ct4  = numeric(), p_value = numeric())
# for (i in 1:100){ sim_df7 <- rbind(sim_df7, simulation_results_5[[i]][[1]])}
sim_df9 <- data.frame(Ct1 = numeric(), Ct2  = numeric(), Ct3  = numeric(), Ct4  = numeric(), p_value = numeric())
for (i in 1:100){ sim_df9 <- rbind(sim_df9, simulation_results_9[[i]][[1]])}

# obs_df_sand_df <- data.frame(Ct1 = numeric(), Ct2  = numeric(), Ct3  = numeric(), Ct4  = numeric(), p_value = numeric())
# for (i in 1:10){ obs_df_sand_df <- rbind(obs_df_sand_df, obs_df_sand[[i]][[1]])}

save(sim_df1, sim_df2, sim_df5, file = "data/script_generated_data/simulation_Ct_fromFit.Rdata")
save(sim_df1_9, sim_df2_9, sim_df5_9, simulation_results_1_9, simulation_results_2_9, simulation_results_5_9, file = "data/script_generated_data/simulation_Ct_fromFit_9taxa.Rdata")
save(sim_df1_7, sim_df2_7, sim_df5_7, simulation_results_1_7, simulation_results_2_7, simulation_results_5_7, file = "data/script_generated_data/simulation_Ct_fromFit_7taxa.Rdata")
save(sim_df1_4, sim_df2_4, sim_df5_4, simulation_results_1_4, simulation_results_2_4, simulation_results_5_4, file = "data/script_generated_data/simulation_Ct_fromFit_4taxa.Rdata")
save(sim_df9_4, sim_df9_7, sim_df9_9, simulation_results_9_4, simulation_results_9_7, simulation_results_9_9, file = "data/script_generated_data/simulation_Ct_fromFit_theta19.Rdata")


source("code/utility/func_plotctnew.R")
# plotCt(simulation_results_5[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx, nsim = 100)

dev.off()
pdf(file = "output/plotct_sim_traits.pdf", height = 8.5, width = 11.33)
layout.matrix <- matrix(c(1:15), nrow = 3, ncol = 5)
layout(mat = layout.matrix, heights = c(1,1,1,1,1),widths = c(1,2,2,2,2))
plotCt.phylo(simulation_results_1_9[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx9, nsim = 100)
plotCt.phylo(simulation_results_1_7[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx7, nsim = 100)
plotCt.phylo(simulation_results_1_4[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx, nsim = 100)

# Optima 1,2
plotCt.new(simulation_results_1_9[[1]][[2]], .ylabel="Phenotypic distance", phy = anilios_tree, focaltaxa = focal_tx9, nsim = 100)
plotCt.new(simulation_results_1_7[[1]][[2]], .ylabel="Phenotypic distance", phy = anilios_tree, focaltaxa = focal_tx7, nsim = 100)
plotCt.new(simulation_results_1_4[[1]][[2]], .ylabel="Phenotypic distance", phy = anilios_tree, focaltaxa = focal_tx4, nsim = 100)

# optima 1,2
plotCt.new(simulation_results_2_9[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx9, nsim = 100)
plotCt.new(simulation_results_2_7[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx7, nsim = 100)
plotCt.new(simulation_results_2_4[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx4, nsim = 100)

# Optima 1,5
plotCt.new(simulation_results_5_9[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx9, nsim = 100)
plotCt.new(simulation_results_5_7[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx7, nsim = 100)
plotCt.new(simulation_results_5_4[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx4, nsim = 100)

# Optima 1,9
plotCt.new(simulation_results_9_4[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx9, nsim = 100)
plotCt.new(simulation_results_9_4[[5]][[2]], phy = anilios_tree, focaltaxa = focal_tx7, nsim = 100)
plotCt.new(simulation_results_9_4[[3]][[2]], phy = anilios_tree, focaltaxa = focal_tx4, nsim = 100)
dev.off()

# PlotCt function
load("data/script_generated_data/simulation_Ct_fromFit_theta19.Rdata")
load("data/script_generated_data/simulation_Ct_fromFit.Rdata")
load("data/script_generated_data/simulation_Ct_fromFit_7taxa.Rdata")
load("data/script_generated_data/simulation_Ct_fromFit_4taxa.Rdata")

pdf(file = "output/plotct_sim_traits_7tx.pdf", height = 11, width = 8)
layout.matrix <- matrix(c(1:15), nrow = 3, ncol = 5)
layout(mat = layout.matrix, heights = c(1,1,1,1,1),widths = c(2,2,2,2,2))

# 9tx
plotCt.phylo(simulation_results_1_9[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx9, nsim = 100)
plotCt.new(simulation_results_1_9[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx9, nsim = 100)
plotCt.new(simulation_results_2_9[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx9, nsim = 100)
plotCt.new(simulation_results_5_9[[5]][[2]], phy = anilios_tree, focaltaxa = focal_tx9, nsim = 100)
plotCt.new(simulation_results_9_9[[6]][[2]], phy = anilios_tree, focaltaxa = focal_tx9, nsim = 100)

# 7tx
plotCt.phylo(simulation_results_1_7[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx7, nsim = 100)
plotCt.new(simulation_results_1_7[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx7, nsim = 100)
plotCt.new(simulation_results_2_7[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx7, nsim = 100)
plotCt.new(simulation_results_5_7[[5]][[2]], phy = anilios_tree, focaltaxa = focal_tx7, nsim = 100)
plotCt.new(simulation_results_9_7[[6]][[2]], phy = anilios_tree, focaltaxa = focal_tx7, nsim = 100)

# 4 tx
plotCt.phylo(simulation_results_1_4[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx, nsim = 100)
plotCt.new(simulation_results_1_7[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx, nsim = 100)
plotCt.new(simulation_results_2_7[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx, nsim = 100)
plotCt.new(simulation_results_5_7[[5]][[2]], phy = anilios_tree, focaltaxa = focal_tx, nsim = 100)
plotCt.new(simulation_results_9_7[[6]][[2]], phy = anilios_tree, focaltaxa = focal_tx, nsim = 100)
dev.off()

library(ggplot2)
library(ggExtra)
library(ggpubr)
# Save the scatter plot in a variable
p1 <- ggplot(sim_df1, aes(x = p_value, y = Ct1)) +  geom_point() + ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m1 <- ggMarginal(p1, type = "histogram", fill = 4, xparams=list(fill = "red"))
p2 <- ggplot(sim_df2, aes(x = p_value, y = Ct1)) +  geom_point() + ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m2 <- ggMarginal(p2, type = "histogram", fill = 4, xparams=list(fill = "red"))
p5 <- ggplot(sim_df5, aes(x = p_value, y = Ct1)) +  geom_point() + ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m5 <- ggMarginal(p5, type = "histogram", fill = 4, xparams=list(fill = "red"))
p9 <- ggplot(sim_df9_9, aes(x = p_value, y = Ct1)) +  geom_point() + ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m9 <- ggMarginal(p9, type = "histogram", fill = 4, xparams=list(fill = "red"))

# 7 focal taxa
# load("data/script_generated_data/simulation_Ct_fromFit_7taxa.Rdata")
# rm(simulation_results_1_7, simulation_results_2_7, simulation_results_5_7)
p1_7 <- ggplot(sim_df1_7, aes(x=p_value, y=Ct1))+geom_point()+ ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m1_7 <- ggMarginal(p1_7, type = "histogram", fill = 4,xparams=list(fill = "red"))
p2_7 <- ggplot(sim_df2_7, aes(x=p_value, y=Ct1))+geom_point()+ ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m2_7 <- ggMarginal(p2_7, type = "histogram", fill = 4,xparams=list(fill = "red"))
p5_7 <- ggplot(sim_df5_7, aes(x=p_value, y=Ct1))+geom_point()+ ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m5_7 <- ggMarginal(p5_7, type = "histogram", fill = 4,xparams=list(fill = "red"))
p9_7 <- ggplot(sim_df9_7, aes(x=p_value, y=Ct1))+geom_point()+ ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m9_7 <- ggMarginal(p9_7, type = "histogram", fill = 4,xparams=list(fill = "red"))

# 4 focal taxa
# load("data/script_generated_data/simulation_Ct_fromFit_4taxa.Rdata")
rm(simulation_results_1_4, simulation_results_2_4, simulation_results_5_4)
p1_4 <- ggplot(sim_df1_4, aes(x = p_value, y = Ct1)) +  geom_point() + ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m1_4 <- ggMarginal(p1_4, type = "histogram", fill = 4,xparams=list(fill = "red"))
p2_4 <- ggplot(sim_df2_4, aes(x = p_value, y = Ct1)) +  geom_point() + ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m2_4 <- ggMarginal(p2_4, type = "histogram", fill = 4,xparams=list(fill = "red"))
p5_4 <- ggplot(sim_df5_4, aes(x = p_value, y = Ct1)) +  geom_point() + ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m5_4 <- ggMarginal(p5_4, type = "histogram", fill = 4,xparams=list(fill = "red"))
p9_4 <- ggplot(sim_df9_4, aes(x = p_value, y = Ct1)) +  geom_point() + ylim(-2,1)+xlim(0,1)+geom_hline(yintercept=0,color="blue")+geom_vline(xintercept=0.05,color="red"); m9_4 <- ggMarginal(p9_4, type = "histogram", fill = 4,xparams=list(fill = "red"))

# Plot
pdf(file = "output/power_test_c_measures.pdf", height = 8, width = 11.33)
ggpubr::ggarrange(m1, m2, m5, m9,
                  m1_7, m2_7, m5_7, m9_7,
                  m1_4, m2_4, m5_4, m9_4,
                  labels = LETTERS[1:12], 
                  nrow = 3, ncol = 4)
dev.off()

# Visualise
plotCt(simulation_results_9_4[[1]][[2]], phy = anilios_tree, focaltaxa = focal_tx, nsim = 100)

# Calculate proportion of significant p-values

which(sim_df1$Ct1>0)
summ_resultr <- function(.x, trt){
  df <- data.frame(trt = trt,
                   median_Ct1 = round(median(.x$Ct1), 2),
                   prop_convg = length(which(.x$Ct1>0)/length(.x$Ct1)*100),
                   percent_sig = length(which(.x$p_value < 0.05))/length(.x$p_value)*100)
  return(df)
}

# 
power_result <- rbind(summ_resultr(sim_df1, "9FT_1_1"),summ_resultr(sim_df2, "9FT_1_2"), summ_resultr(sim_df5, "9FT_1_5"),summ_resultr(sim_df9_9, "9FT_1_9"),
                      summ_resultr(sim_df1_7, "7FT_1_1"), summ_resultr(sim_df2_7, "7FT_1_2"),summ_resultr(sim_df5_7, "7FT_1_5"),summ_resultr(sim_df9_7, "7FT_1_9"),
                      summ_resultr(sim_df1_4, "4FT_1_1"),summ_resultr(sim_df2_4, "4FT_1_2"),summ_resultr(sim_df5_4, "4FT_1_5"),summ_resultr(sim_df9_4, "4FT_1_9")
)

write.csv(power_result, file = "manuscript_word/supplements/sim_power_conv.csv", row.names = F)
