# 01_model_fitting.R
# Sarin Tiatragul 
# March 2023

# Fitting continuos trait evolution models and testing for phylogenetic signal in various traits
# Also testing if variation in morphology among species correlate with their current environments? 

# libraries ---------------------------------------------------------------
# Packages for analyses
library(phytools); library(ape); library(mvMORPH); library(geiger)
library(geomorph); library(phylolm); library(dplyr)
library(MuMIn)
library(R.utils)

'%notin%' <- Negate('%in%')
source('code/utility/func_calc_ses.R')

##
#### MAIN QUESTIONS
## 
# * Which evolutionary model best fit the data?
# * Is there phylogenetic signal?
# * What is the best predictor of body shape variation in blindsnakes? 

# data --------------------------------------------------------------------
# Load full tree, subsetted tree, and subset 
load('data/script_generated_data/subset_mcmctree_shape.rda')
anilios_tree <- ape::drop.tip(sub_phy, tip = c("Ramphotyphlops_multilineatus","Acutotyphlops_subocularis"))
# anilios_tree <- ape::drop.tip(blindsnake_tree, tip = c("Ramphotyphlops_multilineatus","Acutotyphlops_subocularis"))

# Load two D array of y symmetric component of dorsal head shape
load('data/script_generated_data/dorsal_head_shape.rda')

# Load data frame of summarised traits
anilios_data <- read.csv(file = 'data/script_generated_data/anilios_summary_data.csv', row.names = 1)
anilios_tree$tip.label

# Model fitting -----------------------------------------------------------

###
# * Which evolutionary model best fit the phylogeny data?
###

#### Univariate --- maximum likelihood models

fit1 <- phylolm(log(max_svl) ~ 1, data = anilios_data, phy = anilios_tree, model = "BM", measurement_error = TRUE, boot=100)
fit2 <- phylolm(log(max_svl) ~ 1, data = anilios_data, phy = anilios_tree, model = "OUrandomRoot", measurement_error = TRUE, boot=100)
fit3 <- phylolm(log(max_svl) ~ 1, data = anilios_data, phy = anilios_tree, model = "OUfixedRoot", measurement_error = TRUE, boot=100)
fit4 <- phylolm(log(max_svl) ~ 1, data = anilios_data, phy = anilios_tree, model = "EB", measurement_error = TRUE, lower.bound = -10, upper.bound = 10, boot=100)
fit5 <- phylolm(log(max_svl) ~ 1, data = anilios_data, phy = anilios_tree, model = "lambda", boot=100)

AIC(fit1); AIC(fit2); AIC(fit3); AIC(fit4); AIC(fit5) # EB is the best model

#### Using mvMORPH package we fit multivariate data set for body traits and head shape
#### Multivariate body traits include log shape ratios of mbd, mtw, hwe, hda, tl 

Y_body <- as.matrix(anilios_data[, c("sh_mbd", "sh_mtw", "sh_hwe", "sh_hda", "sh_tl")])
body_sh_traits <- list(Y_body=Y_body) # prepared as list for mvMORPH

# Fit the 'BM', 'OU', and 'EB' models to 'Y' using mvMORPH::mvgls
# Using maximum likelihood method
fit_bsh_bm <- mvgls(Y_body~1, data=body_sh_traits, anilios_tree, "BM", method="LL", error = TRUE)
fit_bsh_ou <- mvgls(Y_body~1, data=body_sh_traits, anilios_tree, "OU", method="LL", error = TRUE)
fit_bsh_eb <- mvgls(Y_body~1, data=body_sh_traits, anilios_tree, "EB", method="LL", error = TRUE)
fit_bsh_lambda <- mvgls(Y_body~1, data=body_sh_traits, anilios_tree, "lambda", method="LL")

# Combine model fit in list
bsh_model_fit <- list(fit_bsh_bm, fit_bsh_ou, fit_bsh_eb, fit_bsh_lambda)

# Data table to compare
fit_bod_evol_mods_aic <- lapply(bsh_model_fit, AIC)
names(fit_bod_evol_mods_aic) <- c('BM', 'OU', 'EB', 'lambda')

summary(fit_bsh_bm)
summary(fit_bsh_ou)
summary(fit_bsh_eb)
summary(fit_bsh_lambda)

#### Multivariate head traits from symmetric component
# Array of mean dorsal landmarks from gpa and symmetric analyses from dorsal_headshape.rda
rownames(sp_means_dorsal) <- gsub(x = rownames(sp_means_dorsal), pattern = "^", replacement = "Anilios_")
means_2d_dorsal <- sp_means_dorsal[rownames(sp_means_dorsal) %in% anilios_tree$tip.label,]
means_2d_dorsal <- means_2d_dorsal[match(anilios_tree$tip.label, rownames(means_2d_dorsal)), ]

# Turn head landmark traits 
Y_head <- as.matrix(means_2d_dorsal[-c(1,2)])
head_traits <- list(Y_head=Y_head, c_size = means_2d_dorsal$c_size, body_shape = anilios_data$sh_mbd, svl = anilios_data$log_svl)

## Test for shape-size allometry across species 
# We include c_size here to account for allometric effects. The residuals is size-free shape
# Using ridge quadratic null penalised likelihood LOOCV (recommended for low and high dimensional multivariate traits from in Clavel et al. 2019)
fit_head_bm <- mvgls(Y_head~c_size, data=head_traits, anilios_tree, "BM", penalty="RidgeArch", target = "unitVariance", method="PL-LOOCV", error = TRUE)
fit_head_ou <- mvgls(Y_head~c_size, data=head_traits, anilios_tree, "OU", penalty="RidgeArch", target = "unitVariance", method="PL-LOOCV", error = TRUE)
fit_head_eb <- mvgls(Y_head~c_size, data=head_traits, anilios_tree, "EB", penalty="RidgeArch", target = "unitVariance", method="PL-LOOCV", error = TRUE)
fit_head_lambda <- mvgls(Y_head~c_size, data=head_traits, anilios_tree, "lambda", penalty="RidgeArch", target = "unitVariance", method="PL-LOOCV")

# Combine model fit in list
head_model_fit <- list(fit_head_bm, fit_head_ou, fit_head_eb, fit_head_lambda)

# Data table to compare
fit_head_GIC <- data.frame(models = c("BM", "OU", "EB", "Lambda"), GIC = NA)
for (i in 1:length(head_model_fit)) {
  fit_head_GIC$GIC[i] <- GIC(head_model_fit[[i]])$GIC
}

# GIC Table for body shape
fit_head_GIC[order(fit_head_GIC$GIC),]

# BM and lambda has the lowest GIC score for head traits

summary(fit_head_bm)
summary(fit_head_ou)
summary(fit_head_eb)
summary(fit_head_lambda)
# manova.gls(fit_head_lambda, nperm=999, test="Pillai", type = "II")

fit_head_lambda$coefficients
fit_head_lambda$sigma

# Relationship between headshape and body size ----------------------------
fit_head_svl <- mvgls(Y_head~svl, data=head_traits, anilios_tree, "lambda", penalty="RidgeArch", target = "unitVariance", method="PL-LOOCV")
fit_head_mbd <- mvgls(Y_head~body_shape, data=head_traits, anilios_tree, "lambda", penalty="RidgeArch", target = "unitVariance", method="PL-LOOCV")

summary(fit_head_svl)
summary(fit_head_mbd)

aov_fit_head_svl <- manova.gls(fit_head_svl, nperm=999, test="Pillai", type = "II")
aov_fit_head_mbd <- manova.gls(fit_head_mbd, nperm=999, test="Pillai", type = "II")

## check using procD.pgls to see if results similar

body_width <- head_traits$body_shape
names(body_width) <- rownames(head_traits$Y_head)
fit_head_mbd.pgls <- geomorph::procD.pgls(head_traits$Y_head ~ body_width, phy = anilios_tree)
summary(fit_head_mbd.pgls)

### No correlation between body body width and head shape

# Use 2-block partial least squares test under Brownian motion with geomorph phylo.integration()
two.b.pls <- phylo.integration(A = head_traits$Y_head, 
                  A2 = body_sh_traits$Y_body, 
                  phy = anilios_tree, iter = 999)

two.b.pls
labels_anilios <- gsub(pattern = "Anilios_", replacement = "", x = rownames(head_traits$Y_head))
plot(two.b.pls, label = labels_anilios)

## Weak correlation between between shape ratios and snout shape

###
### * Is there phylogenetic signal
###

summary(fit_bsh_lambda)
summary(fit_head_lambda)

# Phylogenetic signal in traits, closely related species tend to occupy similar morphospace
## Interpretation: Closely related species tend to have similar traits in both body shape and head

# Estimate PCA for traits -------------------------------------------------
## Following approach in Garcia-Porta et al 2022. We calculate PCA based on the variance covariance matrix 
## form the best evolutionary models

# dev.off()
body_pca_results <- mvgls.pca(fit_bsh_lambda, plot=T) 
head_pca_results <- mvgls.pca(fit_head_lambda, plot=T)

rownames(body_pca_results$vectors) <- colnames(fit_bsh_lambda$coefficients)
rownames(body_pca_results$vectors)

body_pca_results$contrib <- round(body_pca_results$values / sum(body_pca_results$values) * 100, 1)
head_pca_results$contrib <- round(head_pca_results$values / sum(head_pca_results$values) * 100, 1)


## Visualise in phoylomorphospace with PCA (not phy PCA)
load('data/script_generated_data/phylomorph_tip_colours.rda')

# xlim = range(body_pca_df$PC1), ylim = range(body_pca_df$PC2),
body_pca_mat <- body_pca_results$scores[,1:2]
# body_pca_mat[,1] <- body_pca_mat[,1]*-1

## Add conventional pca data
load('data/script_generated_data/conventional_pca.rda')
load('data/script_generated_data/linear_measurements/PC_log_ratio_pc_obj.rda')
body_pc_variance <- get_eig(b_mean_PC)
body_pc_variance$variance.percent[1]

rownames(blindsnake_data) <- blindsnake_data$species

# PRINT TO PDF
# pdf(file = 'output/combo_pca_phypca.pdf', width = 11, height = 8.5)
par(mfrow=c(2,2))

body_intree <- blindsnake_data[, c("fil_PC1", "fil_PC2")]
body_intree$fil_PC1 <- body_intree$fil_PC1* -1

# CONVENTIONAL PCA body
phylomorphospace(sub_phy, body_intree[, c("fil_PC1", "fil_PC2")], 
                 label = "horizontal", fsize = 1, ftype = 'i', bty='l',
                 xlab = paste("Body PC1 (", round(body_pc_variance$variance.percent[1], 2), "%)", sep=""),
                 ylab = paste("Body PC2 (", round(body_pc_variance$variance.percent[2], 2), "%)", sep=""),
                 pch = 20, cex = 0, node.size = c(0, 1.4),
                 xlim = c(-6,6), ylim = c(-2, 3),
                 # xlim = c(-6,6), 
                 # ylim = c(-0.4, 0.6),
                 control = list(col.node = tip_cols_cont))
points(y = not_in_phy_df$fil_PC2, x = not_in_phy_df$fil_PC1 * -1, pch = 19, 
       col = not_in_phy_tip_colours, cex = 1.4)
text(y = not_in_phy_df$fil_PC2 - 0.03, x = not_in_phy_df$fil_PC1 * -1,
     labels = rownames(not_in_phy_df), cex = 0.8)

# Plot head PC1 and PC2
head_intree <- blindsnake_data[, c("mean_PC1", "mean_PC2")]
head_intree$mean_PC1 <- head_intree$mean_PC1 *-1

phylomorphospace(blindsnake_tree, head_intree[, c("mean_PC1", "mean_PC2")], 
                 label = "horizontal", fsize = 1, ftype = 'i', bty='l',
                 xlab = paste("Head PC1 (", round(head_pc_variance$variance.percent[1], 2), "%)", sep=""),
                 ylab = paste("Head PC2 (", round(head_pc_variance$variance.percent[2], 2), "%)", sep=""), 
                 pch = 20, cex = 0, node.size = c(0, 1.4),
                 # xlim = c(-0.2019596, 0.2), ylim = range(head_pca_df$head_PC2),
                 # xlim = c(-0.15, 0.15), ylim = c(-0.10, 0.05),
                 # xlim = range(-head_pca_df$head_PC1), ylim = range(head_pca_df$head_PC2),
                 xlim = c(-0.20, 0.20), ylim = c(-.1, 0.1), 
                 control = list(col.node = tip_cols_cont))
points(y = not_in_phy_df$mean_PC2, x = not_in_phy_df$mean_PC1*-1, pch = 19, 
       col = not_in_phy_tip_colours, cex = 1.4)
text(y = not_in_phy_df$mean_PC2 - 0.005, x = not_in_phy_df$mean_PC1*-1,
     labels = rownames(not_in_phy_df), cex = 0.8)


## Phy-PCAs
# Switch the PC vectors around (so that PC has all positive loading)
# Phylogenetic PCA
phylomorphospace(anilios_tree, body_pca_mat, 
                 label = "horizontal", fsize = 1, ftype = 'i', bty='l',
                 xlab = paste("Body PC1 (", body_pca_results$contrib[1], "%)", sep=""),
                 ylab = paste("Body PC2 (", body_pca_results$contrib[2], "%)", sep=""), 
                 xlim = c(-2,2), ylim = c(-0.4, 0.6),
                 # xlim = range(-body_pca_mat[,1]), ylim = range(body_pca_mat[,2]),
                 pch = 20, cex = 0, node.size = c(0, 1.4),
                 control = list(col.node = soil_tip_cols))
# Include legend
# legend_image <- as.raster(matrix(rev(my_palette)), ncol=1)
# text(x=-1.5, y = seq(0, 0.4,l=5), labels = round(seq(min(anilios_data$mean_bulk, na.rm = T), max(anilios_data$mean_bulk, na.rm = T),l=5), 2))
# rasterImage(legend_image, -1.5, 0, -1.25, 0.4)

head_pca_mat <- head_pca_results$scores[,1:2]

dev.off()
par(mfrow=c(2,2))

# Plot head shape
phylomorphospace(anilios_tree, head_pca_mat, 
                 label = "horizontal", fsize = 1, ftype = 'i', bty='l',
                 xlab = paste("Head PC1 (", head_pca_results$contrib[1], "%)", sep=""),
                 ylab = paste("Head PC2 (", head_pca_results$contrib[2], "%)", sep=""), 
                 xlim = c(-0.2,0.2), ylim = c(-0.11, 0.06),
                 # xlim = range(head_pca_mat[,1]), ylim = range(head_pca_mat[,2])
                 xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2),
                 pch = 20, cex = 0, node.size = c(0, 1.4),
                 control = list(col.node = soil_tip_cols))


dev.off()
anilios_data[anilios_tree$tip.label,]




#####
##### * What is the best predictor of body shape variation in blind snakes? 
#####
## With mvMORPH we can fit multivariate traits instead of PC scores. 

# Env predictors including soil compactness, percentage sand, mean temp, and mean prec
soil.bulk <- anilios_data$mean_bulk; names(soil.bulk) <- rownames(anilios_data)
temp_m <- anilios_data$temp_mean; names(temp_m) <- rownames(anilios_data)
arid_m <- anilios_data$ARID; names(arid_m) <- rownames(anilios_data)
arid_c <- anilios_data$arid_cat; names(arid_c) <- rownames(anilios_data)

## Fit univariate maximum svl
fit_max_svl_1 <- phylolm(log(max_svl) ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
fit_max_svl_1 <- phylolm(log(max_svl) ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "EB", lower.bound = -10, upper.bound = 10, measurement_error = TRUE, boot = 100)

summary(fit_max_svl_1)
summary(fit_max_svl_1eb)

cor(anilios_data$mean_bulk, anilios_data$ARID, method = "pearson")

fit_max_svl_2 <- phylolm(log(max_svl) ~ temp_mean, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 100)
fit_max_svl_3 <- phylolm(log(max_svl) ~ mean_bulk, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 100)
fit_max_svl_4 <- phylolm(log(max_svl) ~ ARID, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 100)
fit_max_svl_5 <- phylolm(log(max_svl) ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 100)

sum_fit_svl_2 <- summary(fit_max_svl_2); sum_fit_svl_2 <- round(sum_fit_svl_2$coefficients, 2)
sum_fit_svl_3 <- summary(fit_max_svl_3); sum_fit_svl_3 <- round(sum_fit_svl_3$coefficients, 2)
sum_fit_svl_4 <- summary(fit_max_svl_4); sum_fit_svl_4 <- round(sum_fit_svl_4$coefficients, 2)
sum_fit_svl_5 <- summary(fit_max_svl_5); sum_fit_svl_5 <- round(sum_fit_svl_5$coefficients, 2)

# Step fitting
anilios_data$log_max_svl <- log(anilios_data$max_svl)
step_fit <- phylostep(log_max_svl ~ mean_bulk + temp_mean + ARID, direction = "backward", data = anilios_data, phy = anilios_tree, model = "lambda", k=2)
step_fit

# pdf("output/evolution_4_svl_temp.pdf", height = 7.5, width = 13.33)
plot(log(max_svl) ~ temp_mean, data = anilios_data, bty ="n", pch=20, cex = 1.4,
     xlab = "Temperature", ylab = "Maximum length")
abline(fit_max_svl_2, lty = 5)
text(log(max_svl) ~ temp_mean, data = anilios_data, labels = anilios_data$species)


###
### MULTIVARIATE MODEL FIT
###

# Add the predictors to the previous 'data' list
body_sh_traits_env <- list(Y_body=Y_body, soil = soil.bulk, temp = temp_m, arid = arid_m,
                           PCs = body_pca_results$scores)

# Difference body shape between arid and non arid species? 
## We use Pagel's lambda because this seen as a mixed model [@clavelReliable2020]

## Fit model individually

# Body ~ soil compactness? 
fit_soil <- mvgls(Y_body ~ soil, data=body_sh_traits_env, anilios_tree, model="lambda", error =TRUE)
aov_soil_bd <- manova.gls(fit_soil, nperm=999, test="Pillai", type = "II")

summary(fit_soil)
fit_soil$sigm

# Body ~ temperature
fit_temp <- mvgls(Y_body ~ temp, data=body_sh_traits_env, anilios_tree, model="lambda")
aov_temp_bd <- manova.gls(fit_temp, nperm=999, test="Pillai", type = "II")

# Body ~ aridity
fit_arid <- mvgls(Y_body ~ arid, data=body_sh_traits_env, anilios_tree, model="lambda")
aov_arid_bd <- manova.gls(fit_arid, nperm=999, test="Pillai", type = "II")

## Fitted together
# Body ~ multiple variables
fit_all <- mvgls(Y_body ~ soil + temp + arid, data=body_sh_traits_env, anilios_tree, model="lambda")
# fit_allou <- mvgls(Y_body ~ soil + temp + biome, data=body_sh_traits_env, anilios_tree, model="OU", error = TRUE)
summary(fit_all)
aov_mult_bd <- manova.gls(fit_all, nperm=999, test="Pillai", type = "II")
ses.calc(aov_mult_bd)
effectsize(aov_mult_bd)


# Using parametric test ("LL")
# Body ~ soil compactness? 

# Full model
fit_mult_ll <- mvgls(Y_body ~ soil + temp + arid, data=body_sh_traits_env, anilios_tree, model="lambda", method = "LL")
aov_mult_bd_ll <- manova.gls(fit_mult_ll, type = "III")

# Partial models
fit_soar_ll <- mvgls(Y_body ~ soil + arid, data=body_sh_traits_env, anilios_tree, model="lambda", method = "LL")
aov_soar_bd_ll <- manova.gls(fit_soar_ll, type = "III")

fit_sotm_ll <- mvgls(Y_body ~ soil + temp, data=body_sh_traits_env, anilios_tree, model="lambda", method = "LL")
aov_sotm_bd_ll <- manova.gls(fit_sotm_ll, type = "III")

fit_artm_ll <- mvgls(Y_body ~ arid + temp, data=body_sh_traits_env, anilios_tree, model="lambda", method = "LL")
aov_artm_bd_ll <- manova.gls(fit_artm_ll, type = "III")

# Single models
fit_soil_ll <- mvgls(Y_body ~ soil, data=body_sh_traits_env, anilios_tree, model="lambda", method = "LL")
aov_soil_bd_ll <- manova.gls(fit_soil_ll, type = "III")

# Body ~ temperature
fit_temp_ll <- mvgls(Y_body ~ temp, data=body_sh_traits_env, anilios_tree, model="lambda", method = "LL")
aov_temp_bd_ll <- manova.gls(fit_temp_ll, type = "III")

# Body ~ aridity
fit_arid_ll <- mvgls(Y_body ~ arid, data=body_sh_traits_env, anilios_tree, model="lambda", method = "LL")
aov_arid_bd_ll <- manova.gls(fit_arid_ll, type = "III")

fit_bod_mods <- list(soil = fit_soil_ll, 
                     temp = fit_temp_ll, 
                     arid = fit_arid_ll, 
                     soar = fit_soar_ll, sotm = fit_sotm_ll, artm = fit_artm_ll, mult = fit_mult_ll)

fit_bod_mods_aic <- lapply(fit_bod_mods, AIC)
lapply(fit_bod_mods, "[[", "formula") %>% paste()

aov_mult_bd_ll



###
### Fitting regression models with head traits
###

# In the first approach we regress head shape against env variables + log centroid size
head_traits_env <- list(Y_head = as.matrix(means_2d_dorsal[-c(1,2)]),
                        res_head = fit_head_lambda$residuals,
                        soil = soil.bulk, temp = temp_m, 
                        arid = arid_m, biome = arid_c,
                        PCs = head_pca_results$scores[,1:4],
                        c_size = log(means_2d_dorsal$c_size))

## START with full model

# Head ~ multiple variables
fit_head_mult <- mvgls(Y_head ~ soil + temp + arid_m + c_size + c_size*soil*temp*arid_m, data=head_traits_env, anilios_tree, model="lambda")
# aov_mult_hd <- manova.gls(fit_head_mult, nperm=999, test="Pillai", type = "III")
fit_head_mult2 <- mvgls(Y_head ~ soil + temp + arid_m + c_size, data=head_traits_env, anilios_tree, model="lambda")
# aov_mult_hd2 <- manova.gls(fit_head_mult2, nperm=999, test="Pillai", type = "III")
fit_head_mult3 <- mvgls(Y_head ~ soil + temp + arid_m + c_size + soil*temp*arid_m, data=head_traits_env, anilios_tree, model="lambda")
# aov_mult_hd3 <- manova.gls(fit_head_mult3, nperm=999, test="Pillai", type = "II")


fit_mod_head_mult <- lapply(list(fit_head_mult, fit_head_mult2, fit_head_mult3), EIC)
lapply(fit_mod_head_mult, "[[", "formula") %>% paste()

# EIC scores fit_head_mult2 model without interaction is preferred 

# Reduced terms
fit_head_soar <- mvgls(Y_head ~ soil + arid_m + c_size, data=head_traits_env, anilios_tree, model="lambda")
aov_head_soar <- manova.gls(fit_head_soar, nperm=999, test="Pillai", type = "III")
fit_head_sotm <- mvgls(Y_head ~ soil + temp + c_size, data=head_traits_env, anilios_tree, model="lambda")
aov_head_sotm <- manova.gls(fit_head_sotm, nperm=999, test="Pillai", type = "III")
fit_head_artm <- mvgls(Y_head ~ arid_m + temp + c_size, data=head_traits_env, anilios_tree, model="lambda")
aov_head_artm <- manova.gls(fit_head_artm, nperm=999, test="Pillai", type = "III")

# Individual predictors
# Head ~ soil compactness? 
fit_head_soil <- mvgls(Y_head ~ soil + c_size, data=head_traits_env, anilios_tree, model="lambda")
aov_soil_hd <- manova.gls(fit_head_soil, nperm=999, test="Pillai", type = "III")

# Head ~ temperature
fit_head_temp <- mvgls(Y_head ~ temp + c_size, data=head_traits_env, anilios_tree, model="lambda")
aov_temp_hd <- manova.gls(fit_head_temp, nperm=999, test="Pillai", type = "III")

# Head ~ aridity
fit_head_arid <- mvgls(Y_head ~ arid_m + c_size, data=head_traits_env, anilios_tree, model="lambda")
aov_arid_hd <- manova.gls(fit_head_arid, nperm=999, test="Pillai", type = "III")

# Print results from MANOVA 
aov_soil_hd; aov_arid_hd; aov_temp_hd; aov_mult_hd2

## Calculate EIC
fit_head_mods <- list(fit_head_soil, fit_head_temp, fit_head_arid, fit_head_soar, 
                      fit_head_sotm, fit_head_artm, fit_head_mult2)

fit_head_mods_eic <- lapply(fit_head_mods, EIC)
names(fit_head_mods_eic) <- lapply(fit_head_mods, "[[", "formula") %>% paste()

# Add the predictors to the previous 'data' list
head_traits_env <- list(Y_head = as.matrix(means_2d_dorsal[-c(1,2)]),
                        res_head = fit_head_lambda$residuals,
                        soil = soil.bulk, temp = temp_m, 
                        arid = arid_m, biome = arid_c,
                        PCs = head_pca_results$scores[,1:4],
                        c_size = log(means_2d_dorsal$c_size))

# In another approach
# We accounted for static allometry in the head shape data set by regressing head shape against log centroid size 
# and used the regression residuals to represents size-free variation of head shape in the PCA

 # Head ~ soil compactness? 
fit_head_soil1 <- mvgls(res_head ~ soil, data=head_traits_env, anilios_tree, model="lambda")
summary(fit_head_soil1)
aov_soil_hd1 <- manova.gls(fit_head_soil, nperm=999, test="Pillai", type = "II")

 # Head ~ temperature
fit_head_temp1 <- mvgls(res_head ~ temp, data=head_traits_env, anilios_tree, model="lambda")
summary(fit_head_temp1)
aov_temp_hd1 <- manova.gls(fit_head_temp, nperm=999, test="Pillai", type = "II")

 # Head ~ aridity
fit_head_arid1 <- mvgls(res_head ~ arid_m, data=head_traits_env, anilios_tree, model="lambda")
summary(fit_head_arid1)
aov_arid_hd1 <- manova.gls(fit_head_arid, nperm=999, test="Pillai", type = "II")

 # Head ~ multiple variables
fit_head_all1 <- mvgls(res_head ~ soil + temp + arid_m, data=head_traits_env, anilios_tree, model="lambda")
summary(fit_head_all1)
aov_mult_hd1 <- manova.gls(fit_head_all1, nperm=999, test="Pillai", type = "II")

# Print results from MANOVA and MANCOVA
aov_soil_hd1; aov_arid_hd1; aov_temp_hd1; aov_mult_hd1


# DATA PRESENTATION -------------------------------------------------------

###
### Plotting response together
###

# pdf(file = "output/mvmorph_traits_by_env.pdf", width = 8.5, height = 11)
par(mfrow=c(3,3))
plot(log(max_svl) ~ mean_bulk, data = anilios_data, bty ="n", pch=20, cex = 1.4, col = "steelblue",
     xlab = expression(Soil ~ density ~ (g/cm^3)) , ylab = "Maximum length (mm)")
# text(y = log(anilios_data$max_svl), x = anilios_data$mean_bulk, labels = anilios_data$species)
text(x = 1.2, y = 6.3, labels = paste("beta = ", sum_fit_svl_3[2], "+/-", sum_fit_svl_3[4]), adj = 0, cex = 1.4)
text(x = 1.2, y = 6.25, labels = paste("P = ", sum_fit_svl_3[12]), adj = 0, cex = 1.4)
abline(fit_max_svl_3, lty = 5)

plot(log(max_svl) ~ ARID, data = anilios_data, bty ="n", pch=20, cex = 1.4, col = "steelblue",
     xlab = "Aridity index", ylab = "Maximum length (mm)")
text(x = 0.2, y = 6.3, labels = paste("beta = ", sum_fit_svl_4[2], "+/-", sum_fit_svl_4[4]), adj = 0, cex = 1.4)
text(x = 0.2, y = 6.25, labels = paste("P = ", sum_fit_svl_4[12]), adj = 0, cex = 1.4)
abline(fit_max_svl_4, lty = 5)

plot(log(max_svl) ~ temp_mean, data = anilios_data, bty ="n", pch=20, cex = 1.4, col = "steelblue",
     xlab = "Temperature (°C)", ylab = "Maximum length (mm)");
text(x = 17, y = 6.3, labels = paste("beta = ", sum_fit_svl_2[2], "+/-", sum_fit_svl_2[4]), adj = 0, cex = 1.4)
text(x = 17, y = 6.25, labels = paste("P = ", sum_fit_svl_2[12]), adj = 0, cex = 1.4)
abline(fit_max_svl_2, lty = 5)


### MULTIVARIATE Body shape 1 

plot(soil.bulk, body_sh_traits_env$Y_body[,1], bty ="n", pch=20, cex = 1.4, col = "steelblue", xlab = expression(Soil ~ density ~ (g/cm^3)), ylab = "Log(midbody width / svl)")
text(soil.bulk, body_sh_traits_env$Y_body[,1], labels = rownames(body_sh_traits_env$Y_body))
text(x = 1.2, y = -3.5, labels = paste("beta = ", round(fit_soil$coefficients[2], 2)), adj = 0, cex = 1.4)
text(x = 1.2, y = -3.8, labels = paste("P = ", aov_soil_bd$pvalue), adj = 0, cex = 1.4)
abline(fit_soil, lty = 5)

plot(arid_m, body_sh_traits_env$Y_body[,1], bty ="n", pch=20, cex = 1.4, col = "steelblue", xlab = "Aridity index", ylab = "Log(midbody width / svl)")
text(arid_m, body_sh_traits_env$Y_body[,1], labels = rownames(body_sh_traits_env$Y_body))
text(x = 0.2, y = -3.5, labels = paste("beta = ", round(fit_arid$coefficients[2], 2)), adj = 0, cex = 1.4)
text(x = 0.2, y = -3.8, labels = paste("P = ", aov_arid_bd$pvalue), adj = 0, cex = 1.4)
abline(fit_arid, lty = 5)

plot(temp_m, body_sh_traits_env$Y_body[,1], bty ="n", pch=20, cex = 1.4, col = "steelblue", xlab = "Temperature (°C)", ylab = "Log(midbody width / svl)")
text(x = 16, y = -3.5, labels = paste("beta = ", round(fit_temp$coefficients[2], 2)), adj = 0, cex = 1.4)
text(x = 16, y = -3.8, labels = paste("P = ", aov_temp_bd$pvalue), adj = 0, cex = 1.4)
abline(fit_temp, lty = 5)

### MULTIVARIATE HEAD SHAPE

mv_pc_head_soil <- mvgls.pca(fit_head_soil, plot=F)
mv_pc_head_arid <- mvgls.pca(fit_head_arid, plot=F)
mv_pc_head_temp <- mvgls.pca(fit_head_temp, plot=F)

plot(soil.bulk, mv_pc_head_soil$scores[,1], bty ="n", pch=20, cex = 1.4, ylim = c(-0.15, 0.16), col = "steelblue", xlab = expression(Soil ~ density ~ (g/cm^3)), ylab = "Head shape PC1")
# plot(fit_head_soil, bty ="n", pch=20, cex = 1.4, col = "steelblue")

plot(arid_m, mv_pc_head_arid$scores[,1], bty ="n", pch=20, cex = 1.4, ylim = c(-0.15, 0.16), col = "steelblue", xlab = "Aridity index", ylab = "Head shape PC1")
# plot(fit_head_arid, bty ="n", pch=20, cex = 1.4, col = "steelblue")

plot(temp_m, mv_pc_head_temp$scores[,1], bty ="n", pch=20, cex = 1.4, ylim = c(-0.15, 0.16), col = "steelblue", xlab = "Temperature (°C)", ylab = "Head shape PC1")
# plot(fit_head_temp, bty ="n", pch=20, cex = 1.4, col = "steelblue")

dev.off()

# PLOT SUPPLEMENTARY FIGURE ALL BODY TRAITS AGAINST PREDICTORS
par(mfrow = c(3,5))

# install.packages("purrr")
# install.packages("magick")
# library(purrr); library(magick);

plot(soil.bulk, body_sh_traits_env$Y_body[,1], bty ="n", pch=20, cex = 1.4, xlab = expression(Soil ~ density ~ (g/cm^3)), ylab = "Log(midbody width / svl)")
plot(soil.bulk, body_sh_traits_env$Y_body[,2], bty ="n", pch=20, cex = 1.4, xlab = expression(Soil ~ density ~ (g/cm^3)), ylab = "Log(midtail width / svl)")
plot(soil.bulk, body_sh_traits_env$Y_body[,3], bty ="n", pch=20, cex = 1.4, xlab = expression(Soil ~ density ~ (g/cm^3)), ylab = "Log(head width / svl)")
plot(soil.bulk, body_sh_traits_env$Y_body[,4], bty ="n", pch=20, cex = 1.4, xlab = expression(Soil ~ density ~ (g/cm^3)), ylab = "Log(head depth/ svl)")
plot(soil.bulk, body_sh_traits_env$Y_body[,5], bty ="n", pch=20, cex = 1.4, xlab = expression(Soil ~ density ~ (g/cm^3)), ylab = "Log(tail length / svl)")

plot(arid_m, body_sh_traits_env$Y_body[,1], bty ="n", pch=20, cex = 1.4, xlab = "Aridity index", ylab = "Log(midbody width / svl)")
plot(arid_m, body_sh_traits_env$Y_body[,2], bty ="n", pch=20, cex = 1.4, xlab = "Aridity index", ylab = "Log(midtail width / svl)")
plot(arid_m, body_sh_traits_env$Y_body[,3], bty ="n", pch=20, cex = 1.4, xlab = "Aridity index", ylab = "Log(head width / svl)")
plot(arid_m, body_sh_traits_env$Y_body[,4], bty ="n", pch=20, cex = 1.4, xlab = "Aridity index", ylab = "Log(head depth/ svl)")
plot(arid_m, body_sh_traits_env$Y_body[,5], bty ="n", pch=20, cex = 1.4, xlab = "Aridity index", ylab = "Log(tail length / svl)")

plot(temp_m, body_sh_traits_env$Y_body[,1], bty ="n", pch=20, cex = 1.4, xlab = "Temperature (°C)", ylab = "Log(midbody width / svl)")
plot(temp_m, body_sh_traits_env$Y_body[,2], bty ="n", pch=20, cex = 1.4, xlab = "Temperature (°C)", ylab = "Log(midtail width / svl)")
plot(temp_m, body_sh_traits_env$Y_body[,3], bty ="n", pch=20, cex = 1.4, xlab = "Temperature (°C)", ylab = "Log(head width / svl)")
plot(temp_m, body_sh_traits_env$Y_body[,4], bty ="n", pch=20, cex = 1.4, xlab = "Temperature (°C)", ylab = "Log(head depth/ svl)")
plot(temp_m, body_sh_traits_env$Y_body[,5], bty ="n", pch=20, cex = 1.4, xlab = "Temperature (°C)", ylab = "Log(tail length / svl)")




# VISUAULISATION ----------------------------------------------------------

## Another version of PCA plot with aridity species
# PRINT TO PDF
load('data/script_generated_data/conventional_pca.rda')
arid_breakpoinst <- c(0, 0.03, 0.2, 0.5, 0.65, Inf)
arid_cat_labels <- c("Hyper Arid", "Arid", "Semi-Arid", "Dry sub-humid", "Humid")

blindsnake_data$arid_category <- cut(blindsnake_data$ARID, breaks = arid_breakpoinst, labels = arid_cat_labels, include.lowest = TRUE)
blindsnake_data_full$arid_category <- cut(blindsnake_data_full$ARID, breaks = arid_breakpoinst, labels = arid_cat_labels, include.lowest = TRUE)

stateblindsnake <- as.character(blindsnake_data$arid_category)
names(stateblindsnake) <- blindsnake_data$species
stateblindsnake[which(is.na(stateblindsnake))] <- "Humid"

tips_col_arid <- tip_cols_cont
col_pal <- scico(n = length(unique(stateblindsnake)), palette = "vik")

stateblindsnake <- gsub("Dry sub-humid", replacement = "#001260", stateblindsnake)
stateblindsnake <- gsub("Humid", replacement = "#71A7C4", stateblindsnake)
stateblindsnake <- gsub("Semi-Arid", replacement = "#D29773", stateblindsnake)
stateblindsnake <- gsub("Arid", replacement = "#590007", stateblindsnake)
tips_col_arid[which(sub_phy$tip.label %in% names(stateblindsnake))] <- stateblindsnake
tips_col_arid

tip_col_arid_phy <- soil_tip_cols
stateblindsnake_phy <- as.character(anilios_data$arid_cat)
names(stateblindsnake_phy) <- anilios_data$species
stateblindsnake_phy <- gsub("Dry sub-humid", replacement = "#001260", stateblindsnake_phy)
stateblindsnake_phy <- gsub("Humid", replacement = "#71A7C4", stateblindsnake_phy)
stateblindsnake_phy <- gsub("Semi-Arid", replacement = "#D29773", stateblindsnake_phy)
stateblindsnake_phy <- gsub("Arid", replacement = "#590007", stateblindsnake_phy)
tip_col_arid_phy[which(anilios_tree$tip.label %in% names(stateblindsnake_phy))] <- stateblindsnake_phy

tip_col_arid_phy

notinarid <- blindsnake_data_full[names(not_in_phy_tip_colours),]["arid_category"]
notinarid1 <- as.character(notinarid$arid_category)
statenotin <- notinarid1
names(statenotin) <- rownames(notinarid)
statenotin[which(is.na(statenotin))] <- "Humid"

statenotin <- gsub("Dry sub-humid", replacement = "#001260", statenotin)
statenotin <- gsub("Humid", replacement = "#71A7C4", statenotin)
statenotin <- gsub("Semi-Arid", replacement = "#D29773", statenotin)
statenotin <- gsub("Arid", replacement = "#590007", statenotin)

statenotin

dev.off()
# pdf(file = 'output/combo_pca_phypca_aridity.pdf', width = 11, height = 8.5)
par(mfrow = c(2,2))
# CONVENTIONAL PCA body
# adjust PC so that it mirrors that from phy-PCA
conv_pca_df <- blindsnake_data[, c("fil_PC1", "fil_PC2")]
conv_pca_df$fil_PC1 <- conv_pca_df$fil_PC1 *-1
# conv_pca_df$fil_PC2 <- conv_pca_df$fil_PC2 *-1

phylomorphospace(sub_phy, conv_pca_df, 
                 label = "horizontal", fsize = 1, ftype = 'i', bty='l',
                 xlab = paste("Body PC1 (", round(body_pc_variance$variance.percent[1], 2), "%)", sep=""),
                 ylab = paste("Body PC2 (", round(body_pc_variance$variance.percent[2], 2), "%)", sep=""),
                 pch = 20, cex = 0, node.size = c(0, 1.4),
                 xlim = c(-6,6), ylim = c(-2, 3),
                 # xlim = c(-6,6), 
                 # ylim = c(-0.4, 0.6),
                 control = list(col.node = tips_col_arid))
points(x = -1 * not_in_phy_df$fil_PC1, y = not_in_phy_df$fil_PC2, pch = 19, 
       col = statenotin, cex = 1.4)
text(x = -1 * not_in_phy_df$fil_PC1 - 0.03, y = not_in_phy_df$fil_PC2,
     labels = rownames(not_in_phy_df), cex = 0.8)


# Plot head PC1 and PC2
conv_pca_head_df <- blindsnake_data[, c("mean_PC1", "mean_PC2")]
conv_pca_head_df$mean_PC2 <- conv_pca_head_df$mean_PC2*-1

phylomorphospace(blindsnake_tree, conv_pca_head_df, 
                 label = "horizontal", fsize = 1, ftype = 'i', bty='l',
                 xlab = paste("Head PC1 (", round(head_pc_variance$variance.percent[1], 2), "%)", sep=""),
                 ylab = paste("Head PC2 (", round(head_pc_variance$variance.percent[2], 2), "%)", sep=""), 
                 pch = 20, cex = 0, node.size = c(0, 1.4),
                 # xlim = c(-0.2019596, 0.2), ylim = range(head_pca_df$head_PC2),
                 # xlim = c(-0.15, 0.15), ylim = c(-0.10, 0.05),
                 # xlim = range(-head_pca_df$head_PC1), ylim = range(head_pca_df$head_PC2),
                 xlim = c(-0.2,0.2), ylim = c(-.1, 0.1),
                 control = list(col.node = tips_col_arid))
points(y = not_in_phy_df$mean_PC2*-1, x = not_in_phy_df$mean_PC1, pch = 19, 
       col = statenotin, cex = 1.4)
text(y = not_in_phy_df$mean_PC2 *-1 - 0.005, x = not_in_phy_df$mean_PC1 ,
     labels = rownames(not_in_phy_df), cex = 0.8)


## Phy-PCAs
# Switch the PC vectors around (so that PC has all positive loading)
# Phylogenetic PCA
phylomorphospace(anilios_tree, body_pca_mat, 
                 label = "horizontal", fsize = 1, ftype = 'i', bty='l',
                 xlab = paste("Body PC1 (", body_pca_results$contrib[1], "%)", sep=""),
                 ylab = paste("Body PC2 (", body_pca_results$contrib[2], "%)", sep=""), 
                 xlim = c(-2,2), ylim = c(-0.4, 0.6),
                 # xlim = range(-body_pca_mat[,1]), ylim = range(body_pca_mat[,2]),
                 pch = 20, cex = 0, node.size = c(1, 1.4),
                 control = list(col.node = tip_col_arid_phy))
# Include legend
legend(x="topleft",legend=c("Arid","Semi-arid","Dry sub-humid", "Humid"),
       pch=21,pt.cex=1.5,pt.bg=c("#590007", "#D29773", "#001260", "#71A7C4") ,bty="n")


head_pca_mat <- head_pca_results$scores[,1:2]*-1

# Plot head shape
phylomorphospace(anilios_tree, head_pca_mat, 
                 label = "horizontal", fsize = 1, ftype = 'i', bty='l',
                 xlab = paste("Head PC1 (", head_pca_results$contrib[1], "%)", sep=""),
                 ylab = paste("Head PC2 (", head_pca_results$contrib[2], "%)", sep=""), 
                 xlim = c(-0.2,0.2), ylim = c(-0.11, 0.05),
                 # xlim = range(head_pca_mat[,1]), ylim = range(head_pca_mat[,2])
                 # xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2),
                 pch = 20, cex = 0, node.size = c(1, 1.4),
                 control = list(col.node = tip_col_arid_phy))
points(x = head_pca_mat[,1], y = head_pca_mat[,2], col = stateblindsnake[rownames(head_pca_mat)], pch =20, cex = 1.4)

dev.off()


## Combine plots together

# pdf(file = "output/PC1_pred.pdf", width = 11, height = 8.5)
par(mfrow=c(2,2))


pdf(file = "output/evolution_PC1_pred.pdf", width = 13.33, height = 7.5)
par(mfrow=c(1,3))
plot(body_pca_results$scores[,1] ~ anilios_data$mean_bulk, bty ="n", pch=20, cex = 1.4,
     xlab = "Soil density", ylab = "Body PC1")
# text(body_pca_results$scores[,1] ~ anilios_data$mean_bulk, labels = names(body_pca_results$scores[,1]))

plot(body_pca_results$scores[,1] ~ anilios_data$ARID, bty ="n", pch=20, cex = 1.4,
     xlab = "Aridity", ylab = "Body PC1")
# text(body_pca_results$scores[,1] ~ anilios_data$ARID, labels = names(body_pca_results$scores[,1]))

plot(body_pca_results$scores[,1] ~ anilios_data$temp_mean, bty ="n", pch=20, cex = 1.4,
     xlab = "Temperature", ylab = "Body PC1")
# text(body_pca_results$scores[,1] ~ anilios_data$temp_mean, labels = names(body_pca_results$scores[,1]))
dev.off()

pdf(file = "output/evolution_head_PC1_pred.pdf", width = 13.33, height = 7.5)
par(mfrow=c(1,3))
plot(-head_pca_results$scores[,1] ~ anilios_data$mean_bulk, bty ="n", pch=20, cex = 1.4,
     xlab = "Soil density", ylab = "Head PC1")
# text(-head_pca_results$scores[,1] ~ anilios_data$mean_bulk, labels = names(head_pca_results$scores[,1]))

plot(-head_pca_results$scores[,1] ~ anilios_data$ARID, bty ="n", pch=20, cex = 1.4,
     xlab = "Aridity", ylab = "Head PC1")

plot(-head_pca_results$scores[,1] ~ anilios_data$temp_mean, bty ="n", pch=20, cex = 1.4,
     xlab = "Temperature", ylab = "Head PC1")
# text(-head_pca_results$scores[,1] ~ anilios_data$temp_mean, labels = names(head_pca_results$scores[,1]))
dev.off()

plot(body_pca_results$scores[,1], head_pca_results$scores[,1],
     ylab = "Body PC1", xlab = "Head PC1")


par(mfrow=c(4,3))

# Relative body width by body length
plot(log(max_svl) ~ mean_bulk, data=blindsnake_data, bty="n", pch = 20, cex =1.4)
plot(log(max_svl) ~ ARID, data=blindsnake_data, bty="n", pch = 20, cex =1.4)
plot(log(max_svl) ~ temp_mean, data=blindsnake_data, bty="n", pch = 20, cex =1.4)

plot(sh_mbd ~ mean_bulk, data=blindsnake_data, bty="n", pch = 20, cex =1.4)
plot(sh_mbd ~ ARID, data=blindsnake_data, bty="n", pch = 20, cex =1.4)
plot(sh_mbd ~ temp_mean, data=blindsnake_data, bty="n", pch = 20, cex =1.4)

# Head shape by traits
plot(mean_PC1 ~ mean_bulk, data=blindsnake_data, bty="n", pch = 20, cex =1.4)
plot(mean_PC1 ~ ARID, data=blindsnake_data, bty="n", pch = 20, cex =1.4)
plot(mean_PC1 ~ temp_mean, data=blindsnake_data, bty="n", pch = 20, cex =1.4)

# Head shape by traits
plot(mean_PC2 ~ mean_bulk, data=blindsnake_data, bty="n", pch = 20, cex =1.4)
plot(mean_PC2 ~ ARID, data=blindsnake_data, bty="n", pch = 20, cex =1.4)
plot(mean_PC2 ~ temp_mean, data=blindsnake_data, bty="n", pch = 20, cex =1.4)

# Head shape by body shape
plot(mean_PC1 ~ log_svl, data=blindsnake_data, bty="n", pch = 20, cex =1.4)
plot(mean_PC2 ~ log_svl, data=blindsnake_data, bty="n", pch = 20, cex =1.4)

plot(mean_PC1 ~ sh_mbd, data=blindsnake_data, bty="n", pch = 20, cex =1.4)
plot(mean_PC2 ~ sh_mbd, data=blindsnake_data, bty="n", pch = 20, cex =1.4)
