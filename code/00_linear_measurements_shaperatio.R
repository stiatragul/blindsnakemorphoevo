# 00_linear_measurements_shaperatio.R
# Sarin Tiatragul
# Created: Sept 2021
# Modified: Feb 2023

# Purpose:
# Convert linear body measurements to shape ratios 
# Body proportions calculate by dividing all linear variables by body size (svl)

# Libraries ---------------------------------------------------------------

library(ggplot2); library(lme4); library(stringr)
library(factoextra); library(gridExtra);
library(emmeans); library(dplyr)
'%notin%' <- Negate('%in%')

# Read data ---------------------------------------------------------------
linear_data <- read.csv(file = 'data/script_generated_data/linear_measurements/blindsnake_body_traits.csv')

# Rename column names to make it easier
colnames(linear_data)[18] <- "nostril_2_eye"
colnames(linear_data)[19] <- "eye_2_eye"
colnames(linear_data)[20] <- "forehead_2_snout"
colnames(linear_data)[21] <- "lip_2_snout"

linear_data <- linear_data %>% 
  dplyr::filter(!is.na(head_width)) %>% 
  dplyr::filter(!is.na(head_depth)) %>% 
  dplyr::filter(!reg_no %in% c("R.2202")) %>% 
  dplyr::filter(!stringr::str_starts(sex, "j"))

# Number of specimens 
linear_data %>% nrow()
# Number of species
linear_data %>% filter(genus == "Anilios") %>% distinct(species) %>% arrange(species)
# sample size per species
linear_data %>% dplyr::group_by(species) %>% summarise(n = n())
# Specimens included table
specimens_table <- data.frame(institution = linear_data$institution,
                              catalogue_number = linear_data$reg_no, 
                              species = paste(linear_data$genus, linear_data$species),
                              sex = linear_data$sex,
                              longitude = linear_data$longitude,
                              latitude = linear_data$latitude)
specimens_table <- specimens_table %>% arrange(species)

raw_lin_df <- linear_data[, 12:24]
raw_lin_df$headlength_xray %>% summary()

# Ehmann style length / diameter
linear_data$l_d <- linear_data$svl/linear_data$midbody_diameter


# CONVERT TO LOG SHAPE RATIOS ----------------------------------------------

# For our study, we are interested in the body shape of the snake so we will only 
# look at four variables: midbody diameter, midtail width, head depth, and head width
# All of this will be divided by SVL

# So we first select these traits that we are interested in

body_traits <- data.frame(svl = linear_data$svl,
                          ttl = linear_data$total_length,
                          mbd = linear_data$midbody_diameter,
                          mtw = linear_data$midtail_diameter,
                          tl = linear_data$tail_length,
                          hwe = linear_data$head_width,
                          hda = linear_data$tip_of_mouth_to_head_depth,
                          hd = linear_data$head_depth)

# Divide each trait by size
body_ratios <- body_traits[, c("mbd", "mtw", "hwe", "hda", "hd", "tl")] / body_traits$svl

# Log transform body ratio
# log_ratio_df <- log(body_ratios)
log_ratio_df <- log(body_ratios)

# Change column names to reflect that this is shape ratio (sh)
colnames(log_ratio_df) <- paste("sh", colnames(log_ratio_df), sep = "_")

# Combine with other metadata
body_ratio_df <- cbind(linear_data[, 1:13], log_ratio_df)


# PCA OF BODY SHAPE -------------------------------------------------------

# Assign individual id to each row
rownames(log_ratio_df) <- body_ratio_df$reg_no

# Calculate principal component of log body ratios 
PC_log_body_ratio <- prcomp(log_ratio_df, scale = TRUE, center = TRUE)

save(PC_log_body_ratio, file = "data/script_generated_data/linear_measurements/PC_log_ratio_pc_obj.rda")

# Visualise where individuals land on PC axis to check if there's any outlier
factoextra::fviz_pca_ind(PC_log_body_ratio)

# Calculate variance explained by each PC component
factoextra::fviz_eig(PC_log_body_ratio, "variance")

# Visualise contribution from each trait
pc_img1 <- factoextra::fviz_pca_var(PC_log_body_ratio, 
                         axes = c(1, 2),
                         geom = c("arrow", "text"), pointsize = 2,
                         col.var = "contrib", # Color by contribution of PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE     # Avoid text overlapping
); pc_img1

# Axis 2 and 3
factoextra::fviz_pca_var(PC_log_body_ratio, 
                         axes = c(2, 3),
                         geom = c("point", "text"), pointsize = 2,
                         col.var = "contrib", # Color by contribution of PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE     # Avoid text overlapping
)

# Bind PC1, PC2, and PC3 into data frame 
PC_log_ratio_df <- cbind(body_ratio_df, PC_log_body_ratio$x[,1:3])

# PC SCORES SUMMARISE BY SPECIES -----------------------------------------------

# Append species to PC data
PC_log_ratio_df$species <- paste(PC_log_ratio_df$genus, PC_log_ratio_df$species)

# Summarise data per species. Find mean for latitude and longitude also
PC_body_shape_sum <- aggregate(PC_log_ratio_df[, c("svl", "sh_mbd", "sh_mtw", "sh_hwe", "sh_hda", "sh_hd", "sh_tl", "PC1", "PC2", "PC3")], 
                               list(PC_log_ratio_df$species), 
                               FUN = "mean", na.rm = TRUE, na.action = na.omit) 

# Put species in the column name
colnames(PC_body_shape_sum)[1] <- 'species'

PC_log_ratio_df$log_svl <- log(PC_log_ratio_df$svl)
PC_log_ratio_df$log_tail <- log(PC_log_ratio_df$tail_length)

# TEST FOR DIMORPHISM -----------------------------------------------------
# Is there sexual dimorphism in body ratio? 
# Need to filter specimens that don't have sex data

# Check for bias of measurement
PC_log_ratio_df %>% group_by(sex) %>% summarise(n = n())
## We have more females than males in our dataset

PC_log_ratio_df <- PC_log_ratio_df %>% dplyr::filter(sex %notin% c('j', 'juvenile'))

# Make a subset of data where sex has been identified
dimorph_df <- PC_log_ratio_df %>%
  dplyr::filter(sex %in% c('f', 'm')) %>% 
  dplyr::filter(!species %in% c('Anilios aspina', 'Anilios longissimus', 'Anilios margaretae',
                                'Anilios sp', 'Anilios zonula', 'Ramphotyphlops multilineatus', 
                                'Sundatyphlops polygrammicus', 'Anilios fossor', 'Anilios vagurima',
                                'Anilios chamodracaena'))

# Vector of species where sample size is at least 3 for males and females
include_species <- dimorph_df %>% dplyr::group_by(species, sex) %>% dplyr::summarise(n = n()) %>% 
  dplyr::filter(n >= 3) %>% dplyr::distinct(species) %>% c()

# Subset data to test for sexual dimorphism for each trait using only species that have at least 3 males and 3 females
subset_dimorph <- dimorph_df[dimorph_df$species %in% include_species$species,]

# Fit linear models and check which species have sexual dimorphism
source('code/utility/func_sexual_dimorphism_tester.R')
d_svl <- sex_dimorphic_tester("svl ~ sex + species + sex:species", subset_dimorph)
d_log_svl <- sex_dimorphic_tester("log_svl ~ sex + species + sex:species", subset_dimorph)
d_log_tl  <- sex_dimorphic_tester("log_tail ~ sex + species + sex:species", subset_dimorph)
d_sh_mbd <- sex_dimorphic_tester("sh_mbd ~ sex + species + sex:species", subset_dimorph)
d_sh_mtw <- sex_dimorphic_tester("sh_mtw ~ sex + species + sex:species", subset_dimorph)
d_sh_hwe <- sex_dimorphic_tester("sh_hwe ~ sex + species + sex:species", subset_dimorph)
d_sh_hda <- sex_dimorphic_tester("sh_hda ~ sex + species + sex:species", subset_dimorph)
d_sh_hd  <- sex_dimorphic_tester("sh_hd ~ sex + species + sex:species", subset_dimorph)
d_sh_tl  <- sex_dimorphic_tester("sh_tl ~ sex + species + sex:species", subset_dimorph)

# Filter out males from species that show sexual dimorphism because we have more females
filtered_svl <- PC_log_ratio_df %>% filter(!(species %in% d_svl$dimorph_sp & sex == "m")) 
filtered_log_svl <- PC_log_ratio_df %>% filter(!(species %in% d_log_svl$dimorph_sp & sex == "m")) 
filtered_log_tl <- PC_log_ratio_df %>% filter(!(species %in% d_log_tl$dimorph_sp & sex == "m")) 
filtered_sh_mbd <- PC_log_ratio_df %>% filter(!(species %in% d_sh_mbd$dimorph_sp & sex == "m")) 
filtered_sh_mtw <- PC_log_ratio_df %>% filter(!(species %in% d_sh_mtw$dimorph_sp & sex == "m")) 
filtered_sh_hwe <- PC_log_ratio_df %>% filter(!(species %in% d_sh_hwe$dimorph_sp & sex == "m")) 
filtered_sh_hda <- PC_log_ratio_df %>% filter(!(species %in% d_sh_hda$dimorph_sp & sex == "m")) 
filtered_sh_hd <- PC_log_ratio_df %>% filter(!(species %in% d_sh_hd$dimorph_sp  & sex == "m")) 
filtered_sh_tl <- PC_log_ratio_df %>% filter(!(species %in% d_sh_tl$dimorph_sp  & sex == "m")) 

reg_svl <- filtered_svl %>% dplyr::group_by(species) %>% dplyr::summarise(svl = mean(svl), n=n()) %>% dplyr::mutate(sampling = ifelse(species %in% d_svl$dimorph_sp, "females_only", "full"))
log_svl <- filtered_log_svl %>% dplyr::group_by(species) %>% dplyr::summarise(log_svl = mean(log_svl), n=n()) %>% dplyr::mutate(sampling = ifelse(species %in% d_log_svl$dimorph_sp, "females_only", "full"))
log_tl  <- filtered_log_tl  %>% dplyr::group_by(species) %>% dplyr::summarise(log_tl = mean(log_tail), n=n())%>% dplyr::mutate(sampling = ifelse(species %in% d_log_svl$dimorph_sp, "females_only", "full"))
sh_mbd  <- filtered_sh_mbd  %>% dplyr::group_by(species) %>% dplyr::summarise(sh_mbd = mean(sh_mbd), n=n())%>% dplyr::mutate(sampling = ifelse(species %in% d_sh_mbd$dimorph_sp, "females_only", "full"))
sh_mtw  <- filtered_sh_mtw  %>% dplyr::group_by(species) %>% dplyr::summarise(sh_mtw = mean(sh_mtw), n=n())%>% dplyr::mutate(sampling = ifelse(species %in% d_sh_mtw$dimorph_sp, "females_only", "full"))
sh_hwe  <- filtered_sh_hwe  %>% dplyr::group_by(species) %>% dplyr::summarise(sh_hwe = mean(sh_hwe), n=n())%>% dplyr::mutate(sampling = ifelse(species %in% d_sh_hwe$dimorph_sp, "females_only", "full"))
sh_hda  <- filtered_sh_hda  %>% dplyr::group_by(species) %>% dplyr::summarise(sh_hda = mean(sh_hda), n=n())%>% dplyr::mutate(sampling = ifelse(species %in% d_sh_hda$dimorph_sp, "females_only", "full"))
sh_hd   <- filtered_sh_hd %>% dplyr::group_by(species) %>% dplyr::summarise(sh_hd = mean(sh_hd), n=n())%>% dplyr::mutate(sampling = ifelse(species %in% d_sh_hd$dimorph_sp, "females_only", "full"))
sh_tl   <- filtered_sh_tl %>% dplyr::group_by(species) %>% dplyr::summarise(sh_tl = mean(sh_tl), n=n())%>% dplyr::mutate(sampling = ifelse(species %in% d_sh_tl$dimorph_sp, "females_only", "full"))

## Write each of these to csv as supplementary data
df_list <- list(svl = reg_svl, log_svl=log_svl, log_tl = log_tl, sh_mbd = sh_mbd, sh_mtw = sh_mtw, sh_hwe = sh_hwe, sh_hda = sh_hda, sh_hd = sh_hd, sh_tl = sh_tl)

# apply the write.csv() function to each dataframe in the list
lapply(seq_along(df_list), function(i) {
  filename <- paste0("data/script_generated_data/shape_ratios/supp_data", names(df_list)[i], ".csv")
  write.csv(df_list[[i]], filename, row.names = FALSE)
})

# Replace values in PC_body_shape_sum
# We will disregard the PC scores for now because they did not take in account of our dataset
PC_body_shape_sum$svl <- reg_svl$svl
PC_body_shape_sum$log_svl <- log_svl$log_svl
PC_body_shape_sum$log_tl <- log_tl$log_tl
PC_body_shape_sum$sh_mbd <- sh_mbd$sh_mbd
PC_body_shape_sum$sh_mtw <- sh_mtw$sh_mtw
PC_body_shape_sum$sh_hwe <- sh_hwe$sh_hwe
PC_body_shape_sum$sh_hda <- sh_hda$sh_hda
PC_body_shape_sum$sh_hd  <- sh_hd$sh_hd
PC_body_shape_sum$sh_tl  <- sh_tl$sh_tl

# Redo PCA on filtered values

# PCA of mean values ------------------------------------------------------

rownames(PC_body_shape_sum) <- PC_body_shape_sum$species

b_mean_PC <- prcomp(PC_body_shape_sum[,3:8], center = T, scale = T) 
# Visualise where individuals land on PC axis to check if there's any outlier

factoextra::fviz_pca_ind(b_mean_PC, )

# Calculate variance explained by each PC component
factoextra::fviz_eig(b_mean_PC, "variance")

# Visualise contribution from each trait
factoextra::fviz_pca_var(b_mean_PC, 
                         axes = c(1, 2),
                         geom = c("arrow", "text"), pointsize = 2,
                         col.var = "contrib", # Color by contribution of PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE     # Avoid text overlapping
)

save(b_mean_PC, file = "data/script_generated_data/linear_measurements/PC_log_ratio_pc_obj.rda")

filter_mean_PC_scores <- b_mean_PC$x[,1:3]
colnames(filter_mean_PC_scores) <- c("fil_PC1", "fil_PC2", "fil_PC3")
PC_body_shape_sum <- cbind(PC_body_shape_sum, filter_mean_PC_scores)

# Write data to csv -------------------------------------------------------
write.csv(body_ratio_df, file = 'data/script_generated_data/shape_ratios/blindsnake_logbodyshape_ratio.csv', row.names = FALSE, na = "")
write.csv(PC_body_shape_sum, file = 'data/script_generated_data/shape_ratios/blindsnake_sp_pc_logbodyrat.csv', row.names = FALSE, na = "")
write.csv(PC_log_ratio_df, file = 'data/script_generated_data/shape_ratios/log_body_shape_ratio_pc_all.csv', row.names = FALSE, na = "")

# Make table to show dimorphism -------------------------------------------

list_emm_tables <- list(svl = d_svl$table, log_svl = d_log_svl$table, log_tl = d_log_tl$table, sh_mbd = d_sh_mbd$table, 
     sh_mtw = d_sh_mtw$table, sh_hwe = d_sh_hwe$table, sh_hda = d_sh_hda$table, sh_hd = d_sh_hd$table, sh_tl = d_sh_tl$table)

