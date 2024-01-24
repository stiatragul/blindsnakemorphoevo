# function to check for sexual dimorphism using lm() and emmeans(

library(emmeans)
library(dplyr)

sex_dimorphic_tester <- function(formula, .data) {
  # Fit the model
  
  fit <- lm(as.formula(formula), data = .data)

  # Calculate the pairwise comparisons
  emm_fit <- emmeans::emmeans(object = fit, specs = ~ sex:species)
  emm_fit_contrast <- emmeans::contrast(emm_fit, "pairwise", by = "species", interaction = TRUE, adjust = "bonferroni")
  
  pairwise_df <- as.data.frame(emm_fit_contrast)
  
  contrast_table <- data.frame(species = pairwise_df$species, esimate = pairwise_df$estimate, SE = pairwise_df$SE, 
                                    t.ratio = pairwise_df$t.ratio, p.value = pairwise_df$p.value)
  
  significant_species <- contrast_table %>%
    filter(p.value < 0.05) %>%
    distinct(species) 
  
  vector_to_subset <- significant_species$species
  
  return(list(table = contrast_table, dimorph_sp = as.character(c(vector_to_subset))))

}
# 
# fit_log_svl <- lm(log_svl ~ sex + species + sex:species, data = subset_dimorph)
# summary(fit_log_svl)
# emm1 <- emmeans::emmeans(object = fit_log_svl, specs = ~ sex:species)
# multivariate_test <- emmeans::contrast(emm1, "pairwise", by = "species", 
#                                        interaction = TRUE, adjust = "bonferroni")
# 
# pairwise_df <- as.data.frame(multivariate_test)
# 
# pair_wise_contrasts <- data.frame(species = pairwise_df$species, esimate = pairwise_df$estimate, SE = pairwise_df$SE, 
#                                   t.ratio = pairwise_df$t.ratio, p.value = pairwise_df$p.value)
# 
# significant_species <- pair_wise_contrasts %>%
#   filter(p.value < 0.05) %>%
#   distinct(species)


# 
# # Summarise function ------------------------------------------------------
# 
# filtR_annotatR <- function(.df, .trait, .exclude_vec, name_trait) {
#   
#   .n
#   df <- filtered_svl %>% dplyr::group_by(species) %>% dplyr::summarise(mean_trait = mean(.trait), n=n())
#   
#   dummy_notes <- PC_log_ratio_df %>% dplyr::group_by(species) %>% dplyr::summarise(n = n()) %>% 
#     dplyr::mutate(notes = ifelse(species %in% .exclude_vec, "females_only", no = "both"))
#   
#   df$notes <- dummy_notes$notes  
#   names(df) <- c("species", paste(name_trait), "n", "notes")
#   
# }
# 
# filtR_annotatR(filtered_svl, svl, .exclude_vec = d_svl$dimorph_sp, name_trait = "svl")
# 
# 
# df <- filtered_svl %>% dplyr::group_by(species) %>% dplyr::summarise(mean_trait = mean(.trait), n=n())
# 
# dummy_notes <- PC_log_ratio_df %>% dplyr::group_by(species) %>% dplyr::summarise(n = n()) %>% 
#   dplyr::mutate(notes = ifelse(species %in% .exclude_vec, "females_only", no = "both"))
# 
# df$notes <- dummy_notes$notes  
# 
# 
