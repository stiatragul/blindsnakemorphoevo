# Repository for "Morphological and niche evolution dynamics support a nonadaptive mode of radiation across Australian blindsnakes"

---

These data scripts were used to perform analyses included in the research paper "Morphological and niche evolution dynamics support a nonadaptive mode of radiation across Australian blindsnakes" 

Main questions for the study:

1. What are the main axes of morphological variation?
1. Does variation in morphology among species correlate with their current environments? 
1. Are lineages that occupy ecologically similar habitats morphologically convergent? 
1. Is speciation predominantly allopatric or sympatric? 
1. Do sister species have greater morphological and ecological niche overlap than expected relative to non-sister species pairs?

## Data structure

Contents in the data folder is archived as a zip and can be downloaded from Zenodo (for all versions see https://zenodo.org/doi/10.5281/zenodo.10397830). Once you unzip the zipped files, you will see three folders and some files that are no in any folders. 

**/data/** - files that were manually created and the phylogeny

**/data/script_generated_data/** - A combination of processed data needed to run the analyses

**/data/dorsal/** - photographs of the head from the dorsal view. These photos were used for placing landmarks and semilandmarks. 

**/data/worldclim2_30s/** - cropped and merged annual temperature from WorldClim2 (Fick and Hijmans 2017), soil bulk density from Soil and Landscape Grid of Australia, and Global Aridity Index from Zomer et al. (2022). 

**DREaD/** - contains some files required to replicate DREaD analysis

## Code/Software

All scripts can be run using open source software. Scripts should be run in order to create necessary files that will be saved in /data/script_generated_data/ for further scripts. R is required to run R scripts (.R).

### /Code

  - utility/*.R - scripts for custom functions. These are sourced in other scripts.
  - DREaD/*.R - scripts associated with DREaD analyses
  - 00_linear_measurement_shaperatio.R - script used to account for sexual dimorphism and calculate conventional PCA. Addresses Q1.
  - 01_model_fitting.R - script used to address Q2 and plot visualisations.
  - 02_convergence.R - this script calculates Ct1-4 and C5 scores. Addresses Q3.
  - 02_convergence_model_fitting.R - this script evaluates fit of different evolutionary models to traits. Addresses Q3.
  - 02_convergence_test_simulations.R - simulation studies to show that our phylogeny has sufficient power to detect convergence.
  - 03_niche_enmtools_bias_account.R - calculates ecological niche models (ENMs) for each species using MAXENT. Runs Age-Overlap Correlation tests for geography and ENMs. Partially addresses Q4.
  - 03_DREaD_Blindsnakes_AS.R - script to run DREaD analysis. 
  - 03_morpho_niche_overlap_plots.R - Runs Age-Overlap Correlation tests for body shape and snout shape. Plots AOCs. Partially addresses Q4. 
  - 04_pairwise_distance_test.R - Binomial tests between sister and non-sister pairs for ENMs and Geographic Range. Partially addresses Q5
  - 04_morpho_pairwise.R -  Binomial tests between sister and non-sister pairs for body shape and snout shape. Partially addresses Q5

## Contact

Should you have questions about these scripts or would like to request raw data, please do not hesitate to contact Sarin Tiatragul (contact information can be found in the paper) or on Github (https://github.com/stiatragul)
