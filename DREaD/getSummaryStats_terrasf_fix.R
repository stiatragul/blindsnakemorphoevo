
getSummaryStats <- function(phy, clade, sp_rasters){

  # Standardise tree age 
  root.age.true <- max(findMRCA(phy, type=c("height")))
  phy$edge.length<-phy$edge.length/root.age.true
  root.age<-1
  
  # get sister species
  sisters <- findSisters(phy, solve.polytomies=T)
  
  # find divergence times
  sisters.split.time <- sisterSpeciesSplitTimes(phy, sisters, root.age=1)
  
  #get sister species range overlap
  RO_ss <- sisterSpeciesOverlap(sp_rasters, sisters, root.age=1)
  
  #get sister species range asymmetry
  Asym_ss <- sisterSpeciesRangeAsymmetry(sp_rasters, sisters)
  
  # get sister species range distances
  print(paste("running range distance ... takes a min"))
  RD_ss <- sisterSpeciesRangeDistances(sp_rasters, sisters, aggregate_size = 2)
  
  # get sister-species-outgroup overlap metrics
  TO <- outgroupOverlap(sp_rasters, sisters, sisters.overlap = RO_ss)
  
  #! first batch of summary statistics - mean and SD
  range.overlap.av <- mean(RO_ss, na.rm=T)
  range.overlap.sd <- sd(RO_ss, na.rm=T)
  
  range.asymmetry.av <- mean(Asym_ss, na.rm=T)
  range.asymmetry.sd <- sd(Asym_ss, na.rm=T)
  
  range.distance.av <- mean(RD_ss, na.rm=T)
  range.distance.sd <- sd(RD_ss, na.rm=T)
  
  TO.av <- mean(TO, na.rm=T)
  TO.sd <- sd(TO, na.rm=T)
  
  #! second batch of summary statistics - linear models
  if(length(unique(Asym_ss))==1){
    range.asymmetry.slope <- 0
    range.asymmetry.intercept <- 0
  } else{
    lm.asym    <-  lm(Asym_ss ~ sisters.split.time)
    range.asymmetry.slope<-lm.asym$coefficients[[2]]
    range.asymmetry.intercept<-lm.asym$coefficients[[1]]
  }
  
  if(length(unique(RO_ss))==1){
    range.overlap.slope <- 0
    range.overlap.intercept <- 0
  } else{
    lm.overlap <-  lm(RO_ss ~ sisters.split.time)
    range.overlap.slope <- lm.overlap$coefficients[[2]]
    range.overlap.intercept <- lm.overlap$coefficients[[1]]
  }
  
  if(length(unique(RD_ss))==1){
    range.distance.slope <- 0
    range.distance.intercept <- 0
  } else{
    lm.distance <- lm(RD_ss ~ sisters.split.time)
    range.distance.slope <- lm.distance$coefficients[[2]]
    range.distance.intercept <- lm.distance$coefficients[[1]]
  }

  
  #! third batch of summary statistics - species range sizes
  
  species.range.sizes <- rangeSize(sp_rasters)
  species.range.sizes.stand <- species.range.sizes/max(species.range.sizes)
  
  range.size.skew <- skewness(species.range.sizes.stand)
  range.size.av <- mean(species.range.sizes.stand, na.rm=T)
  range.size.sd <- sd(species.range.sizes.stand, na.rm=T)
  
  #! fourth batch of summary statistics - Phylogenetic balance and treeshape metrics
  
  phy_treeshape <- as.treeshape(phy)
  beta <- maxlik.betasplit(phy_treeshape, up = 10, remove.outgroup = FALSE, confidence.interval = "none", conf.level = 0.95, size.bootstrap = 100)$max_lik
  collessI <- colless(phy_treeshape)
  sackinI <- sackin(phy_treeshape)
  ltt <- ltt(phy, plot=FALSE)
  gamma <-ltt$gamma
  
  #! fifth batch of summary statistics - further overlap summary metrics
  
  symp0 <-   length(which(RO_ss == 0))/   length(RO_ss)
  symp50 <-  length(which(RO_ss >= 0.5))/ length(RO_ss)
  symp75 <-  length(which(RO_ss >= 0.75))/length(RO_ss)
  symp90 <-  length(which(RO_ss >= 0.9))/ length(RO_ss)
  symp100 <- length(which(RO_ss == 1))/   length(RO_ss)
  
  range.overlap.kurt<- kurtosis(RO_ss)
  range.overlap.skew <- skewness(RO_ss)
  
  #! sixth batch of summary statistics - biomodality
  sisterpairs <- length(sisters)
  
  bimodality100 <-((symp0*sisterpairs) *  (symp100*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  bimodality90  <-((symp0*sisterpairs) *  (symp90*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  bimodality75  <-((symp0*sisterpairs) *  (symp75*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  bimodality50  <-((symp0*sisterpairs) *  (symp50*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  
  #! seventh batch of summary statistics - ARC
  print(paste("running ARC ... takes a min"))
  
  ARC.clade <- list()
  
  for (i in 1:length(phy$tip.label)){
    enm.sp <- enmtools.species()
    enm.sp$species.name <- phy$tip.label[[i]]
    enm.sp$range <- sp_rasters[[which(names(sp_rasters)==enm.sp$species.name)]]
    ARC.clade$species[[i]] <- enm.sp
  }
  names(ARC.clade$species) <- phy$tip.label
  ARC.clade$tree <- phy
  
  ARC <- enmtools.aoc(clade = ARC.clade,  nreps = 0, overlap.source = "range")
  
  
  ARCint <- ARC$coefficients[[1,1]]#What's this?
  ARCslope <- ARC$coefficients[[1,2]]#What's this?
  
  # put it together
  
  tmp_df <- data.frame(clade         = clade,
                       ntips         = length(phy$tip.label), 
                       ROmean        = range.overlap.av,
                       ROsd          = range.overlap.sd,
                       ROslope       = range.overlap.slope, 
                       ROintercept   = range.overlap.intercept, 
                       ROskew        = range.overlap.skew, 
                       ROkurtosis    = range.overlap.kurt,
                       RO0           = symp0, 
                       RO50          = symp50, 
                       RO75          = symp75, 
                       RO90          = symp90, 
                       RO100         = symp100,
                       TOmean        = TO.av,
                       TOsd          = TO.sd, 
                       asymmean      = range.asymmetry.av, 
                       asymsd        = range.asymmetry.sd,
                       asymslope     = range.asymmetry.slope, 
                       ASYMintercept = range.asymmetry.intercept,
                       RDmean        = range.distance.av, 
                       RDsd          = range.distance.sd, 
                       RDintercept   = range.distance.intercept, 
                       RDslope       = range.distance.slope,
                       BIMOD50       = bimodality50, 
                       BIMODE75      = bimodality75, 
                       BIMOD90       = bimodality90, 
                       BIMOD100      = bimodality100,
                       RSskew        = range.size.skew, 
                       RSmean        = range.size.av, 
                       RSsd          = range.size.sd, 
                       CI            = collessI, 
                       Beta          = beta, 
                       Gamma         = gamma, 
                       SI            = sackinI,
                       ARCslope      = ARCslope, 
                       ARCint        = ARCint) 
  
  return(tmp_df)
}
