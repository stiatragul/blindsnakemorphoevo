enmtools.aoc <- function (clade, env = NULL, overlap.source, nreps = 100, f = NULL, 
          overlap.matrix = NULL, metric = "D") 
{
  description <- "Age-Overlap Correlation from Monte Carlo Test"
  clade <- check.clade(clade)
  enmtools.aoc.precheck(clade, nreps, overlap.source, env, 
                        f, overlap.matrix, metric)
  if (overlap.source == "bc") {
    empirical.models <- lapply(clade$species, function(x) enmtools.bc(x, 
                                                                      env = env))
  }
  if (overlap.source == "dm") {
    empirical.models <- lapply(clade$species, function(x) enmtools.dm(x, 
                                                                      env = env))
  }
  if (overlap.source == "gam") {
    empirical.models <- lapply(clade$species, function(x) enmtools.gam(x, 
                                                                       env = env, f = f))
  }
  if (overlap.source == "glm") {
    empirical.models <- lapply(clade$species, function(x) enmtools.glm(x, 
                                                                       env = env, f = f))
  }
  if (overlap.source == "mx") {
    empirical.models <- lapply(clade$species, function(x) enmtools.maxent(x, 
                                                                          env = env))
  }
  if (is.na(match(overlap.source, c("range", "points", "matrix")))) {
    if (grepl(pattern = "^env", metric)) {
      overlap <- sapply(empirical.models, function(x) sapply(empirical.models, 
                                                             function(y) env.overlap(x, y, env = env)[metric]))
    }
    else {
      overlap <- sapply(empirical.models, function(x) sapply(empirical.models, 
                                                             function(y) raster.overlap(x, y)[metric]))
    }
  }
  if (overlap.source == "range") {
    overlap <- sapply(clade$species, function(x) sapply(clade$species, 
                                                        function(y) geog.range.overlap(x, y)))
  }
  if (overlap.source == "points") {
    overlap <- sapply(clade$species, function(x) sapply(clade$species, 
                                                        function(y) point.overlap(x, y)))
  }
  if (overlap.source == "matrix") {
    overlap = overlap.matrix
  }
  tree <- clade$tree
  tree$node.label <- NULL
  rownames(overlap) <- colnames(overlap)
  empirical.df <- node.overlap(overlap, tree)
  empirical.model <- lm(empirical.df$overlap ~ empirical.df$age)
  do.rep <- function(inds) {
    tree$tip.label <- tree$tip.label[inds]
    rep.df <- node.overlap(overlap, tree)
    return(list(rep.df = rep.df, rep.lm = lm(rep.df$overlap ~ 
                                               rep.df$age)))
  }
  reps <- list()
  for (i in 1:nreps) {

    this.rep <- sample(nrow(overlap))
    reps[[paste0("rep.", i)]] <- do.rep(sample(length(tree$tip.label)))$rep.lm
  }
  reps.aoc <- rbind(empirical.model$coefficients, do.call(rbind, 
                                                          lapply(reps, function(x) x$coefficients)))
  rownames(reps.aoc) <- c("empirical", paste("rep", 1:nreps))
  p.values <- apply(reps.aoc, 2, function(x) min(rank(x)[1], 
                                                 rank(-x)[1])/length(x))
                                                                                                                   
  output <- list(coefficients = reps.aoc, p.values = p.values, 
                 tree = tree, empirical.overlap = overlap, 
                 empirical.df = empirical.df, empirical.model = empirical.model, 
                 reps = reps)
  class(output) <- "enmtools.aoc"
  return(output)
}

plot.enmtools.aoc <- function(x, ...){
  
  check.packages("ape")
  
  plot(x$tree, no.margin=TRUE, edge.width=2, cex=1)
  ape::nodelabels(format(x$empirical.df$overlap, digits = 1, nsmall = 2))
  
  grid.arrange(x$regressions.plot, x$intercept.plot,
               x$slope.plot, ncol = 2) +
    theme(plot.title = element_text(hjust = 0.5))
  
}

enmtools.aoc.precheck <- function(clade, nreps, overlap.source, env,  model, overlap.matrix, metric){
  
  if(!inherits(clade$tree, "phylo")){
    stop("Tree is not a phylo object!")
  }
  
  if(is.null(clade$tree$edge.length)){
    stop("Tree does not have branch lengths!")
  }
  
  # Check to make sure env data is good
  if(!is.na(match(overlap.source, c("bc", "dm", "mx", "glm", "gam")))){
    if(!inherits(env, c("raster", "RasterLayer", "RasterStack", "RasterBrick"))){
      stop("No environmental rasters were supplied!")
    }
  }
  
  if(overlap.source == "range"){
    if(any(is.na(lapply(clade$species, function(x) x$range)))){
      
      stop(paste("Overlap source set to range, but some species are missing range rasters: ",
                 paste(names(clade$species)[which(is.na(lapply(clade$species, function(x) x$range)))], collapse = ", ")))
    }
  }
  
}