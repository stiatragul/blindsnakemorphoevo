################
# Fitting the model only a single time appears a bit unpredictable
# The final implementation should fit the models multiple times, I'll keep them up here:

# Let's make a function that uses mclapply to fit the model a number of times with different starting params
# start by creating sets of plausible starting parameters from across the surface

library(parallel)

search.surface.DAISIE <- function(model, n.iter = 10, n.proc = 8, nsamples = NULL, datatype="single",
                           g1.params=c(lambda_c=1, mu=0.5, K=30, immigration=0.01, lambda_a=0.1),
                           g2.params=c("lambda_c", "mu", "K", "immigration", "lambda_a"), 
                           results=c("best", "all"), start.params=NULL, 
                           id.parsopt, dd.model, pars.fix, id.parsfix, id.parsnoshift, no.types) {
  beginning <- Sys.time()
  
  if(is.null(nsamples)) stop("please provide the number of taxa in this model as the nsamples value")
  
  if(is.null(start.params)){
    init.params <- lapply(1:n.iter, function(x) {
      c(if(!"lambda_c" %in% g1.params){runif(1, min = 0.001, max = 1)}, # cladogenetic speciation (lambda_c)
        if(!"mu" %in% g1.params){runif(1, min = 0.001, max = 1)}, # extinction (mu)
        if(!"K" %in% g1.params){runif(1, min = length(model[[2]]$branching_times)/2, 
                                                max = length(model[[2]]$branching_times)*2)}, # carrying capacity (K), this could be generalized (e.g. length(branching.times)+10 or something)
        if(!"immigration" %in% g1.params){runif(1, min = 0.0001, max = 1)}, # immigration rate (gamma)
        if(!"lambda_a" %in% g1.params){runif(1, min = 0.0001, max = 1)}, # anagenetic speciation (lambda_a)
        if("lambda_c" %in% g2.params){runif(1, min = 0.001, max = 1)}, # cladogenetic speciation group2 (lambda_c2)
        if("mu" %in% g2.params){runif(1, min = 0.001, max = 1)}, # extinction group2 (mu2)
        if("K" %in% g2.params){runif(1, min = 5, max = 10)}, # carrying capacity group2 (K2)
        if("immigration" %in% g2.params){runif(1, min = 0.0001, max = 1)}, # immigration rate group2 (gamma2)
        if("lambda_a" %in% g2.params){runif(1, min = 0.0001, max = 1)} # anagenetic speciation group2 (lambda_a2)
        )
    })
  } else { init.params <- rep(list(start.params), n.iter)}
  
  for (k in 1:length(init.params)){init.params[[k]] <- init.params[[k]][id.parsopt]}
  
  if(no.types == 1) {
    res.list <- mclapply(1:n.iter, function(x) {
      DAISIE_ML(datalist = model,
                datatype = 'single',
                ddmodel = dd.model,
                initparsopt = init.params[[x]],
                #initparsopt = unlist(eqr_noDD_0[c(1,2,4,5)]), # leaving this here to show that you have to unlist parameters from a previous model fit, annoying!
                idparsopt = id.parsopt,
                parsfix = pars.fix,
                idparsfix = id.parsfix)}, mc.cores = n.proc)
  }
  
  if(no.types == 2) {
    res.list <- mclapply(1:n.iter, function(x) {
      DAISIE_ML(datalist = model,
                datatype = 'single',
                ddmodel = dd.model,
                initparsopt = init.params[[x]],
                idparsopt = id.parsopt,
                parsfix = pars.fix,
                idparsfix = id.parsfix,
                idparsnoshift = id.parsnoshift)}, mc.cores = n.proc)
  }
  
  
  res.list <- res.list[order(sapply(res.list, '[[', 1))] # sort the model fits if you'd like
  
  res.conv <- unlist(lapply(res.list, function(x) x$conv)) # make a vector of the convergence diagnostics, so we can get the index number to remove bad ones
  res.list[which(!res.conv==0)] <- NULL # remove any model fits that didn't converge
  
  if(length(res.list)==0) {stop("No Models Converged!")}
  
  res.values <- unlist(lapply(res.list, function(x) x$loglik)) # make a vector of the values, so we can get the index number of the best
  best.res <- res.list[[which.max(res.values)]]
  
  ## THIS STEP IS GIVING ME ISSUES, AND NEEDS TO BE RESOLVED BEFORE IT CAN BE IMPLEMENTED! ##
  # now use the estimated parameters from the best fit to fit the model again and check for an improvement
#  redo.res <- DAISIE_ML(datalist = model,
#                        ddmodel = dd.model,
#                        initparsopt = unlist(best.res[1:length(id.parsopt)]),
#                        idparsopt = id.parsopt,
#                        parsfix = pars.fix,
#                        idparsfix = id.parsfix)
#  if(redo.res$conv == -1){best.res <- best.res} # remove any model fits that didn't converge
#  else if(redo.res$loglik > best.res$loglik){best.res <- redo.res}
  
  model$nsamples <- nsamples
  
  end <- Sys.time()
  duration <- format(end-beginning)
  print(paste("Computation time to fit the model from", n.iter, "starting points:", duration))
  system("say beep")
  
  if(results=="all"){return(list(input.model = model,
                                 start.params = init.params, 
                                 all.results = res.list, 
                                 best.result = best.res))}
  else if(results=="best"){return(best.res)}
  
}
# The function lets us control a few things via the commands:
# model: just the model of interest, you have to have built it already
# n.iter: the number of model fittings you'd like completed, defaults to 10
# traits: the input traits for model fitting
# n.proc: number of processors. this function will fit the model in parallel.
# no.S: defaults to 1. if your model requires estimating/fitting more than 1 S parameter, say so
# results: would you like the function to report just the best fitting run, or all the results from each fit attempt