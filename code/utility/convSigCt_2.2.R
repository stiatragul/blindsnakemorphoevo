#'Computes and conducts significance tests on Ct-metric scores for putatively convergent tips (or groups of tips) given a set of user provided phenotypic characters and a time calibrated phylogeny.
#'
#'calcConvCt Computes and conducts significance tests on Ct-metric scores for putatively convergent tips (or groups of tips) given a set of user provided phenotypic characters and a time calibrated phylogeny.
#'
#'@param phy The time calibrated phylogeny of interest in phylo format
#'@param traits a matrix of numeric phenotypic traits with rownames matching tip labels of phy
#'@param focaltaxa a vector of tip labels for the putatively convergent taxa to be compared
#'@param groups an optional vector of groups with names matching focaltaxa. Indicates the group identity of all putatively convergent taxa and limits Ct measures to intergroup comparisons only
#'@param user.ace a matrix of user supplied ancestral trait values at internal nodes (formatted as "traits" but with node number as rownames)
#'@param nsim number of simulated (Brownian motion) datasets used to build the null distribution
#'@param ... optional arguments to be passed to calcConvCt. If conservative == TRUE, Dmax.t will be restricted to occurr before the oldest stem lineage of the two groups involved in each pairwise comparison. Stem lineage age for each group is defined as the height of the parent node of the branch subtending the most recent common ancestor of tips within a group. Where groups include a single tip, the parent node of the tip's subtending branch is used. Requires group object to be provided by user.. If VERBOSE is TRUE, model information will be printed during computation, including time limits imposed on Dmax.t if the conservative option is chosen.
#'
#'@details Function incorporates the optimizations introduced by Zelditch et al. (2017), which significantly improve runtimes
#'@details Reconstructions part way along branches are obtained using equation [2] of Felsenstein (1985), following code modified from the phytools (Revell, 2012) function contMap
#'
#'@return A list of the following components:
#'@return pvals a matrix containing Ct1 - Ct4 and p-values from significance tests for each
#'@return meas.Cmat a matrix of Ct values for each pairwise comparison of focaltaxa
#'@return meas.path a list of dataframes, one per pairwise comparison of focaltaxa, each containing information from all timepoint measurements of the two putatively convergent lineages. These provide the nodes at which comparisons were drawn, the evolutionary path along which that node fell (i.e., leading to one of two tips), the node height, reconstructed ancestral states at that node for each phenotypic trait, reconstructed ancestral values for each trait along the opposite path, and the phenotypic distance between the two lineages at that point.
#'@return sim.avg average Ct values from all pairwise comparisons between focaltaxa using simulated Brownian motion traits, number of columns corresponds to the user provided number of simulations
#'@return sim.path a list of dataframes as in meas.path, but obtained using simulated data. Length of object determined by number of pairwise comparisons multiplied by the number of simulated datasets.
#'@return grp.mean a matrix of Ct-metrics summarized for inter-group comparisons, returned only if user defined groups were specified. Provides overall results matching those reported in "mean", results for each unique inter-group comparison, and results averaged with equal weight given to each unique inter-group comparison (i.e., ignoring differences in the number of tips in each group).
#'@return grp.pvals a matrix of p-values associated with Ct values in grp.mean object. Returned only if user defined groups were specified.
#'@return limits a list of tree heights used to constrain Dmax.t calculations for each pairwise comparison in conservative analyses. Only returned if conservative == TRUE.
#'
#'@import geiger MASS phytools
#'
#'@export
#'
#'@references Grossnickle DM, Brightly WH, Weaver LN, Stanchak KE, Roston RA, Pevsner SK, Stayton CT, Polly PD, Law CJ. 2022. A cautionary note on quantitative measures of phenotypic convergence. in revision
#'Zelditch ML, Ye J, Mitchell JS, Swiderski DL. 2017. Rare ecomorphological convergence on a complex adaptive landscape: Body size and diet mediate evolution of jaw shape in squirrels (Sciuridae). Evolution 71: 633-649
#'Stayton CT. 2015. The definition, recognition, and interpretation of convergent evolution and two new measures for quantifying and assessing the significance of convergence. Evolution 69(8): 2140-2153.
#'Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217-223.
#'Felsenstein, J. 1985. Phylogenies and the comparative method. American Naturalist, 125, 1-15.
#'
#'@examples
#'library(phytools)
#'library(geiger)
#'
#'# create time calibrated tree
#'mytree<-rtree(100)
#'mycalibration <- makeChronosCalib(mytree, node="root", age.max=50)
#'phy <- chronos(mytree, lambda = 1, model = "correlated", calibration = mycalibration, control = chronos.control() )
#'class(phy)<-"phylo"
#'
#'# create three normally distributed phenotypic traits
#'traits <- cbind(rnorm(Ntip(phy)),rnorm(Ntip(phy)),rnorm(Ntip(phy)))
#'rownames(traits) <- phy$tip.label
#'
#'#	select two random tips, excluding sister taxa
#'pairs <- apply(combn(phy$tip.label,2),2,function(x) nodepath(phy,which(phy$tip.label == x[1]),which(phy$tip.label == x[2])))
#'nosis <- combn(phy$tip.label,2)[,unlist(lapply(pairs, function(x) length(x) > 3))]
#'focaltaxa <- nosis[,sample(1:ncol(nosis),1)]
#'
#'system.time(run <- calcConvCt(phy, traits, focaltaxa))
#'system.time(run2 <- convSigCt(phy, traits, focaltaxa, nsim=100))
#'
#'plot.C(output = run2,phy = phy,focaltaxa = focaltaxa)

convSigCt <- function (phy, traits, focaltaxa, groups = NULL, user.ace = NULL, nsim = 1000, ...) {	## added user.ace object to supply user defined ancestral states to calcConvCt
    data <- calcConvCt(phy, traits, focaltaxa, groups, ace = user.ace, ...)							## added definition of ace to calcConvCt
    phylMat <- vcv.phylo(phy)
    phylMat2 <- phyl.vcv(traits, phylMat, 0)
    simDat <- sim.char(phy, phylMat2$R, nsim, model = "BM", 
        root = 0)
    simOut <- apply(simDat, 3, calcConvCt, phy = phy, focaltaxa = focaltaxa, 
        groups = groups, VERBOSE = FALSE, ...)
    simAvg <- sapply(simOut, "[[", 1)
    pvals <- sapply(1:4, function(x) length(which(simAvg[x, ] >= 
        data[["mean"]][x])))/nsim
    out <- cbind(data[["mean"]], pvals)
    to_return <- list(pvals = out, meas.Cmat = data$Cmat, meas.path = data$meas.path, 
        sim.avg = simAvg, sim.path = do.call(c, sapply(simOut, 
            "[", 3)))
    if (!is.null(groups)) {
        grpAvg <- lapply(simOut, "[[", 4)
        grp.pvals <- matrix(sapply(1:prod(dim(data$grp.mean)), 
            function(x) length(which(lapply(grpAvg, "[[", 
                x) >= data$grp.mean[[x]]))/nsim), 4, ncol(data$grp.mean))
        rownames(grp.pvals) <- rownames(data$grp.mean)
        colnames(grp.pvals) <- colnames(data$grp.mean)
        to_return[["grp.mean"]] <- data$grp.mean
        to_return[["grp.pvals"]] <- grp.pvals
        print(list(round(to_return$grp.mean, 3), round(to_return$grp.pvals, 
            3)))
    }
    else {
        print(round(to_return$pvals, 3))
    }
	if(!is.null(data[["limits"]])){to_return[["limits"]] <- data[["limits"]]}
    return(to_return)
}
