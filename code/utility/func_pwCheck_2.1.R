#'Calculates the number of measurements that can be made between two lineages for each pairwise comparison within a set of putatively convergent tips (group identity may also be taken into account). Useful for determining which comparisons are not informative, and constructing a group object before running calcConvCt or convSigCt.
#'
#'@param phy The time calibrated phylogeny of interest in phylo format
#'@param focaltaxa a vector of tip labels for the putatively convergent taxa to be compared
#'@param groups an optional vector of groups with names matching focaltaxa. Indicates the group identity of all putatively convergent taxa and limits Ct measures to intergroup comparisons only
#'
#'@return A list of the following components:
#'@return taxa a matrix of uninformative tip comparisons
#'@return path a vector with the number of measurements for all pairwise comparisons - a summary of this is also printed when running the function
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
#'focaltaxa <- sample(phy$tip.label, 5)
#'
#'pwCheck(phy,focaltaxa)

pwCheck <- function(phy, focaltaxa, groups = NULL) {
	Combinations <- combn(focaltaxa, 2)
	if(!is.null(groups)){
	Combinations <- Combinations[,apply(Combinations,2,function(x) length(unique(groups[x]))>1)]
	}
	
	lng_path <-
	apply(Combinations,2,function(y) {
	sum(unlist(lapply(which(phy$tip.label %in% y),function(x) length(nodepath(phy,getMRCA(phy,y),x)) - 2)))
	})
	
	print(summary(lng_path))
	if(!is.null(groups)){
	return(list("taxa" = Combinations[,which(lng_path == 0)],
				"groups" = apply(Combinations[,which(lng_path == 0)],2,function(x) groups[x]),
				"path" = lng_path
				))
	} else { return(list("taxa" = Combinations[,which(lng_path == 0)],
							"path" = lng_path
							))}
}
