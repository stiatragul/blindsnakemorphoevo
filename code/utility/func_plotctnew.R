# Modified version of plotCt from convevol

#'Plots calcConv or convSig output.
#'
#'plotCt Plots calcConv or convSig output.
#'
#'@param output object containing calcConv or convSig output
#'@param phy The time calibrated phylogeny of interest in phylo format
#'@param focaltaxa a vector of tip labels for the putatively convergent taxa to be compared
#'@param nsim number of null simulations to plot
#'@param col vector of colors to use for all unique intergroup comparisons a default option is given usable with up to five groups. If number of groups is 1 or less than length of col, not all colors will be used
#'@param groups an optional vector of groups with names matching focaltaxa, indicating the group identity of all focaltaxa
#'@param ... optional arguments to be passed to tiplabels
#'
#'@details Creates a plot that shows the phenotypic distances between pairs of putatively convergent lineages over time. When these distances decrease, convergence has occurred. When more than two putatively
#'convergent taxa are analyzed, all pairs are plotted.
#'
#'@return A plot identifying putatively convergent taxa in the provided phylogeny and tracking the change in phenotypic distance between taxa since their most recent common ancestor
#'
#'@import phytools stats
#'
#'@importFrom utils combn
#'@importFrom graphics par legend
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
#'
#\donttest{'# create time calibrated tree
#'phy<-rcoal(100)
#'
#'# create three normally distributed phenotypic traits
#'traits <- cbind(rnorm(Ntip(phy)),rnorm(Ntip(phy)),rnorm(Ntip(phy)))
#'rownames(traits) <- phy$tip.label
#'focaltaxa <- sample(phy$tip.label, 5)
#'
#'system.time(run2 <- convSigCt(phy, traits, focaltaxa, nsim=10))
#'
#'plotCt(output = run2,phy = phy,focaltaxa = focaltaxa)}


plotCt.new <- function(output, phy, focaltaxa, nsim = 25, .ylabel=NULL,
                       col = c("black","forest green","dodgerblue2","firebrick1","purple","orange","salmon","goldenrod","springgreen2","plum1"), groups = NULL,...){
  measured<-output$meas.path
  
  if ("meas.Cmat" %in% attributes(output)$names){
    type = "simulated"
  } else {type = "non-simulated"}
  if (nsim == 0) {type = "non-simulated"}
  
  if(!is.null(groups)){
    tipNums <- as.character(1:Ntip(phy))
    names(tipNums) <- phy$tip.label
    names(groups)<-tipNums[names(groups)]
    
    path_col_list<-lapply(measured,function(x) groups[x[x["path"] == "tip","node"]])
    path_col_list<-lapply(path_col_list,function(x) x[order(x)])
    path_col<-unlist(lapply(path_col_list,function(x) paste(x,collapse = "_")))
    
    col<-rep(col,length.out = ncol(combn(unique(groups),2)))
    COL<-setNames(col,unique(path_col))
    clr<-COL[match(path_col,names(COL))]
  } else {clr<-rep(col[1],length(measured))
  COL<-setNames(col[1],"MEASURED")
  }
  
  
  start<-measured[[1]]
  range<-max(na.omit(unlist(sapply(measured,"[[",ncol(start)))))
  
  if(type == "simulated"){
    range<-c(range,max(na.omit(unlist(sapply(output$sim.path,"[[",ncol(start))))))
    xmin<-min(unlist(lapply(output$meas.path,function(x) min(x[!x["path"] == "mrca","height"]-max(nodeHeights(phy))))))
  } else {xmin<-min(unlist(lapply(output$meas.path,function(x) min(x[!x["path"] == "mrca","height"]-max(nodeHeights(phy))))))}
  
  lim<-max(range)
  # par(mfrow = c(1,2))
  
  # edge_paths<- lapply(output$meas.path,function(x) which(phy$edge[,2] %in% x[["node"]]))
  # edge_paths<- lapply(edge_paths,function(x) x[x > min(x)])
  # edge_width<- rep(1,nrow(phy$edge))
  # edge_width[unique(unlist(edge_paths))]<-2.5
  # edge_color<- rep("dark grey",nrow(phy$edge))
  # edge_color[unique(unlist(edge_paths))]<-"black"
  # 
  # plot(phy,cex=0.01, edge.width = edge_width, edge.color = edge_color)
  # tiplabels(pie = rep(1,length(focaltaxa)),tip = which(phy$tip.label %in% focaltaxa),cex=0.25)
  # if(!is.null(groups)){
  #   grp.index<-match(unique(groups),groups)
  #   tiplabels(tip = as.numeric(names(groups[grp.index])),text = groups[grp.index],...)
  # }
  start$height<-start$height - max(nodeHeights(phy))
  plot(anc.diff~height,data=start[order(start$height),],type="l",pch=16,xlab="Time (Ma)", bty = 'n',
       ylab= paste(.ylabel, ""), 
       ylim = c(0,lim),xlim = c(xmin,0))
  
  if(type == "simulated"){
    if (nsim > length(output$sim.path)) {nsim <- length(output$sim.path)}
    sim<-sample(output$sim.path,nsim)
    
    for(i in 1:length(sim)){
      y<-sim[[i]]
      y<-y[order(y$height),]
      y$height<-y$height - max(nodeHeights(phy))
      points(y$height,y$anc.diff,type="l", col = "grey", lty = 3)
    }
  }
  
  for(i in 1:length(measured)){
    y<-measured[[i]]
    y<-y[order(y$height),]
    y$height<-y$height - max(nodeHeights(phy))
    points(y$height,y$anc.diff,type="l",col = clr[i],lwd=2)
    points(y$height,y$anc.diff,pch=16,cex = 0.67,col = clr[i])
  }
  
  # legend("topleft",c(names(COL),"Null Simulations"),col = c(COL,"grey"), lty = c(rep(1,length(COL)),3),lwd=2)
}


plotCt.phylo <- function(output, phy, focaltaxa, nsim = 25, .ylabel=NULL, .edge = "black",
                       col = c("black","forest green","dodgerblue2","firebrick1","purple","orange","salmon","goldenrod","springgreen2","plum1"), groups = NULL,...){
  measured<-output$meas.path
  
  if ("meas.Cmat" %in% attributes(output)$names){
    type = "simulated"
  } else {type = "non-simulated"}
  if (nsim == 0) {type = "non-simulated"}
  
  if(!is.null(groups)){
    tipNums <- as.character(1:Ntip(phy))
    names(tipNums) <- phy$tip.label
    names(groups)<-tipNums[names(groups)]
    
    path_col_list<-lapply(measured,function(x) groups[x[x["path"] == "tip","node"]])
    path_col_list<-lapply(path_col_list,function(x) x[order(x)])
    path_col<-unlist(lapply(path_col_list,function(x) paste(x,collapse = "_")))
    
    col<-rep(col,length.out = ncol(combn(unique(groups),2)))
    COL<-setNames(col,unique(path_col))
    clr<-COL[match(path_col,names(COL))]
  } else {clr<-rep(col[1],length(measured))
  COL<-setNames(col[1],"MEASURED")
  }
  
  
  start<-measured[[1]]
  range<-max(na.omit(unlist(sapply(measured,"[[",ncol(start)))))
  
  if(type == "simulated"){
    range<-c(range,max(na.omit(unlist(sapply(output$sim.path,"[[",ncol(start))))))
    xmin<-min(unlist(lapply(output$meas.path,function(x) min(x[!x["path"] == "mrca","height"]-max(nodeHeights(phy))))))
  } else {xmin<-min(unlist(lapply(output$meas.path,function(x) min(x[!x["path"] == "mrca","height"]-max(nodeHeights(phy))))))}
  
  lim<-max(range)
  # par(mfrow = c(1,2))
  
  edge_paths<- lapply(output$meas.path,function(x) which(phy$edge[,2] %in% x[["node"]]))
  edge_paths<- lapply(edge_paths,function(x) x[x > min(x)])
  edge_width<- rep(1,nrow(phy$edge))
  edge_width[unique(unlist(edge_paths))]<-2.5
  edge_color<- rep("dark grey",nrow(phy$edge))
  edge_color[unique(unlist(edge_paths))]<-"black"

  plot(phy, show.tip.label=FALSE, edge.width = edge_width, edge.color = .edge)
  # tiplabels(pie = rep(1,length(focaltaxa)),tip = which(phy$tip.label %in% focaltaxa),cex=0.25)
  if(!is.null(groups)){
    grp.index<-match(unique(groups),groups)
    tiplabels(tip = as.numeric(names(groups[grp.index])),text = groups[grp.index],...)
  }
  start$height<-start$height - max(nodeHeights(phy))
  # plot(anc.diff~height,data=start[order(start$height),],type="l",pch=16,xlab="Time (Ma)", bty = 'n',
  #      ylab= paste(.ylabel, "distance between lineages"), 
  #      ylim = c(0,lim),xlim = c(xmin,0))
  # 
  # if(type == "simulated"){
  #   if (nsim > length(output$sim.path)) {nsim <- length(output$sim.path)}
  #   sim<-sample(output$sim.path,nsim)
  #   
  #   for(i in 1:length(sim)){
  #     y<-sim[[i]]
  #     y<-y[order(y$height),]
  #     y$height<-y$height - max(nodeHeights(phy))
  #     points(y$height,y$anc.diff,type="l", col = "grey", lty = 3)
  #   }
  # }
  
  # for(i in 1:length(measured)){
  #   y<-measured[[i]]
  #   y<-y[order(y$height),]
  #   y$height<-y$height - max(nodeHeights(phy))
  #   points(y$height,y$anc.diff,type="l",col = clr[i],lwd=2)
  #   points(y$height,y$anc.diff,pch=16,cex = 0.67,col = clr[i])
}
  