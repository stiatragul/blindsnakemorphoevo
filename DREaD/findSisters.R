findSisters <- function(phy,solve.polytomies=FALSE,include.polytomies=TRUE){
  require(phytools)
  n <- length(phy$tip.label)
  nb.node <- phy$Nnode
  if (n < 4)
    stop("tree has fewer than 4 tips")
  if(solve.polytomies){
    phy=multi2di(phy)
    include.polytomies=FALSE
  }
  if(include.polytomies){
    sister.ancestors <- which(tabulate(phy$edge[, 1][phy$edge[, 2] <= n])>=2)
    cherries <- vector("list",length(sister.ancestors))
    for(i in 1:length(cherries)){
      x <- getDescendants(phy, sister.ancestors[i])
      cherries[[i]] <- phy$tip.label[x]
      if(length(which(is.na(cherries[[i]])))>0) cherries[[i]]<-NULL  #remove internal nodes
    }
  }
  if(!include.polytomies){
    sister.ancestors <- which(tabulate(phy$edge[, 1][phy$edge[, 2] <= n])==2)
    cherries <- vector("list",length(sister.ancestors))
    for(i in 1:length(cherries)){
      x <- getDescendants(phy, sister.ancestors[i])
      cherries[[i]] <- c(phy$tip.label[x[1]], phy$tip.label[x[2]])
    }
  }
  cherries[!sapply(cherries,is.null)]
  return(cherries)
}

# function to find the species in the outgroup of a sister species pair
findOG <- function(phy,sistips){
  require(phytools)
  mrca <- mrca(phy)
  og <- list()
  for(i in 1:length(sistips)){
    sistip1 <- sistips[[i]][1]
    sistip2 <- sistips[[i]][2]
    curr.mrca <- mrca[which(rownames(mrca)==sistip1),which(colnames(mrca)==sistip2)]
    curr.anc <- phy$edge[,1][which(phy$edge[,2]==curr.mrca)]
    curr.tips <- phy$tip.label[getDescendants(phy,curr.anc)[getDescendants(phy,curr.anc)<Ntip(phy)+1]]
    curr.tips1 <- c(sistip1,sistip2)
    og[[i]] <- setdiff(curr.tips,curr.tips1)
  }
  return(og)
}