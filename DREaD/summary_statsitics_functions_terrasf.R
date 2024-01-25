
sisterSpeciesSplitTimes <- function(phy, sisters, root.age=1){
  
  sisters.split.time<-c()
  
  for (j in 1:length(sisters)) {
    
    sisters.split.time[j]<-(root.age-findMRCA(phy, tips=paste(sisters[[j]]), type=c("height")))
    
  }
  names(sisters.split.time) <- sisters
  return(sisters.split.time)
}

sisterSpeciesOverlap <- function(sp.rasters, sisters, root.age=1){
  
  sisters.overlap<-list()
  
  sisters.split.time<-c()
  
  for (j in 1:length(sisters)) {
    
    sisters.split.time[j]<-(root.age-findMRCA(phy, tips=paste(sisters[[j]]), type=c("height")))
    
    sis.1 <- sp.rasters[[which(names(sp.rasters) == sisters[[j]][1])]]
    
    sis.2 <- sp.rasters[[which(names(sp.rasters) == sisters[[j]][2])]]
    
    sis.1.area <- sum(values(sis.1), na.rm=T)
    
    sis.2.area <- sum(values(sis.2), na.rm=T)
    
    if(sis.1.area == 0) {print(paste(j," species with no range"))}
    if(sis.2.area == 0) {print(paste(j," species with no range"))}
    
    sis.overlap<-sis.1+sis.2
    
    if(any(values(sis.overlap)==2, na.rm=T)) {
      
      if (sis.1.area <= sis.2.area) {
        
        range.overlap <- length(which(values(sis.overlap) == 2))/ sis.1.area
      }
      
      if (sis.1.area > sis.2.area) {
        
        range.overlap <- length(which(values(sis.overlap) == 2))/ sis.2.area
      }
    } else { range.overlap <- 0 }
    
    sisters.overlap[[length(sisters.overlap) + 1]] <- range.overlap
  }
  
  names(sisters.overlap) <- sisters
  sisters.overlap <- unlist(sisters.overlap)
  return(sisters.overlap)
}

sisterSpeciesRangeAsymmetry <- function(sp.rasters, sisters){
  
  sister.asymmetry <- c()
  
  for (j in 1:length(sisters)) {
    
    sis.1 <- sp.rasters[[which(names(sp.rasters) == sisters[[j]][1])]]
    
    sis.2 <- sp.rasters[[which(names(sp.rasters) == sisters[[j]][2])]]
    
    sis.1.area <- sum(values(sis.1), na.rm=T)
    
    sis.2.area <- sum(values(sis.2), na.rm=T)
    
    if(sis.1.area < sis.2.area) {
      
      sister.asymmetry[j] <-sis.1.area / (sis.1.area + sis.2.area)
      
    } else {
      
      sister.asymmetry[j] <-sis.2.area / (sis.2.area + sis.1.area)
      
    }
    
  }
  
  names(sister.asymmetry) <- sisters
  sister.asymmetry <- unlist(sister.asymmetry)
  
  return(sister.asymmetry)
  
}

sisterSpeciesRangeDistances <- function(sp.rasters, sisters, aggregate_size=2){
  
  range.distance <- c()
  
  range.distance.max <- c()
  
  for (j in 1:length(sisters)) {
    
    sis.1 <- sp.rasters[[which(names(sp.rasters) == sisters[[j]][1])]]
    
    sis.2 <- sp.rasters[[which(names(sp.rasters) == sisters[[j]][2])]]
    
    if(aggregate_size > 0){
      
      sis.1 <- aggregate(sis.1, aggregate_size)
      
      sis.2 <- aggregate(sis.2, aggregate_size)
      
    }
    
    p1 <- data.frame(crds(as.points(sis.1)))
    
    p2 <- data.frame(crds(as.points(sis.2)))

    p1 <- st_as_sf(p1, coords = c("x", "y"), crs = st_crs(4326))
    p2 <- st_as_sf(p2, coords = c("x", "y"), crs = st_crs(4326))
    
    #over the maximum distance between two points in the domain - this is so can standardise the distance as a function of the total possible distance for empirical stuff
    range.distance.test<- try(  min(st_distance(p1, p2)))
    
    while(class(range.distance.test) == "try-error") {
      
      sis.1 <- aggregate(sis.1, aggregate_size)
      
      sis.2 <- aggregate(sis.2, aggregate_size)
      
      p1 <- data.frame(crds(as.points(sis.1)))
  
      p2 <- data.frame(crds(as.points(sis.2)))
      
      p1 <- st_as_sf(p1, coords = c("x", "y"), crs = st_crs(4326))
      
      p2 <- st_as_sf(p2, coords = c("x", "y"), crs = st_crs(4326))
      
      #over the maximum distance between two points in the domain - this is so can standardise the distance as a function of the total possible distance for empirical stuff
      range.distance.test<- try(  min(st_distance(p1, p2)))
    }
    range.distance[j]<- range.distance.test
    range.distance.max[j]<-max(st_distance(p1, p2))
  }
  
  names(range.distance) <- sisters
  
  range.distance <- range.distance/max(range.distance.max, na.rm=T)
  
  range.distance <- unlist(range.distance)
  return(range.distance)
}

outgroupOverlap <- function(sp.rasters, sisters, sisters.overlap=NA){ 
  
  df_OG<-data.frame(model = NA, TO = NA)

  OG.rasters <- list()
  
  OG<-findOG(phy, sisters)
  
  OG.overlap <- list()
  
  outgroup.overlap <- list()
  
  TO <- c()
  
  for (j in 1:length(sisters)) {
    
    OG.overlap.sis2 <- c()
    
    OG.overlap.sis1 <- c()
    
    for (k in 1:length(OG[[j]])) {
      
      sis.1 <- sp.rasters[[which(names(sp.rasters) == sisters[[j]][1])]]
      
      sis.2 <- sp.rasters[[which(names(sp.rasters) == sisters[[j]][2])]]
      
      sis.1.area <- sum(values(sis.1), na.rm=T)
      
      sis.2.area <- sum(values(sis.2), na.rm=T)
      
      og <- sp.rasters[[which(names(sp.rasters) == OG[[j]][k])]]
      
      sis1.overlap<-sis.1+og
      
      sis2.overlap<-sis.2+og
      
      sis.overlap<-sis.1+sis.2
      
      if(any(values(sis1.overlap) == 2, na.rm=T)) {
        
        if (sum(values(sis.1), na.rm=T) < sum(values(og), na.rm=T)) {
          
          range.overlap <- length(which(values(sis1.overlap) == 2))/ sum(values(sis.1), na.rm=T)
        }
        
        if (sum(values(sis.1), na.rm=T) > sum(values(og), na.rm=T)) {
          
          range.overlap <- length(which(values(sis1.overlap) == 2))/ sum(values(og), na.rm=T)
        }
        
      } else { range.overlap <- 0 }
      
      OG.overlap.sis1[k] <- range.overlap
      
      if(any(values(sis2.overlap) == 2, na.rm=T)) {
        
        if (sum(values(sis.2), na.rm=T) <= sum(values(og), na.rm=T)) {
          
          range.overlap <- length(which(values(sis2.overlap) == 2))/ sum(values(sis.2), na.rm=T)
        }
        
        if (sum(values(sis.2), na.rm=T) > sum(values(og), na.rm=T)) {
          
          range.overlap <- length(which(values(sis2.overlap) == 2))/ sum(values(og), na.rm=T)
        }
      } else { range.overlap <- 0 }
      
      OG.overlap.sis2[k] <- range.overlap
      
      og.overlap.list<-list(OG.overlap.sis1, OG.overlap.sis2)
      
    }  
    
    OG.overlap[[j]]<-og.overlap.list
    
    outgroup.overlap[[j]]<-mean(mean(og.overlap.list[[1]], mean(og.overlap.list[[2]])))
    
    TO[j]<- sisters.overlap[[j]] - outgroup.overlap[[j]]
  }

  names(TO) <- sisters

  return(TO)
}

rangeSize <- function(sp.rasters){
  
  species.range.sizes <- c()
  
  for (j in 1:nlyr(sp.rasters)) {
    
    species.range.sizes[j] <- length(which(!is.na(values(sp.rasters[[j]]))))
    
  }
  return(species.range.sizes)
}

calculateARC <- function(phy, sp.rasters){
  
  ARC.clade <- list()
  
  for (i in 1:length(phy$tip.label)){
    enm.sp <- enmtools.species()
    enm.sp$species.name <- phy$tip.label[[i]]
    enm.sp$presence.points <- rasterToPoints(sp.rasters[[which(names(sp.rasters)==enm.sp$species.name)]])[,1:2]
    if (length(enm.sp$presence.points)>2){
      enm.sp$presence.points <- as.data.frame(enm.sp$presence.points)
      colnames(enm.sp$presence.points)[1] <- "Longitude"
      colnames(enm.sp$presence.points)[2] <- "Latitude"
    } else{
      enm.sp$presence.points <- data.frame(Longitude=enm.sp$presence.points[[1]],Latitude=enm.sp$presence.points[[2]])
    }
    enm.sp$range <- sp.rasters[[which(names(sp.rasters)==enm.sp$species.name)]]
    ARC.clade$species[[i]] <- enm.sp
  }
  names(ARC.clade$species) <- phy$tip.label
  ARC.clade$tree <- phy
  
  ARC <- try(enmtools.aoc(clade = ARC.clade,  nreps = 0, overlap.source = "range"))
  
 return(ARC)
  
}
