bias.sample <- function(bias.raster, npoints, replace = TRUE, biased = TRUE){
  
  # Convert the raster to a set of points
  bias.points <- data.frame(rasterToPoints(bias.raster))
  
  if(biased == TRUE){
    # Sampling is biased
    bias.points <- bias.points[sample(nrow(bias.points), 
                                      size = npoints, 
                                      prob = bias.points[,3],
                                      replace = replace),]
  } else {
    # Sampling is unbiased
    bias.points <- bias.points[sample(nrow(bias.points), 
                                      size = npoints, 
                                      replace = replace),]
  }
  
  return(bias.points[,1:2])
}

# bias.sample(biasBlindsnakes, npoints = 1000)
