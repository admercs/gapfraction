dist2can <- function(chm=NA, thresh.val=1.25, degrees=FALSE, from=FALSE) {
  
  require(raster)
  
  canopy <- chm
  values(canopy)[values(canopy) <  thresh.val] <- NA
  values(canopy)[values(canopy) >= thresh.val] <- 1
  
  rws <- nrow(chm)
  cls <- ncol(chm)
  center <- cellFromRowCol(object=canopy, rownr=rws/2, colnr=cls/2)
  
  # Calculate nearest distance to canopy
  dist2c <- distance(x=canopy)
  dist2c <- dist2c[center]
  
  # Calculate direction to canopy
  dir2c <- direction(x=canopy, degrees=degrees, from=from, doEdge=F)
  dir2c <- dir2c[center]
  
  return( c(dist2can=dist2c, dir2can=dir2c) )
}
