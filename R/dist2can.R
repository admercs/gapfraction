dist2can <- function(chm=NA, thresh.val=1.25, degrees=FALSE, from=FALSE) {

  canopy <- chm
  raster::values(canopy)[raster::values(canopy) <  thresh.val] <- NA
  raster::values(canopy)[raster::values(canopy) >= thresh.val] <- 1

  rws <- nrow(chm)
  cls <- ncol(chm)
  center <- raster::cellFromRowCol(object=canopy, rownr=rws/2, colnr=cls/2)

  dist2c <- raster::distance(x=canopy)
  dist2c <- dist2c[center]

  dir2c <- raster::direction(x=canopy, degrees=degrees, from=from, doEdge=F)
  dir2c <- dir2c[center]

  return( c(dist2can=dist2c, dir2can=dir2c) )
}
