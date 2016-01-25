can2pix <- function(chm=NA, thresh.val=1.25, silent=FALSE) {

  #require(raster)

  canopy    <- raster::values(chm)[raster::values(chm) > thresh.val]
  all.cells <- length(raster::values(chm)[!is.na(raster::values(chm))])
  can.cells <- length(raster::values(canopy)[!is.na(raster::values(canopy))])
  output    <- can.cells/all.cells

  col <- myColorRamp(colors=c('blue','green','yellow','red'), values=raster::values(canopy))

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    raster::plot(canopy, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy Ratio')
  }
  return(output)
}
