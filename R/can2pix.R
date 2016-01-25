can2pix <- function(chm=NA, thresh.val=1.25, silent=FALSE) {

  canopy    <- raster::values(chm)[raster::values(chm) > thresh.val]
  all.cells <- length(raster::values(chm)[!is.na(raster::values(chm))])
  can.cells <- length(canopy[!is.na(canopy)])
  output    <- can.cells/all.cells

  col <- myColorRamp(colors=c('blue','green','yellow','red'), values=canopy)

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(chm, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy Ratio')
  }
  return(output)
}
