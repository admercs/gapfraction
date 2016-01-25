can2pix <- function(chm=NA, thresh.val=1.25, silent=FALSE) {

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  canopy    <- raster::values(chm)[raster::values(chm) > thresh.val]
  all.cells <- length(raster::values(chm)[!is.na(raster::values(chm))])
  can.cells <- length(canopy[!is.na(canopy)])
  output    <- can.cells/all.cells

  col <- myColorRamp(colors=c('blue','green','yellow','red'), values=raster::values(chm)[!is.na(raster::values(chm))])

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(chm, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy Ratio')
  }
  return(output)
}
