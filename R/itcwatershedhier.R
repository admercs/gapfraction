itc.watershed.hier <- function(chm.stack=NA, ht2rad=NA, tolerance=0.1, fun=max, res=1, num=TRUE, stacked=FALSE, silent=FALSE, ws.plot=FALSE) {

  #require(raster)
  #require(igraph)

  if(nlayers(chm.stack)==1) stop('Only a single raster layer found; Please input a raster stack')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  val <- seq(from=0, to=max(values(chm.stack)[!is.na(values(chm.stack))]), length.out=length(values(chm.stack)))
  bw  <- myColorRamp(colors=c('black','white'), values=val)

  chm.out   <- stackApply(chm.stack, indices=c(1), fun=max, na.rm=T)
  chm.unstk <- unstack(chm.stack)
  itc.layer <- list()

  for(i in 2:(length(chm.unstk))) {
    chm.in  <- chm.unstk[[i]]
    chm.wat <- itc.watershed(chm=chm.in, tolerance=tolerance, ht2rad=ht2rad, res=res, silent=silent, ws.plot=ws.plot, num=F)
    itc.layer[i] <- chm.wat
  }

  itc.layer  <- itc.layer[!sapply(itc.layer, is.null)]
  itc.stack  <- stack(itc.layer)
  itc.out    <- stackApply(itc.stack, indices=c(1), fun=fun, na.rm=T)
  itc.out[itc.out==0] <- NA
  itc.crowns <- ht2rad(chm.out[!is.na(itc.out)])
  crown.area <- sum(itc.crowns) / (length(values(chm.out)[!is.na(values(chm.out))]) * res^2) * 100
  itc.trees  <- clump(itc.out, directions=8)
  ntrees     <- length(unique(itc.trees))
  message(ntrees,' total trees counted at a tolerance of ',tolerance,' meters')

  if(stacked==FALSE) {
    if(num==FALSE) {
      plot(chm.out, col=gray.colors(256, start=0, end=1, gamma=1, alpha=NULL))
      plot(itc.out, add=TRUE, col='red', legend=FALSE)
      pts <- rasterToPoints(itc.out)
      points(pts, cex=chm.out[!is.na(itc.out)]/4, pch=10, lwd=1)
      return(itc.out)
    } else return(c( trees=ntrees, crown.area=crown.area ))
  } else if(stacked==TRUE) {
    plot(itc.stack, col=col)
    return(itc.stack)
  }
}
