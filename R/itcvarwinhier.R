itc.varwin.hier <- function(chm.stack=NA, ht2rad=NA, min.h=1, type='circle', res=1, fun=max, num=TRUE, stacked=FALSE, silent=FALSE, plots=FALSE, geoTIFF=FALSE) {

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  chm.out <- raster::stackApply(chm.stack, indices=c(1), fun=max, na.rm=T)
  if(max(raster::values(chm.out)[!is.na(raster::values(chm.out))]) <= min.h) return(c(trees=NA, crown.area=NA))
  chm.unstk <- raster::unstack(chm.stack)
  itc.layer <- list()

  for(i in 2:(length(chm.unstk))) {
    chm <- chm.unstk[[i]]
    if(max(raster::values(chm)[!is.na(raster::values(chm))]) <= min.h) next
    itc.layer[i] <- itc.varwin(chm=chm, ht2rad=ht2rad, min.h=min.h, type=type, res=res, silent=silent, plots=plots, geoTIFF=geoTIFF, num=FALSE)
  }

  itc.layrr <- itc.layer[!sapply(itc.layer, is.null)]
  itc.stack <- raster::stack(itc.layrr)
  itc.fun   <- raster::stackApply(itc.stack, indices=c(1), fun=fun, na.rm=T)
  itc.out   <- raster::stackApply(itc.stack, indices=c(1), fun=max, na.rm=T)

  itc.trees <- raster::clump(itc.out, directions=8)
  ntrees    <- length(unique(raster::values(itc.trees)))
  message('There are an estimated ',ntrees,' trees in the plot')

  itc.crowns <- ht2rad(chm.out) * itc.out
  crown.area <- sum(itc.crowns[!is.na(itc.crowns)]) / (length(raster::values(chm.out)[!is.na(raster::values(chm.out))]) * res)

  if(num==TRUE) {
    return( c(trees=ntrees, crown.area=crown.area) )
  } else if(stacked==FALSE) {
    return(itc.fun)
  } else if(stacked==TRUE) {
    return(itc.stack)
  }
}
