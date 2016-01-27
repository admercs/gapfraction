itc.varwin.hier <- function(chm.stack=NA, ht2rad=NA, type='circle', res=1, fun=max, num=TRUE, stacked=FALSE, silent=FALSE, plots=FALSE, geoTIFF=FALSE) {

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  chm.out   <- raster::stackApply(chm.stack, indices=c(1), fun=max, na.rm=T)
  chm.unstk <- raster::unstack(chm.stack)
  itc.layer <- list()

  for(i in 2:(length(chm.unstk))) {
    layer.in  <- chm.unstk[[i]]
    itc.out   <- itc.varwin(chm=layer.in, ht2rad=ht2rad, type=type, res=res, silent=silent, plots=plots, geoTIFF=geoTIFF, num=FALSE)
    layer.ext <- raster::extend(itc.out, c(1,1))
    itc.clump <- raster::clump(layer.ext, directions=8)
    f <- freq(itc.clump)
    f <- as.data.frame(f)
    exclude <- f$value[which(f$count > 1)]
    raster::values(itc.clump)[raster::values(itc.clump) %in% exclude] <- NA
    raster::values(itc.clump)[!raster::values(itc.clump) %in% exclude] <- 1
    itc.layer[i] <- itc.clump
  }

  itc.layrr <- itc.layer[!sapply(itc.layer, is.null)]
  itc.stack <- raster::stack(itc.layrr)
  itc.out   <- raster::stackApply(itc.stack, indices=c(1), fun=fun, na.rm=T)
  raster::values(itc.out)[raster::values(itc.out)==0] <- NA
  itc.trees <- raster::clump(itc.out, directions=8)
  ntrees    <- length(unique(raster::values(itc.trees)))
  message('There are an estimated ',ntrees,' trees in the plot')

  itc.crowns <- ht2rad(raster::values(chm.out)[!is.na(raster::values(itc.trees))])
  crown.area <- sum(itc.crowns[!is.na(itc.crowns)]) / (length(raster::values(chm.out)[!is.na(raster::values(chm.out))]) * res^2) * 100

  if(num==TRUE) {
    return( c(trees=ntrees, crown.area=crown.area) )
  } else if(stacked==FALSE) {
    return(itc.out)
  } else if(stacked==TRUE) {
    return(itc.stack)
  }
}
