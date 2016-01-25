itc.varwin <- function(chm=NA, ht2rad=NA, type='circle', res=1, num=TRUE, plots=FALSE, geoTIFF=FALSE) {

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  run.focal <- function(chm=chm, hts=hts, rd2=rd2, rad, type=type) {
    htt <- hts[rd2==rad | rd2==rad-1]
    x2  <- chm > min(htt) & chm < max(htt)
    x3  <- x2 * chm
    fun <- function(z) ifelse(z[length(z)/2 + 0.5]==max(z), 1, NA)
    wts <- raster::focalWeight(x=x3, d=rad, type=type)
    res <- raster::focal(x=chm, w=wts, fun=fun, na.rm=F, pad=rad, padValue=NA, NAonly=F)
    return(res)
  }

  if(nlayers(chm) > 1) {
    chm <- raster::stackApply(chm, indices=c(1), fun=max, na.rm=T)
  }

  if(is.character(chm)) {
    isPath <- TRUE
    folder <- dirname(chm)
    name   <- strsplit(basename(chm), '\\.')[[1]][1]
    chm    <- raster::raster(chm)
  } else isPath <- FALSE

  hts <- sort(unique(round(chm[chm >= 2 & !is.na(chm)])))
  rd1 <- ht2rad(hts)
  rd2 <- round(rd1)
  rd3 <- sort(unique(rd2))
  for(i in 1:length(rd3)) rd3[i] <- ifelse(rd3[i] %% 2 != 0, rd3[i], rd3[i] + 1)
  if(length(rd3)==0) stop('Trees of no suitable height classes exist for the crown area moving window')
  message('Computing ', length(rd3), ' moving window(s)')

  if(length(rd3==1)) {
    itc.out <- run.focal(chm=chm, hts=hts, rd2=rd2, rad=rd3, type=type)
  } else {
    itc.stk <- raster::stack()
    for(i in 1:length(rd3)) {
      itc.new  <- run.focal(chm=chm, hts=hts, rd2=rd2, rad=rd3[i], type=type)
      itc.stk  <- raster::stack(itc.stk, itc.new)
    }
    itc.out <- raster::stackApply(itc.stk, indices=c(1), fun=max, na.rm=T)
  }

  itc.trees <- raster::clump(itc.out, directions=8)
  ntrees <- length(unique(raster::values(itc.trees)[!is.na(raster::values(itc.trees))]))
  message('There are an estimated ',ntrees,' trees in the plot')
  itc.crowns <- ht2rad(chm) * itc.out
  crown.area <- sum(raster::values(itc.crowns)[!is.na(raster::values(itc.crowns))]) / (length(raster::values(chm)[!is.na(raster::values(chm))]) * res^2) * 100

  val <- seq(from=0, to=max(raster::values(chm)[!is.na(raster::values(chm))]), length.out=length(raster::values(chm)))
  bw  <- myColorRamp(cols=c('black','white'), vals=val)

  if(plots==TRUE) {
    jpeg(file.path(LASfolder, paste(LASname,'_trees.jpg',sep='')), width=8, height=8, units='in', res=300, quality=100)
    par(mfrow=c(1,1), pty='s', xpd=TRUE)
    plot(chm, col=bw)
    plot(itc.out, col='red', add=T, legend=F)
    pts <- raster::rasterToPoints(itc.out)
    points(pts, cex=x[!is.na(itc.out)]/5, pch=10, lwd=1)
    dev.off()
  }
  if(isPath==TRUE & geoTIFF==TRUE) {
    fname <- paste(name,'_itd.tiff',sep='')
    raster::writeRaster(x=itc.out, filename=file.path(folder,fname), format='GTiff', overwrite=T)
  } else if(isPath==FALSE & geoTIFF==TRUE) {
    raster::writeRaster(x=itc.out, filename='itd.tiff', format='GTiff', overwrite=T)
  }
  plot(chm, col=bw)
  plot(itc.out, col='red', add=T, legend=F)
  pts <- raster::rasterToPoints(itc.out)
  points(pts, cex=chm[!is.na(itc.out)]/4, pch=10, lwd=1)

  if(num==TRUE) {
    return( c(trees=ntrees, crown.area=crown.area))
  } else return(itc.out)
}
