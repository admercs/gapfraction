rawchm <- function(las.path=NA, las.proj=NA, las.reproj=NA, breaks=c(2,5,10,15), nx=100, ny=100, stacked=FALSE, plots=FALSE, geoTIFF=FALSE) {
  
  require(raster)
  require(spatstat)
  require(rLiDAR)
  require(rasterVis)
  
  # Raw CHM function
  chm <- function(las=NA, nx=nx, ny=ny, w=chull.all) {
    centers <- gridcentres(w, nx=nx, ny=ny)
    in.win  <- inside.owin(x=centers$x, y=centers$y, w=w)
    xo      <- centers$x[in.win]
    yo      <- centers$y[in.win]
    grd.2d  <- matrix(c(xo,yo), nrow=length(xo), ncol=2)
    las.ext <- extent(grd.2d)
    #las.ext <- extent(las[,1:2])
    las.ras <- raster(las.ext, ncols=nx, nrows=ny)
    las.chm <- rasterize(las[,1:2], las.ras, las[,3], fun=max)
    return(las.chm)
  }
  
  # Color ramp function
  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }
  
  # Read LAS data
  LAS       <- readLAS(las.path, short=FALSE)
  LAS       <- LAS[order(LAS[,3], decreasing=FALSE), ]
  LASfolder <- dirname(las.path)
  LASname   <- strsplit(basename(las.path), '\\.')[[1]][1]
  
  # Reproject LAS data if specified
  if(!is.na(las.reproj)) {
    CRSin  <- las.proj
    CRSout <- las.reproj
    spLAS  <- SpatialPoints(LAS[,c(1:3)], proj4string=CRS(CRSin))
    tmLAS  <- spTransform(spLAS, CRSobj=CRS(CRSout))
    rpLAS  <- coordinates(tmLAS)
    LAS    <- cbind(rpLAS, LAS[,c(4:12)])
  }
  
  if(!is.na(las.reproj)) {
    r.crs <- las.reproj
  } else r.crs <- las.proj
  
  val <- seq(from=0, to=max(LAS[,3]), length.out=length(LAS[,3]))
  col <- myColorRamp(colors=c('blue','green','yellow','red'), values=val)
  
  # Compute convex hull of grid to use as a window
  chull.all <- convexhull.xy(x=LAS[,1], y=LAS[,2])
  
  # Raw canopy height models
  ground   <- chm(las=LAS[LAS[,3] == 0,], nx=nx, ny=ny, w=chull.all)
  chm.all  <- chm(las=LAS[LAS[,5] == 1,], nx=nx, ny=ny, w=chull.all)
  chm.brks <- list()
  
  for(i in 1:length(breaks)) {
    chm.brks[i] <- chm(las=LAS[LAS[,5] == 1 & LAS[,3] >=  breaks[i],], nx=nx, ny=ny, w=chull.all)
  }
  
  chm.brks <- stack(chm.brks)
  chms <- stack(ground, chm.brks, chm.all)
  names(chms) <- c('Ground Returns',breaks,'All First Returns')
  chm  <- stackApply(chms, indices=c(1), fun=max, na.rm=T)
  
  if(plots==TRUE) {
    par(mfrow=c(1,1), pty='s', xpd=TRUE)
    jpeg(file.path(LASfolder, paste(LASname,'_chm_raw_breeaks.jpg',sep='')), width=12, height=8, units='in', res=300, quality=100)
    plot(chm.brks, col=col, box=F)
    dev.off()
    
    jpeg(file.path(LASfolder, paste(LASname,'_chm_raw.jpg',sep='')), width=12, height=8, units='in', res=300, quality=100)
    plot(chm, col=col, box=F)
    dev.off()
    
    jpeg(file.path(LASfolder, paste(LASname,'_chm_raw_level.jpg',sep='')), width=12, height=8, units='in', res=300, quality=100)
    levelplot(chm, contour=T, col.regions=col)
    dev.off()
  }
  if(stacked==FALSE) {
    if(geoTIFF==TRUE) {
      fname <- paste(LASname,'_rawchm.tiff',sep='')
      writeRaster(x=chm, filename=file.path(LASfolder,fname), format='GTiff', overwrite=T)
    }
    plot(chm, col=col)
    return(chm)
  } else if(stacked==TRUE) {
    if(geoTIFF==TRUE) {
      fname <- paste(LASname,'_rawchm_stack.tiff',sep='')
      writeRaster(x=chms, filename=file.path(LASfolder,fname), format='GTiff', overwrite=T)
    }
    plot(chms, col=col)
    return(chms) 
  }
}
