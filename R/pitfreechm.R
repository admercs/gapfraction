pitfreechm <- function(las.path=NA, las.proj=NA, las.reproj=NA, breaks=c(2,5,10,15), nx=100, ny=100, ko=2.5, ku=20, stacked=FALSE, plots=FALSE, geoTIFF=FALSE) {

  if (is.na(LASpath)) stop('Please input a full file path to the LAS file')

  #require(abind)
  #require(magic)
  #require(geometry)
  #require(raster)
  #require(sp)
  #require(spatstat)
  #require(rLiDAR)
  #require(Matrix)
  #require(rasterVis)

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  bary.2d <- function(X=NA, f=NA, Xi=NA, k=NA) {
    dn  <- geometry::delaunayn(X)
    tri <- geometry::tsearch(x=X[,1], y=X[,2], t=dn, xi=Xi[,1], yi=Xi[,2], bary=T)

    active <- dn[tri$idx,]

    cart.pt1 <- matrix(X[active[,1],], ncol=2); colnames(cart.pt1) <- c('x','y')
    cart.pt2 <- matrix(X[active[,2],], ncol=2); colnames(cart.pt2) <- c('x','y')
    cart.pt3 <- matrix(X[active[,3],], ncol=2); colnames(cart.pt3) <- c('x','y')

    tri.dist1 <- raster::pointDistance(p1=cart.pt1, p2=cart.pt2, lonlat=F, allpairs=F)
    tri.dist2 <- raster::pointDistance(p1=cart.pt2, p2=cart.pt3, lonlat=F, allpairs=F)
    tri.dist3 <- raster::pointDistance(p1=cart.pt1, p2=cart.pt3, lonlat=F, allpairs=F)

    active[tri.dist1 > k,] <- NA
    active[tri.dist2 > k,] <- NA
    active[tri.dist3 > k,] <- NA

    active.filt <- active[complete.cases(active[,c(1,2,3)]),]
    Xi.filt <- Xi[complete.cases(active[,c(1,2,3)]),]
    tri <- geometry::tsearch(x=X[,1], y=X[,2], t=dn, xi=Xi.filt[,1], yi=Xi.filt[,2], bary=T)
    M <- Matrix::sparseMatrix(i=rep(1:nrow(Xi.filt),each=3), j=as.numeric(t(active.filt)), x=as.numeric(t(tri$p)), dims=c(nrow(Xi.filt), length(f)))

    result <- as.numeric(M %*% f)
    output <- matrix(c(Xi.filt, result), ncol=3)
    colnames(output) <- c('x','y','z')
    return(output)
  }

  chm <- function(las=NA, nx=nx, ny=ny, w=chull.all) {
    centers <- spatstat::gridcentres(w, nx=nx, ny=ny)
    in.win  <- spatstat::inside.owin(x=centers$x, y=centers$y, w=w)
    xo      <- centers$x[in.win]
    yo      <- centers$y[in.win]
    grd.2d  <- matrix(c(xo,yo), nrow=length(xo), ncol=2)
    las.ext <- raster::extent(grd.2d)
    #las.ext <- extent(las[,1:2])
    las.ras <- raster::raster(las.ext, ncols=nx, nrows=ny)
    las.chm <- raster::rasterize(las[,1:2], las.ras, las[,3], fun=max)
    return(las.chm)
  }

  tin <- function(las=NA, nx=nx, ny=ny, k=NA, w=chull.all) {
    centers <- spatstat::gridcentres(w, nx=nx, ny=ny)
    in.win  <- spatstat::inside.owin(x=centers$x, y=centers$y, w=w)
    xo      <- centers$x[in.win]
    yo      <- centers$y[in.win]
    grd.2d  <- matrix(c(xo,yo), nrow=length(xo), ncol=2)
    tin.ext <- raster::extent(grd.2d)
    tin.ras <- raster::raster(tin.ext, ncols=nx, nrows=ny, crs=CRS(r.crs))
    tin.las <- bary.2d(X=las[,1:2], f=las[,3], Xi=grd.2d, k=k)
    chm.tin <- raster::rasterize(tin.las[,1:2], tin.ras, tin.las[,3], fun=max)
    return(chm.tin)
  }

  # Read LAS data
  LAS       <- rLiDAR::readLAS(las.path, short=FALSE)
  LAS       <- LAS[order(LAS[,3], decreasing=FALSE), ]
  LASfolder <- dirname(las.path)
  LASname   <- strsplit(basename(las.path), '\\.')[[1]][1]

  # Reproject LAS data if specified
  if(!is.na(las.reproj)) {
    CRSin  <- las.proj
    CRSout <- las.reproj
    spLAS  <- sp::SpatialPoints(LAS[,c(1:3)], proj4string=CRS(CRSin))
    tmLAS  <- sp::spTransform(spLAS, CRSobj=CRS(CRSout))
    rpLAS  <- sp::coordinates(tmLAS)
    LAS    <- cbind(rpLAS, LAS[,c(4:12)])
  }

  if(!is.na(las.reproj)) {
    r.crs <- las.reproj
  } else r.crs <- las.proj

  # Create color ramp for plots
  val <- seq(from=0, to=max(LAS[,3]), length.out=length(LAS[,3]))
  col <- myColorRamp(colors=c('blue','green','yellow','red'), values=val)

  # Compute convex hull to use as a window
  chull.all <- spatstat::convexhull.xy(x=LAS[,1], y=LAS[,2])

  # Barycentrically interpolated (effective TIN) canopy height models
  ground  <- tin(las=LAS[LAS[,3] == 0,], nx=nx, ny=ny, k=ku, w=chull.all)
  tin.all <- tin(las=LAS[LAS[,5] == 1,], nx=nx, ny=ny, k=ko)
  tin.all <- raster::stackApply(raster::stack(tin.all,ground), indices=c(1), fun=max, na.rm=T)
  tins    <- list()

  for(i in 1:length(breaks)) {
    tin.break <- tin(las=LAS[LAS[,5] == 1 & LAS[,3] >=  breaks[i],], nx=nx, ny=ny, k=ko, w=chull.all)
    tin.break <- raster::stack(tin.break, ground)
    tins[i]   <- raster::stackApply(tin.break, indices=c(1), fun=max, na.rm=T)
  }

  # Pit-free CHM
  tins <- raster::stack(tins)
  chms <- raster::stack(ground, tins, tin.all)
  names(chms) <- c('Ground Returns',breaks,'All First Returns')
  pitfree <- raster::stackApply(chms, indices=c(1), fun=max, na.rm=T)

  if(plots==TRUE) {
    jpeg(file.path(LASfolder, paste(LASname,'_chm_tin.jpg',sep='')), width=12, height=8, units='in', res=300, quality=100)
    par(mfrow=c(2,3), pty='s', xpd=TRUE)
    raster::plot(ground,  col=col, box=F)
    raster::plot(tin.all, col=col, box=F)
    raster::plot(tins, col=col, box=F)
    dev.off()

    jpeg(file.path(LASfolder, paste(LASname,'_chm_pitfree.jpg',sep='')), width=16, height=16, units='in', res=300, quality=100)
    par(mfrow=c(2,2), pty='s', xpd=TRUE)
    plot(LAS, pch=19, cex=1, col=col)
    raster::plot(chm.all, col=col, box=F)
    raster::plot(tin.all, col=col, box=F)
    raster::plot(pitfree, col=col, box=F)
    dev.off()

    jpeg(file.path(LASfolder, paste(LASname,'_chm_all_vis.jpg',sep='')), width=8, height=8, units='in', res=300, quality=100)
    par(mfrow=c(1,1), pty='s', xpd=TRUE)
    rasterVis::levelplot(chm.all, contour=TRUE, col.regions=col)
    dev.off()

    jpeg(file.path(LASfolder, paste(LASname,'_chm_tin_vis.jpg',sep='')), width=8, height=8, units='in', res=300, quality=100)
    par(mfrow=c(1,1), pty='s', xpd=TRUE)
    rasterVis::levelplot(tin.all, contour=TRUE, col.regions=col)
    dev.off()

    jpeg(file.path(LASfolder, paste(LASname,'_chm_pitfree_vis.jpg',sep='')), width=8, height=8, units='in', res=300, quality=100)
    par(mfrow=c(1,1), pty='s', xpd=TRUE)
    rasterVis::levelplot(pitfree, contour=TRUE, col.regions=col)
    dev.off()
  }
  if(stacked==FALSE) {
    raster::plot(pitfree, col=col)
    if(geoTIFF==TRUE) {
      fname <- paste(LASname,'_pfchm.tiff',sep='')
      raster::writeRaster(x=pitfree, filename=file.path(LASfolder,fname), format='GTiff', overwrite=T)
    }
    return(pitfree)
  } else if(stacked==TRUE) {
    raster::plot(chms,  col=col)
    if(geoTIFF==TRUE) {
      fname <- paste(LASname,'_pfchm_stack.tiff',sep='')
      raster::writeRaster(x=chms, filename=file.path(LASfolder,fname), format='GTiff', overwrite=T)
    }
    return(chms)
  }
}
