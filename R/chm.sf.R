#' Barycentric interpolation-based Spike-free Canopy Height Model
#'
#' This function implements a 2-D barycentric interpolation-based spike-free canopy height model algorithm
#' @param las Path or name of LAS file. Defaults to NA.
#' @param las.proj Proj4 projection string to use for projection. Defaults to NA.
#' @param las.reproj Proj4 projection string to use for reprojection. Defaults to NA.
#' @param breaks List of breaks to use for vertical stratification, either values or percentages of maximum height. Defaults to percentages: 0.10, 0.25, 0.50, 0.75.
#' @param percent Boolean switch for the use of percentages. Defaults to TRUE.
#' @param nx Number of pixels along the x-axis. For a 50m radius plot, nx=100 is 1m resolution. Defaults to 100.
#' @param ny Number of pixels along the y-axis. For a 50m radius plot, ny=100 is 1m resolution. Defaults to 100.
#' @param ko Barycentric triangle edge length filter for overstory layers. Defaults to 2.5.
#' @param ku Barycentric triangle edge length filter for the understory layer. Defaults to 20.
#' @param stacked Boolean switch for the output of a single raster or a raser stack. Defaults to FALSE.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @param plots Boolean switch for the saving of plot files to the las.path folder. Defaults to FALSE.
#' @param geoTIFF Boolean switch for the saving of projected GeoTIFF files to the las.path folder. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://www.sciencedirect.com/science/article/pii/S0303243416300873}
#' @references \url{http://www.riegl.com/uploads/tx_pxpriegldownloads/khosravipour_SilviLaser2013.pdf}
#' @keywords chm, canopy height model, spike-free chm, barycentric interpolation
#' @export
#' @return The results of \code{chm.sf}
#' @examples
#' chm.sf(las='C:/plot.las', las.proj='+init=epsg:26911', las.reproj=NA, breaks=c(0.10,0.25,0.50,0.75), percent=TRUE, nx=100, ny=100, ko=2.5, ku=20, stacked=FALSE, silent=FALSE, plots=FALSE, geoTIFF=FALSE)

chm.sf <- function(las=NA, las.proj=NA, las.reproj=NA, breaks=c(0.10,0.25,0.50,0.75), percent=TRUE, nx=100, ny=100, ko=2.5, ku=20, stacked=FALSE, silent=FALSE, plots=FALSE, geoTIFF=FALSE) {

  if(length(las)==1 & any(is.na(eval(las)))) stop('Please input full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  bary.2d <- function(X=NA, f=NA, Xi=NA, k=NA) {

    dn  <- geometry::delaunayn(X, options='QbB')
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
    Xi.filt     <- Xi[complete.cases(active[,c(1,2,3)]),]
    dims        <- c(nrow(Xi.filt), ncol(Xi.filt), length(f), length(X[,1]), length(X[,2]))

    if(is.null(dim(Xi.filt))) {
      return(NULL)
    } else if(any(dims < 2)) {
      return(NULL)
    } else {
      tri <- geometry::tsearch(x=X[,1], y=X[,2], t=dn, xi=Xi.filt[,1], yi=Xi.filt[,2], bary=T)
      M   <- Matrix::sparseMatrix(i=rep(1:nrow(Xi.filt),each=3), j=as.numeric(t(active.filt)),
                                     x=as.numeric(t(tri$p)), dims=c(nrow(Xi.filt), length(f)))
      result <- as.numeric(M %*% f)
      output <- matrix(c(Xi.filt, result), ncol=3)
      colnames(output) <- c('x','y','z')
      return(output)
    }
  }

  tin <- function(las=NA, nx=nx, ny=ny, k=NA, w=chull.all) {
    centers <- spatstat::gridcentres(w, nx=nx, ny=ny)
    in.win  <- spatstat::inside.owin(x=centers$x, y=centers$y, w=w)
    xo      <- centers$x[in.win]
    yo      <- centers$y[in.win]
    grd.2d  <- matrix(c(xo,yo), nrow=length(xo), ncol=2)
    tin.ext <- raster::extent(grd.2d)
    tin.ras <- raster::raster(tin.ext, ncols=nx, nrows=ny, crs=r.crs)
    tin.las <- bary.2d(X=las[,1:2], f=las[,3], Xi=grd.2d, k=k)
    if(is.null(tin.las)) {
      chm.tin <- NULL
    } else chm.tin <- raster::rasterize(tin.las[,1:2], tin.ras, tin.las[,3], fun=max)
    return(chm.tin)
  }

  if(!exists("las")) {
    LAS       <- lidR::readLAS(las)
    LASfolder <- dirname(las)
    LASname   <- strsplit(basename(las),'\\.')[[1]][1]
  } else LAS  <- las
  LAS  <- LAS[order(LAS[,3], decreasing=FALSE),]

  val <- seq(from=0, to=max(LAS[,3]), length.out=length(LAS[,3]))
  col <- myColorRamp(colors=c('blue','green','yellow','red'), values=val)

  if(percent==TRUE) {
    zmax <- max(LAS[,3])
    for(i in 1:length(breaks)) {
      breaks[i] <- zmax * breaks[i]
    }
  }

  if(!is.na(las.reproj)) {
    CRSin  <- las.proj
    CRSout <- las.reproj
    spLAS  <- sp::SpatialPoints(LAS[,c(1:3)], proj4string=sp::CRS(CRSin))
    tmLAS  <- sp::spTransform(spLAS, CRSobj=sp::CRS(CRSout))
    rpLAS  <- sp::coordinates(tmLAS)
    LAS    <- cbind(rpLAS, LAS[,c(4:12)])
  }

  if(!is.na(las.reproj)) {
    r.crs <- las.reproj
  } else r.crs <- las.proj

  chull.all <- spatstat::convexhull.xy(x=LAS[,1], y=LAS[,2])
  ground    <- tin(las=LAS[LAS[,3] == 0,], nx=nx, ny=ny, k=ku, w=chull.all)
  tin.all   <- tin(las=LAS[LAS[,5] == 1,], nx=nx, ny=ny, k=ko)
  tin.all   <- raster::stackApply(raster::stack(tin.all, ground), indices=c(1), fun=max, na.rm=T)
  tins      <- list()
  pos       <- c()

  for(i in 1:length(breaks)) {
    las.lay <- LAS[LAS[,5] == 1 & LAS[,3] >=  breaks[i],]
    if(is.null(dim(las.lay))) next
    dupl    <- !duplicated.matrix(las.lay[,1:2])
    las     <- las.lay[dupl,]
    if(dim(las)[1] > 1) {
      tin.break <- tin(las=las, nx=nx, ny=ny, k=ko, w=chull.all)
    }
    if(is.null(tin.break)) {
      next
    } else {
      tin.break <- raster::stack(tin.break, ground)
      tins[[i]] <- raster::stackApply(tin.break, indices=c(1), fun=max, na.rm=T)
    }
    pos <- c(pos,i)
  }

  tins <- tins[!sapply(tins, is.null)]
  if(length(tins)==0 | is.null(tins)) {
    if(silent==FALSE) plot(tin.all, col=col)
    return(tin.all)
  }
  tins <- raster::stack(tins)
  chms <- raster::stack(ground, tins, tin.all)
  pitfree <- raster::stackApply(chms, indices=c(1), fun=max, na.rm=T)

  if(plots==TRUE) {
    jpeg(file.path(LASfolder, paste(LASname,'_chm_tin.jpg',sep='')), width=12, height=8, units='in', res=300, quality=100)
    par(mfrow=c(2,3), pty='s', xpd=TRUE)
    plot(ground,  col=col, box=F)
    plot(tin.all, col=col, box=F)
    plot(tins, col=col, box=F)
    dev.off()

    jpeg(file.path(LASfolder, paste(LASname,'_chm_pitfree.jpg',sep='')), width=16, height=16, units='in', res=300, quality=100)
    par(mfrow=c(2,2), pty='s', xpd=TRUE)
    plot(LAS, pch=19, cex=1, col=col)
    plot(chm.all, col=col, box=F)
    plot(tin.all, col=col, box=F)
    plot(pitfree, col=col, box=F)
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
    if(silent==FALSE) plot(pitfree, col=col)
    if(geoTIFF==TRUE) {
      fname <- paste(LASname,'_pfchm.tiff',sep='')
      raster::writeRaster(x=pitfree, filename=file.path(LASfolder,fname), format='GTiff', overwrite=T)
    }
    return(pitfree)
  } else if(stacked==TRUE) {
    if(silent==FALSE) plot(chms,  col=col)
    if(geoTIFF==TRUE) {
      fname <- paste(LASname,'_pfchm_stack.tiff',sep='')
      raster::writeRaster(x=chms, filename=file.path(LASfolder,fname), format='GTiff', overwrite=T)
    }
    return(chms)
  }
}
