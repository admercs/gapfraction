#' Simple Canopy Height Model
#'
#' This function implements a simple grid-based canopy height model
#' @param las Path or name of LAS file. Defaults to NA.
#' @param las.proj Proj4 projection string to use for projection. Defaults to NA.
#' @param las.reproj Proj4 projection string to use for reprojection. Defaults to NA.
#' @param nx Number of pixels along the x-axis. For a 50m radius plot, nx=100 is 1m resolution. Defaults to 100.
#' @param ny Number of pixels along the y-axis. For a 50m radius plot, ny=100 is 1m resolution. Defaults to 100.
#' @param fun Function for the calculate of height values in each cell. Defaults to max.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @param plots Boolean switch for the saving of plot files to the las.path folder. Defaults to FALSE.
#' @param geoTIFF Boolean switch for the saving of projected GeoTIFF files to the las.path folder. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://link.springer.com/book/10.1007/978-94-017-8663-8}
#' @keywords chm, canopy height model
#' @export
#' @return The results of \code{chm}
#' @examples
#' chm(las='C:/plot.las', las.proj='+init=epsg:26911', las.reproj=NA, nx=100, ny=100, fun=max, silent=FALSE, plots=FALSE, geoTIFF=FALSE)
#' chm(las='C:/plot.las', las.proj='+init=epsg:26911', las.reproj=NA, nx=100, ny=100, fun=function(x) quantile(x, 0.95), silent=FALSE, plots=FALSE, geoTIFF=FALSE)

chm <- function(las=NA, las.proj=NA, las.reproj=NA, nx=100, ny=100, fun=max, silent=FALSE, plots=FALSE, geoTIFF=FALSE) {

  if(length(las)==1 & any(is.na(eval(las)))) stop('Please input a full file path to the LAS file')

  chm.grid <- function(las=NA, nx=nx, ny=ny, w=chull, fun=fun) {
    if(is.null(dim(las))) return(NULL)
    centers <- spatstat::gridcentres(w, nx=nx, ny=ny)
    in.win  <- spatstat::inside.owin(x=centers$x, y=centers$y, w=w)
    xo      <- centers$x[in.win]
    yo      <- centers$y[in.win]
    grd.2d  <- matrix(c(xo,yo), nrow=length(xo), ncol=2)
    las.ext <- raster::extent(grd.2d)
    las.ras <- raster::raster(las.ext, ncols=nx, nrows=ny)
    las.chm <- raster::rasterize(las[,1:2], las.ras, las[,3], fun=fun)
    return(las.chm)
  }

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  if(!exists("las")) {
    LAS       <- lidR::readLAS(las)
    LASfolder <- dirname(las)
    LASname   <- strsplit(basename(las),'\\.')[[1]][1]
  } else LAS  <- las
  LAS  <- LAS[order(LAS[,3], decreasing=FALSE),]

  val <- seq(from=0, to=max(LAS[,3]), length.out=length(LAS[,3]))
  col <- myColorRamp(colors=c('blue','green','yellow','red'), values=val)

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

  chull <- spatstat::convexhull.xy(x=LAS[,1], y=LAS[,2])
  chm   <- chm.grid(las=LAS, nx=nx, ny=ny, w=chull, fun=fun)

  if(plots==TRUE) {
    par(mfrow=c(1,1), pty='s', xpd=TRUE)

    jpeg(file.path(LASfolder, paste(LASname,'_chm_raw_breeaks.jpg',sep='')), width=12, height=8, units='in', res=300, quality=100)
    plot(chm.brks, col=col, box=F)
    dev.off()

    jpeg(file.path(LASfolder, paste(LASname,'_chm_raw.jpg',sep='')), width=12, height=8, units='in', res=300, quality=100)
    plot(chm, col=col, box=F)
    dev.off()

    jpeg(file.path(LASfolder, paste(LASname,'_chm_raw_level.jpg',sep='')), width=12, height=8, units='in', res=300, quality=100)
    rasterVis::levelplot(chm, contour=T, col.regions=col)
    dev.off()
  }
  if(geoTIFF==TRUE) {
    fname <- paste(LASname,'_rawchm.tiff',sep='')
    raster::writeRaster(x=chm, filename=file.path(LASfolder,fname), format='GTiff', overwrite=T)
  }
  if(silent==FALSE) plot(chm, col=col)
  return(chm)
}
