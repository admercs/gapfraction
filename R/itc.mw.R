#' Moving Window Individual Tree Crown Detection
#'
#' This function detects individual tree crowns using a moving window function with variable window sizes
#' @param chm Name of the CHM raster object output from a \code{chm} function with \code{stacked=FALSE}. Defaults to NA.
#' @param res Resolution of the input canopy height model in meters, used to adjust the ht2rad function. Defaults to 1.
#' @param ht2rad Function used to convert tree height to crown radius. Defaults to NA.
#' @param min.h Canopy minimum height. Defaults to 1.
#' @param type Shape parameter for the moving window function. Options are circle, Gauss, and rectangle. Defaults to circle.
#' @param num Boolean switch for the output of numeric values instead of a raster object. Defaults to TRUE.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @param plots Boolean switch for the display of plots. Defaults to FALSE.
#' @param geoTIFF Boolean switch for saving GeoTIFF raster of tree positions to a folder, if the chm parameter is a folder. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @keywords itc, moving window
#' @export
#' @return The results of \code{itc.mw}
#' @examples
#' itc.mw(chm=chm, res=1, ht2rad=function(x) 0.15746*x/res, min.h=1, type='circle', num=TRUE, silent=FALSE, plots=FALSE, geoTIFF=FALSE)

itc.mw <- function(chm=NA, res=1, ht2rad=NA, min.h=1, type='circle', num=TRUE, silent=FALSE, plots=FALSE, geoTIFF=FALSE) {

  if(max(raster::values(chm)[!is.na(raster::values(chm))]) <= min.h) return(c(trees=NA, crown.area=NA))

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  run.focal <- function(chm=chm, hts=hts, rd2=rd2, rad, type=type) {
    htt <- hts[rd2==rad | rd2==rad-1]
    x2  <- chm > min(htt) & chm < max(htt)
    x3  <- x2 * chm
    fun <- function(z) ifelse(z[length(z)/2 + 0.5]==suppressWarnings(max(z[!is.na(z) & z > min.h])), 1, NA)
    wts <- raster::focalWeight(x=x3, d=rad, type=type)
    res <- raster::focal(x=chm, w=wts, fun=fun, na.rm=FALSE, pad=rad, padValue=NA, NAonly=FALSE)
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

  hts <- sort(unique(round(chm[!is.na(chm) & chm > min.h])))
  rd1 <- ht2rad(hts)
  rd2 <- round(rd1)
  rd3 <- sort(unique(rd2))
  for(i in 1:length(rd3)) rd3[i] <- ifelse(rd3[i] %% 2 != 0, rd3[i], rd3[i] + 1)
  if(length(rd3)==0) stop('Trees of no suitable height classes exist for the crown area moving window')
  message('Computing ', length(rd3), ' moving window(s)')

  if(length(rd3)==1) {
    itc.out <- run.focal(chm=chm, hts=hts, rd2=rd2, rad=rd3, type=type)
  } else {
    itc.stk <- raster::stack()
    for(i in 1:length(rd3)) {
      itc.new <- run.focal(chm=chm, hts=hts, rd2=rd2, rad=rd3[i], type=type)
      itc.stk <- raster::stack(itc.stk, itc.new)
    }
    itc.out <- raster::stackApply(itc.stk, indices=c(1), fun=max, na.rm=T)
  }

  itc.trees  <- raster::clump(itc.out, directions=8)
  ntrees     <- length(unique(raster::values(itc.trees)[!is.na(raster::values(itc.trees))]))
  message('There are an estimated ',ntrees,' trees in the plot')
  itc.crowns <- ht2rad(chm) * itc.out
  crown.area <- sum(raster::values(itc.crowns)[!is.na(raster::values(itc.crowns))]) / (length(raster::values(chm)[!is.na(raster::values(chm))]) * res)

  val <- seq(from=0, to=max(raster::values(chm)[!is.na(raster::values(chm))]), length.out=length(raster::values(chm)))
  bw  <- myColorRamp(colors=c('black','white'), values=val)

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
  if(silent==FALSE) {
    plot(chm, col=bw)
    plot(itc.out, col='red', add=T, legend=F)
    pts <- raster::rasterToPoints(itc.out)
    points(pts, cex=chm[!is.na(itc.out)]/4, pch=10, lwd=1)
  }
  if(num==TRUE) {
    return( c(trees=ntrees, crown.area=crown.area))
  } else return(itc.out)
}
