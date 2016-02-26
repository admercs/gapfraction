#' Watershed Segmentation Individual Tree Crown Detection
#'
#' This function detects individual tree crowns using watershed segmentation
#' @param chm Name of the CHM raster object output from a \code{chm} function with \code{stacked=FALSE}. Defaults to NA.
#' @param res Resolution of the input canopy height model in meters, used to adjust the ht2rad function. Defaults to 1.
#' @param ht2rad Function used to convert tree height to crown radius. Defaults to NA.
#' @param min.h Canopy minimum height. Defaults to 1.
#' @param tolerance The minimum height between an object's highest point and the point where it contacts another object. If the height is smaller than the tolerance, the object will be combined with its highest neighbor. Defaults to 0.1.
#' @param num Boolean switch for the output of numeric values instead of a raster object. Defaults to TRUE.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @param ws.plot Boolean switch for the display of the watershed segmentation plot. Defaults to FALSE.
#' @param geoTIFF Boolean switch for saving GeoTIFF raster of tree positions to a folder, if the chm parameter is a folder. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://www.mdpi.com/2072-4292/4/4/950/htm}
#' @keywords itc, watershed segmentation
#' @export
#' @return The results of \code{itc.wat}
#' @examples
#' itc.wat(chm=chm, res=1, ht2rad=function(x) 0.15746*x/res, min.h=1, tolerance=0.1, num=TRUE, silent=FALSE, ws.plot=FALSE, geoTIFF=FALSE)

itc.wat <- function(chm=NA, res=1, ht2rad=NA, min.h=1, tolerance=0.1, num=TRUE, silent=FALSE, ws.plot=FALSE, geoTIFF=FALSE) {

  if(max(raster::values(chm)[!is.na(raster::values(chm))]) <= min.h) return(c(trees=NA, crown.area=NA))
  if (!'EBImage' %in% installed.packages()) {
    source('https://bioconductor.org/biocLite.R')
    biocLite('EBImage')
  }
  require(EBImage)

  if(is.character(chm)) {
    isPath <- TRUE
    folder <- dirname(chm)
    name   <- strsplit(basename(chm), '\\.')[[1]][1]
    chm    <- raster::raster(chm)
  } else isPath <- FALSE

  img <- EBImage::flip(EBImage::Image(data=values(chm), dim=c(ncol(chm),nrow(chm),1), colormode='Grayscale'))
  wat <- EBImage::watershed(img, tolerance=tolerance, ext=1)

  wshed <- sort(unique(as.numeric(imageData(wat))))
  trees <- sapply(wshed, function(x) max(img[!is.na(img) & wat==x]))

  itc.trees <- img
  for(i in 1:length(wshed)) {
    itc.trees[!is.na(itc.trees) & itc.trees != 0 & itc.trees > min.h & wat==wshed[i] & itc.trees==trees[i]] <- 9999
  }

  itc.trees[itc.trees != 9999] <- NA
  itc.trees[itc.trees == 9999] <- 1
  ntrees <- sum(itc.trees[!is.na(itc.trees)])
  message(ntrees,' trees counted at a tolerance of ',tolerance,' meters')
  itc.crowns <- ht2rad(img[itc.trees==1])
  crown.area <- sum(itc.crowns[!is.na(itc.crowns)]) / (length(img[!is.na(img)]) * res)

  itc.out <- raster::flip(raster(t(matrix(itc.trees, ncol=ncol(chm), nrow=nrow(chm))), crs=crs(chm)), direction=2)
  extent(itc.out) <- raster::extent(chm)

  if(silent==FALSE) {
    plot(chm, col=gray.colors(256, start=0, end=1, gamma=1, alpha=NULL))
    plot(itc.out, add=T, col='red', legend=F)
    pts <- raster::rasterToPoints(itc.out)
    points(pts, cex=chm[!is.na(itc.out)]/4, pch=10, lwd=1)
  }

  if(ws.plot==TRUE) {
    watcol <- wat
    watcol[is.na(watcol)] <- 9999
    watcol <- EBImage::channel(watcol, mode='gray')
    EBImage::display(colorLabels(watcol))
  }

  if(isPath==TRUE & geoTIFF==TRUE) {
    fname <- paste(name,'_itc.tiff',sep='')
    raster::writeRaster(x=itc.out, filename=file.path(folder,fname), format='GTiff', overwrite=T)
  } else if(isPath==FALSE & geoTIFF==TRUE) {
    raster::writeRaster(x=itc.out, filename='itc.tiff', format='GTiff', overwrite=T)
  }

  if(num==TRUE) {
    return( c(trees=ntrees, crown.area=crown.area) )
  } else {
    return(itc.out)
  }
}
