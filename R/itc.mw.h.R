#' Hierarchical Moving Window Individual Tree Crown Detection
#'
#' This function detects individual tree crowns using a moving window function with variable window diameters
#' @param chm.stack Name of the CHM raster stack output from a \code{chm} function with \code{stacked=TRUE}. Defaults to NA.
#' @param res Resolution of the input canopy height model in meters, used to adjust the ht2rad function. Defaults to 1.
#' @param ht2rad Function used to convert tree height to crown radius. Defaults to NA.
#' @param min.h Canopy minimum height. Defaults to 1.
#' @param type Shape parameter for the moving window function. Options are circle, Gauss, and rectangle. Defaults to circle.
#' @param fun Function used for combining raster layers of crown position. Defaults to max.
#' @param num Boolean switch for the output of numeric values instead of a raster object. Defaults to TRUE.
#' @param stacked Boolean switch for the output of a raster stack of tree crown positions. Defaults to FALSE.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @param plots Boolean switch for the display of plots. Defaults to FALSE.
#' @param geoTIFF Boolean switch for saving GeoTIFF raster of tree positions to a folder, if the chm parameter is a folder. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://www.mdpi.com/2072-4292/4/4/950/htm}
#' @keywords itc, moving window
#' @export
#' @return The results of \code{itc.mw.h}
#' @examples
#' itc.mw.h(chm.stack=chms, res=1, ht2rad=function(x) 0.15746*x/res, min.h=1, type='circle', fun=max, num=TRUE, stacked=FALSE, silent=FALSE, plots=FALSE, geoTIFF=FALSE)

itc.mw.h <- function(chm.stack=NA, res=1, ht2rad=NA, min.h=1, type='circle', fun=max, num=TRUE, stacked=FALSE, silent=FALSE, plots=FALSE, geoTIFF=FALSE) {

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
    itc.layer[[i]] <- itc.mw(chm=chm, ht2rad=ht2rad, min.h=min.h, type=type, res=res, silent=silent, plots=plots, geoTIFF=geoTIFF, num=FALSE)
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
