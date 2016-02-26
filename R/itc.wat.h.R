#' Hierarchical Watershed Segmentation Individual Tree Crown Detection
#'
#' This function detects individual tree crowns using hierarchical watershed segmentation
#' @param chm.stack Name of the CHM raster stack object output from a \code{chm} function with \code{stacked=TRUE}. Defaults to NA.
#' @param res Resolution of the input canopy height model in meters, used to adjust the ht2rad function. Defaults to 1.
#' @param ht2rad Function used to convert tree height to crown radius. Defaults to NA.
#' @param min.h Canopy minimum height. Defaults to 1.
#' @param tolerance The minimum height between an object's highest point and the point where it contacts another object. If the height is smaller than the tolerance, the object will be combined with its highest neighbor. Defaults to 0.1.
#' @param fun Function used to combine rasters of tree crown position. Defaults to max.
#' @param num Boolean switch for the output of numeric values instead of a raster object. Defaults to TRUE.
#' @param stacked Boolean switch for the output of a raster stack. Defaults to FALSE.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @param ws.plot Boolean switch for the display of the watershed segmentation plot. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @keywords itc, watershed segmentation
#' @export
#' @return The results of \code{itc.wat.h}
#' @examples
#' itc.wat.h(chm.stack=chms, res=1, ht2rad=function(x) 0.15746*x/res, min.h=1, tolerance=0.1, fun=max, num=TRUE, stacked=FALSE, silent=FALSE, ws.plot=FALSE)
#' itc.wat.h(chm.stack=chms, res=1, ht2rad=function(x) 0.15746*x/res, min.h=1, tolerance=0.1, fun=length, num=TRUE, stacked=FALSE, silent=FALSE, ws.plot=FALSE)

itc.wat.h <- function(chm.stack=NA, res=1, ht2rad=NA, min.h=1, tolerance=0.1, fun=max, num=TRUE, stacked=FALSE, silent=FALSE, ws.plot=FALSE) {

  if(!nlayers(chm.stack) > 1) return(c(trees=NULL, crown.area=NULL))

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  val <- seq(from=0, to=max(values(chm.stack)[!is.na(values(chm.stack))]), length.out=length(values(chm.stack)))
  bw  <- myColorRamp(colors=c('black','white'), values=val)

  chm.out <- raster::stackApply(chm.stack, indices=c(1), fun=max, na.rm=T)
  if(max(raster::values(chm.out)[!is.na(raster::values(chm.out))]) <= min.h) return(c(trees=NA, crown.area=NA))
  chm.unstk <- raster::unstack(chm.stack)
  itc.layer <- list()

  for(i in 2:(length(chm.unstk))) {
    chm <- chm.unstk[[i]]
    if(max(raster::values(chm)[!is.na(raster::values(chm))]) <= min.h) next
    itc.layer[i] <- try(itc.wat(chm=chm, tolerance=tolerance, ht2rad=ht2rad, min.h=min.h, res=res, silent=silent, ws.plot=ws.plot, num=F))
  }

  itc.layer  <- itc.layer[!sapply(itc.layer, is.null)]
  itc.stack  <- raster::stack(itc.layer)
  itc.out    <- raster::stackApply(itc.stack, indices=c(1), fun=fun, na.rm=T)
  itc.out[itc.out==0] <- NA
  itc.crowns <- ht2rad(chm.out[!is.na(itc.out)])
  crown.area <- sum(itc.crowns) / (length(raster::values(chm.out)[!is.na(raster::values(chm.out))]) * res)
  itc.trees  <- raster::clump(itc.out, directions=8)
  ntrees     <- length(unique(raster::values(itc.trees)))
  message(ntrees,' total trees counted at a tolerance of ',tolerance,' meters')

  if(stacked==FALSE) {
    if(num==FALSE) {
      plot(chm.out, col=gray.colors(256, start=0, end=1, gamma=1, alpha=NULL))
      plot(itc.out, add=TRUE, col='red', legend=FALSE)
      pts <- raster::rasterToPoints(itc.out)
      points(pts, cex=chm.out[!is.na(itc.out)]/4, pch=10, lwd=1)
      return(itc.out)
    } else return(c( trees=ntrees, crown.area=crown.area ))
  } else if(stacked==TRUE) {
    plot(itc.stack, col=bw)
    return(itc.stack)
  }
}
