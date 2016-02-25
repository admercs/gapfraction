#' Canopy-to-total-return ratio Fractional Cover
#'
#' This function calculates fractional cover as the canopy-to-total-return ratio above a height threshold
#' @param las.path Path of LAS file. Defaults to NA.
#' @param thresh.val Specifies the value to use for canopy height thresholding. Defaults to 1.25.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @keywords gap fraction, voronoi, thiessen
#' @export
#' @return The results of \code{fc.rr}
#' @examples
#' fc.r(las.path='C:/plot.las', thresh.val=1.25, silent=FALSE)

fc.r <- function(las.path=NA, thresh.val=1.25, silent=FALSE) {

  if (is.na(las.path)) stop('Please input a full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  LAS    <- rLiDAR::readLAS(las.path, short=FALSE)
  LAS    <- LAS[order(LAS[,3], decreasing=FALSE), ]
  col    <- myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,3])
  result <- length(LAS[,3][LAS[,3] >= thresh.val]) / length(LAS[,3])

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LAS[,1], LAS[,2], pch=19, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy to Return Ratio')
  }
  return(result)
}
