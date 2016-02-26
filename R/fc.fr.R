#' Canopy-to-first-return Ratio Fractional Cover
#'
#' This function calculates fractional cover as the canopy-to-first-return ratio
#' @param las.path Path of LAS file. Defaults to NA.
#' @param thresh.val Specifies the value to use for canopy height thresholding. Defaults to 1.25.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @keywords fractional canopy cover, fractional cover, canopy cover
#' @export
#' @return The results of \code{fc.fr}
#' @examples
#' fc.fr(las.path='C:/plot.las', thresh.val=1.25, silent=FALSE)

fc.fr <- function(las.path=NA, thresh.val=1.25, silent=FALSE) {

  if (is.na(las.path)) stop('Please input a full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  LAS <- rLiDAR::readLAS(las.path, short=FALSE)
  LAS <- LAS[order(LAS[,'ReturnNumber'], decreasing=FALSE), ]
  col <- myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,'ReturnNumber'])

  result <- nrow(LAS[LAS[,'Z'] >= thresh.val & LAS[,'ReturnNumber'] == 1,]) / nrow(LAS[LAS[,'ReturnNumber'] == 1,])

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LAS[,1], LAS[,2], pch=19, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy to Return Ratio')
  }
  return(result)
}
