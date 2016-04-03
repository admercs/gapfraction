#' Canopy-to-total-return Ratio Fractional Cover
#'
#' This function calculates fractional cover per the canopy-to-total-return ratio
#' @param las.path Path of LAS file. Defaults to NA.
#' @param thresh.val Specifies the value to use for canopy height thresholding. Defaults to 1.25.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://link.springer.com/chapter/10.1007\%2F978-94-017-8663-8_20}
#' @keywords fractional canopy cover, fractional cover, canopy cover
#' @export
#' @return The results of \code{fc.r}
#' @examples
#' fc.r(las.path='C:/plot.las', thresh.val=1.25, silent=FALSE)

fc.r <- function(las=NA, thresh.val=1.25, silent=FALSE) {

  if (is.na(las)) stop('Please input a full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  if(!exists("las")) {
    LAS       <- rLiDAR::readLAS(las, short=FALSE)
    LASfolder <- dirname(las)
    LASname   <- strsplit(basename(las),'\\.')[[1]][1]
  } else LAS  <- las

  LAS <- LAS[order(LAS[,'Classification'], decreasing=FALSE), ]
  col <- myColorRamp(colors=c('blue','orange','yellow','green','red'), values=LAS[,'Classification'])

  result <- nrow(LAS[LAS[,'Z'] >= thresh.val,]) / nrow(LAS)

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LAS[,1], LAS[,2], pch=19, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy to Return Ratio')
  }
  return(result)
}
