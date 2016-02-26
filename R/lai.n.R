#' Contact Frequency and Fractional Cover-based LAI
#'
#' This function calculates LAI as a function of the effective LAI and fractional cover per contact frequency
#' @param las.path Path of LAS file. Defaults to NA.
#' @param thresh.val Specifies the value to use for canopy height thresholding. Defaults to 1.25.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://www.sciencedirect.com/science/article/pii/S0034425706001751}
#' @keywords effective LAI, LAI, leaf area index
#' @export
#' @return The results of \code{lai.n}
#' @examples
#' lai.n(las.path='C:/plot.las', thresh.val=1.25, silent=FALSE)

lai.n <- function(las.path=NA, thresh.val=1.25, silent=FALSE) {

  if (is.na(las.path)) stop('Please input a full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  LAS <- rLiDAR::readLAS(las.path, short=FALSE)
  LAS <- LAS[order(LAS[,'Intensity'], decreasing=FALSE), ]
  col <- myColorRamp(colors=c('brown','red','orange','yellow'), values=LAS[,'Intensity'])

  n.tot <- nrow(LAS)
  n.veg <- nrow(LAS[LAS[,'Z'] >= thresh.val,])
  fc    <- n.veg/n.tot

  n.first  <- nrow(LAS[LAS[,'ReturnNumber']==1 & LAS[,'Z'] >= thresh.val, ])
  n.last   <- nrow(LAS[LAS[,'ReturnNumber']==LAS[,'NumberOfReturns'] & LAS[,'Z'] >= thresh.val, ])
  n.single <- nrow(LAS[LAS[,'ReturnNumber']==1 & LAS[,'NumberOfReturns']==1 & LAS[,'Z'] >= thresh.val, ])
  lai      <- n.first/(n.last+n.single)

  result <- lai*fc

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LAS[,1], LAS[,2], pch=19, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy to Return Ratio')
  }
  return(result)
}
