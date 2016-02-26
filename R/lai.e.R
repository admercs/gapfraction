#' Simple Effective LAI Estimation with Beer's Law, the Ground-to-total-returns Ratio, and a Spherical Leaf Angle Distribution
#'
#' This function calculates effective LAI using Beer's Law with the ground-to-total-returns ratio and a spherical leaf angle distribution
#' @param las.path Path of LAS file. Defaults to NA.
#' @param k Specifies the leaf angle distribution to use. Defaults to 0.5 for spherical.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://www.sciencedirect.com/science/article/pii/S0168192309000409}
#' @keywords effective LAI, LAI, leaf area index
#' @export
#' @return The results of \code{lai.e}
#' @examples
#' lai.e(las.path='C:/plot.las', k=0.5, silent=FALSE)

lai.e <- function(las.path=NA, k=0.5, silent=FALSE) {

  if (is.na(las.path)) stop('Please input a full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  LAS <- rLiDAR::readLAS(las.path, short=FALSE)
  LAS <- LAS[order(LAS[,'Classification'], decreasing=FALSE), ]
  col <- myColorRamp(colors=c('blue','red','orange','green','purple'), values=LAS[,'Classification'])

  n.total  <- nrow(LAS)
  n.ground <- nrow(LAS[LAS[,'Classification'] == 2,])
  result   <- (-1/k)*log(n.ground/n.total)

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LAS[,1], LAS[,2], pch=19, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy to Return Ratio')
  }
  return(result)
}
