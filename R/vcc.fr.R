#' Canopy-to-first-return Ratio Vertical Canopy Cover
#'
#' This function calculates canopy cover as the Canopy-to-first-return Ratio
#' @param las Path or name of LAS file. Defaults to NA.
#' @param thresh.val Specifies the value to use for canopy height thresholding. Defaults to 1.25.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://www.sciencedirect.com/science/article/pii/S0034425708003106}
#' @references \url{http://link.springer.com/chapter/10.1007\%2F978-94-017-8663-8_20}
#' @keywords vertical canopy cover, fractional canpopy cover, canopy cover
#' @export
#' @return The results of \code{vcc.fr}
#' @examples
#' vcc.fr(las='C:/plot.las', thresh.val=1.25, silent=FALSE)

vcc.fr <- function(las=NA, thresh.val=1.25, silent=FALSE) {

  if(length(las)==1 & any(is.na(eval(las)))) stop('Please input a full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    if(any(v == -Inf | is.na(v)==T | is.null(v)==T | v < 0)) return('blue')
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  if(!exists("las")) {
    LAS       <- rLiDAR::readLAS(las, short=FALSE)
    LASfolder <- dirname(las)
    LASname   <- strsplit(basename(las),'\\.')[[1]][1]
  } else LAS  <- las

  LAS <- LAS[order(LAS[,'ReturnNumber'], decreasing=FALSE), ]
  col <- myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,'ReturnNumber'])

  if (length(LAS[LAS[,'Z'] >= thresh.val]) < 1) { return(0) }

  result <- nrow(LAS[LAS[,'Z'] >= thresh.val & LAS[,'ReturnNumber'] == 1,]) / nrow(LAS[LAS[,'ReturnNumber'] == 1,])

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LAS[,1], LAS[,2], pch=19, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy to Return Ratio')
  }
  return(result)
}
