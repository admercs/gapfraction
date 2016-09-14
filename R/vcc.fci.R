#' First-echo Cover Index of Vertical Canopy Cover
#'
#' This function calculates fractional cover per the First-echo Cover Index
#' @param las Path or name of LAS file. Defaults to NA.
#' @param thresh.val Specifies the value to use for canopy height thresholding. Defaults to 1.25.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://geography.swan.ac.uk/silvilaser/papers/oral_papers/Forestry\%20Applications\%20\%26\%20Inventory/Holmgren.pdf}
#' @references \url{http://link.springer.com/chapter/10.1007\%2F978-94-017-8663-8_20}
#' @keywords vertical canopy cover, fractional canopy cover, canopy cover
#' @export
#' @return The results of \code{vcc.fci}
#' @examples
#' vcc.fci(las='C:/plot.las', thresh.val=1.25, silent=FALSE)

vcc.fci <- function(las=NA, thresh.val=1.25, silent=FALSE) {

  if(length(las)==1 & any(is.na(eval(las)))) stop('Please input a full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    if(any(v == -Inf | is.na(v)==T | is.null(v)==T | v < 0)) return('brown')
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  if(!exists("las")) {
    LAS       <- rLiDAR::readLAS(las, short=FALSE)
    LASfolder <- dirname(las)
    LASname   <- strsplit(basename(las),'\\.')[[1]][1]
  } else LAS  <- las

  LAS <- LAS[order(LAS[,'ReturnNumber'], decreasing=FALSE), ]
  col <- myColorRamp(colors=c('brown','red','orange','yellow'), values=LAS[,'ReturnNumber'])

  if (length(LAS[LAS[,'Z'] >= thresh.val]) < 1) { return(0) }

  all.first  <- nrow(LAS[LAS[,'ReturnNumber']==1, ])
  all.single <- nrow(LAS[LAS[,'ReturnNumber']==1 & LAS[,'NumberOfReturns']==1, ])
  can.first  <- nrow(LAS[LAS[,'ReturnNumber']==1 & LAS[,'Z'] >= thresh.val, ])
  can.single <- nrow(LAS[LAS[,'ReturnNumber']==1 & LAS[,'NumberOfReturns']==1 & LAS[,'Z'] >= thresh.val, ])

  result <- (can.single+can.first)/(all.single+all.first)

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LAS[,1], LAS[,2], pch=19, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy to Return Ratio')
  }
  return(result)
}
