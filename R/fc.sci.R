#' Solberg's Cover Index of Fractional Cover
#'
#' This function calculates fraction cover per Solberg's Cover Index
#' @param las.path Path of LAS file. Defaults to NA.
#' @param thresh.val Specifies the value to use for canopy height thresholding. Defaults to 1.25.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @keywords fractional canopy cover, fractional cover, canopy cover
#' @export
#' @return The results of \code{fc.sci}
#' @examples
#' fc.sci(las.path='C:/plot.las', thresh.val=1.25, silent=FALSE)

fc.sci <- function(las.path=NA, thresh.val=1.25, silent=FALSE) {

  if (is.na(las.path)) stop('Please input a full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  LAS <- rLiDAR::readLAS(las.path, short=FALSE)
  LAS <- LAS[order(LAS[,'Intensity'], decreasing=FALSE), ]
  col <- myColorRamp(colors=c('brown','red','orange','yellow'), values=LAS[,'Intensity'])

  all.single <- nrow(LAS[LAS[,'ReturnNumber']==1 & LAS[,'NumberOfReturns']==1, ])
  all.first  <- nrow(LAS[LAS[,'ReturnNumber']==1, ])
  all.last   <- nrow(LAS[LAS[,'ReturnNumber']==LAS[,'NumberOfReturns'], ])

  can.single <- nrow(LAS[LAS[,'ReturnNumber']==1 & LAS[,'NumberOfReturns']==1 & LAS[,'Z'] >= thresh.val, ])
  can.first  <- nrow(LAS[LAS[,'ReturnNumber']==1 & LAS[,'Z'] >= thresh.val, ])
  can.last   <- nrow(LAS[LAS[,'ReturnNumber']==LAS[,'NumberOfReturns'] & LAS[,'Z'] >= thresh.val, ])

  top.terms <- can.single+0.5*(can.first+can.last)
  bot.terms <- all.single+0.5*(all.first+all.last)
  result    <- top.terms/bot.terms

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LAS[,1], LAS[,2], pch=19, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy to Return Ratio')
  }
  return(result)
}
