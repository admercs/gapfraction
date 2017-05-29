#' Beer's Law-modified Intensity-return Ratio of Vertical Canopy Cover
#'
#' This function calculates canopy cover per the Beer's Law-modified Intensity-return Ratio, also known as the Intensity Cover Index
#' @param las Path or name of LAS file. Defaults to NA.
#' @param thresh.val Specifies the value to use for canopy height thresholding. Defaults to 1.25.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://www.sciencedirect.com/science/article/pii/S0034425708003106}
#' @references \url{http://link.springer.com/chapter/10.1007\%2F978-94-017-8663-8_20}
#' @keywords vertical canopy cover, fractional canopy cover, canopy cover
#' @export
#' @return The results of \code{vcc.bl}
#' @examples
#' vcc.bl(las='C:/plot.las', thresh.val=1.25, silent=FALSE)

vcc.bl <- function(las=NA, thresh.val=1.25, silent=FALSE) {

  if(length(las)==1 & any(is.na(eval(las)))) stop('Please input a full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  if(!exists("las")) {
    LAS       <- lidR::readLAS(las)
    LASfolder <- dirname(las)
    LASname   <- strsplit(basename(las),'\\.')[[1]][1]
  } else LAS  <- las

  LAS <- LAS[order(LAS[,'Intensity'], decreasing=FALSE), ]
  col <- myColorRamp(colors=c('brown','red','orange','yellow'), values=LAS[,'Intensity'])

  i.total    <- sum(LAS[, 'Intensity'])
  i.ground   <- sum(LAS[LAS[,'Classification']==2, 'Intensity'])
  i.first    <- sum(LAS[LAS[,'ReturnNumber']==1, 'Intensity'])
  i.last     <- sum(LAS[LAS[,'ReturnNumber']==LAS[,'NumberOfReturns'], 'Intensity'])
  i.inter    <- sum(LAS[LAS[,'ReturnNumber']!=1 & LAS[,'ReturnNumber']!=LAS[,'NumberOfReturns'], 'Intensity'])
  i.single   <- sum(LAS[LAS[,'ReturnNumber']==1 & LAS[,'NumberOfReturns']==1, 'Intensity'])
  i.g.single <- sum(LAS[LAS[,'ReturnNumber']==1 & LAS[,'NumberOfReturns']==1 & LAS[,'Classification']==2, 'Intensity'])
  i.g.last   <- sum(LAS[LAS[,'ReturnNumber']==LAS[,'NumberOfReturns'] & LAS[,'Classification']==2, 'Intensity'])

  top.terms <- (i.g.single/i.total) + sqrt(i.g.last/i.total)
  bot.terms <- ((i.first+i.single)/i.total) + sqrt((i.inter+i.last)/i.total)
  result    <- 1-(top.terms/bot.terms)

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LAS[,1], LAS[,2], pch=19, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy to Return Ratio')
  }
  return(result)
}
