#' Canopy-to-total-pixel Ratio Vertical Canopy Cover
#'
#' This function calculates canopy cover as the ratio of canopy pixels to non-canopy pixels
#' @param chm Name of the CHM raster object output from a CHM function with stacked=FALSE. Defaults to NA.
#' @param thresh.val Threshold value used for minimum canopy height. Defaults to 1.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://www.sciencedirect.com/science/article/pii/S0924271614001646}
#' @references \url{http://treephys.oxfordjournals.org/content/early/2009/06/24/treephys.tpp042.full}
#' @keywords vertical canopy cover, fractional canopy cover, canopy cover
#' @export
#' @return The results of \code{vcc.p}
#' @examples
#' vcc.p(chm=chm, thresh.val=1.25, silent=FALSE)

vcc.p <- function(chm=NA, thresh.val=1.25, silent=FALSE) {

  if(length(las)==1 & any(is.na(eval(las)))) stop('Please input a CHM')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  canopy    <- raster::values(chm)[raster::values(chm) >= thresh.val]
  all.cells <- length(raster::values(chm)[!is.na(raster::values(chm))])
  can.cells <- length(canopy[!is.na(canopy)])
  result    <- can.cells/all.cells

  val <- seq(from=0, to=max(raster::values(chm)[!is.na(raster::values(chm))]), length.out=length(raster::values(chm)[!is.na(raster::values(chm))]))
  col <- myColorRamp(colors=c('blue','green','yellow','red'), values=val)

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(chm, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy Ratio')
  }
  return(result)
}
