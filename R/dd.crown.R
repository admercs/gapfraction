#' Nearest Crown Distance and Direction from Plot Center
#'
#' This function calculates the distance and direction to the nearest crown pixel from the plot center
#' @param crowns Raster object of tree crown positions. Defaults to NA.
#' @param thresh.val Threshold value used for tree crown selection. Defaults to 1.
#' @param degrees Boolean switch for the output of direction values in degrees rather than radians. Defaults to FALSE.
#' @param from Boolean switch for the output of direction values from nearest crowns rather than to nearest crowns. Defaults to FALSE.
#' @keywords crown, distance, direction
#' @export
#' @return The results of \code{dd.crown}
#' @examples
#' dd.crown(crowns=itc.crowns, thresh.val=1, degrees=FALSE, from=FALSE)

dd.crown <- function(crowns=NA, thresh.val=1, degrees=FALSE, from=FALSE) {

  if(class(crowns)=='logical') return(c(c.dist=NA, c.dir=NA))
  if(length(raster::values(crowns)[!is.na(raster::values(crowns))]) < 1) return(c(c.dist=NA, c.dir=NA))

  raster::values(crowns)[raster::values(crowns) <  thresh.val] <- NA
  raster::values(crowns)[raster::values(crowns) >= thresh.val] <- 1

  rws <- nrow(crowns)
  cls <- ncol(crowns)

  center <- raster::cellFromRowCol(object=crowns, rownr=rws/2, colnr=cls/2)

  c.dist <- raster::distance(x=crowns)
  c.dist <- c.dist[center]
  c.dir  <- raster::direction(x=crowns, degrees=degrees, from=from, doEdge=FALSE)
  c.dir  <- c.dir[center]

  return(c(c.dist=c.dist, c.dir=c.dir))
}
