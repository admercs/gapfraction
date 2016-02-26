#' Nearest Canopy Pixel Distance and Direction from Plot Center
#'
#' This function calculates the distance and direction to the nearest crown from the plot center
#' @param chm Name of the CHM raster object output from a CHM function with stacked=FALSE. Defaults to NA.
#' @param thresh.val Threshold value used for minimum canopy height. Defaults to 1.
#' @param degrees Boolean switch for the output of direction values in degrees rather than radians. Defaults to FALSE.
#' @param from Boolean switch for the output of direction values from nearest crowns rather than to nearest crowns. Defaults to FALSE.
#' @keywords canopy, distance, direction
#' @export
#' @return The results of \code{dd.canopy}
#' @examples
#' dd.canopy(chm=chm, thresh.val=1.25, degrees=FALSE, from=FALSE)

dd.canopy <- function(chm=NA, thresh.val=1.25, degrees=FALSE, from=FALSE) {

  if(max(raster::values(chm)[!is.na(raster::values(chm))]) < thresh.val) return(c(can.dist=NA, can.dir=NA))

  canopy <- chm
  raster::values(canopy)[raster::values(canopy) <  thresh.val] <- NA
  raster::values(canopy)[raster::values(canopy) >= thresh.val] <- 1

  rws <- nrow(chm)
  cls <- ncol(chm)
  center <- raster::cellFromRowCol(object=canopy, rownr=rws/2, colnr=cls/2)

  can.dist <- raster::distance(x=canopy)
  can.dist <- can.dist[center]
  can.dir  <- raster::direction(x=canopy, degrees=degrees, from=from, doEdge=F)
  can.dir  <- can.dir[center]

  return(c(can.dist=can.dist, can.dir=can.dir))
}
