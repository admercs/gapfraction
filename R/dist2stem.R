dist2stem <- function(stems=NA, thresh.val=1, degrees=FALSE, from=FALSE) {

  if(class(stems)=='logical') return(c(dist2stem=NA, dir2stem=NA))
  if(length(raster::values(stems)[!is.na(raster::values(stems))]) < 1) return(c(dist2stem=NA, dir2stem=NA))

  raster::values(stems)[raster::values(stems) <  thresh.val] <- NA
  raster::values(stems)[raster::values(stems) >= thresh.val] <- 1

  rws <- nrow(stems)
  cls <- ncol(stems)

  center <- raster::cellFromRowCol(object=stems, rownr=rws/2, colnr=cls/2)

  dist2s <- raster::distance(x=stems)
  dist2s <- dist2s[center]

  dir2s <- raster::direction(x=stems, degrees=degrees, from=from, doEdge=F)
  dir2s <- dir2s[center]

  return(c( dist2stem=dist2s, dir2stem=dir2s) )
}
