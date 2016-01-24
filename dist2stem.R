dist2stem <- function(stems=NA, thresh.val=1, degrees=FALSE, from=FALSE) {
  
  require(raster)
  
  values(stems)[values(stems) <  thresh.val] <- NA
  values(stems)[values(stems) >= thresh.val] <- 1
  
  rws <- nrow(stems)
  cls <- ncol(stems)
  center <- cellFromRowCol(object=stems, rownr=rws/2, colnr=cls/2)
  
  # Calculate nearest distance to stem
  dist2s <- distance(x=stems)
  dist2s <- dist2s[center]
  
  # Calculate direction to nearest stem
  dir2s <- direction(x=stems, degrees=degrees, from=from, doEdge=F)
  dir2s <- dir2s[center]
  
  return(c( dist2stem=dist2s, dir2stem=dir2s) )
}
