can2pix <- function(chm=NA, thresh.val=1.25) {
  
  require(raster)
  
  canopy    <- values(chm)[values(chm) > thresh.val]
  all.cells <- length(values(chm)[!is.na(values(chm))])
  can.cells <- length(canopy[!is.na(canopy)])
  output    <- can.cells/all.cells
  
  return(output)
}
