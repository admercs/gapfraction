itc.watershed <- function(chm=NA, ht2rad=NA, tolerance=0.1, res=1, num=TRUE, silent=FALSE, ws.plot=FALSE) {
  
  require(raster)
  if (!'EBImage' %in% installed.packages()) {
    source('https://bioconductor.org/biocLite.R')
    biocLite('EBImage')
  }
  require(EBImage)
  
  img <- flip(Image(data=values(chm), dim=c(ncol(chm),nrow(chm),1), colormode='Grayscale'))
  wat <- watershed(img, tolerance=0.1, ext=1)
  
  wshed <- sort(unique(as.numeric(imageData(wat))))
  trees <- sapply(wshed, function(x) max(img[!is.na(img) & wat==x]))
  
  itc.trees <- img
  for(i in 1:length(wshed)) {
    itc.trees[!is.na(itc.trees) & img!=0 & wat==wshed[i] & itc.trees==trees[i]] <- 9999
  }
  itc.trees[itc.trees != 9999] <- NA
  itc.trees[itc.trees == 9999] <- 1
  ntrees     <- length(itc.trees[!is.na(itc.trees) & itc.trees==1])
  message(ntrees,' trees counted at a tolerance of ',tolerance,' meters')
  itc.crowns <- ht2rad(img[itc.trees==1])
  crown.area <- sum(itc.crowns[!is.na(itc.crowns)]) / (length(img[!is.na(img)]) * res^2) * 100
  
  crown.ras <- raster::flip(raster(t(matrix(itc.trees, ncol=ncol(chm), nrow=nrow(chm))), crs=crs(chm)), direction=2)
  extent(crown.ras) <- extent(chm)
  
  if(silent==FALSE) {
    plot(chm, col=gray.colors(256, start=0, end=1, gamma=1, alpha=NULL))
    plot(crown.ras, add=T, col='red', legend=F)
    pts <- rasterToPoints(crown.ras)
    points(pts, cex=chm[!is.na(crown.ras)]/4, pch=10, lwd=1)
  }
  
  if(ws.plot==TRUE) {
    watcol <- wat
    watcol[is.na(watcol)] <- 9999
    watcol <- channel(watcol, mode='gray')
    display(colorLabels(watcol))
  }
  
  if(num==TRUE) {
    return( c(trees=ntrees, crown.area=crown.area) )
  } else {
    return(crown.ras)
  }
}
