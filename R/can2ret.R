can2ret <- function(LASpath=NA, reprojection=NA, thresh.val=1.25, silent=FALSE) {

  if (is.na(LASpath)) stop('Please input a full file path to the LAS file')

  #require(sp)
  #require(rLiDAR)

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  LAS       <- readLAS(LASpath, short=FALSE)
  LAS       <- LAS[order(LAS[,3], decreasing=FALSE), ]
  LASfolder <- dirname(LASpath)
  LASname   <- strsplit(basename(LASpath), '\\.')[[1]][1]

  if(!is.na(reprojection)) {
    try(if(length(reprojection) != 2) stop('Please supply input and output CRS strings in the form: c(input,output)'))
    CRSin  <- reprojection[1]
    CRSout <- reprojection[2]
    spLAS  <- SpatialPoints(LAS[,c(1:3)], proj4string=CRS(CRSin))
    tmLAS  <- spTransform(spLAS, CRSobj=CRS(CRSout))
    rpLAS  <- coordinates(tmLAS)
    LAS    <- cbind(rpLAS, LAS[,c(4:12)])
  }

  col      <- myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,9])
  can2ret  <- length(LAS[,3][LAS[,3] >= thresh.val & LAS[,9] != 2]) / length(LAS[,3])

  if(silent==FALSE) {
    par(mfrow=c(1,1), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LAS[,1], LAS[,2], pch=19, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir Canopy to Return Ratio')
  }
  return(can2ret)
}
