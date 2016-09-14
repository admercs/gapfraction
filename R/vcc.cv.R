#' Cartesian-Voronoi Vertical Canopy Cover
#'
#' This function implements the Cartesian-Voronoi vertical canopy cover algorithm of Alexander (2013)
#' @param las Path or name of LAS file. Defaults to NA.
#' @param reprojection Proj4 projection string to use for reprojection. Defaults to NA.
#' @param col Specifies the LiDAR metric to use to color points of first plot in display. Options include height, intensity, nreturn, and class. Defaults to height.
#' @param col2 Specifies the LiDAR metric to use to color points of second plot in display. Defaults to NA.
#' @param thresh.var Specifies the LiDAR metric to use for thresholding canopy points. Options include height, intensity, nreturn, and class. Defaults to height.
#' @param thresh.val Specifies the value to use for thresholding. Defaults to 1.25.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @param plots Boolean switch for the saving of plot files to the las.path folder. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references \url{http://www.sciencedirect.com/science/article/pii/S0034425713000771}
#' @keywords vertical canopy cover, fractional canopy cover, canopy cover
#' @export
#' @return The results of \code{vcc.cv}
#' @examples
#' vcc.cv(las='C:/plot.las', reprojection=NA, col='height', col2='intensity', thresh.var='height', thresh.val=1.25, silent=FALSE, plots=FALSE)

vcc.cv <- function(las=NA, reprojection=NA, col='height', col2=NA, thresh.var='height', thresh.val=1.25, silent=TRUE, plots=FALSE) {

  if(length(las)==1 & any(is.na(eval(las)))) stop('Please input a full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  if(!exists("las")) {
    LAS       <- rLiDAR::readLAS(las, short=FALSE)
    LASfolder <- dirname(las)
    LASname   <- strsplit(basename(las),'\\.')[[1]][1]
  } else LAS  <- las
  LAS  <- LAS[order(LAS[,3], decreasing=FALSE),]

  if(!is.na(reprojection)) {
    try(if(length(reprojection) != 2) stop('Please supply input and output CRS strings in the form: c(input,output)'))
    CRSin  <- reprojection[1]
    CRSout <- reprojection[2]
    spLAS  <- sp::SpatialPoints(LAS[,c(1:3)], proj4string=CRS(CRSin))
    tmLAS  <- sp::spTransform(spLAS, CRSobj=CRS(CRSout))
    rpLAS  <- sp::coordinates(tmLAS)
    LAS    <- cbind(rpLAS, LAS[,c(4:12)])
  }

  x    <- (LAS[,1]-min(LAS[,1])) - (diff(range(LAS[,1]))/2)
  y    <- (LAS[,2]-min(LAS[,2])) - (diff(range(LAS[,2]))/2)
  z    <-  LAS[,3]

  dupl <- !deldir::duplicatedxy(x=x, y=y)

  x <- x[dupl]
  y <- y[dupl]
  z <- z[dupl]

  nrows <- length(x)
  col   <- rep(col,  nrows)
  col2  <- rep(col2, nrows)

  col <- ifelse(col == 'height',    myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,3]),
         ifelse(col == 'intensity', myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,4]),
         ifelse(col == 'nreturn',   myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,5]),
         ifelse(col == 'class',     myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,9]),
         'forestgreen'))))

  col2 <- ifelse(col2 == 'height',    myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,3]),
          ifelse(col2 == 'intensity', myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,4]),
          ifelse(col2 == 'nreturn',   myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,5]),
          ifelse(col2 == 'class',     myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,9]),
          'forestgreen'))))

  xy <- matrix(c(x,y), nrow=length(x), ncol=2, dimnames=list(c(1:length(x)),c('x','y')))
  md <- geometry::delaunayn(xy, full=F)

  thresh.var <- rep(thresh.var, nrows)

  thresh.var <- ifelse(thresh.var == 'height',    LAS[,3],
                ifelse(thresh.var == 'intensity', LAS[,4],
                ifelse(thresh.var == 'nreturn',   LAS[,5],
                ifelse(thresh.var == 'class',     LAS[,9],
                z))))

  canopy <- ifelse(thresh.var >= thresh.val, 1, 0)
  mv     <- deldir::deldir(x=x, y=y, z=canopy, rw=NULL, eps=1e-09, plotit=FALSE, suppressMsge=TRUE)
  result <- sum(mv$summary$dir.area*mv$summary$z) / mv$dir.area

  cvex  <- spatstat::convexhull.xy(matrix(c(x=x,y=y),ncol=2))
  clipp <- cvex$bdry[[1]]

  point.symbols <- 19

  if (plots == TRUE) {

    jpeg(file.path(LASfolder, paste(LASname,'_closure.jpg',sep='')), width=8, height=8, units='in', res=300, quality=100)
    par(mfrow=c(2,2), mar=c(2,2,3,2), pty='s', xpd=TRUE)

    plot(x, y, pch=point.symbols, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir')
    plot(x, y, pch=point.symbols, col=col2, bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir')

    plot(x, y, pch=point.symbols, col=col, bty='n', xlab='Latitude', ylab='Longitude', main='Delaunay Cartesian Nadir')
    geometry::trimesh(md, xy, add=TRUE)

    plot(x, y, pch=point.symbols, col=col, bty='n', xlab='Latitude', ylab='Longitude', main='Voronoi Cartesian Nadir')
    fillcol <- ifelse(thresh.var >= thresh.val & LAS[,9] != 2, col, NA)
    deldir::plot.tile.list(deldir::tile.list(mv), verbose=FALSE, close=TRUE, pch=NA, fillcol=fillcol, col.pts=NA, border='white',
                           showpoints=FALSE, add=TRUE, asp=1, clipp=clipp, alpha=0.5)
    dev.off()
  }
  if(silent==FALSE) {
    par(mfrow=c(2,2), mar=c(2,2,3,2), pty='s', xpd=TRUE)

    plot(x, y, pch=point.symbols, col=col,  bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir')
    plot(x, y, pch=point.symbols, col=col2, bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir')

    plot(x, y, pch=point.symbols, col=col, bty='n', xlab='Latitude', ylab='Longitude', main='Delaunay Cartesian Nadir')
    geometry::trimesh(md, xy, add=TRUE)

    plot(x, y, pch=point.symbols, col=col, bty='n', xlab='Latitude', ylab='Longitude', main='Voronoi Cartesian Nadir')
    fillcol <- ifelse(thresh.var >= thresh.val & LAS[,9] != 2, col, NA)
    deldir::plot.tile.list(deldir::tile.list(mv), verbose=FALSE, close=TRUE, pch=NA, fillcol=fillcol, col.pts=NA, border='white',
                   showpoints=FALSE, add=TRUE, asp=1, clipp=clipp, alpha=0.5)

    par(mfrow=c(1,1))
  }
  return(result)
}
