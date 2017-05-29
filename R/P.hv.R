#' Hemispherical-Voronoi Gap Fraction
#'
#' This function implements Erickson's hemispherical-Voronoi gap fraction algorithm with four common lens geometries: equi-distant, equi-angular, stereographic, and orthographic
#' @param las Path or name of LAS file. Defaults to NA.
#' @param model Hemispherical lens geometry model to use. Options include equi-distant (\code{"equidist"}), equi-angular (\code{"equiangle"}), stereographic (\code{"stereo"}), and orthographic (\code{"ortho"}). Defaults to \code{"equidist"}.
#' @param thresh.val Specifies the value to use for thresholding. Defaults to 1.25.
#' @param thresh.var Specifies the LiDAR metric to use for thresholding canopy points. Options include height, intensity, nreturn, and class. Defaults to height.
#' @param reprojection Proj4 projection string to use for reprojection. Defaults to NA.
#' @param pol.deg Specifies the polar resolution for the radial plot lines. Defaults to 5.
#' @param azi.deg Specifies the azimuthal resolution for the radial plot lines. Defaults to 45.
#' @param col Specifies the LiDAR metric to use to color points of first plot in display. Options include height, intensity, nreturn, and class. Defaults to height.
#' @param plots Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @param plots.each Boolean switch for displaying individual of plots. Defaults to FALSE.
#' @param plots.save Boolean switch for the saving of plot files to the las.path folder. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @references Forthcoming
#' @keywords gap fraction, voronoi, thiessen
#' @export
#' @return The results of \code{P.hv}
#' @examples
#' P.hv(las='C:/plot.las', model='equidist', thresh.val=1.25, thresh.var='height', reprojection=NA, pol.deg=5, azi.deg=45, col='height', plots=TRUE, plots.each=FALSE, plots.save=FALSE)

P.hv <- function(las=NA, model='equidist', thresh.val=1.25, thresh.var='height', reprojection=NA, pol.deg=5, azi.deg=45, col='height', plots=TRUE, plots.each=FALSE, plots.save=FALSE) {

  if(length(las)==1 & any(is.na(eval(las)))) stop('Please input a full file path to the LAS file')

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  crt2str <- function(m) {
    x     <- m[,1]
    y     <- m[,2]
    z     <- m[,3]
    xi    <- x/(sqrt(x^2+y^2+z^2)+z)
    yi    <- y/(sqrt(x^2+y^2+z^2)+z)
    nrows <- length(x)
    out   <- matrix(c(xi,yi), nrow=nrows, ncol=2, dimnames=list(c(1:nrows),c('x','y')))
    return(out)
  }

  crt2eqd <- function(m) {
    x     <- m[,1]
    y     <- m[,2]
    z     <- m[,3]
    xi    <- (x/(sqrt(x^2+y^2)))*atan2(sqrt(x^2+y^2),z)
    yi    <- (y/(sqrt(x^2+y^2)))*atan2(sqrt(x^2+y^2),z)
    nrows <- length(x)
    out   <- matrix(c(xi,yi), nrow=nrows, ncol=2, dimnames=list(c(1:nrows),c('x','y')))
    return(out)
  }

  crt2ort <- function(m) {
    x     <- m[,1]
    y     <- m[,2]
    z     <- m[,3]
    xi    <- x/sqrt(x^2+y^2+z^2)
    yi    <- y/sqrt(x^2+y^2+z^2)
    nrows <- length(x)
    out   <- matrix(c(xi,yi), nrow=nrows, ncol=2, dimnames=list(c(1:nrows),c('x','y')))
    return(out)
  }

  crt2esa <- function(m) {
    x     <- m[,1]
    y     <- m[,2]
    z     <- m[,3]
    xi    <- x/sqrt(2*(x^2+y^2))*sqrt(1-(z/sqrt(x^2+y^2+z^2)))
    yi    <- y/sqrt(2*(x^2+y^2))*sqrt(1-(z/sqrt(x^2+y^2+z^2)))
    nrows <- length(x)
    out   <- matrix(c(xi,yi), nrow=nrows, ncol=2, dimnames=list(c(1:nrows),c('x','y')))
    return(out)
  }

  deg2rad <- function(x) return(x*(pi/180))

  if(!exists("las")) {
    LAS       <- lidR::readLAS(las)
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

  x <- (LAS[,1]-min(LAS[,1])) - (diff(range(LAS[,1]))/2)
  y <- (LAS[,2]-min(LAS[,2])) - (diff(range(LAS[,2]))/2)
  z <-  LAS[,3]

  nrows <- length(x)
  m     <- matrix(c(x,y,z), nrow=nrows, ncol=3, dimnames=list(c(1:nrows),c('X','Y','Z')))

  col <- rep(col, nrows)
  if(!any(col %in% c('height','intensity','nreturn','class'))) col <- 'forestgreen'
  if(any(col == 'height'))    col <- myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,3])
  if(any(col == 'intensity')) col <- myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,4])
  if(any(col == 'nreturn'))   col <- myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,5])
  if(any(col == 'class'))     col <- myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,9])

  thresh.var <- rep(thresh.var, nrows)
  if(!any(thresh.var %in% c('height','intensity','nreturn','class'))) thresh.var <- z
  if(any(thresh.var == 'height'))    thresh.var <- LAS[,3]
  if(any(thresh.var == 'intensity')) thresh.var <- LAS[,4]
  if(any(thresh.var == 'nreturn'))   thresh.var <- LAS[,5]
  if(any(thresh.var == 'class'))     thresh.var <- LAS[,9]

  LASord <- LAS[order(LAS[,3], decreasing=FALSE),]
  LAScol <- myColorRamp(colors=c('blue','green','yellow','red'), values=LASord[,3])

  if(model == 'stereo')    xy <- crt2str(m)
  if(model == 'ortho')     xy <- crt2ort(m)
  if(model == 'equidist')  xy <- crt2eqd(m)
  if(model == 'equiangle') xy <- crt2esa(m)

  md    <- geometry::delaunayn(xy, full=F)
  cvex  <- spatstat::convexhull.xy(matrix(c(x=xy[,1], y=xy[,2]), ncol=2))
  clipp <- cvex$bdry[[1]]

  canopy  <- ifelse(thresh.var >= thresh.val, 1, 0)
  mv      <- deldir::deldir(x=xy[,1], y=xy[,2], z=canopy, rw=NULL, eps=1e-09, digits=6, plotit=FALSE, suppressMsge=TRUE)
  if(!is.null(mv)) {
    result <- 1-(sum(mv$summary$dir.area*mv$summary$z) / mv$dir.area)
  } else result <- 0

  pol.res <- deg2rad(pol.deg)
  azi.res <- deg2rad(azi.deg)

  start      <- pi/2
  rp.type    <- 's'
  point.symbols <- 19
  point.col  <- col
  clockwise  <- TRUE
  show.grid.labels <- FALSE
  cex        <- 1
  nsets      <- 1
  maxlength  <- diff(range(xy))/2
  radial.lim <- c(0, maxlength)
  ngpos      <- (pi/2)/pol.res
  grid.pos   <- seq((pi/2)*(1/ngpos), pi/2, length.out=ngpos)
  nlpos      <- (2*pi)/azi.res
  label.pos  <- seq(0, pi*(2-2/nlpos), length.out=nlpos)

  if (plots.save == TRUE) {
    if(!exists("LASfolder")) stop("You can only save plots when a LAS file path is used")
    pdf(file.path(LASfolder, paste(LASname,'_closure.pdf',sep='')), width=8, height=8, units='in', res=300)
    par(mfrow=c(2,2), mar=c(2,2,3,2), pty='s', xpd=TRUE)

    plot(LASord[,1],LASord[,2], pch=point.symbols, col=LAScol, bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir')

    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar')
    points(xy[,1], xy[,2], pch=point.symbols, col=point.col)
    radial.grid.hemi(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE, model=model, r = NA)

    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar Delaunay')
    points(xy[,1], xy[,2], pch=point.symbols, col=point.col)
    geometry::trimesh(md, xy, add=TRUE)
    radial.grid.hemi(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE, model=model, r = NA)

    fillcol <- ifelse(thresh.var >= thresh.val, col, NA)
    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar Voronoi')
    deldir::plot.tile.list(deldir::tile.list(mv), verbose=FALSE, close=TRUE, pch=NA, fillcol=fillcol, col.pts=NA, border='white',
                   showpoints=FALSE, add=TRUE, asp=1, clipp=clipp, alpha=0.5)
    radial.grid.hemi(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE, model=model, r = NA)
    dev.off()
  }

  if(plots == TRUE & plots.each == FALSE) {
    par(mfrow=c(2,2), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LASord[,1],LASord[,2], pch=point.symbols, col=LAScol, bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir')

    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar')
    points(xy[,1], xy[,2], pch=point.symbols, col=point.col)
    radial.grid.hemi(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE, model=model, r = NA)

    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar Delaunay')
    points(xy[,1], xy[,2], pch=point.symbols, col=point.col)
    geometry::trimesh(md, xy, add=TRUE)
    radial.grid.hemi(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE, model=model, r = NA)

    fillcol <- ifelse(thresh.var >= thresh.val, col, NA)
    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar Voronoi')
    deldir::plot.tile.list(deldir::tile.list(mv), verbose=FALSE, close=TRUE, pch=NA, fillcol=fillcol, col.pts=NA, border='white',
                   showpoints=FALSE, add=TRUE, asp=1, clipp=clipp, alpha=0.5)
    radial.grid.hemi(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE, model=model, r = NA)
    par(mfrow=c(1,1))
  }

  if(plots == TRUE & plots.each == TRUE) {
    par(mfrow=c(1,1))
    plot(LASord[,1],LASord[,2], pch=point.symbols, col=LAScol, bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir')

    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar')
    points(xy[,1], xy[,2], pch=point.symbols, col=point.col)
    radial.grid.hemi(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                     clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                     grid.bg = 'transparent', show.radial.grid = TRUE, model=model, r = NA)

    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar Delaunay')
    points(xy[,1], xy[,2], pch=point.symbols, col=point.col)
    geometry::trimesh(md, xy, add=TRUE)
    radial.grid.hemi(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                     clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                     grid.bg = 'transparent', show.radial.grid = TRUE, model=model, r = NA)

    fillcol <- ifelse(thresh.var >= thresh.val, col, NA)
    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar Voronoi')
    deldir::plot.tile.list(deldir::tile.list(mv), verbose=FALSE, close=TRUE, pch=NA, fillcol=fillcol, col.pts=NA, border='white',
                           showpoints=FALSE, add=TRUE, asp=1, clipp=clipp, alpha=0.5)
    radial.grid.hemi(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                     clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                     grid.bg = 'transparent', show.radial.grid = TRUE, model=model, r = NA)
  }
  return(result)
}
