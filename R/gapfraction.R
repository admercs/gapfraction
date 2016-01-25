gapfraction <- function(LASpath=NA, model='equidist', reprojection=NA, col='height', thresh.var='height', thresh.val=1.25, silent=TRUE, plots=FALSE) {

  if(is.na(LASpath)) stop('Please input a full file path to the LAS file')

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

  LAS       <- rLiDAR::readLAS(LASpath, short=FALSE)
  LAS       <- LAS[order(LAS[,3], decreasing=FALSE), ]
  LASfolder <- dirname(LASpath)
  LASname   <- strsplit(basename(LASpath), '\\.')[[1]][1]

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

  col <- ifelse(col == 'height',    myColorRamp(colors=c('blue','green','yellow','red'), values=z),
         ifelse(col == 'theta',     myColorRamp(colors=c('blue','green','yellow','red'), values=a[,1]),
         ifelse(col == 'phi',       myColorRamp(colors=c('blue','green','yellow','red'), values=a[,2]),
         ifelse(col == 'r',         myColorRamp(colors=c('blue','green','yellow','red'), values=a[,3]),
         ifelse(col == 'intensity', myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,4]),
         ifelse(col == 'nreturn',   myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,5]),
         ifelse(col == 'class',     myColorRamp(colors=c('blue','green','yellow','red'), values=LAS[,9]),
        'forestgreen')))))))

  LASord <- LAS[order(LAS[,3], decreasing=FALSE),]
  LAScol <- myColorRamp(colors=c('blue','green','yellow','red'), values=LASord[,3])

  if(model == 'stereo')    xy <- crt2str(m)
  if(model == 'ortho')     xy <- crt2ort(m)
  if(model == 'equidist')  xy <- crt2eqd(m)
  if(model == 'equiangle') xy <- crt2esa(m)

  md <- geometry::delaunayn(xy, full=F)

  thresh.var <- rep(thresh.var, nrows)

  thresh.var <- ifelse(thresh.var == 'height',    z,
                ifelse(thresh.var == 'theta',     a[,1],
                ifelse(thresh.var == 'phi',       a[,2],
                ifelse(thresh.var == 'r',         a[,3],
                ifelse(thresh.var == 'intensity', LAS[,4],
                ifelse(thresh.var == 'nreturn',   LAS[,5],
                ifelse(thresh.var == 'class',     LAS[,9],
                z)))))))

  canopy <- ifelse(thresh.var >= thresh.val & LAS[,9] != 2, 1, 0)
  mv     <- deldir::deldir(x=xy[,1], y=xy[,2], z=canopy, rw=NULL, eps=1e-09, plotit=FALSE, suppressMsge=TRUE)
  gf     <- ( sum(mv$summary$dir.area * mv$summary$z) / mv$del.area )

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
  grid.pos   <- c(maxlength*(1/4), maxlength*(1/2), maxlength*(3/4), maxlength)
  angles     <- seq(0, 1.96*pi, by=0.04*pi)
  nlpos      <- 8
  label.pos  <- seq(0, pi*(2-2/nlpos), length.out=nlpos)
  cvex       <- spatstat::convexhull.xy(matrix(c(x=xy[,1], y=xy[,2]), ncol=2))
  clipp      <- cvex$bdry[[1]]

  if (plots == TRUE) {
    jpeg(file.path(LASfolder, paste(LASname,'_closure.jpg',sep='')), width=8, height=8, units='in', res=300, quality=100)
    par(mfrow=c(2,2), mar=c(2,2,3,2), pty='s', xpd=TRUE)

    plot(LASord[,1],LASord[,2], pch=19, col=LAScol, bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir')

    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar')
    points(xy[,1], xy[,2], pch=point.symbols, col=point.col)
    plotrix::radial.grid(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE)

    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar Delaunay')
    points(xy[,1], xy[,2], pch=point.symbols, col=point.col)
    geometry::trimesh(md, xy, add=TRUE)
    plotrix::radial.grid(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE)

    fillcol <- ifelse(thresh.var >= thresh.val & LAS[,9] != 2, col, NA)
    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar Voronoi')
    deldir::plot.tile.list(deldir::tile.list(mv), verbose=FALSE, close=TRUE, pch=NA, fillcol=fillcol, col.pts=NA, border='white',
                   showpoints=FALSE, add=TRUE, asp=1, clipp=clipp, alpha=0.5)
    plotrix::radial.grid(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE)
    dev.off()
  }
  if(silent==FALSE) {
    par(mfrow=c(2,2), mar=c(2,2,3,2), pty='s', xpd=TRUE)
    plot(LASord[,1],LASord[,2], pch=10, col=LAScol, bty='n', xlab='Latitude', ylab='Longitude', main='Cartesian Nadir')

    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar')
    points(xy[,1], xy[,2], pch=point.symbols, col=point.col)
    plotrix::radial.grid(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE)

    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar Delaunay')
    points(xy[,1], xy[,2], pch=point.symbols, col=point.col)
    geometry::trimesh(md, xy, add=TRUE)
    plotrix::radial.grid(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE)

    fillcol <- ifelse(thresh.var >= thresh.val & LAS[,9] != 2, col, NA)
    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type='n', axes=FALSE, xlab=NA, ylab=NA, main='Polar Voronoi')
    deldir::plot.tile.list(deldir::tile.list(mv), verbose=FALSE, close=TRUE, pch=NA, fillcol=fillcol, col.pts=NA, border='white',
                   showpoints=FALSE, add=TRUE, asp=1, clipp=clipp, alpha=0.5)
    plotrix::radial.grid(labels = NA, label.pos = label.pos, radlab = FALSE, radial.lim = radial.lim, start = start,
                clockwise = clockwise, label.prop = 1.1, grid.pos = grid.pos, grid.col = 'gray',
                grid.bg = 'transparent', show.radial.grid = TRUE)
    par(mfrow=c(1,1))
  }
  return(gf)
}
