lai <- function(LASpath=NA, pol.deg=5, azi.deg=45, reprojection=NA, silent=FALSE, plots=FALSE) {

  if (is.na(LASpath)) stop('Please input a full file path to the LAS file')

  #require(sp)
  #require(rLiDAR)
  #require(plotrix)

  myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue=255)
  }

  crt2sph <- function(m) {
    x     <- m[,1]
    y     <- m[,2]
    z     <- m[,3]
    r     <- sqrt(x^2 + y^2 + z^2)
    theta <- acos(z/r)
    phi   <- atan2(y=y, x=x)
    nrows <- length(x)
    out   <- matrix(c(theta,phi,r), nrow=nrows, ncol=3, dimnames=list(c(1:nrows),c('theta','phi','r')))
    return(out)
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

  x    <- (LAS[,1]-min(LAS[,1])) - (diff(range(LAS[,1]))/2)
  y    <- (LAS[,2]-min(LAS[,2])) - (diff(range(LAS[,2]))/2)
  z    <-  LAS[,3]

  nrows <- length(x)
  m     <- matrix(c(x,y,z), nrow=nrows, ncol=3, dimnames=list(c(1:nrows),c('X','Y','Z')))
  a     <- crt2sph(m=m)

  deg2rad <- function(x) return(x*(pi/180))
  rad2deg <- function(x) return(x*(180/pi))

  pol.res <- deg2rad(pol.deg)
  azi.res <- deg2rad(azi.deg)

  theta <- a[,1]
  phi   <- a[,2]
  r     <- a[,3]

  R      <- max(r)
  n.pts  <- length(a[,1])
  n.pol  <- (pi/2)/pol.res
  n.azi  <- (pi*2)/azi.res
  n.win  <- n.pol * n.azi

  pol <- c(0, seq(from=pol.res, to=(pi/2), by=pol.res))
  azi <- c(0, seq(from=azi.res, to=(pi*2), by=azi.res))

  pt.density <- n.pts / (pi*R^2)
  gapfrac.n  <- matrix(nrow=n.pol, ncol=n.azi)

  for(i in 1:(length(pol)-1)) {
    for(j in 1:(length(azi)-1)) {
      gapfrac.n[i,j] <- (length(r[which(findInterval(theta,c(pol[i],pol[i+1]))==1 & findInterval(phi,c(azi[j],azi[j+1]))==1)]) / pt.density)
    }
  }
  message('Estimated ground plane points: ',round((1-(sum(gapfrac.n)/(n.pts/pt.density)))*100,2),'%')
  message('Canopy-to-total-return ratio: ',round((sum(gapfrac.n)/(n.pts/pt.density))*100,2),'%')

  quad.area <- matrix(nrow=n.pol, ncol=n.azi)

  for(i in 1:(length(pol)-1)) {
    for(j in 1:(length(azi)-1)) {
      quad.area[i,j] <- ( (R^2) * (sin(pol[i+1])-sin(pol[i])) * (azi[j+1]-azi[j]) ) / ((4*pi*R^2)/2)
    }
  }
  if(sum(quad.area) != 1) message('Caution: Quadrangles do not sum to calculated hemisphere area')

  gapfrac <- quad.area / gapfrac.n
  gapfrac[!is.finite(gapfrac)] <- 0
  message('Total gap fraction: ',round(sum(gapfrac)*100,2),'%')

  e.lai <- c()
  for(i in (1:length(pol)-1)[-1]) {
    log.gap  <- log(sum(gapfrac[i,]))
    e.lai[i] <- -2 * ifelse(is.finite(log.gap),log.gap,0) * cos(pol[i+1]) * sin(pol[i+1])
  }
  message('Mean effective LAI: ',round(mean(e.lai),2))

  aci <- c()
  for(i in (1:length(pol)-1)[-1]) {
    mlog <- ifelse(is.finite(log(gapfrac[i,])), log(gapfrac[i,]), 0)
    aci[i] <- 1-0.5*((var(gapfrac[i,])/mean(gapfrac[i,])^2)*cos(pol[i+1])*sin(pol[i+1]) /
        (-mean(mlog)*cos(pol[i+1])*sin(pol[i+1])))
  }
  aci[!is.finite(aci)] <- 0

  pol.sum <- apply(gapfrac, 1, sum)
  names(pol.sum) <- rad2deg(pol[-1])

  azi.sum <- apply(gapfrac, 2, sum)
  names(azi.sum) <- rad2deg(azi[-1])

  result <- c(gapfraction=sum(gapfrac), e.lai=mean(e.lai), aci=mean(aci))

  if (plots==TRUE) {

    jpeg(file.path(LASfolder, paste(LASname,'_gf_lai_aci.jpg',sep='')), width=8, height=8, units='in', res=300, quality=100)
    par(mfrow=c(2,2), mar=c(2,2,3,2), pty='s', xpd=TRUE)

    plot(pol.sum, type='b', xaxt='n', xlab='Zenith Angle', ylab='Gap Fraction', lwd=2, main='Gap Fraction by Zenith Angle')
    axis(1, at=1:length(pol.sum), labels=names(pol.sum))

    plot(azi.sum, type='b', xaxt='n', xlab='Azimuth Angle', ylab='Gap Fraction', lwd=2, main='Gap Fraction by Azimuth Angle')
    axis(1, at=1:length(azi.sum), labels=names(azi.sum))

    plot(e.lai, type='b', xaxt='n', xlab='Zenith Angle', ylab='Effective LAI', lwd=2, main='Effective LAI by Zenith Angle')
    axis(1, at=1:length(pol.sum), labels=names(pol.sum))

    plot(aci, type='b', xaxt='n', xlab='Zenith Angle', ylab='Apparent Clumping Index', lwd=2, main='Apparent Clumping Index by Zenith Angle')
    axis(1, at=1:length(pol.sum), labels=names(pol.sum))

    dev.off()
  }
  if(silent==FALSE) {
    par(mfrow=c(2,2), mar=c(2,2,3,2), pty='s', xpd=TRUE)

    plot(pol.sum, type='b', xaxt='n', xlab='Zenith Angle', ylab='Gap Fraction', lwd=2, main='Gap Fraction by Zenith Angle')
    axis(1, at=1:length(pol.sum), labels=names(pol.sum))

    plot(azi.sum, type='b', xaxt='n', xlab='Azimuth Angle', ylab='Gap Fraction', lwd=2, main='Gap Fraction by Azimuth Angle')
    axis(1, at=1:length(azi.sum), labels=names(azi.sum))

    plot(e.lai, type='b', xaxt='n', xlab='Zenith Angle', ylab='Effective LAI', lwd=2, main='Effective LAI by Zenith Angle')
    axis(1, at=1:length(pol.sum), labels=names(pol.sum))

    plot(aci, type='b', xaxt='n', xlab='Zenith Angle', ylab='Apparent Clumping Index', lwd=2, main='Apparent Clumping Index by Zenith Angle')
    axis(1, at=1:length(pol.sum), labels=names(pol.sum))

    par(mfrow=c(1,1))
  }
  return(result)
}
