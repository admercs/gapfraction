#' Parallel Hemispherical-Voronoi Gap Fraction
#'
#' This function allows you to run \code{gf.hv} on multiple plots in parallel using multi-core CPUs with \code{snow}
#' @param las.files String vector or list of LAS files. Defaults to NA.
#' @param models String vector or list of hemispherical lens models to use. Options include equi-distant (\code{"equidist"}), equi-angular (\code{"equiangle"}), stereographic (\code{"stereo"}), orthographic (\code{"ortho"}), or \code{"all"}. Defaults to \code{"all"}.
#' @param thresh.vals String vector or list of height thresholds to use. Defaults to 2.
#' @param thresh.var Specifies the LiDAR metric to use for thresholding canopy points. Options include height, intensity, nreturn, and class. Defaults to height.
#' @param reprojection Proj4 projection string to use for reprojection. Defaults to NA.
#' @param pol.deg Specifies the polar resolution for the radial plot lines. Defaults to 5.
#' @param azi.deg Specifies the azimuthal resolution for the radial plot lines. Defaults to 45.
#' @param col Specifies the LiDAR metric to use to color points of first plot in display. Options include height, intensity, nreturn, and class. Defaults to height.
#' @param silent Boolean switch for the interactive display of plots. Defaults to FALSE.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @keywords gap fraction, parallel, hemispherical Voronoi, tesselation
#' @references Forthcoming
#' @export
#' @return The results of \code{gf.hv.par}
#' @examples
#' gf.hv.par(las.files='C:/plot.las', models=c('equidist','stereo'), threshs=seq(1,4,0.5))
#' gf.hv.par(las.files=las.list, models='all', thresh.vals=1.25, thresh.var='height', reprojection=NA, pol.deg=5, azi.deg=45, col='height', silent=TRUE, plots=FALSE)

gf.hv.par <- function(las.files=NA, models='all', thresh.vals=seq(1,4,0.5), thresh.var='height', reprojection=NA, pol.deg=5, azi.deg=45, col='height', silent=TRUE, plots=FALSE) {

  require(foreach)

  if(any(models=='all')) models <- c('equidist','equiangle','ortho','stereo')

  ncores   <- parallel::detectCores()-1
  clust    <- snow::makeCluster(ncores, type='SOCK')
  doSNOW::registerDoSNOW(clust)

  ntasks   <- (length(las.files)*length(models)*length(thresh.vals))
  pb       <- utils::txtProgressBar(max=ntasks, style=3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts     <- list(progress=progress)

  results  <-
    foreach::foreach(i = models,     .combine='cbind') %:%
    foreach::foreach(j = thresh.vals,.combine='cbind') %:%
    foreach::foreach(k = las.files,  .combine='rbind', .packages=c('gapfraction'), .errorhandling='pass', .options.snow=opts) %dopar% {
      Sys.sleep(1)
      return(gf.hv(las.path=k, model=i, thresh.val=j, thresh.var=thresh.var, reprojection=reprojection, pol.deg=pol.deg, azi.deg=azi.deg, col=col, silent=silent, plots=plots))
    }
  snow::stopCluster(clust)
  return(results)
}
