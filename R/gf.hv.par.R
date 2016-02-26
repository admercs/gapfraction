#' Parallel Hemispherical-Voronoi Gap Fraction
#'
#' This function allows you to run \code{gf.hv} on multiple plots in parallel using multi-core CPUs with \code{snow}
#' @param las.files String vector or list of LAS files. Defaults to NA.
#' @param models String vector or list of hemispherical lens models to use. Options include equi-distant (\code{"equidist"}), equi-angular (\code{"equiangle"}), stereographic (\code{"stereo"}), and orthographic (\code{"ortho"}). Defaults to \code{"equidist"}.
#' @param threshs String vector or list of height thresholds to use. Defaults to 2.
#' @author Adam Erickson, \email{adam.erickson@@ubc.ca}
#' @keywords gap fraction, parallel, hemispherical Voronoi, tesselation
#' @references Forthcoming
#' @export
#' @return The results of \code{gf.hv.par}
#' @examples
#' gf.hv.par(las.files='C:/plot.las', models=c('equidist','stereo'), threshs=seq(1,3))
#' gf.hv.par(las.files=las.list, models='equidist', threshs=2)

gf.hv.par <- function(las.files=NA, models=NA, threshs=NA) {

  require(foreach)

  ncores <- parallel::detectCores()-1
  clust  <- snow::makeCluster(ncores, type='SOCK')
  doSNOW::registerDoSNOW(clust)

  ntasks   <- (length(models)*length(threshs)*length(las.files))
  pb       <- utils::txtProgressBar(max=ntasks, style=3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts     <- list(progress=progress)

  results <-
    foreach::foreach(i = models,   .combine='rbind', .packages=c('gapfraction')) %:%
    foreach::foreach(j = threshs,  .combine='rbind') %:%
    foreach::foreach(k = las.files,.combine='cbind', .options.snow=opts) %dopar% {
      Sys.sleep(0.1)
      gf.hv(las.path=k, model=i, thresh.val=j, thresh.var='height', silent=TRUE, plots=FALSE)
    }
  snow::stopCluster(clust)
  results <- t(results)
  return(results)
}
