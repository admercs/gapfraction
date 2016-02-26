#' Parallel Hemispherical-Voronoi Gap Fraction
#'
#' This function allows you to run gapfraction in parallel on multi-core CPUs
#' @param las.files List of LAS files. Defaults to NA.
#' @param models List of hemispherical lens models to use. Defaults to equidist.
#' @param threshs List of height thresholds to use. Defaults to 2.
#' @keywords gapfraction
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
      gapfraction(las.path=k, model=i, thresh.val=j, thresh.var='height', silent=TRUE, plots=FALSE)
    }
  snow::stopCluster(clust)
  results <- t(results)
  return(results)
}
