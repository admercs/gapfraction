gapfraction.par <- function(models=models, threshs=threshs, las.files=las.files) {
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
