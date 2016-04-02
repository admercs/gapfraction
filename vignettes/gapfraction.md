---
title: "gapfraction: functions for processing LiDAR scans of forests"
author: "Adam Erickson"
date: "April 1, 2016"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{"gapfraction: functions for processing LiDAR scans of forests"}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

The `gapfraction` package is designed for processing airborne laser scanning (ALS) light-detection-and-ranging (LiDAR) data of forests. A version of the package was implemented in the final chapter of my doctoral dissertation at University of British Columbia^[Erickson, A. (2016) Forecasting Brown Bear (Ursus Arctos) Habitat Through the Integration of Remote Sensing, a Process-based Tree Establishment Model, and a Forest Landscape Model. University of British Columbia.]. The package is designed to be used with LiDAR data pre-processed with LAStools^[LAStools: http://rapidlasso.com/lastools/] or USDA Fusion^[USDA Fusion: http://forsys.cfr.washington.edu/fusion/fusionlatest.html]. The main input to functions in the `gapfraction` package are LAS format height-normalized point clouds, typically LiDAR plots corresponding to field plots. The compressed LAZ format is not yet supported, as I rely on the `rLiDAR` `readLAS` function for importing data^[Silva, Crookston, Hudak, and Vierling (2015). rLiDAR: LiDAR Data Processing and Visualization. R package version 0.1. https://CRAN.R-project.org/package=rLiDAR].

## Data Pre-processing

For LiDAR data without ground point classifications, height-normalized point clouds can be produced either with two LAStools command line functions, `lasground` and `lasheight`, or with three functions in Fusion, `GroundFilter`, `GridSurfaceCreate`, and `CanopyModel`. If the ground points are already classified then you only need to use the `lasheight` function of LAStools, while the process for Fusion still requires the same three functions. Hence, I recommend that you use LAStools. Besides, the LAStools ground point classification algorithm is superior to that of Fusion, producing more accurate height-normalized point clouds. That because Fusion uses the Kraus and Pfeifer (1998) algorithm^[Kraus and Pfeifer (1998) Determination of terrain models in wooded areas with airborne laser scanner data. http://www.sciencedirect.com/science/article/pii/S0924271698000094], while LAStools implements an optimized version of the Axelsson (1999) algorithm^[Axelsson (1999) Processing of laser scanner data—algorithms and applications. http://www.sciencedirect.com/science/article/pii/S0924271699000088]. For more information, read Maguya, Junttila, and Kauranne (2014)^[Maguya, Junttila, and Kauranne (2014) http://www.mdpi.com/2072-4292/6/7/6524]. An example application of `lasground` and `lasheight` is provided below. Very simple!

```
lasground -i lidar.las -o lidar_g.las 
lasheight -i lidar_g.las -o lidar_n.las -replace_z 
```

To loop through LAS files stored in a folder, simple use the `system` function included in base R. The syntax follows something like the following pseudo-code:

```
folder <- 'C:/lidar'
files  <- list.files(folder, pattern="\\.las$", full.names=TRUE)

for (i in 1:length(files)) {
  file   <- files[i]
  basenm <- basename(file)
  filenm <- strsplit(basenm,'.',fixed=TRUE)[[1]][1]
  ground <- paste(folder,filenm,'_ground.las',sep='')
  htnorm <- paste(folder,filenm,'_norm.las',sep='')
  system(paste('lasground -i ',file,' -o ',ground, sep=''))
  while (!file.exists(ground)) { Sys.sleep(1) }
  system(paste('lasheight -i ',ground,' -o ',htnorm,' -replace_z', sep=''))
  while (!file.exists(htnorm)) { Sys.sleep(1) }
}
```
What this `for` loop does is read in each LAS file path, extract the name of the file without extension, create the filenames of the ground and height-normalized outputs, execute `lasground`, wait for the output, execute `lasheight` using the ground file as the input, waits for the output, then proceeds to the next iteration. The code would be simple to parallelize using the `foreach` package to speed the operation.

## Functions Included

The `gapfraction` package implements my new fast-pit-free canopy height model (CHM) algorithm based on Khosravipour et al. (2013)^[Khosravipour et al. (2013) Development of an algorithm to generate a LiDAR pit-free canopy height model. http://www.riegl.com/uploads/tx_pxpriegldownloads/khosravipour_SilviLaser2013.pdf], two new LiDAR metrics of canopy gap fraction ($P_o$) and angular canopy closure ($ACC$), several recent individual tree crown (ITC) detection methods, canopy distance and direction metrics, effective leaf area index ($L_e$) and apparent clumping index ($ACI$) estimation methods, as well as four mathematical fisheye (hemispherical) lens models: equi-angular, equi-distant, stereographic, and orthographic. An alphabetical list of functions in the `gapfraction` package is provided below.

- `chm` Simple canopy height model
- `chm.pf` Fast-pit-free canopy height model
- `dd.canopy` Euclidean distance and direction to nearest canopy pixel from plot center
- `dd.crown` Euclidean distance and direction to nearest tree crown from plot centers
- `fc.aci` Above-height cover index of fractional canopy cover
- `fc.bl` Beer-Lambert-Law-modified intensity-return ratio of fractional canopy cover
- `fc.cv` 2-D Cartesian Voronoi fractional canopy cover
- `fc.fci` First-echo cover index of fractional canopy cover
- `fc.fr` Canopy-to-first-return ration of fractional canopy cover
- `fc.ir` Intesity-return ratio of fractional canopy cover
- `fc.p` Canopy-to-total-pixel ratio of fractional canopy cover
- `fc.r` Canopy-to-total-return ratio of fractional canopy cover
- `fc.sci` Solberg's cover index of fractional canopy cover
- `gf.hv` Hemipsherical Voronoi canopy gap fraction
- `gf.hv.par` Parallel hemispherical Voronoi canopy gap fraction with SOCKS
- `gf.laie.aci` Point-density-normalized canopy gap fraction, effective LAI, and ACI
- `itc.mw` Variable-window individual tree crown detection
- `itc.mw.h` Hierarchical variable-window individual tree crown detection
- `itc.wat` Watershed segmentation individual tree crown detection
- `itc.wat.h` Hierarchical watershed segmentation individual tree crown detection
- `lai.e` Ground-to-total-return ratio with a spherical leaf angle distribution
- `lai.n` Contact frequency and fractional canopy cover-based effective LAI
- `radial.grid.hemi` Modified radial.grid function supporting hemispherical lens geometries
- `sun.path` Modified solar position plots of Thomas Steiner

## Vignette Info

Note the various macros within the `vignette` section of the metadata block above. These are required in order to instruct R how to build the vignette. Note that you should change the `title` field and the `\VignetteIndexEntry` to match the title of your vignette.

## Styles

The `html_vignette` template includes a basic CSS theme. To override this theme you can specify your own CSS in the document metadata as follows:

    output: 
      rmarkdown::html_vignette:
        css: mystyles.css

## Figures

The figure sizes have been customised so that you can easily put two images side-by-side. 

```{r, fig.show='hold'}
plot(1:10)
plot(10:1)
```

You can enable figure captions by `fig_caption: yes` in YAML:

    output:
      rmarkdown::html_vignette:
        fig_caption: yes

Then you can use the chunk option `fig.cap = "Your figure caption."` in **knitr**.

## More Examples

You can write math expressions, e.g. $Y = X\beta + \epsilon$, and tables, e.g. using `knitr::kable()`.

```{r, echo=FALSE, results='asis'}
knitr::kable(head(mtcars, 10))
```

Also a quote using `>`:

> "He who gives up [code] safety for [code] speed deserves neither."
([via](https://twitter.com/hadleywickham/status/504368538874703872))
