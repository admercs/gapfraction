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

The `gapfraction` package is designed for processing airborne and terrestrial laser scanning (ALS and TLS) light-detection-and-ranging (LiDAR) data of forests. The gapfraction package is based on the final chapter of my doctoral dissertation at University of British Columbia^[Erickson, A. (2016) Forecasting Brown Bear (Ursus Arctos) Habitat Through the Integration of Remote Sensing, a Process-based Tree Establishment Model, and a Forest Landscape Model. University of British Columbia.]. The package is designed to be used with LiDAR data pre-processed with LAStools^[LAStools: http://rapidlasso.com/lastools/] or USDA Fusion^[USDA Fusion: http://forsys.cfr.washington.edu/fusion/fusionlatest.html].

The main input to functions in the `gapfraction` package are normalized point clouds, typically LiDAR plots corresponding to field plots. For LiDAR data without ground point classifications, normalized point clouds can be produced either with two LAStools command line functions, `lasground` and `lasheight`, or with three functions in Fusion, `GroundFilter`, `GridSurfaceCreate`, and `CanopyModel`. If the ground points are already classified then you only need to use the `lasheight` function of LAStools, while the process for Fusion still requires the same three functions. Hence, I recommend that you use LAStools. Besides, the LAStools ground point classification algorithm is superior to that of Fusion, producing superior height normalization results. That because Fusion uses the Kraus and Pfeifer (1998) algorithm^[Kraus and Pfeifer (1998) Determination of terrain models in wooded areas with airborne laser scanner data. http://www.sciencedirect.com/science/article/pii/S0924271698000094], while LAStools implements an optimized version of Axelsson (1999)^[Axelsson (1999) Processing of laser scanner data—algorithms and applications. http://www.sciencedirect.com/science/article/pii/S0924271699000088]. For more information, read Maguya, Junttila, Kauranne (2014)^[Maguya, Junttila, and Kauranne (2014) http://www.mdpi.com/2072-4292/6/7/6524].

The package implements my new fast-pit-free canopy height model (CHM) algorithm based on Khosravipour et al. (2013)^[Khosravipour et al. (2013) Development of an algorithm to generate a LiDAR pit-free canopy height model. http://www.riegl.com/uploads/tx_pxpriegldownloads/khosravipour_SilviLaser2013.pdf], two new LiDAR metrics of canopy gap fraction ($P_o$) and angular canopy closure ($ACC$), several recent individual tree crown (ITC) detection methods, canopy distance and direction metrics, effective leaf area index ($LAI_e$) and apparent clumping index ($ACI$) estimation methods, as well as four mathematical fisheye (hemispherical) lens models: equi-angular, equi-distant, stereographic, and orthographic.

- Never uses retina figures
- Has a smaller default figure size
- Uses a custom CSS stylesheet instead of the default Twitter Bootstrap style

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
