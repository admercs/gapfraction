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

The 'gapfraction' package is designed for processing airborne and terrestrial laser scanning (ALS and TLS) light-detection-and-ranging (LiDAR) data of forests. The package is based on the final chapter of my doctoral dissertation at University of British Columbia footnotes^[Erickson, A. (2016) Forecasting Brown Bear (Ursus Arctos) Habitat Through the Integration of Remote Sensing, a Process-Based Tree Establishment Model, and a Forest Landscape Model. University of British Columbia.]. The package implements my fast-pit-free canopy height model (CHM) algorithm based on Khosravipour et al. (2013) footnotes^[Khosravipour et al. (2013) Development of an algorithm to generate a LiDAR pit-free canopy height model. http://www.riegl.com/uploads/tx_pxpriegldownloads/khosravipour_SilviLaser2013.pdf]. 

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

You can write math expressions, e.g. $Y = X\beta + \epsilon$, footnotes^[Something cool!], and tables, e.g. using `knitr::kable()`.

```{r, echo=FALSE, results='asis'}
knitr::kable(head(mtcars, 10))
```

Also a quote using `>`:

> "He who gives up [code] safety for [code] speed deserves neither."
([via](https://twitter.com/hadleywickham/status/504368538874703872))
