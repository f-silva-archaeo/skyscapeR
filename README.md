[![cran version](http://www.r-pkg.org/badges/version/skyscapeR)](https://cran.rstudio.com/web/packages/skyscapeR) 
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/skyscapeR?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/skyscapeR?color=82b4e8)](https://github.com/metacran/cranlogs.app)

# skyscapeR
_skyscapeR_ is a open source R package for data reduction, visualization and analysis in skyscape archaeology, archaeoastronomy and cultural astronomy. It is intended to become a fully-fledged, transparent and peer-reviewed package offering a robust set of quantitative methods while retaining simplicity of use.

For information on how to use _skyscapeR_ download [the official vignette](https://github.com/f-silva-archaeo/skyscapeR/blob/master/doc/vignette.html). This is slightly out of date though, so watch this space.


## Release Notes
### v0.3.0.9000 notes
This version begins the _plotly_ implementation process for the plotting functions, which has other side-effects of changing other functions. It's not considerable stable yet, though it run hrough the vignette without a hitch. Some new functions, however, are not finalised yet (avoid anything in the plotting_plotly.R file for now).


### v0.2.9 notes
This version has abandoned _astrolibR_ ephemeris completely and, in its stead, has begun using the Swiss Ephemeris version of the JPL DE431 dataset. This is implemented via package _swephR_ which is still in development and not on CRAN. You can find the offical _swephR_ GitHub [here](https://github.com/rstub/swephR).
 
However, as of v0.3.0, when you install _skyscapeR_ from GitHub it will automatically install _swephR_. Just do:

```r
if(!requireNamespace("devtools", quietly = TRUE)) { install.packages("devtools") }
devtools::install_github('f-silva-archaeo/skyscapeR')
```

### v0.2.2 CRAN release
The latest release version (v0.2.2) is available on CRAN. This Git contains the latest development version which has several bug fixes and additional tools (some of which might not have been fully tested).
