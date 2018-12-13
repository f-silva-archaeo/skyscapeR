[![cran version](http://www.r-pkg.org/badges/version/skyscapeR)](https://cran.rstudio.com/web/packages/skyscapeR) 
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/skyscapeR?)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/skyscapeR?color=82b4e8)](https://github.com/metacran/cranlogs.app)

# skyscapeR
_skyscapeR_ is a open source R package for data reduction, visualization and analysis in skyscape archaeology, archaeoastronomy and cultural astronomy. It is intended to become a fully-fledged, transparent and peer-reviewed package offering a robust set of quantitative methods while retaining simplicity of use.





Currently it includes functions to:

* transform horizontal (Az/Alt) to equatorial (Dec/RA) coordinates;
* automatically convert theodolite or compass field measurements;
* create or download horizon profiles and plot them;
* overlay horizon profiles with visible paths of common celestial objects/events;
* easily construct azimuth polar plots;
* easily construct declination curvigrams;
* run null hypothesis significance test on declination curvigrams;
* estimate dates of stellar phases/events and their seasonality.

Future versions will add further functionality, as well as an easy-to-use Graphical User Interface.

## Release Notes
The latest release version (v0.2.2) is available on CRAN. This Git contains the latest development version which has several bug fixes and additional tools (some of which might not have been fully tested).

For information on how to use _skyscapeR_ see [the official vignette](https://cran.rstudio.com/web/packages/skyscapeR/vignettes/skyscapeR.html).

## v0.2.9 release notes
This version has abandoned _astrolibR_ ephemeris completely and, in its stead, has begun using the Swiss Ephemeris version of the JPL DE431 dataset. This is implemented via package _swephR_ which is still in development and not on CRAN. You can find the offical _swephR_ GitHub [here](https://github.com/rstub/swephR).

It is therefore essential to first install _swephR_ and _swephRdata_ by doing:
```r
if (!requireNamespace("drat", quietly = TRUE)) install.packages("drat")
drat::addRepo("rstub")
install.packages("swephR")
install.packages("swephRdata")
```

After this, the development version of _skyscapeR_ can be installed as follows:
```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github('f-silva-archaeo/skyscapeR')
```



