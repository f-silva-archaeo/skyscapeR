# skyscapeR
_skyscapeR_ is a open source R package for data reduction, visualization and analysis in skyscape archaeology, archaeoastronomy and cultural astronomy. It is intended to become a fully-fledged, transparent and peer-reviewed package offering a robust set of quantitative methods while retaining simplicity of use.


Currently it includes functions to transform horizontal (Az/Alt) to equatorial (Dec/RA) coordinates, create or download horizon profiles, plot them and overlay them with visible paths of common celestial objects/events for prehistoric or historic periods. It also includes the ability to easily construct azimuth polar plots and declination curvigrams. Future versions will add data reduction, significance testing and model selection facilities, as well as an easy-to-use Graphical User Interface.


## how to install
The package requires that the latest version of R is installed first. See the [R Project website](https://www.r-project.org/) for details. ALso suggested is the installation of [RStudio](https://www.rstudio.com/). With R installed, make sure package _devtools_ is also installed by running:

```
install.packages('devtools')
```

_skyscapeR_, and associated requirements, can then be installed by running:
```
install_github('skyscapeR')
```

Check the manual pages for help on individual functions.
