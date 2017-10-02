# skyscapeR
_skyscapeR_ is a open source R package for data analysis in skyscape archaeology, archaeoastronomy and cultural astronomy, currently in alpha. It is intended to become a fully-fledged, transparent package offering robustness of quantitative methods yet simplicity of use.


Currently it includes functions to transform horizontal to equatorial coordinates, create or download horizon profiles, plot them and overlay them with visible paths of common celestial objects/events for most prehistoric and all historic epochs. It also includs the ability to easily construct azimuth polar plots and curvigrams of declinations. Future versions will add further data analysis, significance testing and model selection facilities, as well as an easy-to-use Graphical User Interface.


## how to install
The package requires that the latest version of R is installed first. See the [R Project website](https://www.r-project.org/) for details. ALso suggested is the installation of [RStudio](https://www.rstudio.com/).
With R installed, make sure package _devtools_ is also installed by doing:

```
install.packages('devtools')
```

With _devtools_ installed, you can install _skyscapeR_, and associated requirements, by doing:
```
install_github('skyscapeR')
```

Check the man pages for help on individual functions.
