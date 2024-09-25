skyscapeR v1.1 (Release date: 2024-09-25)
==============

Changes:

* added dates optional facility to findTargets()
* fixed issues with sun and moon reporting on findtargets()
* Added modelTest()
* Added all stars up to magV 6, i.e. all visible stars
* Added planetary.extremes()
* Fixed issue with reporting of sunrise/set dates in findTargets()
* Fixed issue with az2dec()
* Fixed manual page for randomTest()
* Fixed plotting issue with plotAzimuth()
* Fixed typos in vignette
* Added parameter col to plot.skyscapeR.sigtest()
* coordtrans() can now transform skyscapeR.spd objects
* coordtrans() has been sped up through vectorisation
* Changed randomTest() to use SNOW, now allowing for progress bar
* randomTest() can now perform the test on multiple measurements at the same site/horizon much faster
* Small fixes to riseset()
* Corrected issues with exportHor()
* Changed star() to use internal version of stellar ephemeris and naming (replacing swephR's version)
* Complete rewrite of the engine behind star.phases() to allow for seasonality calculations outside of swephR's time interval
* Fixed issue in star.phases() to do with alt.hor and circumpolar stars
* Added facility to plot multiple star.phases objects
* Added hor2min.dec() and hor2max.dec() to calculate min/max decs for given horizon profile or location
* Added sky colour option to plot.skyscapeR.horizon()
* Fixed issue with jd2time where day was not correct for change of date
* Fixed issue with riseset where date was not displayed accurately when the day changed
