skyscapeR v1.0.0 (Release date: 2021-10-29)
==============

Changes:

* Fixed all check problems that led to archival of previous version
* Fixed links in README.md
* Replaced curvigram() and sigTest() with az.pdf(), coordtrans(), spd() and randomTest() based on Silva 2020's approach
* Added new celestial events - spatial.equinox(), EFM() - and analytical tools - findTargets(), 
* Added global variables accessible via skyscapeR.vars()
* Fixed issues with star.phases()
* Implemented swephR as new ephemeris 
* Added error messages when trying to model sky outside of swephR range
* Added bernoulli.trial() for discrete approach
* Added a vignette
