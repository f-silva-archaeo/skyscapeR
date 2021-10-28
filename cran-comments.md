## Resubmission
This is a resubmission of a package that was archived on CRAN in 2020 due to check problems.

In this version I have:

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

## Test environments
* local OS X install, R 4.1.1
* Rhub
* win-builder

## R CMD check results
0 errors ✓ | 0 warnings ✓ | 0 notes ✓

