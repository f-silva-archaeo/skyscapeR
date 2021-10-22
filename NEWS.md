skyscapeR v1.0.0 (Release date: 2021-10-28)
==============

Changes:

* Replaced curvigram() and sigTest() with az.pdf(), coordtrans(), spd() and randomTest() based on Silva 2020's approach
* Added new celestial events - spatial.equinox(), EFM() - and analytical tools - findTargets(), 
* Added global variables accessible via skyscapeR.vars()
* Fixed issues with star.phases()
* Added error messages when trying to model sky outside of swephR range
* Added bernoulli.trial() for discrete approach
* Fixed several minor bugs related with upgrade to R > 4.0
