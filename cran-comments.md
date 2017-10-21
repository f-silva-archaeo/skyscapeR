## Resubmission
This is a resubmission. In this version I have:

* Fixed mis-spelled words in DESCRIPTION and throughout.

* Replaced 'An R package' in description by 'A toolset'
  
I was also asked about:

* Rather, is there some reference about the method you can 
add in the Description field in the form Authors (year) <doi:.....>? 

The package contains a variety of tools/methods, which are referenced 
in the appropriate .Rd files. I have now added extra references where
appropriate.

* Package has a FOSS license but eventually depends on the following 
package which may restrict use: palinsol 

I've implemented the necessary changes so that this package no longer 
depends on package palinsol.



## Test environments
* local OS X install, R 3.4.2
* ubuntu 12.04 (on travis-ci), R 3.4.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

* There are no downstream dependencies.

* All revdep maintainers were notified of the release on 20 Oct 2017.
