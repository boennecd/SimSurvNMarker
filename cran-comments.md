## Test environments
* Ubuntu 18.04 LTS with gcc 8.3.0
  R version 3.6.3
* win-builder (devel, release, and oldrelease)
* `rhub::check_for_cran()`
  
## R CMD check results
There were no ERRORs, or notes.

There is a WARNING stating that the package have been archived.

I have removed the C++14 depends. This should solve the issue with the 
compiler used by the old R version.
