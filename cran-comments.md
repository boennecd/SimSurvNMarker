## Test environments
* Ubuntu 18.04 LTS with gcc 8.3.0
  R version 3.6.3
* win-builder (devel, release, and oldrelease)
* `rhub::check_for_cran()`
  
## R CMD check results
There were no ERRORs, or notes.

There is a WARNING stating that the package have been archived.

I have removed the `assert`s that would be tricked if `-DNDEBUG` is not defined
during compilation.
