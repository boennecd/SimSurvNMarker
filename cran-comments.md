## Test environments
* Ubuntu 18.04 LTS with gcc 8.3.0
  R version 3.6.3
* Ubuntu 18.04 LTS with gcc 8.3.0
  R version 3.6.3 with `--use-valgrind`  
* Ubuntu 16.04 LTS (on travis-ci)
  R version 4.0.0
* Ubuntu 18.04 LTS with clang 6.0.0 with ASAN and 
  UBSAN checks
  R devel (2020-08-19 r79050)
* win-builder (devel and release)
* `rhub::check_for_cran()`
  
## R CMD check results
There were no ERRORs, WARNINGs, or notes.

## Resubmission
This is a resubmission. In this version I have:

* removed the re-setting of the seed which changed `.GlobalEnv`. 
* reset `par()` in the examples.
