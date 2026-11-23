## New submission

This is the first submission of SimSurvNMarker 0.1.4. CRAN archived
version 0.1.3 on 2026-02-15. Version 0.1.4 removes
the C++14 requests from DESCRIPTION and src/Makevars. R now selects the
C++ standard. The examples and test fixtures use the current `dim`
argument to `structure()`.

## Check results

On Debian 13 with GCC 14.2, `R CMD check --as-cran` gave the same result
with R 4.6.2 and R-devel (r90959). Both checks finished with 0 errors,
0 warnings, and 1 NOTE. The incoming feasibility NOTE identifies the
package as a new submission because CRAN archived the previous version.
