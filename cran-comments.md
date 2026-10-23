## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Test environments

* local macOS (aarch64), R 4.5.2
* win-builder (R-devel)  -- to be run
* macOS builder (R-release) -- to be run

## Notes for the reviewer

* 'Python' with 'PyRadiomics' is optional. It is used only by
  pyradiomics_features(); its example, test and vignette chunk run only when
  the environment variable OSIRAD_PYTHON gives a 'Python' path, so no
  'Python' is used on CRAN.
* The example data in inst/extdata are a synthetic phantom (no patient data).
