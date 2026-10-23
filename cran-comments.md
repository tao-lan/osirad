## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* Possibly misspelled words in DESCRIPTION (DICOM, Griethuysen, ROI, ROIs,
  Radiomics, et, al): these are correct. DICOM is the medical image format
  (written out in the Description), ROI is region of interest, radiomics is the
  field of quantitative image features, and Griethuysen et al. is the
  reference for 'PyRadiomics'.

## Test environments

* local macOS (aarch64), R 4.5.2
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)
* win-builder: Windows, R-devel (2026-10-22 r90748 ucrt) -- 1 NOTE (above)
* macOS builder: macOS Sonoma (aarch64), R 4.6.1 Patched -- Status OK

## Notes for the reviewer

* 'Python' with 'PyRadiomics' is optional. It is used only by
  pyradiomics_features(); its example, test and vignette chunk run only when
  the environment variable OSIRAD_PYTHON gives a 'Python' path, so no
  'Python' is used on CRAN.
* The example data in inst/extdata are a synthetic phantom (no patient data).
