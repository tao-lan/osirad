#' Radiomics features of the ROI voxels with 'PyRadiomics'
#'
#' Computes radiomics features of the ROI voxels with the 'Python' library
#' 'PyRadiomics' (van Griethuysen et al. 2017), called through 'reticulate'.
#' The PET series is read in 'Python' and converted to SUV with the same SUV
#' factor as in R; the ROI voxels form one mask (label 1). Only the original
#' image is used (no filters).
#'
#' @param x Result of [extract_roi_voxels()].
#' @param pet_dir Folder with the DICOM files of the PET series (the same as in
#'   [extract_roi_voxels()]).
#' @param settings Named list of 'PyRadiomics' settings, see
#'   <https://pyradiomics.readthedocs.io/en/latest/customization.html>. The
#'   default is a fixed bin width of 0.25 SUV and no intensity normalisation;
#'   no resampling.
#' @param feature_classes 'PyRadiomics' feature classes to compute; all
#'   features of each class are computed.
#' @param python Optional path to the 'Python' executable that has
#'   'PyRadiomics' installed. It must be set before 'Python' is first used in
#'   the R session.
#'
#' @return A list with:
#'   \describe{
#'     \item{features}{Data frame, one row per feature: `feature` (name) and
#'       `value`.}
#'     \item{diagnostics}{Data frame, one row per 'PyRadiomics' diagnostics
#'       item (versions, settings, image and mask information): `item` and
#'       `value` (as text).}
#'     \item{check}{Named vector: ROI voxels in `x`, voxels in the 'Python'
#'       mask, largest distance (mm) between a voxel and its slice, and
#'       largest difference between the SUV in 'Python' and in R.}
#'   }
#'
#' @references van Griethuysen JJM, Fedorov A, Parmar C, et al. (2017).
#'   Computational radiomics system to decode the radiographic phenotype.
#'   Cancer Research, 77(21), e104-e107. \doi{10.1158/0008-5472.CAN-17-0339}
#'
#' @examplesIf nzchar(Sys.getenv("OSIRAD_PYTHON"))
#' # Needs 'Python' with 'PyRadiomics'; the path is read from the
#' # environment variable OSIRAD_PYTHON
#' xml_file = system.file("extdata", "phantom", "ROI.xml", package = "osirad")
#' pet_dir = system.file("extdata", "phantom", "pet", package = "osirad")
#' x = extract_roi_voxels(xml_file, pet_dir, roi_name = "Tumour")
#' rad = pyradiomics_features(x, pet_dir, python = Sys.getenv("OSIRAD_PYTHON"))
#' head(rad$features)
#' rad$check
#' @export
pyradiomics_features = function(x, pet_dir,
                                settings = list(binWidth = 0.25, normalize = FALSE),
                                feature_classes = c("shape", "firstorder", "glcm", "glrlm",
                                                    "glszm", "gldm", "ngtdm"),
                                python = NULL) {
  if(!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is needed: install.packages(\"reticulate\")")
  }
  if(!is.null(python)) reticulate::use_python(python, required = TRUE)
  if(!reticulate::py_module_available("radiomics")) {
    stop("'Python' module 'radiomics' (PyRadiomics) not found in ",
         reticulate::py_config()$python, ". See vignette(\"osirad\").")
  }
  py = reticulate::import_from_path("osirad_pyradiomics",
                                    path = system.file("python", package = "osirad"))
  v = x$voxels
  res = py$extract(normalizePath(pet_dir), x$suv_factor, as.integer(v$x_px),
                   as.integer(v$y_px), v$z_mm, v$suv, settings, feature_classes)
  list(
    features = data.frame(feature = unlist(res$feature), value = unlist(res$value)),
    diagnostics = data.frame(item = unlist(res$diagnostic), value = unlist(res$diagnostic_value)),
    check = c(roi_voxels = nrow(v), mask_voxels = res$n_mask,
              max_z_distance = res$max_z_distance, max_suv_difference = res$max_suv_difference)
  )
}
