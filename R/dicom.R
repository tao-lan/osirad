# one DICOM header value (oro.dicom header table); NA if the tag is missing
dicom_tag = function(hdr, name) {
  hdr$value[hdr$name == name][1]
}

# multi-valued DICOM tag as numbers
dicom_num = function(hdr, name) {
  as.numeric(strsplit(trimws(dicom_tag(hdr, name)), "[ \\\\]+")[[1]])
}

# DICOM time "hhmmss.ff" -> seconds
dicom_seconds = function(x) {
  x = as.numeric(x)
  h = x %/% 10000
  m = (x %/% 100) %% 100
  3600 * h + 60 * m + (x - 10000 * h - 100 * m)
}

#' SUV factor from a PET DICOM header
#'
#' Computes the factor that converts PET activity concentration (Bq/ml) to
#' body-weight standardised uptake value (SUV):
#' \deqn{SUV = activity \times weight / dose_{decayed},}{SUV = activity * weight / decayed dose,}
#' with the patient weight in g and the injected dose (Bq) decayed from the
#' injection time (`RadiopharmaceuticalStartTime`) to the series start
#' (`SeriesTime`).
#'
#' @param hdr DICOM header table of one slice of the PET series, as returned in
#'   the `hdr` element of [oro.dicom::readDICOMFile()].
#'
#' @return The SUV factor, a single number (g/Bq). Multiply the rescaled pixel
#'   values (Bq/ml) by it to get SUV.
#'
#' @details The image values must be in Bq/ml (`Units` = `BQML`) and decay
#'   corrected to the series start (`DecayCorrection` = `START`); otherwise the
#'   function stops.
#'
#' @examples
#' pet_dir = system.file("extdata", "phantom", "pet", package = "osirad")
#' f = list.files(pet_dir, full.names = TRUE)[1]
#' hdr = oro.dicom::readDICOMFile(f, pixelData = FALSE)$hdr
#' suv_factor_from_header(hdr)
#' @export
suv_factor_from_header = function(hdr) {
  units = dicom_tag(hdr, "Units")
  if(!is.na(units) && trimws(units) != "BQML") {
    stop("PET Units are '", units, "'; the SUV factor needs Bq/ml (BQML)")
  }
  decay = dicom_tag(hdr, "DecayCorrection")
  if(!is.na(decay) && trimws(decay) != "START") {
    stop("PET DecayCorrection is '", decay, "'; the SUV factor needs START")
  }
  weight_g = 1000 * as.numeric(dicom_tag(hdr, "PatientsWeight"))
  dose_bq = as.numeric(dicom_tag(hdr, "RadionuclideTotalDose"))
  half_life = as.numeric(dicom_tag(hdr, "RadionuclideHalfLife"))
  dt = dicom_seconds(dicom_tag(hdr, "SeriesTime")) -
    dicom_seconds(dicom_tag(hdr, "RadiopharmaceuticalStartTime"))
  f = weight_g / (dose_bq * exp(-log(2) * dt / half_life))
  if(!is.finite(f)) {
    stop("SUV factor cannot be computed: PatientWeight, RadionuclideTotalDose, ",
         "RadionuclideHalfLife, SeriesTime or RadiopharmaceuticalStartTime is missing")
  }
  f
}

#' Read one PET slice in SUV
#'
#' Reads one DICOM file of a PET series and converts the pixel values to SUV
#' with the slice's own `RescaleSlope` and `RescaleIntercept` and the SUV
#' factor.
#'
#' @param file Path to the DICOM file.
#' @param suv_factor SUV factor, from [suv_factor_from_header()].
#'
#' @return A list with:
#'   \describe{
#'     \item{suv}{Matrix of SUV, `[row, column]` = `[y_px + 1, x_px + 1]`, in
#'       the DICOM orientation (row 1 at the top of the image).}
#'     \item{ipp}{`ImagePositionPatient`: position (mm) of the first pixel.}
#'     \item{iop}{`ImageOrientationPatient`: direction cosines of the rows and
#'       columns.}
#'     \item{spacing}{`PixelSpacing` (mm): row spacing, column spacing.}
#'   }
#'
#' @examples
#' pet_dir = system.file("extdata", "phantom", "pet", package = "osirad")
#' f = list.files(pet_dir, full.names = TRUE)
#' hdr = oro.dicom::readDICOMFile(f[1], pixelData = FALSE)$hdr
#' sl = read_pet_slice(f[4], suv_factor_from_header(hdr))
#' range(sl$suv)
#' @export
read_pet_slice = function(file, suv_factor) {
  d = readDICOMFile(file, flipud = FALSE)
  slope = as.numeric(dicom_tag(d$hdr, "RescaleSlope"))
  intercept = as.numeric(dicom_tag(d$hdr, "RescaleIntercept"))
  if(is.na(slope)) slope = 1
  if(is.na(intercept)) intercept = 0
  list(
    suv = (d$img * slope + intercept) * suv_factor,
    ipp = dicom_num(d$hdr, "ImagePositionPatient"),
    iop = dicom_num(d$hdr, "ImageOrientationPatient"),
    spacing = dicom_num(d$hdr, "PixelSpacing")
  )
}

# pixel position (x_px, y_px) -> patient coordinates in mm
px_to_mm = function(x, y, sl) {
  mm = outer(x, sl$spacing[2] * sl$iop[1:3]) + outer(y, sl$spacing[1] * sl$iop[4:6])
  mm = sweep(mm, 2, sl$ipp, "+")
  colnames(mm) = c("x_mm", "y_mm", "z_mm")
  mm
}
