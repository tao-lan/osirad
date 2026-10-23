#' Extract the voxels of 'OsiriX' ROIs from a PET DICOM series
#'
#' Finds the position and the SUV of every voxel inside the regions of
#' interest (ROIs) of an 'OsiriX' ROI XML file. For each ROI slice, the PET
#' DICOM file is read and converted to SUV, and the pixels whose SUV equals the
#' voxel values stored in the XML file are found with [match_roi_voxels()].
#'
#' @param xml_file Path to the ROI XML file exported from 'OsiriX'.
#' @param pet_dir Folder with the DICOM files of the PET series the ROIs were
#'   drawn on (one file per slice, no other files). The files are sorted by
#'   name; file number `ImageIndex + 1` is the ROI slice.
#' @param roi_name Optional character vector of ROI names to keep. `NULL`
#'   (default) keeps all ROIs.
#'
#' @return An object of class `osirad_roi`, a list with:
#'   \describe{
#'     \item{voxels}{Data frame, one row per ROI voxel: `roi_id`, `roi_name`,
#'       `image_index`, pixel position (`x_px` = column, `y_px` = row, both
#'       starting at 0), position in mm (`x_mm`, `y_mm`, `z_mm`, patient
#'       coordinates) and `suv`.}
#'     \item{match}{Data frame, one row per ROI: the counts of
#'       [match_roi_voxels()] and the number of voxels found.}
#'     \item{suv_factor}{SUV factor from [suv_factor_from_header()].}
#'     \item{slices}{List of SUV slices from [read_pet_slice()], one per ROI
#'       slice, named by `image_index`.}
#'   }
#'
#' @examples
#' xml_file = system.file("extdata", "phantom", "ROI.xml", package = "osirad")
#' pet_dir = system.file("extdata", "phantom", "pet", package = "osirad")
#' x = extract_roi_voxels(xml_file, pet_dir)
#' x
#' head(x$voxels)
#' @export
extract_roi_voxels = function(xml_file, pet_dir, roi_name = NULL) {
  rois = read_osirix_rois(xml_file, roi_name)
  if(length(rois) == 0) {
    stop("No ROI found in ", xml_file)
  }

  pet_files = list.files(pet_dir, full.names = TRUE)
  pet_files = sort(pet_files[!file.info(pet_files)$isdir], method = "radix")
  n_total = rois[[1]]$image_total_num
  if(!is.na(n_total) && length(pet_files) != n_total) {
    stop("pet_dir has ", length(pet_files), " files but the ROIs were drawn on a series of ",
         n_total, " slices")
  }
  hdr = readDICOMFile(pet_files[1], pixelData = FALSE)$hdr
  series_uid = dicom_tag(hdr, "SeriesInstanceUID")
  roi_uid = unique(stats::na.omit(vapply(rois, `[[`, "", "series_uid")))
  if(!is.na(series_uid) && length(roi_uid) > 0 && !all(trimws(roi_uid) == trimws(series_uid))) {
    stop("The ROIs were drawn on another series (SeriesInstanceUID ", roi_uid[1],
         ") than the one in pet_dir (", series_uid, ")")
  }
  suv_factor = suv_factor_from_header(hdr)

  # PET slices with ROIs, in SUV
  slices = list()
  for(k in sort(unique(vapply(rois, `[[`, 0, "k")))) {
    slices[[as.character(k - 1)]] = read_pet_slice(pet_files[k], suv_factor)
  }

  # ROI voxels, ROI by ROI
  voxels = list()
  match = list()
  for(i in seq_along(rois)) {
    r = rois[[i]]
    sl = slices[[as.character(r$image_index)]]
    m = match_roi_voxels(r$point_values, sl$suv, r$points_px)
    rc = m$rc
    voxels[[i]] = data.frame(
      roi_id = rep(i, nrow(rc)), roi_name = rep(r$roi_name, nrow(rc)),
      image_index = rep(r$image_index, nrow(rc)),
      x_px = unname(rc[, "col"]) - 1, y_px = unname(rc[, "row"]) - 1,
      px_to_mm(rc[, "col"] - 1, rc[, "row"] - 1, sl),
      suv = sl$suv[rc]
    )
    match[[i]] = data.frame(roi_id = i, roi_name = r$roi_name, image_index = r$image_index,
                            t(m$counts), voxels = nrow(rc))
  }
  voxels = do.call(rbind, voxels)
  match = do.call(rbind, match)
  rownames(voxels) = NULL

  if(sum(match$not_found) > 0) {
    warning(sum(match$not_found), " of ", sum(match$xml_values),
            " ROI voxel values were not found in the PET slices; ",
            "check that pet_dir is the series the ROIs were drawn on")
  }

  structure(list(voxels = voxels, match = match, suv_factor = suv_factor, slices = slices),
            class = "osirad_roi")
}

#' @export
print.osirad_roi = function(x, ...) {
  m = x$match
  cat("'OsiriX' ROI voxels:", length(unique(m$roi_id)), "ROI(s) on",
      length(unique(m$image_index)), "slice(s),", nrow(x$voxels), "voxels\n")
  cat("ROI names:", paste(unique(m$roi_name), collapse = ", "), "\n")
  cat("SUV factor:", signif(x$suv_factor, 6), "\n")
  cat("XML values:", sum(m$xml_values),
      "| clear match:", sum(m$clear),
      "| listed twice:", sum(m$listed_twice),
      "| tie by neighbours:", sum(m$tie_by_neighbours),
      "| tie by distance:", sum(m$tie_by_distance),
      "| not found:", sum(m$not_found), "\n")
  shared = sum(duplicated(x$voxels[, c("image_index", "x_px", "y_px")]))
  if(shared > 0) cat("Voxels in more than one ROI:", shared, "\n")
  invisible(x)
}
