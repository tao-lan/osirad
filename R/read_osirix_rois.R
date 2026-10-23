# "(x, y)" or "(x, y, z)" strings -> numeric matrix
parse_osirix_points = function(x) {
  x = gsub("[()]", "", x)
  m = do.call(rbind, strsplit(x, ","))
  storage.mode(m) = "numeric"
  colnames(m) = c("x", "y", "z")[seq_len(ncol(m))]
  m
}

# text of the element that follows <key>key</key> in a plist dict
plist_value = function(node, key) {
  xml_text(xml_find_first(node, paste0("./key[.='", key, "']/following-sibling::*[1]")))
}

# text of all elements in the array that follows <key>key</key>
plist_array = function(node, key, type) {
  xml_text(xml_find_all(node, paste0("./key[.='", key, "']/following-sibling::array[1]/", type)))
}

#' Read regions of interest exported from 'OsiriX'
#'
#' Reads the XML file written by the 'OsiriX' menu
#' "Plugins > ROI Tools > Export ROIs" (XML format). The file stores, for each
#' region of interest (ROI) on each slice, the polygon corner points, the value
#' of every voxel inside the ROI and summary statistics. It does not store the
#' voxel positions; [extract_roi_voxels()] finds them.
#'
#' @param xml_file Path to the ROI XML file.
#' @param roi_name Optional character vector. Only the ROIs with these names
#'   are kept. `NULL` (default) keeps all ROIs.
#'
#' @return A list with one element per ROI per slice. Each element is a list
#'   with:
#'   \describe{
#'     \item{image_index}{Slice number in 'OsiriX' (`ImageIndex`, starts at 0).}
#'     \item{k}{`image_index + 1`, the file number in the sorted DICOM series.}
#'     \item{image_width, image_height, image_total_num}{Slice size in pixels
#'       and number of slices in the series.}
#'     \item{roi_name}{ROI name.}
#'     \item{series_uid}{`SeriesInstanceUID` of the image series the ROI was
#'       drawn on.}
#'     \item{points_px, points_mm}{Polygon corner points, in pixels (columns
#'       `x`, `y`) and in mm (columns `x`, `y`, `z`).}
#'     \item{point_values}{Value of every voxel inside the ROI
#'       (`Point_value`).}
#'     \item{n_voxels_osirix}{Length of `point_values`.}
#'     \item{area_cm2, area_pix2, roi_mean, roi_min, roi_max, roi_total}{ROI
#'       statistics computed by 'OsiriX'.}
#'   }
#'
#' @examples
#' xml_file = system.file("extdata", "phantom", "ROI.xml", package = "osirad")
#' rois = read_osirix_rois(xml_file)
#' length(rois)
#' str(rois[[1]])
#' @export
read_osirix_rois = function(xml_file, roi_name = NULL) {
  if(!file.exists(xml_file)) {
    stop("XML file not found: ", xml_file)
  }
  doc = read_xml(xml_file)
  images = xml_find_all(doc, ".//key[.='Images']/following-sibling::array[1]/dict")
  out = list()

  for(img_node in images) {
    image_index = as.integer(plist_value(img_node, "ImageIndex"))
    rois = xml_find_all(img_node, "./key[.='ROIs']/following-sibling::array[1]/dict")

    for(roi_node in rois) {
      name = plist_value(roi_node, "Name")
      if(!is.null(roi_name) && !(name %in% roi_name)) next

      pts_px_txt = plist_array(roi_node, "Point_px", "string")
      if(length(pts_px_txt) < 3) next
      pts_mm_txt = plist_array(roi_node, "Point_mm", "string")
      point_values = as.numeric(plist_array(roi_node, "Point_value", "real"))

      out[[length(out) + 1]] = list(
        image_index = image_index,
        k = image_index + 1,
        image_width = as.integer(plist_value(img_node, "ImageWidth")),
        image_height = as.integer(plist_value(img_node, "ImageHeight")),
        image_total_num = as.integer(plist_value(img_node, "ImageTotalNum")),
        roi_name = name,
        series_uid = plist_value(roi_node, "SeriesInstanceUID"),
        points_px = parse_osirix_points(pts_px_txt),
        points_mm = if(length(pts_mm_txt) > 0) parse_osirix_points(pts_mm_txt) else NULL,
        point_values = point_values,
        n_voxels_osirix = length(point_values),
        area_cm2 = as.numeric(plist_value(roi_node, "AreaCm2")),
        area_pix2 = as.numeric(plist_value(roi_node, "AreaPix2")),
        roi_mean = as.numeric(plist_value(roi_node, "Mean")),
        roi_min = as.numeric(plist_value(roi_node, "Min")),
        roi_max = as.numeric(plist_value(roi_node, "Max")),
        roi_total = as.numeric(plist_value(roi_node, "Total"))
      )
    }
  }

  out
}
