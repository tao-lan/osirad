#' Find the voxels of one ROI on one PET slice
#'
#' Finds the pixels whose SUV equals the voxel values that 'OsiriX' stored for
#' one region of interest (ROI). The search area is the range of the ROI
#' corner points plus 2 pixels on each side.
#'
#' @param point_values Voxel values of the ROI (`point_values` from
#'   [read_osirix_rois()]).
#' @param suv SUV matrix of the slice (`suv` from [read_pet_slice()]).
#' @param points_px ROI corner points in pixels (`points_px` from
#'   [read_osirix_rois()]).
#' @param tol Relative tolerance for two values to be equal. The default
#'   allows for 'OsiriX' storing the values as 32-bit floating point numbers.
#'
#' @details For each distinct value that 'OsiriX' lists \eqn{m} times and that
#'   is found at \eqn{c} pixels of the search area:
#'   \itemize{
#'     \item \eqn{c = m}: clear match.
#'     \item \eqn{c < m}: 'OsiriX' listed the same pixel more than once; the
#'       pixels found are kept.
#'     \item \eqn{c > m}: tie (the same SUV at several pixels). The \eqn{m}
#'       pixels with the most clear matches of this ROI among their 8
#'       neighbours are kept (an ROI is one connected region); if this is still
#'       equal, the pixels nearest the centre of the corner points are kept.
#'     \item \eqn{c = 0}: not found.
#'   }
#'   Each ROI is searched with its own values, so several ROIs on one slice
#'   are handled one by one.
#'
#' @return A list with:
#'   \describe{
#'     \item{rc}{Two-column matrix (`row`, `col`) of the ROI pixels in `suv`.}
#'     \item{counts}{Named vector: number of XML values (`xml_values`) and how
#'       many of them were a clear match (`clear`), listed twice
#'       (`listed_twice`), a tie decided by neighbours (`tie_by_neighbours`), a
#'       tie decided by distance (`tie_by_distance`) or not found
#'       (`not_found`).}
#'   }
#'
#' @examples
#' xml_file = system.file("extdata", "phantom", "ROI.xml", package = "osirad")
#' pet_dir = system.file("extdata", "phantom", "pet", package = "osirad")
#' rois = read_osirix_rois(xml_file)
#' r = rois[[1]]
#' f = list.files(pet_dir, full.names = TRUE)
#' hdr = oro.dicom::readDICOMFile(f[1], pixelData = FALSE)$hdr
#' sl = read_pet_slice(f[r$k], suv_factor_from_header(hdr))
#' m = match_roi_voxels(r$point_values, sl$suv, r$points_px)
#' m$counts
#' @export
match_roi_voxels = function(point_values, suv, points_px, tol = 1e-5) {
  xr = range(points_px[, "x"]) + c(-2, 2)
  yr = range(points_px[, "y"]) + c(-2, 2)
  rc = which(!is.na(suv), arr.ind = TRUE)
  rc = rc[rc[, "col"] - 1 >= xr[1] & rc[, "col"] - 1 <= xr[2] &
          rc[, "row"] - 1 >= yr[1] & rc[, "row"] - 1 <= yr[2], , drop = FALSE]
  box_suv = suv[rc]
  values = unique(point_values)
  listed = vapply(values, function(v) sum(point_values == v), numeric(1))
  hits = lapply(values, function(v) which(abs(box_suv - v) <= tol * max(abs(v), 1e-6)))
  found = lengths(hits)

  keep = unlist(hits[found <= listed])
  keep_key = paste(rc[keep, "row"], rc[keep, "col"])
  centre = colMeans(points_px[, c("x", "y"), drop = FALSE])
  n_tie_neighbours = 0
  n_tie_distance = 0
  for(j in which(found > listed)) {
    h = hits[[j]]
    nb = vapply(h, function(q) {
      sum(outer(rc[q, "row"] + -1:1, rc[q, "col"] + -1:1, paste) %in% keep_key)
    }, numeric(1))
    d = (rc[h, "col"] - 1 - centre[1])^2 + (rc[h, "row"] - 1 - centre[2])^2
    ord = order(-nb, d)
    m = listed[j]
    if(nb[ord[m]] > nb[ord[m + 1]]) {
      n_tie_neighbours = n_tie_neighbours + m
    } else {
      n_tie_distance = n_tie_distance + m
    }
    keep = c(keep, h[ord[1:m]])
  }

  list(
    rc = rc[keep, , drop = FALSE],
    counts = c(
      xml_values = length(point_values),
      clear = sum(listed[found == listed]),
      listed_twice = sum(listed[found > 0 & found < listed]),
      tie_by_neighbours = n_tie_neighbours,
      tie_by_distance = n_tie_distance,
      not_found = sum(listed[found == 0])
    )
  )
}
