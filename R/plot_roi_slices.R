#' Plot PET slices with the ROI voxels
#'
#' Draws each ROI slice in grey scale (SUV) with the ROI voxels as red squares,
#' four slices per page. Row 1 of the image is at the top, as in 'OsiriX'.
#'
#' @param x Result of [extract_roi_voxels()].
#' @param image_index Slices to draw (`ImageIndex`). `NULL` (default) draws
#'   all ROI slices.
#' @param ncol,nrow Number of panels per page.
#'
#' @return `x`, invisibly. Called for the plot.
#'
#' @examples
#' xml_file = system.file("extdata", "phantom", "ROI.xml", package = "osirad")
#' pet_dir = system.file("extdata", "phantom", "pet", package = "osirad")
#' x = extract_roi_voxels(xml_file, pet_dir)
#' plot_roi_slices(x, image_index = 4)
#' @export
plot_roi_slices = function(x, image_index = NULL, ncol = 2, nrow = 2) {
  ks = sort(as.integer(names(x$slices)))
  if(!is.null(image_index)) ks = ks[ks %in% image_index]
  n = min(length(ks), ncol * nrow)
  old = par(mfrow = c(ceiling(n / ncol), min(n, ncol)), mar = rep(0.15, 4), oma = rep(0, 4))
  on.exit(par(old))
  for(k in ks) {
    s = x$slices[[as.character(k)]]$suv
    n_row = nrow(s)
    n_col = ncol(s)
    image(x = 0:(n_col - 1), y = 0:(n_row - 1), z = t(s[n_row:1, ]),
          col = gray(0:255 / 255), asp = 1, axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
    v = x$voxels[x$voxels$image_index == k, ]
    y_plot = n_row - 1 - v$y_px
    rect(v$x_px - 0.5, y_plot - 0.5, v$x_px + 0.5, y_plot + 0.5, border = "red", lwd = 0.5)
    text(0, n_row - 1, adj = c(0, 1), col = "yellow", cex = 0.8,
         paste("ImageIndex", k, "| ROI", paste(unique(v$roi_id), collapse = ", "),
               "|", nrow(v), "voxels"))
  }
  invisible(x)
}
