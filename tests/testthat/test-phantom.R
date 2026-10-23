xml_file = system.file("extdata", "phantom", "ROI.xml", package = "osirad")
pet_dir = system.file("extdata", "phantom", "pet", package = "osirad")
truth = read.csv(test_path("phantom_truth.csv"))

test_that("read_osirix_rois() reads all ROIs and their values", {
  rois = read_osirix_rois(xml_file)
  expect_length(rois, 8)
  expect_setequal(vapply(rois, `[[`, "", "roi_name"), c("Tumour", "Node"))
  expect_equal(sum(vapply(rois, `[[`, 0, "n_voxels_osirix")), nrow(truth))
  expect_equal(colnames(rois[[1]]$points_px), c("x", "y"))
  expect_equal(colnames(rois[[1]]$points_mm), c("x", "y", "z"))
  expect_length(read_osirix_rois(xml_file, roi_name = "Node"), 2)
})

test_that("suv_factor_from_header() gives the phantom SUV factor", {
  hdr = oro.dicom::readDICOMFile(list.files(pet_dir, full.names = TRUE)[1], pixelData = FALSE)$hdr
  # 70 kg, 3.7e8 Bq, half-life 6588 s, 3600 s from injection to series start
  expect_equal(suv_factor_from_header(hdr), 70000 / (3.7e8 * exp(-log(2) * 3600 / 6588)))
})

test_that("extract_roi_voxels() recovers the true ROI voxels", {
  x = extract_roi_voxels(xml_file, pet_dir)
  expect_s3_class(x, "osirad_roi")
  v = x$voxels
  key = function(d) sort(paste(d$roi_name, d$image_index, d$x_px, d$y_px))
  expect_equal(key(v), key(truth))
  expect_equal(sum(x$match$not_found), 0)
  expect_equal(sum(x$match$tie_by_neighbours), 1)
  # SUV in R equals the 32-bit values stored by OsiriX
  m = merge(v, truth, by = c("roi_name", "image_index", "x_px", "y_px"))
  expect_equal(m$suv.x, m$suv.y, tolerance = 1e-6)
  # mm position: x_mm = -62 + 4 * x_px in the phantom
  expect_equal(v$x_mm, -62 + 4 * v$x_px)
  expect_equal(v$y_mm, -62 + 4 * v$y_px)
})

test_that("roi_name keeps only the chosen ROI", {
  x = extract_roi_voxels(xml_file, pet_dir, roi_name = "Node")
  expect_equal(unique(x$voxels$roi_name), "Node")
  expect_equal(nrow(x$voxels), sum(truth$roi_name == "Node"))
})

test_that("a wrong PET folder is detected", {
  tmp = file.path(tempdir(), "osirad_pet")
  dir.create(tmp, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE))
  file.copy(list.files(pet_dir, full.names = TRUE)[1:3], tmp)
  expect_error(extract_roi_voxels(xml_file, tmp), "slices")
  expect_error(read_osirix_rois(file.path(tmp, "none.xml")), "not found")
})

test_that("print and plot methods run", {
  x = extract_roi_voxels(xml_file, pet_dir)
  expect_output(print(x), "not found: 0")
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot_roi_slices(x))
})
