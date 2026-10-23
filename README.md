# osirad

**osirad** finds the exact voxels of regions of interest (ROIs) drawn in the
'OsiriX' DICOM viewer on PET images, and passes them to 'PyRadiomics' for
radiomics features.

'OsiriX' exports ROIs as an XML file with the polygon corner points and the
value of every voxel inside each ROI, but not the voxel positions. osirad
converts the PET DICOM slices to SUV and finds the pixels whose SUV equals the
values in the XML file.

## Installation

Development version:

```r
# install.packages("remotes")
remotes::install_github("tao-lan/osirad")
```

## Example

```r
library(osirad)
xml_file = system.file("extdata", "phantom", "ROI.xml", package = "osirad")
pet_dir = system.file("extdata", "phantom", "pet", package = "osirad")

x = extract_roi_voxels(xml_file, pet_dir)
x                     # match report
head(x$voxels)        # one row per ROI voxel: position (px, mm) and SUV
plot_roi_slices(x)    # ROI voxels on the PET slices

# radiomics features ('Python' with pyradiomics needed)
rad = pyradiomics_features(x, pet_dir, python = "/path/to/python")
```

See `vignette("osirad")` for the full workflow, including the 'Python' setup.
