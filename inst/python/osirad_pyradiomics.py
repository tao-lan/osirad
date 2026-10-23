# Called by osirad::pyradiomics_features() through 'reticulate'.
# Needs the Python modules numpy, SimpleITK and radiomics (pyradiomics).

import logging
import numpy as np
import SimpleITK as sitk
import radiomics
from radiomics import featureextractor


def extract(pet_dir, suv_factor, x_px, y_px, z_mm, suv, settings, feature_classes):
    radiomics.setVerbosity(logging.ERROR)

    # PET volume in SUV. SimpleITK sorts the slices by position and applies the
    # RescaleSlope of each slice, so the image is in Bq/ml.
    reader = sitk.ImageSeriesReader()
    files = reader.GetGDCMSeriesFileNames(pet_dir)
    reader.SetFileNames(files)
    reader.MetaDataDictionaryArrayUpdateOn()
    bqml = reader.Execute()
    suv_img = sitk.Cast(bqml, sitk.sitkFloat64) * float(suv_factor)

    # Mask (1 = ROI) on the PET grid. The slice of each voxel is the one whose
    # z position (ImagePositionPatient) is nearest to z_mm: the slice spacing
    # is not always constant.
    x_px = np.atleast_1d(np.asarray(x_px, dtype=int))
    y_px = np.atleast_1d(np.asarray(y_px, dtype=int))
    z_mm = np.atleast_1d(np.asarray(z_mm, dtype=float))
    suv = np.atleast_1d(np.asarray(suv, dtype=float))
    slice_z = np.array([float(reader.GetMetaData(i, "0020|0032").split("\\")[2])
                        for i in range(len(files))])
    k = np.abs(z_mm[:, None] - slice_z[None, :]).argmin(axis=1)
    mask_arr = np.zeros(sitk.GetArrayViewFromImage(suv_img).shape, dtype=np.uint8)  # [z, y, x]
    mask_arr[k, y_px, x_px] = 1
    mask = sitk.GetImageFromArray(mask_arr)
    mask.CopyInformation(suv_img)

    # Check values: the mask voxels must have the SUV computed in R
    suv_py = sitk.GetArrayViewFromImage(suv_img)[k, y_px, x_px]

    params = {
        "setting": dict(settings, label=1),
        "imageType": {"Original": {}},
        "featureClass": {c: [] for c in np.atleast_1d(feature_classes)},  # [] = all features
    }
    extractor = featureextractor.RadiomicsFeatureExtractor(params)
    result = extractor.execute(suv_img, mask)

    names = list(result.keys())
    is_diag = [n.startswith("diagnostics_") for n in names]
    return {
        "feature": [n for n, d in zip(names, is_diag) if not d],
        "value": [float(result[n]) for n, d in zip(names, is_diag) if not d],
        "diagnostic": [n for n, d in zip(names, is_diag) if d],
        "diagnostic_value": [str(result[n]) for n, d in zip(names, is_diag) if d],
        "n_mask": int(mask_arr.sum()),
        "max_z_distance": float(np.abs(slice_z[k] - z_mm).max()),
        "max_suv_difference": float(np.abs(suv_py - suv).max()),
    }
