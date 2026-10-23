# Synthetic PET phantom for the osirad examples, tests and vignette.
#
# Writes to inst/extdata/phantom/:
#   pet/IM001.dcm ... IM008.dcm  PET DICOM series, 32 x 32 pixels, 8 slices, Bq/ml
#   ROI.xml                      OsiriX-style ROI export (plist), 2 ROIs: "Tumour", "Node"
# and to tests/testthat/:
#   phantom_truth.csv            the true ROI voxels (pixel centre inside the polygon)
#
# The phantom copies the layout of a real OsiriX export:
#   - file order = ImageIndex order, z decreases with the file number
#   - slice spacing 3.27 mm with one 3.35 mm gap (a bed-position join)
#   - each slice has its own RescaleSlope
#   - Point_value = SUV of each ROI voxel as 32-bit float
#   - one background pixel next to the tumour on ImageIndex 2 has the same SUV
#     as one tumour voxel (a tie that match_roi_voxels() must solve)
#
# Run from the package root:
#   python data-raw/make_phantom.py
# Needs numpy and pydicom.

import os
import numpy as np
import pydicom
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

rng = np.random.default_rng(1)
out_dir = "inst/extdata/phantom"
pet_dir = os.path.join(out_dir, "pet")
os.makedirs(pet_dir, exist_ok=True)

# -------------- geometry and SUV factor --------------

nx = ny = 32
nz = 8
spacing = 4.0
x0 = y0 = -62.0
gaps = [3.27, 3.27, 3.27, 3.35, 3.27, 3.27, 3.27]
z = 60.0 - np.concatenate([[0.0], np.cumsum(gaps)])

weight_kg = 70.0
dose_bq = 3.7e8
half_life = 6588.0
start_time = "100000.00"
series_time = "110000"
dt = 3600.0
suv_factor = weight_kg * 1000 / (dose_bq * np.exp(-np.log(2) * dt / half_life))

# -------------- SUV volume [slice, row, col] --------------

rows, cols = np.mgrid[0:ny, 0:nx]
suv = np.empty((nz, ny, nx))
for k in range(nz):
    tumour = 9.0 * np.exp(-((cols - 15.3) ** 2 + (rows - 13.7) ** 2) / (2 * 3.5 ** 2)
                          - (k - 3.5) ** 2 / (2 * 1.3 ** 2))
    node = 5.0 * np.exp(-((cols - 24.2) ** 2 + (rows - 23.5) ** 2) / (2 * 1.5 ** 2)
                        - (k - 4.5) ** 2 / (2 * 0.6 ** 2))
    suv[k] = 1.0 + tumour + node + rng.normal(0, 0.15, (ny, nx))
suv = np.clip(suv, 0.05, None)

# stored integers with one RescaleSlope per slice
stored = np.empty(suv.shape, dtype=np.int16)
slopes = []
for k in range(nz):
    bqml = suv[k] / suv_factor
    slope = float("%.6g" % (bqml.max() / 32000))
    stored[k] = np.round(bqml / slope).astype(np.int16)
    slopes.append(slope)
# SUV as computed from the DICOM files
suv_dicom = stored * np.array(slopes)[:, None, None] * suv_factor

# -------------- ROI polygons --------------

def polygon(cx, cy, r, n=12, phase=0.3):
    a = phase + 2 * np.pi * np.arange(n) / n
    rr = r * (1 + 0.06 * np.sin(3 * a))
    return np.column_stack([cx + rr * np.cos(a), cy + rr * np.sin(a)])

def inside(px, py, poly):
    # ray casting: pixel centre (px, py) inside the closed polygon
    x, y = poly[:, 0], poly[:, 1]
    xj, yj = np.roll(x, 1), np.roll(y, 1)
    res = np.zeros(px.shape, dtype=bool)
    for a, b, c, d in zip(x, y, xj, yj):
        cross = ((b > py) != (d > py)) & (px < (c - a) * (py - b) / (d - b + 1e-300) + a)
        res ^= cross
    return res

rois = []  # (image_index, name, polygon)
for k, r in zip(range(1, 7), [2.2, 3.4, 4.3, 4.3, 3.4, 2.2]):
    rois.append((k, "Tumour", polygon(15.3, 13.7, r)))
for k in (4, 5):
    rois.append((k, "Node", polygon(24.2, 23.5, 1.8, phase=0.1)))

# tie on ImageIndex 2: a background pixel inside the search box gets the stored
# value of the tumour voxel nearest the centre
k_tie = 2
poly = rois[1][2]
mask = inside(cols, rows, poly)
d = (cols - 15.3) ** 2 + (rows - 13.7) ** 2
d[~mask] = np.inf
r_in, c_in = np.unravel_index(np.argmin(d), d.shape)
r_out, c_out = int(np.floor(poly[:, 1].min())) - 1, int(np.floor(poly[:, 0].min())) - 1
assert not mask[r_out, c_out]
stored[k_tie, r_out, c_out] = stored[k_tie, r_in, c_in]
suv_dicom = stored * np.array(slopes)[:, None, None] * suv_factor

# -------------- DICOM files --------------

PET_SOP = "1.2.840.10008.5.1.4.1.1.128"
uid = lambda *s: generate_uid(entropy_srcs=["osirad-phantom", *map(str, s)])
study_uid, series_uid, frame_uid = uid("study"), uid("series"), uid("frame")

for k in range(nz):
    sop_uid = uid("image", k)
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = PET_SOP
    meta.MediaStorageSOPInstanceUID = sop_uid
    meta.TransferSyntaxUID = ExplicitVRLittleEndian
    meta.ImplementationClassUID = pydicom.uid.PYDICOM_IMPLEMENTATION_UID
    fname = os.path.join(pet_dir, "IM%03d.dcm" % (k + 1))
    ds = FileDataset(fname, {}, file_meta=meta, preamble=b"\0" * 128)
    ds.is_little_endian = True
    ds.is_implicit_VR = False

    ds.SOPClassUID = PET_SOP
    ds.SOPInstanceUID = sop_uid
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = series_uid
    ds.FrameOfReferenceUID = frame_uid
    ds.PatientName = "Phantom^Synthetic"
    ds.PatientID = "PHANTOM"
    ds.PatientWeight = "70"
    ds.StudyDate = ds.SeriesDate = ds.AcquisitionDate = "20260101"
    ds.SeriesTime = ds.AcquisitionTime = series_time
    ds.Modality = "PT"
    ds.SeriesDescription = "Synthetic PET phantom"
    ds.SeriesNumber = 1
    ds.InstanceNumber = k + 1
    ds.ImageType = ["ORIGINAL", "PRIMARY"]
    ds.Units = "BQML"
    ds.DecayCorrection = "START"
    ds.CorrectedImage = ["DECY", "ATTN"]

    rp = Dataset()
    rp.Radiopharmaceutical = "FDG -- fluorodeoxyglucose"
    rp.RadiopharmaceuticalStartTime = start_time
    rp.RadionuclideTotalDose = "%.1f" % dose_bq
    rp.RadionuclideHalfLife = "%.1f" % half_life
    rp.RadionuclidePositronFraction = "0.97"
    ds.RadiopharmaceuticalInformationSequence = Sequence([rp])

    ds.ImagePositionPatient = ["%.6f" % x0, "%.6f" % y0, "%.6f" % z[k]]
    ds.ImageOrientationPatient = ["1", "0", "0", "0", "1", "0"]
    ds.SliceLocation = "%.6f" % z[k]
    ds.PixelSpacing = ["%.1f" % spacing, "%.1f" % spacing]
    ds.SliceThickness = "3.27"
    ds.Rows, ds.Columns = ny, nx
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated, ds.BitsStored, ds.HighBit = 16, 16, 15
    ds.PixelRepresentation = 1
    ds.RescaleSlope = "%.6g" % slopes[k]
    ds.RescaleIntercept = "0"
    ds.PixelData = stored[k].astype("<i2").tobytes()
    ds.save_as(fname, write_like_original=False)

# -------------- ROI.xml and true voxels --------------

def real(v):
    return "<real>%s</real>" % repr(float(v))

def px_mm(p, k):
    return x0 + spacing * p[:, 0], y0 + spacing * p[:, 1], np.full(len(p), z[k])

truth = []
images = []
for k in sorted(set(r[0] for r in rois)):
    roi_xml = []
    for j, (kk, name, poly) in enumerate([r for r in rois if r[0] == k]):
        m = inside(cols, rows, poly)
        rr, cc = np.nonzero(m)  # row-major order
        vals = suv_dicom[k][rr, cc].astype(np.float32).astype(np.float64)
        for a, b, v in zip(cc, rr, vals):
            truth.append((name, k, a, b, v))
        xm, ym, zm = px_mm(poly, k)
        area = 0.5 * abs(np.dot(poly[:, 0], np.roll(poly[:, 1], 1)) - np.dot(poly[:, 1], np.roll(poly[:, 0], 1)))
        length = np.sum(np.hypot(*(poly - np.roll(poly, 1, axis=0)).T))
        sd = np.std(vals, ddof=1)
        uid_roi = uid("sop", k)
        roi_xml.append("\n".join([
            "\t\t\t\t<dict>",
            "\t\t\t\t\t<key>AreaCm2</key>", "\t\t\t\t\t" + real(area * spacing ** 2 / 100),
            "\t\t\t\t\t<key>AreaPix2</key>", "\t\t\t\t\t" + real(area),
            "\t\t\t\t\t<key>Center</key>",
            "\t\t\t\t\t<string>(%f, %f, %f)</string>" % (xm.mean(), ym.mean(), z[k]),
            "\t\t\t\t\t<key>Dev</key>", "\t\t\t\t\t" + real(sd),
            "\t\t\t\t\t<key>IndexInImage</key>", "\t\t\t\t\t<integer>%d</integer>" % j,
            "\t\t\t\t\t<key>LengthCm</key>", "\t\t\t\t\t" + real(length * spacing / 10),
            "\t\t\t\t\t<key>LengthPix</key>", "\t\t\t\t\t" + real(length),
            "\t\t\t\t\t<key>Max</key>", "\t\t\t\t\t" + real(vals.max()),
            "\t\t\t\t\t<key>Mean</key>", "\t\t\t\t\t" + real(vals.mean()),
            "\t\t\t\t\t<key>Median</key>", "\t\t\t\t\t" + real(np.median(vals)),
            "\t\t\t\t\t<key>Min</key>", "\t\t\t\t\t" + real(vals.min()),
            "\t\t\t\t\t<key>Name</key>", "\t\t\t\t\t<string>%s</string>" % name,
            "\t\t\t\t\t<key>NumberOfPoints</key>", "\t\t\t\t\t<integer>%d</integer>" % len(poly),
            "\t\t\t\t\t<key>Point_mm</key>", "\t\t\t\t\t<array>",
            *["\t\t\t\t\t\t<string>(%f, %f, %f)</string>" % t for t in zip(xm, ym, zm)],
            "\t\t\t\t\t</array>",
            "\t\t\t\t\t<key>Point_px</key>", "\t\t\t\t\t<array>",
            *["\t\t\t\t\t\t<string>(%f, %f)</string>" % (a, b) for a, b in poly],
            "\t\t\t\t\t</array>",
            "\t\t\t\t\t<key>Point_value</key>", "\t\t\t\t\t<array>",
            *["\t\t\t\t\t\t" + real(v) for v in vals],
            "\t\t\t\t\t</array>",
            "\t\t\t\t\t<key>RadiusHeightCM</key>", "\t\t\t\t\t<real>0.0</real>",
            "\t\t\t\t\t<key>RadiusHeightPIX</key>", "\t\t\t\t\t<real>0.0</real>",
            "\t\t\t\t\t<key>RadiusWidthCM</key>", "\t\t\t\t\t<real>0.0</real>",
            "\t\t\t\t\t<key>RadiusWidthPIX</key>", "\t\t\t\t\t<real>0.0</real>",
            "\t\t\t\t\t<key>SOPInstanceUID</key>", "\t\t\t\t\t<string>%s</string>" % uid_roi,
            "\t\t\t\t\t<key>SeriesInstanceUID</key>", "\t\t\t\t\t<string>%s</string>" % series_uid,
            "\t\t\t\t\t<key>StudyInstanceUID</key>", "\t\t\t\t\t<string>%s</string>" % study_uid,
            "\t\t\t\t\t<key>Total</key>", "\t\t\t\t\t" + real(vals.sum()),
            "\t\t\t\t\t<key>Type</key>", "\t\t\t\t\t<integer>11</integer>",
            "\t\t\t\t</dict>"]))
    images.append("\n".join([
        "\t\t<dict>",
        "\t\t\t<key>ImageHeight</key>", "\t\t\t<integer>%d</integer>" % ny,
        "\t\t\t<key>ImageIndex</key>", "\t\t\t<integer>%d</integer>" % k,
        "\t\t\t<key>ImageTotalNum</key>", "\t\t\t<integer>%d</integer>" % nz,
        "\t\t\t<key>ImageWidth</key>", "\t\t\t<integer>%d</integer>" % nx,
        "\t\t\t<key>NumberOfROIs</key>", "\t\t\t<integer>%d</integer>" % len(roi_xml),
        "\t\t\t<key>ROIs</key>", "\t\t\t<array>", *roi_xml, "\t\t\t</array>",
        "\t\t</dict>"]))

with open(os.path.join(out_dir, "ROI.xml"), "w") as f:
    f.write("\n".join([
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">',
        '<plist version="1.0">', "<dict>", "\t<key>Images</key>", "\t<array>",
        *images, "\t</array>", "</dict>", "</plist>", ""]))

with open("tests/testthat/phantom_truth.csv", "w") as f:
    f.write("roi_name,image_index,x_px,y_px,suv\n")
    for t in truth:
        f.write("%s,%d,%d,%d,%r\n" % (t[0], t[1], t[2], t[3], float(t[4])))

print("SUV factor:", suv_factor)
print("ROI voxels:", len(truth), "| tie pixel (row, col):", (r_out, c_out), "copies", (r_in, c_in))
