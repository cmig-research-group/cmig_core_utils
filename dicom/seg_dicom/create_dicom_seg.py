import sys
from pathlib import Path
import os
import scipy.io as sio
import pydicom
from pydicom.sr.codedict import codes
import highdicom as hd

path_mat = sys.argv[1]
path_output = sys.argv[2]

mat_contents = sio.loadmat(path_mat)


# Mask volume
mask = mat_contents['seg_fields'][0,0]['data']['imgs'][0,0]
mask = mask.astype('bool')
if len(mask.shape) == 3:
    mask = mask.transpose(2, 1, 0)
else:
    mask = mask.transpose(2, 1, 0, 3)


# Reference DICOMs
path_ref_ims = mat_contents['seg_fields'][0,0]['source_images'][0]
path_ref_ims = Path(path_ref_ims)
image_files = path_ref_ims.glob('*')
image_datasets = [pydicom.dcmread(str(f)) for f in image_files]
for im in image_datasets:
    if not hasattr(im, 'PatientBirthDate'):
        im.PatientBirthDate = None
    if not hasattr(im, 'AccessionNumber'):
        im.AccessionNumber = '69'
instance_numbers = [im.InstanceNumber for im in image_datasets]
inds = sorted(range(len(instance_numbers)), key=lambda k: instance_numbers[k])
image_datasets = [image_datasets[k] for k in inds]


# DICOM SEG requirements --------------------------------------------------------
# Algorithm
algo_name = mat_contents['seg_fields'][0,0]['algorithm'][0,0]['name'][0]
algo_version = mat_contents['seg_fields'][0,0]['algorithm'][0,0]['version'][0]
algo_family = mat_contents['seg_fields'][0,0]['algorithm'][0,0]['family'][0]
if algo_family == 'ArtificialIntelligence':
    algo_family_code = codes.cid7162.ArtificialIntelligence
    algo_type = hd.seg.SegmentAlgorithmTypeValues.AUTOMATIC

algorithm_identification = hd.AlgorithmIdentificationSequence(
    name=algo_name,
    version=algo_version,
    family=algo_family_code
)

# Descriptions
descriptions = []
fname = ''
for i in range(mat_contents['seg_fields'][0,0]['description'].size):
    seg_label = mat_contents['seg_fields'][0,0]['description'][0,i]['label'][0]
    seg_id = mat_contents['seg_fields'][0,0]['description'][0,i]['tracking_id'][0]
    seg_type = mat_contents['seg_fields'][0,0]['description'][0,i]['type'][0]
    if seg_type == 'Organ':
        seg_category = codes.cid7150.AnatomicalStructure
        seg_type = codes.cid7166.Organ
    elif seg_type == 'Abnormal':
        seg_category = codes.cid7150.MorphologicallyAbnormalStructure
        seg_type = codes.cid242.Abnormal
    tmp = hd.seg.SegmentDescription(
        segment_number=i+1,
        segment_label=seg_label,
        segmented_property_category=seg_category,
        segmented_property_type=seg_type,
        algorithm_type=algo_type,
        algorithm_identification=algorithm_identification,
        tracking_uid=hd.UID(),
        tracking_id=seg_id)
    descriptions.append(tmp)
    fname = fname + seg_label + '_'

# Seg instance
seg_series_desc = mat_contents['seg_fields'][0,0]['metadata'][0,0]['series_description'][0]
seg_series_num = mat_contents['seg_fields'][0,0]['metadata'][0,0]['series_number'][0][0]
seg_instance_num = mat_contents['seg_fields'][0,0]['metadata'][0,0]['instance_number'][0][0]
seg_manufacturer = mat_contents['seg_fields'][0,0]['metadata'][0,0]['manufacturer'][0]
seg_model = mat_contents['seg_fields'][0,0]['metadata'][0,0]['manufacturer_model_name'][0]
seg_software_version = mat_contents['seg_fields'][0,0]['metadata'][0,0]['software_versions'][0]
seg_device_serial = mat_contents['seg_fields'][0,0]['metadata'][0,0]['device_serial_number'][0]

seg_dataset = hd.seg.Segmentation(
    source_images=image_datasets,
    pixel_array=mask,
    segmentation_type=hd.seg.SegmentationTypeValues.BINARY,
    segment_descriptions=descriptions,
    series_instance_uid=hd.UID(),
    series_number=seg_series_num,
    sop_instance_uid=hd.UID(),
    instance_number=seg_instance_num,
    manufacturer=seg_manufacturer,
    manufacturer_model_name=seg_model,
    software_versions=seg_software_version,
    device_serial_number=seg_device_serial,
    series_description=seg_series_desc
)

seg_dataset.save_as(path_output + '/' + fname + 'SEG.dcm')
