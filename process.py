#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 16:47:21 2024

@author: gabriel
"""

import SimpleITK as sitk
import numpy as np
import sys
import os

bold_file = os.path.abspath(sys.argv[1])
TR = float(sys.argv[2])
cutoff = sys.argv[3]
output_folder = os.path.abspath(sys.argv[4])

from rabies.preprocess_pkg.utils import convert_to_RAS
bold_file = convert_to_RAS(img_file=bold_file, out_dir=output_folder)

img = sitk.ReadImage(bold_file)
if cutoff=='none':
    img_subset = img
else:
    cutoff_idx = cutoff.split(',')
    if not len(cutoff_idx)==2:
        raise ValueError(f"{cutoff} syntax is wrong.")
    img_subset = img[:,:,:,int(cutoff_idx[0]):int(cutoff_idx[1])]
    print(f"{img_subset.GetSize()[3]} timepoints were selected.")
    

#sitk.WriteImage(img_2min, '2min.nii.gz')

arr = sitk.GetArrayFromImage(img_subset)

# create a median
median = np.median(arr, axis=0)
median_img = sitk.GetImageFromArray(median)
from rabies.utils import copyInfo_3DImage
median_img = copyInfo_3DImage(image_3d=median_img, ref_3d=img_subset)
sitk.WriteImage(median_img, f'{output_folder}/median.nii.gz')


arr = sitk.GetArrayFromImage(img_subset)
arr_2d = arr.reshape(arr.shape[0],-1)

# apply highpass
from rabies.confound_correction_pkg.utils import butterworth, smooth_image
filtered = butterworth(arr_2d, TR=TR,
                        high_pass=0.01, low_pass=None)

# convert to 4D sitk image
reshaped = filtered.reshape(arr.shape)
reshaped_img = sitk.GetImageFromArray(reshaped, isVector=False)
from rabies.utils import copyInfo_4DImage
reshaped_img = copyInfo_4DImage(image_4d=reshaped_img, ref_3d=img_subset, ref_4d=img_subset)

# apply smoothing
mask = np.ones(arr.shape[1:]).astype(int) # create a mask of the entire mask, because the smoothing function needs one
mask_img = sitk.GetImageFromArray(mask)
mask_img = copyInfo_3DImage(image_3d=mask_img, ref_3d=img_subset)

import nibabel as nb
affine = nb.load(bold_file).affine[:3,:3] # still not sure how to match nibabel's affine reliably
smoothed_img = smooth_image(reshaped_img, affine, fwhm=0.5, mask_img=mask_img)
sitk.WriteImage(smoothed_img, f'{output_folder}/processed_img.nii.gz')

