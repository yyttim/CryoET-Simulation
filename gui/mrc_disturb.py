# -*- coding: utf-8 -*-
# @Time    : 2024/9/3 7:55
import mrcfile
import numpy as np

input_file_path = '../data/templates/mb_mrcs_10A/5ide.mrc'
output_file_path = '../data/templates/mb_mrcs_10A/5ide_perturbed_noise.mrc'

with mrcfile.open(input_file_path, mode='r') as mrc:
    # Get the original data and header
    original_data = mrc.data
    header = mrc.header

    # Perturb the data in the MRC file
    perturbed_data = original_data + np.random.normal(scale=0.5, size=original_data.shape)
    perturbed_data = perturbed_data.astype(np.float32)  # Convert to the appropriate data type

    # Create a new MRC file with the perturbed data and transfer the header information
    with mrcfile.new(output_file_path, overwrite=True) as new_mrc:
        new_mrc.set_data(perturbed_data)
        for key in header.dtype.names:
            new_mrc.header[key] = header[key]