# -*- coding: utf-8 -*-
# @Time    : 2024/9/3 7:55
import mrcfile
import numpy as np

input_file_path = '../data/templates/mb_mrcs_10A/5ide.mrc'
output_file_path = '../data/templates/mb_mrcs_10A/5ide_perturbed.mrc'

with mrcfile.open(input_file_path, mode='r') as mrc:
    original_data = mrc.data
    header = mrc.header

    mask = original_data > 0.1

    perturbed_data = original_data.copy()
    perturbed_data[mask] += np.random.normal(scale=0.5, size=perturbed_data[mask].shape)
    perturbed_data = perturbed_data.astype(np.float32)

    # Create a new MRC
    with mrcfile.new(output_file_path, overwrite=True) as new_mrc:
        new_mrc.set_data(perturbed_data)
        for key in header.dtype.names:
            new_mrc.header[key] = header[key]