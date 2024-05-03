"""
This script takes text data and stores it with .h5 format (binary) to reduce file size. 
Levels of compression between 0 and 9 can be used using compression_opts option.

Note: This script requires h5py, and is run outside of Abaqus using my own python environment.
"""

import numpy as np
import h5py

# Load the data from the .dat file
dir = "C:/Users/lsalo/abaqus/wilmington/output/plastic_initialization/extracted/S/ENTIRE_DOMAIN/"
fn = "Geostatic"
data = np.loadtxt(dir + fn + ".dat", skiprows=1, delimiter=",")  # skip the header row

# Open a new HDF5 file
with h5py.File(dir + fn + ".h5", 'w') as f:
    # Create a dataset for each column
    element_ids = f.create_dataset('ELEMENT_ID', data=data[:, 0], dtype='int16', compression="gzip", compression_opts=2)
    Sxx = f.create_dataset('Sxx', data=data[:, 1], dtype='float32', compression="gzip", compression_opts=2)
    Syy = f.create_dataset('Syy', data=data[:, 2], dtype='float32', compression="gzip", compression_opts=2)
    Szz = f.create_dataset('Szz', data=data[:, 3], dtype='float32', compression="gzip", compression_opts=2)
    Sxy = f.create_dataset('Sxy', data=data[:, 4], dtype='float32', compression="gzip", compression_opts=2)
    Sxz = f.create_dataset('Sxz', data=data[:, 5], dtype='float32', compression="gzip", compression_opts=2)
    Syz = f.create_dataset('Syz', data=data[:, 6], dtype='float32', compression="gzip", compression_opts=2)
