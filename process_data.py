#!/usr/bin/env python
import os
import h5py
import numpy as np
import scipy.sparse as sp
from scipy.io import mmwrite
import gzip
import shutil


h5_file = os.path.join(os.environ['CELLBENDER_OUT'], "output_file_filtered.h5")
out_dir = os.path.join(os.environ['CELLBENDER_OUT'], "sc_out")
os.makedirs(out_dir, exist_ok=True)

def compress_file(filepath):
    with open(filepath, 'rb') as f_in, gzip.open(filepath + '.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(filepath)  # Optionally remove the uncompressed file

with h5py.File(h5_file, 'r') as f:
    mat_group = f['matrix']
    
    # Decode barcodes using numpy's vectorized string functions:
    barcodes = np.char.decode(mat_group['barcodes'][:], 'utf-8')
    
    features_group = mat_group['features']
    # Decode the features (id, name, feature_type) if they're stored as bytes:
    ids = np.char.decode(features_group['id'][:], 'utf-8')
    names = np.char.decode(features_group['name'][:], 'utf-8')
    types = np.char.decode(features_group['feature_type'][:], 'utf-8')
    features = np.column_stack((ids, names, types))
    
    data = mat_group['data'][:]
    indices = mat_group['indices'][:]
    indptr = mat_group['indptr'][:]
    shape = tuple(mat_group['shape'][:])

    expression_matrix = sp.csc_matrix((data, indices, indptr), shape=shape)

# Write the outputs into the 'sc_out' directory
mmwrite(os.path.join(out_dir, 'matrix.mtx'), expression_matrix)
np.savetxt(os.path.join(out_dir, 'barcodes.tsv'), barcodes, fmt='%s', delimiter='\t')
np.savetxt(os.path.join(out_dir, 'features.tsv'), features, fmt='%s', delimiter='\t')
# Compress each file
compress_file(os.path.join(out_dir, 'barcodes.tsv'))
compress_file(os.path.join(out_dir, 'features.tsv'))
compress_file(os.path.join(out_dir, 'matrix.mtx'))
