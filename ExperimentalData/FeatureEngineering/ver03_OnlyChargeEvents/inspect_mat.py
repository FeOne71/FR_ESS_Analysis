import scipy.io
import pandas as pd
import os

# Path to the mat file
mat_file_path = r'd:\JCW\Projects\KEPCO_ESS_Local\ExperimentalData\FeatureEngineering\ver03_OnlyChargeEvents\Results\min10_std2_rng0_Both\Feature_Summary_Table_Both.mat'

try:
    mat = scipy.io.loadmat(mat_file_path)
    print("Keys in mat file:", mat.keys())
    
    # Inspect specific keys if found
    for key in mat.keys():
        if not key.startswith('__'):
            val = mat[key]
            print(f"Key: {key}, Type: {type(val)}")
            if hasattr(val, 'shape'):
                print(f"  Shape: {val.shape}")
            # Try to convert to dataframe if it looks tabular
            # ...
            
except Exception as e:
    print(f"Error loading mat file: {e}")
