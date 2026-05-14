# This script generates an ensemble of parameter files using
# Latin Hypercube Sampling. Beware, this will generate a
# lot of files!
# It uses the scipy Quasi Monte Carlo (QMC) library


import numpy as np
from scipy.stats import qmc


# USER SPECIFIES THIS SECTION
# Index of zero means that the parameter value should apply to all indices
# as copies (same value)
# KM_DECOMP_NH4 = 0.18
# KM_DECOMP_NO3 = 0.41
# KM_NIT = 0.76
# KM_DEN = 0.11

default_config  = {
    "files_input": ["file1.json","file2.nc"],
    "files_output_pref":["fatesparams_lh","elmparams_lh"],
    "n_samples": 20,
    "parameters":{
        "fates_cnp_vmax_nh4": {"file": 1, "dim_id": [0], "min":2.0e-10, "max":2.0e-6},
        "fates_cnp_vmax_no3": {"file": 1, "dim_id": [0], "min":2.0e-10, "max":2.0e-6},
        "fates_cnp_vmax_p": {"file": 1, "dim_id": [0], "min":2.0e-11, "max":2.0e-7},
        "fates_cnp_eca_km_nh4": {"file": 1, "dim_id": [0], "min":0.01, "max":10},
        "fates_cnp_eca_km_no3": {"file": 1, "dim_id": [0], "min":0.01, "max":10},
        "fates_cnp_eca_km_p":   {"file": 1, "dim_id": [0], "min":0.01, "max":10},
        "fates_cnp_eca_vmax_ptase":   {"file": 1, "dim_id": [0], "min":0.01, "max":10},
        "fates_cnp_eca_km_ptase":   {"file": 1, "dim_id": [0], "min":0.01, "max":1},
        "fates_cnp_nfix1": {"file": 1, "dim_id": [0], "min":0.01, "max":0.4},
        "fates_turnover_fnrt": {"file": 1, "dim_id": [0], "min":0.5, "max":3},
        "fates_allom_l2fr": {"file": 1, "dim_id": [0], "min":0.33, "max":3.0},
        "KM_NIT":{"file": 2, "dim_id": 0, "min":0.01, "max":10.0},
        "KM_DEN":{"file": 2, "dim_id": 0, "min":0.01, "max":10.0},
        "KM_DECOMP_NH4":{"file": 2, "dim_id": 0, "min":0.01, "max":10.0},
        "KM_DECOMP_NO3":{"file": 2, "dim_id": 0, "min":0.01, "max":10.0},
        "KM_DECOMP_P":{"file": 2, "dim_id": 0, "min":0.01, "max":10.0},
    }
}


config = default_config

n_dims = len(config['parameters'])
n_samples = config['n_samples']

l_bounds = []
u_bounds = []

for pname, pobj in config['parameters'].items():
    l_bounds.append(float(pobj["min"]))
    u_bounds.append(float(pobj["max"]))

sampler = qmc.LatinHypercube(d=n_dims)
lh_norm_sample = sampler.random(n=n_samples)
lh_sample = qmc.scale(lh_norm_sample, l_bounds, u_bounds)

print(f"LH Simulations: {len(lh_sample)}")
print(f"Sampler Shape: {lh_sample.shape}")
#print(lh_sample)


# Step 3
# Create copies of the original file, and replace parameters with sample

for i in range(n_samples):

    for j,file_in in enumerate(config['files_input']):

        file_out_base = config['files_output_pref'][j]
        print(f"{file_in}, {file_out_base}")
