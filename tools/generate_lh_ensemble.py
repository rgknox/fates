# This script generates an ensemble of parameter files using
# Latin Hypercube Sampling. Beware, this will generate a
# lot of files!
# It uses the scipy Quasi Monte Carlo (QMC) library

import os
import argparse
import sys
import datetime
import time
import math
import code  # For development: code.interact(local=dict(globals(), **locals()))
import json
import numpy as np
import write_json
import subprocess

##from scipy.stats import qmc


# USER SPECIFIES THIS SECTION
# Index of zero means that the parameter value should apply to all indices
# as copies (same value)
# KM_DECOMP_NH4 = 0.18
# KM_DECOMP_NO3 = 0.41
# KM_NIT = 0.76
# KM_DEN = 0.11
# "/home/rgknox/Models/InputDatasets/e3sm_input_datasets//lnd/clm2/paramdata/clm_params.cbgc.c07292018.nc"],

default_config  = {
    "files_input": ["../parameter_files/fates_params_default.json"],
    "files_output_pref":["fatesparams_lh","elmparams_lh"],
    "n_samples": 2,
    "parameters":{
        "fates_stoich_nitr":  {"file": 1, "dim_id": [0,1], "min":2.0e-10, "max":2.0e-6},
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
        "fates_canopy_closure_thresh":  {"file": 1, "dim_id": 0, "min":0.33, "max":3.0},
        #"KM_NIT":{"file": 2, "dim_id": 0, "min":0.01, "max":10.0},
        #"KM_DEN":{"file": 2, "dim_id": 0, "min":0.01, "max":10.0},
        #"KM_DECOMP_NH4":{"file": 2, "dim_id": 0, "min":0.01, "max":10.0},
        #"KM_DECOMP_NO3":{"file": 2, "dim_id": 0, "min":0.01, "max":10.0},
        #"KM_DECOMP_P":{"file": 2, "dim_id": 0, "min":0.01, "max":10.0},
    }
}

def main():


    # Parse argument(s)
    # Only one optional argument, the config file
    # -------------------------------------------
    
    time_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    change_log = time_str+': generate_lh_ensemble.py'+' '.join(sys.argv[1:])
    
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')

    parser.add_argument('--config', dest='input_config', \
                        type=str, help="Input filename.  Required.", required=False)

    args = parser.parse_args()


    # Open and load the configuration file
    # ------------------------------------
    
    if(args.input_config is None):
        config = default_config
        print(f"\nUsing Default config file from script\n")
    else:
        print(f"\nOpening config file:{args.input_config}\n")
        with open(args.input_config, 'r') as file:
            config = json.load(file)


    # Define the LH array, then scale by the bounds in the config file
    # ----------------------------------------------------------------
            
    n_dims = len(config['parameters'])
    n_samples = config['n_samples']

    l_bounds = []
    u_bounds = []

    for pname, pobj in config['parameters'].items():
        l_bounds.append(float(pobj["min"]))
        u_bounds.append(float(pobj["max"]))

    ##sampler = qmc.LatinHypercube(d=n_dims)
    ##lh_norm_sample = sampler.random(n=n_samples)
    ##lh_sample = qmc.scale(lh_norm_sample, l_bounds, u_bounds)
    ##print(f"LH Simulations: {len(lh_sample)}")
    ##print(f"Sampler Shape: {lh_sample.shape}")


    # Parameter file handling.
    # For each sample:
    # 1) Open the base file(s)
    # 2) Either load it to dictionary (json) or
    #    create a copy of it (nc)
    # 3) replace parameters with new values
    # 4) Write out finalized dictionary (json)
    # -----------------------------------------
    
    for i in range(n_samples):

        for j,file_in in enumerate(config['files_input']):

            file_out_base = config['files_output_pref'][j]

            if file_in.endswith('.nc'):
                nc_file_in  = file_in
                nc_file_out = file_out_base+f"_{i+1:04d}"+".nc"
                
            elif file_in.endswith('.json'):
                json_file_in  = file_in
                json_file_out = file_out_base+f"_{i+1:04d}"+".json"
                try:
                    with open(file_in, 'r') as file:
                        json_obj = json.load(file)
                except FileNotFoundError:
                    print("The file path is incorrect or the file does not exist.")
                    exit(2)
                except json.JSONDecodeError:
                    print(f"Error: {file_in} contains invalid JSON formatting.")
                    exit(2)
                
            else:
                print("File extension not recognized.")
                exit(2)

        # Loop through all parameters
        for pname, pobj in config['parameters'].items():

            if(pobj['file']==1):
                try:
                    file_pobj = json_obj['parameters'][pname]
                except KeyError:
                    print(f"Could not find parameter {pname}")
                    print(f"in the base json file: {json_file_in}")
                    print(f"exiting")
                    exit(0)
                    
                dim_names = json_obj['parameters'][pname]['dims']
                
            else:
                # TO DO: INCLUDE NC ALTERNATIVE
                print("DNE")
                exit(0)
                
            dim_sizes = []
            for idname in dim_names:
                if(idname!='scalar'):
                    dim_sizes.append( json_obj['dimensions'][idname] )
            
            if isinstance(pobj['dim_id'], list):

                if(len(dim_names)!=len(pobj['dim_id'])):
                    print(f"The dimensionality of the parameter")
                    print(f"in your config did not match dimensionality")
                    print(f"in the base parameter file")
                    print(f"parameter: {pname}")
                    print(f"num dims in file: {len(dim_names)}")
                    print(f"file dims: {dim_names}")
                    len_dim =  len(pobj['dim_id'])
                    print(f"num dims in config: {len_dim}")
                    exit(2)
                
                if(len(pobj['dim_id'])>1):
                    for i in range(dim_sizes[0]):
                        if(i==(pobj['dim_id'][0]-1) or pobj['dim_id'][0]==0):
                            for j in range(dim_sizes[1]):
                                if(j==(pobj['dim_id'][1]-1) or pobj['dim_id'][1]==0):
                                    json_obj['parameters'][pname]['data'][i][j] = 42
                elif(len(pobj['dim_id'])==1):
                    for i in range(dim_sizes[0]):
                        if(i==(pobj['dim_id'][0]-1) or pobj['dim_id'][0]==0):
                            json_obj['parameters'][pname]['data'][i] = 42
            else:

                json_obj['parameters'][pname]['data'][0] = 42

        # Parameters have been replaced
        # if JSON, now we push the dictionary to a file
        # if NC, we updated in place we are done

        with open(json_file_out, 'w') as outfile:
            write_json.traverse_data(outfile,json_obj)


# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main()
