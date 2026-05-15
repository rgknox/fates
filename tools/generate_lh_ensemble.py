# =======================================================================================
#  This script generates an ensemble of parameter files using
#  Latin Hypercube Sampling. Beware, this could generate a lot of files!
#  It uses the scipy Quasi Monte Carlo (QMC) library.
# =======================================================================================

import os
import argparse
import sys
import datetime
import time
import math
import code  # For development: code.interact(local=dict(globals(), **locals()))
import json
import shutil
import numpy as np
import write_json
import subprocess
import netCDF4 as nc
from scipy.stats import qmc


# USER SPECIFIES THIS SECTION
# Index of zero means that the parameter value should apply to all indices
# as copies (same value)
# KM_DECOMP_NH4 = 0.18
# KM_DECOMP_NO3 = 0.41
# KM_NIT = 0.76
# KM_DEN = 0.11
# "/home/rgknox/Models/InputDatasets/e3sm_input_datasets//lnd/clm2/paramdata/clm_params.cbgc.c07292018.nc"],

# Constants
nc_type=0
json_type=1
true=True
false=False

default_config  = {
    "files_input": ["../parameter_files/fates_params_default.json","clm_params.cbgc.c07292018.nc"],
    "files_output_pref":["fatesparams_lh","elmparams_lh"],
    "n_samples": 2,
    "notes":["Parameter metadata 'file' is indexed from 1 to 2, not 0 to 1, should match the file listed in 'files_input' above.",
             "Parameter metadata 'log' can be either '10','natural' or 'none'",
             "The value for min and max, if it is log scaled, is the ACTUAL VALUE, not the exponent!!!",
             "2D example 'fates_stoich_nitr':{'file': 1, 'dim_id': [0,1], 'log': true, 'min':2.0e-10, 'max':2.0e-6}"],
    "parameters":{
        "fates_cnp_vmax_nh4":       {"file": 1, "dim_id": [0],   "log": true, "min":2.0e-10, "max":2.0e-6},
        "fates_cnp_vmax_no3":       {"file": 1, "dim_id": [0],   "log": true, "min":2.0e-10, "max":2.0e-6},
        "fates_cnp_vmax_p":         {"file": 1, "dim_id": [0],   "log": true, "min":2.0e-11, "max":2.0e-7},
        "fates_cnp_eca_vmax_ptase": {"file": 1, "dim_id": [0],   "log": true, "min":2.0e-11, "max":2.0e-6},
        "fates_cnp_eca_km_nh4":     {"file": 1, "dim_id": [0],   "log": true, "min":0.01, "max":10},
        "fates_cnp_eca_km_no3":     {"file": 1, "dim_id": [0],   "log": true, "min":0.01, "max":10},
        "fates_cnp_eca_km_p":       {"file": 1, "dim_id": [0],   "log": true, "min":0.01, "max":10},
        "fates_cnp_eca_km_ptase":   {"file": 1, "dim_id": [0],   "log": true, "min":0.01, "max":10},
        "fates_cnp_nfix1":          {"file": 1, "dim_id": [0],   "log": false,"min":0.01, "max":0.4},
        "fates_turnover_fnrt":      {"file": 1, "dim_id": [0],   "log": true, "min":0.5,  "max":3.0},
        "fates_allom_l2fr":         {"file": 1, "dim_id": [0],   "log": true, "min":0.33, "max":3.0},
        "KM_NIT":                   {"file": 2, "dim_id": 0,     "log": true, "min":0.01, "max":10.0},
        "KM_DEN":                   {"file": 2, "dim_id": 0,     "log": true, "min":0.01, "max":10.0},
        "KM_DECOMP_NH4":            {"file": 2, "dim_id": 0,     "log": true, "min":0.01, "max":10.0},
        "KM_DECOMP_NO3":            {"file": 2, "dim_id": 0,     "log": true, "min":0.01, "max":10.0},
        "KM_DECOMP_P":              {"file": 2, "dim_id": 0,     "log": true, "min":0.01, "max":10.0}
    }
}

def main():


    # Parse argument(s)
    # Only one optional argument, the config file
    # -------------------------------------------
    
    time_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    change_log = time_str+': generate_lh_ensemble.py'+' '.join(sys.argv[1:])
    
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')

    parser.add_argument('--config', dest='config_file', type=str, required=False, \
                        help="JSON conifuration file, not required if using the in-line configuration found at the head of script")

    args = parser.parse_args()


    # Open and load the configuration file
    # ------------------------------------
    
    if(args.config_file is None):
        config = default_config
        print(f"\n\nUsing Default config file from script\n\n")
    else:
        print(f"\n\nOpening config file:{args.config_file}\n\n")
        with open(args.config_file, 'r') as file:
            config = json.load(file)


    # Define the LH array, then scale by the bounds in the config file
    # ----------------------------------------------------------------
            
    n_dims = len(config['parameters'])
    n_samples = config['n_samples']

    l_bounds = []
    u_bounds = []
    for pname, pobj in config['parameters'].items():
        if(pobj['log']):
            l_bounds.append(np.log10(float(pobj["min"])))
            u_bounds.append(np.log10(float(pobj["max"])))
        else:
            l_bounds.append(float(pobj["min"]))
            u_bounds.append(float(pobj["max"]))

    sampler = qmc.LatinHypercube(d=n_dims)
    lh_norm_sample = sampler.random(n=n_samples)
    lh_sample = qmc.scale(lh_norm_sample, l_bounds, u_bounds)
    print(f"Creating a Latin Hypercube")
    print(f"Number of samples {len(lh_sample)}")
    print(f"Number of parameters: {n_dims}")
    print(f"Resulting parameter array shape: {lh_sample.shape}\n\n")


    # Parameter file handling.  For each sample:
    # 1) Open the base file(s)
    # 2) Either load it to dictionary (json) or
    #    create a copy of it (nc)
    # 3) replace parameters with new values
    # 4) Write out finalized dictionary (json)
    # -----------------------------------------
    
    for k in range(n_samples):

        file_type=[]
        print(f"Sample: {k}")
        for j,file_in in enumerate(config['files_input']):

            file_out_base = config['files_output_pref'][j]

            if file_in.endswith('.nc'):
                file_type.append(nc_type)
                nc_file_in  = file_in
                nc_file_out = file_out_base+f"_{k+1:04d}"+".nc"
                try:
                    shutil.copy2(nc_file_in, nc_file_out)
                except IOError as e:
                    print(f"Unable to copy netcdf file. {e}")
                    exit(1)
                nc_obj = nc.Dataset(nc_file_out, 'r+')
                
            elif file_in.endswith('.json'):
                file_type.append(json_type)
                json_file_in  = file_in
                json_file_out = file_out_base+f"_{k+1:04d}"+".json"
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
        ip=-1
        for pname, pobj in config['parameters'].items():
            ip=ip+1

            if(pobj['log']):
                pvalue = 10.0**(lh_sample[k][ip])
            else:
                pvalue = lh_sample[k][ip]
            
            if(file_type[pobj['file']-1]==json_type):
                try:
                    file_pobj = json_obj['parameters'][pname]
                except KeyError:
                    print(f"Could not find parameter {pname}")
                    print(f"in the base json file: {json_file_in}")
                    print(f"exiting")
                    exit(0)
                    
                dim_names = json_obj['parameters'][pname]['dims']
                dim_sizes = []
                for idname in dim_names:
                    if(idname!='scalar'):
                        dim_sizes.append( json_obj['dimensions'][idname] )
                        
            elif(file_type[pobj['file']-1]==nc_type):

                nc_var = nc_obj.variables[pname]
                dim_names = list(nc_var.dimensions)
                dim_sizes = list(nc_var.shape)

                
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
                                    if(file_type[pobj['file']-1]==json_type):
                                        json_obj['parameters'][pname]['data'][i][j] = pvalue
                                    else:
                                        nc_var[i,j] = pvalue

                                    
                elif(len(pobj['dim_id'])==1):
                    for i in range(dim_sizes[0]):
                        if(i==(pobj['dim_id'][0]-1) or pobj['dim_id'][0]==0):
                            if(file_type[pobj['file']-1]==json_type):
                                json_obj['parameters'][pname]['data'][i] = pvalue
                            else:
                                nc_var[i] = pvalue
            else:
                if(file_type[pobj['file']-1]==json_type):
                    json_obj['parameters'][pname]['data'][0] = pvalue
                else:
                    nc_var[0] = pvalue

        # Parameters have been replaced
        # if JSON, now we push the dictionary to a file
        # if NC, we updated in place we are done
        
        for j,file_in in enumerate(config['files_input']):
            if(file_type[j]==json_type):
                with open(json_file_out, 'w') as outfile:
                    write_json.traverse_data(outfile,json_obj)
                print(f"wrote: {json_file_out}")
            else:
                nc_obj.sync()
                nc_obj.close()
                print(f"wrote: {nc_file_out}")
                
        print(f"\n")
                

# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main()
