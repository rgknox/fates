import xarray as xr
import json
import numpy as np
import argparse
import datetime
import sys
import code  # For development: code.interact(local=dict(globals(), **locals()))
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)
import write_json


def clean_nested_list(data):
    """
    Recursively traverses lists to replace values > 1e+30 or NaNs with None.
        This preserves the nested [[...], [...]] structure.
    """
    threshold = 1e+30
    if isinstance(data, list):
        return [clean_nested_list(item) for item in data]
    else:
        # Check for threshold and NaN
        if data is None or np.isnan(data) or data >= threshold:
            return None
        return data

def get_structured_data(ds, var_name):
    """Converts DataArray to nested list and cleans values."""
    # Convert the specific variable to a numpy array first
    arr = ds[var_name].values
        
    # .tolist() on a 2D array creates the [[row1], [row2]] structure automatically
    nested_list = arr.tolist()
        
    return clean_nested_list(nested_list)
    
def compare_netcdf_files(file1_path, file2_path, diff_dic, json_path):
    
    # Load datasets
    ds1 = xr.open_dataset(file1_path)
    ds2 = xr.open_dataset(file2_path)
    
    output_data = {}
    
    vars1 = set(ds1.data_vars)
    vars2 = set(ds2.data_vars)

    # Initialize the two JSON objects
    result = {
        "new_parameters": {},
        "changed_parameters": {}
    }
    
    # 1. Find variables in file 2 but NOT in file 1
    unique_to_2 = vars2 - vars1

    # 2. Variables in both that have different data values
    common_vars = vars2.intersection(vars1)
    different_data = [v for v in common_vars if not ds1[v].equals(ds2[v])]

    #code.interact(local=dict(globals(), **locals()))

    numpft_ds1 = ds1.dims["fates_pft"]
    numpft_ds2 = ds2.dims["fates_pft"]

    if(numpft_ds1!=numpft_ds2):
        print("The two parameter sets must have the same number of PFTS")
        print("aborting")
        exit(2)

    if(numpft_ds1!=14):
        print("WARNING: You are looking at diffs between two files")
        print("The intention is to record these diffs, and then apply")
        print("them to the default fates parameter file. This had 14")
        print("PFTs last time we checked. So it is unclear how you")
        print("would apply these diffs, since the PFT set is different")
        print("If you are ok with this, override this exit")
        exit(2)

    diff_dic["pft_trim_list"] = list(range(1, numpft_ds1 + 1))
    str_pft_list = ",".join(map(str, diff_dic["pft_trim_list"]))

    diff_dic["parameters"] = {}
    diff_dic["parameters"]["pft_parameters"] = {}
    diff_dic["parameters"]["pft_parameters"][str_pft_list] = {}
    diff_dic["parameters"]["non_pft_parameters"] = {}

    print("\n\nThe following parameters were not in the base file")
    print(" and will not be added to the patch file.")
    print(" Feel free to manually add these to your patch if")
    print(" you believe they are in the file you will apply")
    print(" this patch to.")
    for var in unique_to_2:
        print(f"New Parameter: {var}")
        print(get_structured_data(ds2, var))
        
    for var in different_data:
        if("fates_pft" in ds2[var].dims):
            diff_dic["parameters"]["pft_parameters"][str_pft_list][var] = get_structured_data(ds2, var)
        else:
            diff_dic["parameters"]["non_pft_parameters"][var] = get_structured_data(ds2,var)


    print("\nParsing differences complete, dictionary built, writing patch to file")
    with open(json_path, 'w') as outfile:
        write_json.traverse_data(outfile,diff_dic)

    print("\nComplete, exiting")

def main():
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    #
    parser.add_argument('--basefile',dest='basefile', type=str, help="base netcdf file")
    parser.add_argument('--newfile',dest='newfile',type=str,help="new netdf file")
    parser.add_argument('--jsonout',dest='jsonout',type=str,help="json text output file")
    
    args = parser.parse_args()

    time_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    change_log = time_str+': modify_fates_paramfile_json.py'+' '.join(sys.argv[1:])

    # This is the root of the new dictionary we are creating. This will be dumped into
    # a json file.
    diff_dic = {}
    diff_dic["notes"] = "This patch file was created via:"+change_log

    diff_dic["usage"] = f"This file should be passed as input to the batch patch python script in fates, from the parameter_files directory:, like so: (cd parameter_files;python tools/batch_patch_params.py --fin {args.jsonout})"

    diff_dic["base_file"] = "UNSPECIFIED_BASEFILE.json"
    diff_dic["new_file"]  = "UNSPECIFIED_NEWFILE.json"
    
    compare_netcdf_files(args.basefile, args.newfile, diff_dic, args.jsonout)
    

if __name__ == "__main__":
    main()
