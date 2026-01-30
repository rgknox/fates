import xarray as xr
import json
import numpy as np
import argparse
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
    
def compare_netcdf_files(file1_path, file2_path, json_path):
    
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
    

    for var in unique_to_2:
        result["new_parameters"][var] = get_structured_data(ds2, var)
        
    for var in different_data:
        result["changed_parameters"][var] = get_structured_data(ds2, var)


        
    # Print as formatted JSON
    #print(json.dumps(output_data, indent=4))
    with open(json_path, 'w') as outfile:

        write_json.traverse_data(outfile,result)
        

def main():
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    #
    parser.add_argument('--basefile',dest='basefile', type=str, help="base netcdf file")
    parser.add_argument('--newfile',dest='newfile',type=str,help="new netdf file")
    parser.add_argument('--jsonout',dest='jsonout',type=str,help="json text output file")
    
    args = parser.parse_args()

    
    compare_netcdf_files(args.basefile, args.newfile, args.jsonout)
    

if __name__ == "__main__":
    main()
