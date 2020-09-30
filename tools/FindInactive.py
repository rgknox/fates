#!/usr/bin/env python
#
# This script scans the FatesHistoryInterfaceMod.F90 file
# to list out all the variables that are default inactive.
# We use this to re-populate the AllVars regression test
# primarily.  Note flags to filter-in variables
# that are included in hydro, nitrogen or phosphorus
# active runs



import sys
import os
import argparse
import code  # For development: code.interact(local=locals())

def main():
    
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')

    parser.add_argument('--f', '--input', dest='fnamein', type=str, help="Path to FatesHistoryInterfaceMod.F90", required=True)
    parser.add_argument('-hydro-active', action='store_true')
    parser.add_argument('-nitr-active', action='store_true')
    parser.add_argument('-phos-active', action='store_true')
    parser.add_argument('--check-name-length', dest='max_namelen', type=int, help="Maximum number of characters allowed in names",required=False)

    # As of 9/30/2020, the maximum number of charachters for history output is 32 characters, remove 6 for "FATES_" and that makes 26

    
    args = parser.parse_args()

    #code.interact(local=locals())

   

    
    # Load up the pft names if they are there
    inactive_list = []
    file_obj = open(args.fnamein, "r")
    contents = file_obj.read()
    cont_split = contents.split('\n')

    hydro_excl = False
    nitr_excl = False
    phos_excl = False

    if((args.max_namelen is None)):
        print('\n\nINACTIVE HISTORY VARIABLES: \n-------------------------------\n')
    else:
        print('\n\nLONG HISTORY VARIABLE NAMES: \n-------------------------------\n')
    
    for line_num,line_str in enumerate(cont_split):


        if(not (args.max_namelen is None)):
            if(('vname' in line_str)and('call' in line_str)and('set_history_var' in line_str)):
                #                code.interact(local=locals())
                varname_str = line_str.split('\'')[1]
                if(len(varname_str)>args.max_namelen):
                    print('variable: {} is {} characters'.format(varname_str,len(varname_str)))

        else:
                    

            # Check to see if we encountered an exlusion flag

            if(not(args.hydro_active) and ('hydro_active_if' in line_str) ):
                if(hydro_excl):
                    hydro_excl=False
                else:
                    hydro_excl=True

            if(not(args.nitr_active) and ('nitrogen_active_if' in line_str) ):
                if(nitr_excl):
                    nitr_excl=False
                else:
                    nitr_excl=True

            if(not(args.phos_active) and ('phosphorus_active_if' in line_str) ):
                if(phos_excl):
                    phos_excl=False
                else:
                    phos_excl=True


            if(not(hydro_excl or nitr_excl or phos_excl)):
                
                
                if ('inactive' in line_str):
                    # Work backwards until we find "vname"
                    not_found = True
                    srch_num = line_num
                    count = 0
                    while(not_found):
                        srch_num = srch_num-1
                        count = count+1
                        if('vname' in cont_split[srch_num]):
                            #print(cont_split[srch_num])
                            print('\'{}\','.format(cont_split[srch_num].split('\'')[1]))
                            not_found = False
                    
                        if(count>3):
                            not_found = False

    
    print('\n-------------------------------\n')
    
# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main()
