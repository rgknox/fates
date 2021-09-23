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

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-active', action='store_true', help="active variables only")
    group.add_argument('-inactive', action='store_true', help="inactive variables only, this is default")

    parser.add_argument('--f', '--input', dest='fnamein', type=str, help="Path to FatesHistoryInterfaceMod.F90", required=True)
    parser.add_argument('-hydro-active', action='store_true', help="include hydro variables")
    parser.add_argument('-nitr-active', action='store_true', help="include nitrogen variables")
    parser.add_argument('-phos-active', action='store_true', help="include phosphorous variables")

    parser.add_argument('-multiplex-exclude', action='store_true', help="exclude multiplexed variables")

    args = parser.parse_args()

    # Load up the pft names if they are there
    inactive_list = []
    file_obj = open(args.fnamein, "r")
    contents = file_obj.read()
    cont_split = contents.split('\n')

    hydro_excl = False
    nitr_excl = False
    phos_excl = False

    # Default to inactive even if flag not used
    use_default = 'inactive'

    if(args.active):
        use_default = 'active'

    elif(args.inactive):
        use_default = 'inactive'

    # Multiplexed list
    # ! scpf = size class x PFT
    # ! cacpf = cohort age class x PFT
    # ! cnlf = canopy layer x leaf layer
    # ! cnlfpft = canopy layer x leaf layer x PFT
    # ! scag = size class bin x age bin
    # ! scagpft = size class bin x age bin x PFT
    # ! agepft  = age bin x PFT
    # ! agefuel = age bin x fuel size class
    multip_list = ['cwdsc','elem','scls','scpf','cacls','cacpf','cnlf','cnlfpft','scag','scagpft','age','agepft','agefuel']

    # Initialize global counts
    total = 0

    print('\nListHist.py call inputs:')
    print('------------------------')
    print('use_default: {}'.format(use_default))
    print('args.hydro_active: {}'.format(args.hydro_active))
    print('args.nitr_active: {}'.format(args.nitr_active))
    print('args.phos_active: {}'.format(args.phos_active))

    print('\n\nHISTORY VARIABLES: \n-------------------------------\n')

    for line_num,line_str in enumerate(cont_split):

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

            if('use_default' in line_str):
                if(use_default in line_str):
                    # Work backwards until we find "vname"
                    not_found_vname = True
                    srch_num = line_num
                    count = 0
                    while(not_found_vname):
                        srch_num = srch_num-1
                        count += 1
                        if('vname=' in cont_split[srch_num]):
                            vname_string = cont_split[srch_num].split('\'')[1]
                            # print(cont_split[srch_num])
                            # print(vname_string)
                            not_found_vname = False
                        if(count>3):
                            not_found_vname = False
                    # Work forwards to see if index is multiplexed
                    not_found_index = True
                    while(not_found_index):
                        srch_num += 1
                        if('index =' in cont_split[srch_num]):
                            suffix = ((cont_split[srch_num].split('_')[-1]).split(')')[0]).strip()
                            # print(cont_split[srch_num])
                            # print(suffix)
                            not_found_index = False
                            if(not(suffix != 'si' and args.multiplex_exclude)):
                                print('\'{}\','.format(vname_string))
                                # print('line num: {}'.format(line_num))
                                total += 1

    print('\n-------------------------------\n')
    print('Variables listed: {}'.format(total))

# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main()
