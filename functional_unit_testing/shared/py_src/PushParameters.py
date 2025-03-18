import ctypes
from ctypes import *
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

# Radiation constants, number of broadbands, and their indices (visible and NIR)
n_bands = 2
visb = 1
nirb = 2

def PushParameters(f90,params,dims):


    numpft = dims['fates_pft']
    
    # Allocate and push  photosynthesis parameters
    # -------------------------------------------------------------------------------------
    print('Allocating parameter space for {} pfts'.format(numpft))
    iret = f90.alloc_leaf_param_sub(ci(numpft))

    # These are photosynthesis PFT parameters that also need to be pushed to the fortran objects
    f90_photo_pft_params = ['fates_leaf_stomatal_btran_model','fates_leaf_agross_btran_model', \
                            'fates_leaf_c3psn','fates_leaf_stomatal_slope_ballberry', \
                            'fates_leaf_stomatal_slope_medlyn','fates_leaf_fnps', \
                            'fates_leaf_stomatal_intercept','fates_maintresp_reduction_curvature', \
                            'fates_maintresp_reduction_intercept','fates_maintresp_reduction_upthresh', \
                            'fates_maintresp_leaf_atkin2017_baserate','fates_maintresp_leaf_ryan1991_baserate', \
                            'fates_leaf_vcmaxha','fates_leaf_jmaxha', \
                            'fates_leaf_vcmaxhd','fates_leaf_jmaxhd', \
                            'fates_leaf_vcmaxse','fates_leaf_jmaxse']

    
    # Push Parameter File values to the fortran objects, also save some values in local lists
    for param_name in f90_photo_pft_params:
        for pft in range(numpft):
            iret = f90.set_leaf_param_sub(c8(float(params[param_name].data[pft])),ci(pft+1),*ccharnb(param_name))

            
    # Allocate and push radiation parameters
    # -----------------------------------------------------------------------------------
    iret = f90.alloc_radparams_sub(ci(numpft),ci(n_bands))
    for ft in range(numpft):
        pft = ft+1
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_rhovis'].data[ft])), c_int(pft),c_int(visb),*ccharnb("rhol"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_rhonir'].data[ft])), c_int(pft),c_int(nirb),*ccharnb("rhol"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_tauvis'].data[ft])), c_int(pft),c_int(visb),*ccharnb("taul"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_taunir'].data[ft])), c_int(pft),c_int(nirb),*ccharnb("taul"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_stem_rhovis'].data[ft])), c_int(pft),c_int(visb),*ccharnb("rhos"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_stem_rhonir'].data[ft])), c_int(pft),c_int(nirb),*ccharnb("rhos"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_stem_tauvis'].data[ft])), c_int(pft),c_int(visb),*ccharnb("taus"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_stem_taunir'].data[ft])), c_int(pft),c_int(nirb),*ccharnb("taus"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_xl'].data[ft])), c_int(pft), c_int(0),*ccharnb("xl"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_clumping_index'].data[ft])), c_int(pft),c_int(0),*ccharnb("clumping_index"))


    
        
