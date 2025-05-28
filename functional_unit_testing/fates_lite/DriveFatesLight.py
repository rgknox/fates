# =======================================================================================
#
# For usage: $python LeafBiophysDriver.py --help
#
# This script runs unit tests on the leaf biophysics functions
#
#
# =======================================================================================

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colormaps
from datetime import datetime
import argparse
#from matplotlib.backends.backend_pdf import PdfPages
import platform
import xml.etree.ElementTree as et
import numpy as np
import matplotlib
import os
import sys
import getopt
import code  # For development: code.interact(local=dict(globals(), **locals()))
import time
import importlib
import csv
import subprocess
import re
import pandas as pd
import glob
#import streamlit as st
import ctypes
from ctypes import *
from operator import add
sys.path.append('../shared/py_src')
sys.path.append('../leaf_biophys')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb
from SetupCanopy import SetupCanopyDiags
from MetDrivers import GetCosZ
from MetDrivers import met_driver
from SuppLeafPhysics import GetJmaxKp25Top
from PushParameters import PushParameters
from PushParameters import PushXMLPhotoParameters
from PushParameters import PushXMLRadParameters
from PushParameters import GetParamFromAttrib
from PushParameters import GetParamList
from PushParameters import XMLToDic
from PushParameters import GetStrVecFromTag
from PushParameters import PushDictPhotoParameters
from PushParameters import PushDictRadParameters

from CDLRead import CDLParse
import CtypesInit

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 11}

matplotlib.rc('font', **font)

plt.rcParams.update({'font.size': 12})

plt.set_cmap('Dark2')

# Global constants to use in all Leaf Biophysics unit testing
# =======================================================================================

# Freezing point of water in Kelvin (at standard atmosphere)
tfrz_1atm = 273.15

# 25 degrees C in Kelvin (used because T25 functions)
leaf_tempk25 = tfrz_1atm + 25.0
 
# Simple conversion, number of micro-moles in a mole
umol_per_mol = 1.e6
mol_per_umol = 1.e-6

# Radiation constants, number of broadbands, and their indices (visible and NIR)
n_bands = 2
visb = 1
nirb = 2

# When there is hydrualic limitation on photosynthesis
# (via Vcmax reductions), then the btran factor is 1
btran_nolimit = 1.0

# Daylight limitations can be imposed on Vcmax, a value of
# 1 means daylight length is at its maximum
# WE DO NOT YET SCALE BY DAYL FACTOR...
dayl_factor_full = 1.0

# If Kumerathunga respiration is used, it requires moving averages
# of leaf temperature
t_growth_kum = -999
t_home_kum = -999

normalized_boundary = 1
absolute_boundary = 2

# 1 standard atmosphere in [Pa]
can_press_1atm = 101325.0

# Atmospheric CO2 partial pressure [Pa] at 400 ppm
co2_ppress_400ppm = 0.0004*can_press_1atm

# Atmospheric O2 partial pressure [Pa] %29.5 of atmosphere
o2_ppress_209kppm = 0.2095*can_press_1atm

# 70% of atmospheric CO2 is a reasonablish guess for
# intercellular CO2 concentration during primary production
# We can use this to test the gross assimilation routines
# directly, without having to solve for the equilibrium
# intercellular CO2
ci_ppress_static = 0.7*co2_ppress_400ppm

# When there is hydrualic limitation on photosynthesis
# (via Vcmax reductions), then the btran factor is 1
btran_nolimit = 1.0

# Respiration scaler at canopy top
rdark_scaler_top = 1.0

# Conversion factor for radiant energy in [Watts/m2] to
# photon flux density [umol/m2/s]
wm2_to_umolm2s = 4.6

# This signals that we perform a solution for a unit raditation calculation
# and then afterwards, scale by the magnitude of the downwelling flux
normalized_boundary = 1

# Initialize all of the fortran shared objects, and create aliases
# for the various subroutines and functions
# =======================================================================================

exec(open("../shared/py_src/CtypesInit.py").read())




# Subroutines
# =======================================================================================







        
# ========================================================================

def main(argv):

    # This class call instanteates all the fortran shared objects
    # and creates aliases for their functions and subroutines
    f90 = f90_modules('../shared/bld/')

    
    #numpft = int(xmlroot.find('numpft').text.strip())

    #cdlfile = "../../parameter_files/fates_params_default.cdl"
    #print('Loading Base parameters from: {}'.format(cdlfile))
    #params,dims = CDLParse(cdlfile,verbose=False)

    # Load controls
    xml_param_file = "fateslite_controls.xml"
    xmlroot = et.parse(xml_param_file).getroot()

    numpft = int(xmlroot.find('numpft').text.strip())
    
    # Set some of the major module switches (argument to subroutine is float)
    # -----------------------------------------------------------------------------------
    scalar_root = xmlroot.find('f90_params').find('scalar_dim')
    pft_root = xmlroot.find('f90_params').find('pft_dim')
    
    # Daylength factor: 1) scale vcmax and jmax,  0) do not scale
    daylength_switch = float(GetParamFromAttrib(scalar_root,'fates_daylength_factor_switch')[0])
    iret = f90.set_leaf_param_sub(c8(daylength_switch),ci(0),*ccharnb('fates_daylength_factor_switch'))

    # Stomatal Model: 1) Ball-Berry, 2) for Medlyn
    stomatal_switch = float(GetParamFromAttrib(scalar_root,'fates_leaf_stomatal_model')[0])
    iret = f90.set_leaf_param_sub(c8(stomatal_switch),ci(0),*ccharnb('fates_leaf_stomatal_model')) 

    # Stomatal assimilation model, ie 1) net or 2) gross will be used to estimate leaf co2 partial pressure
    stoma_assim_switch = float(GetParamFromAttrib(scalar_root,'fates_leaf_stomatal_assim_model')[0])
    iret = f90.set_leaf_param_sub(c8(stoma_assim_switch),ci(0),*ccharnb('fates_leaf_stomatal_assim_model'))

    # Photosynthesis Temperature acclimation 0) no acclimation, 1) kumarathunge et al 2019
    temp_acclim_switch = float(GetParamFromAttrib(scalar_root,'fates_leaf_photo_tempsens_model')[0])
    iret = f90.set_leaf_param_sub(c8(temp_acclim_switch),ci(0),*ccharnb('fates_leaf_photo_tempsens_model'))

    # Electron transport model  1) FvCB1980 2) JohnsonBerry2021
    electron_transp_switch = float(GetParamFromAttrib(scalar_root,'fates_electron_transport_model')[0])
    iret = f90.set_leaf_param_sub(c8(electron_transp_switch),ci(0),*ccharnb('fates_electron_transport_model'))

    # Maintenance respiration model # 1=Ryan (1991), 2=Atkin et al (2017)
    maintresp_leaf_model = float(GetParamFromAttrib(scalar_root,'fates_maintresp_leaf_model')[0])

    
    # Call this external to push the default parameters to the F90 objects

    # Testing an alternative read method. This populates
    # a dictionary for the parameters
    scalar_params,pft_params = XMLToDic(xmlroot)
    #lat = GetStrVecFromTag(xmlroot,'lat')
    PushDictPhotoParameters(f90,scalar_params,pft_params)
    PushDictRadParameters(f90,pft_params)
    
    # PushXMLPhotoParameters(f90,xmlroot)
    # PushXMLRadParameters(f90,xmlroot)
    
    # Set some plant trait data
    # *Note that the photosynthesis scheme takes vcmax25_top, jmax25_top and kp24_top
    # as ARGUMENTS and not pft traits... This is because FATES allows vcmax25_top to
    # change as a function of leaf age, and therefore this is not a parameter constant
    # at the pft level. In this unit framework we do not factor in leaf age, so
    # it is essentially a parameter constant
    # -----------------------------------------------------------------------------------

    vcmax25_top = np.zeros(numpft)
    jmax25_top = np.zeros(numpft)
    kp25_top = np.zeros(numpft)
    leaf_slatop = np.zeros(numpft)
    lnc_top = np.zeros(numpft)
    for ft in range(numpft):
        vcmax25_top[ft] = float(GetParamFromAttrib(pft_root,'fates_leaf_vcmax25top')[ft])

        # Jmax and Kp (co2 curve initial slope) are scaled off of Vcmax
        jmax25_top[ft],kp25_top[ft] =  GetJmaxKp25Top(vcmax25_top[ft])
        
        # Leaf Nitrogen Concentration at the canopy top (N:C ratio) (n/m2]
        leaf_nc_ratio = float(GetParamFromAttrib(pft_root,'fates_stoich_nitr')[ft])

        # Specific leaf area
        leaf_slatop[ft]   = float(GetParamFromAttrib(pft_root,'fates_leaf_slatop')[ft])
    
        # Leaf N Conc at the crown top [gN/m2]
        lnc_top[ft]       = leaf_nc_ratio/leaf_slatop[ft]

        
    # Lets setup the site information and the driver data
    # location and time-zone are for caluclating zenith angles from
    # local time info
    # ------------------------------------------------------------------------------------

    met_driver_csvfile = "bci_met_data/BCI_met_drivers_2003_2016.csv"
    site_root = xmlroot.find('site_info')
    latitude  = GetParamList(site_root,'lat','float')[0]
    longitude = GetParamList(site_root,'lon','float')[0]
    tzone     = GetParamList(site_root,'tzone','integer')[0]
    dvai      = GetParamList(site_root,'dvai','float')[0]
    
    if(longitude>180.):
        longitude = longitude - 360.

    met = met_driver(met_driver_csvfile,latitude,longitude,tzone,f90)
    ctrl_root = xmlroot.find('analysis_controls')
    met.FilterTimes(GetParamList(ctrl_root,'met_filter','string')[0])

    if(GetParamList(ctrl_root,'view_met_forcing','logical')[0]):
        met.EvalMetForcing()

    
    ntimes = met.ndata

    n_can,n_col,elem_diags,site_diags = SetupCanopyDiags(xmlroot,ntimes,dvai,f90)

    can_root   = xmlroot.find('canopy_structure')
    ground_vis_albedo = GetParamList(can_root,'ground_vis_albedo','float')
    ground_nir_albedo = GetParamList(can_root,'ground_nir_albedo','float')
    frac_snow = GetParamList(can_root,'frac_snow','float')[0]
    n_layer = GetParamList(can_root,'nlayers','integer')[0]
    btran_avg = GetParamList(can_root,'btran_avg','float')[0]
    
    # Site level scattering parameters
    iret = f90.grndsnow_albedo_sub(c_int(visb),c_double(ground_vis_albedo[1]),*ccharnb('albedo_grnd_diff'))
    iret = f90.grndsnow_albedo_sub(c_int(visb),c_double(ground_vis_albedo[0]),*ccharnb('albedo_grnd_beam'))
    iret = f90.grndsnow_albedo_sub(c_int(nirb),c_double(ground_nir_albedo[1]),*ccharnb('albedo_grnd_diff'))
    iret = f90.grndsnow_albedo_sub(c_int(nirb),c_double(ground_nir_albedo[0]),*ccharnb('albedo_grnd_beam'))
    iret = f90.canopy_prep_sub(c8(frac_snow))
    
    # Initialize variables that are output from the f90 modules
    # these need to be defined as a cytype, and are floating point,
    # integer or strings
    
    qsat_f          = c_double(-9);    esat_f          = c_double(-9);    qsdt_f          = c_double(-9)
    esdt_f          = c_double(-9);   
    albedo_beam_f   = c_double(-9)
    albedo_diff_f   = c_double(-9);    canabs_beam_f   = c_double(-9);    canabs_diff_f   = c_double(-9)
    ffbeam_beam_f   = c_double(-9);    ffdiff_beam_f   = c_double(-9);    ffdiff_diff_f   = c_double(-9)
    

    # Start the main model time loop
    # -----------------------------------------------------------------------------------

    pmod = int(ntimes/100)
    do_time_print = True
    for it in range(ntimes):

        ft = 0
        pft = 1
        
        tvegk = met.data['t_veg'][it]

        cosz  = met.data['cosz'][it]

        if( do_time_print and it>pmod and np.mod(it,pmod) == 0):
            print('{:04d}-{:02d}-{:02d}, {:.1f}% complete'.format(met.data['yr'][it], \
                                                      met.data['mon'][it], \
                                                      met.data['day'][it], \
                                                      100*float(it)/float(ntimes)))
        
        # Update scattering parameters based on zenith angle
        # optical depth, backscatter
        iret = f90.zenith_prep_sub(c8(cosz))

        # Perform a normalized solve of the scattering environment
        iret = f90.solver_sub(ci(visb),ci(normalized_boundary),c8(1.0),c8(1.0), \
                       byref(albedo_beam_f),byref(albedo_diff_f), \
                       byref(canabs_beam_f),byref(canabs_diff_f), \
                       byref(ffbeam_beam_f),byref(ffdiff_beam_f),byref(ffdiff_diff_f))

        # Scale the normalized solution by the actual upper boundary conditions
        iret = f90.setdown_sub(ci(visb),c8(met.data['visbdn'][it]),c8(met.data['visddn'][it]))
    
        # Lets walk through the various canopy elements and perform some
        # physics operations. These elements may be the crowns
        # of a single cohort, or could hypothetically be a mixed of plants
        # with the same scattering coefficients over depth. 
        # 1) retrieve the absorbed radiation
        # 2) perform photosynthesis and
        # 3) maintenance respiration
        # -------------------------------------------------------------------------------

        for ican,icol in site_diags.elem_veg:

            pft = elem_diags[ican][icol].pft
            ft  = pft - 1

            # Initialize values that communicate with fortran modules
            solve_iter_f    = c_int(-9);    
            rd_abs_leaf_f   = c_double(-9);    rb_abs_leaf_f   = c_double(-9);    r_abs_stem_f    = c_double(-9)
            r_abs_snow_f    = c_double(-9);    leaf_sun_frac_f = c_double(-9);    r_diff_dn_f     = c_double(-9)
            rd_abs_f        = c_double(-9);    rb_abs_f        = c_double(-9);    r_diff_up_f     = c_double(-9)
            r_beam_f        = c_double(-9);    
            vcmax_f         = c_double(-9);    jmax_f          = c_double(-9)
            kp_f            = c_double(-9);    agross_f        = c_double(-9);    gstoma_f        = c_double(-9)
            anet_f          = c_double(-9);    lmr_f           = c_double(-9);    c13_f           = c_double(-9)
            ac_f            = c_double(-9);    aj_f            = c_double(-9);    ap_f            = c_double(-9)
            co2_interc_f    = c_double(-9);    mm_kco2_f       = c_double(-9);    mm_ko2_f        = c_double(-9)
            co2_cpoint_f    = c_double(-9);    gs0_f           = c_double(-9);    gs1_f           = c_double(-9)
            gs2_f           = c_double(-9);

            # Vegetation Temperature [K]
            tvegk = met.data['t_veg'][it]
    
            # Set canopy gas parameters, such as the MM coefficients and the CO2 compensation point
            iret = f90.cangas_sub(c8(met.data['can_press'][it]), \
                                  c8(o2_ppress_209kppm), \
                                  c8(tvegk), \
                                  byref(mm_kco2_f), \
                                  byref(mm_ko2_f), \
                                  byref(co2_cpoint_f))
    
            # Determine the nitrogen attenuation rate in the canopy
            # -------------------------------------------------------------------------------
            # Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
            # kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al 
            # (2010) Biogeosciences, 7, 1833-1859
            pft_root = xmlroot.find('f90_params').find('pft_dim')
            kn = f90.decaycoeffvcmax_fun(c8(vcmax25_top), \
                                         c8(float(GetParamFromAttrib(pft_root,'fates_leafn_vert_scaler_coeff1')[ft])), \
                                         c8(float(GetParamFromAttrib(pft_root,'fates_leafn_vert_scaler_coeff2')[ft])))
    
            par_abs_check = 0
            for il in range(elem_diags[ican][icol].n_layer):

                # These area indices are WRT the element's top
                # --------------------------------------------
                vai_top = elem_diags[ican][icol].vai_top[il]
                vai_bot = elem_diags[ican][icol].vai_bot[il]
                dlai    = elem_diags[ican][icol].dlai
                dvai    = vai_bot - vai_top
                vai_mid = 0.5*(vai_bot+vai_top)
                sil     = elem_diags[ican][icol].sil[il]
                
                # Query the solver for the absorbed radiation over this interval
                iret = f90.getabsrad_sub(ci(ican+1),ci(icol+1),ci(visb),c8(vai_top),c8(vai_bot), \
                                         byref(rd_abs_f),byref(rb_abs_f), \
                                         byref(rd_abs_leaf_f),byref(rb_abs_leaf_f),byref(r_abs_stem_f), \
                                         byref(r_abs_snow_f),byref(leaf_sun_frac_f))

        
                # Query the solver for the radiative intensities at the mid-point of this interval
                iret = f90.getintens_sub(ci(ican+1), ci(icol+1), ci(visb), c8(vai_mid), \
                                         byref(r_diff_dn_f), byref(r_diff_up_f), byref(r_beam_f))
        
        
                # Nitrogen attenuation at this canopy depth
                leaf_frac = elem_diags[ican][icol].total_lai/(elem_diags[ican][icol].total_lai+elem_diags[ican][icol].total_sai)
                nscaler = np.exp(-kn*elem_diags[ican][icol].canopy_lai_mid[il])

                # Scale down N and biophysical rates
                iret = f90.biophysrate_sub(ci(pft), c8(vcmax25_top), \
                                           c8(jmax25_top), c8(kp25_top), \
                                           c8(nscaler), c8(tvegk), c8(met.data['dayl_factor'][it]), \
                                           c8(t_growth_kum),c8(t_home_kum),c8(btran_avg), \
                                           byref(vcmax_f), byref(jmax_f), byref(kp_f), \
                                           byref(gs0_f), byref(gs1_f), byref(gs2_f))

                # Leaf Maintenance Respiration (temp and pft dependent)
                if(maintresp_leaf_model==1):
                    iret = f90.lmr_ryan_sub(c8(lnc_top), c8(nscaler), ci(pft), c8(tvegk), byref(lmr_f))
                elif(maintresp_leaf_model==2):
                    iret = f90.lmr_atkin_sub(c8(lnc_top), c8(rdark_scaler), c8(tvegk), \
                                             c8(atkin_mean_leaf_tempk), byref(lmr_f) )
                else:
                    print('unknown leaf respiration model')
                    exit(1)

        
                par_abs_check = par_abs_check + (r_abs_stem_f.value+r_abs_snow_f.value)
                for ipar in [0,1]:
                    if(ipar==0):
                        lit_frac = leaf_sun_frac_f.value
                        par_abs_leaf_umol = (rb_abs_leaf_f.value/(dlai*lit_frac) + \
                                             rd_abs_leaf_f.value/dlai)*wm2_to_umolm2s
                    
                    else:
                        lit_frac = 1.- leaf_sun_frac_f.value
                        par_abs_leaf_umol = (rd_abs_leaf_f.value/dlai)*wm2_to_umolm2s

                    # Energy conservation check (diff absorbed by leaf, stem and snow)
                    par_abs_check = par_abs_check + par_abs_leaf_umol*lit_frac*dlai/wm2_to_umolm2s

                    
                    # Main call to FATES photosynthesis
            
                    iret = f90.leaflayerphoto_sub( c8(par_abs_leaf_umol),  \
                                                   ci(pft),   \
                                                   c8(vcmax_f.value),   \
                                                   c8(jmax_f.value),    \
                                                   c8(kp_f.value),      \
                                                   c8(gs0_f.value), \
                                                   c8(gs1_f.value), \
                                                   c8(gs2_f.value), \
                                                   c8(tvegk), \
                                                   c8(met.data['can_press'][it]), \
                                                   c8(co2_ppress_400ppm), \
                                                   c8(o2_ppress_209kppm), \
                                                   c8(met.data['vpress_sat'][it]), \
                                                   c8(met.data['g_b_umol'][it]), \
                                                   c8(met.data['vpress'][it]), \
                                                   c8(mm_kco2_f.value), \
                                                   c8(mm_ko2_f.value), \
                                                   c8(co2_cpoint_f.value), \
                                                   c8(lmr_f.value), \
                                                   c8(0.5), \
                                                   byref(agross_f), \
                                                   byref(gstoma_f), \
                                                   byref(anet_f), \
                                                   byref(c13_f), \
                                                   byref(co2_interc_f), \
                                                   byref(solve_iter_f))


                    agross_rubisco = f90.agross_rubiscoc3_fun(c8(vcmax_f.value), c8(co2_interc_f.value), \
                                                              c8(o2_ppress_209kppm), c8(co2_cpoint_f.value), \
                                                              c8(mm_kco2_f.value),c8(mm_ko2_f.value) )
                
                    agross_rubp  = f90.agross_rubpc3_fun(c8(par_abs_leaf_umol), c8(jmax_f.value), \
                                                         c8(float(GetParamFromAttrib(pft_root,'fates_leaf_fnps')[ft])), \
                                                         c8(co2_interc_f.value),c8(co2_cpoint_f.value))

                    # Diagnostics
                    # ----------------------------------------------------------------------
                
                    # We don't present the element diagnostic updater here to reduce verbosity in the main routine
                    elem_diags[ican][icol].Update(il,lit_frac,r_diff_dn_f.value,r_diff_up_f.value,r_beam_f.value,rd_abs_leaf_f.value,\
                                                  rb_abs_leaf_f.value, r_abs_stem_f.value,leaf_sun_frac_f.value, \
                                                  lmr_f.value,agross_f.value,gstoma_f.value,co2_interc_f.value )                
                
                    # Fraction that this leaf portion contributes to averages amongst the column
                    leaf_per_ground    = lit_frac*dlai*elem_diags[ican][icol].crown_area_frac

                    if(site_diags.lai_ground[sil]<0.000000001):
                        print(sil,site_diags.lai_ground[sil],il,ican,icol,site_diags.lai_ground[1])
                        exit(0)
            
                    # Fraction that this leaf portion contributes to averages amongst the layer
                    # Note its possible this solver is being called for a layer with only
                    # stems...
                    if(dlai>0.):
                        leaf_site_frac = leaf_per_ground/site_diags.lai_ground[sil]
                    else:
                        leaf_site_frac = 0.
                
                    site_diags.lmr[it]        = site_diags.lmr[it] + leaf_per_ground*lmr_f.value
                    site_diags.agross[it]     = site_diags.agross[it] + leaf_per_ground*agross_f.value
                    site_diags.gstoma[it]     = site_diags.gstoma[it] + leaf_per_ground*gstoma_f.value
                    site_diags.anet[it]       = site_diags.anet[it] + leaf_per_ground*anet_f.value
                    site_diags.r_abs_leaf[it] = site_diags.r_abs_leaf[it] + leaf_per_ground*par_abs_leaf_umol/wm2_to_umolm2s
                    
                    site_diags.sunfrac_vl[it,sil]    = site_diags.sunfrac_vl[it,sil] + leaf_site_frac*leaf_sun_frac_f.value
                    site_diags.lmr_vl[it,sil]        = site_diags.lmr_vl[it,sil] + leaf_site_frac*lmr_f.value
                    site_diags.agross_vl[it,sil]     = site_diags.agross_vl[it,sil] + leaf_site_frac*agross_f.value
                    site_diags.gstoma_vl[it,sil]     = site_diags.gstoma_vl[it,sil] + leaf_site_frac*gstoma_f.value
                    site_diags.anet_vl[it,sil]       = site_diags.anet_vl[it,sil] + leaf_site_frac*anet_f.value
                    site_diags.r_abs_leaf_vl[it,sil] = site_diags.r_abs_leaf_vl[it,sil] + leaf_site_frac*par_abs_leaf_umol/wm2_to_umolm2s
                    site_diags.vcmax_vl[it,sil]      = site_diags.vcmax_vl[it,sil] + leaf_site_frac*vcmax_f.value
                    site_diags.ag_rubisco_vl[it,sil] = site_diags.ag_rubisco_vl[it,sil] + leaf_site_frac*agross_rubisco
                    site_diags.ag_rubp_vl[it,sil]    = site_diags.ag_rubp_vl[it,sil] + leaf_site_frac*agross_rubp
                    site_diags.co2_interc_vl[it,sil] = site_diags.co2_interc_vl[it,sil] + leaf_site_frac*co2_interc_f.value
                
                    ssleaf_site_frac =  dlai*elem_diags[ican][icol].crown_area_frac/site_diags.lai_ground[sil]
                    if(ipar==0):
                        site_diags.r_abs_leaf_sunl[it,sil] = site_diags.r_abs_leaf_sunl[it,sil] + ssleaf_site_frac * par_abs_leaf_umol/wm2_to_umolm2s
                        site_diags.ag_rubisco_sunl[it,sil] = site_diags.ag_rubisco_sunl[it,sil] + ssleaf_site_frac * agross_rubisco
                        site_diags.ag_rubp_sunl[it,sil]    = site_diags.ag_rubp_sunl[it,sil] + ssleaf_site_frac * agross_rubp
                    else:
                        site_diags.r_abs_leaf_shal[it,sil] = site_diags.r_abs_leaf_shal[it,sil] + ssleaf_site_frac * par_abs_leaf_umol/wm2_to_umolm2s
                        site_diags.ag_rubisco_shal[it,sil] = site_diags.ag_rubisco_shal[it,sil] + ssleaf_site_frac * agross_rubisco
                        site_diags.ag_rubp_shal[it,sil]    = site_diags.ag_rubp_shal[it,sil] + ssleaf_site_frac * agross_rubp
                    
            # Display instantaneous profiles
            # Complete leaf layer and sun/shade loops first
            do_profiles = GetParamList(xmlroot.find('analysis_controls'),'view_element_inst_profs','logical')[0]
            if(do_profiles):
                elem_diags[ican][icol].PlotInstProfiles()

                
        # Complete element x leaf-layer x sun/shade loops
        # ------------------------------------------------------------------------------
        
        # Display instantaneous profiles
        # Complete leaf layer and sun/shade loops first
        ctrl_root = xmlroot.find('analysis_controls')
        do_profiles = GetParamList(ctrl_root,'view_site_inst_profs','logical')[0]
        if(do_profiles):
            fig_dp1,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(6.0,6.5))
            y  = site_diags.vai_top
            x  = site_diags.r_abs_leaf_vl[it,:]
            ax1.plot(x,y)
            ax1.invert_yaxis()
            ax1.set_ylabel('VAI')
            ax1.set_xlabel('Absorbed PAR umol/m2/s')
            ax1.grid('on')
            x = site_diags.vcmax_vl[it,:]
            ax2.plot(x,y)
            ax2.invert_yaxis()
            ax2.set_ylabel('VAI')
            ax2.set_xlabel('Vc,max umol/m2/s')
            ax2.grid('on')
            x = site_diags.agross_rubisco_vl[it,:]
            ax3.plot(x,y)
            ax3.invert_yaxis()
            ax3.set_ylabel('VAI')
            ax3.set_xlabel('Ag (rubisco) umol/m2/s')
            ax3.grid('on')
            x = site_diags.agross_rubp_vl[it,:]
            ax4.plot(x,y)
            ax4.invert_yaxis()
            ax4.set_ylabel('VAI')
            ax4.set_xlabel('Ag (RUBP) umol/m2/s')
            ax4.grid('on')
            plt.tight_layout()
            plt.show()

                    
    # Perform Diagnostics
    # ------------------------------------------------------------------------

    # Look at sun-shade fractions on the first element
    fig22,((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(8.5,6.5))
    #fig22,((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8)) = plt.subplots(2,4,figsize=(8.5,6.5))
    
    cmap = mpl.colormaps['Dark2']
    n_colors = 4
    colors = cmap(np.linspace(0, 1, n_colors))
    #code.interact(local=dict(globals(), **locals()))
    
    cosz_24 = met.GetDiurnalMean(met.data['cosz'])
    vai  = site_diags.vai_top

    icc=0
    sunfrac_24_vl = met.GetDiurnalMean(site_diags.sunfrac_vl)

    ics = []
    for ic,cosz in enumerate(cosz_24):
        if(cosz>0.01 and ((ic+1) % 3 == 0)):
            ics.append(ic)
            
    for icc,ic in enumerate(ics):
        ax1.plot(sunfrac_24_vl[ic,:],vai,label=f'{ic}:00',color=colors[icc])
    ax1.invert_yaxis()
    ax1.set_ylabel('VAI [m2/m2]')
    ax1.set_xlabel('Mean Sunlit Fraction')
    ax1.grid('on')
    ax1.legend()

    r_abs_leaf_24_vl = met.GetDiurnalMean(site_diags.r_abs_leaf_vl)
    for icc,ic in enumerate(ics):
        ax2.plot(r_abs_leaf_24_vl[ic,:],vai,color=colors[icc])
    ax2.invert_yaxis()
    ax2.set_xlabel('Absorbed Vis [W/m2]')
    ax2.grid('on')
    
    ag_rubisco_24_vl =  met.GetDiurnalMean(site_diags.ag_rubisco_vl)
    ag_rubp_24_vl =  met.GetDiurnalMean(site_diags.ag_rubp_vl)
    for icc,ic in enumerate(ics):
        ax3.plot(ag_rubisco_24_vl[ic,:],vai,color=colors[icc])
    for icc,ic in enumerate(ics):
        ax3.plot(ag_rubp_24_vl[ic,:],vai,color=colors[icc],linestyle = '--')
    ax3.invert_yaxis()
    ax3.set_xlabel('Ag (Rubisco & RuBP)\n   [umol/m2/s]')
    ax3.grid('on')

    co2_interc_24_vl =  met.GetDiurnalMean(site_diags.co2_interc_vl)
    for icc,ic in enumerate(ics):
        ax4.plot(co2_interc_24_vl[ic,:],vai,color=colors[icc])
    ax4.invert_yaxis()
    ax4.set_xlabel('co2_c [Pa]')
    ax4.grid('on')
    ax4.set_ylabel('VAI [m2/m2]')
    ax4.axvline(x=co2_ppress_400ppm, color='k', linestyle='solid') 

    
    vcmax_24_vl = met.GetDiurnalMean(site_diags.vcmax_vl)
    for icc,ic in enumerate(ics):
        ax5.plot(vcmax_24_vl[ic,:],vai,color=colors[icc])
    ax5.invert_yaxis()
    ax5.set_xlabel('Vcmax [umol/m2/s]')
    ax5.grid('on')
   
    
    lmr_24_vl = met.GetDiurnalMean(site_diags.lmr_vl)
    for icc,ic in enumerate(ics):
        ax6.plot(lmr_24_vl[ic,:],vai,color=colors[icc])
    ax6.invert_yaxis()
    ax6.set_xlabel('Rleaf [umol/m2/s]')
    ax6.grid('on')


    #rabs_leaf_24_sunl = met.GetDiurnalMean(site_diags.r_abs_leaf_sunl)
    #for icc,ic in enumerate(ics):
    #    ax2.plot(rabs_leaf_24_sunl[ic,:],vai,color=colors[icc])
    #ax2.invert_yaxis()
    #ax2.set_xlabel('Rabs (sunlit) [umol/m2/s]')
    #ax2.grid('on')

    #ag_rubisco_24_sunl = met.GetDiurnalMean(site_diags.ag_rubisco_sunl)
    #ag_rubp_24_sunl = met.GetDiurnalMean(site_diags.ag_rubp_sunl)
    #for icc,ic in enumerate(ics):
    #    ax3.plot(ag_rubisco_24_sunl[ic,:],vai,color=colors[icc])
    #for icc,ic in enumerate(ics):
    #    ax3.plot(ag_rubp_24_sunl[ic,:],vai,color=colors[icc],linestyle = '--')
    #ax3.invert_yaxis()
    #ax3.set_xlabel('Sunlit Ag (Rubisco & RuBP)\n   [umol/m2/s]')
    #ax3.grid('on')
    
    
    plt.tight_layout()
    plt.show()

    
    fig23, ax1= plt.subplots(figsize=(7.5,7.5))

    nhbins = 20
    boff   = 0
    
    hist,bins = np.histogram(site_diags.aglimit_apar, bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax1.plot(binsc[boff:],hist[boff:]/len(site_diags.aglimit_apar), color = 'k')

    
    fig22,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(7.5,7.5))

    nhbins = 20
    boff   = 0
    
    hist,bins = np.histogram(site_diags.anet[:], bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax1.plot(binsc[boff:],hist[boff:]/ntimes, color = 'k')

    ax1.set_ylabel('Probability')
    ax1.set_xlabel('Anet [umol/m2/s]')
    ax1.grid('on')

    hist,bins = np.histogram(site_diags.agross[:], bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax2.plot(binsc[boff:],hist[boff:]/ntimes, color = 'k')
    ax2.set_xlabel('Ag [umol/m2/s]')
    ax2.grid('on')

    hist,bins = np.histogram(site_diags.gstoma[:]*1.e-6, bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax3.plot(binsc[boff:],hist[boff:]/ntimes, color = 'k')

    ax3.set_ylabel('Probability')
    ax3.set_xlabel('gs [mol/m2/s]')
    ax3.grid('on')

    hist,bins = np.histogram(site_diags.lmr[:], bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax4.plot(binsc[boff:],hist[boff:]/ntimes, color = 'k')

    ax4.set_xlabel('Rl [umol/m2/s]')
    ax4.grid('on')
    
    plt.tight_layout()
    plt.show()


    
    print('Deallocating parameter space')
    iret = f90.dealloc_leaf_param_sub()

    exit(0)


    
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
