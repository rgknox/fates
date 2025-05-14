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
from CrownPhysics import CanopyElementPhysics
from SuppLeafPhysics import GetJmaxKp25Top
from PushParameters import PushParameters
from PushParameters import PushXMLPhotoParameters
from PushParameters import PushXMLRadParameters
from PushParameters import GetParamFromAttrib
from PushParameters import GetParamList

from CDLRead import CDLParse
import CtypesInit

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 11}

matplotlib.rc('font', **font)

plt.rcParams.update({'font.size': 12})


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
    iret = f90.set_leaf_param_sub(c8(1),ci(0),*ccharnb('fates_electron_transport_model'))

    # Maintenance respiration model # 1=Ryan (1991), 2=Atkin et al (2017)
    maintresp_leaf_model = float(GetParamFromAttrib(scalar_root,'fates_maintresp_leaf_model')[0])


    # Call this external to push the default parameters to the F90 objects
    PushXMLPhotoParameters(f90,xmlroot)
    PushXMLRadParameters(f90,xmlroot)
    
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

    site_root   = xmlroot.find('site_structure')
    ground_vis_albedo = GetParamList(site_root,'ground_vis_albedo','float')
    ground_nir_albedo = GetParamList(site_root,'ground_nir_albedo','float')
    frac_snow = GetParamList(site_root,'frac_snow','float')[0]
    n_layer = GetParamList(site_root,'nlayers','integer')[0]
    
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
        for ican in range(n_can):
            for icol in range(n_col):
                pft = elem_diags[ican][icol].pft
                if(pft>0):
                    site_diags,elem_diags[ican][icol] = \
                        CanopyElementPhysics(ican,icol,f90,met,it,xmlroot,pft, \
                                         vcmax25_top[ft],jmax25_top[ft],kp25_top[ft],lnc_top[ft], \
                                         co2_ppress_400ppm, o2_ppress_209kppm, btran_nolimit, \
                                         maintresp_leaf_model,site_diags,elem_diags[ican][icol])
        
                    
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
            x = site_diags.agross_rubpc3_vl[it,:]
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
    fig22,(ax1,ax2) = plt.subplots(1,2,figsize=(8.0,4.5))

    sunfrac_24_vl = met.GetDiurnalMean(site_diags.sunfrac_vl)
    y  = site_diags.vai_top
    ax1.plot(sunfrac_24_vl,y)
    ax1.invert_yaxis()
    ax1.set_ylabel('VAI')
    ax1.set_xlabel('Mean Sunlit Fraction')
    ax1.grid('on')
    

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
