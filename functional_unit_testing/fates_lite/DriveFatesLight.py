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

from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

from MetDrivers import GetCosZ
from MetDrivers import met_driver

from PushParameters import PushParameters
from PushParameters import PushXMLParameters
from PushParameters import GetParamFromAttrib

from CDLRead import CDLParse
import CtypesInit

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 11}

matplotlib.rc('font', **font)

plt.rcParams.update({'font.size': 14})


# Global constants to use in all Leaf Biophysics unit testing
# =======================================================================================

# Freezing point of water in Kelvin (at standard atmosphere)
tfrz_1atm = 273.15

# 25 degrees C in Kelvin (used because T25 functions)
leaf_tempk25 = tfrz_1atm + 25.0

# Daylight limitations can be imposed on Vcmax, a value of
# 1 means daylight length is at its maximum
dayl_factor_full = 1.0

# If Kumerathunga respiration is used, it requires moving averages
# of leaf temperature
t_growth_kum = -999
t_home_kum = -999

# Conversion factor for radiant energy in [Watts/m2] to
# photon flux density [umol/m2/s]
wm2_to_umolm2s = 4.6
 
# Simple conversion, number of micro-moles in a mole
umol_per_mol = 1.e6
mol_per_umol = 1.e-6

# Radiation constants, number of broadbands, and their indices (visible and NIR)
n_bands = 2
visb = 1
nirb = 2

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

# Nitrogen scaler at canopy top
nscaler_top = 1.0

# This signals that we perform a solution for a unit raditation calculation
# and then afterwards, scale by the magnitude of the downwelling flux
normalized_boundary = 1

# Initialize all of the fortran shared objects, and create aliases
# for the various subroutines and functions
# =======================================================================================

exec(open("../shared/py_src/CtypesInit.py").read())


# Subroutines
# =======================================================================================




def GetJmaxKp25Top(vcmax25_top):

    # Calculate Jmax and Kp at the canopy top at 25C
    # they scale off of vcmax
    #
    # jmax25_top:  Canopy top maximum electron transport
    #              rate at 25C (umol electrons/m**2/s)
    #
    # kp25top      Canopy top initial slope of CO2 response
    #              curve (C4 plants) at 25C
    
    jmax25_top = 1.67   * vcmax25_top
    kp25_top   = 20000.  * vcmax25_top
    
    # q10 response of product limited psn.
    # co2_rcurve_islope = co2_rcurve_islope25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
    
    return jmax25_top, kp25_top



    
        
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


    # Call this external to push the default parameters to the F90 objects
    #PushParameters(f90,params,dims)
    PushXMLParameters(f90,xmlroot)

    # code.interact(local=dict(globals(), **locals()))
    # Lets create a synthetic met driver based on nearest neighbor data?
    # ------------------------------------------------------------------------------------

    met_driver_csvfile = "bci_met_data/BCI_met_drivers_2003_2016.csv"

    met = met_driver(met_driver_csvfile)
    met.FilterTimes('daytime')

    
    ntimes = met.ndata

    # For cosine..
    zone = 12-5
    doys = [0,31,28,31,30,31,30,31,31,30,31,30]
    cdoys = np.cumsum(doys)
    latitude=9.153
    longitude=280.1539-360.0

    for it in range(ntimes):
        
        tvegk = met.data['t_veg'][it]

        
        doy = cdoys[met.data['mon'][it]-1] + met.data['day'][it]
        lhour = met.data['hod'][it]-5.0
        if(lhour<0):
            lhour = lhour-24
            doy   = doy-1

        met.data['cosz'][it] = GetCosZ(zone,doy,latitude,longitude,lhour)

    # Simple check to see if the cosz follows the diurnal signal
    #met.EvalCosz()
        
    # User specifications
    # -----------------------------------------------------------------------------------
    
    pft                   = 1

    fates_maintresp_leaf_model = 1     # 1=Ryan (1991), 2=Atkin et al (2017)
    
    ft                  = pft-1 # Python array pft index (starts with 0)

    total_lai             = 5.0
    total_sai             = 0.1*total_lai
    n_elem                = 1
    n_col                 = 1
    area                  = 1.0   # Assume a completely closed canopy
    ground_albedo_diff    = 0.3
    ground_albedo_beam    = 0.3
    frac_snow             = 0.0

    n_layer = 10 # number of discrete canopy layers

    avai = np.zeros(n_layer)  # Accumulated VAI (top of bin)
    
    rd_abs_leaf = np.zeros(n_layer)
    rb_abs_leaf = np.zeros(n_layer)
    r_abs_stem  = np.zeros(n_layer)
    rd_dn = np.zeros(n_layer)
    rd_up = np.zeros(n_layer)
    rbeam =  np.zeros(n_layer)
    sunfrac = np.zeros(n_layer)
    sunfrac_v2 = np.zeros(n_layer)
    for il in range(n_layer):
        avai[il] = (total_lai+total_sai)*float(il)/float(n_layer)
    davai =  (total_lai+total_sai)/float(n_layer)
    dalai =  total_lai/float(n_layer)
    vfrac = 1./float(n_layer)
    
    # The number of PFTs we actually simulate (do not change)
    nusepft = 1
    
    # Set some plant trait data
    # -----------------------------------------------------------------------------------
    # We calculate solutions for every independant output from the FATES model
    # so initialize output vectors with the same size (use r_b or any variable)

    # params['fates_leaf_vcmax25top'].data.shape = (1, 14) Fist index is "age"
    #vcmax25_top = float(params['fates_leaf_vcmax25top'].data[0,ft])
    vcmax25_top = float(GetParamFromAttrib(pft_root,'fates_leaf_vcmax25top')[ft])
    
    jmax25_top,kp25_top =  GetJmaxKp25Top(vcmax25_top)

    # Leaf Nitrogen Concentration at the canopy top (N:C ratio) (n/m2]
    # params['fates_stoich_nitr'].data.shape = (4, 14) First index is organ (leaf)
    # params['fates_stoich_nitr'].meta['units'] = 'gN/gC'
    #leaf_nc_ratio = params['fates_stoich_nitr'].data[0,ft]
    leaf_nc_ratio = float(GetParamFromAttrib(pft_root,'fates_stoich_nitr')[ft])

    #params['fates_leaf_slatop'].meta['units'] = 'm^2/gC'
    #leaf_slatop   = params['fates_leaf_slatop'].data[ft]
    leaf_slatop   = float(GetParamFromAttrib(pft_root,'fates_leaf_slatop')[ft])
    
    # Leaf N Conc at the canopy top [gN/m2]
    lnc_top       = leaf_nc_ratio/leaf_slatop


    
    
    # Create a vegetation canopy and initialize scattering elements
    # -----------------------------------------------------------------------------------

    iret = f90.alloc_twostream_sub(ci(n_elem),ci(n_col))
    iret = f90.param_prep_sub(ci(nusepft))
    iret = f90.setup_canopy_sub(c_int(1),c_int(1),c_int(pft),c_double(area),c_double(total_lai),c_double(total_sai))
    iret = f90.grndsnow_albedo_sub(c_int(visb),c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = f90.grndsnow_albedo_sub(c_int(visb),c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))
    iret = f90.grndsnow_albedo_sub(c_int(nirb),c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = f90.grndsnow_albedo_sub(c_int(nirb),c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))
    iret = f90.canopy_prep_sub(c8(frac_snow))

    

    
    # Initialize variables that are output from the f90 modules
    # these need to be defined as a cytype, and are floating point,
    # integer or strings
    
    qsat_f          = c_double(-9);    esat_f          = c_double(-9);    qsdt_f          = c_double(-9)
    esdt_f          = c_double(-9);    vcmax_f         = c_double(-9);    jmax_f          = c_double(-9)
    kp_f            = c_double(-9);    agross_f        = c_double(-9);    gstoma_f        = c_double(-9)
    anet_f          = c_double(-9);    lmr_f           = c_double(-9);    c13_f           = c_double(-9)
    ac_f            = c_double(-9);    aj_f            = c_double(-9);    ap_f            = c_double(-9)
    co2_interc_f    = c_double(-9);    mm_kco2_f       = c_double(-9);    mm_ko2_f        = c_double(-9)
    co2_cpoint_f    = c_double(-9);    gs0_f           = c_double(-9);    gs1_f           = c_double(-9)
    gs2_f           = c_double(-9);    solve_inter_f   = c_double(-9);    albedo_beam_f   = c_double(-9)
    albedo_diff_f   = c_double(-9);    canabs_beam_f   = c_double(-9);    canabs_diff_f   = c_double(-9)
    ffbeam_beam_f   = c_double(-9);    ffdiff_beam_f   = c_double(-9);    ffdiff_diff_f   = c_double(-9)
    rd_abs_leaf_f   = c_double(-9);    rb_abs_leaf_f   = c_double(-9);    r_abs_stem_f    = c_double(-9)
    r_abs_snow_f    = c_double(-9);    leaf_sun_frac_f = c_double(-9);    r_diff_dn_f     = c_double(-9)
    rd_abs_f        = c_double(-9);    rb_abs_f        = c_double(-9);    r_diff_up_f     = c_double(-9)
    r_beam_f        = c_double(-9);    leaf_sun_frac_v2_f = c_double(-9);
    
    # Initialize output arrays
    # -----------------------------------------------------------------------------------

    cosz = np.zeros(ntimes)
    lmr = np.zeros(ntimes)
    agross = np.zeros(ntimes)
    gstoma = np.zeros(ntimes)
    anet = np.zeros(ntimes)
    g_b_umol = np.zeros(ntimes)
    ag_limit = np.zeros([n_layer,2])
    ag_sslimit = np.zeros([n_layer,2,2])


    aglimit_which = []
    aglimit_apar  = []
    aglimit_temp  = []
    
    # Start the main model time loop
    # -----------------------------------------------------------------------------------
    
    for it in range(ntimes):
        
        tvegk = met.data['t_veg'][it]
        
        doy = cdoys[met.data['mon'][it]-1] + met.data['day'][it]
        lhour = met.data['hod'][it]-5.0
        if(lhour<0):
            lhour = lhour-24
            doy   = doy-1

        cosz[it] = GetCosZ(zone,doy,latitude,longitude,lhour)

        # Update scattering parameters based on zenith angle
        # optical depth, backscatter
        iret = f90.zenith_prep_sub(c8(cosz[it]))

        # Perform a normalized solve of the scattering environment
        iret = f90.solver_sub(ci(visb),ci(normalized_boundary),c8(1.0),c8(1.0), \
                       byref(albedo_beam_f),byref(albedo_diff_f), \
                       byref(canabs_beam_f),byref(canabs_diff_f), \
                       byref(ffbeam_beam_f),byref(ffdiff_beam_f),byref(ffdiff_diff_f))

        # Scale the normalized solution by the actual upper boundary conditions
        iret = f90.setdown_sub(ci(visb),c8(met.data['visbdn'][it]),c8(met.data['visddn'][it]))

        # Convert leaf boundary layer resistance to its reciprocal, conductance, and then
        # convert it from velocity units to micromoles/m2/s
        g_b_umol[it] = f90.velotomolarcf_fun(c8(met.data['can_press'][it]),c8(met.data['t_can'][it]))/met.data['r_b'][it]
        
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

        kn = f90.decaycoeffvcmax_fun(c8(float(GetParamFromAttrib(pft_root,'fates_leaf_vcmax25top')[ft])), \
                                     c8(float(GetParamFromAttrib(pft_root,'fates_leafn_vert_scaler_coeff1')[ft])), \
                                     c8(float(GetParamFromAttrib(pft_root,'fates_leafn_vert_scaler_coeff2')[ft])))
        
        par_abs_check = 0
        for il in range(n_layer):

            # Query the solver for the absorbed radiation over this interval
            iret = f90.getabsrad_sub(ci(1),ci(1),ci(visb),c8(avai[il]),c8(avai[il]+davai), \
                                     byref(rd_abs_f),byref(rb_abs_f), \
                                     byref(rd_abs_leaf_f),byref(rb_abs_leaf_f),byref(r_abs_stem_f), \
                                     byref(r_abs_snow_f),byref(leaf_sun_frac_f),byref(leaf_sun_frac_v2_f))

            rd_abs_leaf[il] = rd_abs_leaf_f.value
            rb_abs_leaf[il] = rb_abs_leaf_f.value
            r_abs_stem[il]  = r_abs_stem_f.value
            sunfrac[il]     = leaf_sun_frac_f.value
            sunfrac_v2[il]  = leaf_sun_frac_v2_f.value
            
            # Query the solver for the radiative intensities at the mid-point of this interval
            iret = f90.getintens_sub(ci(1), ci(1), ci(visb), c8(avai[il]+0.5*davai), \
                                     byref(r_diff_dn_f), byref(r_diff_up_f), byref(r_beam_f))

            rd_dn[il] = r_diff_dn_f.value
            rd_up[il] = r_diff_up_f.value
            rbeam[il] = r_beam_f.value
            
            

            # Nitrogen attenuation at this canopy depth
            nscaler = np.exp(-kn*(avai[il]+0.5*davai)*total_lai/(total_lai+total_sai))
            
            #par_abs_umol_m2 = (rd_abs_leaf[il] + rb_abs_leaf[il])*wm2_to_umolm2s/dalai

            # Scale down N and biophysical rates
            #nscaler = 
            iret = f90.biophysrate_sub(ci(pft), c8(vcmax25_top), \
                                       c8(jmax25_top), c8(kp25_top), \
                                       c8(nscaler_top), c8(tvegk), c8(dayl_factor_full), \
                                       c8(t_growth_kum),c8(t_home_kum),c8(btran_nolimit), \
                                       byref(vcmax_f), byref(jmax_f), byref(kp_f), \
                                       byref(gs0_f), byref(gs1_f), byref(gs2_f))

            # Leaf Maintenance Respiration (temp and pft dependent)
            if(fates_maintresp_leaf_model==1):
                iret = f90.lmr_ryan_sub(c8(lnc_top),c8(nscaler_top), ci(pft), c8(tvegk), byref(lmr_f))
            elif(fates_maintresp_leaf_model==2):
                iret = f90.lmr_atkin_sub(c8(lnc_top),c8(rdark_scaler),c8(tvegk),c8(atkin_mean_leaf_tempk),byref(lmr_f) )
            else:
                print('unknown leaf respiration model')
                exit(1)

            par_abs_check = par_abs_check + (r_abs_stem_f.value+r_abs_snow_f.value)
            #par_abs_check = par_abs_check + (rd_abs_f.value+rb_abs_f.value)
            
            for ipar in [0,1]:
                if(ipar==0):
                    areafrac = sunfrac_v2[il]
                    par_abs_leaf_umol = (rb_abs_leaf[il]/(dalai*areafrac) + rd_abs_leaf[il]/dalai)*wm2_to_umolm2s
                    
                else:
                    areafrac = 1.-sunfrac_v2[il]
                    par_abs_leaf_umol = (rd_abs_leaf[il]/dalai)*wm2_to_umolm2s

                # Energy conservation check (diff absorbed by leaf, stem and snow)
                par_abs_check = par_abs_check + par_abs_leaf_umol*areafrac*dalai/wm2_to_umolm2s
                
                    
                    
                if(par_abs_leaf_umol>400000):
                    print('high par_abs_leaf_umol:',ipar,il,met.data['visbdn'][it]*wm2_to_umolm2s, \
                          met.data['visddn'][it]*wm2_to_umolm2s, \
                          rb_abs_leaf[il]*wm2_to_umolm2s,dalai,areafrac,rd_abs_leaf[il]*wm2_to_umolm2s)
                    exit(0)

                    
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
                                               c8(g_b_umol[it]), \
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
                                               byref(solve_inter_f))
                 
                lmr[it]        = lmr[it] + areafrac*vfrac*lmr_f.value
                agross[it]     = agross[it] + areafrac*vfrac*agross_f.value
                gstoma[it]     = gstoma[it] + areafrac*vfrac*gstoma_f.value
                anet[it]       = anet[it] + areafrac*vfrac*anet_f.value

                agross_rubisco = f90.agross_rubiscoc3_fun(c8(vcmax_f.value), c8(co2_interc_f.value), \
                                                          c8(o2_ppress_209kppm), c8(co2_cpoint_f.value), \
                                                          c8(mm_kco2_f.value),c8(mm_ko2_f.value) )
                
                agross_rubpc3  = f90.agross_rubpc3_fun(c8(par_abs_leaf_umol), c8(jmax_f.value), \
                                                       c8(float(GetParamFromAttrib(pft_root,'fates_leaf_fnps')[ft])), \
                                                       c8(co2_interc_f.value),c8(co2_cpoint_f.value))
                aglimit_apar.append(par_abs_leaf_umol)
                aglimit_temp.append(tvegk)
                
                #code.interact(local=dict(globals(), **locals()))
                if(agross_rubisco<agross_rubpc3):
                    ag_limit[il,0] = ag_limit[il,0] + areafrac
                    ag_sslimit[il,ipar,0] = ag_sslimit[il,ipar,0] + areafrac
                    aglimit_which.append(1.)
                                         
                else:
                    ag_limit[il,1] = ag_limit[il,1] + areafrac
                    ag_sslimit[il,ipar,1] = ag_sslimit[il,ipar,1] + areafrac
                    aglimit_which.append(0.)

        plot_vert_prof = True
        if(plot_vert_prof):
            fig00, axs = plt.subplots(ncols=3,nrows=3,figsize=(7,7))
            ax1s = axs.reshape(-1)

            ax = ax1s[0]
            ax.plot(rd_abs_leaf,avai)
            ax.set_ylim([0,total_lai+total_sai])
            ax.invert_yaxis()
            ax.set_xlabel('Absorbed Par [W/m2]')
            ax.grid('on')
            
            ax = ax1s[1]
            ax.plot(sunfrac,avai,label='v1')
            ax.plot(sunfrac_v2,avai,label='v2')
            ax.set_ylim([0,total_lai+total_sai])
            ax.invert_yaxis()
            ax.set_xlabel('Sunlit Fraction [/]')
            ax.grid('on')
            ax.legend()
            
            ax = ax1s[2]
            ax.plot(rd_dn,avai+0.5*davai,label='Rd(dn)')
            ax.plot(rd_up,avai+0.5*davai,label='Rd(up)')
            ax.plot(rbeam,avai+0.5*davai,label='Rb')
            ax.set_ylim([0,total_lai+total_sai])
            ax.invert_yaxis()
            ax.set_xlabel('Par Flux Rate [W/m2]')
            ax.legend()
            ax.grid('on')
            
            ax = ax1s[3]
            ax.axis('off')
            #ax.text(0.25,0.5,'cosz: {} \nvisbdn: {5.2f} \nvisddn: {5.2f}'.format(cosz[it],met.data['visbdn'][it],met.data['visddn'][it]))
            ax.text(0.25,0.5,f"cosz: {cosz[it]:5.2f} \nvisbdn: {met.data['visbdn'][it]:5.2f} \nvisddn: {met.data['visddn'][it]:5.2f}")
            plt.show()

        # Check absorbed PAR
        # ----------------------------------------------------------------------------------
        err = met.data['visbdn'][it]*canabs_beam_f.value + met.data['visddn'][it]*canabs_diff_f.value - par_abs_check
        if(err>1.e-4):
            print('Energy Conservation Check on Absorbed Par FAILED, error: {} [W/m2]\n'.format(err))
            print(met.data['visbdn'][it],canabs_beam_f.value,met.data['visddn'][it],canabs_diff_f.value,par_abs_check)
            exit(2)
        

    fig7,(ax1,ax2) = plt.subplots(2,1,figsize=(5.5,7.5))

    ax1.scatter(aglimit_apar,aglimit_which)
    ax1.set_ylabel('1 = Rubisco, 0 = RuBP')
    ax1.set_xlabel('Apar [umol/m2/s]')
    ax1.grid('on')
    
    ax2.scatter(aglimit_temp,aglimit_which)
    ax2.set_ylabel('1 = Rubisco, 0 = RuBP')
    ax2.set_xlabel('Temperature [K]')
    ax2.grid('on')
           
    fig55, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(6.5,6.5))

    ax3.plot(ag_limit[:,1]/(ag_limit[:,0]+ag_limit[:,1]),avai[:])
    ax3.invert_yaxis()
    ax3.set_ylabel('LAI')
    ax3.set_xlabel('Fraction RuBP Limited')
    ax3.set_title('All')
    ax3.grid('on')
    
    ax1.plot(ag_sslimit[:,0,1]/(ag_sslimit[:,0,0]+ag_sslimit[:,0,1]),avai[:])
    ax1.invert_yaxis()
    ax1.set_ylabel('LAI')
    ax1.set_xlabel('Fraction RuBP Limited')
    ax1.set_title('Sunlit')
    ax1.grid('on')

    ax2.plot(ag_sslimit[:,1,1]/(ag_sslimit[:,1,0]+ag_sslimit[:,1,1]),avai[:])
    ax2.invert_yaxis()
    ax2.set_ylabel('LAI')
    ax2.set_xlabel('Fraction RuBP Limited')
    ax2.set_title('Shaded')
    ax2.grid('on')

    fig23, ax1= plt.subplots(figsize=(7.5,7.5))

    nhbins = 20
    boff   = 0
    
    hist,bins = np.histogram(aglimit_apar, bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax1.plot(binsc[boff:],hist[boff:]/len(aglimit_apar), color = 'k')

    
    fig22,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(7.5,7.5))

    nhbins = 20
    boff   = 0
    
    hist,bins = np.histogram(anet[:], bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax1.plot(binsc[boff:],hist[boff:]/ntimes, color = 'k')

    ax1.set_ylabel('Probability')
    ax1.set_xlabel('Anet [umol/m2/s]')
    ax1.grid('on')

    hist,bins = np.histogram(agross[:], bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax2.plot(binsc[boff:],hist[boff:]/ntimes, color = 'k')
    ax2.set_xlabel('Ag [umol/m2/s]')
    ax2.grid('on')

    hist,bins = np.histogram(gstoma[:]*1.e-6, bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax3.plot(binsc[boff:],hist[boff:]/ntimes, color = 'k')

    ax3.set_ylabel('Probability')
    ax3.set_xlabel('gs [mol/m2/s]')
    ax3.grid('on')

    hist,bins = np.histogram(lmr[:], bins=nhbins)
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
