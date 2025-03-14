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
import xml.etree.ElementTree as ET
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
import streamlit as st
import ctypes
from ctypes import *
from operator import add
sys.path.append('../shared/py_src')

from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

from GetCosZ import GetCosZ
from MetDrivers import met_driver


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

    cdlfile = "../../parameter_files/fates_params_default.cdl"

    print('Loading Base parameters from: {}'.format(cdlfile))
    params,dims = CDLParse(cdlfile,verbose=False)

    # Set some of the major module switches (argument to subroutine is float)
    # -----------------------------------------------------------------------------------

    # Daylength factor: 1) scale vcmax and jmax,  0) do not scale
    iret = f90.set_leaf_param_sub(c8(1),ci(0),*ccharnb('fates_daylength_factor_switch'))

    # Stomatal Model: 1) Ball-Berry, 2) for Medlyn
    iret = f90.set_leaf_param_sub(c8(2),ci(0),*ccharnb('fates_leaf_stomatal_model')) 

    # Stomatal assimilation model, ie 1) net or 2) gross will be used to estimate leaf co2 partial pressure
    iret = f90.set_leaf_param_sub(c8(1),ci(0),*ccharnb('fates_leaf_stomatal_assim_model'))

    # Photosynthesis Temperature acclimation 0) no acclimation, 1) kumarathunge et al 2019
    iret = f90.set_leaf_param_sub(c8(0),ci(0),*ccharnb('fates_leaf_photo_tempsens_model'))

    # Electron transport model  1) FvCB1980 2) JohnsonBerry2021
    iret = f90.set_leaf_param_sub(c8(1),ci(0),*ccharnb('fates_electron_transport_model'))

    numpft = dims['fates_pft']
    
    # Allocating parameters
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
    leaf_c3psn = []
    leaf_stomatal_intercept = []
    for param_name in f90_photo_pft_params:
        for pft in range(numpft):
            iret = f90.set_leaf_param_sub(c8(float(params[param_name].data[pft])),ci(pft+1),*ccharnb(param_name))
            if(param_name == 'fates_leaf_c3psn'):
                leaf_c3psn.append(int(params[param_name].data[pft]))
            if(param_name == 'fates_leaf_stomatal_intercept'):
                leaf_stomatal_intercept.append(int(params[param_name].data[pft]))

    if(len(leaf_c3psn) != numpft):
        print('Did not find fates_leaf_c3psn')
        exit(1)

    
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


    # code.interact(local=dict(globals(), **locals()))
    # Lets create a synthetic met driver based on nearest neighbor data?
    # ------------------------------------------------------------------------------------
    
    met = met_driver("/home/rgknox/Downloads/bci_met.txt")
    met.FilterTimes('daytime')

    
    ntime = met.ndata

    # For cosine..
    zone = 12-5
    doys = [0,31,28,31,30,31,30,31,31,30,31,30]
    cdoys = np.cumsum(doys)
    latitude=9.153
    longitude=280.1539-360.0



    
    # User specifications
    # -----------------------------------------------------------------------------------
    
    pft                   = 1

    fates_maintresp_leaf_model = 1     # 1=Ryan (1991), 2=Atkin et al (2017)
    
    paft                  = pft-1 # Python array pft index (starts with 0)

    total_lai             = 5.0
    total_sai             = 0.1*total_lai
    n_elem                = 1
    n_col                 = 1
    area                  = 1.0   # Assume a completely closed canopy
    ground_albedo_diff    = 0.3
    ground_albedo_beam    = 0.3
    frac_snow             = 0.0

    
    # The number of PFTs we actually simulate (do not change)
    nusepft = 1
    
    # Set some plant trait data
    # -----------------------------------------------------------------------------------
    # We calculate solutions for every independant output from the FATES model
    # so initialize output vectors with the same size (use r_b or any variable)

    # params['fates_leaf_vcmax25top'].data.shape = (1, 14) Fist index is "age"
    vcmax25_top = float(params['fates_leaf_vcmax25top'].data[0,paft])
    
    jmax25_top,kp25_top =  GetJmaxKp25Top(vcmax25_top)

    # Leaf Nitrogen Concentration at the canopy top (N:C ratio) (n/m2]
    # params['fates_stoich_nitr'].data.shape = (4, 14) First index is organ (leaf)
    # params['fates_stoich_nitr'].meta['units'] = 'gN/gC'
    leaf_nc_ratio = params['fates_stoich_nitr'].data[0,paft]

    #params['fates_leaf_slatop'].meta['units'] = 'm^2/gC'
    leaf_slatop   = params['fates_leaf_slatop'].data[paft]

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

    # number of discrete canopy layers
    n_layer = 50
    avai = np.zeros(n_layer)
    rd_abs_leaf = np.zeros(n_layer)
    rb_abs_leaf = np.zeros(n_layer)
    sunfrac = np.zeros(n_layer)
    for il in range(n_layer):
        avai[il] = (total_lai+total_sai)*float(il)/float(n_layer)
    davai =  (total_lai+total_sai)/float(n_layer)
    dalai =  total_lai/float(n_layer)
    vfrac = 1./float(n_layer)


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
    r_abs_snow_f    = c_double(-9);    leaf_sun_frac_f = c_double(-9)

    
    # Initialize output arrays
    # -----------------------------------------------------------------------------------

    cosz = np.zeros(ntime)
    lmr = np.zeros(ntime)
    agross = np.zeros(ntime)
    gstoma = np.zeros(ntime)
    anet = np.zeros(ntime)
    g_b_umol = np.zeros(ntime)
    ag_limit = np.zeros([n_layer,2])
    ag_sslimit = np.zeros([n_layer,2,2])


    aglimit_which = []
    aglimit_apar  = []
    aglimit_temp  = []
    
    # Start the main model time loop
    # -----------------------------------------------------------------------------------
    
    for it in range(ntime):
        
        tvegk = met.data['t_veg'][it]
        
        doy = cdoys[met.data['mon'][it]-1] + met.data['day'][it]
        lhour = met.data['hod'][it]-5.0
        if(lhour<0):
            lhour = lhour-24
            doy   = doy-1

        cosz[it] = GetCosZ(zone,doy,latitude,longitude,lhour)

        #print("cosz: {}, Rbeam: {}, Rdiff: {}".format(cosz[it],met.data['Rbeam'][it],met.data['Rdiff'][it]))
        iret = f90.zenith_prep_sub(c8(cosz[it]))
        iret = f90.solver_sub(ci(visb),ci(normalized_boundary),c8(1.0),c8(1.0), \
                       byref(albedo_beam_f),byref(albedo_diff_f), \
                       byref(canabs_beam_f),byref(canabs_diff_f), \
                       byref(ffbeam_beam_f),byref(ffdiff_beam_f),byref(ffdiff_diff_f))
        iret = f90.setdown_sub(ci(visb),c8(met.data['Rbeam'][it]),c8(met.data['Rdiff'][it]))

        g_b_umol[it] = f90.velotomolarcf_fun(c8(met.data['can_press'][it]),c8(met.data['t_can'][it]))/met.data['r_b'][it]
        
        # Get Humidity
        iret = f90.qsat_sub(c8(met.data['t_can'][it]),c8(met.data['can_press'][it]), \
                            byref(qsat_f),byref(esat_f), \
                            byref(qsdt_f),byref(esdt_f))
        
        iret = f90.cangas_sub(c8(met.data['can_press'][it]), \
                              c8(o2_ppress_209kppm), \
                              c8(tvegk), \
                              byref(mm_kco2_f), \
                              byref(mm_ko2_f), \
                              byref(co2_cpoint_f))

        for il in range(n_layer):
            
            iret = f90.getabsrad_sub(ci(1),ci(1),ci(visb),c8(avai[il]),c8(avai[il]+davai), \
                                     byref(rd_abs_leaf_f),byref(rb_abs_leaf_f),byref(r_abs_stem_f), \
                                     byref(r_abs_snow_f),byref(leaf_sun_frac_f))
            
            rd_abs_leaf[il] = rd_abs_leaf_f.value
            rb_abs_leaf[il] = rb_abs_leaf_f.value
            sunfrac[il]     = leaf_sun_frac_f.value
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

            for ipar in [0,1]:
                if(ipar==0):
                    areafrac = sunfrac[il]
                    par_abs  = (rb_abs_leaf[il]/(dalai*areafrac) + rd_abs_leaf[il]/dalai)*wm2_to_umolm2s
                else:
                    areafrac = 1.-sunfrac[il]
                    par_abs  = (rd_abs_leaf[il]/dalai)*wm2_to_umolm2s
                    
                iret = f90.leaflayerphoto_sub( c8(par_abs),  \
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
                                               c8(esat_f.value), \
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
                agross_rubpc3  = f90.agross_rubpc3_fun(c8(par_abs), c8(jmax_f.value), \
                                                       c8(params['fates_leaf_fnps'].data[paft]), \
                                                       c8(co2_interc_f.value),c8(co2_cpoint_f.value))
                aglimit_apar.append(par_abs)
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
                    
            plot_vert_prof = False
            if(plot_vert_prof):
                fig00, axs = plt.subplots(ncols=1,nrows=2,figsize=(5,8))
                ax1s = axs.reshape(-1)

                ax = ax1s[0]
                ax.plot(rd_abs_leaf,avai)
                ax.set_ylim([0,total_lai+total_sai])
                ax.invert_yaxis()
                ax.set_xlabel('[W/m2 ground]')
                ax.set_title('Absorbed Par')
                
                ax = ax1s[1]
                ax.plot(np.cumsum(rd_abs_leaf),avai)
                ax.set_ylim([0,total_lai+total_sai])
                ax.invert_yaxis()
                ax.set_xlabel('[W/m2]')
                ax.set_title('Integrated Absorbed Par')
                
                plt.show()

        

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
    
    fig22,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(7.5,7.5))

    nhbins = 20
    boff   = 0
    
    hist,bins = np.histogram(anet[:], bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax1.plot(binsc[boff:],hist[boff:]/ntime, color = 'k')

    ax1.set_ylabel('Probability')
    ax1.set_xlabel('Anet [umol/m2/s]')
    ax1.grid('on')

    hist,bins = np.histogram(agross[:], bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax2.plot(binsc[boff:],hist[boff:]/ntime, color = 'k')
    ax2.set_xlabel('Ag [umol/m2/s]')
    ax2.grid('on')

    hist,bins = np.histogram(gstoma[:]*1.e-6, bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax3.plot(binsc[boff:],hist[boff:]/ntime, color = 'k')

    ax3.set_ylabel('Probability')
    ax3.set_xlabel('gs [mol/m2/s]')
    ax3.grid('on')

    hist,bins = np.histogram(lmr[:], bins=nhbins)
    binsc =  [0.5*(bins[i+1]+bins[i]) for i in range(len(bins)-1) ]
    ax4.plot(binsc[boff:],hist[boff:]/ntime, color = 'k')

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
