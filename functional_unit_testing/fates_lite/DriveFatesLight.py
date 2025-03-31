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

from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

from MetDrivers import GetCosZ
from MetDrivers import met_driver
from CrownPhysics import CanopyElementPhysics
 
from PushParameters import PushParameters
from PushParameters import PushXMLParameters
from PushParameters import GetParamFromAttrib
from PushParameters import GetParamList

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

# Nitrogen scaler at canopy top
nscaler_top = 1.0

# This signals that we perform a solution for a unit raditation calculation
# and then afterwards, scale by the magnitude of the downwelling flux
normalized_boundary = 1

# Initialize all of the fortran shared objects, and create aliases
# for the various subroutines and functions
# =======================================================================================

exec(open("../shared/py_src/CtypesInit.py").read())

n_layer = 10

# Subroutines
# =======================================================================================

# Radiation diagnostics
class elem_diags_type:

    def __init__(self,total_lai,total_sai,n_layer,ntimes):

        self.total_lai = total_lai
        self.total_sai = total_sai
        self.n_layer = n_layer
        self.avai = np.zeros(n_layer)  # Accumulated VAI (top of bin)
        self.rd_abs_leaf = np.zeros(n_layer)
        self.rb_abs_leaf = np.zeros(n_layer)
        self.r_abs_stem  = np.zeros(n_layer)
        self.rd_dn = np.zeros(n_layer)
        self.rd_up = np.zeros(n_layer)
        self.rbeam =  np.zeros(n_layer)
        self.sunfrac = np.zeros(n_layer)
        self.sunfrac_v2 = np.zeros(n_layer)
        for il in range(n_layer):
            self.avai[il] = (total_lai+total_sai)*float(il)/float(n_layer)
        self.davai =  (total_lai+total_sai)/float(n_layer)
        self.dalai =  total_lai/float(n_layer)
        self.vfrac = 1./float(n_layer)
        self.ag_limit = np.zeros([n_layer,2])
        self.ag_sslimit = np.zeros([n_layer,2,2])
        return
    
    def ZeroDiag(self):

        # Zero out the instantaneous diagnostics that
        # are not indexed by time.
        self.rd_abs_leaf = np.zeros(n_layer)
        self.rb_abs_leaf = np.zeros(n_layer)
        self.r_abs_stem  = np.zeros(n_layer)
        self.rd_dn = np.zeros(n_layer)
        self.rd_up = np.zeros(n_layer)
        self.rbeam =  np.zeros(n_layer)
        self.sunfrac = np.zeros(n_layer)
        self.sunfrac_v2 = np.zeros(n_layer)
        self.ag_limit = np.zeros([n_layer,2])
        self.ag_sslimit = np.zeros([n_layer,2,2])

        return

class site_diags_type:

    def __init__(self,ntimes):

        self.ntimes = ntimes
        self.aglimit_which = []
        self.aglimit_apar  = []
        self.aglimit_temp  = []
        self.lmr = np.zeros(ntimes)
        self.agross = np.zeros(ntimes)
        self.gstoma = np.zeros(ntimes)
        self.anet = np.zeros(ntimes)
        
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

    # Maintenance respiration model # 1=Ryan (1991), 2=Atkin et al (2017)
    maintresp_leaf_model = float(GetParamFromAttrib(scalar_root,'fates_maintresp_leaf_model')[0])


    # Call this external to push the default parameters to the F90 objects
    #PushParameters(f90,params,dims)
    PushXMLParameters(f90,xmlroot)

    # Set some plant trait data
    # -----------------------------------------------------------------------------------
    # We calculate solutions for every independant output from the FATES model
    # so initialize output vectors with the same size (use r_b or any variable)

    # params['fates_leaf_vcmax25top'].data.shape = (1, 14) Fist index is "age"
    #vcmax25_top = float(params['fates_leaf_vcmax25top'].data[0,ft])
    vcmax25_top = np.zeros(numpft)
    jmax25_top = np.zeros(numpft)
    kp25_top = np.zeros(numpft)
    leaf_slatop = np.zeros(numpft)
    lnc_top = np.zeros(numpft)
    for ft in range(numpft):
        vcmax25_top[ft] = float(GetParamFromAttrib(pft_root,'fates_leaf_vcmax25top')[ft])
    
        jmax25_top[ft],kp25_top[ft] =  GetJmaxKp25Top(vcmax25_top[ft])
        
        # Leaf Nitrogen Concentration at the canopy top (N:C ratio) (n/m2]
        # params['fates_stoich_nitr'].data.shape = (4, 14) First index is organ (leaf)
        # params['fates_stoich_nitr'].meta['units'] = 'gN/gC'
        #leaf_nc_ratio = params['fates_stoich_nitr'].data[0,ft]
        leaf_nc_ratio = float(GetParamFromAttrib(pft_root,'fates_stoich_nitr')[ft])

        #params['fates_leaf_slatop'].meta['units'] = 'm^2/gC'
        #leaf_slatop   = params['fates_leaf_slatop'].data[ft]
        leaf_slatop[ft]   = float(GetParamFromAttrib(pft_root,'fates_leaf_slatop')[ft])
    
        # Leaf N Conc at the crown top [gN/m2]
        lnc_top[ft]       = leaf_nc_ratio/leaf_slatop[ft]


    
    # code.interact(local=dict(globals(), **locals()))
    # Lets create a synthetic met driver based on nearest neighbor data?
    # ------------------------------------------------------------------------------------

    met_driver_csvfile = "bci_met_data/BCI_met_drivers_2003_2016.csv"

    met = met_driver(met_driver_csvfile)
    met.FilterTimes('daytime')

    
    ntimes = met.ndata

    site_diags = site_diags_type(ntimes)
    
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


    # Create a canopy based on input from the control file
    # -----------------------------------------------------------------------------------
    
    cstruct_root = xmlroot.find('canopy_structure')
    cohort_pft  = GetParamList(cstruct_root,'cohort_pft','integer')
    cohort_area = GetParamList(cstruct_root,'cohort_area_frac','float')
    cohort_lai  = GetParamList(cstruct_root,'cohort_lai','float')
    cohort_sai  = GetParamList(cstruct_root,'cohort_sai','float')
    ncohorts    = len(cohort_pft)
    site_root   = xmlroot.find('site_structure')
    ground_vis_albedo = GetParamList(site_root,'ground_vis_albedo','float')
    ground_nir_albedo = GetParamList(site_root,'ground_nir_albedo','float')
    frac_snow = GetParamList(site_root,'frac_snow','float')[0]
    n_layer = GetParamList(site_root,'nlayers','integer')[0] 

    
    # Lets walk through the list and assign the layer of each
    cohort_can    = list(range(ncohorts)) # The canopy layer associated with each cohort (0 indexed)
    cohort_col    = list(range(ncohorts)) # Column index of the cohort                   (0 indexed)
    layer_ncohort = [0]                   # Number of cohorts (elements) in each layer
    area_current  = 0
    ican = 0
    icol = 0
    for ico in range(ncohorts):

        if( cohort_area[ico]-1.0 < 0.001 ):
            cohort_area[ico] = np.min([1.0,cohort_area[ico]])
        else:
            print('\n\nYou specified a cohort that takes up more area than 1.0\n')
            print('Dont do that. Exiting.\n')
            exit(1)

        if(cohort_area[ico]+area_current>1.0):
            layer_ncohort.append(0)
            ican = ican+1
            area_current  = 0.
            icol = 0
        else:
            area_current = area_current + cohort_area[ico]
            icol = icol + 1
                
        cohort_can[ico]     = ican
        layer_ncohort[ican] = layer_ncohort[ican]+1
        cohort_col[ico]     = icol-1

    # Determine the size of the canopy scattering matrix
    n_can  = len(layer_ncohort)
    n_col  = np.max(layer_ncohort)

    # Setup the cohort index matrix, initialize non-cohorts
    # as air, with pft -1
    elem_cohort = np.zeros([n_can,n_col],dtype=int)-1
    for ico in range(ncohorts):
        ican  = cohort_can[ico]
        icol  = cohort_col[ico]
        print("cohort ican: {} icol: {}".format(ican,icol))
        elem_cohort[ican,icol] = ico

    

    # Initialize scattering elements

    elem_diags = []
    iret = f90.alloc_twostream_sub(ci(n_can),ci(n_col))
    iret = f90.param_prep_sub()  # This routine creates parameters that are derived from others
    for ican in range(n_can):
        veg_area = 0.0
        n_veg    = 0

        elem_diags.append(list())
        for icol in range(n_col):
            ico = elem_cohort[ican,icol]
            print("column: ",ico)
            if(ico>=0):
                print("adding cohort")
                veg_area = veg_area+cohort_area[ico]
                n_veg = n_veg + 1
                iret = f90.setup_canopy_sub(c_int(ican+1),c_int(icol+1), \
                                            c_int(cohort_pft[ico]), c_double(cohort_area[ico]), \
                                            c_double(cohort_lai[ico]), c_double(cohort_sai[ico]))
                elem_diags[ican].append(elem_diags_type(cohort_lai[ico],cohort_sai[ico],n_layer,ntimes))
            else:
                air_pft  = 0
                air_area = (1.-veg_area)/float(n_col-n_veg)
                air_lai  = 1.0 # notional
                air_sai  = 1.0 # notional
                iret = f90.setup_canopy_sub(c_int(ican+1),c_int(icol+1),c_int(air_pft), \
                                            c_double(air_area),c_double(air_lai),c_double(air_sai))

    # Site level scattering parameters
    iret = f90.grndsnow_albedo_sub(c_int(visb),c_double(ground_vis_albedo[1]),*ccharnb('albedo_grnd_diff'))
    iret = f90.grndsnow_albedo_sub(c_int(visb),c_double(ground_vis_albedo[0]),*ccharnb('albedo_grnd_beam'))
    iret = f90.grndsnow_albedo_sub(c_int(nirb),c_double(ground_nir_albedo[1]),*ccharnb('albedo_grnd_diff'))
    iret = f90.grndsnow_albedo_sub(c_int(nirb),c_double(ground_nir_albedo[0]),*ccharnb('albedo_grnd_beam'))
    iret = f90.canopy_prep_sub(c8(frac_snow))

    visualize_elements = False
    if(visualize_elements):
        fig10, ax = plt.subplots(ncols=1,nrows=1,figsize=(7,7))
        maxvai = 0.
        vai0   = 0.
        for ican in range(n_can):
            veg_area = 0.0
            vai0     = vai0+maxvai
            maxvai   = 0.0
            for icol in range(n_col):
                ico = elem_cohort[ican,icol]
                if(ico>=0):
                    vai = cohort_lai[ico]+cohort_sai[ico]
                    # xy, width height
                    rect = patches.Rectangle((veg_area, vai0), cohort_area[ico], \
                                             vai, linewidth=1, edgecolor='k', facecolor='g')
                    ax.add_patch(rect)
                    ax.text(veg_area+0.5*cohort_area[ico], vai0+0.5*vai, \
                            'cohort {}'.format(ico+1), horizontalalignment='center', \
                            verticalalignment='center')
                    maxvai = np.max([maxvai,vai])
                    veg_area = veg_area + cohort_area[ico]
                    
                    
        plot_vai = 1.2*(vai0+maxvai)
        ax.set_xlabel('Area Fraction [m2/m2]')
        ax.set_ylabel('Integrated LAI+SAI [m2/m2]')
        ax.set_ylim([0,plot_vai])
        ax.set_xlim([0,1])
        ax.invert_yaxis()
        plt.show()

    # Initialize variables that are output from the f90 modules
    # these need to be defined as a cytype, and are floating point,
    # integer or strings
    
    qsat_f          = c_double(-9);    esat_f          = c_double(-9);    qsdt_f          = c_double(-9)
    esdt_f          = c_double(-9);   
    albedo_beam_f   = c_double(-9)
    albedo_diff_f   = c_double(-9);    canabs_beam_f   = c_double(-9);    canabs_diff_f   = c_double(-9)
    ffbeam_beam_f   = c_double(-9);    ffdiff_beam_f   = c_double(-9);    ffdiff_diff_f   = c_double(-9)
    
    # Initialize output arrays
    # -----------------------------------------------------------------------------------

   

    cosz = np.zeros(ntimes)
    g_b_umol = np.zeros(ntimes)
    
    
    # Start the main model time loop
    # -----------------------------------------------------------------------------------
    
    for it in range(ntimes):

        ft = 0
        pft = 1
        
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
                ico = elem_cohort[ican,icol]
                if(ico>=0):
                    CanopyElementPhysics(ican,icol,ico,f90,met,it,xmlroot,cohort_pft[ico], \
                                         vcmax25_top[ft],jmax25_top[ft],kp25_top[ft],lnc_top[ft], \
                                         co2_ppress_400ppm, o2_ppress_209kppm, btran_nolimit, dayl_factor_full, \
                                         g_b_umol,maintresp_leaf_model,site_diags,elem_diags[ican][icol])
        

    # Per Element

    for ican in range(n_can):
            for icol in range(n_col):
                ico = elem_cohort[ican,icol]
                if(ico>=0):
                    ft = cohort_pft[ico]
                    fig55, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(6.5,6.5))
    
                    ax3.plot(elem_diags[ican][icol].ag_limit[:,1]/(elem_diags[ican][icol].ag_limit[:,0]+elem_diags[ican][icol].ag_limit[:,1]),elem_diags[ican][icol].avai[:])
                    ax3.invert_yaxis()
                    ax3.set_ylabel('LAI')
                    ax3.set_xlabel('Fraction RuBP Limited')
                    ax3.set_title('All')
                    ax3.grid('on')
    
                    ax1.plot(elem_diags[ican][icol].ag_sslimit[:,0,1]/(elem_diags[ican][icol].ag_sslimit[:,0,0]+elem_diags[ican][icol].ag_sslimit[:,0,1]),elem_diags[ican][icol].avai[:])
                    ax1.invert_yaxis()
                    ax1.set_ylabel('LAI')
                    ax1.set_xlabel('Fraction RuBP Limited')
                    ax1.set_title('Sunlit')
                    ax1.grid('on')
                    
                    ax2.plot(elem_diags[ican][icol].ag_sslimit[:,1,1]/(elem_diags[ican][icol].ag_sslimit[:,1,0]+elem_diags[ican][icol].ag_sslimit[:,1,1]),elem_diags[ican][icol].avai[:])
                    ax2.invert_yaxis()
                    ax2.set_ylabel('LAI')
                    ax2.set_xlabel('Fraction RuBP Limited')
                    ax2.set_title('Shaded')
                    ax2.grid('on')
                    plt.show()

                        
    fig7,(ax1,ax2) = plt.subplots(2,1,figsize=(5.5,7.5))

    ax1.scatter(site_diags.aglimit_apar,site_diags.aglimit_which)
    ax1.set_ylabel('1 = Rubisco, 0 = RuBP')
    ax1.set_xlabel('Apar [umol/m2/s]')
    ax1.grid('on')
    
    ax2.scatter(site_diags.aglimit_temp,site_diags.aglimit_which)
    ax2.set_ylabel('1 = Rubisco, 0 = RuBP')
    ax2.set_xlabel('Temperature [K]')
    ax2.grid('on')
           
    

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
