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
import CtypesLeafBiophys

import ctypes
from ctypes import *
from operator import add
sys.path.append('../shared/py_src')

from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 11}

matplotlib.rc('font', **font)


# Global constants to use in all Leaf Biophysics unit testing
# =======================================================================================

do_comparekitajimacond = False


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

wm2_to_umolm2s = 4.6
 
# Simple conversion, number of micro-moles in a mole
umol_per_mol = 1.e6
mol_per_umol = 1.e-6

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


# Respiration does affect conductance, which also affects
# photosynthesis, but for now we use 0
zero_lmr = 0.0

# Set Leaf water potential to zero for some calcs
zero_lwp = 0.0

# When there is hydrualic limitation on photosynthesis
# (via Vcmax reductions), then the btran factor is 1
btran_nolimit = 1.0

# Respiration scaler at canopy top
rdark_scaler_top = 1.0

# Nitrogen scaler at canopy top
nscaler_top = 1.0

normalized_boundary = 1

leaf_rhonir = [0.46, 0.41, 0.39, 0.46, 0.41, 0.41, 0.46, 0.41, 0.41, 0.28, 0.28, 0.28 ]
leaf_rhovis = [0.11, 0.09, 0.08, 0.11, 0.08, 0.08, 0.11, 0.08, 0.08, 0.05, 0.05, 0.05 ]
leaf_taunir = [0.33, 0.32, 0.42, 0.33, 0.43, 0.43, 0.33, 0.43, 0.43, 0.4,  0.4,  0.4 ]
leaf_tauvis = [0.06, 0.04, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05]
leaf_xl     = [0.32, 0.01, 0.01, 0.32, 0.2, 0.59, 0.32, 0.59, 0.59, -0.23, -0.23, -0.23]
leaf_clumping_index = [0.85, 0.85, 0.8, 0.85, 0.85, 0.9, 0.85, 0.9, 0.9, 0.75, 0.75, 0.75]
stem_rhonir = [0.49, 0.36, 0.36, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.53, 0.53, 0.53]
stem_rhovis = [0.21, 0.12, 0.12, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.31, 0.31, 0.31]
stem_taunir = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.25, 0.25, 0.25]
stem_tauvis = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.12, 0.12, 0.12]


# Create aliases for the ctype Fortran objects
# =======================================================================================

exec(open("CtypesLeafBiophys.py").read())


# Subroutines
# =======================================================================================

def GetCosZ(ZONE,JULIAN,LAT,LON,LHOUR):

    #SUBROUTINE COSZENITH (LON,LATD,LHOUR,ZONE,JULIAN,CZENITH)
    #
    #     The purpose is to calculate the following:
    #        1)  Day angle (GAMMA)
    #        2)  Solar DEClination
    #        3)  Equation of time
    #        4)  Local apparent time
    #        5)  Hour angle
    #        6)  Cosine of zenith angle
    #
    #     All equations come from "An Introduction to
    #     Solar Radition" By Muhammad Iqbal, 1983.
    #
    # Integers:
    #     ZONE,                ! time zone (1-24) GMT=12
    #     JULIAN               ! julian day
    # Real
    #     CZENITH,             ! cosine of zenith angle (radians)
    #     DECd,                ! solar declination (degrees)
    #     DEC,                 ! solar declination (radians)
    #     et,                  ! equation of time (minutes)
    #     GAMMA,               ! day angle (radians)
    #     LATime,              ! local apparent time
    #     LCORR,               ! longitudical correction
    #     LHOUR,               ! local standard time
    #     LON,                 ! local longitude (deg)
    #     LLAT,                ! local latitude in radians
    #     LATD ,               ! local latitude in degrees
    #     LS,                  ! standard longitude (deg)
    #     OMEGAD,              ! omega in degrees
    #     OMEGA ,              ! omega in radians
    #     PI,                  ! universal constant PI [-]
    #     ZENITH,              ! zenith angle(radians)
    #     ZEND                 ! zenith angle(degress)
    #
    #     Neither ZENITH nor ZEND are necessary for this program.
    #     I originally used them as checks, and left them here in
    #     case anyone else had a use for them.
    #
    #     1)  Day angle GAMMA (radians) page 3

    PI= 3.141592              # universal constant PI
    GAMMA=2.*PI*(JULIAN-1)/365.

    #     2) Solar declination (assumed constant for a 24 hour period)  page 7
    #     in radians
    #
    DEC=(0.006918-0.399912*np.cos(GAMMA)+0.070257*np.sin(GAMMA) \
         -0.006758*np.cos(2.*GAMMA)+0.000907*np.sin(2.*GAMMA) \
         -0.002697*np.cos(3.*GAMMA)+0.00148*np.sin(3.*GAMMA))
    DECd=DEC*(180./PI)
    #
    #     maximum error 0.0006 rad (<3'), leads to error of less than 1/2 degree
    #     in ZENITH angle
    #     3)  Equation of time  page 11

    et=(0.000075+0.001868*np.cos(GAMMA)-0.032077*np.sin(GAMMA) \
        -0.014615*np.cos(2*GAMMA)-0.04089*np.sin(2*GAMMA))*229.18
    #
    #     4) Local apparent time  page 13
    #
    #     LS     standard longitude (nearest 15 degree meridian)
    #     LON     local longitude
    #     LHOUR  local standard time
    #     LATIME local apparent time
    #     LCORR  longitudunal correction (minutes)
    #
    
    LS=((ZONE-1)*15)-180.
    LCORR=4.*(LS-LON)*(-1.)
    LATIME=LHOUR+LCORR/60.+et/60.
    if (LATIME<0.):
        LATIME=LATIME+24.
    if (LATIME>24.):
        LATIME=LATIME-24.

    #     5) Hour angle OMEGA  page 15
    #
    #     hour angle is zero at noon, postive in the morning
    #     It ranges from 180 to -180
    #
    #     OMEGAD is OMEGA in degrees, OMEGA is in radians
    #
    OMEGAD=(LATIME-12.)*(-15.)
    OMEGA=OMEGAD*PI/180.
    #
    #     6)  Zenith angle page 15
    #
    #     CZENITH cosine of zenith angle (radians)
    #     LAT=local latitude in degrees
    #     LLAT=local latitude in radians
    #
    LLAT=LAT*PI/180.
    CZENITH=np.sin(DEC)*np.sin(LLAT)+np.cos(DEC)*np.cos(LLAT)*np.cos(OMEGA)
    CZENITH=np.max([0.,CZENITH])
    ZENITH=np.arcsin(CZENITH)
    ZEND=180.*ZENITH/PI
    
    return CZENITH


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


# Plot support routines
# ========================================================================

def ReadModoutput(filepath):

    print('Experiment 1: Evaluating Photosynthesis and Conductance For Site Specific Drivers')


    #df = pd.read_csv(filepath, delimiter=",", names=["yr,mon,day,tod,r_b,u_ref,forc_q,vpress,t_veg,can_press,forc_t,thm,t_can,Rbeam,Rdiff"], header=None)
    df = pd.read_csv(filepath, delimiter=",", header=None)
    #code.interact(local=dict(globals(), **locals()))

    yr = df[0].values
    mon = df[1].values
    day = df[2].values
    tod = df[3].values
    r_b = df[4].values
    u_ref = df[5].values
    forc_q = df[6].values
    vpress = df[7].values
    t_veg = df[8].values
    can_press = df[9].values
    forc_t = df[10].values
    thm  = df[11].values
    t_can = df[12].values
    Rbeam  = df[13].values
    Rdiff  = df[14].values

    return yr,mon,day,tod,r_b,u_ref,forc_q,vpress,t_veg,can_press,forc_t,thm,t_can,Rbeam,Rdiff

# =======================================================================================

def main(argv):


    # Load the xml control file
    xmlfile = "leaf_biophys_controls.xml"
    xmlroot = ET.parse(xmlfile).getroot()
    
    numpft = int(xmlroot.find('numpft').text.strip())

    # Allocating parameters
    print('Allocating parameter space for {} pfts'.format(numpft))
    iret = f90_alloc_leaf_param_sub(ci(numpft))

    # Push scalar parameters
    print('Pushing parameters from the xml file to the f90 lb_params datastructure')
    scalar_root = xmlroot.find('f90_params').find('scalar_dim')
    for param in scalar_root.iter('param'):
        iret = f90_set_leaf_param_sub(c8(float(param.text.split(',')[0])),ci(0),*ccharnb(param.attrib['name'].strip()))

    # Push pft parameters to fortran instantiations
    pft_root = xmlroot.find('f90_params').find('pft_dim')
    leaf_c3psn = []
    leaf_stomatal_intercept = []
    for param in pft_root.iter('param'):
        for pft in range(numpft):
            iret = f90_set_leaf_param_sub(c8(float(param.text.split(',')[pft])),ci(pft+1),*ccharnb(param.attrib['name'].strip()))
            if(param.attrib['name'].strip() == 'fates_leaf_c3psn'):
                leaf_c3psn.append(int(param.text.split(',')[pft]))
            if(param.attrib['name'].strip() == 'fates_leaf_stomatal_intercept'):
                leaf_stomatal_intercept.append(int(param.text.split(',')[pft]))
                
    # Dump parameters
    #iret = f90_dump_param_sub()
    
    # Read in non-fortran parameters from the xml file
    leafn_vert_scaler_coeff1 = []
    leafn_vert_scaler_coeff2 = []
    fates_leaf_vcmax25top    = []
    fates_stoich_nitr = []
    fates_leaf_slatop = []
    
    print('Reading non-fortran pft parameters')
    py_pft_root = xmlroot.find('py_params').find('pft_dim')
    for param in py_pft_root.iter('param'):
        for pft in range(numpft):
            if (param.attrib['name']=='fates_leafn_vert_scaler_coeff1'):
                leafn_vert_scaler_coeff1.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_leafn_vert_scaler_coeff2'):
                leafn_vert_scaler_coeff2.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_leaf_vcmax25top'):
                fates_leaf_vcmax25top.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_stoich_nitr'):
                fates_stoich_nitr.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_leaf_slatop'):
                fates_leaf_slatop.append(np.float64(param.text.split(',')[pft]))
                
    print('Reading non-fortran scalar parameters')
    py_scalar_root = xmlroot.find('py_params').find('scalar_dim')
    for param in py_scalar_root.iter('param'):
        if (param.attrib['name']=='fates_maintresp_leaf_model'):
            fates_maintresp_leaf_model = int(param.text.split(',')[0])

            
    # Drive conductance with Amax values from Kitajima et al for BTEs
    # Compare with stomatal conductance from Kitajima et al.
    if(do_comparekitajimacond):
        CompareKitajimaCond(0,fates_stoich_nitr[0],fates_leaf_slatop[0], \
                            leaf_stomatal_intercept[0],fates_maintresp_leaf_model)



        
    pft = 0

    [yr,mon,day,tod,r_b,u_ref,forc_q,vpress,t_veg,can_press,forc_t,thm,t_can,Rbeam,Rdiff] =  ReadModoutput("/home/rgknox/Downloads/bci_met.txt")

    hod = tod/3600.0

    # Create masks for wet season, dry season, morning and afternoon
    #
    # From: https://gml.noaa.gov/grad/solcalc/
    # Solar noon Nov 1 2011:  17:02:57 GMT
    #            Dec 31 2011: 17:22:06 GMT
    #            Feb 1 2013:  17:32:58 GMT
    #            Mar 31 2013: 17:23:33 GMT

    # Sunrise: 11:17

    mod_wet_ssn  = [ yr[iyr] == 2011 and (mon[iyr]==11 or mon[iyr]==12) for iyr,year in enumerate(yr)]
    mod_dry_ssn  = [ yr[iyr] == 2013 and (mon[iyr]==2 or mon[iyr]==3) for iyr,year in enumerate(yr)]

    morning = [ (hod[iyr] <= 17.25 and hod[iyr] > 8) and (Rbeam[iyr]>0. or Rdiff[iyr]>0.) for iyr,year in enumerate(yr)]
    afternoon = [ (hod[iyr] > 17.25 or hod[iyr] < 8) and (Rbeam[iyr]>0. or Rdiff[iyr]>0.) for iyr,year in enumerate(yr)]

    mod_wet_morn = [ (yr[iyr] == 2011 and (mon[iyr]==11 or mon[iyr]==12)) and morning[iyr] for iyr,year in enumerate(yr)]
    mod_wet_aftn = [ (yr[iyr] == 2011 and (mon[iyr]==11 or mon[iyr]==12)) and afternoon[iyr] for iyr,year in enumerate(yr)]

    mod_dry_morn  = [ (yr[iyr] == 2013 and (mon[iyr]==2 or mon[iyr]==3)) and morning[iyr] for iyr,year in enumerate(yr)]
    mod_dry_aftn  = [ (yr[iyr] == 2013 and (mon[iyr]==2 or mon[iyr]==3)) and afternoon[iyr] for iyr,year in enumerate(yr)]

    
    
    numpts = np.sum(mod_dry_morn)
    
    yr = yr[mod_dry_morn]
    hod = hod[mod_dry_morn]
    mon = mon[mod_dry_morn]
    day = day[mod_dry_morn]
    t_can = t_can[mod_dry_morn]
    t_veg = t_veg[mod_dry_morn]
    can_press = can_press[mod_dry_morn]
    vpress = vpress[mod_dry_morn]
    Rbeam = Rbeam[mod_dry_morn]
    Rdiff = Rdiff[mod_dry_morn]
    r_b = r_b[mod_dry_morn]

   

    
    qsat_f        = c_double(-9)
    esat_f        = c_double(-9)
    qsdt_f        = c_double(-9)
    esdt_f        = c_double(-9)
    vcmax_f       = c_double(-9)
    jmax_f        = c_double(-9)
    kp_f          = c_double(-9)
    agross_f      = c_double(-9)
    gstoma_f      = c_double(-9)
    anet_f        = c_double(-9)
    lmr_f         = c_double(-9)
    c13_f         = c_double(-9)
    ac_f          = c_double(-9)
    aj_f          = c_double(-9)
    ap_f          = c_double(-9)
    co2_interc_f  = c_double(-9)
    mm_kco2_f     = c_double(-9)
    mm_ko2_f      = c_double(-9)
    co2_cpoint_f  = c_double(-9)
    gs0_f         = c_double(-9)
    gs1_f         = c_double(-9)
    gs2_f         = c_double(-9)
    solve_inter_f = c_double(-9)

    cd_albedo_beam = c_double(-9.0)
    cd_albedo_diff = c_double(-9.0)
    cd_canabs_beam = c_double(-9.0)
    cd_canabs_diff = c_double(-9.0)
    cd_ffbeam_beam = c_double(-9.0)
    cd_ffdiff_beam = c_double(-9.0)
    cd_ffdiff_diff = c_double(-9.0)
    cd_rd_abs_leaf = c_double(-9.0)
    cd_rb_abs_leaf = c_double(-9.0)
    cd_r_abs_stem  = c_double(-9.0)
    cd_r_abs_snow  = c_double(-9.0)
    cd_leaf_sun_frac = c_double(-9.0)
    
    # We calculate solutions for every independant output from the FATES model
    # so initialize output vectors with the same size (use r_b or any variable)

    #numpts = r_b.size
    
    g_b_umol   = np.zeros(numpts)
    qsat       = np.zeros(numpts)
    vpress     = np.zeros(numpts)
    lmr        = np.zeros(numpts)
    agross     = np.zeros(numpts)
    gstoma     = np.zeros(numpts)
    anet       = np.zeros(numpts)
    ac         = np.zeros(numpts)
    aj         = np.zeros(numpts)
    ap         = np.zeros(numpts)
    co2_interc = np.zeros(numpts)
    vcmax      = np.zeros(numpts)
    jmax       = np.zeros(numpts)
    kp         = np.zeros(numpts)
    rh         = np.zeros(numpts)
    par_abs    = np.zeros(numpts)
    cosz       = np.zeros(numpts)

    jmax25_top,kp25_top =  GetJmaxKp25Top(fates_leaf_vcmax25top[pft])


    
    # Leaf Nitrogen Concentration at the canopy top
    lnc_top  = fates_stoich_nitr[pft]/fates_leaf_slatop[pft]


    # Create a vegetation canopy and initialize scattering elements

    total_lai = 5.0
    total_sai = 0.1*total_lai
    n_elem   = 1
    n_col     = 1
    area      = 1.0  # Assume a completely closed canopy
    ground_albedo_diff = 0.3
    ground_albedo_beam = 0.3
    frac_snow  = 0.0

    iret = alloc_twostream_call(ci(n_elem),ci(n_col))

    
    pft   = 1
    ft    = 0
    n_pft = 1
    n_bands = 2
    iret = alloc_radparams_call(ci(n_pft),ci(n_bands))
    # rho (vis+nir)
    iret = set_radparams_call(c_double(leaf_rhovis[ft]),c_int(pft),c_int(visb),*ccharnb("rhol"))
    iret = set_radparams_call(c_double(leaf_rhonir[ft]),c_int(pft),c_int(nirb),*ccharnb("rhol"))
    iret = set_radparams_call(c_double(stem_rhovis[ft]),c_int(pft),c_int(visb),*ccharnb("rhos"))
    iret = set_radparams_call(c_double(stem_rhonir[ft]),c_int(pft),c_int(nirb),*ccharnb("rhos"))
    # tau (vis+nir)
    iret = set_radparams_call(c_double(leaf_tauvis[ft]),c_int(pft),c_int(visb),*ccharnb("taul"))
    iret = set_radparams_call(c_double(leaf_taunir[ft]),c_int(pft),c_int(nirb),*ccharnb("taul"))
    iret = set_radparams_call(c_double(stem_tauvis[ft]),c_int(pft),c_int(visb),*ccharnb("taus"))
    iret = set_radparams_call(c_double(stem_taunir[ft]),c_int(pft),c_int(nirb),*ccharnb("taus"))
    # orientations
    iret = set_radparams_call(c_double(leaf_xl[ft]),c_int(pft),c_int(0),*ccharnb("xl"))
    iret = set_radparams_call(c_double(leaf_clumping_index[ft]),c_int(pft),c_int(0),*ccharnb("clumping_index"))

    iret = param_prep_call(ci(n_pft))
    
    iret = setup_canopy_call(c_int(1),c_int(1),c_int(pft),c_double(area),c_double(total_lai),c_double(total_sai))
    iret = grndsnow_albedo_call(c_int(visb),c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = grndsnow_albedo_call(c_int(visb),c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))
    iret = grndsnow_albedo_call(c_int(nirb),c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = grndsnow_albedo_call(c_int(nirb),c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))
    iret = canopy_prep_call(c8(frac_snow))

    print('canopy prep done')

    # variability of leaf temp, dry-morn: 0.33776200082982527
    # variability of leaf temp, dry-aftn: 0.16511968767353083

    # For cosine..
    zone = 12-5
    doys = [0,31,28,31,30,31,30,31,31,30,31,30]
    cdoys = np.cumsum(doys)
    latitude=9.153
    longitude=280.1539-360.0

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
    
    for it, tvegk in enumerate(t_veg):

        doy = cdoys[mon[it]-1] + day[it]
        lhour = hod[it]-5.0
        if(lhour<0):
            lhour = lhour-24
            doy   = doy-1
            
        cosz[it] = GetCosZ(zone,doy,latitude,longitude,lhour)

        print("cosz: {}, Rbeam: {}, Rdiff: {}".format(cosz[it],Rbeam[it],Rdiff[it]))
        iret = zenith_prep_call(c8(cosz[it]))
        iret = solver_call(ci(visb),ci(normalized_boundary),c8(1.0),c8(1.0), \
                       byref(cd_albedo_beam),byref(cd_albedo_diff), \
                       byref(cd_canabs_beam),byref(cd_canabs_diff), \
                       byref(cd_ffbeam_beam),byref(cd_ffdiff_beam),byref(cd_ffdiff_diff))
        iret = setdown_call(ci(visb),c8(Rbeam[it]),c8(Rdiff[it]))




        
        g_b_umol[it] = f90_velotomolarcf_sub(c8(can_press[it]),c8(tvegk))/r_b[it]
        
        # Get Humidity
        iret = f90_qsat_sub(c8(t_can[it]),c8(can_press[it]), \
                            byref(qsat_f),byref(esat_f), \
                            byref(qsdt_f),byref(esdt_f))

        #print("humidity done")
        
        qsat[it] = qsat_f.value
        
        #vpress[it] = can_press[it]/ ( eps/q[it] + 1 - eps)
        #vpress[it] = can_press[it]*q[it]/eps
        
        rh[it] = 100.*vpress[it]/esat_f.value

        iret = f90_cangas_sub(c8(can_press[it]), \
                              c8(o2_ppress_209kppm), \
                              c8(t_can[it]), \
                              byref(mm_kco2_f), \
                              byref(mm_ko2_f), \
                              byref(co2_cpoint_f))

        #print("mm and co2point done")

        for il in range(n_layer):
            #iret = getintens_call(ci(1),ci(1),ci(visb),c8(avai[il]),byref(cd_r_diff_dn),byref(cd_r_diff_up),byref(cd_r_beam))
            iret = getabsrad_call(ci(1),ci(1),ci(visb),c8(avai[il]),c8(avai[il]+davai), \
                                  byref(cd_rd_abs_leaf),byref(cd_rb_abs_leaf),byref(cd_r_abs_stem), \
                                  byref(cd_r_abs_snow),byref(cd_leaf_sun_frac))
            rd_abs_leaf[il] = cd_rd_abs_leaf.value
            rb_abs_leaf[il] = cd_rb_abs_leaf.value
            sunfrac[il]     = cd_leaf_sun_frac.value
            par_abs_umol_m2 = (rd_abs_leaf[il] + rb_abs_leaf[il])*wm2_to_umolm2s/dalai

            # Scale down N and biophysical rates
            #nscaler = 
            iret = f90_biophysrate_sub(ci(pft+1), c8(fates_leaf_vcmax25top[pft]), \
                                       c8(jmax25_top), c8(kp25_top), \
                                       c8(nscaler_top), c8(tvegk), c8(dayl_factor_full), \
                                       c8(t_growth_kum),c8(t_home_kum),c8(btran_nolimit), \
                                       byref(vcmax_f), byref(jmax_f), byref(kp_f), \
                                       byref(gs0_f), byref(gs1_f), byref(gs2_f))

            # Leaf Maintenance Respiration (temp and pft dependent)
            if(fates_maintresp_leaf_model==1):
                iret = f90_lmr_ryan_sub(c8(lnc_top),c8(nscaler_top), ci(pft+1), c8(tvegk), byref(lmr_f))
            elif(fates_maintresp_leaf_model==2):
                iret = f90_lmr_atkin_sub(c8(lnc_top),c8(rdark_scaler),c8(tvegk),c8(atkin_mean_leaf_tempk),byref(lmr_f) )
            else:
                print('unknown leaf respiration model')
                exit(1)

            
            iret = f90_leaflayerphoto_sub(c8(par_abs_umol_m2), \
                                          c8(par_abs_umol_m2),  \
                                          c8(1.0),     \
                                          ci(pft+1),   \
                                          c8(vcmax_f.value),   \
                                          c8(jmax_f.value),    \
                                          c8(kp_f.value),      \
                                          c8(gs0_f.value), \
                                          c8(gs1_f.value), \
                                          c8(gs2_f.value), \
                                          c8(tvegk), \
                                          c8(can_press[it]), \
                                          c8(co2_ppress_400ppm), \
                                          c8(o2_ppress_209kppm), \
                                          c8(g_b_umol[it]), \
                                          c8(vpress[it]), \
                                          c8(mm_kco2_f.value), \
                                          c8(mm_ko2_f.value), \
                                          c8(co2_cpoint_f.value), \
                                          c8(lmr_f.value), \
                                          c8(0.5), \
                                          byref(agross_f), \
                                          byref(gstoma_f), \
                                          byref(anet_f), \
                                          byref(c13_f), \
                                          byref(ac_f), \
                                          byref(aj_f), \
                                          byref(ap_f), \
                                          byref(co2_interc_f), \
                                          byref(solve_inter_f))
            
            agross[it]     = agross[it]+vfrac*agross_f.value
            gstoma[it]     = gstoma[it]+vfrac*gstoma_f.value
            anet[it]       = anet[it]+vfrac*anet_f.value
            ac[it]         = ac[it]+vfrac* ac_f.value
            aj[it]         = aj[it]+vfrac*aj_f.value
            ap[it]         = ap[it]+vfrac*ap_f.value
            co2_interc[it] = co2_interc[it]+vfrac*co2_interc_f.value
        
            
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


        

        #print("photo done")
        
        
        
    fig,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(6.5,8.5))

    
    ax1.scatter(t_veg-tfrz_1atm,agross)
    ax1.set_xlabel('Tveg [C]')
    ax1.set_ylabel('Ag [umol m-2 s-1]')
    ax1.grid(True)

    ax2.scatter(t_veg-tfrz_1atm,gstoma*mol_per_umol)
    ax2.set_xlabel('Tveg [C]')
    ax2.set_ylabel('gs [mol m-2 s-1]')
   
    ax2.grid(True)
    
    ax3.scatter(g_b_umol*mol_per_umol,agross)
    ax3.set_xlabel('gb [mol m-2 s-1]')
    ax3.set_ylabel('Ag [umol m-2 s-1]')
    ax3.grid(True)

    ax4.scatter(g_b_umol*mol_per_umol,gstoma*mol_per_umol)
    ax4.set_xlabel('gb [mol m-2 s-1]')
    ax4.set_ylabel('gs [mol m-2 s-1]')
    ax4.grid(True)

    ax5.scatter(par_abs,agross)
    ax5.set_xlabel('Par Abs. [W/m2]')
    ax5.set_ylabel('Ag [umol m-2 s-1]')
    ax5.grid(True)

    ax6.scatter(hod,cosz)
    ax6.set_xlabel('hour')
    ax6.set_ylabel('cosz')
    ax6.grid(True)

    
    # Plot histograms of the driver data
    fig,((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8)) = plt.subplots(4,2,figsize=(5.5,9.5))

    boxprops = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # Histogram of g_b_mol

    ax1.hist(g_b_umol*mol_per_umol, bins = 50)
    ax1.set_ylabel('count')
    ax1.set_xlabel('gb [mol H20 m-2 s-1]')
    ax1.grid(True)
    ax1.text(0.25, 0.95, 'min(gb)=%5.3f'%(np.min(g_b_umol*mol_per_umol)), transform=ax1.transAxes, fontsize=12,
             verticalalignment='top', bbox=boxprops)

    # Histogram of vegetation temperature

    ax2.hist(t_veg-tfrz_1atm, bins = 50)
    ax2.set_ylabel('count')
    ax2.set_xlabel('Tv [C]')
    ax2.grid(True)
    
    # Histogram of wind-speed

    ax3.hist(u_ref, bins = 50)
    ax3.set_ylabel('count')
    ax3.set_xlabel('U ref [m/s]')
    ax3.grid(True)
    # place a text box in upper left in axes coords
    ax3.text(0.25, 0.95, 'min(Uref)=%3.1f'%(np.min(u_ref)), transform=ax3.transAxes, fontsize=12,
             verticalalignment='top', bbox=boxprops)

    # Histogram of relative humidity
    ax4.hist(rh, bins = 50)
    ax4.set_ylabel('count')
    ax4.set_xlabel('RH [%]')
    ax4.grid(True)
    ax4.text(0.25, 0.95, 'max(RH)=%3.1f'%(np.max(rh)), transform=ax4.transAxes, fontsize=12,
             verticalalignment='top', bbox=boxprops)

    
    ax5.hist(vcmax, bins = 50)
    ax5.set_ylabel('count')
    ax5.set_xlabel('Vcmax [umol/m2/s]')
    ax5.grid(True)
    
    ax6.hist(jmax, bins = 50)
    ax6.set_ylabel('count')
    ax6.set_xlabel('Jmax [umol/m2/s]')
    ax6.grid(True)
    
    ax7.scatter(t_veg-tfrz_1atm,vcmax)
    ax7.set_xlabel('Tv [C]')
    ax7.set_ylabel('Vcmax [umol/m2/s]')
    ax7.grid(True)
    
    ax8.scatter(t_veg-tfrz_1atm,jmax)
    ax8.set_xlabel('Tv [C]')
    ax8.set_ylabel('Jmax [umol/m2/s]')
    ax8.grid(True)
    #plt.tight_layout()
    
    plt.show()

    
    print('Deallocating parameter space')
    iret = f90_dealloc_leaf_param_sub()
    
    print('Functional Unit Testing Complete')
    exit(0)


    
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
