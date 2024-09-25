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

import ctypes
from ctypes import *
from operator import add
sys.path.append('../shared/py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 11}

matplotlib.rc('font', **font)


def GetModSymbol(mod_path,symbol):

    subprocess.run(["nm", mod_path],capture_output=True)
    
    list = str(subprocess.run(["nm",mod_path],capture_output=True)).split()
    p = re.compile(symbol, re.IGNORECASE)
    modlist = [ s for s in list if p.search(s)]

    # If nothing came up, thats bad
    if (len(modlist)==0):
        print('Failed to find the right module symbol for:{} in module {}'.format(symbol,mod_path))
        print(list)
        print('Exiting')
        exit(2)

    # Its possible a routine name could also be part of another routine name,
    # such as alloc_such_and_such  vs dealloc_such_and_such, we want the shorter
    mod_symbols = []
    slen = 10000
    for item in modlist:
        ms_str = item.split('\\')[0]
        if len(ms_str)<slen:
            mod_symbol = ms_str
            slen = len(ms_str)

        
    return mod_symbol

# ========================================================================

def GetJmaxKp25Top(vcmax25_top):

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
        

def MolarToVeloCF(press,tempk):

    umol_per_kmol = 1000.0
    rgas = 8314.4598
    
    cf = press/(rgas * tempk )*umol_per_kmol

    return cf


    
# Instantiate the F90 modules

f90_const_obj = ctypes.CDLL('bld/FatesConstantsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_shr_obj = ctypes.CDLL('bld/WrapShrMod.o',mode=ctypes.RTLD_GLOBAL)
f90_fatesutils_obj = ctypes.CDLL('bld/FatesUtilsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_leaf_biophys_obj = ctypes.CDLL('bld/LeafBiophysicsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_leaf_biophys_supp_obj = ctypes.CDLL('bld/LeafBiophysSuppMod.o',mode=ctypes.RTLD_GLOBAL)

# Identify subroutine objects in the helper code LeafBiophysSuppMod
f90_set_leaf_param_sub = getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','setleafparam'))
f90_alloc_leaf_param_sub = getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','allocleafparam'))
f90_dealloc_leaf_param_sub = getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','deallocleafparam'))
f90_dump_param_sub =  getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','DumpParams'))
f90_set_leaf_param_sub.argtypes = [POINTER(c_double),POINTER(c_int),c_char_p,c_long]
f90_biophysrate_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerBiophysicalRate'))
f90_biophysrate_sub.argtypes = [POINTER(c_int), \
                                POINTER(c_double),POINTER(c_double),POINTER(c_double), \
                                POINTER(c_double),POINTER(c_double),POINTER(c_double), \
                                POINTER(c_double),POINTER(c_double),POINTER(c_double), \
                                POINTER(c_double),POINTER(c_double),POINTER(c_double)]

f90_leaflayerphoto_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerPhotosynthesis'))
f90_leaflayerphoto_sub.argtypes = [POINTER(c_double),POINTER(c_double),POINTER(c_double), \
                                   POINTER(c_int), \
                                   POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double), \
                                   POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double), \
                                   POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double), \
                                   POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double), \
                                   POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double), \
                                   POINTER(c_double),POINTER(c_double)]

# Example of setting a function return call
#dftcdpsi_from_psi.restype = c_double

    
# For debugging
dump_parameters = False



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
    for param in pft_root.iter('param'):
        for pft in range(numpft):
            iret = f90_set_leaf_param_sub(c8(float(param.text.split(',')[pft])),ci(pft+1),*ccharnb(param.attrib['name'].strip()))
            if(param.attrib['name'].strip() == 'fates_leaf_c3psn'):
                leaf_c3psn.append(int(param.text.split(',')[pft]))
                
    # Dump parameters
    if(dump_parameters):
        iret = f90_dump_param_sub()

    
    # Read in non-fortran parameters from the xml file
    leafn_vert_scaler_coeff1 = []
    leafn_vert_scaler_coeff2 = []
    fates_leaf_vcmax25top    = []
    print('Reading non-fortran parameters')
    py_pft_root = xmlroot.find('py_params').find('pft_dim')
    for param in py_pft_root.iter('param'):
        for pft in range(numpft):
            if (param.attrib['name']=='fates_leafn_vert_scaler_coeff1'):
                leafn_vert_scaler_coeff1.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_leafn_vert_scaler_coeff2'):
                leafn_vert_scaler_coeff2.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_leaf_vcmax25top'):
                fates_leaf_vcmax25top.append(np.float64(param.text.split(',')[pft]))
                


    # Call the BiophysicalRates() Routine, which calculates
    # the actual vmax, jmax and rcurve-islope terms for the
    # leaf-scale conditions
    nscaler = 1.0
    tfrz_1atm = 273.15
    leaf_tempk = tfrz_1atm + 25.0  # Use 25 to start
    dayl_factor = 1.0
    t_growth_kum = -999
    t_home_kum = -999
    btran = 1.0

    # Leaf temperature ranges [C]
    leaf_tempk_min = -30.0
    leaf_tempk_max = 60.0
    leaf_tempk_n = 100
    leaf_tempk_vec = np.linspace(leaf_tempk_min,leaf_tempk_max,num=leaf_tempk_n) + tfrz_1atm

    # Absorbed PAR ranges [W/m2]
    par_abs_min = 0.01
    par_abs_max = 200
    par_abs_n  = 100
    par_abs_vec = np.linspace(par_abs_min,par_abs_max,num=par_abs_n)

    # Relative Humidity Ranges
    rh_max = 0.99
    rh_min = 0.01
    rh_n   = 100
    rh_vec = np.linspace(rh_min,rh_max,num=rh_n)

    # Generic boundary layer conductance from Kimura et al ~1.2 cm/s
    # reasonable ranges for green houses are 0-3 cm/s
    # perhaps storms would be >10cm/s?
    bl_cond = MolarToVeloCF(10100,298.0)*1.2
    
    
    # Initialize fortran return values
    vcmax_f  = c_double(-9.0)
    jmax_f   = c_double(-9.0)
    kp_f     = c_double(-9.0)
    agross_f = c_double(-9.0)
    gstoma_f = c_double(-9.0)
    anet_f   = c_double(-9.0)
    ac_f     = c_double(-9.0)
    aj_f     = c_double(-9.0)
    ap_f     = c_double(-9.0)
    veg_qs_f = c_double(-9.0)
    veg_es_f = c_double(-9.0)
    mm_kco2_f = c_double(-9.0)
    mm_ko2_f  = c_double(-9.0)
    co2_cpoint_f = c_double(-9.0)

    can_press_1atm = 101000.0
    co2_ppress_400ppm = 0.0004*can_press_1atm
    o2_ppress_209kppm = 0.29*can_press_1atm
    btran_nolim = 1.0
    
    for pft in range(numpft):

        jmax25_top,kp25_top =  GetJmaxKp25Top(fates_leaf_vcmax25top[pft])
        vcmax = np.zeros([leaf_tempk_n])
        jmax  = np.zeros([leaf_tempk_n])
        kp    = np.zeros([leaf_tempk_n])

        agross = np.zeros([leaf_tempk_n,par_abs_n])
        gstoma = np.zeros([leaf_tempk_n,par_abs_n])
        anet   = np.zeros([leaf_tempk_n,par_abs_n])
        ac     = np.zeros([leaf_tempk_n,par_abs_n])
        aj     = np.zeros([leaf_tempk_n,par_abs_n])
        ap     = np.zeros([leaf_tempk_n,par_abs_n])
        
        for it, leaf_tempk in enumerate(leaf_tempk_vec):
        
            iret = f90_biophysrate_sub(ci(pft+1), c8(fates_leaf_vcmax25top[pft]), c8(jmax25_top), c8(kp25_top), \
                                       c8(nscaler), c8(leaf_tempk), c8(dayl_factor), \
                                       c8(t_growth_kum),c8(t_home_kum),c8(btran), \
                                       vcmax_f, jmax_f, kp_f)
            vcmax[it] = vcmax_f.value
            jmax[it]  = jmax_f.value
            kp[it]    = kp_f.value

            iret = f90_qsat_sub(c8(leaf_tempk_f),c8(can_press_1atm),veg_qs_f,veg_es_f)

            iret = f90_cangas_params(c8(can_press_1atm), \
                                     c8(o2_ppress_209kppm), \
                                     c8(leaf_tempk_f), \
                                     mm_kco2_f, \
                                     mm_ko2_f, \
                                     co2_cpoint_f)

            for ir, rh in enumerate(rh_vec):
            
                vpress = rh * veg_es_f.value
            
                for ip, par_abs in enumerate(par_abs_vec):
            
                    iret = f90_leaflayerphoto_sub(c8(par_abs),  \
                                                  c8(1.0),     \
                                                  ci(pft+1),   \
                                                  vcmax[it],   \
                                                  jmax[it],    \
                                                  kp[it],      \
                                                  c8(leaf_tempk), \
                                                  c8(can_press_1atm), \
                                                  c8(co2_ppress_400ppm), \
                                                  c8(o2_ppress_209kppm), \
                                                  c8(btran_nolim), \
                                                  c8(bl_cond), \
                                                  c8(vpress), \
       mm_kco2,           &  ! in
       mm_ko2,            &  ! in
       co2_cpoint,        &  ! in
       lmr,               &  ! in
       leaf_psi,          &  ! in   (currently dummy)
       psn_out,           &  ! out
       gs_out,            &  ! out
       anet_out,          &  ! out
       c13disc_out,       &  ! out
       ac,                &  ! out
       aj,                &  ! out
       ap)                   ! out


            
        
        fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8.5,7.5))

        ap = ax1.plot(leaf_tempk_vec-tfrz_1atm,vcmax)
        
        ax1.set_ylabel('Vcmax [umol/m2/s]')
        ax1.set_title('PFT: {}'.format(pft+1))
        ax1.grid(True)

        ap = ax2.plot(leaf_tempk_vec-tfrz_1atm,jmax)
        ax2.set_xlabel('Leaf Temperature [C]')
        ax2.set_ylabel('jmax [umol/m2/s]')
        ax2.set_title('PFT: {}'.format(pft+1))
        ax2.grid(True)

        print(pft,leaf_c3psn[pft])
        
        if(leaf_c3psn[pft] == 0):
            ap = ax3.semilogy(leaf_tempk_vec-tfrz_1atm,kp)
            ax3.set_xlabel('Leaf Temperature [C]')
            ax3.set_ylabel('Kp (C4 only)')
            ax3.grid(True)
            ax3.set_xlabel('Leaf Temperature [C]')
        else:
            ax1.set_xlabel('Leaf Temperature [C]')
            ax3.axis("off")

        ax4.axis("off")



        




        
    plt.show()
        
    code.interact(local=dict(globals(), **locals()))

    # Environmental Conditions
    #   Absorbed PAR [w/m2 leaf]
    #   Leaf Temperature [K]
    #   Air Temperature  [K]
    #   Air Saturation Vapor Pressure at leaf surface []
    #   Air Vapor Pressure at leaf surface []

        
    #kn = decay_coeff_vcmax(currentCohort%vcmax25top, &
    #                       prt_params%leafn_vert_scaler_coeff1(ft), &
    #                       prt_params%leafn_vert_scaler_coeff2(ft))

    
    #StomatalCondMedlyn(anet,ft,veg_esat,can_vpress,stomatal_intercept_btran,leaf_co2_ppress,can_press,gb,gs)
        

        
    
    print('Deallocating parameter space')
    iret = f90_dealloc_leaf_param_sub()
    
    print('Functional Unit Testing Complete')
    exit(0)

    
# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
