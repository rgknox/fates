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


# Global constants to use in all Leaf Biophysics unit testing
# =======================================================================================


# For debugging
dump_parameters = False

# Should we evaluate vcmax, jmax and kp actual?
do_evalvjkbytemp = True

# Should we run the comparison against Kitajima's plots?
do_comparekitajimacond = False

# Should we run the parameter sensitivity experiment
do_paramsense = True

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

# Simple conversion, number of micro-moles in a mole
umol_per_mol = 1.e6

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
    
# Subroutines
# =======================================================================================

def GetModSymbol(mod_path,symbol):

    # This routine uses some minor string magic and a system call using
    # "nm" to list the symbols inside a fortran object.  Symbols are simply
    # the text strings the define the routines and data structures inside
    # the compiled objects. When symbols are created inside a compiled object,
    # different compilers do different things, like force to lower or upper
    # case, or attach trailing or preceding underscores.  This routine
    # helps identify the exact symbol name, from the name of the fortran
    # routine (as it appears in code) that we want to match with.
    
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


# Plot support routines
# ========================================================================

def EvalVJKByTemp(numpft,fates_leaf_vcmax25top,leaf_c3psn):

    print('\n')
    print('Experiment 1: Evaluating Vcmax,Jmax,Kp by Temperature')

    # Plot out vcmax, jmax and kp as a function of temperature
    # Assumes canopy top position, and is dependent on the
    # PFT base rate
        
    # Leaf temperature ranges [C]
    leaf_tempc_min = -20.0
    leaf_tempc_max = 50.0
    leaf_tempc_n = 100
    leaf_tempc_vec = np.linspace(leaf_tempc_min,leaf_tempc_max,num=leaf_tempc_n)

    vcmax_f = c_double(-9)
    jmax_f  = c_double(-9)
    kp_f    = c_double(-9)
    
    for pft in range(numpft):

        print('Evaluating PFT {}'.format(pft+1))
        jmax25_top,kp25_top =  GetJmaxKp25Top(fates_leaf_vcmax25top[pft])
        vcmax = np.zeros([leaf_tempc_n])
        jmax  = np.zeros([leaf_tempc_n])
        kp    = np.zeros([leaf_tempc_n])

        for it, leaf_tempc in enumerate(leaf_tempc_vec):

            leaf_tempk = leaf_tempc + tfrz_1atm
            
            iret = f90_biophysrate_sub(ci(pft+1), c8(fates_leaf_vcmax25top[pft]), \
                                       c8(jmax25_top), c8(kp25_top), \
                                       c8(nscaler_top), c8(leaf_tempk), c8(dayl_factor_full), \
                                       c8(t_growth_kum),c8(t_home_kum),c8(btran_nolimit), \
                                       byref(vcmax_f), byref(jmax_f), byref(kp_f))
            
            vcmax[it] = vcmax_f.value
            jmax[it]  = jmax_f.value
            kp[it]    = kp_f.value

            
        if(leaf_c3psn[pft] == 0):
            fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8.5,7.5))
        else:
            fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8.5,5.5))
        
        ax1.plot(leaf_tempc_vec,vcmax)
        ax1.set_ylabel('Vcmax [umol/m2/s]')
        ax1.set_title('PFT: {}'.format(pft+1))
        ax1.grid(True)

        ax2.plot(leaf_tempc_vec,jmax)
        ax2.set_xlabel('Leaf Temperature [C]')
        ax2.set_ylabel('Jmax [umol/m2/s]')
        ax2.set_title('PFT: {}'.format(pft+1))
        ax2.grid(True)

        if(leaf_c3psn[pft] == 0):
            ax3.plot(leaf_tempc_vec,kp/umol_per_mol)
            ax3.set_xlabel('Leaf Temperature [C]')
            ax3.set_ylabel('Kp [mol/m2/s]')
            ax3.grid(True)
            ax3.set_xlabel('Leaf Temperature [C]')
            ax4.axis("off")
        else:
            ax1.set_xlabel('Leaf Temperature [C]')


            

# ========================================================================

    
def LinePlotY3dM1(ax,x1,x2,x3,y3d,str_x2,str_x3,add_labels):

    # This takes a 3d array and plots that array over
    # the first of its 3 dimensions.  It does this
    # by evaluating 3 percentiles in each of the other
    # two dimensions, such as the 0, 50% and 99%, and plotting
    # the nine different combinations. The plot axis
    # is passed in, and there is minimal stylizing of the
    # axis and plot formatting, which is left to be done
    # outside of this, and in the calling script
    
    # Find the indices on the second and third dimensions
    # that are the 10th,50th and 90th percentiles
    ix2_10 = int(len(x2)*0.0)
    ix2_50 = int(len(x2)*0.5)
    ix2_90 = int(len(x2)*0.99)
    ix3_10 = int(len(x3)*0.0)
    ix3_50 = int(len(x3)*0.5)
    ix3_90 = int(len(x3)*0.99)


    ix2 = ix2_10;ix3=ix3_10;
    fmt_str = '%s = %4.1f %s = %4.1f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dotted',color='red',label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dotted',color='red')
        
    ix2 = ix2_10;ix3=ix3_50;
    fmt_str = '%s = %4.1f %s = %4.1f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='solid',color='red',label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='solid',color='red')
        
    ix2 = ix2_10;ix3=ix3_90;
    fmt_str = '%s = %4.1f %s = %4.1f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dashed',color='red',label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dashed',color='red')
        
    ix2 = ix2_50;ix3=ix3_10;
    fmt_str = '%s = %4.1f %s = %4.1f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dotted',color='black',label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dotted',color='black')
        
    ix2 = ix2_50;ix3=ix3_50;
    fmt_str = '%s = %4.1f %s = %4.1f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='solid',color='black',label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='solid',color='black')
        
    ix2 = ix2_50;ix3=ix3_90;
    fmt_str = '%s = %4.1f %s = %4.1f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dashed',color='black',label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dashed',color='black')
        
    ix2 = ix2_90;ix3=ix3_10;
    fmt_str = '%s = %4.1f %s = %4.1f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dotted',color='blue',label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dotted',color='blue')
        
    ix2 = ix2_90;ix3=ix3_50;
    fmt_str = '%s = %4.1f %s = %4.1f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='solid',color='blue',label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='solid',color='blue')
        
    ix2 = ix2_90;ix3=ix3_90;
    fmt_str = '%s = %4.1f %s = %4.1f'%(str_x2,x2[ix2],str_x3,x3[ix3])
    if(add_labels):
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dashed',color='blue',label=fmt_str)
    else:
        ap1 = ax.plot(x1,y3d[:,ix2,ix3],linestyle='dashed',color='blue')
        
    ax.grid(True)
    
    
# =======================================================================================


def CompareKitajimaCond(pft,fates_stoich_nitr,fates_leaf_slatop,leaf_stomatal_intercept,fates_maintresp_leaf_model):
    
    # Compare conductances to the paired Amax conductance
    # measurements performed by Kitajima et al. 2005
    # These are for early-late, early and pioneer tropical trees at BCI

    # Photosynthesis (umol m-2 s-1) [early-late, early-late, early, early, pioneer]
    amax_sun_mean_kitajima = [7.3,10.5,10.8,13.7,25.5]
    amax_sun_sd_kitajima   = [2.8, 3.1, 3.0,2.2, 3.6]
    amax_sha_mean_kitajima = [5.3, 8.0, 5.5, 4.5 ]
    amax_sha_sd_kitajima   = [3.6, 2.5, 2.8, 1.5 ]
    
    # Conductance (mol h2o m-2 s-1) [early-late, early-late, early, early, pioneer]
    gs_sun_mean_kitajima = [0.213,0.505,0.412,0.276,0.8]
    gs_sun_sd_kitajima   = [0.098,0.283,0.182,0.083,0.345]
    gs_sha_mean_kitajima = [0.164,0.217,0.264,0.139]
    gs_sha_sd_kitajima   = [0.128,0.069,0.192,0.084]
    
    # Leaf temperature ranges [C]
    leaf_tempc_min = 15.0
    leaf_tempc_max = 40.0
    leaf_tempc_n = 100
    leaf_tempc_vec = np.linspace(leaf_tempc_min,leaf_tempc_max,num=leaf_tempc_n)
    
    # Relative Humidity Ranges
    rh_max = 0.99
    rh_min = 0.70
    rh_n   = 100
    rh_vec = np.linspace(rh_min,rh_max,num=rh_n)

    gs_max_sun_k = np.zeros([leaf_tempc_n,rh_n,3,5])
    gs_max_sha_k = np.zeros([leaf_tempc_n,rh_n,3,4])
    veg_esat_vec = np.zeros([leaf_tempc_n])

    # Generic boundary layer conductance from Kimura et al ~1.2 cm/s
    # reasonable ranges for green houses are 0-3 cm/s
    # perhaps storms would be >10cm/s?
    # Convert to molar form using 1 standard atm at 25C
    # Units:  umol/m2/s
    bl_cond = 100000.0*1.2*0.01*f90_velotomolarcf_sub(c8(can_press_1atm),c8(leaf_tempk25))
    
    
    # Leaf Nitrogen Concentration at the top
    lnc_top  = fates_stoich_nitr/fates_leaf_slatop

    # Range of Amax's between the 
    leaf_tempc_vec = np.linspace(leaf_tempc_min,leaf_tempc_max,num=leaf_tempc_n)
    
    veg_qs_f = c_double(-9.0)
    veg_es_f = c_double(-9.0)
    qsdt_dummy_f = c_double(-9.0)
    esdt_dummy_f = c_double(-9.0)
    lmr_f = c_double(-9.0)
    gstoma_f = c_double(-9.0)
    
    rdark_scaler = 1.0
    nscaler = 1.0
    
    for it, leaf_tempc in enumerate(leaf_tempc_vec):

        leaf_tempk = leaf_tempc + tfrz_1atm
        iret = f90_qsat_sub(c8(leaf_tempk),c8(can_press_1atm), \
                            byref(veg_qs_f),byref(veg_es_f), \
                            byref(qsdt_dummy_f),byref(esdt_dummy_f))
        
        veg_esat_vec[it] = veg_es_f.value
        
        # Leaf Maintenance Respiration (temp and pft dependent)
        if(fates_maintresp_leaf_model==1):
            iret = f90_lmr_ryan_sub(c8(lnc_top),c8(nscaler), ci(pft+1), c8(leaf_tempk), byref(lmr_f))
        elif(fates_maintresp_leaf_model==2):
            iret = f90_lmr_atkin_sub(c8(lnc_top),c8(rdark_scaler),c8(leaf_tempk),c8(atkin_mean_leaf_tempk),byref(lmr_f) )
        else:
            print('unknown leaf respiration model')
            exit(1)

        anet_k = np.zeros([3])
        
        # Calculate conductance consistent with Amax in Kitajima
        for ir, rh in enumerate(rh_vec):
            vpress = rh * veg_esat_vec[it]

            for isp, amax in enumerate(amax_sun_mean_kitajima):

                anet_k[0] = (amax-amax_sun_sd_kitajima[isp]) - lmr_f.value
                anet_k[1] = amax - lmr_f.value
                anet_k[2] = (amax+amax_sun_sd_kitajima[isp]) - lmr_f.value

                for i in range(len(anet_k)):
                    # Photosynthetic rate at light saturation (Amax)
                    iret = f90_gs_medlyn(c8(anet_k[i]),ci(pft+1), c8(veg_esat_vec[it]), c8(vpress), \
                                         c8(btran_nolimit*leaf_stomatal_intercept), c8(co2_ppress_400ppm),c8(can_press_1atm), \
                                         c8(bl_cond),byref(gstoma_f))
                
                    gs_max_sun_k[it,ir,i,isp] = gstoma_f.value

                
    fig4,ax1 = plt.subplots(1,1,figsize=(6.5,5.5))

    ap = ax1.hist(gs_max_sun_k.reshape(-1)/umol_per_mol, bins = 30)
    
    ax1.set_ylabel('count')
    ax1.set_xlabel('gs [mol H20 m-2 s-1]')
    ap = ax1.axvline(x=gs_sun_mean_kitajima[0], linewidth=2, color=[0.3,0.3,0.3])
    ap = ax1.axvline(x=gs_sun_mean_kitajima[1], linewidth=2, color=[0.3,0.3,0.3])
    ap = ax1.axvline(x=gs_sun_mean_kitajima[2], linewidth=2, color=[0.3,0.3,0.3])
    ap = ax1.axvline(x=gs_sun_mean_kitajima[3], linewidth=2, color=[0.3,0.3,0.3])
    ap = ax1.axvline(x=gs_sun_mean_kitajima[4], linewidth=2, color=[0.3,0.3,0.3])
    ax1.set_title('Conductance Comparison with Kitajima et. al. 2005')


# =======================================================================================

# Instantiate the F90 modules

f90_const_obj = ctypes.CDLL('bld/FatesConstantsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_shr_obj = ctypes.CDLL('bld/WrapShrMod.o',mode=ctypes.RTLD_GLOBAL)
f90_fatesutils_obj = ctypes.CDLL('bld/FatesUtilsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_leaf_biophys_obj = ctypes.CDLL('bld/LeafBiophysicsMod.o',mode=ctypes.RTLD_GLOBAL)
f90_leaf_biophys_supp_obj = ctypes.CDLL('bld/LeafBiophysSuppMod.o',mode=ctypes.RTLD_GLOBAL)

# Identify subroutine objects, so we can call them
f90_set_leaf_param_sub = getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','setleafparam'))
f90_alloc_leaf_param_sub = getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','allocleafparam'))
f90_dealloc_leaf_param_sub = getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','deallocleafparam'))
f90_dump_param_sub =  getattr(f90_leaf_biophys_supp_obj, GetModSymbol('bld/LeafBiophysSuppMod.o','DumpParams'))
f90_set_leaf_param_sub.argtypes = [POINTER(c_double),POINTER(c_int),c_char_p,c_long]
f90_biophysrate_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerBiophysicalRate'))
f90_leaflayerphoto_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerPhotosynthesis'))
f90_qsat_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','QSat'))
f90_cangas_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','GetCanopyGasParameters'))
f90_lmr_ryan_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerMaintenanceRespiration_Ryan_1991'))
f90_lmr_atkin_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerMaintenanceRespiration_Atkin_etal_2017'))
f90_agross_rubiscoc3  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossRubiscoC3'))
f90_agross_rubpc3  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossRuBPC3'))
f90_agross_rubpc4  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossRuBPC4'))
f90_agross_pepc4  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossPEPC4'))

f90_gs_medlyn = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','StomatalCondMedlyn'))
f90_gs_ballberry = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','StomatalCondBallBerry'))
f90_velotomolarcf_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','VeloToMolarCF'))

# For functions, define the return value
f90_agross_rubiscoc3.restype = c_double
f90_agross_rubpc3.restype = c_double
f90_agross_rubpc4.restype = c_double
f90_agross_pepc4.restype = c_double
f90_velotomolarcf_sub.restype = c_double
    
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
    if(dump_parameters):
        iret = f90_dump_param_sub()

    
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

    if(do_evalvjkbytemp):
        EvalVJKByTemp(numpft,fates_leaf_vcmax25top,leaf_c3psn)
            
    # Drive conductance with Amax values from Kitajima et al for BTEs
    # Compare with stomatal conductance from Kitajima et al.
    if(do_comparekitajimacond):
        CompareKitajimaCond(0,fates_stoich_nitr[0],fates_leaf_slatop[0], \
                            leaf_stomatal_intercept[0],fates_maintresp_leaf_model)


    #if(do_paramsense):
    #    ParamSense()
        
    # Leaf temperature ranges [C]
    leaf_tempc_min = -20.0
    leaf_tempc_max = 40.0
    leaf_tempc_n = 100
    leaf_tempc_vec = np.linspace(leaf_tempc_min,leaf_tempc_max,num=leaf_tempc_n)
    
    # These variables are mostly dependent on leaf temperature,
    # and weakly dependent on atmospheric pressure
    veg_qsat_vec = np.zeros([leaf_tempc_n])
    veg_esat_vec = np.zeros([leaf_tempc_n])
    mm_kco2_vec = np.zeros([leaf_tempc_n])
    mm_ko2_vec = np.zeros([leaf_tempc_n])
    co2_cpoint_vec = np.zeros([leaf_tempc_n])
    
    # Absorbed PAR ranges [W/m2]
    par_abs_min = 0.01
    par_abs_max = 100
    par_abs_n  = 100
    par_abs_vec = np.linspace(par_abs_min,par_abs_max,num=par_abs_n)

    # Relative Humidity Ranges
    rh_max = 0.99
    rh_min = 0.20
    rh_n   = 100
    rh_vec = np.linspace(rh_min,rh_max,num=rh_n)
    
    # Generic boundary layer conductance from Kimura et al ~1.2 cm/s
    # reasonable ranges for green houses are 0-3 cm/s
    # perhaps storms would be >10cm/s?
    # Convert to molar form using 1 standard atm at 25C
    # Units:  umol/m2/s
    bl_cond = 100000.*1.2*0.01*f90_velotomolarcf_sub(c8(can_press_1atm),c8(leaf_tempk25))

    # Lets look at canopy top
    # kn = DecayCoeffVcmax(currentCohort%vcmax25top, &
    #                     EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff1(ft), &
    #                     EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff2(ft))
    # rdark_scaler = exp(-kn * cumulative_lai)

    rdark_scaler = 1.0

    # kn = DecayCoeffVcmax(currentCohort%vcmax25top, &
    #                      prt_params%leafn_vert_scaler_coeff1(ft), &
    #                      prt_params%leafn_vert_scaler_coeff2(ft))

    nscaler = 1.0

    # Initialize fortran return values
    # these are all temps and doubles
    vcmax_f      = c_double(-9.0)
    jmax_f       = c_double(-9.0)
    kp_f         = c_double(-9.0)
    agross_f     = c_double(-9.0)
    gstoma_f     = c_double(-9.0)
    anet_f       = c_double(-9.0)
    lmr_f        = c_double(-9.0)
    c13_f        = c_double(-9.0)
    ac_f         = c_double(-9.0)
    aj_f         = c_double(-9.0)
    ap_f         = c_double(-9.0)
    co2_interc_f = c_double(-9.0)
    veg_qs_f     = c_double(-9.0)
    veg_es_f     = c_double(-9.0)
    mm_kco2_f    = c_double(-9.0)
    mm_ko2_f     = c_double(-9.0)
    co2_cpoint_f = c_double(-9.0)
    qsdt_dummy_f = c_double(-9.0)
    esdt_dummy_f = c_double(-9.0)

    print('Prepping Canopy Gas Parameters')
    
    for it, leaf_tempc in enumerate(leaf_tempc_vec):

        leaf_tempk = leaf_tempc + tfrz_1atm
            
        iret = f90_qsat_sub(c8(leaf_tempk),c8(can_press_1atm), \
                            byref(veg_qs_f),byref(veg_es_f), \
                            byref(qsdt_dummy_f),byref(esdt_dummy_f))

        
        veg_qsat_vec[it] = veg_qs_f.value
        veg_esat_vec[it] = veg_es_f.value
        
        iret = f90_cangas_sub(c8(can_press_1atm), \
                              c8(o2_ppress_209kppm), \
                              c8(leaf_tempk), \
                              byref(mm_kco2_f), \
                              byref(mm_ko2_f), \
                              byref(co2_cpoint_f))

        mm_kco2_vec[it] = mm_kco2_f.value
        mm_ko2_vec[it] = mm_ko2_f.value
        co2_cpoint_vec[it] = co2_cpoint_f.value

    
        

    print('\n')
    print('Experiment 1: Evaluating Photosynthesis Equations by pft/Tl/RH/PR')
    
    for pft in [0,11]: #range(numpft):

        print('Evaluating PFT {}'.format(pft+1))
        
        jmax25_top,kp25_top =  GetJmaxKp25Top(fates_leaf_vcmax25top[pft])
        vcmax = np.zeros([leaf_tempc_n])
        jmax  = np.zeros([leaf_tempc_n])
        kp    = np.zeros([leaf_tempc_n])
        lmr    = np.zeros([leaf_tempc_n])
        agross = np.zeros([leaf_tempc_n,rh_n,par_abs_n])
        gstoma = np.zeros([leaf_tempc_n,rh_n,par_abs_n])
        anet   = np.zeros([leaf_tempc_n,rh_n,par_abs_n])
        ac     = np.zeros([leaf_tempc_n,rh_n,par_abs_n])
        aj     = np.zeros([leaf_tempc_n,rh_n,par_abs_n])
        ap     = np.zeros([leaf_tempc_n,rh_n,par_abs_n])
        co2_interc = np.zeros([leaf_tempc_n,rh_n,par_abs_n])

        gs_max_sun_k = np.zeros([leaf_tempc_n,rh_n])
        
        # When calling component limitations exclusively
        # using an approximation of interstitial co2 as
        # 0.7*canopy_co2
        
        ac2    = np.zeros([leaf_tempc_n])
        aj2    = np.zeros([leaf_tempc_n,par_abs_n])
        ap2    = np.zeros([leaf_tempc_n])
        
        # Leaf Nitrogen Concentration at the top
        lnc_top  = fates_stoich_nitr[pft]/fates_leaf_slatop[pft]
        
        for it, leaf_tempc in enumerate(leaf_tempc_vec):

            leaf_tempk = leaf_tempc + tfrz_1atm
            
            iret = f90_biophysrate_sub(ci(pft+1), c8(fates_leaf_vcmax25top[pft]), \
                                       c8(jmax25_top), c8(kp25_top), \
                                       c8(nscaler), c8(leaf_tempk), c8(dayl_factor_full), \
                                       c8(t_growth_kum),c8(t_home_kum),c8(btran_nolimit), \
                                       byref(vcmax_f), byref(jmax_f), byref(kp_f))
            
            vcmax[it] = vcmax_f.value
            jmax[it]  = jmax_f.value
            kp[it]    = kp_f.value

            if(leaf_c3psn[pft] == 1):
                ap2[it] = 0.0
            else:
                ap2[it] = f90_agross_pepc4(c8(co2_ppress_400ppm),c8(kp[it]),c8(can_press_1atm))
                
            # Leaf Maintenance Respiration (temp and pft dependent)
            if(fates_maintresp_leaf_model==1):
                iret = f90_lmr_ryan_sub(c8(lnc_top),c8(nscaler), ci(pft+1), c8(leaf_tempk), byref(lmr_f))
            elif(fates_maintresp_leaf_model==2):
                iret = f90_lmr_atkin_sub(c8(lnc_top),c8(rdark_scaler),c8(leaf_tempk),c8(atkin_mean_leaf_tempk),byref(lmr_f) )
            else:
                print('unknown leaf respiration model')
                exit(1)

            lmr[it] = lmr_f.value

            if(leaf_c3psn[pft] == 1):
                
                ac2[it] = f90_agross_rubiscoc3(c8(vcmax[it]),c8(ci_ppress_static),c8(o2_ppress_209kppm), \
                                             c8(co2_cpoint_vec[it]),c8(mm_kco2_vec[it]),c8(mm_ko2_vec[it]))
            else:
                ac2[it] = vcmax[it]

                 
            for ip, par_abs in enumerate(par_abs_vec):

                if(leaf_c3psn[pft] == 1):
                    aj2[it,ip] = f90_agross_rubpc3(c8(par_abs),c8(jmax[it]),c8(ci_ppress_static),c8(co2_cpoint_vec[it]))
                else:
                    aj2[it,ip] = f90_agross_rubpc4(c8(par_abs))
                    
                for ir, rh in enumerate(rh_vec):

                    vpress = rh * veg_esat_vec[it]
            
                    iret = f90_leaflayerphoto_sub(c8(par_abs),  \
                                                  c8(1.0),     \
                                                  ci(pft+1),   \
                                                  c8(vcmax[it]),   \
                                                  c8(jmax[it]),    \
                                                  c8(kp[it]),      \
                                                  c8(leaf_tempk), \
                                                  c8(can_press_1atm), \
                                                  c8(co2_ppress_400ppm), \
                                                  c8(o2_ppress_209kppm), \
                                                  c8(btran_nolimit), \
                                                  c8(bl_cond), \
                                                  c8(vpress), \
                                                  c8(mm_kco2_vec[it]), \
                                                  c8(mm_ko2_vec[it]), \
                                                  c8(co2_cpoint_vec[it]), \
                                                  c8(lmr[it]), \
                                                  c8(zero_lwp), \
                                                  byref(agross_f), \
                                                  byref(gstoma_f), \
                                                  byref(anet_f), \
                                                  byref(c13_f), \
                                                  byref(ac_f), \
                                                  byref(aj_f), \
                                                  byref(ap_f), \
                                                  byref(co2_interc_f))

                    agross[it,ir,ip] = agross_f.value
                    gstoma[it,ir,ip] = gstoma_f.value
                    anet[it,ir,ip] = anet_f.value
                    ac[it,ir,ip] = ac_f.value
                    aj[it,ir,ip] = aj_f.value
                    ap[it,ir,ip] = ap_f.value
                    co2_interc[it,ir,ip] = co2_interc_f.value

                    
        


        # Plot out component gross assimilation rates
        # by temperature with constant Ci, and by Ci with
        # constant temperature
            
        fig2,ax1 = plt.subplots(1,1,figsize=(6.5,5.5))

        ax1.plot(leaf_tempc_vec,ac2,label='Ac')
        ix2_10 = int(par_abs_n*0.1)
        ix2_50 = int(par_abs_n*0.5)
        ix2_90 = int(par_abs_n*0.9)
        ax1.plot(leaf_tempc_vec,aj2[:,ix2_10],color=[0.5,0.5,0.5],linestyle='dotted',label='Aj apar=%4.1f'%(par_abs_vec[ix2_10]))
        ax1.plot(leaf_tempc_vec,aj2[:,ix2_50],color=[0.5,0.5,0.5],linestyle='dashed',label='Aj apar=%4.1f'%(par_abs_vec[ix2_50]))
        ax1.plot(leaf_tempc_vec,aj2[:,ix2_90],color=[0.5,0.5,0.5],linestyle='solid',label='Aj apar=%4.1f'%(par_abs_vec[ix2_90]))
        #if(leaf_c3psn[pft] == 0):
        #    ax1.axhline(y=ap2,linewidth=2, color=[0.3,0.3,0.3])
        ax1.plot(leaf_tempc_vec,ap2[:],color='orange',linestyle='solid',label='Ap')
        
        ax1.set_ylabel('[umol/m2/s]')
        ax1.set_xlabel('Leaf Temperature [C]')
        ax1.set_title('PFT: %3i, Vcmax25: %4.1f, Jmax25: %4.1f, Ci: %4.1f'%(pft+1,fates_leaf_vcmax25top[pft],jmax25_top,ci_ppress_static))
        ax1.grid(True)
        fig2.legend(loc='upper right')

        # Lets plot metrics by temperature, using the
        # 10th, 50th and 90th percentiles of both RH and PAR
        # Agross, Anet, Gstoma, Ac, Aj, Ap
        
        fig3, ((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8)) = plt.subplots(4,2,figsize=(7.,9.))

        LinePlotY3dM1(ax1,leaf_tempc_vec,rh_vec,par_abs_vec,agross,'RH','APAR',True)
        ax1.set_ylabel('Agross \n [umol/m2/s]')
        ax1.set_xticklabels([])
        
        LinePlotY3dM1(ax2,leaf_tempc_vec,rh_vec,par_abs_vec,gstoma*1.e-6,'RH','APAR',False)
        ax2.set_ylabel('Gs \n [mol/m2/s]')
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()
        ax2.axhline(y=bl_cond*1.e-6,linewidth=2, color=[0.3,0.3,0.3])
        ax2.set_xticklabels([])
        
        LinePlotY3dM1(ax3,leaf_tempc_vec,rh_vec,par_abs_vec,ac,'RH','APAR',False)
        ax3.set_ylabel('Ac (Rubisco) \n [umol/m2/s]')
        ax3.set_xticklabels([])
        
        LinePlotY3dM1(ax4,leaf_tempc_vec,rh_vec,par_abs_vec,aj,'RH','APAR',False)
        ax4.set_ylabel('Aj (RuBP) \n [umol/m2/s]')
        ax4.yaxis.set_label_position("right")
        ax4.yaxis.tick_right()
        ax4.set_xticklabels([])
        
        LinePlotY3dM1(ax5,leaf_tempc_vec,rh_vec,par_abs_vec,co2_interc,'RH','APAR',False)
        ax5.set_ylabel('Ci \n [Pa]')
        ax5.axhline(y=co2_ppress_400ppm,linewidth=2, color=[0.3,0.3,0.3])
        
        ax6.plot(leaf_tempc_vec, lmr)
        ax6.set_ylabel('LMR \n [umol/m2/s]')
        ax6.grid(True)
        ax6.yaxis.set_label_position("right")
        ax6.yaxis.tick_right()
        
        ax5.set_xlabel('Leaf Temperature [C]')
        ax6.set_xlabel('Leaf Temperature [C]')

        ax7.axis("off")
        ax8.axis("off")
        
        
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.02, hspace=0.03)
        fig3.legend(loc='lower center',labelspacing = 0.2)
        
        
        plt.show()

    
    print('Deallocating parameter space')
    iret = f90_dealloc_leaf_param_sub()
    
    print('Functional Unit Testing Complete')
    exit(0)

    
# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
