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

# ========================================================================

def PlotVJKByTemp(leaf_tempc_vec,vcmax,jmax,kp,pft,leaf_c3psn):

                    
    fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8.5,7.5))

    ap = ax1.plot(leaf_tempc_vec,vcmax)
    
    ax1.set_ylabel('Vcmax [umol/m2/s]')
    ax1.set_title('PFT: {}'.format(pft+1))
    ax1.grid(True)

    ap = ax2.plot(leaf_tempc_vec,jmax)
    ax2.set_xlabel('Leaf Temperature [C]')
    ax2.set_ylabel('jmax [umol/m2/s]')
    ax2.set_title('PFT: {}'.format(pft+1))
    ax2.grid(True)

    if(leaf_c3psn == 0):
        ap = ax3.semilogy(leaf_tempc_vec,kp)
        ax3.set_xlabel('Leaf Temperature [C]')
        ax3.set_ylabel('Kp (C4 only)')
        ax3.grid(True)
        ax3.set_xlabel('Leaf Temperature [C]')
    else:
        ax1.set_xlabel('Leaf Temperature [C]')
        ax3.axis("off")

    ax4.axis("off")

    
def LinePlotY3dM1(ax,x1,x2,x3,y3d,str_x2,str_x3,name_str,add_labels,x_line):

    # Find the indices on the second and third dimensions
    # that are the 10th,50th and 90th percentiles
    ix2_10 = int(len(x2)*0.1)
    ix2_50 = int(len(x2)*0.5)
    ix2_90 = int(len(x2)*0.9)
    ix3_10 = int(len(x3)*0.1)
    ix3_50 = int(len(x3)*0.5)
    ix3_90 = int(len(x3)*0.9)


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
        
    ax.set_ylabel(name_str)

    if x_line is not None:
        ax.hlines(y=x_line, xmin=np.min(x1), xmax=np.max(x1),linewidth=2, color=[0.3,0.3,0.3])


    
    
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
#f90_biophysrate_sub.argtypes = [POINTER(c_int), \
#                                POINTER(c_double),POINTER(c_double),POINTER(c_double), \
#                                POINTER(c_double),POINTER(c_double),POINTER(c_double), \
#                                POINTER(c_double),POINTER(c_double),POINTER(c_double), \
#                                POINTER(c_double),POINTER(c_double),POINTER(c_double)]

f90_leaflayerphoto_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerPhotosynthesis'))
f90_qsat_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','QSat'))
f90_cangas_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','GetCanopyGasParameters'))
f90_lmr_ryan_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerMaintenanceRespiration_Ryan_1991'))
f90_lmr_atkin_sub = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','LeafLayerMaintenanceRespiration_Atkin_etal_2017'))
f90_agross_rubiscoc3  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossRubiscoC3'))
f90_agross_rubpc3  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossRuBPC3'))
f90_agross_rubpc4  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossRuBPC4'))
f90_agross_pepc4  = getattr(f90_leaf_biophys_obj,GetModSymbol('bld/LeafBiophysicsMod.o','AgrossPEPC4'))
f90_agross_rubiscoc3.restype = c_double
f90_agross_rubpc3.restype = c_double
f90_agross_rubpc4.restype = c_double
f90_agross_pepc4.restype = c_double


# Example of setting a function return call
#dftcdpsi_from_psi.restype = c_double

    
# For debugging
dump_parameters = False

# Plot Control
do_plot_vjk = True


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

                


    # Call the BiophysicalRates() Routine, which calculates
    # the actual vmax, jmax and rcurve-islope terms for the
    # leaf-scale conditions
    tfrz_1atm = 273.15
    leaf_tempk25 = tfrz_1atm + 25.0  # Use 25 to start
    dayl_factor = 1.0
    t_growth_kum = -999
    t_home_kum = -999
    btran = 1.0

    # Leaf temperature ranges [C]
    leaf_tempc_min = 00.0
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
    rh_min = 0.01
    rh_n   = 100
    rh_vec = np.linspace(rh_min,rh_max,num=rh_n)

    can_press_1atm = 101325.0
    co2_ppress_400ppm = 0.0004*can_press_1atm
    o2_ppress_209kppm = 0.29*can_press_1atm
    btran_nolim = 1.0

    ci_ppress_static = 0.7*co2_ppress_400ppm
    
    # Generic boundary layer conductance from Kimura et al ~1.2 cm/s
    # reasonable ranges for green houses are 0-3 cm/s
    # perhaps storms would be >10cm/s?
    # Convert to molar form using 1 standard atm at 25C
    bl_cond = MolarToVeloCF(can_press_1atm,298.0)*1.2*0.01

    # Respiration does affect conductance, which also affects
    # photosynthesis, but for now we use 0
    zero_lmr = 0.0

    # Set Leaf water potential to zero for some calcs
    zero_lwp = 0.0


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
    vcmax_f  = c_double(-9.0)
    jmax_f   = c_double(-9.0)
    kp_f     = c_double(-9.0)
    agross_f = c_double(-9.0)
    gstoma_f = c_double(-9.0)
    anet_f   = c_double(-9.0)
    lmr_f    = c_double(-9.0)
    c13_f    = c_double(-9.0)
    ac_f     = c_double(-9.0)
    aj_f     = c_double(-9.0)
    ap_f     = c_double(-9.0)
    co2_interc_f = c_double(-9.0)
    veg_qs_f = c_double(-9.0)
    veg_es_f = c_double(-9.0)
    mm_kco2_f = c_double(-9.0)
    mm_ko2_f  = c_double(-9.0)
    co2_cpoint_f = c_double(-9.0)
    qsdt_dummy_f = c_double(-9.0)
    esdt_dummy_f = c_double(-9.0)

    aj2_f = c_double(-9.)


    

    print('Prepping Canopy Gas Parameters')
    
    for it, leaf_tempc in enumerate(leaf_tempc_vec):

        leaf_tempk = leaf_tempc + tfrz_1atm
            
        iret = f90_qsat_sub(c8(leaf_tempk),c8(can_press_1atm), \
                            byref(veg_qs_f),byref(veg_es_f), \
                            byref(qsdt_dummy_f),byref(esdt_dummy_f))


        
        
        veg_qsat_vec[it] = veg_qs_f.value

        print(veg_qsat_vec[it], (veg_qs_f.value * can_press_1atm) / (0.622 + veg_qs_f.value), veg_es_f.value)
        
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
    
    for pft in range(numpft):

        print('Evaluating PFT {}'.format(pft))
        
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
                                       c8(nscaler), c8(leaf_tempk), c8(dayl_factor), \
                                       c8(t_growth_kum),c8(t_home_kum),c8(btran), \
                                       byref(vcmax_f), byref(jmax_f), byref(kp_f))
            
            vcmax[it] = vcmax_f.value
            jmax[it]  = jmax_f.value
            kp[it]    = kp_f.value

            if(leaf_c3psn[pft] == 1):
                ap2[it] = 0.0
            else:
                ap2_f = f90_agross_pepc4(c8(co2_ppress_400ppm),c8(kp[it]),c8(can_press_1atm))
                ap2[it] = ap2_f.value
                
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
                
                ac2_f = f90_agross_rubiscoc3(c8(vcmax[it]),c8(ci_ppress_static),c8(o2_ppress_209kppm), \
                                             c8(co2_cpoint_vec[it]),c8(mm_kco2_vec[it]),c8(mm_ko2_vec[it]))
                ac2[it] = ac2_f#.value
            else:
                print('Add c4 for Ac')
                exit(0)
                
            for ip, par_abs in enumerate(par_abs_vec):

                if(leaf_c3psn[pft] == 1):
                    aj2_f = f90_agross_rubpc3(c8(par_abs),c8(jmax[it]),c8(ci_ppress_static),c8(co2_cpoint_vec[it]))
                    aj2[it,ip] = aj2_f#.value
                else:
                    print('Add c4 for Aj')
                    exit(0)
                    
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
                                                  c8(btran_nolim), \
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

                    #print("par %4.1f Tl  %4.1f P     %7.1f Pco2 %4.1f"%(par_abs,leaf_tempk,can_press_1atm,co2_ppress_400ppm))
                    #print("          Po2 %6.1f btran %2.1f Gb   %5.1f"%(o2_ppress_209kppm,btran_nolim,bl_cond))
                    #print("          Ec  %5.1f Kco2  %5.1f Ko2  %5.1f"%(vpress,mm_kco2_vec[it],mm_ko2_vec[it]))
                    #print("          Pcc %4.1f Rl    %6.2f"%(co2_cpoint_vec[it],lmr[it]))

                    
        if(do_plot_vjk):
            PlotVJKByTemp(leaf_tempc_vec,vcmax,jmax,kp,pft,leaf_c3psn)

        
        fig2,ax1 = plt.subplots(1,1,figsize=(6.5,5.5))

        ap = ax1.plot(leaf_tempc_vec,ac2,label='Ac')
        ix2_10 = int(par_abs_n*0.1)
        ix2_50 = int(par_abs_n*0.5)
        ix2_90 = int(par_abs_n*0.9)
        ap = ax1.plot(leaf_tempc_vec,aj2[:,ix2_10],color=[0.5,0.5,0.5],linestyle='dotted',label='Aj apar=%4.1f'%(par_abs_vec[ix2_10]))
        ap = ax1.plot(leaf_tempc_vec,aj2[:,ix2_50],color=[0.5,0.5,0.5],linestyle='dashed',label='Aj apar=%4.1f'%(par_abs_vec[ix2_50]))
        ap = ax1.plot(leaf_tempc_vec,aj2[:,ix2_90],color=[0.5,0.5,0.5],linestyle='solid',label='Aj apar=%4.1f'%(par_abs_vec[ix2_90]))
        if(leaf_c3psn[pft] == 0):
            ap = ax.hlines(y=ap2, xmin=np.min(leaf_tempc_vec), xmax=np.max(leaf_tempc_vec),linewidth=2, color=[0.3,0.3,0.3])
            
        ax1.set_ylabel('[umol/m2/s]')
        ax1.set_xlabel('Leaf Temperature [C]')
        ax1.set_title('PFT: %3i, Vcmax25: %4.1f, Jmax25: %4.1f, Ci: %4.1f'%(pft+1,fates_leaf_vcmax25top[pft],jmax25_top,ci_ppress_static))
        ax1.grid(True)
        fig2.legend(loc='upper right')
        
        # Plot out inputs

        # Lets plot metrics by temperature, using the
        # 10th, 50th and 90th percentiles of both RH and PAR
        # Agross, Anet, Gstoma, Ac, Aj, Ap
        
        fig3, ((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8)) = plt.subplots(4,2,figsize=(8.5,6.5))

        LinePlotY3dM1(ax1,leaf_tempc_vec,rh_vec,par_abs_vec,agross,'RH','APAR','Agross [umol/m2/s]',True,None)
        LinePlotY3dM1(ax2,leaf_tempc_vec,rh_vec,par_abs_vec,gstoma,'RH','APAR','Gs [umol/m2/s]',False,bl_cond)
        LinePlotY3dM1(ax3,leaf_tempc_vec,rh_vec,par_abs_vec,ac,'RH','APAR','Ac (Rubisco) [umol/m2/s]',False,None)
        LinePlotY3dM1(ax4,leaf_tempc_vec,rh_vec,par_abs_vec,aj,'RH','APAR','Aj (RuBP) [umol/m2/s]',False,None)
        LinePlotY3dM1(ax5,leaf_tempc_vec,rh_vec,par_abs_vec,co2_interc,'RH','APAR','CO2 Interc [Pa]',False,None)
        

        ax6.plot(leaf_tempc_vec, lmr)
        ax6.set_ylabel('LMR [umol/m2/s]')
        ax6.grid(True)

        
        #ax6.plot(leaf_tempc_vec,mm_kco2_vec)
        
        ax8.plot(leaf_tempc_vec,mm_ko2_vec)
        ax8.set_ylabel('MM Ko2 [Pa]')
        ax8.grid(True)


        
        ax7.plot(leaf_tempc_vec,co2_cpoint_vec)
        ax7.set_ylabel('Co2 cpoint [Pa]')
        ax7.grid(True)
        
        #ax6.set_ylabel('Qsat [g/kg]')
        
        #@if(
        #LinePlotY3dM1(ax2,leaf_tempc_vec,rh_vec,par_abs_vec,ac,'RH','APAR','Aj (RuBP) [umol/m2/s]',False,None)
        
       
        

        # PLOT IDEA
        # CONTOURS OF THE AC/AJ EQUIVALENCY FOR TEMP AND LIGHT FOR DIFFERENT
        # VCMAXES, DIFFERENT PLOTS FOR DIFFERENT RH AND RBlayer

        #ax6.axis("off")

        fig3.legend(loc='lower right')
        plt.tight_layout()
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
