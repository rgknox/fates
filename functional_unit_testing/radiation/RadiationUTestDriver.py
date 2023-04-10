# =======================================================================================
#
# For usage: $python RadiationUTestDriver.py --help
#
# This script runs unit tests on the two-stream functions.
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
import numpy as np
import matplotlib
import os
import sys
import getopt
#import code  # For development: code.interact(local=dict(globals(), **locals()))
import code  # For development: code.interact(local=locals())  code.interact(local=dict(globals(), **locals()))
import time
import importlib
import ctypes
from ctypes import *
from operator import add
sys.path.append('../shared/py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 12}

matplotlib.rc('font', **font)


# Instantiate the F90 modules
f90_twostr_obj = ctypes.CDLL('bld/TwoStreamPPAMod.o',mode=ctypes.RTLD_GLOBAL)
f90_wrap_obj = ctypes.CDLL('bld/RadiationWrapMod.o',mode=ctypes.RTLD_GLOBAL)

# Create aliases for the calls and define arguments if it helps with clarity
alloc_twostream_call =  f90_wrap_obj.__radiationwrapmod_MOD_initallocate
dealloc_twostream_call = f90_wrap_obj.__radiationwrapmod_MOD_dealloc
alloc_radparams_call = f90_twostr_obj.__twostreamppamod_MOD_allocateradparams
set_radparams_call   = f90_wrap_obj.__radiationwrapmod_MOD_setradparam
set_radparams_call.argtypes = [POINTER(c_double),POINTER(c_int),POINTER(c_int),c_char_p,c_long]
param_prep_call = f90_twostr_obj.__twostreamppamod_MOD_paramprep

setup_canopy_call = f90_wrap_obj.__radiationwrapmod_MOD_setupcanopy
setup_canopy_call.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_double),POINTER(c_double),POINTER(c_double)]

grndsnow_albedo_call = f90_wrap_obj.__radiationwrapmod_MOD_setgroundsnow
grndsnow_albedo_call.argtypes = [POINTER(c_double),c_char_p,c_long]

canopy_prep_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapcanopyprep
zenith_prep_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapzenithprep
solver_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapsolve

getintens_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapgetintensity
getparams_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapgetparams
forceparam_call = f90_wrap_obj.__radiationwrapmod_MOD_wrapforceparams
forceparam_call.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_double),c_char_p,c_long]

leaf_rhonir = [0.46, 0.41, 0.39, 0.46, 0.41, 0.41, 0.46, 0.41, 0.41, 0.28, 0.28, 0.28 ]
leaf_rhovis = [0.11, 0.09, 0.08, 0.11, 0.08, 0.08, 0.11, 0.08, 0.08, 0.05, 0.05, 0.05 ]
leaf_taunir = [0.33, 0.32, 0.42, 0.33, 0.43, 0.43, 0.33, 0.43, 0.43, 0.4,  0.4,  0.4 ]
leaf_tauvis = [0.06, 0.04, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05]
leaf_xl = [0.32, 0.01, 0.01, 0.32, 0.2, 0.59, 0.32, 0.59, 0.59, -0.23, -0.23, -0.23]
leaf_clumping_index = [0.85, 0.85, 0.8, 0.85, 0.85, 0.9, 0.85, 0.9, 0.9, 0.75, 0.75, 0.75]
stem_rhonir = [0.49, 0.36, 0.36, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.53, 0.53, 0.53]
stem_rhovis = [0.21, 0.12, 0.12, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.31, 0.31, 0.31]
stem_taunir = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.25, 0.25, 0.25]
stem_tauvis = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.12, 0.12, 0.12]

visb = 1
nirb = 2




class elem_type:
    def __init__(self,n_ll):

        self.area = -9.0
        self.lai  = -9.0
        self.sai  = -9.0
        
        self.avai = np.zeros([n_ll])
        self.r_dn = np.zeros([n_ll])
        self.r_up = np.zeros([n_ll])
        self.r_b  = np.zeros([n_ll])
        self.r_abs = np.zeros([n_ll-1])

        
def main(argv):

    # All tests will use 2 bands 1=vis, 2=nir

    # Initialize radiation parameters
    n_bands = 2
    n_pft   = 12

    iret = alloc_radparams_call(ci(n_pft),ci(n_bands))

    for ft in range(n_pft):

        pft=ft+1
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
        
        # Process the core 2Stream parameters from parameters in file
        iret = param_prep_call(ci(n_pft),ci(visb))
        iret = param_prep_call(ci(n_pft),ci(nirb))

    
    # Test 1, a single element, visible

    if(False):
        SingleElementPerturbTest()

    if(True):
        SerialParallelCanopyTest()


    
def SerialParallelCanopyTest():


    # Lets first construct a bunch of cohorts, 5 cohorts
    # equal area, but folding by 2 in LAI

    cohort_lai  = np.array([0.25,0.5,1.0,2.0,4.0])
    
    cohort_area = np.array([0.2,0.2,0.2,0.2,0.2])
    sai_frac = 0.1
    
    pft = 1
    
    # Serial approach: 5 layers with veg and ghost
    n_col = 2
    n_layer = 5
    iret = alloc_twostream_call(ci(n_layer),ci(n_col))

    elems = []
    elems.append([])
    elems.append([])
    n_vai = 100
    
    
    for i in range(n_layer):
        ican = i+1

        elems[0].append(elem_type(n_vai))
        elems[1].append(elem_type(n_vai))
        
        icol = 1
        area = np.sum(cohort_area[i:])
        if(i==0):
            lai = cohort_lai[i]
        else:
            lai = cohort_lai[i]-cohort_lai[i-1]
        
        sai  = lai*sai_frac
        elems[0][-1].lai  = lai
        elems[0][-1].sai  = sai
        elems[0][-1].area = area
        elems[0][-1].avai = np.linspace(0,lai+sai,num=n_vai)
        iret = setup_canopy_call(c_int(ican),c_int(icol),c_int(pft),c_double(area),c_double(lai),c_double(sai))

        icol = 2
        area = np.sum(cohort_area[i:])
        elems[0][-1].lai  = 0.0
        elems[0][-1].sai  = 0.0
        elems[0][-1].area = area
        lai  = 0.0
        sai  = 0.0
        air_pft = 0
        iret = setup_canopy_call(c_int(ican),c_int(icol),c_int(air_pft),c_double(area),c_double(lai),c_double(sai))
        
    # Decide on a band:
    ib = visb
    
    cd_r_beam = c_double(-9.0)
    cd_r_diff_up = c_double(-9.0)
    cd_r_diff_dn = c_double(-9.0)
    cd_kb = c_double(-9.0)
    cd_kd = c_double(-9.0)
    cd_om = c_double(-9.0)
    cd_betad = c_double(-9.0)
    cd_betab = c_double(-9.0)

    R_beam = 100.
    R_diff = 100.
    cosz   = np.cos(0.0)

    ground_albedo_diff = 0.3
    ground_albedo_beam = 0.3
    
    iret = grndsnow_albedo_call(c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = grndsnow_albedo_call(c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))
    iret = canopy_prep_call(ci(ib))
    iret = zenith_prep_call(c8(cosz))
    iret = solver_call(c8(R_beam),c8(R_diff))

    for i in range(n_layer):
        
        ican = i+1
        icol = 1
        for iv in range(n_vai):
            iret = getintens_call(ci(ican),ci(icol),c8(elems[0][i].avai[iv]),byref(cd_r_diff_dn),byref(cd_r_diff_up),byref(cd_r_beam))
            elems[0][i].r_dn[iv] = cd_r_diff_dn.value
            elems[0][i].r_up[iv] = cd_r_diff_up.value
            elems[0][i].r_b[iv] = cd_r_beam.value
            print(elems[0][i].r_up[iv])
            if(iv>0):
                elems[0][i].r_abs[iv-1] = (elems[0][i].r_dn[iv]-elems[0][i].r_dn[iv-1]) + \
                    (elems[0][i].r_up[iv-1]-elems[0][i].r_up[iv]) \
                    (elems[0][i].r_b[iv]-elems[0][i].r_b[iv-1])


    #

    fig, axs = plt.subplots(ncols=2,nrows=n_layer,figsize=(8,8))
    ax1s = axs.reshape(-1)
    
    for i in range(n_layer):

        #ax = plt.subplot2grid( shape, loc,rowspan,colspan)
        #ax = plt.subplot2grid((
        #ax = plt.subplot2grid((2, 2), (0, 0))


        #elems[0][i].r_dn[iv]
        ax = ax1s[ic]
        ap = ax.plot(elems[0][i].r_dn,elems[0][i].avai)
        ax.invert_yaxis()
        ax.set_xlabel('')
        ax.set_ylabel('Integrated VAI [m2/m2]')
        ax.set_title('Beam Intensity [W/m2]')
        ax.grid(True)

        ic=ic+1
        ax = ax1s[ic]
        ap = ax.plot(elems[1][i].r_dn,elems[1][i].avai)
        ax.invert_yaxis()
        ax.set_xlabel('')
        ax.set_ylabel('Integrated VAI [m2/m2]')
        ax.set_title('Beam Intensity [W/m2]')
        ax.grid(True)

        ic=ic+1
   
    plt.show()
    

                
def SingleElementPerturbTest():

    
    # ===================================================================================
    n_col    = 1
    n_layer  = 1
    iret = alloc_twostream_call(ci(n_layer),ci(n_col))

    ican = 1 # Single canopy layer
    icol = 1 # Single PFT
    pft  = 1 # Use PFT number 1
    area = 0.9  # Assume only 90% of the ground is covered
    lai  = 2.0  # LAI
    sai  = 0.5  # SAI
    vai = lai+sai
    iret = setup_canopy_call(c_int(ican),c_int(icol),c_int(pft),c_double(area),c_double(lai),c_double(sai))


    # Decide on a band:

    ib = visb
    
    cd_r_beam = c_double(-9.0)
    cd_r_diff_up = c_double(-9.0)
    cd_r_diff_dn = c_double(-9.0)
    cd_kb = c_double(-9.0)
    cd_kd = c_double(-9.0)
    cd_om = c_double(-9.0)
    cd_betad = c_double(-9.0)
    cd_betab = c_double(-9.0)
    

    # Make parameter pertubations, bump up 50%
    pp_dict = {}
    pp_dict['Kb'] = 0.74*1.5
    pp_dict['Kd'] = 1.03*1.5
    pp_dict['om'] = 0.18*1.5
    pp_dict['betab'] = 0.45*1.5
    pp_dict['betad'] = 0.6*1.5

    R_beam = 100.
    R_diff = 0.
    cosz   = np.cos(0.0)
    n_vai  = 100
    vai_a  = np.linspace(0,vai,num=n_vai)

    dv = vai/n_vai
    
    r_diff_up = np.zeros(n_vai)   
    r_diff_dn = np.zeros(n_vai)
    r_beam    = np.zeros(n_vai)

    drdv_diff_up = np.zeros(n_vai-1) # Delta
    drdv_diff_dn = np.zeros(n_vai-1) # Delta
    drdv_ubeam    = np.zeros(n_vai-1) # Delta
    drdv_dbeam    = np.zeros(n_vai-1) # Delta
    
    p_r_diff_up = np.zeros([n_vai,len(pp_dict)])
    p_r_diff_dn = np.zeros([n_vai,len(pp_dict)])
    p_r_beam    = np.zeros([n_vai,len(pp_dict)])
    p_drdv_diff_up = np.zeros([n_vai-1,len(pp_dict)])
    p_drdv_diff_dn = np.zeros([n_vai-1,len(pp_dict)])
    p_drdv_ubeam    = np.zeros([n_vai-1,len(pp_dict)])
    p_drdv_dbeam    = np.zeros([n_vai-1,len(pp_dict)])

    ground_albedo_diff = 0.3
    ground_albedo_beam = 0.3
    
    iret = grndsnow_albedo_call(c_double(ground_albedo_diff),*ccharnb('albedo_grnd_diff'))
    iret = grndsnow_albedo_call(c_double(ground_albedo_beam),*ccharnb('albedo_grnd_beam'))
    iret = canopy_prep_call(ci(ib))
    iret = zenith_prep_call(c8(cosz))
    iret = solver_call(c8(R_beam),c8(R_diff))
    iret = getparams_call(ci(ican),ci(icol),byref(cd_kb),byref(cd_kd),byref(cd_om),byref(cd_betad),byref(cd_betab))

    for iv in range(n_vai):
        iret = getintens_call(ci(ican),ci(icol),c8(vai_a[iv]),byref(cd_r_diff_dn),byref(cd_r_diff_up),byref(cd_r_beam))
        r_beam[iv] = cd_r_beam.value
        r_diff_up[iv] = cd_r_diff_up.value
        r_diff_dn[iv] = cd_r_diff_dn.value

        if(iv>0):
            drdv_ubeam[iv-1] = -cd_om.value*cd_betab.value*(r_beam[iv]-r_beam[iv-1])/dv
            drdv_dbeam[iv-1] = -cd_om.value*(1.-cd_betab.value)*(r_beam[iv]-r_beam[iv-1])/dv
            drdv_diff_dn[iv-1] = -(r_diff_dn[iv]-r_diff_dn[iv-1])/dv
            drdv_diff_up[iv-1] = (r_diff_up[iv]-r_diff_up[iv-1])/dv
            
    # Redo the scattering with perturbations
    i = -1
    for key,val in pp_dict.items():
        i=i+1
        iret = canopy_prep_call(ci(ib))
        iret = zenith_prep_call(c8(cosz))
        iret = forceparam_call(c_int(ican),c_int(icol),c_double(val),*ccharnb(key))
        iret = solver_call(c8(R_beam),c8(R_diff))
        
        for iv in range(n_vai):
            iret = getintens_call(ci(ican),ci(icol),c8(vai_a[iv]),byref(cd_r_diff_dn),byref(cd_r_diff_up),byref(cd_r_beam))
            #print(iv,i,cd_r_beam.value)
            p_r_beam[iv,i] = cd_r_beam.value
            p_r_diff_up[iv,i] = cd_r_diff_up.value
            p_r_diff_dn[iv,i] = cd_r_diff_dn.value

            if(iv>0):
                p_drdv_ubeam[iv-1] = -cd_om.value*cd_betab.value*(p_r_beam[iv]-p_r_beam[iv-1])/dv
                p_drdv_dbeam[iv-1] = -cd_om.value*(1.-cd_betab.value)*(p_r_beam[iv]-p_r_beam[iv-1])/dv
                p_drdv_diff_dn[iv-1] = -(p_r_diff_dn[iv]-p_r_diff_dn[iv-1])/dv
                p_drdv_diff_up[iv-1] = (p_r_diff_up[iv]-p_r_diff_up[iv-1])/dv


        # Derivative Plot
        if(False):
            fig0, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(9,7))

            ap = ax1.plot(drdv_dbeam,vai_a[1:],p_drdv_dbeam[:,i],vai_a[1:])
            last_color = ap[-1].get_color()
            ax1.invert_yaxis()
            ax1.set_xlabel('')
            ax1.set_ylabel('Integrated VAI [m2/m2]')
            ax1.set_title('dIb/dv Beam Intensity [W/m2/m2]')
            ax1.grid(True)
            
            ax2.plot(drdv_diff_dn,vai_a[1:],p_drdv_diff_dn[:,i],vai_a[1:])
            ax2.invert_yaxis()
            ax2.set_xlabel('')
            ax2.set_yticklabels('')
            ax2.set_ylabel('')
            ax2.set_title('dIdn/dv Down Diffuse Intensity [W/m2/m2] ')
            ax2.grid(True)

            ax3.plot(drdv_ubeam,vai_a[1:],p_drdv_ubeam[:,i],vai_a[1:])
            ax3.invert_yaxis()
            ax3.set_xlabel('')
            ax3.set_ylabel('Integrated VAI [m2/m2]')
            ax3.set_title('dIb/dv Beam Intensity [W/m2/m2]')
            ax3.grid(True)
            
            ax4.plot(drdv_diff_up,vai_a[1:],p_drdv_diff_up[:,i],vai_a[1:])
            ax4.invert_yaxis()
            ax4.set_xlabel('')
            ax4.set_ylabel('Integrated VAI [m2/m2]')
            ax4.set_title('dIup/dv Up Diffuse Intensity [W/m2/m2]')
            ax4.grid(True)
        
                
        fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(9,7))

        ap = ax1.plot(r_beam,vai_a,p_r_beam[:,i],vai_a)
        first_color = ap[0].get_color()
        last_color = ap[-1].get_color()
        ax1.invert_yaxis()
        ax1.set_xlabel('')
        ax1.set_ylabel('Integrated VAI [m2/m2]')
        ax1.set_title('Beam Intensity [W/m2]')
        ax1.grid(True)

        ax2.plot(r_diff_dn,vai_a,p_r_diff_dn[:,i],vai_a)
        ax2.invert_yaxis()
        ax2.set_xlabel('')
        ax2.set_yticklabels('')
        ax2.set_ylabel('')
        ax2.set_title('Down Diffuse Intensity [W/m2] ')
        ax2.grid(True)
        
        ax3.plot(r_diff_up,vai_a,p_r_diff_up[:,i],vai_a)
        ax3.invert_yaxis()
        ax3.set_xlabel('')
        ax3.set_ylabel('Integrated VAI [m2/m2]')
        ax3.set_title('Up Diffuse Intensity [W/m2]')
        ax3.grid(True)

        ax4.axis("off")
        ax4.set_axis_off()

        if(ib==visb):
            band_name = "Visible"
        elif(ib==nirb):
            band_name = "Near Infrared"
        else:
            print("Unknown band")
            exit(2)
        
            
        param_str = r"""In-element Scattering Profiles

Broad band: {0}
$R_{{b,atm}} = ${1:.0f}
$R_{{d,atm}} = ${2:.0f}
$cos(\phi) = ${3:.2f}
$K_b = ${4:.2f}
$K_d = ${5:.2f} 
$\omega = ${6:.2f} 
$\beta_b = ${7:.2f}
$\beta_d = ${8:.2f}
$\alpha_{{gd}} = ${9:.2f}
$\alpha_{{gb}} = ${9:.2f}""".format(band_name,R_beam,R_diff,cosz,cd_kb.value,cd_kd.value,cd_om.value,cd_betab.value,cd_betad.value,ground_albedo_diff,ground_albedo_beam)    
        ax4.text(0.1, 0.5, param_str, horizontalalignment='left', \
                 verticalalignment='center', transform=ax4.transAxes,backgroundcolor=[1.0,1.0,1.0],fontsize=12,color=first_color)
        ax4.text(0.5,0.5,r"{0}={1:.2f}".format(key,val),color=last_color)
        plt.subplots_adjust(wspace=0.1, hspace=0.25)
        plt.show()


    dealloc_twostream_call()
    
    
# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
