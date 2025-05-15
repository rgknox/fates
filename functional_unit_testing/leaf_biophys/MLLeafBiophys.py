# =======================================================================================
#
# For usage: $python LeafBiophysDriver.py --help
#
# This script runs unit tests on the leaf biophysics functions
#
#
# =======================================================================================

#
# Note on conda environment
#
# conda env list
# 
# conda environments:
#
# base                     /home/rgknox/local/anaconda3
# py38                  *  /home/rgknox/local/anaconda3/envs/py38
# pytorch-py38             /home/rgknox/local/anaconda3/envs/pytorch-py38

# conda activate pytorch-py38



import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime
import argparse
from matplotlib.backends.backend_pdf import PdfPages
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
import random
import ctypes
from ctypes import *
from operator import add
sys.path.append('../shared/py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb
import CtypesInit
from PushParameters import PushParameters
from PushParameters import PushXMLPhotoParameters
from PushParameters import GetParamFromAttrib
from PushParameters import GetParamList

import torch
from torch import nn
from torch.utils.data import DataLoader
from torchvision import datasets
from torchvision.transforms import ToTensor

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 12}

matplotlib.rc('font', **font)


# Training a model to predict
# 1) g_s Stomatal Conductance
# 2) a_g Gross assimilation
#
# We train one model per set of model switches
#
#

device = "cpu"

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

# Simple conversion, number of micro-moles in a mole
umol_per_mol = 1.e6
mol_per_umol = 1.e-6

# 1 standard atmosphere in [Pa]
can_press_1atm = 101325.0

# Atmospheric CO2 partial pressure [Pa] at 400 ppm
co2_ppress_400ppm = 0.0004*can_press_1atm

# Atmospheric O2 partial pressure [Pa] %29.5 of atmosphere
o2_ppress_209kppm = 0.2095*can_press_1atm

# Respiration scaler at canopy top
rdark_scaler_top = 1.0

# Nitrogen scaler at canopy top
nscaler_top = 1.0


# Create aliases for the ctype Fortran objects
# =======================================================================================

exec(open("../shared/py_src/CtypesInit.py").read())

# Subroutines
# =======================================================================================
def nan_hook(self, inp, output):
    if not isinstance(output, tuple):
        outputs = [output]
    else:
        outputs = output

    for i, out in enumerate(outputs):
        nan_mask = torch.isnan(out)
        if nan_mask.any():
            print("In", self.__class__.__name__)
            raise RuntimeError(f"Found NAN in output {i} at indices: ", nan_mask.nonzero(), "where:", out[nan_mask.nonzero()[:, 0].unique(sorted=True)])

        
# Define the model
class Net(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, output_size)

    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
        return out

    
def Normalize2DTensor(x,dim=0,x_mean=None,x_std=None):    

    # If this is called for validation data
    # we want to use the mean and std from the training
    # set
    
    if(x_mean is None):
        x_mean = x.mean(dim=0)
    if(x_std is None):
        x_std  = x.std(dim=0)
        
    x_norm = (x - x_mean) / x_std
    return x_mean,x_std,x_norm

def DeNormalize2DTensor(x_norm,x_mean,x_std):

    #x_norm = (x - x_mean) / x_std
    #x = x_norm*x_std+x_mean
    
    return x_norm*x_std+x_mean


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



def main(argv):

    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--fin', dest='xmlfile', type=str, help="XML control file, required.", required=True)
    parser.add_argument('--smoketest', action='store_true')
    args = parser.parse_args()

    # This class call instanteates all the fortran shared objects
    # and creates aliases for their functions and subroutines
    f90 = f90_modules('../shared/bld/')

    
    # Load the xml control file
    # -----------------------------------------------------------------------------------
    xmlroot = et.parse(args.xmlfile).getroot()
    
    # We will allocate 1 token pft to hold data, we will change the values as needed
    numpft = 1

    # Push scalar parameters
    print('Pushing parameters from the xml file to the f90 lb_params datastructure')
    scalar_root = xmlroot.find('f90_params').find('scalar_dim')
    for param in scalar_root.iter('param'):
        print('pushing: '+param.attrib['name'].strip())
        iret = f90.set_leaf_param_sub(c8(float(param.text.split(',')[0])),ci(0),*ccharnb(param.attrib['name'].strip()))

        
    PushXMLPhotoParameters(f90,xmlroot)
    
    # Push pft parameters to fortran instantiations
    pft_root = xmlroot.find('f90_params').find('pft_dim')
    leaf_c3psn = []
    leaf_stomatal_intercept = []
    for param in pft_root.iter('param'):
        print('pushing: '+param.attrib['name'].strip())
        for pft in range(numpft):
            iret = f90.set_leaf_param_sub(c8(float(param.text.split(',')[pft])),ci(pft+1),*ccharnb(param.attrib['name'].strip()))
            if(param.attrib['name'].strip() == 'fates_leaf_c3psn'):
                leaf_c3psn.append(int(param.text.split(',')[pft]))
            if(param.attrib['name'].strip() == 'fates_leaf_stomatal_intercept'):
                leaf_stomatal_intercept.append(int(param.text.split(',')[pft]))
                
    # Read in non-fortran parameters from the xml file
    fates_leaf_vcmax25top    = []
    fates_stoich_nitr = []
    fates_leaf_slatop = []
    
    print('Reading non-fortran pft parameters')
    py_pft_root = xmlroot.find('py_params').find('pft_dim')
    for param in py_pft_root.iter('param'):
        for pft in range(numpft):
            if (param.attrib['name']=='fates_leaf_vcmax25top'):
                fates_leaf_vcmax25top.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_leaf_slatop'):
                fates_leaf_slatop.append(np.float64(param.text.split(',')[pft]))
            if (param.attrib['name']=='fates_stoich_nitr'):
                fates_stoich_nitr.append(np.float64(param.text.split(',')[pft]))
                
    print('Reading non-fortran scalar parameters')
    py_scalar_root = xmlroot.find('py_params').find('scalar_dim')
    for param in py_scalar_root.iter('param'):
        if (param.attrib['name']=='fates_maintresp_leaf_model'):
            fates_maintresp_leaf_model = int(param.text.split(',')[0])


    # Axes to test:
    #   Temperature
    #   PAR
    #   Humidity
    #   Boundary Layer conductance
    #   BTRAN
    #     vcmax (via btran)
    #     stomatal intercept (via btran)
    #     stomatal slope 1   (via btran)
    #     stomatal slope 2   (via btran/medlyn)
    
    # Axes implicitly tested
    # daylength factors (scale vcmax, but covered by btran scaling)
    # 

    # Leaf temperature ranges [C]
    leaf_tempc_min = -50.0
    leaf_tempc_max = 60.0
    leaf_tempc_n = 5
    leaf_tempc_vec = np.linspace(leaf_tempc_min,leaf_tempc_max,num=leaf_tempc_n)

    # Relative Humidity Ranges
    rh_max = 1.00
    rh_min = 0.1
    rh_n   = 5
    rh_vec = np.linspace(rh_min,rh_max,num=rh_n)
    
    # Absorbed PAR ranges [W/m2]
    par_abs_min = 10.0
    par_abs_max = 250
    par_abs_n  = 5
    par_abs_vec = np.linspace(par_abs_min,par_abs_max,num=par_abs_n)

    # Boundary Conductance ranges [umol/m2/s]
    gb_min =  500000.0            # Lower limit imposed by CLM/ELM 0.5 mol/m2/s
    gb_max = 5000000.0            # 50% larger than  Roughly largestthe largest values seen at BCI (which are 2.5mol/m2/s)
    gb_n  = 5
    gb_vec = np.linspace(gb_min,gb_max,num=gb_n)

    # btran ranges

    btran_n   = 5
    btran_min = -2
    btran_max = 0
    btran_vec = np.logspace(btran_min,btran_max,num=btran_n)

    # vcmax25top ranges
    vcmax25t_n = 15
    vcmax25t_min = 2
    vcmax25t_max = 125
    vcmax25t_vec = np.linspace(vcmax25t_min,vcmax25t_max,num=vcmax25t_n)

    # air pressure

    air_press_max = 1.1*can_press_1atm
    air_press_min = 0.85*can_press_1atm
    air_press_n   = 5
    air_press_vec = np.linspace(air_press_min,air_press_max,air_press_n)
    

    # Set convergence tolerance
    ci_tol = 0.1
    
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
    solve_iter_f = c_int(-9)
    gs0_f        = c_double(-9.0)
    gs1_f        = c_double(-9.0)
    gs2_f        = c_double(-9.0)

    pfails = 0
    ip = 0
    printfail = True
    
    # unit conversion of W/m2 to umol photons/m^2/s
    wm2_to_umolm2s = 4.6
  
    print('\nStarting Model training')
    print(' {} leaf temperature values [C] from {} to {}'.format(leaf_tempc_n,leaf_tempc_min,leaf_tempc_max))
    print(' {} RH values [fraction] from {} to {}'.format(rh_n,rh_min,rh_max))
    print(' {} PAR Abs [W/m2] values from {} to {}'.format(par_abs_n,par_abs_min,par_abs_max))
    print(' {} BL conductance (gb) [umol/m2/s] values from {} to {}'.format(gb_n,gb_min,gb_max))
    print(' {} BTRAN values [fraction] from {} to {}'.format(btran_n,np.exp(btran_min),np.exp(btran_max)))
    print(' {} Vcmax 25 top values [umol/m2/s] from {} to {}'.format(vcmax25t_n,vcmax25t_min,vcmax25t_max))
    

    c3_path_index = 1
    iret = f90.set_leaf_param_sub(c8(float(c3_path_index)),ci(1),*ccharnb('fates_leaf_c3psn'))
    
    btran_on_gs_gs2  = 4       # apply btran to the whole non-intercept portion
    iret = f90.set_leaf_param_sub(c8(float(btran_on_gs_gs2)),ci(1),*ccharnb('fates_leaf_stomatal_btran_model'))

    btran_on_ag_vcmax_jmax = 2 # apply btran to vcmax and jmax
    iret = f90.set_leaf_param_sub(c8(float(btran_on_ag_vcmax_jmax)),ci(1),*ccharnb('fates_leaf_agross_btran_model'))

    medlyn_model = 2
    ballberry_model = 1
    iret = f90.set_leaf_param_sub(c8(float(medlyn_model)),ci(1),*ccharnb('fates_leaf_stomatal_model'))

    
    n_model_runs = len(vcmax25t_vec)*len(leaf_tempc_vec)*len(btran_vec)*len(gb_vec)*len(rh_vec)*len(par_abs_vec)

    
    # Report every 5%
    ntestmod = int(n_model_runs/20)

    n_features_in  = 12   # number of controls on photosynthesis
    n_features_out = 2    # number of predictions (ie agross and gstoma)

    model_in  = np.zeros([n_model_runs,n_features_in])
    model_out = np.zeros([n_model_runs,n_features_out])
    
    print('\nRunning a total of {} tests: \n'.format(n_model_runs))
    time0 = time.process_time()

    for vcmax25_top in vcmax25t_vec:
                    
        jmax25_top,kp25_top =  GetJmaxKp25Top(vcmax25_top)
                    
        for leaf_tempc in leaf_tempc_vec:
                        
            leaf_tempk = leaf_tempc + tfrz_1atm
            
            iret = f90.qsat_sub(c8(leaf_tempk),c8(can_press_1atm), \
                                byref(veg_qs_f),byref(veg_es_f), \
                                byref(qsdt_dummy_f),byref(esdt_dummy_f))
            
            iret = f90.cangas_sub(c8(can_press_1atm), \
                                  c8(o2_ppress_209kppm), \
                                  c8(leaf_tempk), \
                                  byref(mm_kco2_f), \
                                  byref(mm_ko2_f), \
                                  byref(co2_cpoint_f))
        
            # Leaf Nitrogen Concentration at the top
            lnc_top  = fates_stoich_nitr[0]/fates_leaf_slatop[0]

            # Leaf Maintenance Respiration (temp and pft dependent)
            if(fates_maintresp_leaf_model==1):
                iret = f90.lmr_ryan_sub(c8(lnc_top),c8(nscaler_top), ci(1), c8(leaf_tempk), byref(lmr_f))
            elif(fates_maintresp_leaf_model==2):
                iret = f90.lmr_atkin_sub(c8(lnc_top),c8(rdark_scaler_top),c8(leaf_tempk),c8(atkin_mean_leaf_tempk),byref(lmr_f) )
            else:
                print('unknown leaf respiration model')
                exit(1)

            for btran in btran_vec:
                iret = f90.biophysrate_sub(ci(1), \
                                           c8(vcmax25_top), c8(jmax25_top), c8(kp25_top), \
                                           c8(nscaler_top), c8(leaf_tempk), c8(dayl_factor_full), \
                                           c8(t_growth_kum), c8(t_home_kum), c8(btran), \
                                           byref(vcmax_f), byref(jmax_f), byref(kp_f), byref(gs0_f), byref(gs1_f), byref(gs2_f))

                for gb in gb_vec:
                    for par_abs in par_abs_vec:
                        par_abs_umol = par_abs*wm2_to_umolm2s
                            
                        for rh in rh_vec:
                            vpress = rh * veg_es_f.value
                            

                            # We don't need to add the pft to the training
                            # because we assume this model is c3, has fnps=0.15
                            # and the typical metabolic rates, everything else
                            # is an argument.
                                
                            model_in[ip,0] = par_abs_umol
                            model_in[ip,1] = vcmax_f.value
                            model_in[ip,2] = jmax_f.value
                            #model_in[ip,3] = kp_f.value                 #
                            #model_in[ip,4] = gs0_f.value                #
                            model_in[ip,3] = gs2_f.value                
                            model_in[ip,4] = leaf_tempk
                            #model_in[ip,7] = can_press_1atm             #
                            #model_in[ip,8] = co2_ppress_400ppm          #
                            model_in[ip,5] = veg_es_f.value
                            model_in[ip,6] = gb
                            model_in[ip,7] = vpress
                            model_in[ip,8] = mm_kco2_f.value
                            model_in[ip,9] = mm_ko2_f.value
                            model_in[ip,10] = co2_cpoint_f.value
                            model_in[ip,11] = lmr_f.value
                                
                            try:
                                iret = f90.leaflayerphoto_sub(c8(par_abs_umol),  \
                                                              ci(1),   \
                                                              c8(vcmax_f.value),   \
                                                              c8(jmax_f.value),    \
                                                              c8(kp_f.value),      \
                                                              c8(gs0_f.value), \
                                                              c8(gs1_f.value), \
                                                              c8(gs2_f.value), \
                                                              c8(leaf_tempk), \
                                                              c8(can_press_1atm), \
                                                              c8(co2_ppress_400ppm), \
                                                              c8(o2_ppress_209kppm), \
                                                              c8(veg_es_f.value), \
                                                              c8(gb), \
                                                              c8(vpress), \
                                                              c8(mm_kco2_f.value), \
                                                              c8(mm_ko2_f.value), \
                                                              c8(co2_cpoint_f.value), \
                                                              c8(lmr_f.value), \
                                                              c8(ci_tol), \
                                                              byref(agross_f), \
                                                              byref(gstoma_f), \
                                                              byref(anet_f), \
                                                              byref(c13_f), \
                                                              byref(co2_interc_f), \
                                                              byref(solve_iter_f) )
                                
                                model_out[ip,0] = agross_f.value
                                model_out[ip,1] = gstoma_f.value
                                    
                            except:
                                pfails = pfails+1
                                printfail=True
                                                
                            if (np.mod(ip,ntestmod)==0):
                                print('Completed {} tests -- {} percent complete'.format(ip,100*float(ip)/float(n_model_runs)))

                            if (pfails>0 and np.mod(pfails,100)==0 and printfail):
                                printfail=False
                                print('\n{} fails so far\n'.format(pfails))
                                
                            ip = ip + 1             

                                                
    print("{} Failures out of {} Encountered; {}% of Tests\n".format(pfails,ip,float(pfails)/float(ip)))

    print('Deallocating parameter space')
    iret = f90.dealloc_leaf_param_sub()

    # Train an NN model
    # ------------------------------------------------------------------------------------------

    print('\n\nStarting model training process')
    print('Spliting model data into training and validation sets')
    
    # 1) Split the model in and output into a training set and a validation set
    # Remove validation points from the training sample

    valid_size = 5000
    ivs = [random.randint(0, n_model_runs-1) for _ in range(valid_size)]
    ivs.sort(reverse=True)
    valid_in  = np.zeros([valid_size,n_features_in])
    valid_out = np.zeros([valid_size,n_features_out])
    valid_in  = model_in[ivs,:]
    valid_out = model_out[ivs,:]
    
    train_size = n_model_runs-valid_size
    train_in  = np.zeros([train_size,n_features_in])
    train_out = np.zeros([train_size,n_features_out])
    ios = list(range(n_model_runs))
    for i in ivs:
        del ios[i]

    #for i,ii in enumerate(ios):
    train_in  = model_in[ios,:]
    train_out = model_out[ios,:]
    
    # 2) Generate input/output tensors, this includes normalization
    #    of data

    print('Creating NN tensors')
    # Generate pytorch tensors from photosynthesis code inputs and outputs
    x = torch.from_numpy(train_in.astype(np.float32))
    y = torch.from_numpy(train_out.astype(np.float32))

    x_mean,x_std,x_norm = Normalize2DTensor(x,dim=0)
    y_mean,y_std,y_norm = Normalize2DTensor(y,dim=0)
    
    print('Checking training data for nans')
    if(x_norm.isnan().any() or y_norm.isnan().any()):
        print('Found Nans in training data')
        code.interact(local=dict(globals(), **locals()))
    
    # Construct a model
    print('Constructing model')
    n_hidden  = 20   # number of perceptrons/neurons in the hidden layers
    #    nn.ReLU(),
    #    nn.Sigmoid,
    #    nn.Tanh
    model = nn.Sequential(
        nn.Linear(n_features_in,n_hidden),
        nn.Linear(n_hidden,n_hidden),
        nn.ReLU(),
        nn.Linear(n_hidden,n_hidden),
        nn.Linear(n_hidden,n_features_out),
    ).to(device)

    mod_pattern = 'LLRLL'  # This is a label for the above model
                           # used in plots
    
    learning_rate = 0.001
    criterion = nn.MSELoss()
    optimizer = torch.optim.SGD(model.parameters(),lr=learning_rate)
    num_epochs = 300000
    batch_size = 1000
    min_mse = 0.002
    print('Starting model passes, calculating losses and back-propogating')
    for epoch in range(num_epochs):

        ibs = [random.randint(0, train_size-1) for _ in range(batch_size)]
        x_batch = x_norm[ibs,:]
        y_batch = y_norm[ibs,:]
                    
        # forward pass and loss calculation
        loss = criterion(model(x_batch), y_batch)
        
        optimizer.zero_grad()
        
        # backward pass
        loss.backward()
        
        # update
        optimizer.step()

        if(epoch+1) % 1000 == 0:
            print(f'epoch: {epoch+1},loss = {loss.item():.4f}')

        if(loss.item()<=min_mse):
            break
            

    print('Re-running model with validation set and evaluating with original output')
    # Compare the validation set to model output
    x_valid = torch.from_numpy(valid_in.astype(np.float32))

    # Normalize input data based on training set mean, std
    xv_mean,xv_std,xv_norm = Normalize2DTensor(x_valid,dim=0,x_mean=x_mean,x_std=x_std)

    # Run the validation set through the model
    yv_norm = model(xv_norm)

    # Denormalize the validation set
    yv_pred = DeNormalize2DTensor(yv_norm,y_mean,y_std).detach().numpy()


    # Generate some scatter plots    
    
    fig, ((ax1,ax2)) = plt.subplots(1,2,figsize=(8.5,4.))
    ax1.scatter(valid_out[:,0],yv_pred[:,0])
    
    ax1.set_xlabel('FATES Ag [umol/m2/s]')
    ax1.set_ylabel('NN Ag [umol/m2/s]')
    minax = np.min([yv_pred[:,0],valid_out[:,0]])
    maxax = np.max([yv_pred[:,0],valid_out[:,0]])
    rngax = maxax-minax
    minax = minax-0.1*rngax
    maxax = maxax+0.1*rngax
    ax1.set_xlim([minax,maxax])
    ax1.set_ylim([minax,maxax])
    ax1.text(minax+0.05*(maxax-minax), \
             minax+0.80*(maxax-minax), \
             f'epoch: {epoch+1}\nlr: {learning_rate:.4f}\nN: {n_hidden}\nmodel: {mod_pattern}' , \
             bbox=dict(facecolor=[0.95,0.95,0.95], edgecolor='black'))
    ax1.grid(True)
    
    ax2.scatter(1.e-6*valid_out[:,1],1.e-6*yv_pred[:,1])
    ax2.set_xlabel('FATES gs [mol/m2/s]')
    ax2.set_ylabel('NN gs [mol/m2/s]')
    minax = np.min([1.e-6*yv_pred[:,1],1.e-6*valid_out[:,1]])
    maxax = np.max([1.e-6*yv_pred[:,1],1.e-6*valid_out[:,1]])
    rngax = maxax-minax
    minax = minax-0.1*rngax
    maxax = maxax+0.1*rngax
    ax2.grid(True)
    
    plt.tight_layout()
    plt.show()
    code.interact(local=dict(globals(), **locals()))
    
    print('Functional Unit Testing Complete')
    exit(0)

    
# =======================================================================================
# This is the actual call to main

if __name__ == "__main__":
    main(sys.argv)
