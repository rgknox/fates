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
import importlib.util

# THIS SCRIPT SHOULD BE RUN WITH TORCHRUN!
# e.g.:
# torchrun --nproc_per_node=12 python PhotoTrain.py

current_path = os.getcwd()
fates_path=current_path.split('fates')[0]+'fates'
sys.path.insert(0,fates_path+'/functional_unit_testing/shared/py_src')

from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb
import CtypesInit
from PushParameters import PushParameters
from PushParameters import PushXMLPhotoParameters
from PushParameters import GetParamFromAttrib
from PushParameters import GetParamList

import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset
import torch.distributed as dist
import torch.optim as optim
import torch.multiprocessing as mp
from torch.nn.parallel import DistributedDataParallel as DDP
    
#from torchvision import datasets
#from torchvision.transforms import ToTensor
#from torch.optim.lr_scheduler import ReduceLROnPlateau

import pandas as pd
from ipywidgets import widgets
from IPython.display import display


#! (cd /home/rgknox/Models/CTSM/src/fates/functional_unit_testing/shared/;./build_fates_objs.sh)

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 12}

matplotlib.rc('font', **font)

print("IS CUDA AVAILABLE?:{}".format(torch.cuda.is_available()))

#code.interact(local=dict(globals(), **locals()))

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

# Atmospheric O2 partial pressure [Pa] %29.5 of atmosphere
o2_ppress_209kppm = 0.2095*can_press_1atm


# Respiration scaler at canopy top
rdark_scaler_top = 1.0

# Nitrogen scaler at canopy top
nscaler_top = 1.0

# unit conversion of W/m2 to umol photons/m^2/s
wm2_to_umolm2s = 4.6

# Set convergence tolerance
ci_tol = 0.1

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

if(dist.get_rank()==1):
    
    # Create aliases for the ctype Fortran objects
    # =======================================================================================

    exec(open(fates_path+"/functional_unit_testing/shared/py_src/CtypesInit.py").read())


    # This class call instantiates all the fortran shared objects
    # and creates aliases for their functions and subroutines
    f90 = f90_modules(fates_path+'/functional_unit_testing/shared/bld/')


    controls = {
        "fates_daylength_factor_switch": [1],
        "fates_leaf_stomatal_btran_model": [4],
        "fates_leaf_agross_btran_model": [2],
        "fates_leaf_c3psn": [1],
        "fates_maintresp_leaf_model": [1] }

    parameter_constants = {
        "fates_leaf_fnps": [0.15],
        "fates_maintresp_leaf_atkin2017_baserate": [1.756],
        "fates_maintresp_leaf_ryan1991_baserate": [2.525e-06],
        "fates_leaf_vcmaxha": [65330],
        "fates_leaf_jmaxha": [43540],
        "fates_leaf_vcmaxhd": [149250],
        "fates_leaf_jmaxhd": [149250],
        "fates_leaf_vcmaxse": [485],
        "fates_leaf_jmaxse": [495],
        "fates_leafn_vert_scaler_coeff1": [0.00963],
        "fates_leafn_vert_scaler_coeff2": [2.43],
        "fates_stoich_nitr": [0.033],
        "fates_leaf_slatop": [0.012]}



    # Load the xml control file
    # -----------------------------------------------------------------------------------
    xmlroot = et.parse(fates_path+"/functional_unit_testing/leaf_biophys/leaf_biophys_controls.xml").getroot()

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

    c3_path_index = 1
    iret = f90.set_leaf_param_sub(c8(float(c3_path_index)),ci(1),*ccharnb('fates_leaf_c3psn'))

    btran_on_gs_gs2  = 4       # apply btran to the whole non-intercept portion
    iret = f90.set_leaf_param_sub(c8(float(btran_on_gs_gs2)),ci(1),*ccharnb('fates_leaf_stomatal_btran_model'))

    btran_on_ag_vcmax_jmax = 2 # apply btran to vcmax and jmax
    iret = f90.set_leaf_param_sub(c8(float(btran_on_ag_vcmax_jmax)),ci(1),*ccharnb('fates_leaf_agross_btran_model'))

    medlyn_model = 2
    ballberry_model = 1
    iret = f90.set_leaf_param_sub(c8(float(medlyn_model)),ci(1),*ccharnb('fates_leaf_stomatal_model'))

    n_small = 5
    n_large = 8

    # Leaf temperature ranges [C]
    leaf_tempc_min = 15.0
    leaf_tempc_max = 50.0
    leaf_tempc_n = n_large
    leaf_tempc_vec = np.linspace(leaf_tempc_min,leaf_tempc_max,num=leaf_tempc_n)
    
    # Relative Humidity Ranges
    rh_max = 1.00
    rh_min = 0.1
    rh_n   = n_large
    rh_vec = np.linspace(rh_min,rh_max,num=rh_n)

    # CO2 concentration ranges (ppm)
    co2_max = 450.0
    co2_min = 250.0
    co2_n   = n_small
    co2_vec = np.linspace(co2_min,co2_max,num=co2_n)

    # Atmospheric Pressure ranges
    can_press_min = 090000.0
    can_press_max = 110000.0
    can_press_n   = n_small
    can_press_vec = np.linspace(can_press_min,can_press_max,can_press_n)

    # Absorbed PAR ranges [W/m2]
    par_abs_min = 0.0
    par_abs_max = 300
    par_abs_n  = n_large
    par_abs_vec = np.linspace(par_abs_min,par_abs_max,num=par_abs_n)

    # Boundary Conductance ranges [umol/m2/s]
    gb_min =  500000.0            # Lower limit imposed by CLM/ELM 0.5 mol/m2/s
    gb_max = 5000000.0            # 50% larger than  Roughly largestthe largest values seen at BCI (which are 2.5mol/m2/s)
    gb_n  = n_large
    gb_vec = np.linspace(gb_min,gb_max,num=gb_n)

    # btran ranges
    btran_n   = n_large
    btran_min = 0.01
    btran_max = 1
    btran_vec = np.linspace(btran_min,btran_max,num=btran_n)
    
    # vcmax25top ranges
    vcmax25t_n = n_large
    vcmax25t_min = 0.1
    vcmax25t_max = 60
    vcmax25t_vec = np.linspace(vcmax25t_min,vcmax25t_max,num=vcmax25t_n)

    print(' {} leaf temperature values [C] from {} to {}'.format(leaf_tempc_n,leaf_tempc_min,leaf_tempc_max))
    print(' {} RH values [fraction] from {} to {}'.format(rh_n,rh_min,rh_max))
    print(' {} PAR Abs [W/m2] values from {} to {}'.format(par_abs_n,par_abs_min,par_abs_max))
    print(' {} BL conductance (gb) [umol/m2/s] values from {} to {}'.format(gb_n,gb_min,gb_max))
    print(' {} BTRAN values [fraction] from {} to {}'.format(btran_n,np.exp(btran_min),np.exp(btran_max)))
    print(' {} Vcmax 25 top values [umol/m2/s] from {} to {}'.format(vcmax25t_n,vcmax25t_min,vcmax25t_max))

    # Initialize fortran return values
    # these are all temps and doubles
    vcmax_f      = c_double(-9.0);jmax_f       = c_double(-9.0)
    kp_f         = c_double(-9.0);agross_f     = c_double(-9.0)
    gstoma_f     = c_double(-9.0);anet_f       = c_double(-9.0)
    lmr_f        = c_double(-9.0);c13_f        = c_double(-9.0)
    co2_interc_f = c_double(-9.0);veg_qs_f     = c_double(-9.0)
    veg_es_f     = c_double(-9.0);mm_kco2_f    = c_double(-9.0)
    mm_ko2_f     = c_double(-9.0);co2_cpoint_f = c_double(-9.0)
    qsdt_dummy_f = c_double(-9.0);esdt_dummy_f = c_double(-9.0)
    solve_iter_f = c_int(-9);     gs0_f        = c_double(-9.0)
    gs1_f        = c_double(-9.0);gs2_f        = c_double(-9.0)
    
    n_model_runs = len(vcmax25t_vec)*len(leaf_tempc_vec)*len(btran_vec)*len(gb_vec)*len(rh_vec)*len(par_abs_vec)*len(co2_vec)*len(can_press_vec)

    print('\nRunning mechanistic photosynthesis model for {} combinations'.format(n_model_runs))

    # Report every 5%
    ntestmod = int(n_model_runs/20)

    n_features_in  = 14   # number of controls on photosynthesis
    n_features_out = 2    # number of predictions (ie agross and gstoma)

    model_in  = torch.from_numpy(np.zeros([n_model_runs,n_features_in]).astype(np.float32))
    model_out = torch.from_numpy(np.zeros([n_model_runs,n_features_out]).astype(np.float32))

    print('\nRunning a total of {} tests: \n'.format(n_model_runs))
    time0 = time.process_time()

    lnc_top  = fates_stoich_nitr[0]/fates_leaf_slatop[0]

    ip = 0
    for vcmax25_top in vcmax25t_vec:

  jmax25_top,kp25_top =  GetJmaxKp25Top(vcmax25_top)

  for leaf_tempc in leaf_tempc_vec:
    leaf_tempk = leaf_tempc + tfrz_1atm
    for can_press in can_press_vec:

      o2_ppress = 0.2095*can_press
        
      iret = f90.qsat_sub(c8(leaf_tempk),c8(can_press), \
                          byref(veg_qs_f),byref(veg_es_f), \
                          byref(qsdt_dummy_f),byref(esdt_dummy_f))

      iret = f90.cangas_sub(c8(can_press), \
                            c8(o2_ppress), \
                            c8(leaf_tempk), \
                            byref(mm_kco2_f), \
                            byref(mm_ko2_f), \
                            byref(co2_cpoint_f))

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
              for co2_ppm in co2_vec:
                co2_ppress = co2_ppm * 1.e-6 * can_press

                # We don't need to add the pft to the training
                # because we assume this model is c3, has fnps=0.15
                # and the typical metabolic rates, everything else
                # is an argument.

                model_in[ip,0] = par_abs_umol
                model_in[ip,1] = vcmax_f.value
                model_in[ip,2] = jmax_f.value
                model_in[ip,3] = gs2_f.value
                model_in[ip,4] = leaf_tempk
                model_in[ip,5] = can_press
                model_in[ip,6] = co2_ppress          #
                model_in[ip,7] = veg_es_f.value
                model_in[ip,8] = gb
                model_in[ip,9] = vpress
                model_in[ip,10] = mm_kco2_f.value
                model_in[ip,11] = mm_ko2_f.value
                model_in[ip,12] = co2_cpoint_f.value
                model_in[ip,13] = lmr_f.value


                try:
                  # Call the FATES photosynthesis subroutine:
                  # https://github.com/NGEET/fates/blob/main/biogeophys/LeafBiophysicsMod.F90#L1232
                  iret = f90.leaflayerphoto_sub(c8(par_abs_umol),  \
                                                ci(1),   \
                                                c8(vcmax_f.value),   \
                                                c8(jmax_f.value),    \
                                                c8(kp_f.value),      \
                                                c8(gs0_f.value), \
                                                c8(gs1_f.value), \
                                                c8(gs2_f.value), \
                                                c8(leaf_tempk), \
                                                c8(can_press), \
                                                c8(co2_ppress), \
                                                c8(o2_ppress), \
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
                  print('Photosynthesis model could not find a solution')
                  exit(1)

                if (np.mod(ip,ntestmod)==0):
                    print('Completed {} tests -- {} percent complete'.format(ip,100*float(ip)/float(n_model_runs)))

                ip = ip + 1

# Visualize the distributions of the input data

# 14 variables, 4x4 should work

if(False):
    fig, ((ax1,ax2,ax3,ax4), (ax5,ax6,ax7,ax8), \
          (ax9,ax10,ax11,ax12),(ax13,ax14,ax15,ax16)) = plt.subplots(4,4,figsize=(9.,9.))


    ax1.hist( model_in[:,0] ,bins=n_large)
    ax1.set_title('PAR (umol/m2/s)')
    ax2.hist(model_in[:,1], bins=n_large)
    ax2.set_title('Vcmax (umol/m2/s)')
    ax3.hist(model_in[:,2], bins=n_large)
    ax3.set_title('Jmax (umol/m2/s)')
    ax4.hist(model_in[:,3], bins=n_large)
    ax4.set_title('Gs2 (-)')
    ax5.hist(model_in[:,4], bins=n_large)
    ax5.set_title('Tl (K)')
    ax6.hist(model_in[:,5], bins=n_large)
    ax6.set_title('P_atm (Pa)')
    ax7.hist(model_in[:,6], bins=n_large)
    ax7.set_title('P_co2 (Pa)')
    ax8.hist(model_in[:,7], bins=n_large)
    ax8.set_title('E_sat (Pa)')
    
    ax9.hist(model_in[:,8], bins=n_large)
    ax9.set_title('g_b (umol/m2/s)')
    ax10.hist(model_in[:,9], bins=n_large)
    ax10.set_title('E (Pa)')
    ax11.hist(model_in[:,10], bins=n_large)
    ax11.set_title('Kco2')
    ax12.hist(model_in[:,11], bins=n_large)
    ax12.set_title('Ko2')
    ax13.hist(model_in[:,12], bins=n_large)
    ax13.set_title('Co2 Cpoint (Pa)')
    ax14.hist(model_in[:,13], bins=n_large)
    ax14.set_title('LMR (umol/m2/s')
    
    plt.tight_layout()
    plt.show()


# Lets export this to a dataset, the dataloader
# requires these methods:

class CustomDataset(Dataset):
    def __init__(self, inputs, outputs):
        self.inputs = inputs
        self.outputs = outputs
        self.in_mean = self.inputs.mean(dim=0)
        self.in_std = self.inputs.std(dim=0)
        self.out_mean = self.outputs.mean(dim=0)
        self.out_std = self.outputs.std(dim=0)
        
    def __len__(self):
        return len(self.inputs)

    def __getitem__(self, idx):

        # Save this if you need to do normalization
        #x = (self.inputs[idx] - self.in_mean) / self.in_std
        #y = (self.outputs[idx] - self.out_mean) / self.out_std
        x = self.inputs[idx]
        y = self.outputs[idx]
        return x, y


batch_size = 1024
dataset = CustomDataset(model_in, model_out)

model_out_std = model_out.std(dim=0)

data_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

class PhotoNeuralNetwork(nn.Module):
    def __init__(self):
        super(PhotoNeuralNetwork, self).__init__()
        self.fc1   = nn.Linear(14, 64)
        self.relu1 = nn.ReLU()
        self.fc2   = nn.Linear(64, 32)
        self.relu2 = nn.ReLU()
        self.fc3   = nn.Linear(32, 2)
        #self.relu3 = nn.ReLU()
        #self.fc4   = nn.Linear(8, 2)
        
    def forward(self, x):
        x = self.relu1(self.fc1(x))
        x = self.relu2(self.fc2(x))
        x = self.fc3(x)
        #x = self.relu3(x)
        #x = self.fc4(x)
        return x

    #def set_normalization_params(self, mean, std):
    #    # Method to update the mean and std after calculating from training data
    #    self.data_mean = mean
    #    self.data_std = std

mod_pattern = '13-L64-Re-L32-Re-2'  # This is a label for model architecture in plotting

model = PhotoNeuralNetwork().to(device)


def setup():
    os.environ['MASTER_ADDR'] = '127.0.0.1'
    os.environ['MASTER_PORT'] = '29500'

    dist.init_process_group(backend='gloo', init_method='env://')
    # initialize the process group
    #dist.init_process_group("gloo", rank=rank, world_size=world_size)

def cleanup():
    dist.destroy_process_group()
    

setup()
ddp_model = DDP(model, device_ids=[dist.get_rank()])


#model = Net()
# Mean Absolute Fractional Error (scale invariant)

def mafe_loss(pred, target):
    eps = 1e-8
    return torch.mean(torch.abs((pred - target) / (torch.abs(target) + eps)))

def norm_mse_loss(pred,target,std):
    return torch.mean( ((pred-target)/std)**2 )

criterion = norm_mse_loss #nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)

# Train the network
loss_history = []
for epoch in range(1000):  
    running_loss = 0.0
    for batch_idx, (inputs, data_out) in enumerate(data_loader):
        optimizer.zero_grad()
        model_out = ddp_model(inputs)
        loss = criterion(model_out, data_out, model_out_std )
        loss.backward()
        optimizer.step()
        running_loss += loss.item()
        loss_history.append(loss.item())

    print(f'Rank {dist.get_rank()} Train Epoch: {epoch} [{batch_idx*len(data)}/{len(train_loader.dataset)} ({100.*batch_idx/len(train_loader):.0f}%)]\tLoss: {loss.item():.6f}')

    #print('Epoch %d, loss = %.3f' %
    #(epoch + 1, loss.item() ))


cleanup()


if(dist.get_rank()==1):
    # give the model run a unique string
    datestr = datetime.now().strftime("%Y%m%d-%H%M")

    script_module = torch.jit.script(model)
    script_module.save("./c3psn_modelsd_szv2_i13_{}_c{}.pt".format(mod_pattern,datestr))


#traced_model = torch.jit.trace(model, dummy_input)
#    frozen_model = torch.jit.freeze(traced_model)
#    ## print(frozen_model.graph)
#    ## print(frozen_model.code)
#    frozen_model.save(filename)

    
#script_model = torch.jit.trace(model)

# Save the ScriptModule
#script_model.save("/content/c3psn_model_i13_{}_c{}.pt".format(mod_pattern,datestr))


