import matplotlib
import matplotlib.pyplot as plt
import torch
import torch.multiprocessing as mp
import torch.distributed as dist
from torch.utils.data import DataLoader, Dataset
from torch.nn.parallel import DistributedDataParallel as DDP
import numpy as np
import argparse
import pandas as pd
import time
import os
import sys
import xml.etree.ElementTree as et
import code  # For development: code.interact(local=dict(globals(), **locals()))
from datetime import datetime
# Ray Imports
import ray
from ray import tune
from ray.train.torch import TorchTrainer
from ray.tune.schedulers import ASHAScheduler

current_path = os.getcwd()
fates_path=current_path.split('fates')[0]+'fates'
sys.path.insert(0,fates_path+'/functional_unit_testing/shared/py_src')
exec(open(fates_path+"/functional_unit_testing/shared/py_src/CtypesInit.py").read())

from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb
import CtypesInit
from PushParameters import PushParameters
from PushParameters import PushXMLPhotoParameters
from PushParameters import GetParamFromAttrib
from PushParameters import GetParamList

# Number of shared parameters
n_shared = 7

# Number of dedicated parameter (par/nscaler)
n_per_leaflayer = 2

# Right now we are only predicing Ag, not gs
n_outfeatures = 1

STOP_SIGNAL = torch.zeros(1, dtype=torch.int32)

        
# LOSS FUNCTIONS
def MafeLoss(pred, target):
    eps = 1e-8
    return torch.mean(torch.abs((pred - target) / (torch.abs(target) + eps)))

def NormMSELoss(pred,target,std):
    return torch.mean( ((pred-target)/std)**2.0 )

def NormHuberLoss(pred,target,std):
    delta = 1.0
    a = torch.abs((pred-target)/std)
    quadratic_mask = a <= delta
    quadratic_loss = 0.5 * a**2

    # 2. Linear Region: |a| > delta
    # Loss is delta * (|a| - 0.5 * delta)
    linear_mask = a > delta
    linear_loss = delta * (a - 0.5 * delta)
    
    # Combine the losses based on the mask
    loss = torch.where(quadratic_mask, quadratic_loss, linear_loss)
    
    # Return the mean loss across the entire batch
    return torch.mean(loss)


def DefineMask(mask,n_shared,n_output,m_in,m_out):

    print_mask = False
    
    mask[:,:] = False
    
    for i in range(n_shared):
        mask[:,i] = True

    for i in range(n_output):
        j0 = n_shared+i*(m_in)
        j1 = j0+m_in
        for j in range(j0,j1):
            i0 = i*m_out
            i1 = i0+m_out
            mask[i0:i1,j] = True

    if(print_mask):
        print(f"Printing mask")
        for i in range(mask.size(0)):
            print(f"{mask[i,:]}")
        code.interact(local=dict(globals(), **locals()))
        
    return mask


def ReplicateModel(model_in,n_reps,n_shared):

    # Extend the model to handle a vector of inputs
    # Tell it how many replicates, and the number of inputs
    # that are replicated. It will assume the last inputs
    # are those to be replicated

    
    
    n_input  = model_in.fc1.weight.data.size(1)
    n_output = model_in.fc4.weight.data.size(0)
    n_hidden1 = model_in.fc1.weight.data.size(0)
    n_hidden2 = model_in.fc2.weight.data.size(0)
    n_hidden3 = model_in.fc3.weight.data.size(0)
    
    n_in_per_rep = n_input - n_shared

    model_out = PhotoNeuralNetwork( n_shared + n_reps*n_in_per_rep, \
                                    n_output*n_reps, \
                                    [n_reps*n_hidden1,n_reps*n_hidden2,n_reps*n_hidden3])

    model_out.fc1 = TransferWeightsBias(model_out.fc1,model_in.fc1,n_shared,n_reps)
    model_out.fc2 = TransferWeightsBias(model_out.fc2,model_in.fc2,0,n_reps)
    model_out.fc3 = TransferWeightsBias(model_out.fc3,model_in.fc3,0,n_reps)
    model_out.fc4 = TransferWeightsBias(model_out.fc4,model_in.fc4,0,n_reps)

    return model_out

    
def TransferWeightsBias(fc_out,fc_in,n_shared,n_reps):

    #
    # FC is being changed from [n_out x n_in] -> [nreps*n_out x n_shared+nreps*n_nonshared]
    #
    #
    #   Example 4 out, 3 in (2 shared) with 2 replicates
    #   x x x       x x x o 
    #   x x x       x x x o
    #   x x x  ->   x x x o
    #   x x x       x x x o
    #               x x o x
    #               x x o x
    #               x x o x
    #               x x o x

    #code.interact(local=dict(globals(), **locals()))
    
    fc_out.weight.data[:,:] = 0.
    fc_out.bias.data[:] = 0.

    n_input  = fc_in.weight.data.size(1)
    n_nonshared = n_input - n_shared

    n_output = fc_in.weight.data.size(0)

    # Transfer Biases
    for i in range(n_reps):
        i0 = i*n_output
        i1 = i0+n_output
        fc_out.bias.data[i0:i1] = fc_in.bias.data[:]

    # Transfer weights
    for i in range(n_reps):
        i0 = i*n_output
        i1 = i0+n_output
        for j in range(n_shared):
            fc_out.weight.data[i0:i1,j] = fc_in.weight.data[0:n_output,j]

        for jj in range(n_nonshared):
            j_in  = jj+n_shared
            j_out = jj+n_shared+i*n_nonshared
            fc_out.weight.data[i0:i1,j_out] = fc_in.weight.data[0:n_output,j_in]
        

    return fc_out
            
        
class PhotoNeuralNetwork(torch.nn.Module):
    def __init__(self,n_input,n_output,n_hidden):
        super(PhotoNeuralNetwork, self).__init__()
        self.nametag = 'tuning'
        self.fc1   = torch.nn.Linear(n_input, n_hidden[0])
        self.relu1 = torch.nn.ReLU()
        self.fc2   = torch.nn.Linear(n_hidden[0],n_hidden[1])
        self.fc3   = torch.nn.Linear(n_hidden[1],n_hidden[2])
        self.fc4   = torch.nn.Linear(n_hidden[2],n_output)
        
    def forward(self, x):
        x = self.relu1(self.fc1(x))
        x = self.fc2(x)
        x = self.fc3(x)
        x = self.fc4(x)
        return x

    def NormScale(self,in_mean,in_std,out_mean,out_std):

        # We apply this after the model has been
        # trained on normalized input data. It
        # essentially reverses the normalization
        # process, so inference does not need
        # to perform that step
        
        with torch.no_grad():

            # Lets scale the model's first linear
            # weights and biases by the normalization factor
            W = self.fc1.weight.data  # W is typically shape [output_size, input_size]
            B = self.fc1.bias.data    # B is typically shape [output_size]
            
            # 1. Update Weights: W_new = W / sigma
            # The unsqueeze(0) ensures 'stds' has shape [1, n_input_features] for correct division
            W_new = W.div(in_std.unsqueeze(0)) 
    
            # 2. Update Bias: B_new = B - Sum(W_i * mu_i / sigma_i)
            # The term to subtract is the vector-matrix product: W_i * (mu_i / sigma_i)
    
            # Calculate mu_i / sigma_i for each feature
            scale_term = in_mean.div(in_std) # shape [3]

            # Calculate the adjustment for the bias: Sum over the input features
            # W_new @ scale_term is a matrix multiplication that sums over the input features
            bias_adjustment = W.matmul(scale_term) # shape [output_size]

            B_new = B.sub(bias_adjustment)
    
            # 3. Overwrite the trained model parameters
            self.fc1.weight.data.copy_(W_new)
            self.fc1.bias.data.copy_(B_new)

            
            # Extract the trained tensors
            W = self.fc4.weight.data 
            B = self.fc4.bias.data

            W_new = W * out_std.unsqueeze(1) 
            B_new = (B * out_std) + out_mean

            self.fc4.weight.data.copy_(W_new)
            self.fc4.bias.data.copy_(B_new)
            
    
class CustomDataset(Dataset):
    def __init__(self, inputs, outputs):
        self.inputs = inputs
        self.outputs = outputs
        
    def __len__(self):
        return len(self.inputs)

    def __getitem__(self, idx):

        # Save this if you need to do normalization
        x = self.inputs[idx]
        y = self.outputs[idx]
        return x, y


def GetLeafTempc(n_samp):
    return np.random.normal(loc=302-273.4, scale=5, size=n_samp)
    # SCALE E2

def GetRH(n_samp):
    return np.random.uniform(low=0.1, high=1.0, size=n_samp)
    # SCALE E1

def GetCO2(n_samp):
    return np.random.uniform(low=390.0,high=410., size=n_samp)

def GetPress(n_samp):
    return np.random.normal(loc=100000., scale=1000, size=n_samp)

def GetALTempDiff(n_samp):
    # This is T_air - T_leaf
    #return np.random.normal(loc=0, scale = 10, size = n_samp)
    return np.random.uniform(low=-2.0, high=0., size=n_samp)

def GetNScaler(n_samp):
    return np.random.uniform(low=0.05, high=1.0, size = n_samp)

def GetBTran(n_samp):
    return np.random.normal(loc=0.6,scale=0.6,size=n_samp).clip(0.05,1.0)

def GetGB(n_samp):
    return np.random.normal(loc=1.e6,scale=0.5e6,size=n_samp).clip(0.2e6,None)

def GetPARAbsWatts(n_samp):
    # umol/m2/s
    return np.random.uniform(low=0,high=300.0, size = n_samp)

def GetVcmax25Top(n_samp):
    return np.random.normal(loc=60,scale=20,size=n_samp).clip(5.0,200.)

def GetLNCTop(n_samp):
    fates_stoich_nitr = [0.033, 0.029, 0.04, 0.033, 0.04, 0.04, \
                         0.033, 0.04, 0.04, 0.04, 0.04, 0.04]
    fates_leaf_slatop = [0.012, 0.005, 0.024, 0.009, 0.03, 0.03, \
                         0.012, 0.03, 0.03, 0.03, 0.03, 0.03]    
    lnc = np.array([fates_stoich_nitr[i]/fates_leaf_slatop[i] for i in range(len(fates_stoich_nitr))])
    return np.random.normal(loc=lnc.mean(),scale=lnc.std(),size=n_samp).clip(0.1,None)
    

def RankPrepInput(rank, chunk, fates_path, shared_inputs, shared_outputs, numleaf):

    print(f"RANK:[{rank}/{world_size}] PHASE 1: Starting input prep.")

    # Constants
    dayl_factor_full = 1.0
    kumgrowth_tempk = -999.9
    kumhome_tempk   = -999.9
    ci_tol = 0.1
    
    # This class call instantiates all the fortran shared objects
    # and creates aliases for their functions and subroutines
    f90 = f90_modules(fates_path+'/functional_unit_testing/shared/bld/')

    # Load the xml control file
    # -----------------------------------------------------------------------------------
    xmlroot = et.parse(fates_path+"/functional_unit_testing/leaf_biophys/leaf_biophys_controls.xml").getroot()

    numpft = 12
    iret = f90.alloc_leaf_param_sub(ci(numpft))

    photo_parameters = {
        'fates_daylength_factor_switch':1,
        'fates_leaf_stomatal_model':2,
        'fates_leaf_stomatal_assim_model': 1,
        'fates_leaf_photo_tempsens_model': 0,
        'fates_electron_transport_model': 1}   # ,  <!--FvCB1980 = 1, JohnsonBerry2021 = 2 -->
    
    photo_pft_parameters = {
        'fates_leaf_stomatal_btran_model': [4,4,4,4,4,4,4,4,4,4,4,4],
        'fates_leaf_agross_btran_model': [2,2,2,2,2,2,2,2,2,2,2,2 ],
        'fates_leaf_c3psn': [1,1,1,1,1,1,1,1,1,1,1,0],
        'fates_leaf_stomatal_slope_ballberry': [ 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 ],
        'fates_leaf_stomatal_slope_medlyn': [ 4.1, 2.3, 2.3, 4.1, 4.4, 4.4, 4.7, 4.7, 4.7, 2.2, 5.3, 1.6 ],
        'fates_leaf_fnps': [ 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15 ],
        'fates_leaf_stomatal_intercept': [ 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 40000 ],
        'fates_maintresp_leaf_atkin2017_baserate': [  1.756, 1.4995, 1.4995, 1.756, 1.756, 1.756, 2.0749, 2.0749, 2.0749, 2.1956, 2.1956, 2.1956 ],
        'fates_maintresp_leaf_ryan1991_baserate': [ 2.525e-06, 2.525e-06, 2.525e-06, 2.525e-06, 2.525e-06, 2.525e-06, 2.525e-06, 2.525e-06, 2.525e-06, 2.525e-06, 2.525e-06, 2.525e-06 ],
        'fates_leaf_vcmaxha': [ 65330, 65330, 65330, 65330, 65330, 65330, 65330, 65330, 65330, 65330, 65330, 65330 ],
        'fates_leaf_jmaxha': [ 43540, 43540, 43540, 43540, 43540, 43540, 43540, 43540, 43540, 43540, 43540, 43540 ],
        'fates_leaf_vcmaxhd': [  149250, 149250, 149250, 149250, 149250, 149250, 149250, 149250, 149250, 149250, 149250, 149250 ],
        'fates_leaf_jmaxhd': [ 149250, 149250, 149250, 149250, 149250, 149250, 149250, 149250, 149250, 149250, 149250, 149250 ],
        'fates_leaf_vcmaxse': [ 485, 485, 485, 485, 485, 485, 485, 485, 485, 485, 485, 485 ],
        'fates_leaf_jmaxse': [ 495, 495, 495, 495, 495, 495, 495, 495, 495, 495, 495, 495 ]}

    # These arrays are not pushed to the F90 data structures
    fates_leafn_vert_scaler_coeff1 = [0.00963, 0.00963, 0.00963, 0.00963, 0.00963, 0.00963, \
                                      0.00963, 0.00963, 0.00963, 0.00963, 0.00963, 0.00963]
    fates_leafn_vert_scaler_coeff2 = [2.43, 2.43, 2.43, 2.43, 2.43, 2.43, \
                                      2.43, 2.43, 2.43, 2.43, 2.43, 2.43]
    fates_leaf_vcmax25top = [50, 62, 39, 61, 58, 58, \
                             62, 54, 54, 78, 78, 78]
    fates_stoich_nitr = [0.033, 0.029, 0.04, 0.033, 0.04, 0.04, \
                         0.033, 0.04, 0.04, 0.04, 0.04, 0.04]
    fates_leaf_slatop = [0.012, 0.005, 0.024, 0.009, 0.03, 0.03, \
                         0.012, 0.03, 0.03, 0.03, 0.03, 0.03]

    
    # Push scalar parameters
    for key, value in photo_parameters.items():
        #print(f"Key: {key}, Value: {value}")
        iret = f90.set_leaf_param_sub(c8(float(value)),ci(0),*ccharnb(key))
        
    for key, values in photo_pft_parameters.items():
        #print(f"Key: {key}")
        for ft,val in enumerate(values):
            #print(f"Value: {val}")
            iret = f90.set_leaf_param_sub(c8(float(val)),ci(ft+1),*ccharnb(key))
        
    np.random.seed(rank)
    
    n_samp         = len(chunk)
    leaf_tempc_vec = GetLeafTempc(n_samp)
    rh_vec         = GetRH(n_samp)
    co2_ppm_vec    = GetCO2(n_samp)
    altdiff_vec    = GetALTempDiff(n_samp)
    btran_vec      = GetBTran(n_samp)
    gb_vec         = GetGB(n_samp)
    press_vec      = GetPress(n_samp)
    
    pft = 1
    ft  = pft-1  # Python index associated with PFT of interest

    vcmax25top = fates_leaf_vcmax25top[ft]
    jmax25top = 1.67 * vcmax25top
    kp25top   = 20000. * vcmax25top
    lnctop    = fates_stoich_nitr[ft]/fates_leaf_slatop[ft]
    
        
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
    
    for i,ip in enumerate(chunk):

        leaf_tempk = leaf_tempc_vec[i]+273.14
        air_tempk  = leaf_tempk + altdiff_vec[i]
        o2_ppress = 0.2095*press_vec[i]
        co2_ppress = co2_ppm_vec[i] * 1.e-6 * press_vec[i]
        
        iret = f90.qsat_sub(c8(air_tempk),c8(press_vec[i]), \
                            byref(veg_qs_f),byref(veg_es_f), \
                            byref(qsdt_dummy_f),byref(esdt_dummy_f))

        vpress_sat = veg_es_f.value
        vpress = rh_vec[i] * veg_es_f.value
        
        iret = f90.cangas_sub(c8(press_vec[i]), \
                              c8(o2_ppress), \
                              c8(leaf_tempk), \
                              byref(mm_kco2_f), \
                              byref(mm_ko2_f), \
                              byref(co2_cpoint_f))

        nscaler_vec    = GetNScaler(numleaf)
        parabs_vec     = GetPARAbsWatts(numleaf)*4.6

        shared_inputs[ip,0] = btran_vec[i]
        shared_inputs[ip,1] = leaf_tempc_vec[i]
        shared_inputs[ip,2] = press_vec[i]
        shared_inputs[ip,3] = co2_ppress
        shared_inputs[ip,4] = veg_es_f.value
        shared_inputs[ip,5] = gb_vec[i]
        shared_inputs[ip,6] = vpress
            
        for j in range(numleaf):
        
            iret = f90.lmr_ryan_sub(c8(lnctop),c8(nscaler_vec[j]), ci(pft), \
                                    c8(leaf_tempk), byref(lmr_f))

        
            iret = f90.biophysrate_sub(ci(pft), \
                                       c8(vcmax25top), c8(jmax25top), c8(kp25top), \
                                       c8(nscaler_vec[j]), c8(leaf_tempk), c8(dayl_factor_full), \
                                       c8(kumgrowth_tempk), c8(kumhome_tempk), c8(btran_vec[i]), \
                                       byref(vcmax_f), byref(jmax_f), byref(kp_f), byref(gs0_f), byref(gs1_f), byref(gs2_f))
        
            ip0 = 7+(j*2)
            shared_inputs[ip,ip0]   = nscaler_vec[j]
            shared_inputs[ip,ip0+1] = parabs_vec[j]
        
            try:
                # Call the FATES photosynthesis subroutine:
                # https://github.com/NGEET/fates/blob/main/biogeophys/LeafBiophysicsMod.F90#L1232
                iret = f90.leaflayerphoto_sub(c8(parabs_vec[j]),  \
                                              ci(pft),   \
                                              c8(vcmax_f.value),   \
                                              c8(jmax_f.value),    \
                                              c8(kp_f.value),      \
                                              c8(gs0_f.value), \
                                              c8(gs1_f.value), \
                                              c8(gs2_f.value), \
                                              c8(leaf_tempk), \
                                              c8(press_vec[i]), \
                                              c8(co2_ppress), \
                                              c8(o2_ppress), \
                                              c8(veg_es_f.value), \
                                              c8(gb_vec[i]), \
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

                shared_outputs[ip,j] = agross_f.value

            
            except:
                print('Photosynthesis model could not find a solution')
                exit(1)
        
    print(f"RANK:[{rank}/{world_size}] PHASE 1: Ending input prep.")


# 2. Define the main training function for a single trial
# This function is run on *each* distributed worker by Ray Train

def RayTrainPhoto(config):

    shard = ray.train.get_dataset_shard("train_key").materialize()
    
    # Get all column names in the shard. This triggers a necessary metadata fetch.
    all_columns = shard.columns() 
    
    # Use list comprehension and regex to filter for feature columns
    feature_columns = [col for col in all_columns if re.match(r"input_\d+", col)]
    target_columns = [col for col in all_columns if re.match(r"target_\d+", col)]
                       
    # 3. Convert to PyTorch DataLoader using the discovered column names
    #torch_loader = shard.to_torch(
    #    label_columns=target_columns,
    #    feature_columns=feature_columns,
    #    batch_size=config["batch_size"],
    #    label_column_dtypes=torch.float32 
    #)
    data_iterator = shard.iter_torch_batches(
        batch_size=config["batch_size"],
        # Specify all columns and their data types here
        dtypes={col: torch.float32 for col in feature_columns + target_columns} 
    )
    

    # 4. Model Setup (dynamically uses the discovered feature count)
    #train_loader = ray.train.torch.prepare_data_loader(torch_loader)               
    
    model = PhotoNeuralNetwork(len(feature_columns),len(target_columns), \
                               [config["hidden_size1"],config["hidden_size2"],config["hidden_size3"]])
    
    
    # Wrap the model for DistributedDataParallel (DDP)
    # Ray Train sets up the distributed process group, so we just wrap it
    model = ray.train.torch.prepare_model(model) 

    # Loss function and optimizer (tuned hyperparameters)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=config["lr"])

    # --- Training Loop ---
    for epoch in range(config["epochs"]):
        model.train()
        if hasattr(train_loader.sampler, "set_epoch"):             
            train_loader.sampler.set_epoch(epoch) # Important for DDP shuffling

        epoch_loss = 0.0
        for batch in data_iterator:
        #for X_batch, Y_batch in train_loader:
            targets = torch.cat(
                [torch.unsqueeze(batch[f], dim=1) for f in target_columns], dim=1
            )
            inputs = torch.cat(
                [torch.unsqueeze(batch[f], dim=1) for f in feature_columns], dim=1
            )
            
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        avg_loss = epoch_loss / len(dataloader)
        
        # Report the metric back to Ray Tune
        # Ray Train syncs this metric from all workers to the rank 0 process
        # before reporting to Tune.
        ray.train.report({"loss": avg_loss})

    # Optional: Save a final checkpoint
    ray.train.save_checkpoint(
        model_state_dict=model.state_dict(),
        epoch=config["epochs"]
    )


# ==================================================================================================

    
if __name__ == '__main__':

    
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--numproc',dest='numproc', type=int, \
                        help="Define how many processes (world_size) to run.", required=True)
    
    args = parser.parse_args()
    world_size = args.numproc
    
    print(f"\n------- Training Initiated --------\n\n")
    print(f"World size (process count): {world_size}")
    
    # Create the shared tensor in the main process.
    # The '.share_memory_()' method makes its memory accessible to all child processes.
    #shared_tensor = torch.zeros(1, 4, dtype=torch.float32).share_memory_()

    #n_trainset = int(2**20)  # ~1M
    n_trainset = int(2**16)
    
    n_validset = int(2**12)  # Validation set (vectorized version)
    
    #n_trainset = int(2**15)   # ~32k
    
    # LETS DO 10 layers with sunlit/shaded 

    # How many inputs do we have that are specific to the leaf layer?
    # 1) Par Absorbed, 2) Nitrogen content (or would LMR work just as well? It should..
    #                     because it is based off of N content and the other variables
    #                     test this...

    # We may need to predict two outputs because if we also want the conductance,
    # we need to get the conductance slope terms. Or we can run that algorithm
    # to feed them into the model? Instead of btran? The calculations of gs0,gs1 and gs2
    # are only dependent on the parameter constants (pft) and btran.
    
    n_infeatures = n_shared + n_per_leaflayer
    shared_inputs = torch.zeros([n_trainset,n_infeatures],dtype=torch.float32).share_memory_()
    shared_outputs = torch.zeros([n_trainset,n_outfeatures],dtype=torch.float32).share_memory_()

    # Create a list to hold the process objects.
    processes = []
    for rank in range(world_size):
        chunk_size = int(np.floor(n_trainset/world_size))
        ida = int(rank*chunk_size)
        if(rank==(world_size-1)):
            idz = int(n_trainset)
        else:
            idz = int((rank+1)*chunk_size)
        chunk = list(range(ida,idz))
        p = mp.Process(target=RankPrepInput, args=(rank, chunk, fates_path,shared_inputs,shared_outputs,1))
        p.start()
        processes.append(p)
    
    # Wait for all processes to complete their work.
    for p in processes:
        p.join()


    if(False):
        print("About to make plot")
        ViewTrainingData(shared_inputs,shared_outputs)

    # Lets normalize
    shared_inputs_mean = shared_inputs.mean(dim=0)
    shared_inputs_std  = shared_inputs.std(dim=0)
    shared_inputs  = (shared_inputs-shared_inputs_mean)/shared_inputs_std

    shared_outputs_mean = shared_outputs.mean(dim=0)
    shared_outputs_std  = shared_outputs.std(dim=0)
    shared_outputs  = (shared_outputs-shared_outputs_mean)/shared_outputs_std
    
    print(f"Final shared tensor after all processes have written: {shared_inputs[:,0]}")

    # Initialize Ray
    ray.init()

    # Ray wants a dictionary. Intererstingly, it also likes it when
    # each feature has its own entry

    train_data_dict = {
        f"input_{i}": shared_inputs[:, i].numpy() for i in range(shared_inputs.shape[1])}

    for i in range(shared_outputs.shape[1]):
        train_data_dict[f"target_{i}"] = shared_outputs[:,i].numpy()

    #code.interact(local=dict(globals(), **locals()))
    # Ray wants things in a dataframe
    train_dataset = ray.data.from_pandas(pd.DataFrame(train_data_dict))
    #code.interact(local=dict(globals(), **locals()))
    #shard = ray.train.get_dataset_shard("train_key")

    # Get all column names in the shard. This triggers a necessary metadata fetch.
    #    all_columns = shard.columns()

    
                               
    # Define the hyperparameter search space
    search_space = {
        "hidden_size1": tune.choice([32, 16, 8]),
        "hidden_size2": tune.choice([32, 16, 8]),
        "hidden_size3": tune.choice([32, 16, 8]),
        "lr": tune.loguniform(1e-4, 1e-2),
        "batch_size": tune.choice([32, 128, 512, 1024, 2048]),
        "epochs": 20, # Small number for example speed
    }

    #search_space = {
        # The 'config' key holds all the parameters that will be passed
        # to your actual training function or model builder.
    #    "config": model_hparams,
        
        # You might also have trainer-specific arguments here, like
        # "scaling_config": tune.grid_search([...])
    #}

    # Use Ray Train to manage the distributed training environment (DDP)
    # num_workers specifies the number of DDP processes (GPUs or CPUs)
    trainer = TorchTrainer(
        train_loop_per_worker = RayTrainPhoto,
        scaling_config=ray.train.ScalingConfig(num_workers=world_size, use_gpu=False),
        run_config=ray.air.RunConfig(name="ddp_tune_regression"),
        datasets={"train_key":train_dataset},
    )

    # Configure the hyperparameter search with an ASHA scheduler
    scheduler = ASHAScheduler(
        metric="loss",
        mode="min",
        max_t=search_space["epochs"],
        grace_period=1,
    )

    tuner = tune.Tuner(
        trainer,
        param_space={"train_loop_config": search_space},
    )

    # Start the hyperparameter search
    results = tuner.fit()
    
    best_result = results.get_best_result("loss", "min")
    #code.interact(local=dict(globals(), **locals()))
    print(f"\nBest config found: {best_result.config}")
    #print(f"Best validation loss: {best_result.metrics['loss']:.4f}")

    
    
    print("\n--- All Execution Complete ---")


    

    

