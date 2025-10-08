import matplotlib
import matplotlib.pyplot as plt
import torch
import torch.multiprocessing as mp
import torch.distributed as dist
from torch.utils.data import DataLoader, Dataset
from torch.nn.parallel import DistributedDataParallel as DDP
import numpy as np
import argparse
import time
import os
import sys
import xml.etree.ElementTree as et
import code  # For development: code.interact(local=dict(globals(), **locals()))
from datetime import datetime

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

# Simple class to hold hyper-parameters
class hyper_parameters:
    def __init__(self, learning_rate, batch_size, loss_threshold):
        self.learning_rate = learning_rate
        self.batch_size = batch_size
        self.loss_threshold = loss_threshold
        
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


class PhotoNeuralNetwork(torch.nn.Module):
    def __init__(self,n_input,n_output,n_hidden):
        super(PhotoNeuralNetwork, self).__init__()
        self.nametag = 'its complicated'
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
    
def DDPRankTrain(rank, world_size, hyper_params, shared_inputs, shared_outputs, \
                 shared_inputs_mean, shared_inputs_std, \
                 shared_outputs_mean, shared_outputs_std, numleaf):

    """
    Worker function to initialize DDP and train the model using the shared dataset.
    """

    # 1. Setup Environment and Initialize Process Group
    os.environ['MASTER_ADDR'] = 'localhost'
    os.environ['MASTER_PORT'] = '29500'

    if (torch.cuda.is_available()):
        torch.cuda.set_device(rank)
        dist.init_process_group(
            backend="nccl",
            rank=rank,
            world_size=world_size)
        use_gpu = True
        device = torch.device("cuda:0")
        gpu_id = rank
    else:
        dist.init_process_group(
        backend="gloo",
        rank=rank,
        world_size=world_size)
        use_gpu = False
        device = "cpu"
        gpu_id = None
        
    print(f"RANK:[{rank}/{world_size}] PHASE 2: DDP group initialized.")

    # We don't need to normalize the output data because we normalize during the loss function
    
    #shared_outputs = (shared_outputs-shared_outputs.mean(dim=0))/shared_outputs.std(dim=0)
    
    dataset = CustomDataset(shared_inputs, shared_outputs)
    model_out_std = shared_outputs.std(dim=0)

    data_loader = DataLoader(dataset, batch_size=hyper_params.batch_size, shuffle=False, pin_memory=True, \
                             sampler=torch.utils.data.distributed.DistributedSampler(dataset))
    
    n_input = shared_inputs.size(1)
    n_output = shared_outputs.size(1)

    # We define the size of the hidden layers
    # as multiples of the number of independent leaf layers
    
    n_mult1 = 32
    n_mult2 = 16
    n_mult3 = 8

    n_hidden1 = numleaf*n_mult1
    n_hidden2 = numleaf*n_mult2
    n_hidden3 = numleaf*n_mult3
    
    mask1 = torch.zeros(n_hidden1,  n_input,   dtype=torch.bool)
    mask2 = torch.zeros(n_hidden2,  n_hidden1, dtype=torch.bool)
    mask3 = torch.zeros(n_hidden3,  n_hidden2, dtype=torch.bool)
    mask4 = torch.zeros(n_output,   n_hidden3, dtype=torch.bool)

    mask1 = DefineMask(mask1,n_shared,numleaf,n_per_leaflayer,n_mult1).float().to(device)
    mask2 = DefineMask(mask2,0,       numleaf,n_mult1,n_mult2).float().to(device)
    mask3 = DefineMask(mask3,0,       numleaf,n_mult2,n_mult3).float().to(device)
    mask4 = DefineMask(mask4,0,       numleaf,n_mult3,1).float().to(device)

    model = PhotoNeuralNetwork(n_input,n_output,[n_hidden1,n_hidden2,n_hidden3]).to(device)

    if(use_gpu): 
        ddp_model = torch.nn.parallel.DistributedDataParallel(model, device_ids=[gpu_id])
    else:
        ddp_model = torch.nn.parallel.DistributedDataParallel(model, device_ids=[])

    assert mask1.shape == ddp_model.module.fc1.weight.data.shape, \
        f"Shape mismatch: {mask1.shape} != {ddp_model.module.fc1.weight.data.shape}"
    assert mask2.shape == ddp_model.module.fc2.weight.data.shape, \
        f"Shape mismatch: {mask2.shape} != {ddp_model.module.fc2.weight.data.shape}"
    assert mask3.shape == ddp_model.module.fc3.weight.data.shape, \
        f"Shape mismatch: {mask3.shape} != {ddp_model.module.fc3.weight.data.shape}"
    assert mask4.shape == ddp_model.module.fc4.weight.data.shape, \
        f"Shape mismatch: {mask4.shape} != {ddp_model.module.fc4.weight.data.shape}"
    

    #criterion = norm_mse_loss #nn.MSELoss()
    #criterion = norm_huber_loss
    criterion = torch.nn.HuberLoss(reduction='mean', delta=1.0)
    #criterion = torch.nn.MSELoss()
    optimizer = torch.optim.AdamW(ddp_model.parameters(), lr=hyper_params.learning_rate)

    # Train the network
    loss_history = []
    print_interval = 10
    save_interval = 50

    start_time = time.time()
    for epoch in range(500000):
        global STOP_SIGNAL
        STOP_SIGNAL.zero_() 
        STOP_SIGNAL = STOP_SIGNAL.to(device)
        data_loader.sampler.set_epoch(epoch)
        for batch_idx, (batch_inputs, batch_outputs) in enumerate(data_loader):
            optimizer.zero_grad()

            ddp_model.module.fc1.weight.data.mul_(mask1)
            ddp_model.module.fc2.weight.data.mul_(mask2)
            ddp_model.module.fc3.weight.data.mul_(mask3)
            ddp_model.module.fc4.weight.data.mul_(mask4)
            
            pred_outputs = ddp_model(batch_inputs)
            #loss = criterion(pred_outputs, batch_outputs, model_out_std)
            loss = criterion(pred_outputs, batch_outputs)
            loss.backward()
            optimizer.step()

            ddp_model.module.fc1.weight.data.mul_(mask1)
            ddp_model.module.fc2.weight.data.mul_(mask2)
            ddp_model.module.fc3.weight.data.mul_(mask3)
            ddp_model.module.fc4.weight.data.mul_(mask4)
            
        loss_history.append(loss.item())

        if(rank==0 and loss.item() < hyper_params.loss_threshold):
            print(f"Rank {rank}: Loss {loss.item():.4f} is below threshold {hyper_params.loss_threshold}. Sending STOP signal.")
            # Set the signal to 1 (True)
            STOP_SIGNAL[0] = 1 
                
        # 6. Broadcast the decision to all ranks
        # This is the CRITICAL DDP synchronization step
        dist.broadcast(tensor=STOP_SIGNAL, src=0)

        if STOP_SIGNAL.item() == 1:
            print(f"Rank {rank}: Received STOP due to acceptable loss, epoch: {epoch+1}.")
            break # Exit the inner (batch/iteration) loop
        
        if(rank==0 and epoch % print_interval == 0):
            elapsed_time = time.time() - start_time
            print(f'Epoch: {epoch} Loss: {loss.item():.6f}, elapse: {elapsed_time:.4f}')
            start_time = time.time()
            
        if(rank==0 and epoch % save_interval == 0):
            ckp = ddp_model.module.state_dict()
            datestr = datetime.now().strftime("%Y%m%d-%H")
            torch.save(ckp,"checkpoint_{}.pt".format(datestr))
            
            
    if(rank==0):

        # Scale the weights and biases of the first layer
        model.NormScale(shared_inputs_mean,shared_inputs_std,shared_outputs_mean,shared_outputs_std)
        
        # give the model run a unique string
        datestr = datetime.now().strftime("%Y%m%d-%H%M")

        script_module = torch.jit.script(model)
        mod_pattern = model.nametag#'11-L64-Re-L32-Re-2'  # This is a label for model architecture in plotting
        script_module.save("./c3psn_modelsd_szv2_i13_{}_c{}.pt".format(mod_pattern,datestr))

    print(f"RANK:[{rank}/{world_size}] PHASE 2: Training complete")

    # Generate some scatter plots    

    if(rank==0):

        # Lets un-normalize the input data

        # new  = (old-mean)/std
        # old = new*std+mean
        shared_inputs = shared_inputs*shared_inputs_std + shared_inputs_mean
        shared_outputs = shared_outputs*shared_outputs_std + shared_outputs_mean
        
        indices_of_all_data = torch.randperm(shared_inputs.size(0))
        n_plot = 10000
        
        plot_mod  = model(shared_inputs[indices_of_all_data[0:n_plot]]).detach().numpy()
        plot_data = shared_outputs[indices_of_all_data[0:n_plot]].detach().numpy()
        
        #plot_data = batch_outputs.detach().numpy()
        #plot_mod  = pred_outputs.detach().numpy()
        
        
        fig, ((ax1,ax2)) = plt.subplots(1,2,figsize=(8.5,4.))
        ax1.scatter(plot_data[:,0],plot_mod[:,0])
        ax1.set_xlabel('FATES Ag(1) [umol/m2/s]')
        ax1.set_ylabel('NN Ag(1) [umol/m2/s]')
        minax = np.min([plot_data[:,0],plot_mod[:,0]])
        maxax = np.max([plot_data[:,0],plot_mod[:,0]])
        rngax = maxax-minax
        minax = minax-0.1*rngax
        maxax = maxax+0.1*rngax
        ax1.set_xlim([minax,maxax])
        ax1.set_ylim([minax,maxax])
        ax1.text(minax+0.05*(maxax-minax), \
                 minax+0.80*(maxax-minax), \
                 f'epoch: {epoch+1}\nlr: {hyper_params.learning_rate:.4f}\nmodel: {mod_pattern}' , \
                 bbox=dict(facecolor=[0.95,0.95,0.95], edgecolor='black'))
        ax1.grid(True)

        if(True):#plot_data.size(1)>1):
            plot_data[:,1] = plot_data[:,1] #*1.e-6  # Convert from umol to mol
            plot_mod[:,1] = plot_mod[:,1]    #*1.e-6
            ax2.scatter(plot_data[:,1],plot_mod[:,1])
            ax2.set_xlabel('FATES Ag(2) [mol/m2/s]')
            ax2.set_ylabel('NN Ag(2) [mol/m2/s]')
            minax = np.min([plot_data[:,1],plot_mod[:,1]])
            maxax = np.max([plot_data[:,1],plot_mod[:,1]])
            rngax = maxax-minax
            minax = minax-0.1*rngax
            maxax = maxax+0.1*rngax
            ax2.set_xlim([minax,maxax])
            ax2.set_ylim([minax,maxax])
            ax2.grid(True)
            
        plt.tight_layout()
        plt.show()

def PrepBar(bar_inputs):
    n_large = 20
    counts, bin_edges = np.histogram(bar_inputs, bins=n_large)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = (bin_edges[1] - bin_edges[0])*0.9
    
    
    return bin_centers, counts, bin_width

def ViewTrainingData(shared_inputs,shared_outputs):

    #matplotlib.use('Agg')
    fig, ((ax1,ax2,ax3,ax4), (ax5,ax6,ax7,ax8), \
          (ax9,ax10,ax11,ax12),(ax13,ax14,ax15,ax16)) = plt.subplots(4,4,figsize=(9.,9.))

    print("Generating fig")
    n_large = 20

    bin_centers, counts, bin_width = PrepBar(shared_inputs[:,0])
    ax1.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax1.set_title('PAR abs (umol/m2/s)')

    ax2.hist(shared_inputs[:,1], bins=n_large)
    bin_centers, counts, bin_width = PrepBar(shared_inputs[:,1])
    ax2.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax2.set_title('Vcmax25top (umol/m2/s)')

    bin_centers, counts, bin_width = PrepBar(shared_inputs[:,2])
    ax3.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax3.set_title('BTRAN (-)')

    bin_centers, counts, bin_width = PrepBar(shared_inputs[:,3])
    ax4.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax4.set_title('Tleaf (K)')

    bin_centers, counts, bin_width = PrepBar(shared_inputs[:,4])
    ax5.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax5.set_title('Patm (Pa)')

    bin_centers, counts, bin_width = PrepBar(shared_inputs[:,5])
    ax6.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax6.set_title('P_co2 (Pa)')

    bin_centers, counts, bin_width = PrepBar(shared_inputs[:,6])
    ax7.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax7.set_title('Esat (Pa)')

    bin_centers, counts, bin_width = PrepBar(shared_inputs[:,7])
    ax8.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax8.set_title('gb (umol/m2/s)')

    bin_centers, counts, bin_width = PrepBar(shared_inputs[:,8])
    ax9.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax9.set_title('E (Pa)')
    
    bin_centers, counts, bin_width = PrepBar(shared_inputs[:,9])
    ax10.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax10.set_title('Nscaler (-)')

    bin_centers, counts, bin_width = PrepBar(shared_inputs[:,10])
    ax11.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax11.set_title('LMR (umol/m2/s)')

    bin_centers, counts, bin_width = PrepBar(shared_outputs[:,0])
    ax12.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax12.set_title('Ag (umol/m2/s)')
    
    bin_centers, counts, bin_width = PrepBar(shared_outputs[:,1])
    ax13.bar(bin_centers, counts, width=bin_width, color='skyblue', edgecolor='black')
    ax13.set_title('gs (umol/m2/s)')

    plt.tight_layout()
    plt.show()
    
# ==================================================================================================



    
if __name__ == '__main__':

    
    default_learning_rate = 0.01
    default_loss_threshold = 0.001
    
    
    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--numproc',dest='numproc', type=int, \
                        help="Define how many processes (world_size) to run.", required=True)
    parser.add_argument('--learning_rate',dest='learning_rate', type=float, \
                        help="The learning rate (scale of training data)", required=False, default=default_learning_rate)
    parser.add_argument('--global_batch_size',dest='global_batch_size', type=int, \
                        help="The total batch size across all processes", required=True)
    parser.add_argument('--loss_threshold',dest='loss_threshold', type=float, \
                        help="The acceptable loss threshold to converge", required=False, default=default_loss_threshold)
    parser.add_argument('--numleaf',dest='numleaf', type=int, \
                        help="How many leaf layers to solve simultaneously?",required=True)

    
    args = parser.parse_args()
    world_size = args.numproc
    global_batch_size = args.global_batch_size
    learning_rate = args.learning_rate
    local_batch_size = int(float(global_batch_size)/float(world_size))
    loss_threshold = args.loss_threshold
    numleaf = args.numleaf

    hyper_params = hyper_parameters(learning_rate,local_batch_size,loss_threshold)
    
    
    print(f"\n------- Training Initiated --------\n\n")
    print(f"Global batch size: {global_batch_size}")
    print(f"World size (process count): {world_size}")
    print(f"Local (per process) batch size: {local_batch_size}")
    print(f"Acceptable loss threshold: {loss_threshold}")
    print(f"Learning rate: {learning_rate}\n\n")
    
    # Create the shared tensor in the main process.
    # The '.share_memory_()' method makes its memory accessible to all child processes.
    #shared_tensor = torch.zeros(1, 4, dtype=torch.float32).share_memory_()

    #n_trainset = int(2**20)  # ~1M

    n_trainset = int(2**15)   # ~32k
    
    # LETS DO 10 layers with sunlit/shaded 

    # How many inputs do we have that are specific to the leaf layer?
    # 1) Par Absorbed, 2) Nitrogen content (or would LMR work just as well? It should..
    #                     because it is based off of N content and the other variables
    #                     test this...

    # We may need to predict two outputs because if we also want the conductance,
    # we need to get the conductance slope terms. Or we can run that algorithm
    # to feed them into the model? Instead of btran? The calculations of gs0,gs1 and gs2
    # are only dependent on the parameter constants (pft) and btran.
    
    
    
    n_infeatures = n_shared + numleaf*n_per_leaflayer
    
    shared_inputs = torch.zeros([n_trainset,n_infeatures],dtype=torch.float32).share_memory_()
    shared_outputs = torch.zeros([n_trainset,n_outfeatures*numleaf],dtype=torch.float32).share_memory_()

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
        p = mp.Process(target=RankPrepInput, args=(rank, chunk, fates_path,shared_inputs,shared_outputs,numleaf))
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
    
    # The shared tensor is now fully populated.

    print(f"Final shared tensor after all processes have written: {shared_inputs[:,0]}")

    # Training

    mp.spawn(
        DDPRankTrain, 
        args=(world_size, hyper_params, shared_inputs, shared_outputs, \
              shared_inputs_mean, shared_inputs_std, \
              shared_outputs_mean, shared_outputs_std, numleaf),
        nprocs=world_size)


    print("\n--- All Execution Complete ---")


    

    

