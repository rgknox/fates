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

STOP_SIGNAL = torch.zeros(1, dtype=torch.int32)

# LOSS FUNCTIONS
def mafe_loss(pred, target):
    eps = 1e-8
    return torch.mean(torch.abs((pred - target) / (torch.abs(target) + eps)))

def norm_mse_loss(pred,target,std):
    return torch.mean( ((pred-target)/std)**2 )


class PhotoNeuralNetwork(torch.nn.Module):
    def __init__(self):
        super(PhotoNeuralNetwork, self).__init__()
        self.fc1   = torch.nn.Linear(11, 64)
        self.relu1 = torch.nn.ReLU()
        self.fc2   = torch.nn.Linear(64, 32)
        self.relu2 = torch.nn.ReLU()
        self.fc3   = torch.nn.Linear(32, 2)
        #self.relu3 = torch.nn.ReLU()
        #self.fc4   = torch.nn.Linear(8, 2)
        
    def forward(self, x):
        x = self.relu1(self.fc1(x))
        x = self.relu2(self.fc2(x))
        x = self.fc3(x)
        #x = self.relu3(x)
        #x = self.fc4(x)
        return x

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


    
def GetLeafTempc(n_samp):
    return np.random.normal(loc=302-273.4, scale=5, size=n_samp)
    
def GetRH(n_samp):
    return np.random.uniform(low=0.0, high=1.0, size=n_samp)

def GetCO2(n_samp):
    return np.random.uniform(low=250.0,high=450., size=n_samp)

def GetPress(n_samp):
    return np.random.normal(loc=100000., scale=1000, size=n_samp)

def GetALTempDiff(n_samp):
    # This is T_air - T_leaf
    return np.random.normal(loc=0, scale = 10, size = n_samp)

def GetNScaler(n_samp):
    return np.random.uniform(low=0.01,high=1.0, size = n_samp)

def GetBTran(n_samp):
    return np.random.normal(loc=0.6,scale=0.6,size=n_samp).clip(0.0,1.0)

def GetGB(n_samp):
    return np.random.normal(loc=1.e6,scale=0.5e6,size=n_samp).clip(0.2e6,None)

def GetPARAbsUmol(n_samp):
    # umol/m2/s
    return np.random.uniform(low=0,high=300.0*4.6, size = n_samp)

def GetVcmax25Top(n_samp):
    return np.random.normal(loc=60,scale=20,size=n_samp).clip(0.,200.)

def GetLNCTop(n_samp):
    fates_stoich_nitr = [0.033, 0.029, 0.04, 0.033, 0.04, 0.04, \
                         0.033, 0.04, 0.04, 0.04, 0.04, 0.04]
    fates_leaf_slatop = [0.012, 0.005, 0.024, 0.009, 0.03, 0.03, \
                         0.012, 0.03, 0.03, 0.03, 0.03, 0.03]    
    lnc = np.array([fates_stoich_nitr[i]/fates_leaf_slatop[i] for i in range(len(fates_stoich_nitr))])
    return np.random.normal(loc=lnc.mean(),scale=lnc.std(),size=n_samp).clip(0.1,None)
    

def RankPrepInput(rank, chunk, fates_path, shared_inputs, shared_outputs):

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
    nscaler_vec    = GetNScaler(n_samp)
    btran_vec      = GetBTran(n_samp)
    gb_vec         = GetGB(n_samp)
    parabs_vec     = GetPARAbsUmol(n_samp)
    press_vec      = GetPress(n_samp)
    vcmax25top_vec = GetVcmax25Top(n_samp)
    lnctop_vec     = GetLNCTop(n_samp)
    pft = 1
    ft  = pft-1  # Python index associated with PFT of interest
    
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
        vcmax25top = vcmax25top_vec[i]
        jmax25top = 1.67 * vcmax25top
        kp25top   = 20000. * vcmax25top
        
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

        iret = f90.lmr_ryan_sub(c8(lnctop_vec[i]),c8(nscaler_vec[i]), ci(pft), \
                                c8(leaf_tempk), byref(lmr_f))

        
        iret = f90.biophysrate_sub(ci(pft), \
                                   c8(vcmax25top), c8(jmax25top), c8(kp25top), \
                                   c8(nscaler_vec[i]), c8(leaf_tempk), c8(dayl_factor_full), \
                                   c8(kumgrowth_tempk), c8(kumhome_tempk), c8(btran_vec[i]), \
                                   byref(vcmax_f), byref(jmax_f), byref(kp_f), byref(gs0_f), byref(gs1_f), byref(gs2_f))

        shared_inputs[ip,0] = parabs_vec[i]
        shared_inputs[ip,1] = vcmax25top_vec[i]
        shared_inputs[ip,2] = btran_vec[i]
        shared_inputs[ip,3] = leaf_tempk
        shared_inputs[ip,4] = press_vec[i]
        shared_inputs[ip,5] = co2_ppress
        shared_inputs[ip,6] = veg_es_f.value
        shared_inputs[ip,7] = gb_vec[i]
        shared_inputs[ip,8] = vpress
        shared_inputs[ip,9] = nscaler_vec[i]
        shared_inputs[ip,10] = lmr_f.value

        try:
            # Call the FATES photosynthesis subroutine:
            # https://github.com/NGEET/fates/blob/main/biogeophys/LeafBiophysicsMod.F90#L1232
            iret = f90.leaflayerphoto_sub(c8(parabs_vec[i]),  \
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

            shared_outputs[ip,0] = agross_f.value
            shared_outputs[ip,1] = gstoma_f.value
        except:
            print('Photosynthesis model could not find a solution')
            exit(1)
        
    print(f"RANK:[{rank}/{world_size}] PHASE 1: Ending input prep.")
    
def DDPRankTrain(rank, world_size, shared_inputs, shared_outputs):

    """
    Worker function to initialize DDP and train the model using the shared dataset.
    """

    #rank = dist.get_rank()
    batch_size = 1024
    
    # 1. Setup Environment and Initialize Process Group
    os.environ['MASTER_ADDR'] = 'localhost'
    os.environ['MASTER_PORT'] = '29500'
    
    dist.init_process_group(
        backend="gloo",
        rank=rank,
        world_size=world_size
    )
    print(f"RANK:[{rank}/{world_size}] PHASE 2: DDP group initialized.")

    dataset = CustomDataset(shared_inputs, shared_outputs)
    model_out_std = shared_outputs.std(dim=0)

    #sampler_ext = torch.utils.data.distributed.DistributedSampler(
    #    dataset, 
    #    num_replicas=world_size, 
    #    rank=rank, 
    #    shuffle=True
    #)
    
    data_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, \
                             sampler=torch.utils.data.distributed.DistributedSampler(dataset))

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    
    model = PhotoNeuralNetwork().to(device)

    #ddp_model = torch.nn.parallel.DistributedDataParallel(model, device_ids=[rank])
    # Use this for CPUs which don't have to provide a phyisocal device ID
    ddp_model = torch.nn.parallel.DistributedDataParallel(model, device_ids=[], output_device=None)

    learning_rate = 0.0001
    criterion = norm_mse_loss #nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

    LOSS_THRESHOLD = 0.1
    
    # Train the network
    loss_history = []
    print_interval = 100
    
    for epoch in range(50000):
        global STOP_SIGNAL
        STOP_SIGNAL.zero_() 
        STOP_SIGNAL = STOP_SIGNAL.to(device)
        running_loss = 0.0
        for batch_idx, (inputs, data_out) in enumerate(data_loader):
            optimizer.zero_grad()
            model_out = ddp_model(inputs)
            loss = criterion(model_out, data_out, model_out_std )
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
            loss_history.append(loss.item())

        if(rank==0 and loss.item() < LOSS_THRESHOLD):
            print(f"Rank {rank}: Loss {loss.item():.4f} is below threshold {LOSS_THRESHOLD}. Sending STOP signal.")
            # Set the signal to 1 (True)
            STOP_SIGNAL[0] = 1 
                
        # 6. Broadcast the decision to all ranks
        # This is the CRITICAL DDP synchronization step
        dist.broadcast(tensor=STOP_SIGNAL, src=0)

        if STOP_SIGNAL.item() == 1:
            print(f"Rank {rank}: Received STOP signal. Canceling epoch {epoch+1} at batch {batch_idx+1}.")
            break # Exit the inner (batch/iteration) loop
        
        if(rank==0 and epoch % print_interval == 0):
            print(f'Epoch: {epoch} Loss: {loss.item():.6f}')

        if(rank==0):
            ckp = ddp_model.module.state_dict()
            datestr = datetime.now().strftime("%Y%m%d-%H")
            torch.save(ckp,"checkpoint_{}.pt".format(datestr))

            
    if(rank==0):
        # give the model run a unique string
        datestr = datetime.now().strftime("%Y%m%d-%H%M")

        script_module = torch.jit.script(model)
        mod_pattern = '11-L64-Re-L32-Re-2'  # This is a label for model architecture in plotting
        script_module.save("./c3psn_modelsd_szv2_i13_{}_c{}.pt".format(mod_pattern,datestr))

    print(f"RANK:[{rank}/{world_size}] PHASE 2: Training complete")

    # Generate some scatter plots    

    if(rank==0):
        plot_data = data_out.detach().numpy()
        plot_mod  = model_out.detach().numpy()
        plot_data[:,1] = plot_data[:,1]*1.e-6
        plot_mod[:,1] = plot_mod[:,1]*1.e-6
        
        fig, ((ax1,ax2)) = plt.subplots(1,2,figsize=(8.5,4.))
        ax1.scatter(plot_data[:,0],plot_mod[:,0])
        ax1.set_xlabel('FATES Ag [umol/m2/s]')
        ax1.set_ylabel('NN Ag [umol/m2/s]')
        minax = np.min([plot_data[:,0],plot_mod[:,0]])
        maxax = np.max([plot_data[:,0],plot_mod[:,0]])
        rngax = maxax-minax
        minax = minax-0.1*rngax
        maxax = maxax+0.1*rngax
        ax1.set_xlim([minax,maxax])
        ax1.set_ylim([minax,maxax])
        ax1.text(minax+0.05*(maxax-minax), \
                 minax+0.80*(maxax-minax), \
                 f'epoch: {epoch+1}\nlr: {learning_rate:.4f}\nmodel: {mod_pattern}' , \
                 bbox=dict(facecolor=[0.95,0.95,0.95], edgecolor='black'))
        ax1.grid(True)

       
        
        ax2.scatter(plot_data[:,1],plot_mod[:,1])
        ax2.set_xlabel('FATES gs [mol/m2/s]')
        ax2.set_ylabel('NN gs [mol/m2/s]')
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


    

# ==================================================================================================

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Parse command line arguments to this script.')
    parser.add_argument('--numproc',dest='nproc', type=int, \
                        help="Define how many processes (world_size) to run.", required=True)
    args = parser.parse_args()

    world_size = args.nproc

    print("World size: {}".format(world_size))
    
    # Create the shared tensor in the main process.
    # The '.share_memory_()' method makes its memory accessible to all child processes.
    #shared_tensor = torch.zeros(1, 4, dtype=torch.float32).share_memory_()

    n_trainset = 50000
    n_infeatures = 11
    n_outfeatures = 2
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
        p = mp.Process(target=RankPrepInput, args=(rank, chunk, fates_path,shared_inputs,shared_outputs))
        p.start()
        processes.append(p)
    
    # Wait for all processes to complete their work.
    for p in processes:
        p.join()

    # Child processes are complete
        
    # The shared tensor is now fully populated.
    print(f"Final shared tensor after all processes have written: {shared_inputs[:,0]}")

    
    # Training

    mp.spawn(
        DDPRankTrain, 
        args=(world_size, shared_inputs, shared_outputs),
        nprocs=world_size,
        join=True)  # Blocks until all DDP workers have terminated


    print("\n--- All Execution Complete ---")


    

    

