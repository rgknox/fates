import torch
import torch.multiprocessing as mp
import torch.distributed as dist
import numpy as np
import argparse
import time
        
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
    return np.random.normal(loc=1.e6,scale=0.5e6,size=n_samp).clip(0,None)

def GetPARAbsUmol(n_samp):
    # umol/m2/s
    return np.random.uniform(low=0,high=300.0*4.6, size = n_samp)

def GetRH(n_samp):
    return np.random.uniform(low=0.1,high=1.0, size= n_samp)


def TestDistWorker(rank,test_tensor):
    test_tensor[0,rank] = rank
    return test_tensor
    
def RankPrepInput(rank, chunk, shared_inputs):

    # Every process needs to sample from a different

    print("Rank: {}, ida: {}, idz: {}".format(rank,chunk[0],chunk[-1]))

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
    
    #shared_inputs[ip,0] = rank

    for i,ip in enumerate(chunk):

        leaf_tempk = leaf_tempc_vec[i]+273.14
        
        iret = f90.qsat_sub(c8(leaf_tempk),c8(press_vec[i]), \
                            byref(veg_qs_f),byref(veg_es_f), \
                            byref(qsdt_dummy_f),byref(esdt_dummy_f))

        vpress_sat = veg_es_f.value
        
        iret = f90.cangas_sub(c8(press_vec[i]), \
                              c8(o2_ppress), \
                              c8(leaf_tempk), \
                              byref(mm_kco2_f), \
                              byref(mm_ko2_f), \
                              byref(co2_cpoint_f))

        iret = f90.lmr_ryan_sub(c8(lnc_top),c8(nscaler_top), ci(1), \
                                c8(leaf_tempk), byref(lmr_f))

        iret = f90.biophysrate_sub(ci(1), \
                                   c8(vcmax25_top), c8(jmax25_top), c8(kp25_top), \
                                   c8(nscaler_top), c8(leaf_tempk), c8(dayl_factor_full), \
                                   c8(t_growth_kum), c8(t_home_kum), c8(btran), \
                                   byref(vcmax_f), byref(jmax_f), byref(kp_f), byref(gs0_f), byref(gs1_f), byref(gs2_f))
        
        shared_inputs[ip,0] = par_abs_umol[i]
        shared_inputs[ip,1] = vcmax_f.value
        shared_inputs[ip,2] = jmax_f.value
        shared_inputs[ip,3] = gs2_f.value
        shared_inputs[ip,4] = leaf_tempk
        shared_inputs[ip,5] = can_press
        shared_inputs[ip,6] = co2_ppress          #
        shared_inputs[ip,7] = veg_es_f.value
        shared_inputs[ip,8] = gb
        shared_inputs[ip,9] = vpress
        shared_inputs[ip,10] = mm_kco2_f.value
        shared_inputs[ip,11] = mm_ko2_f.value
        shared_inputs[ip,12] = co2_cpoint_f.value
        shared_inputs[ip,13] = lmr_f.value

        
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

    n_trainset = 50
    n_infeatures = 14
    shared_inputs = torch.zeros([n_trainset,n_infeatures],dtype=torch.float32).share_memory_()

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
        p = mp.Process(target=RankPrepInput, args=(rank, chunk, shared_inputs))
        p.start()
        processes.append(p)
    
    # Wait for all processes to complete their work.
    for p in processes:
        p.join()
    
    # The shared tensor is now fully populated.
    print(f"Final shared tensor after all processes have written: {shared_inputs[:,0]}")

    # We invoke with torchrun because we use the distributed methods
    # to parallelize the training process. However, we use multiprocessing
    # for steps prior. So we use the former to identify the world size for
    # the latter
    ##dist.init_process_group(backend='gloo', init_method='env://')
    ##world_size = dist.get_world_size()
    ##dist.destroy_process_group()

    

    

