import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime
import argparse
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


g_per_kg = 1000.0
kg_per_g = 0.001

def DecayCoeffVcmax(vcmax25top,p_slope,p_int):

    kn = np.exp(p_slope * vcmax25top - p_int)

    return kn

def GetSLA(slatop,lai_0,kn,lai):

    sla = slatop / np.exp(-kn*(lai_0+lai))

    return sla

def GetLeafCSLAMax(slatop,slamax,lai_0,kn):

    leafc_slamax = (slatop - slamax * np.exp(-kn * lai_0)) / (-kn * slatop * slamax)

    return leafc_slamax


def LaiFromLeafC(slatop,lai_0,kn,leafc):
        
    lai = (np.log(np.exp( -kn * lai_0) - kn * slatop * leafc ) + (kn * lai_0)) / (-kn)

    return lai


def LeafCFromLai(slatop,lai_0,kn,lai):

    leafc = (np.exp( -kn * lai_0) - np.exp(-kn*(lai+lai_0)))/(kn * slatop)
    
    return leafc

# fates_leafn_vert_scaler_coeff1:long_name = "Coefficient one for decrease in leaf nitrogen through the canopy, from Lloyd et al. 2010." ;
# fates_leafn_vert_scaler_coeff2:long_name = "Coefficient two for decrease in leaf nitrogen through the canopy, from Lloyd et al. 2010." ;
# fates_leaf_vcmax25top:long_name = "maximum carboxylation rate of Rub. at 25C, canopy top" ;
# fates_leaf_slatop:long_name = "Specific Leaf Area (SLA) at top of canopy, projected area basis" "m^2/gC" ;
# fates_leaf_slamax:long_name = "Maximum Specific Leaf Area (SLA), even if under a dense canopy"   "m^2/gC" ;


fates_leafn_vert_scaler_coeff1 = [0.00963, 0.00963, 0.00963, 0.00963, 0.00963, 0.00963, 0.00963, \
                                  0.00963, 0.00963, 0.00963, 0.00963, 0.00963, 0.00963, 0.00963]

fates_leafn_vert_scaler_coeff2 = [2.43, 2.43, 2.43, 2.43, 2.43, 2.43, 2.43, \
                                  2.43, 2.43, 2.43, 2.43, 2.43, 2.43, 2.43]

fates_leaf_vcmax25top = [50, 62, 39, 61, 58, 58, 62, 54, 54, 38, 54, 86, 78, 78]

fates_leaf_slamax = np.array([0.0954, 0.0954, 0.0954, 0.0954, 0.0954, 0.0954, 0.012, \
                     0.03, 0.03, 0.012, 0.032, 0.05, 0.05, 0.05]) * g_per_kg

fates_leaf_slatop = np.array([0.012, 0.005, 0.024, 0.009, 0.03, 0.03, 0.012, 0.03, \
                     0.03, 0.01, 0.032, 0.027, 0.05, 0.05]) * g_per_kg

kn = []

npft = len(fates_leafn_vert_scaler_coeff1)


# Preprocess the decay coefficient
for ft in range(npft):
    kn.append(DecayCoeffVcmax(fates_leaf_vcmax25top[ft],fates_leafn_vert_scaler_coeff1[ft],fates_leafn_vert_scaler_coeff2[ft]))
    print(kn[ft])


n_lai = 50
max_lai = 5.0

lai_vec = np.linspace(0.1,max_lai,n_lai)

lai_0 = 1.0

ft = 0

bleaf = max_lai/fates_leaf_slatop[ft]


dbleaf = 0.5*bleaf/n_lai

sla_vec = np.zeros(lai_vec.shape)
sla_vec2 = np.zeros(lai_vec.shape)
int_lai = np.zeros(lai_vec.shape)
an_lai  = np.zeros(lai_vec.shape)
leafc   = np.zeros(lai_vec.shape)
leafc0  = np.zeros(lai_vec.shape)

#leafc = integral of e^(-kn(x + canopy_lai_above)) / slatop

for il,lai in enumerate(lai_vec):

    sla_vec[il] = GetSLA(fates_leaf_slatop[ft],lai_0,kn[ft],lai)
    
    if(il>0):
        sla_vec2[il] = GetSLA(fates_leaf_slatop[ft],lai_0,kn[ft],int_lai[il-1])
        int_lai[il] = int_lai[il-1] + dbleaf*sla_vec2[il]
        an_lai[il] = LaiFromLeafC(fates_leaf_slatop[ft],lai_0,kn[ft], float(il)*dbleaf)
        
    else:
        sla_vec2[il]=sla_vec[il]


for il,lai in enumerate(lai_vec):

    leafc0[il] = float(il)*dbleaf
    lai = LaiFromLeafC(fates_leaf_slatop[ft],lai_0,kn[ft], leafc0[il])
    leafc[il]  = LeafCFromLai(fates_leaf_slatop[ft],lai_0,kn[ft],lai)


    
        
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(8.5,7.5))
        
ax1.plot(sla_vec,lai_vec)
ax1.plot(sla_vec2,int_lai)
ax1.plot(sla_vec2,an_lai)
ax1.set_ylabel('LAI [m2/m2]')
ax1.set_xlabel('SLA [m2/kg]')
ax1.invert_yaxis()
ax1.grid(True)


ax2.scatter(leafc0,leafc)
ax2.set_ylabel('Leaf C (Back Calculated)')
ax2.set_xlabel('Leaf C (coordinate)')
ax2.grid(True)

plt.show()
