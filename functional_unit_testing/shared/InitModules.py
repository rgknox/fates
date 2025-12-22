import numpy as np
import matplotlib.pyplot as plt
import json
import sys
sys.path.insert(0,'py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb
from CtypesInit import f90_modules # Directly import f90_modules from CtypesInit
from ctypes import *
import code  # For development: code.interact(local=dict(globals(), **locals()))


f90 = f90_modules('bld/')


parameterfile = '../../parameter_files/fates_params_default.json'

with open(parameterfile, 'r') as file:
    params = json.load(file)
    
#code.interact(local=dict(globals(), **locals()))
        
iret = f90.json_setinval(c8(9.9e32+10))
iret = f90.json_setloginit(ci(6))
iret = f90.wrapjson_read(*ccharnb(parameterfile.strip()))

iret = f90.getparams()
iret = f90.alloc_allomparam_sub(ci(14))

dbh = np.linspace(0.2,200.0,1000)

# Carbon biomass 
cfr = np.zeros(np.shape(dbh))
cleaf = np.zeros(np.shape(dbh))
csap  = np.zeros(np.shape(dbh))
cdead = np.zeros(np.shape(dbh))
cstore = np.zeros(np.shape(dbh))

# Nitrogen biomass
nfr = np.zeros(np.shape(dbh))
nleaf = np.zeros(np.shape(dbh))
nsap  = np.zeros(np.shape(dbh))
ndead = np.zeros(np.shape(dbh))
nstore = np.zeros(np.shape(dbh))

allom_nc = params['parameters']['fates_stoich_nitr']['data']

leafid=0
fnrtid=1
sapwid=2
deadid=3

pft = 1
pft0 = pft - 1
no_trim = 1.0
l2fr = 1.0
elong_fnrt = 1.0

# fortran return values
cd_bfr = c_double(-9.0)
cd_dbfrdd = c_double(-9.0)

for i, d in enumerate(dbh):
    iret = f90.bfineroot_sub(c8(d),ci(pft),c8(no_trim),c8(l2fr),c8(elong_fnrt),byref(cd_bfr),byref(cd_dbfrdd))
    cfr[i] = cd_bfr.value
    nfr[i] = cfr[i] * allom_nc[fnrtid][pft0]
    

# fates_stoich_nitr
# "dtype": "float",
# "dims": ["fates_plant_organs", "fates_pft"],
# "long_name": "target nitrogen concentration (ratio with carbon) of organs",
# "units": "gN/gC",
# "data": ["leaf", "fine root", "sapwood", "structure"]

                                                     
fig, axs = plt.subplots(ncols=2,nrows=2,figsize=(6,3))
ax1s = axs.reshape(-1)
ax=ax1s[0]
ap = ax.plot(dbh,cfr)
ax.set_ylabel('Fineroot Carbon [kgC]')
ax.set_xlabel('Stem Diameter [cm]')
ax.grid(True)

ax=ax1s[1]
ap = ax.plot(dbh,nfr)
ax.set_ylabel('Fineroot Nitrogen [kgN]')
ax.set_xlabel('Stem Diameter [cm]')
ax.grid(True)
plt.show()
