import numpy as np
import matplotlib.pyplot as plt
import json
import sys
sys.path.insert(0,'../shared/py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb
from CtypesInit import f90_modules # Directly import f90_modules from CtypesInit
from ctypes import *
import code  # For development: code.interact(local=dict(globals(), **locals()))


f90 = f90_modules('../shared/bld/')


parameterfile = '../../parameter_files/fates_params_default.json'

with open(parameterfile, 'r') as file:
    params = json.load(file)
    
#code.interact(local=dict(globals(), **locals()))
        
iret = f90.json_setinval(c8(9.9e32+10))
iret = f90.json_setloginit(ci(6))
iret = f90.wrapjson_read(*ccharnb(parameterfile.strip()))

iret = f90.getparams()
#iret = f90.alloc_allomparam_sub(ci(14))

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

# Nitrogen Flux
nfr_flux = np.zeros(np.shape(dbh))
nleaf_flux = np.zeros(np.shape(dbh))
nsap_flux  = np.zeros(np.shape(dbh))
ndead_flux = np.zeros(np.shape(dbh))
nstore_flux = np.zeros(np.shape(dbh))


leafid=0
fnrtid=1
sapwid=2
deadid=3
pft = 1
pft0 = pft - 1
no_trim = 1.0
l2fr = 1.0
elong_fnrt = 1.0
elong_leaf = 1.0
no_damage = 1
age0 = 0   # index for age-class 1

allom_nc = params['parameters']['fates_stoich_nitr']['data'] #fraction - ["fates_plant_organs", "fates_pft']
retrans_n = params['parameters']['fates_cnp_turnover_nitr_retrans']['data']  # fraction ["fates_plant_organs", "fates_pft']
turnover_fnrt = params['parameters']['fates_turnover_fnrt']['data'] # years ["fates_pft"]
turnover_branch = params['parameters']['fates_turnover_branch']['data']
turnover_leaf_canopy = params['parameters']['fates_turnover_leaf_canopy']['data']
turnover_leaf_ustory = params['parameters']['fates_turnover_leaf_ustory']['data']

#      "dtype": "float",
#      "dims": ["fates_leafage_class", "fates_pft"],
#"fates_turnover_leaf_ustory": {
#      "dtype": "float",
#      "dims": ["fates_leafage_class", "fates_pft"],
    
# fortran return values
cd_bfr = c_double(-9.0)
cd_dbfrdd = c_double(-9.0)
cd_leaf = c_double(-9.0)
cd_dleafdd = c_double(-9.0)
cd_asap = c_double(-9.0)
cd_sap = c_double(-9.0)
cd_dsapdd = c_double(-9.0)

for i, d in enumerate(dbh):

    # All units in kg/plant and kg/plant/yr
    # Fine-root
    iret = f90.bfineroot_sub(c8(d),ci(pft),c8(no_trim), \
                             c8(l2fr),c8(elong_fnrt),byref(cd_bfr),byref(cd_dbfrdd))
    cfr[i] = cd_bfr.value
    nfr[i] = cfr[i] * allom_nc[fnrtid][pft0]
    nfr_flux[i] = nfr[i]*(1./turnover_fnrt[pft0])*(1.-retrans_n[fnrtid][pft0])

    # Leaf
    iret = f90.bleaf_sub(c8(d),ci(pft),ci(no_damage),c8(no_trim),c8(elong_leaf), \
                        byref(cd_leaf),byref(cd_dleafdd))
    cleaf[i] = cd_leaf.value
    nleaf[i] = cleaf[i]*allom_nc[leafid][pft0]
    nleaf_flux[i] = nleaf[i]*(1./turnover_leaf_canopy[age0][pft0])*(1.-retrans_n[leafid][pft0])

    # Sapwood
    #iret = f90.bsap_allom_sub(cd(d),ci(pft),c8(no_damage),c8(no_trim),c8(elong_stem), \
    #                          byref(asap),byref(cd_sap),byref(cd_dsapdd))

    # Dead wood
    #iret = f90.bdead_allom_sub()
    
    

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

ax=ax1s[2]
ap = ax.plot(dbh,nfr_flux)
ax.set_ylabel('Fineroot Nitrogen Flux [kgN/plant/year]')
ax.set_xlabel('Stem Diameter [cm]')
ax.grid(True)
plt.show()

ax=ax1s[3]
ap = ax.plot(dbh,nfr)
ax.set_ylabel('Fineroot Nitrogen [kgN]')
ax.set_xlabel('Stem Diameter [cm]')
ax.grid(True)
plt.show()

