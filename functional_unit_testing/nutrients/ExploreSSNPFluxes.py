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
iret = f90.init_derived_sub()

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
elong_stem = 1.0
no_damage = 1
age0 = 0   # index for age-class 1

allom_nc = params['parameters']['fates_stoich_nitr']['data'] #fraction - ["fates_plant_organs", "fates_pft']
retrans_n = params['parameters']['fates_cnp_turnover_nitr_retrans']['data']  # fraction ["fates_plant_organs", "fates_pft']
turnover_fnrt = params['parameters']['fates_turnover_fnrt']['data'] # years ["fates_pft"]
turnover_branch = params['parameters']['fates_turnover_branch']['data']
turnover_leaf_canopy = params['parameters']['fates_turnover_leaf_canopy']['data']
turnover_leaf_ustory = params['parameters']['fates_turnover_leaf_ustory']['data']
nitr_store_ratio = params['parameters']['fates_cnp_nitr_store_ratio']['data']

        
#      "dtype": "float",
#      "dims": ["fates_leafage_class", "fates_pft"],
#"fates_turnover_leaf_ustory": {
#      "dtype": "float",
#      "dims": ["fates_leafage_class", "fates_pft"],
    
# fortran return values
cd_bfr = c_double(-9.0)
cd_dbfrdd = c_double(-9.0)
cd_bleaf = c_double(-9.0)
cd_dbleafdd = c_double(-9.0)
cd_asap = c_double(-9.0)
cd_bsap = c_double(-9.0)
cd_dbsapdd = c_double(-9.0)
cd_bdead = c_double(-9.0)
cd_dbdeaddd = c_double(-9.0)
cd_bstore = c_double(-9.0)
cd_dbstoredd = c_double(-9.0)
cd_bagw = c_double(-9.0)
cd_dbagwdd = c_double(-9.0)
cd_bbgw = c_double(-9.0)
cd_dbbgwdd = c_double(-9.0)

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
                        byref(cd_bleaf),byref(cd_dbleafdd))
    cleaf[i] = cd_bleaf.value
    nleaf[i] = cleaf[i]*allom_nc[leafid][pft0]
    nleaf_flux[i] = nleaf[i]*(1./turnover_leaf_canopy[age0][pft0])*(1.-retrans_n[leafid][pft0])

    # Sapwood
    iret = f90.bsap_allom_sub(c8(d),ci(pft),c8(no_damage),c8(no_trim),c8(elong_stem), \
                              byref(cd_asap),byref(cd_bsap),byref(cd_dbsapdd))
    csap[i] = cd_bsap.value
    nsap[i] = csap[i]*allom_nc[sapwid][pft0]
    nsap_flux[i] = nsap[i]*(1./turnover_branch[pft0])*(1.-retrans_n[sapwid][pft0])
    
    # Dead wood
    iret = f90.bagw_allom_sub(c8(d),ci(pft),c8(no_damage),c8(elong_stem),byref(cd_bagw),byref(cd_dbagwdd))
    iret = f90.bbgw_allom_sub(c8(d),ci(pft),c8(elong_stem),byref(cd_bbgw),byref(cd_dbbgwdd))
    iret = f90.bdead_allom_sub(c8(cd_bagw.value),c8(cd_bbgw.value),c8(cd_bsap.value),ci(pft), \
                               byref(cd_bdead),c8(cd_dbagwdd.value),c8(cd_dbbgwdd.value),c8(cd_dbsapdd.value),\
                               byref(cd_dbdeaddd))
    cdead[i] = cd_bdead.value
    ndead[i] = cdead[i]*allom_nc[deadid][pft0]
    ndead_flux[i] = ndead[i]*(1./turnover_branch[pft0])*(1.-retrans_n[deadid][pft0])

    # Storage
    iret = f90.bstore_allom_sub(c8(d),ci(pft),c8(no_damage),c8(no_trim), \
                                byref(cd_bstore),byref(cd_dbstoredd))
    cstore[i] = cd_bstore.value
    # Default hypothesis is the storage N is a factor of leaf N
    # Also, storage does not retranslocate on branchfall
    nstore[i] = nleaf[i] * nitr_store_ratio[pft0]
    nstore_flux[i] = nstore[i]*(1./turnover_branch[pft0])


    
# fates_stoich_nitr
# "dtype": "float",
# "dims": ["fates_plant_organs", "fates_pft"],
# "long_name": "target nitrogen concentration (ratio with carbon) of organs",
# "units": "gN/gC",
# "data": ["leaf", "fine root", "sapwood", "structure"]
#nfr_flux = np.zeros(np.shape(dbh))
#nleaf_flux = np.zeros(np.shape(dbh))
#nsap_flux  = np.zeros(np.shape(dbh))
#ndead_flux = np.zeros(np.shape(dbh))
#nstore_flux = np.zeros(np.shape(dbh))

fig,axs = plt.subplots(ncols=1,nrows=1,figsize=(6,5))
plt.stackplot(dbh, nleaf_flux,nfr_flux,nsap_flux,ndead_flux,nstore_flux,
              labels=['leaf', 'fine-root', 'sapwood', 'dead','storage'],
              alpha=0.8)
plt.legend(loc='upper left')
plt.grid(True, linestyle='--', alpha=0.6)


fig,axs = plt.subplots(ncols=1,nrows=1,figsize=(6,5))
plt.stackplot(dbh, np.log(nleaf_flux),np.log(nfr_flux),np.log(nsap_flux),np.log(ndead_flux),np.log(nstore_flux),
              labels=['leaf', 'fine-root', 'sapwood', 'dead','storage'],
              alpha=0.8)
plt.legend(loc='upper left')
plt.grid(True, linestyle='--', alpha=0.6)


fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(6,5))
ax.plot(dbh,(nleaf_flux+nfr_flux+nsap_flux+ndead_flux+nstore_flux)/cfr)
ax.set_xlabel('Stem Diameter [cm]')
ax.set_ylabel('Steady State N Uptake/Cfr [gN/yr/gC]')
ax.grid(True)


fig, axs = plt.subplots(ncols=2,nrows=2,figsize=(8,7))
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

ax=ax1s[2]
ap = ax.plot(dbh,nfr_flux)
ax.set_ylabel('Fineroot Nitrogen Flux [kgN/plant/year]')
ax.set_xlabel('Stem Diameter [cm]')
ax.grid(True)

ax=ax1s[3]
ap = ax.plot(dbh,nfr)
ax.set_ylabel('Fineroot Nitrogen [kgN]')
ax.set_xlabel('Stem Diameter [cm]')
ax.grid(True)
plt.show()

# Lets look at some of the soil dynamics
# ------------------------------------------------
# smin_nh4: minearlized nh4 in soil [g m-3 soil]
# Bulk density of dry soil material [kg m-3 soil]
# real(r8), parameter :: adsorp_nh4_eff = 2.76_r8
# real(r8), parameter :: m3_per_liter = 1.e-3_r8   ! m3 per liter
# Vol. Soil Water in each layer [m3 H2O/m3 soil]
# solution_conc: concentration of available mineralized nutrient
#                [gN / m3 H2O]


solution_conc = smin_nh4 / (bd(j)*adsorp_nh4_eff*m3_per_liter + h2osoi_vol(j))

compet_decomp = solution_conc / (km_decomp_nh4 * (1. + solution_conc/km_decomp_nh4 + e_km))

compet_nit    = solution_conc / (km_nit * (1. + solution_conc/km_nit + e_km))

compet_plant  = solution_conc / ( km_nh4_plant(ft) * (1._r8 + solution_conc/km_nh4_plant(ft) + e_km))

sum_nh4_demand_scaled = sum_plant_nh4demand + potential_immob*compet_decomp + pot_f_nit*compet_nit

actual_immob_nh4 = min(potential_immob,(smin_nh4/dt)* \
                       (potential_immob*compet_decomp / sum_nh4_demand_scaled))
