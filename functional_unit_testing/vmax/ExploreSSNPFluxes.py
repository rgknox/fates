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

dbh = np.linspace(0.5,200.0,1000)

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
nfr_ss_flux = np.zeros(np.shape(dbh))
nleaf_ss_flux = np.zeros(np.shape(dbh))
nsap_ss_flux  = np.zeros(np.shape(dbh))
ndead_ss_flux = np.zeros(np.shape(dbh))
nstore_ss_flux = np.zeros(np.shape(dbh))

nfr_gr1_flux = np.zeros(np.shape(dbh))
nleaf_gr1_flux = np.zeros(np.shape(dbh))
nsap_gr1_flux  = np.zeros(np.shape(dbh))
ndead_gr1_flux = np.zeros(np.shape(dbh))
nstore_gr1_flux = np.zeros(np.shape(dbh))
repro_gr1_flux =  np.zeros(np.shape(dbh))
nrep_gr1_flux =  np.zeros(np.shape(dbh))

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

secperyear=86400.*365.
g_per_kg = 1000.


allom_nc = params['parameters']['fates_stoich_nitr']['data'] #fraction - ["fates_plant_organs", "fates_pft']
retrans_n = params['parameters']['fates_cnp_turnover_nitr_retrans']['data']  # fraction ["fates_plant_organs", "fates_pft']
turnover_fnrt = params['parameters']['fates_turnover_fnrt']['data'] # years ["fates_pft"]
turnover_branch = params['parameters']['fates_turnover_branch']['data']
turnover_leaf_canopy = params['parameters']['fates_turnover_leaf_canopy']['data']
turnover_leaf_ustory = params['parameters']['fates_turnover_leaf_ustory']['data']
nitr_store_ratio = params['parameters']['fates_cnp_nitr_store_ratio']['data']
vmax_nh4 = params['parameters']['fates_cnp_vmax_nh4']['data']
vmax_no3 = params['parameters']['fates_cnp_vmax_no3']['data']

# fortran return values
cd_bfr = c_double(-9.0)
cd_dbfrdd = c_double(-9.0)
cd_bleaf = c_double(-9.0)
cd_dbleafdd = c_double(-9.0)
cd_asap = c_double(-9.0)
cd_bsap = c_double(-9.0)
cd_dbsapwdd = c_double(-9.0)
cd_bdead = c_double(-9.0)
cd_dbdeaddd = c_double(-9.0)
cd_bstore = c_double(-9.0)
cd_dbstoredd = c_double(-9.0)
cd_bagw = c_double(-9.0)
cd_dbagwdd = c_double(-9.0)
cd_bbgw = c_double(-9.0)
cd_dbbgwdd = c_double(-9.0)

cm_per_yr1 = 1.0

# Calculate the N:C ratio of a recruit
# We use that to determine the N:C ratio of seed
# ---------------------------------------------------------------------------------------

recr_hgt = params['parameters']['fates_recruit_height_min']['data']
recr_seed_alloc = params['parameters']['fates_recruit_seed_alloc']['data']
recr_seed_alloc_mature = params['parameters']['fates_recruit_seed_alloc_mature']['data']
recr_dbh_mature = params['parameters']['fates_recruit_seed_dbh_repro_threshold']['data']
cd_dmin = c_double(-9.0)
cd_ddmindh = c_double(-9)
iret = f90.h2d_allom_sub(c8(recr_hgt[pft0]),ci(pft),byref(cd_dmin),byref(cd_ddmindh))
drec = cd_dmin.value
iret = f90.bfineroot_sub(c8(drec),ci(pft),c8(no_trim),c8(l2fr),c8(elong_fnrt),byref(cd_bfr),byref(cd_dbfrdd))
crec = cd_bfr.value
nrec = cd_bfr.value*allom_nc[fnrtid][pft0]
iret = f90.bleaf_sub(c8(drec),ci(pft),ci(no_damage),c8(no_trim),c8(elong_leaf),byref(cd_bleaf),byref(cd_dbleafdd))
crec = crec + cd_bleaf.value
nrec = nrec + cd_bleaf.value*allom_nc[leafid][pft0]
iret = f90.bsap_allom_sub(c8(drec),ci(pft),c8(no_damage),c8(no_trim),c8(elong_stem),byref(cd_asap),byref(cd_bsap),byref(cd_dbsapwdd))
crec = crec + cd_bsap.value
nrec = nrec + cd_bsap.value*allom_nc[sapwid][pft0]
iret = f90.bagw_allom_sub(c8(drec),ci(pft),c8(no_damage),c8(elong_stem),byref(cd_bagw),byref(cd_dbagwdd))
iret = f90.bbgw_allom_sub(c8(drec),ci(pft),c8(elong_stem),byref(cd_bbgw),byref(cd_dbbgwdd))
iret = f90.bdead_allom_sub(c8(cd_bagw.value),c8(cd_bbgw.value),c8(cd_bsap.value),ci(pft), \
                               byref(cd_bdead),c8(cd_dbagwdd.value),c8(cd_dbbgwdd.value),c8(cd_dbsapwdd.value),\
                               byref(cd_dbdeaddd))
crec = crec + cd_bdead.value
nrec = nrec + cd_bdead.value*allom_nc[deadid][pft0]
iret = f90.bstore_allom_sub(c8(drec),ci(pft),c8(no_damage),c8(no_trim),byref(cd_bstore),byref(cd_dbstoredd))
crec = crec + cd_bstore.value
nrec = nrec + cd_bleaf.value*allom_nc[leafid][pft0]*nitr_store_ratio[pft0]
nc_recr = nrec/crec
# ------------------------------------------------------------------------------------


for i, d in enumerate(dbh):

    # All units in kg/plant and kg/plant/yr
    # Fine-root
    iret = f90.bfineroot_sub(c8(d),ci(pft),c8(no_trim), \
                             c8(l2fr),c8(elong_fnrt),byref(cd_bfr),byref(cd_dbfrdd))
    cfr[i] = cd_bfr.value
    nfr[i] = cfr[i] * allom_nc[fnrtid][pft0]
    nfr_ss_flux[i] = nfr[i]*(1./turnover_fnrt[pft0])*(1.-retrans_n[fnrtid][pft0])*g_per_kg
    nfr_gr1_flux[i] = cd_dbfrdd.value*cm_per_yr1*allom_nc[fnrtid][pft0]*g_per_kg
    
    # Leaf
    iret = f90.bleaf_sub(c8(d),ci(pft),ci(no_damage),c8(no_trim),c8(elong_leaf), \
                        byref(cd_bleaf),byref(cd_dbleafdd))
    cleaf[i] = cd_bleaf.value
    nleaf[i] = cleaf[i]*allom_nc[leafid][pft0]
    nleaf_ss_flux[i] = nleaf[i]*(1./turnover_leaf_canopy[age0][pft0])*(1.-retrans_n[leafid][pft0])*g_per_kg
    nleaf_gr1_flux[i] = cd_dbleafdd.value*cm_per_yr1*allom_nc[leafid][pft0]*g_per_kg

    # Sapwood
    iret = f90.bsap_allom_sub(c8(d),ci(pft),c8(no_damage),c8(no_trim),c8(elong_stem), \
                              byref(cd_asap),byref(cd_bsap),byref(cd_dbsapwdd))
    csap[i] = cd_bsap.value
    nsap[i] = csap[i]*allom_nc[sapwid][pft0]
    nsap_ss_flux[i] = nsap[i]*(1./turnover_branch[pft0])*(1.-retrans_n[sapwid][pft0])*g_per_kg
    nsap_gr1_flux[i] = cd_dbsapwdd.value*cm_per_yr1*allom_nc[sapwid][pft0]*g_per_kg

    # Dead wood
    iret = f90.bagw_allom_sub(c8(d),ci(pft),c8(no_damage),c8(elong_stem),byref(cd_bagw),byref(cd_dbagwdd))
    iret = f90.bbgw_allom_sub(c8(d),ci(pft),c8(elong_stem),byref(cd_bbgw),byref(cd_dbbgwdd))
    iret = f90.bdead_allom_sub(c8(cd_bagw.value),c8(cd_bbgw.value),c8(cd_bsap.value),ci(pft), \
                               byref(cd_bdead),c8(cd_dbagwdd.value),c8(cd_dbbgwdd.value),c8(cd_dbsapwdd.value),\
                               byref(cd_dbdeaddd))
    cdead[i] = cd_bdead.value
    ndead[i] = cdead[i]*allom_nc[deadid][pft0]
    ndead_ss_flux[i] = ndead[i]*(1./turnover_branch[pft0])*(1.-retrans_n[deadid][pft0])*g_per_kg
    ndead_gr1_flux[i] = cd_dbdeaddd.value*cm_per_yr1*allom_nc[deadid][pft0]*g_per_kg
    
    # Storage
    iret = f90.bstore_allom_sub(c8(d),ci(pft),c8(no_damage),c8(no_trim), \
                                byref(cd_bstore),byref(cd_dbstoredd))
    cstore[i] = cd_bstore.value
    # Default hypothesis is the storage N is a factor of leaf N
    # Also, storage does not retranslocate on branchfall
    nstore[i] = nleaf[i] * nitr_store_ratio[pft0]
    nstore_ss_flux[i] = nstore[i]*(1./turnover_branch[pft0])*g_per_kg
    nstore_gr1_flux[i] = nleaf_gr1_flux[i] * nitr_store_ratio[pft0]

    # Recruitment
    if(d>recr_dbh_mature[pft0]):
        recr_frac = recr_seed_alloc[pft0] + recr_seed_alloc_mature[pft0]
    else:
        recr_frac = recr_seed_alloc[pft0]

    gr1_dbdd = cd_dbstoredd.value+cd_dbdeaddd.value+cd_dbsapwdd.value+cd_dbleafdd.value+cd_dbfrdd.value
    nrep_gr1_flux[i] = gr1_dbdd*cm_per_yr1*nc_recr*recr_frac*g_per_kg
    
    
# fates_stoich_nitr
# "dtype": "float",
# "dims": ["fates_plant_organs", "fates_pft"],
# "long_name": "target nitrogen concentration (ratio with carbon) of organs",
# "units": "gN/gC",
# "data": ["leaf", "fine root", "sapwood", "structure"]
#nfr_ss_flux = np.zeros(np.shape(dbh))
#nleaf_ss_flux = np.zeros(np.shape(dbh))
#nsap_ss_flux  = np.zeros(np.shape(dbh))
#ndead_ss_flux = np.zeros(np.shape(dbh))
#nstore_ss_flux = np.zeros(np.shape(dbh))

colors = ['green','orange','blue','purple','cyan','green','orange','blue','purple','cyan',[0.3,0.3,0.5]]
patterns = ['','','','','','///','///','///','///','///','...']

fig,axs = plt.subplots(ncols=1,nrows=1,figsize=(6,5))
stacks = plt.stackplot(dbh, nleaf_ss_flux,nfr_ss_flux,nsap_ss_flux,ndead_ss_flux,nstore_ss_flux,nleaf_gr1_flux,nfr_gr1_flux,nsap_gr1_flux,ndead_gr1_flux,nstore_gr1_flux,nrep_gr1_flux,
                       labels=['ss-leaf', 'ss-fnrt', 'ss-sapw', 'ss-dead','ss-store','gr-leaf', 'gr-fnrt', 'gr-sapw', 'gr-dead','gr-store','recruit'],
                       colors=colors,
                       alpha=0.8)
plt.legend(loc='upper left')
plt.xlabel('Stem Diameter [cm]')
plt.ylabel('Steady State N Demand [g/yr]')
plt.grid(True, linestyle='--', alpha=0.6)

for stack, pattern in zip(stacks, patterns):
    stack.set_hatch(pattern)

ss_flux = nleaf_ss_flux+nfr_ss_flux+nsap_ss_flux+ndead_ss_flux+nstore_ss_flux
gr1_flux = nleaf_gr1_flux+nfr_gr1_flux+nsap_gr1_flux+ndead_gr1_flux+nstore_gr1_flux+nrep_gr1_flux

fig,ax = plt.subplots(ncols=1,nrows=1,figsize=(6,5))
ax.plot(dbh,ss_flux/(cfr*g_per_kg))
ax.plot(dbh,(ss_flux+gr1_flux)/(cfr*g_per_kg))
ax.set_xlabel('Stem Diameter [cm]')
ax.set_ylabel('1cm Growth Uptake/Cfr [gN/yr/gC]')
ax.axhline((vmax_nh4[pft0]+vmax_no3[pft0])*secperyear,label='vmax',linestyle='--',color='k')
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
ap = ax.plot(dbh,nfr_ss_flux)
ax.set_ylabel('Fineroot Nitrogen Flux [kgN/plant/year]')
ax.set_xlabel('Stem Diameter [cm]')
ax.grid(True)

ax=ax1s[3]
ap = ax.plot(dbh,nfr)
ax.set_ylabel('Fineroot Nitrogen [kgN]')
ax.set_xlabel('Stem Diameter [cm]')
ax.grid(True)
plt.show()

