import numpy as np
import matplotlib
import os
import sys
import getopt
import code  # For development: code.interact(local=dict(globals(), **locals()))
from PushParameters import GetParamFromAttrib
from PushParameters import GetParamList
from ctypes import *
sys.path.append('../shared/py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

# Diagnostics that are grouped by radiation scattering element
class elem_diags_type:

    def __init__(self,vai_above,lai_above,total_lai,total_sai,dvai,crown_area_frac,pft,ntimes):

        # Static diagnostics
        self.pft       = pft
        self.vai_above = vai_above   # Maximum depth integrated VAI above this element 
        self.lai_above = lai_above   # Mean LAI in canopy layers above this one (/m2 ground)
        self.total_lai = total_lai   # LAI of this element (/m2 crown)
        self.total_sai = total_sai   # SAI of this element (/m2 crown)
        self.total_vai = total_lai+total_sai
        self.leaf_frac = total_lai/(total_lai+total_sai)
        self.n_layer = int(np.ceil((total_lai+total_sai)/dvai))  # Number of layers we discretize
        self.dlai = dvai*self.leaf_frac
        self.crown_area_frac = crown_area_frac

        self.vai_top = np.zeros(self.n_layer)
        self.vai_bot = np.zeros(self.n_layer)
        self.canopy_vai_mid = np.zeros(self.n_layer)
        self.canopy_lai_mid = np.zeros(self.n_layer)
        self.sil     = np.zeros(self.n_layer,dtype=int)
        for il in range(self.n_layer):
            self.vai_top[il] = dvai*float(il)
            if(il == (self.n_layer-1)):
                self.vai_bot[il] =  self.total_vai
            else:
                self.vai_bot[il] =  self.vai_top[il]+dvai
        
        # Instantaneous diagnostics
        self.rd_abs_leaf = np.zeros(self.n_layer)
        self.rb_abs_leaf = np.zeros(self.n_layer)
        self.r_abs_stem  = np.zeros(self.n_layer)
        self.rd_dn = np.zeros(self.n_layer)
        self.rd_up = np.zeros(self.n_layer)
        self.rbeam =  np.zeros(self.n_layer)
        self.sunfrac = np.zeros(self.n_layer)
        self.sunfrac_v2 = np.zeros(self.n_layer)
        self.lmr  = np.zeros(self.n_layer)
        self.agross = np.zeros(self.n_layer)
        self.anet = np.zeros(self.n_layer)
        self.gstoma = np.zeros(self.n_layer)
        self.co2_interc = np.zeros(self.n_layer)
        self.vcmax = np.zeros(self.n_layer)
        self.jmax = np.zeros(self.n_layer)
        self.kp = np.zeros(self.n_layer)


        # Mean diagnostics
        self.ag_limit = np.zeros([self.n_layer,2])
        self.ag_sslimit = np.zeros([self.n_layer,2,2])
        self.n_zen          = 5          # Number of zenith bins we use for diagnostics
        self.zen_bins       = np.linspace(0,1.-1/self.n_zen,self.n_zen)
        self.sunfrac_zen_ll = np.zeros([self.n_layer,self.n_zen])
        self.sunfrac2_zen_ll= np.zeros([self.n_layer,self.n_zen])
        self.count_zen_ll   = np.zeros([self.n_layer,self.n_zen])
        return
    
    def ZeroDiag(self):

        # Zero out the instantaneous diagnostics that
        # are not indexed by time.
        self.rd_abs_leaf = np.zeros(self.n_layer)
        self.rb_abs_leaf = np.zeros(self.n_layer)
        self.r_abs_stem  = np.zeros(self.n_layer)
        self.rd_dn = np.zeros(self.n_layer)
        self.rd_up = np.zeros(self.n_layer)
        self.rbeam =  np.zeros(self.n_layer)
        self.sunfrac = np.zeros(self.n_layer)
        self.sunfrac_v2 = np.zeros(self.n_layer)
        self.lmr  = np.zeros(self.n_layer)
        self.agross = np.zeros(self.n_layer)
        self.anet = np.zeros(self.n_layer)
        self.gstoma = np.zeros(self.n_layer)
        self.co2_interc = np.zeros(self.n_layer)
        self.vcmax = np.zeros(self.n_layer)
        self.jmax = np.zeros(self.n_layer)
        self.kp = np.zeros(self.n_layer)
        return

class site_diags_type:

    def __init__(self,ntimes,dvai,total_vai):

        self.ntimes        = ntimes
        self.n_layer       = int(np.ceil(total_vai/dvai))
        self.vai_top       = np.zeros(self.n_layer) # top down integrated 'in-canopy"
                                                    # vegetation area index (stem+leaf)
                                                    # m2 vegetation / m2 crowns
                                                    # assessed at the top of the layer
        self.lai_ground    = np.zeros(self.n_layer) # total amount of leaf area
                                                    # in the present layer per ground area
                                                    # m2 leaf / m2 ground 
        for il in range(self.n_layer):
            self.vai_top[il] = dvai*float(il)
            
        self.aglimit_which = []
        self.aglimit_apar  = []
        self.aglimit_temp  = []
        self.lmr           = np.zeros(ntimes) # [umol/m2/s ground]
        self.agross        = np.zeros(ntimes) # [umol/m2/s ground]
        self.gstoma        = np.zeros(ntimes) # [umol/m2/s ground]
        self.anet          = np.zeros(ntimes) # [umol/m2/s ground]
        self.r_abs_leaf    = np.zeros(ntimes) # [W/m2 ground]
        self.g_b_umol      = np.zeros(ntimes) # [umol/m2/s]
        self.mm_kco2       = np.zeros(ntimes)
        self.mm_ko2        = np.zeros(ntimes)
        self.co2_cpoint    = np.zeros(ntimes) # [Pa CO2]
        self.gs0           = np.zeros(ntimes) # [0-1]
        self.gs1           = np.zeros(ntimes) # [0-1]
        self.gs2           = np.zeros(ntimes) # [0-1]
        self.solve_iter    = np.zeros(ntimes) # [count]
        self.co2_interc    = np.zeros(ntimes) # [Pa CO2]
       
        # Time x VAI depth diagnostics
        
        self.lmr_vl            = np.zeros([ntimes,self.n_layer]) # [umol/m2 leaf /s]
        self.agross_vl         = np.zeros([ntimes,self.n_layer]) # [umol/m2 leaf/s ]
        self.gstoma_vl         = np.zeros([ntimes,self.n_layer]) # [umol/m2 leaf/s ]
        self.anet_vl           = np.zeros([ntimes,self.n_layer]) # [umol/m2 leaf/s ]
        self.r_abs_leaf_vl     = np.zeros([ntimes,self.n_layer]) # [W/m2 leaf ]
        self.vcmax_vl          = np.zeros([ntimes,self.n_layer]) # [umol/m2 leaf/s ]
        self.agross_rubisco_vl = np.zeros([ntimes,self.n_layer]) # [umol/m2 leaf/s ]
        self.agross_rubpc3_vl  = np.zeros([ntimes,self.n_layer]) # [umol/m2 leaf/s ]
        self.sunfrac_vl        = np.zeros([ntimes,self.n_layer])

        
def GetElemLayerLAIShare(elem_diags,site_diags):

    for il in range(elem_diags.n_layer):

        vai_mid = 0.5*(elem_diags.vai_bot[il]+elem_diags.vai_top[il])
        
        # These area indices are WRT the canopy top
        # -----------------------------------------
        canopy_vai_mid = elem_diags.vai_above + vai_mid
        
        elem_diags.canopy_vai_mid[il] = canopy_vai_mid
        elem_diags.canopy_lai_mid[il] = elem_diags.lai_above + vai_mid*elem_diags.leaf_frac
        
        # Index for the whole-canopy vertical array (0 for first index)
        sil = np.sum(site_diags.vai_top < canopy_vai_mid)-1

        elem_diags.sil[il] = sil
        
        # Fraction of the site-layer's leaf area  occupied by this element layer's leaves
        site_diags.lai_ground[sil] = site_diags.lai_ground[sil] + elem_diags.crown_area_frac*elem_diags.dlai
        
    return elem_diags,site_diags



# Interpret XML parameters to create and element-based canopy structure
# This also sets up the diagnostic data types

def SetupCanopyDiags(xmlroot,ntimes,dvai,f90):
        
    # Create a canopy based on input from the control file
    # -----------------------------------------------------------------------------------
    
    cstruct_root = xmlroot.find('canopy_structure')
    cohort_pft  = GetParamList(cstruct_root,'cohort_pft','integer')
    cohort_area = GetParamList(cstruct_root,'cohort_area_frac','float')
    cohort_lai  = GetParamList(cstruct_root,'cohort_lai','float')
    cohort_sai  = GetParamList(cstruct_root,'cohort_sai','float')
    ncohorts    = len(cohort_pft)
    site_root   = xmlroot.find('site_structure')
    ground_vis_albedo = GetParamList(site_root,'ground_vis_albedo','float')
    ground_nir_albedo = GetParamList(site_root,'ground_nir_albedo','float')
    frac_snow = GetParamList(site_root,'frac_snow','float')[0]
    n_layer = GetParamList(site_root,'nlayers','integer')[0] 

    
    # Lets walk through the list and assign the layer of each
    cohort_can    = list(range(ncohorts)) # The canopy layer associated with each cohort (0 indexed)
    cohort_col    = list(range(ncohorts)) # Column index of the cohort                   (0 indexed)
    layer_ncohort = [0]                   # Number of cohorts (elements) in each layer
    area_current  = 0
    ican = 0
    icol = 0
    for ico in range(ncohorts):

        if( cohort_area[ico]-1.0 < 0.001 ):
            cohort_area[ico] = np.min([1.0,cohort_area[ico]])
        else:
            print('\n\nYou specified a cohort that takes up more area than 1.0\n')
            print('Dont do that. Exiting.\n')
            exit(1)

        if(cohort_area[ico]+area_current>1.0):
            layer_ncohort.append(0)
            ican = ican+1
            area_current  = 0.
            icol = 0
        else:
            area_current = area_current + cohort_area[ico]
            icol = icol + 1
                
        cohort_can[ico]     = ican
        layer_ncohort[ican] = layer_ncohort[ican]+1
        cohort_col[ico]     = icol-1

    # Determine the size of the canopy scattering matrix
    n_can  = len(layer_ncohort)
    n_col  = np.max(layer_ncohort)

    # Setup the cohort index matrix, initialize non-cohorts
    # as air, with pft -1
    elem_cohort = np.zeros([n_can,n_col],dtype=int)-1
    for ico in range(ncohorts):
        ican  = cohort_can[ico]
        icol  = cohort_col[ico]
        print("cohort ican: {} icol: {}".format(ican,icol))
        elem_cohort[ican,icol] = ico
    

    # Initialize scattering elements
    # -----------------------------------------------------------------------------

    elem_diags = []
    iret = f90.alloc_twostream_sub(ci(n_can),ci(n_col))
    iret = f90.param_prep_sub()  # This routine creates parameters that are derived from others
    lai_above = 0.               # Mean LAI (m2 ground) in canopy above current (for N decay coeffs)
    max_vai_above   = 0.               # Maximum VAI in crown (for assigning vertical position)
    for ican in range(n_can):
        veg_area = 0.0
        n_veg    = 0
        lai_current = 0    # Mean LAI (m2 crown) in current element
        vai_max_layer = 0.
        elem_diags.append(list())
        for icol in range(n_col):
            ico = elem_cohort[ican,icol]
            if(ico>=0):
                veg_area = veg_area + cohort_area[ico]
                lai_current = lai_current + cohort_area[ico]*cohort_lai[ico]
                vai_max_layer = np.max([vai_max_layer,cohort_lai[ico]+cohort_sai[ico]])
                n_veg = n_veg + 1
                iret = f90.setup_canopy_sub(c_int(ican+1),c_int(icol+1), \
                                            c_int(cohort_pft[ico]), c_double(cohort_area[ico]), \
                                            c_double(cohort_lai[ico]), c_double(cohort_sai[ico]))
                
                elem_diags[ican].append(elem_diags_type(max_vai_above, lai_above, cohort_lai[ico], cohort_sai[ico], \
                                                        dvai, cohort_area[ico], cohort_pft[ico], ntimes))
            else:
                air_pft  = 0
                air_area = (1.-veg_area)/float(n_col-n_veg)
                air_lai  = 1.0 # notional
                air_sai  = 1.0 # notional
                iret = f90.setup_canopy_sub(c_int(ican+1),c_int(icol+1),c_int(air_pft), \
                                            c_double(air_area),c_double(air_lai),c_double(air_sai))

        # Set the current LAI and SAI to the layer above
        max_vai_above = max_vai_above + vai_max_layer
        lai_above = lai_above + lai_current



    # Initialize the site-level diagnostics structure
    site_diags = site_diags_type(ntimes,dvai,max_vai_above)
    for ican in range(n_can):
        for icol in range(n_col):
            ico = elem_cohort[ican,icol]
            if(ico>=0):
                elem_diags[ican][icol],site_diags = GetElemLayerLAIShare(elem_diags[ican][icol],site_diags)
    
    

    visualize_elements = False
    if(visualize_elements):
        fig10, ax = plt.subplots(ncols=1,nrows=1,figsize=(7,7))
        maxvai = 0.
        vai0   = 0.
        for ican in range(n_can):
            veg_area = 0.0
            vai0     = vai0+maxvai
            maxvai   = 0.0
            for icol in range(n_col):
                ico = elem_cohort[ican,icol]
                if(ico>=0):
                    vai = cohort_lai[ico]+cohort_sai[ico]
                    # xy, width height
                    rect = patches.Rectangle((veg_area, vai0), cohort_area[ico], \
                                             vai, linewidth=1, edgecolor='k', facecolor='g')
                    ax.add_patch(rect)
                    ax.text(veg_area+0.5*cohort_area[ico], vai0+0.5*vai, \
                            'cohort {}'.format(ico+1), horizontalalignment='center', \
                            verticalalignment='center')
                    maxvai = np.max([maxvai,vai])
                    veg_area = veg_area + cohort_area[ico]
                    
                    
        plot_vai = 1.2*(vai0+maxvai)
        ax.set_xlabel('Area Fraction [m2/m2]')
        ax.set_ylabel('Integrated LAI+SAI [m2/m2]')
        ax.set_ylim([0,plot_vai])
        ax.set_xlim([0,1])
        ax.invert_yaxis()
        plt.show()


    return n_can,n_col,elem_diags,site_diags
