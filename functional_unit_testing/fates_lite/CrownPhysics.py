import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
import ctypes
import numpy as np
from ctypes import *
from operator import add
sys.path.append('../shared/py_src')
import code  # For development: code.interact(local=dict(globals(), **locals()))

from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb
from PushParameters import GetParamFromAttrib
from PushParameters import GetParamList


# Radiation constants, number of broadbands, and their indices (visible and NIR)
n_bands = 2
visb = 1
nirb = 2


# If Kumerathunga respiration is used, it requires moving averages
# of leaf temperature
t_growth_kum = -999
t_home_kum = -999

dayl_factor_full = 1.0

# Conversion factor for radiant energy in [Watts/m2] to
# photon flux density [umol/m2/s]
wm2_to_umolm2s = 4.6

def CanopyElementPhysics(ican,icol,f90,met,it,xmlroot,pft, \
                         vcmax25_top,jmax25_top,kp25_top,lnc_top, \
                         co2_ppress, o2_ppress, btran, \
                         maintresp_leaf_model,site_diags,elem_diags):

    # This routine assumes the radiation solve has already been performed for the canopy.

    ft = pft - 1   # Python uses zero as the base index, fortran uses 1

    # Vegetation Temperature [K]
    tvegk = met.data['t_veg'][it]
    
    solve_inter_f   = c_double(-9);    
    rd_abs_leaf_f   = c_double(-9);    rb_abs_leaf_f   = c_double(-9);    r_abs_stem_f    = c_double(-9)
    r_abs_snow_f    = c_double(-9);    leaf_sun_frac_f = c_double(-9);    r_diff_dn_f     = c_double(-9)
    rd_abs_f        = c_double(-9);    rb_abs_f        = c_double(-9);    r_diff_up_f     = c_double(-9)
    r_beam_f        = c_double(-9);    

    vcmax_f         = c_double(-9);    jmax_f          = c_double(-9)
    kp_f            = c_double(-9);    agross_f        = c_double(-9);    gstoma_f        = c_double(-9)
    anet_f          = c_double(-9);    lmr_f           = c_double(-9);    c13_f           = c_double(-9)
    ac_f            = c_double(-9);    aj_f            = c_double(-9);    ap_f            = c_double(-9)
    co2_interc_f    = c_double(-9);    mm_kco2_f       = c_double(-9);    mm_ko2_f        = c_double(-9)
    co2_cpoint_f    = c_double(-9);    gs0_f           = c_double(-9);    gs1_f           = c_double(-9)
    gs2_f           = c_double(-9);

    # Set canopy gas parameters, such as the MM coefficients and the CO2 compensation point
    iret = f90.cangas_sub(c8(met.data['can_press'][it]), \
                          c8(o2_ppress), \
                          c8(tvegk), \
                          byref(mm_kco2_f), \
                          byref(mm_ko2_f), \
                          byref(co2_cpoint_f))
    
    # Determine the nitrogen attenuation rate in the canopy
    # -------------------------------------------------------------------------------
    # Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
    # kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al 
    # (2010) Biogeosciences, 7, 1833-1859
    pft_root = xmlroot.find('f90_params').find('pft_dim')
    kn = f90.decaycoeffvcmax_fun(c8(vcmax25_top), \
                                 c8(float(GetParamFromAttrib(pft_root,'fates_leafn_vert_scaler_coeff1')[ft])), \
                                 c8(float(GetParamFromAttrib(pft_root,'fates_leafn_vert_scaler_coeff2')[ft])))

    leaf_frac = elem_diags.total_lai/(elem_diags.total_lai+elem_diags.total_sai)

    
    par_abs_check = 0
    for il in range(elem_diags.n_layer):

        # These area indices are WRT the element's top
        # --------------------------------------------
        vai_top = elem_diags.vai_top[il]
        vai_bot = elem_diags.vai_bot[il]
        dlai    = elem_diags.dlai
        dvai    = vai_bot - vai_top
        vai_mid = 0.5*(vai_bot+vai_top)
        sil     = elem_diags.sil[il]
        
        # Query the solver for the absorbed radiation over this interval
        iret = f90.getabsrad_sub(ci(ican+1),ci(icol+1),ci(visb),c8(vai_top),c8(vai_bot), \
                                 byref(rd_abs_f),byref(rb_abs_f), \
                                 byref(rd_abs_leaf_f),byref(rb_abs_leaf_f),byref(r_abs_stem_f), \
                                 byref(r_abs_snow_f),byref(leaf_sun_frac_f))

        
        # Query the solver for the radiative intensities at the mid-point of this interval
        iret = f90.getintens_sub(ci(ican+1), ci(icol+1), ci(visb), c8(vai_mid), \
                                 byref(r_diff_dn_f), byref(r_diff_up_f), byref(r_beam_f))
        
        
        
        # Nitrogen attenuation at this canopy depth
        nscaler = np.exp(-kn*elem_diags.canopy_lai_mid[il]*leaf_frac)

        # Scale down N and biophysical rates
        iret = f90.biophysrate_sub(ci(pft), c8(vcmax25_top), \
                                   c8(jmax25_top), c8(kp25_top), \
                                   c8(nscaler), c8(tvegk), c8(met.data['dayl_factor'][it]), \
                                   c8(t_growth_kum),c8(t_home_kum),c8(btran), \
                                   byref(vcmax_f), byref(jmax_f), byref(kp_f), \
                                   byref(gs0_f), byref(gs1_f), byref(gs2_f))

        # Leaf Maintenance Respiration (temp and pft dependent)
        if(maintresp_leaf_model==1):
            iret = f90.lmr_ryan_sub(c8(lnc_top), c8(nscaler), ci(pft), c8(tvegk), byref(lmr_f))
        elif(maintresp_leaf_model==2):
            iret = f90.lmr_atkin_sub(c8(lnc_top), c8(rdark_scaler), c8(tvegk), \
                                     c8(atkin_mean_leaf_tempk), byref(lmr_f) )
        else:
            print('unknown leaf respiration model')
            exit(1)

        
        par_abs_check = par_abs_check + (r_abs_stem_f.value+r_abs_snow_f.value)
        for ipar in [0,1]:
            if(ipar==0):
                lit_frac = leaf_sun_frac_f.value
                par_abs_leaf_umol = (rb_abs_leaf_f.value/(dlai*lit_frac) + \
                                     rd_abs_leaf_f.value/dlai)*wm2_to_umolm2s
                    
            else:
                lit_frac = 1.- leaf_sun_frac_f.value
                par_abs_leaf_umol = (rd_abs_leaf_f.value/dlai)*wm2_to_umolm2s

            # Energy conservation check (diff absorbed by leaf, stem and snow)
            par_abs_check = par_abs_check + par_abs_leaf_umol*lit_frac*dlai/wm2_to_umolm2s


            # Main call to FATES photosynthesis
            
            iret = f90.leaflayerphoto_sub( c8(par_abs_leaf_umol),  \
                                           ci(pft),   \
                                           c8(vcmax_f.value),   \
                                           c8(jmax_f.value),    \
                                           c8(kp_f.value),      \
                                           c8(gs0_f.value), \
                                           c8(gs1_f.value), \
                                           c8(gs2_f.value), \
                                           c8(tvegk), \
                                           c8(met.data['can_press'][it]), \
                                           c8(co2_ppress), \
                                           c8(o2_ppress), \
                                           c8(met.data['vpress_sat'][it]), \
                                           c8(met.data['g_b_umol'][it]), \
                                           c8(met.data['vpress'][it]), \
                                           c8(mm_kco2_f.value), \
                                           c8(mm_ko2_f.value), \
                                           c8(co2_cpoint_f.value), \
                                           c8(lmr_f.value), \
                                           c8(0.5), \
                                           byref(agross_f), \
                                           byref(gstoma_f), \
                                           byref(anet_f), \
                                           byref(c13_f), \
                                           byref(co2_interc_f), \
                                           byref(solve_inter_f))


            agross_rubisco = f90.agross_rubiscoc3_fun(c8(vcmax_f.value), c8(co2_interc_f.value), \
                                                      c8(o2_ppress), c8(co2_cpoint_f.value), \
                                                      c8(mm_kco2_f.value),c8(mm_ko2_f.value) )
                
            agross_rubp  = f90.agross_rubpc3_fun(c8(par_abs_leaf_umol), c8(jmax_f.value), \
                                                   c8(float(GetParamFromAttrib(pft_root,'fates_leaf_fnps')[ft])), \
                                                   c8(co2_interc_f.value),c8(co2_cpoint_f.value))


            # We don't present the element diagnostic updater here to reduce verbosity in the main routine
            elem_diags.Update(il,lit_frac,r_diff_dn_f.value,r_diff_up_f.value,r_beam_f.value,rd_abs_leaf_f.value,\
                              rb_abs_leaf_f.value, r_abs_stem_f.value,leaf_sun_frac_f.value, \
                              lmr_f.value,agross_f.value,gstoma_f.value,co2_interc_f.value )


            # Track some diagnostics
            
            # Fraction that this leaf portion contributes to averages amongst the column
            leaf_per_ground    = lit_frac*dlai*elem_diags.crown_area_frac

            # 
            if(site_diags.lai_ground[sil]<0.000000001):
                print(sil,site_diags.lai_ground[sil],il,ican,icol,site_diags.lai_ground[1])
                exit(0)

            
            # Fraction that this leaf portion contributes to averages amongst the layer
            if(dlai>0.):
                leaf_site_frac = leaf_per_ground/site_diags.lai_ground[sil]
            else:
                leaf_site_frac = 0.
                
            site_diags.lmr[it]        = site_diags.lmr[it] + leaf_per_ground*lmr_f.value
            site_diags.agross[it]     = site_diags.agross[it] + leaf_per_ground*agross_f.value
            site_diags.gstoma[it]     = site_diags.gstoma[it] + leaf_per_ground*gstoma_f.value
            site_diags.anet[it]       = site_diags.anet[it] + leaf_per_ground*anet_f.value
            site_diags.r_abs_leaf[it] = site_diags.r_abs_leaf[it] + leaf_per_ground*par_abs_leaf_umol/wm2_to_umolm2s

            site_diags.sunfrac_vl[it,sil]    = site_diags.sunfrac_vl[it,sil] + leaf_site_frac*leaf_sun_frac_f.value
            site_diags.lmr_vl[it,sil]        = site_diags.lmr_vl[it,sil] + leaf_site_frac*lmr_f.value
            site_diags.agross_vl[it,sil]     = site_diags.agross_vl[it,sil] + leaf_site_frac*agross_f.value
            site_diags.gstoma_vl[it,sil]     = site_diags.gstoma_vl[it,sil] + leaf_site_frac*gstoma_f.value
            site_diags.anet_vl[it,sil]       = site_diags.anet_vl[it,sil] + leaf_site_frac*anet_f.value
            site_diags.r_abs_leaf_vl[it,sil] = site_diags.r_abs_leaf_vl[it,sil] + leaf_site_frac*par_abs_leaf_umol/wm2_to_umolm2s
            
            #site_diags.r_abs_leaf_vl[it,sil] = site_diags.r_abs_leaf_vl[it,sil] + leaf_site_frac*(rd_abs_leaf_f.value+rb_abs_leaf_f.value)/dlai
            
            site_diags.vcmax_vl[it,sil]      = site_diags.vcmax_vl[it,sil] + leaf_site_frac*vcmax_f.value
            site_diags.ag_rubisco_vl[it,sil] = site_diags.ag_rubisco_vl[it,sil] + leaf_site_frac*agross_rubisco
            site_diags.ag_rubp_vl[it,sil]    = site_diags.ag_rubp_vl[it,sil] + leaf_site_frac*agross_rubp
            site_diags.co2_interc_vl[it,sil] = site_diags.co2_interc_vl[it,sil] + leaf_site_frac*co2_interc_f.value

            
            

    # Display instantaneous profiles
    # Complete leaf layer and sun/shade loops first
    do_profiles = GetParamList(xmlroot.find('analysis_controls'),'view_element_inst_profs','logical')[0]
    if(do_profiles):
        elem_diags.PlotInstProfiles()

        
                

    return site_diags,elem_diags
