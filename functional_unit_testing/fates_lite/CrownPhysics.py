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

    izen = np.sum(czbin < met.data['cosz'][it] for czbin in elem_diags.zen_bins)-1
    
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
        
        elem_diags.rd_abs_leaf[il] = rd_abs_leaf_f.value
        elem_diags.rb_abs_leaf[il] = rb_abs_leaf_f.value
        elem_diags.r_abs_stem[il]  = r_abs_stem_f.value
        elem_diags.sunfrac[il]     = leaf_sun_frac_f.value
        
        elem_diags.sunfrac_zen_ll[il,izen] = elem_diags.sunfrac_zen_ll[il,izen] + leaf_sun_frac_f.value
        elem_diags.count_zen_ll[il,izen] = elem_diags.count_zen_ll[il,izen] + 1
        
        # Query the solver for the radiative intensities at the mid-point of this interval
        iret = f90.getintens_sub(ci(ican+1), ci(icol+1), ci(visb), c8(vai_mid), \
                                 byref(r_diff_dn_f), byref(r_diff_up_f), byref(r_beam_f))

        elem_diags.rd_dn[il] = r_diff_dn_f.value
        elem_diags.rd_up[il] = r_diff_up_f.value
        elem_diags.rbeam[il] = r_beam_f.value

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

        site_diags.sunfrac_vl[it,sil] = site_diags.sunfrac_vl[it,sil] + \
            elem_diags.crown_area_frac*dlai/site_diags.lai_ground[sil]*1. #leaf_sun_frac_f.value CHECK TO SEE IF SUMS TO 1
        
        par_abs_check = par_abs_check + (r_abs_stem_f.value+r_abs_snow_f.value)
        for ipar in [0,1]:
            if(ipar==0):
                lit_frac = leaf_sun_frac_f.value
                par_abs_leaf_umol = (elem_diags.rb_abs_leaf[il]/(dlai*lit_frac) + \
                                     elem_diags.rd_abs_leaf[il]/dlai)*wm2_to_umolm2s
                    
            else:
                lit_frac = 1.- leaf_sun_frac_f.value
                par_abs_leaf_umol = (elem_diags.rd_abs_leaf[il]/dlai)*wm2_to_umolm2s

            # Energy conservation check (diff absorbed by leaf, stem and snow)
            par_abs_check = par_abs_check + par_abs_leaf_umol*lit_frac*dlai/wm2_to_umolm2s
                    
            if(par_abs_leaf_umol>400000):
                print('high par_abs_leaf_umol:',ipar,il,met.data['visbdn'][it]*wm2_to_umolm2s, \
                      met.data['visddn'][it]*wm2_to_umolm2s, \
                      rb_abs_leaf[il]*wm2_to_umolm2s,dalai,lit_frac,rd_abs_leaf[il]*wm2_to_umolm2s)

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

            
            leaf_per_ground           = lit_frac*dlai*elem_diags.crown_area_frac
            site_diags.lmr[it]        = site_diags.lmr[it] + leaf_per_ground*lmr_f.value
            site_diags.agross[it]     = site_diags.agross[it] + leaf_per_ground*agross_f.value
            site_diags.gstoma[it]     = site_diags.gstoma[it] + leaf_per_ground*gstoma_f.value
            site_diags.anet[it]       = site_diags.anet[it] + leaf_per_ground*anet_f.value
            site_diags.r_abs_leaf[it] = site_diags.r_abs_leaf[it] + leaf_per_ground*par_abs_leaf_umol/wm2_to_umolm2s

            # Fraction of the site-layer's leaf area  occupied by this element layer's leaves
            if(site_diags.lai_ground[sil]<0.000000001):
                print(sil,site_diags.lai_ground[sil],il,ican,icol,site_diags.lai_ground[1])
                exit(0)
            
            if(dlai>0.):
                leaf_site_frac = lit_frac*elem_diags.crown_area_frac*dlai/site_diags.lai_ground[sil]
            else:
                leaf_site_frac = 0.
                
            site_diags.lmr_vl[it,sil] = site_diags.lmr_vl[it,sil] + leaf_site_frac*lmr_f.value
            
            site_diags.agross_vl[it,sil] = site_diags.agross_vl[it,sil] + leaf_site_frac*agross_f.value
            site_diags.gstoma_vl[it,sil] = site_diags.gstoma_vl[it,sil] + leaf_site_frac*gstoma_f.value
            site_diags.anet_vl[it,sil] = site_diags.anet_vl[it,sil] + leaf_site_frac*anet_f.value
            site_diags.r_abs_leaf_vl[it,sil] = site_diags.r_abs_leaf_vl[it,sil] + leaf_site_frac*par_abs_leaf_umol/wm2_to_umolm2s
            site_diags.vcmax_vl[it,sil] = site_diags.vcmax_vl[it,sil] + leaf_site_frac*vcmax_f.value
                
            agross_rubisco = f90.agross_rubiscoc3_fun(c8(vcmax_f.value), c8(co2_interc_f.value), \
                                                      c8(o2_ppress), c8(co2_cpoint_f.value), \
                                                      c8(mm_kco2_f.value),c8(mm_ko2_f.value) )
                
            agross_rubpc3  = f90.agross_rubpc3_fun(c8(par_abs_leaf_umol), c8(jmax_f.value), \
                                                   c8(float(GetParamFromAttrib(pft_root,'fates_leaf_fnps')[ft])), \
                                                   c8(co2_interc_f.value),c8(co2_cpoint_f.value))
            
            site_diags.agross_rubisco_vl[it,sil] =  site_diags.agross_rubisco_vl[it,sil] + leaf_site_frac*agross_rubisco
            site_diags.agross_rubpc3_vl[it,sil] =  site_diags.agross_rubpc3_vl[it,sil] + leaf_site_frac*agross_rubpc3
            
            if(agross_rubisco<agross_rubpc3):
                elem_diags.ag_limit[il,0] = elem_diags.ag_limit[il,0] + lit_frac
                elem_diags.ag_sslimit[il,ipar,0] = elem_diags.ag_sslimit[il,ipar,0] + lit_frac

            else:
                elem_diags.ag_limit[il,1] = elem_diags.ag_limit[il,1] + lit_frac
                elem_diags.ag_sslimit[il,ipar,1] = elem_diags.ag_sslimit[il,ipar,1] + lit_frac


    # Display instantaneous profiles
    # Complete leaf layer and sun/shade loops first
    ctrl_root = xmlroot.find('analysis_controls')
    do_profiles = GetParamList(ctrl_root,'view_element_inst_profs','logical')[0]
    if(do_profiles):
        fig_edp1,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(6.0,6.5))
        y  = elem_diags.vai_top
        x  = elem_diags.rd_abs_leaf
        ax1.plot(x,y)
        ax1.invert_yaxis()
        ax1.set_ylabel('VAI')
        ax1.set_xlabel('Absorbed Rd umol/m2/s')
        ax1.grid('on')
        x = elem_diags.rb_abs_leaf
        ax2.plot(x,y)
        ax2.invert_yaxis()
        ax2.set_ylabel('VAI')
        ax2.set_xlabel('Absorbed Rb umol/m2/s')
        ax2.grid('on')
        x = elem_diags.sunfrac_v2
        ax3.plot(x,y)
        ax3.invert_yaxis()
        ax3.set_ylabel('VAI')
        ax3.set_xlabel('Sunlit Fraction umol/m2/s')
        ax3.grid('on')

        txt_str = "Time: {:04d}-{:02d}-{:02d}\nCan: {:2d}\nCol: {:2d}\nCos(z) = {:4.3f}\nT_v = {:3.1f}\n".format(met.data['yr'][it],met.data['mon'][it],met.data['day'][it],ican,icol,    met.data['cosz'][it],\
                                                             met.data['t_veg'][it])
        ax4.text(0.2,0.2,txt_str)
        ax4.axis('off')
        plt.tight_layout()
        plt.show()                
        
                

    return site_diags,elem_diags
