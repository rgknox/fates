import ctypes
from ctypes import *
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb
import code

# Radiation constants, number of broadbands, and their indices (visible and NIR)
n_bands = 2
visb = 1
nirb = 2

def GetParamList(noderoot,node_name,vtype):

    # Read in a list of values from an xml node
    # Convert that string list to the datatype specified
    # and return a list of the new datatype
    
    str_list = noderoot.find(node_name).text.strip().split(',')
    if(vtype=='integer'):
        val_list = [int(str) for str in str_list]
    elif(vtype=='float'):
        val_list = [float(str) for str in str_list]
    elif(vtype=='logical'):
        val_list = [str.lower()=='true' for str in str_list]
    elif(vtype=='string'):
        val_list = str_list
    else:
        print('passed an unknown datatype to GetParamList: {}'.format(vtype))
        exit(1)

    return val_list


def GetParamFromAttrib(noderoot,attr_str):

    # This returns a list of text strings
    
    foundnode = False
    for node in noderoot.iter('param'):
        #print(node.attrib.get('name').strip(),attr_str.strip())
        if (node.attrib.get('name').strip()==attr_str.strip()):
            param_val = node.text.strip().split(',')
            foundnode = True
    if (not foundnode):
        print('Could not find an xml element')
        exit(1)

    return param_val


def PushXMLPhotoParameters(f90,xmlroot):

    numpft = int(xmlroot.find('numpft').text.strip())
    
    # Allocate and push  photosynthesis parameters
    # -------------------------------------------------------------------------------------
    print('Allocating parameter space for {} pfts'.format(numpft))
    iret = f90.alloc_leaf_param_sub(ci(numpft))

    # These are photosynthesis PFT parameters that also need to be pushed to the fortran objects
    f90_photo_pft_params = ['fates_leaf_stomatal_btran_model','fates_leaf_agross_btran_model', \
                            'fates_leaf_c3psn','fates_leaf_stomatal_slope_ballberry', \
                            'fates_leaf_stomatal_slope_medlyn','fates_leaf_fnps', \
                            'fates_leaf_stomatal_intercept','fates_maintresp_reduction_curvature', \
                            'fates_maintresp_reduction_intercept','fates_maintresp_reduction_upthresh', \
                            'fates_maintresp_leaf_atkin2017_baserate','fates_maintresp_leaf_ryan1991_baserate', \
                            'fates_leaf_vcmaxha','fates_leaf_jmaxha', \
                            'fates_leaf_vcmaxhd','fates_leaf_jmaxhd', \
                            'fates_leaf_vcmaxse','fates_leaf_jmaxse']

    
    # Push Parameter File values to the fortran objects, also save some values in local lists
    pft_root = xmlroot.find('f90_params').find('pft_dim')
    for param_name in f90_photo_pft_params:
        print('Pushing parameter: '+param_name)
        for ft in range(numpft):
            param_val = float(GetParamFromAttrib(pft_root,param_name)[ft])
            iret = f90.set_leaf_param_sub(c8(param_val),ci(ft+1),*ccharnb(param_name))

            
def PushXMLRadParameters(f90,xmlroot):

    numpft = int(xmlroot.find('numpft').text.strip())
    
    # Allocate and push radiation parameters
    # -----------------------------------------------------------------------------------
    iret = f90.alloc_radparams_sub(ci(numpft),ci(n_bands))
    for ft in range(numpft):
        pft = ft+1
        #code.interact(local=dict(globals(), **locals()))
        param_val = float(GetParamFromAttrib(pft_root,'fates_rad_leaf_rhovis')[ft])
        
        iret = f90.set_radparams_sub(c8(float(GetParamFromAttrib(pft_root,'fates_rad_leaf_rhovis')[ft]) ), c_int(pft),c_int(visb),*ccharnb("rhol"))
        iret = f90.set_radparams_sub(c8(float(GetParamFromAttrib(pft_root,'fates_rad_leaf_rhonir')[ft]) ), c_int(pft),c_int(nirb),*ccharnb("rhol"))
        iret = f90.set_radparams_sub(c8(float(GetParamFromAttrib(pft_root,'fates_rad_leaf_tauvis')[ft]) ), c_int(pft),c_int(visb),*ccharnb("taul"))
        iret = f90.set_radparams_sub(c8(float(GetParamFromAttrib(pft_root,'fates_rad_leaf_taunir')[ft]) ), c_int(pft),c_int(nirb),*ccharnb("taul"))
        iret = f90.set_radparams_sub(c8(float(GetParamFromAttrib(pft_root,'fates_rad_stem_rhovis')[ft]) ), c_int(pft),c_int(visb),*ccharnb("rhos"))
        iret = f90.set_radparams_sub(c8(float(GetParamFromAttrib(pft_root,'fates_rad_stem_rhonir')[ft]) ), c_int(pft),c_int(nirb),*ccharnb("rhos"))
        iret = f90.set_radparams_sub(c8(float(GetParamFromAttrib(pft_root,'fates_rad_stem_tauvis')[ft]) ), c_int(pft),c_int(visb),*ccharnb("taus"))
        iret = f90.set_radparams_sub(c8(float(GetParamFromAttrib(pft_root,'fates_rad_stem_taunir')[ft]) ), c_int(pft),c_int(nirb),*ccharnb("taus"))
        iret = f90.set_radparams_sub(c8(float(GetParamFromAttrib(pft_root,'fates_rad_leaf_xl')[ft]) ), c_int(pft), c_int(0),*ccharnb("xl"))
        iret = f90.set_radparams_sub(c8(float(GetParamFromAttrib(pft_root,'fates_rad_leaf_clumping_index')[ft]) ), c_int(pft),c_int(0),*ccharnb("clumping_index"))
    


def PushParameters(f90,params,dims):


    numpft = dims['fates_pft']
    
    # Allocate and push  photosynthesis parameters
    # -------------------------------------------------------------------------------------
    print('Allocating parameter space for {} pfts'.format(numpft))
    iret = f90.alloc_leaf_param_sub(ci(numpft))

    # These are photosynthesis PFT parameters that also need to be pushed to the fortran objects
    f90_photo_pft_params = ['fates_leaf_stomatal_btran_model','fates_leaf_agross_btran_model', \
                            'fates_leaf_c3psn','fates_leaf_stomatal_slope_ballberry', \
                            'fates_leaf_stomatal_slope_medlyn','fates_leaf_fnps', \
                            'fates_leaf_stomatal_intercept','fates_maintresp_reduction_curvature', \
                            'fates_maintresp_reduction_intercept','fates_maintresp_reduction_upthresh', \
                            'fates_maintresp_leaf_atkin2017_baserate','fates_maintresp_leaf_ryan1991_baserate', \
                            'fates_leaf_vcmaxha','fates_leaf_jmaxha', \
                            'fates_leaf_vcmaxhd','fates_leaf_jmaxhd', \
                            'fates_leaf_vcmaxse','fates_leaf_jmaxse']

    
    # Push Parameter File values to the fortran objects, also save some values in local lists
    for param_name in f90_photo_pft_params:
        for pft in range(numpft):
            iret = f90.set_leaf_param_sub(c8(float(params[param_name].data[pft])),ci(pft+1),*ccharnb(param_name))

            
    # Allocate and push radiation parameters
    # -----------------------------------------------------------------------------------
    iret = f90.alloc_radparams_sub(ci(numpft),ci(n_bands))
    for ft in range(numpft):
        pft = ft+1
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_rhovis'].data[ft])), c_int(pft),c_int(visb),*ccharnb("rhol"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_rhonir'].data[ft])), c_int(pft),c_int(nirb),*ccharnb("rhol"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_tauvis'].data[ft])), c_int(pft),c_int(visb),*ccharnb("taul"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_taunir'].data[ft])), c_int(pft),c_int(nirb),*ccharnb("taul"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_stem_rhovis'].data[ft])), c_int(pft),c_int(visb),*ccharnb("rhos"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_stem_rhonir'].data[ft])), c_int(pft),c_int(nirb),*ccharnb("rhos"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_stem_tauvis'].data[ft])), c_int(pft),c_int(visb),*ccharnb("taus"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_stem_taunir'].data[ft])), c_int(pft),c_int(nirb),*ccharnb("taus"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_xl'].data[ft])), c_int(pft), c_int(0),*ccharnb("xl"))
        iret = f90.set_radparams_sub(c8(float(params['fates_rad_leaf_clumping_index'].data[ft])), c_int(pft),c_int(0),*ccharnb("clumping_index"))


    
        
