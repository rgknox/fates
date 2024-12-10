import re
import subprocess
import sys
import ctypes
from ctypes import *
from operator import add
sys.path.append('../shared/py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb

# ======================================================================================
#
# This module performs various ctypes related prep work
# for all modules associated with this style of FATES unit tests
# 
#
# 1) it instantiates all of the pre-compiled shared objects
# 2) it associates these objects with an alias
# 3) for some functions, it defines the argument types which makes
#    calling them a little easier
#
# There is a helper function here called GetModSYmbol, which uses
# the system command "nm" to look through the shared object file
# and find the correct symbol.
#
#
# Subroutines
# =======================================================================================

def GetModSymbol(mod_path,symbol):

    # This routine uses some minor string magic and a system call using
    # "nm" to list the symbols inside a fortran object.  Symbols are simply
    # the text strings the define the routines and data structures inside
    # the compiled objects. When symbols are created inside a compiled object,
    # different compilers do different things, like force to lower or upper
    # case, or attach trailing or preceding underscores.  This routine
    # helps identify the exact symbol name, from the name of the fortran
    # routine (as it appears in code) that we want to match with.
    
    subprocess.run(["nm", mod_path],capture_output=True)
    
    list = str(subprocess.run(["nm",mod_path],capture_output=True)).split()
    p = re.compile(symbol, re.IGNORECASE)
    modlist = [ s for s in list if p.search(s)]

    # If nothing came up, thats bad
    if (len(modlist)==0):
        print('Failed to find the right module symbol for:{} in module {}'.format(symbol,mod_path))
        print(list)
        print('Exiting')
        exit(2)

    # Its possible a routine name could also be part of another routine name,
    # such as alloc_such_and_such  vs dealloc_such_and_such, we want the shorter
    mod_symbols = []
    slen = 10000
    for item in modlist:
        ms_str = item.split('\\')[0]
        if len(ms_str)<slen:
            mod_symbol = ms_str
            slen = len(ms_str)

        
    return mod_symbol


# =======================================================================================


def InstantiateF90Modules(abs_fates_path):

    # abs_fates_path is the absolute path to the fates directory

    mod_path = fates_path+'functional_unit_testing/shared/bld/'
    
    # Instantiate DGESV from lapack
    # ---------------------------------------------------------------------------------------
    f_dgesv_obj = ctypes.CDLL(mod_path+'dgesvMod.o',mode=ctypes.RTLD_GLOBAL)


    # LEAF BIOPHYSICS
    #----------------------------------------------------------------------------------------

    # Instantiate the F90 modules
    f90_const_obj = ctypes.CDLL(mod_path+'FatesConstantsMod.o',mode=ctypes.RTLD_GLOBAL)
    f90_shr_obj = ctypes.CDLL(mod_path+'WrapShrMod.o',mode=ctypes.RTLD_GLOBAL)
    f90_fatesutils_obj = ctypes.CDLL(mod_path+'FatesUtilsMod.o',mode=ctypes.RTLD_GLOBAL)
    f90_leaf_biophys_obj = ctypes.CDLL(mod_path+'LeafBiophysicsMod.o',mode=ctypes.RTLD_GLOBAL)
    f90_leaf_biophys_supp_obj = ctypes.CDLL(mod_path+'LeafBiophysSuppMod.o',mode=ctypes.RTLD_GLOBAL)

    # Identify subroutine objects, so we can call them
    fsub_set_leaf_param = getattr(f90_leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','setleafparam'))
    fsub_alloc_leaf_param = getattr(f90_leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','allocleafparam'))
    fsub_dealloc_leaf_param = getattr(f90_leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','deallocleafparam'))
    fsub_dump_param =  getattr(f90_leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','DumpParams'))
    fsub_set_leaf_param.argtypes = [POINTER(c_double),POINTER(c_int),c_char_p,c_long]
    fsub_biophysrate = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','LeafLayerBiophysicalRate'))
    fsub_leaflayerphoto = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','LeafLayerPhotosynthesis'))
    fsub_qsat = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','QSat'))
    fsub_cangas = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','GetCanopyGasParameters'))
    fsub_lmr_ryan = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','LeafLayerMaintenanceRespiration_Ryan_1991'))
    fsub_lmr_atkin = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','LeafLayerMaintenanceRespiration_Atkin_etal_2017'))
    ffun_agross_rubiscoc3  = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossRubiscoC3'))
    ffun_agross_rubpc3  = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossRuBPC3'))
    ffun_agross_rubpc4  = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossRuBPC4'))
    ffun_agross_pepc4  = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossPEPC4'))
    
    fsub_gs_medlyn = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','StomatalCondMedlyn'))
    fsub_gs_ballberry = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','StomatalCondBallBerry'))
    fsub_velotomolarcf = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','VeloToMolarCF'))
    fsub_cifunc = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','CiFunc'))
    fsub_cibisection = getattr(f90_leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','CiBisection'))

    # For functions, define the return value
    ffun_agross_rubiscoc3.restype = c_double
    ffun_agross_rubpc3.restype = c_double
    ffun_agross_rubpc4.restype = c_double
    ffun_agross_pepc4.restype = c_double
    ffun_velotomolarcf.restype = c_double

    # RADIATION
    #----------------------------------------------------------------------------------------
    
    # Instantiate the F90 modules
    f90_shr_obj = ctypes.CDLL(mod_path+'WrapShrMod.o',mode=ctypes.RTLD_GLOBAL)
    f90_mem_obj = ctypes.CDLL(mod_path+'FatesRadiationMemMod.o',mode=ctypes.RTLD_GLOBAL)
    f90_twostr_obj = ctypes.CDLL(mod_path+'TwoStreamMLPEMod.o',mode=ctypes.RTLD_GLOBAL)
    f90_wrap_obj = ctypes.CDLL(mod_path+'RadiationWrapMod.o',mode=ctypes.RTLD_GLOBAL)

    # Create aliases for the calls and define arguments if it helps with clarity
    fsub_alloc_twostream =  getattr(f90_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','initallocate')
    fsub_dealloc_twostream = getattr(f90_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','dealloc')
    fsub_alloc_radparams = getattr(f90_twostr_obj,GetModSymbol(mod_path+'TwoStreamMLPEMod.o','allocateradparams')
    ffun_set_radparams   = f90_wrap_obj.__radiationwrapmod_MOD_setradparam
    ffun_set_radparams.argtypes = [POINTER(c_double),POINTER(c_int),POINTER(c_int),c_char_p,c_long]

    fsub_param_prep = f90_twostr_obj.__twostreammlpemod_MOD_radparamprep
    fsub_setup_canopy = f90_wrap_obj.__radiationwrapmod_MOD_setupcanopy
    fsub_setup_canopy_call.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int), \
                                  POINTER(c_double),POINTER(c_double),POINTER(c_double)]

    fsub_grndsnow_albedo = f90_wrap_obj.__radiationwrapmod_MOD_setgroundsnow
    fsub_grndsnow_albedo.argtypes = [POINTER(c_int),POINTER(c_double),c_char_p,c_long]

    fsub_canopy_prep = f90_wrap_obj.__radiationwrapmod_MOD_wrapcanopyprep
    fsub_zenith_prep = f90_wrap_obj.__radiationwrapmod_MOD_wrapzenithprep
    fsub_solver = f90_wrap_obj.__radiationwrapmod_MOD_wrapsolve
    fsub_setdown = f90_wrap_obj.__radiationwrapmod_MOD_wrapsetdownwelling
    
    fsub_getintens = f90_wrap_obj.__radiationwrapmod_MOD_wrapgetintensity
    fsub_getabsrad = f90_wrap_obj.__radiationwrapmod_MOD_wrapgetabsrad
    fsub_getparams = f90_wrap_obj.__radiationwrapmod_MOD_wrapgetparams
    fsub_forceparam = f90_wrap_obj.__radiationwrapmod_MOD_wrapforceparams
    fsub_forceparam.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_double),c_char_p,c_long]
