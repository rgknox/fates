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


class f90_modules:

    def __init__(self,mod_path):

        # Instantiate DGESV from lapack
        # ---------------------------------------------------------------------------------------

        self.dgesv_obj = ctypes.CDLL(mod_path+'dgesvMod.o',mode=ctypes.RTLD_GLOBAL)

        
        # LEAF BIOPHYSICS
        #----------------------------------------------------------------------------------------

        # Instantiate the F90 modules
        self.const_obj = ctypes.CDLL(mod_path+'FatesConstantsMod.o',mode=ctypes.RTLD_GLOBAL)
        self.shr_obj = ctypes.CDLL(mod_path+'WrapShrMod.o',mode=ctypes.RTLD_GLOBAL)
        self.fatesutils_obj = ctypes.CDLL(mod_path+'FatesUtilsMod.o',mode=ctypes.RTLD_GLOBAL)
        self.leaf_biophys_obj = ctypes.CDLL(mod_path+'LeafBiophysicsMod.o',mode=ctypes.RTLD_GLOBAL)
        self.leaf_biophys_supp_obj = ctypes.CDLL(mod_path+'LeafBiophysSuppMod.o',mode=ctypes.RTLD_GLOBAL)

        # Identify subroutine objects, so we can call them
        self.set_leaf_param_sub = getattr(self.leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','setleafparam'))
        self.alloc_leaf_param_sub = getattr(self.leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','allocleafparam'))
        self.dealloc_leaf_param_sub = getattr(self.leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','deallocleafparam'))
        self.dump_param_sub =  getattr(self.leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','DumpParams'))
        self.biophysrate_sub = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','LeafLayerBiophysicalRate'))
        self.leaflayerphoto_sub = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','LeafLayerPhotosynthesis'))
        self.qsat_sub = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','QSat'))
        self.cangas_sub = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','GetCanopyGasParameters'))
        self.lmr_ryan_sub = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','LeafLayerMaintenanceRespiration_Ryan_1991'))
        self.lmr_atkin_sub = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','LeafLayerMaintenanceRespiration_Atkin_etal_2017'))
        self.gs_medlyn_sub = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','StomatalCondMedlyn'))
        self.gs_ballberry_sub = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','StomatalCondBallBerry'))
        self.cifunc_sub = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','CiFunc'))
        self.cibisection_sub = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','CiBisection'))
        
        self.velotomolarcf_fun = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','VeloToMolarCF'))
        self.agross_rubiscoc3_fun  = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossRubiscoC3'))
        self.agross_rubpc3_fun  = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossRuBPC3'))
        self.agross_rubpc4_fun  = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossRuBPC4'))
        self.agross_pepc4_fun  = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossPEPC4'))
    
        # For functions, define the return value
        self.set_leaf_param_sub.argtypes = [POINTER(c_double),POINTER(c_int),c_char_p,c_long]
        self.agross_rubiscoc3_fun.restype = c_double
        self.agross_rubpc3_fun.restype = c_double
        self.agross_rubpc4_fun.restype = c_double
        self.agross_pepc4_fun.restype = c_double
        self.velotomolarcf_fun.restype = c_double

        # RADIATION
        #----------------------------------------------------------------------------------------
    
        # Instantiate the F90 modules
        self.shr_obj = ctypes.CDLL(mod_path+'WrapShrMod.o',mode=ctypes.RTLD_GLOBAL)
        self.mem_obj = ctypes.CDLL(mod_path+'FatesRadiationMemMod.o',mode=ctypes.RTLD_GLOBAL)
        self.twostr_obj = ctypes.CDLL(mod_path+'TwoStreamMLPEMod.o',mode=ctypes.RTLD_GLOBAL)
        self.rad_wrap_obj = ctypes.CDLL(mod_path+'RadiationWrapMod.o',mode=ctypes.RTLD_GLOBAL)
        
        # Create aliases for the calls and define arguments if it helps with clarity
        self.alloc_twostream_sub =  getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','initallocate'))
        self.dealloc_twostream_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','dealloc'))
        self.alloc_radparams_sub = getattr(self.twostr_obj,GetModSymbol(mod_path+'TwoStreamMLPEMod.o','allocateradparams'))
        self.set_radparams_sub   = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','setradparam'))
        self.param_prep_sub = getattr(self.twostr_obj,GetModSymbol(mod_path+'TwoStreamMLPEMod.o','radparamprep'))
        self.setup_canopy_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','setupcanopy'))
        self.grndsnow_albedo_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','setgroundsnow'))
        self.canopy_prep_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','wrapcanopyprep'))
        self.zenith_prep_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','wrapzenithprep'))
        self.solver_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','wrapsolve'))
        self.setdown_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','wrapsetdownwelling'))
        self.getintens_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','wrapgetintensity'))
        self.getabsrad_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','wrapgetabsrad'))
        self.getparams_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','wrapgetparams'))
        self.forceparam_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','wrapforceparams'))
                              
        self.forceparam_sub.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_double),c_char_p,c_long]
        self.setup_canopy_sub.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int), \
                                          POINTER(c_double),POINTER(c_double),POINTER(c_double)]
        self.set_radparams_sub.argtypes = [POINTER(c_double),POINTER(c_int),POINTER(c_int),c_char_p,c_long]
        self.grndsnow_albedo_sub.argtypes = [POINTER(c_int),POINTER(c_double),c_char_p,c_long]

