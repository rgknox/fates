import re
import subprocess
import sys
import ctypes
import gc
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


        json_modstr = mod_path+'JSONParameterUtilsMod.o'
        wrapjson_modstr = mod_path+'WrapJSONParameterUtilsMod.o'
        paramint_modstr = mod_path+'FatesParametersInterface.o'
        
        # Instantiate DGESV from lapack
        # ---------------------------------------------------------------------------------------

        self.dgesv_obj = ctypes.CDLL(mod_path+'dgesvMod.o',mode=ctypes.RTLD_GLOBAL)

       

        # Instantiate the F90 modules
        self.const_obj = ctypes.CDLL(mod_path+'FatesConstantsMod.o',mode=ctypes.RTLD_GLOBAL)
        self.shr_obj = ctypes.CDLL(mod_path+'WrapShrMod.o',mode=ctypes.RTLD_GLOBAL)
        self.twostr_obj = ctypes.CDLL(mod_path+'TwoStreamMLPEMod.o',mode=ctypes.RTLD_GLOBAL)
        self.glob_obj = ctypes.CDLL(mod_path+'FatesGlobals.o',mode=ctypes.RTLD_GLOBAL)
        self.json_obj = ctypes.CDLL(json_modstr,mode=ctypes.RTLD_GLOBAL)
        self.paramint_obj = ctypes.CDLL(paramint_modstr,mode=ctypes.RTLD_GLOBAL)
        self.wrapjson_obj = ctypes.CDLL(wrapjson_modstr,mode=ctypes.RTLD_GLOBAL)

        self.intface_obj = ctypes.CDLL(mod_path+'FatesInterfaceTypesMod.o',mode=ctypes.RTLD_GLOBAL)
        self.prt_parameters_obj = ctypes.CDLL(mod_path+'PRTParametersMod.o',mode=ctypes.RTLD_GLOBAL)
        self.fatesutils_obj = ctypes.CDLL(mod_path+'FatesUtilsMod.o',mode=ctypes.RTLD_GLOBAL)
        self.leaf_biophys_obj = ctypes.CDLL(mod_path+'LeafBiophysicsMod.o',mode=ctypes.RTLD_GLOBAL)
        self.leaf_biophys_supp_obj = ctypes.CDLL(mod_path+'LeafBiophysSuppMod.o',mode=ctypes.RTLD_GLOBAL)
        self.edparams_obj       = ctypes.CDLL(mod_path+'EDParamsMod.o',mode=ctypes.RTLD_GLOBAL)
        self.pftparm_obj        = ctypes.CDLL(mod_path+'EDPftvarcon.o',mode=ctypes.RTLD_GLOBAL)
        self.sfparam_obj        = ctypes.CDLL(mod_path+'SFParamsMod.o',mode=ctypes.RTLD_GLOBAL)
        self.fatesderived_obj = ctypes.CDLL(mod_path+'FatesParameterDerivedMod.o',mode=ctypes.RTLD_GLOBAL)
        self.paramsderived_obj = ctypes.CDLL(mod_path+'EDParamsDerivedSupp.o',mode=ctypes.RTLD_GLOBAL)
        self.damage_obj = ctypes.CDLL(mod_path+'DamageMainMod.o',mode=ctypes.RTLD_GLOBAL)
        self.allometry_obj = ctypes.CDLL(mod_path+'FatesAllometryMod.o',mode=ctypes.RTLD_GLOBAL)
        self.kinetics_obj = ctypes.CDLL(mod_path+'FatesKineticsMod.o',mode=ctypes.RTLD_GLOBAL)

        # Kinetics
        self.aqueous_uptake_mm_sub = getattr(self.kinetics_obj,GetModSymbol(mod_path+'FatesKineticsMod.o','AqueousUptakeMM'))
        self.solveconc_nr_mm_sub = getattr(self.kinetics_obj,GetModSymbol(mod_path+'FatesKineticsMod.o','SolveConcNewtRaphMM'))
        
        # LEAF BIOPHYSICS
        #----------------------------------------------------------------------------------------
        
        # Identify subroutine objects, so we can call them
        self.set_leaf_param_sub = getattr(self.leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','setleafparam'))
        self.alloc_leaf_param_sub = getattr(self.leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','allocleafparam'))
        self.dealloc_leaf_param_sub = getattr(self.leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','deallocleafparam'))
        self.dump_param_sub =  getattr(self.leaf_biophys_supp_obj, GetModSymbol(mod_path+'LeafBiophysSuppMod.o','DumpParams'))
        self.isalloc_leaf_param_fun =  getattr(self.leaf_biophys_supp_obj,GetModSymbol(mod_path+'LeafBiophysSuppMod.o','IsLeafParamAllocated'))
        self.isalloc_leaf_param_fun.restype = c_bool
        
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

        self.decaycoeffvcmax_fun = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','DecayCoeffVcmax'))
        self.velotomolarcf_fun = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','VeloToMolarCF'))
        self.agross_rubiscoc3_fun  = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossRubiscoC3'))
        self.agross_rubpc3_fun  = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossRuBPC3'))
        self.agross_rubpc4_fun  = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossRuBPC4'))
        self.agross_pepc4_fun  = getattr(self.leaf_biophys_obj,GetModSymbol(mod_path+'LeafBiophysicsMod.o','AgrossPEPC4'))
       
        # For functions, define the return value
        self.decaycoeffvcmax_fun.restype = c_double
        self.set_leaf_param_sub.argtypes = [POINTER(c_double),POINTER(c_int),c_char_p,c_long]
        self.agross_rubiscoc3_fun.restype = c_double
        self.agross_rubpc3_fun.restype = c_double
        self.agross_rubpc4_fun.restype = c_double
        self.agross_pepc4_fun.restype = c_double
        self.velotomolarcf_fun.restype = c_double

        # JSON
        self.json_setinval = getattr(self.json_obj,GetModSymbol(json_modstr,'JSONSetInvalid'))
        self.json_setloginit = getattr(self.json_obj,GetModSymbol(json_modstr,'JSONSetLogInit'))
        self.json_read  = getattr(self.json_obj,GetModSymbol(json_modstr,'JSONRead'))
        self.json_dumpparam = getattr(self.json_obj,GetModSymbol(json_modstr,'JSONDumpParameter'))

        # WrapJSON
        self.wrapjson_read     = getattr(self.wrapjson_obj,GetModSymbol(wrapjson_modstr,'WrapJSONRead'))
        self.wrapjson_setparam = getattr(self.wrapjson_obj,GetModSymbol(wrapjson_modstr,'WrapJSONSetParameter'))

        self.getparams = getattr(self.paramint_obj,GetModSymbol(paramint_modstr,'GetParameterIndices'))
        
        # RADIATION
        #----------------------------------------------------------------------------------------
    
        # Instantiate the F90 modules
        self.shr_obj = ctypes.CDLL(mod_path+'WrapShrMod.o',mode=ctypes.RTLD_GLOBAL)
        self.mem_obj = ctypes.CDLL(mod_path+'FatesRadiationMemMod.o',mode=ctypes.RTLD_GLOBAL)
       
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
        self.isalloc_radparam_fun =  getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','isradparamallocated'))
        self.isalloc_radparam_fun.restype = c_bool
        self.isalloc_radtype_fun = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','isradtypeallocated'))
        self.isalloc_radtype_fun.restype = c_bool
        self.forceparam_sub.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_double),c_char_p,c_long]
        self.setup_canopy_sub.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int), \
                                          POINTER(c_double),POINTER(c_double),POINTER(c_double)]
        self.set_radparams_sub.argtypes = [POINTER(c_double),POINTER(c_int),POINTER(c_int),c_char_p,c_long]
        self.grndsnow_albedo_sub.argtypes = [POINTER(c_int),POINTER(c_double),c_char_p,c_long]

        # Allometry (parameters via parteh)

       

        # Alias routines
        self.forceparam_sub = getattr(self.rad_wrap_obj,GetModSymbol(mod_path+'RadiationWrapMod.o','wrapforceparams'))
        self.alloc_derived_sub = getattr(self.fatesderived_obj,GetModSymbol(mod_path+'FatesParameterDerivedMod.o','InitAllocate'))
        self.init_derived_sub = getattr(self.paramsderived_obj,GetModSymbol(mod_path+'EDParamsDerivedSupp.o','InitDerived'))

        self.setderived_sub = getattr(self.paramsderived_obj,GetModSymbol(mod_path+'EDParamsDerivedSupp.o','SetDerivedParam'))

        self.h2d_allom_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','h2d_allom'))
        self.h_allom_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','h_allom'))
        self.bagw_allom_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','bagw_allom'))
        self.blmax_allom_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','blmax_allom'))
        self.bleaf_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','bleaf'))
        self.storage_fraction_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','storage_fraction_of_target'))
        self.bsap_allom_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','bsap_allom'))
        self.bbgw_allom_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','bbgw_allom'))
        self.bfineroot_sub  = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','bfineroot'))
        self.bdead_allom_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','bdead_allom'))
        self.carea_allom_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','carea_allom'))
        self.bstore_allom_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','bstore_allom'))
        self.ForceDBH_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','ForceDBH'))
        self.CrownDepth_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','CrownDepth'))
        self.set_root_fraction_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','set_root_fraction'))
        self.leafc_from_treelai_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','leafc_from_treelai'))
        self.tree_lai_sai_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','leafc_from_treelai'))
        self.VegAreaLayer_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','VegAreaLayer'))
        self.tree_lai_sai_sub = getattr(self.allometry_obj,GetModSymbol(mod_path+'FatesAllometryMod.o','tree_lai_sai'))

    def Release(self):

        del self.dgesv_obj
        del self.const_obj
        del self.shr_obj
        del self.twostr_obj
        del self.glob_obj
        del self.json_obj
        del self.pint_obj
        del self.wrapjson_obj
        del self.intface_obj
        del self.prt_parameters_obj
        del self.fatesutils_obj
        del self.leaf_biophys_obj
        del self.leaf_biophys_supp_obj
        del self.edparams_obj
        del self.pftparm_obj
        del self.sfparam_obj
        del self.fatesderived_obj
        del self.paramsderived_obj
        del self.damage_obj
        del self.allometry_obj

        gc.collect()
        
