# ======================================================================
# This module is a supplement to physical processes that dont have a
# dedicated function or routine in the FATES fortran code
#
# Perhaps they should, but they live here for now
#=======================================================================

def GetJmaxKp25Top(vcmax25_top: float):

    # Calculate Jmax and Kp at the canopy top at 25C
    # they scale off of vcmax
    #
    # jmax25_top:  Canopy top maximum electron transport
    #              rate at 25C (umol electrons/m**2/s)
    #
    # kp25top      Canopy top initial slope of CO2 response
    #              curve (C4 plants) at 25C
    
    jmax25_top: float = 1.67   * vcmax25_top
    kp25_top: float   = 20000.  * vcmax25_top
    
    # q10 response of product limited psn.
    # co2_rcurve_islope = co2_rcurve_islope25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
    
    return jmax25_top, kp25_top
