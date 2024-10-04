 ustar                  => frictionvel_inst%ustar_patch                 , & ! Output: [real(r8) (:)   ]  friction velocity [m/s]
         um                     => frictionvel_inst%um_patch                    , & ! Output: [real(r8) (:)   ]  wind speed including the stablity effect [m/s]
         uaf                    => frictionvel_inst%uaf_patch                   , & ! Output: [real(r8) (:)   ]  canopy air speed [m/s]


def MolarToVeloCF(press,tempk):

    umol_per_kmol = 1000.0
    rgas = 8314.4598
    
    cf = press/(rgas * tempk )*umol_per_kmol

    return cf




def GetRbFromUcanDleaf(u_can,dleaf):

    # From the CTSM Technical Manual 2.5.125:
    # The leaf boundary layer resistance [s/m]
    # is a function of the turbulent transfer coefficient (C_v) between
    # the canopy surface and the canopy air 0.01 [m s^(-0.5)]
    # The leaf diameter d_leaf [m] and U_av is "the wind speed
    # incident on the leaves" (which equals the friction velocity u*)

    c_v = 0.01
    
    r_b_ms = (1./C_v)*np.sqrt(u_can / d_leaf)

    g_b_mol = 
    

def MoninObukIni( ur, thv, dthv, zldis, z0m, um, obu7 ):

    
def UCanFromUstar(u_m,u_star,ram1):


    u_can = u_m*np.sqrt( 1./(ram1*u_m) )



