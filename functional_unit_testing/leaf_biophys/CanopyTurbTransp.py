# ustar                  => frictionvel_inst%ustar_patch    friction velocity [m/s]
# um                     => frictionvel_inst%um_patch       wind speed including the stablity effect [m/s]
# uaf                    => frictionvel_inst%uaf_patch            canopy air speed [m/s]


def MolarToVeloCF(press,tempk):

    umol_per_kmol = 1000.0
    rgas = 8314.4598
    
    cf = press/(rgas * tempk )*umol_per_kmol

    return cf




def GetRbFromUcanDleaf(u_can,dleaf,press,tempk):

    # From the CTSM Technical Manual 2.5.125:
    # The leaf boundary layer resistance [s/m]
    # is a function of the turbulent transfer coefficient (C_v) between
    # the canopy surface and the canopy air 0.01 [m s^(-0.5)]
    # The leaf diameter d_leaf [m] and U_av is "the wind speed
    # incident on the leaves" (which equals the friction velocity u*)

    c_v = 0.01
    
    r_b_ms = (1./C_v)*np.sqrt(u_can / d_leaf)

    g_b_mol = (1./r_b_ms)/MolarToVeloCF(press,tempk)
    

def MoninObukIni( ur, thv, dthv, zldis, z0m, um, obu ):

    # Initial values of u* and convective velocity
    # ur    ! wind speed at reference height [m/s]
    # thv   ! virtual potential temperature (kelvin)
    # dthv  ! diff of vir. poten. temp. between ref. height and surface
    # zldis ! reference height "minus" zero displacement heght [m]
    # z0m   ! roughness length, momentum [m]
    # um    ! wind speed including the stability effect [m/s]

    # obu   ! monin-obukhov length (m)
    # wc    ! convective velocity [m/s]
    # rib   ! bulk Richardson number
    # zeta  ! dimensionless height used in Monin-Obukhov theory
    # ustar ! friction velocity [m/s]
    
    ustar = 0.06
    wc = 0.5
    
    um=np.max(ur,0.1)

    rib=grav*zldis*dthv/(thv*um*um)

    if (rib >= 0._r8) then      ! neutral or stable
       zeta = rib*log(zldis/z0m)/(1._r8-5._r8*min(rib,0.19_r8))
       zeta = min(this%zetamaxstable,max(zeta,0.01_r8 ))
    else                     ! unstable
       zeta=rib*log(zldis/z0m)
       zeta = max(-100._r8,min(zeta,-0.01_r8 ))
    endif

    obu=zldis/zeta
    
    
    
def UCanFromUstar(u_m,u_star,ram1):


    u_can = u_m*np.sqrt( 1./(ram1*u_m) )



