module LeafBiophysicsMod

  !-------------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! 
  ! This module contains routines for leaf-level biophyics.  These
  ! routines are all called with primitive arguments to facilitate
  ! use accross models, with the exception of an internally defined
  ! set of constants associated with plant functional type.
  !
  ! ------------------------------------------------------------------------------------

  use shr_log_mod,       only : errMsg => shr_log_errMsg
  use shr_sys_mod,       only : shr_sys_abort
  use FatesConstantsMod, only : r8 => fates_r8
  use shr_infnan_mod,    only : shr_infnan_isnan
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesGlobals,      only : fates_log
  use FatesGlobals,      only : FatesWarn,N2S,A2S,I2S
  use FatesConstantsMod, only : itrue
  use FatesConstantsMod, only : nearzero
  use FatesConstantsMod, only : molar_mass_ratio_vapdry
  use FatesConstantsMod, only : molar_mass_water
  use FatesConstantsMod, only : rgas_J_K_mol
  use FatesConstantsMod, only : fates_unset_r8
  use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
  use FatesConstantsMod, only : wm2_to_umolm2s
  use FatesConstantsMod, only : nocomp_bareground
  use FatesConstantsMod, only : photosynth_acclim_model_none
  use FatesConstantsMod, only : photosynth_acclim_model_kumarathunge_etal_2019
  use FatesConstantsMod, only : lmrmodel_ryan_1991
  use FatesConstantsMod, only : lmrmodel_atkin_etal_2017
  use FatesConstantsMod, only : kpa_per_pa 
  use FatesUtilsMod,     only : QuadraticRoots => QuadraticRootsSridharachary
  
  
  implicit none
  private

  public :: LeafLayerPhotosynthesis

  character(len=*), parameter, private :: sourcefile = &
       __FILE__


  character(len=1024) :: warn_msg   ! for defining a warning message

  !-------------------------------------------------------------------------------------

  ! maximum stomatal resistance [s/m] (used across several procedures)
  real(r8),parameter :: rsmax0 =  2.e8_r8
  
  ! minimum Leaf area to solve, too little has shown instability
  real(r8), parameter :: min_la_to_solve = 0.0000000001_r8

  ! Set this to true to perform debugging
  logical,parameter   ::  debug = .false.

  ! Set this to true to remove leaf boundary layer effects on calculation
  ! of stomatal conductance and net assimilation
  logical, parameter :: zero_bl_resist = .true.
  
  ! Ratio of H2O/CO2 gas diffusion in stomatal airspace (approximate)
  real(r8),parameter :: h2o_co2_stoma_diffuse_ratio = 1.6_r8

  
  ! Ratio of H2O/CO2 gass diffusion in the leaf boundary layer (approximate)
  real(r8),parameter :: h2o_co2_bl_diffuse_ratio = 1.4_r8


  ! Constants used to define C3 versus C4 photosynth pathways
  integer, parameter :: c3_path_index = 1
  integer, parameter :: c4_path_index = 0


  ! Constants used to define conductance models
  integer, parameter :: medlyn_model = 2
  integer, parameter :: ballberry_model = 1

  ! Alternatively, Gross Assimilation can be used to estimate
  ! leaf co2 partial pressure and therefore conductance. The default
  ! is to use anet
  integer, parameter :: net_assim_model = 1
  integer, parameter :: gross_assim_model = 2



  ! These are mutable parameter constants, some are differentiated by PFT, others are not
  ! -------------------------------------------------------------------------------------
  
  type, public :: leafbiophys_params_type

     integer, allocatable :: photopath(:)                         ! Photosynthetic pathway index (C3,C4,etc)
     real(r8),allocatable :: medlyn_slope(:)                      ! Stomatal Slope, Medlyn, e.g. g1 [-]
     real(r8),allocatable :: bb_slope(:)                          ! Stomatal Slope, Ball-Berry, e.g. g1 [-]
     real(r8),allocatable :: stomatal_intercept(:)                ! Stomatal int, BB or Medlyn, e.g. g0, [-]
     real(r8),allocatable :: maintresp_leaf_ryan1991_baserate(:)  ! Base maintenance resp rate M.Ryan 1991 [gC gN-1 s-1]
     real(r8),allocatable :: maintresp_leaf_atkin2017_baserate(:) ! Base maintenance resp rate Atkin 2017 [umol CO2 m-2 s-1]
     real(r8),allocatable :: vcmaxha(:)                           ! activation energy for vcmax (J/mol)
     real(r8),allocatable :: jmaxha(:)                            ! activation energy for jmax (J/mol)
     real(r8),allocatable :: vcmaxhd(:)                           ! deactivation energy for vcmax (J/mol)
     real(r8),allocatable :: jmaxhd(:)                            ! deactivation energy for jmax (J/mol)
     real(r8),allocatable :: vcmaxse(:)                           ! entropy term for vcmax (J/mol/K)
     real(r8),allocatable :: jmaxse(:)                            ! entropy term for jmax (J/mol/K)
       
  end type leafbiophys_params_type

  type(leafbiophys_params_type),public :: lb_params


  
  ! A possible sequence of calls for leaf biophysics is as follows:
  ! 1) determine any gas quantities or parameters that are derived and
  !    are applicable to a super-set of leaf-layers (like MM, and compensation points)
  ! 2) loop over discrete portions of leaves, perhaps differentiated by spatial position and pft
  ! 3) determine if this particular leaf is present and has been solved
  ! 4) solve for leaf maintenance respiration (e.g. dark respiration) 
  ! 5) update leaf-level rates, such as vcmax, jmax
  ! 6) solve for photosynthesis
  ! 7) solve for maintenance respiration of other tissues (other library?)
  
contains

  subroutine StomatalCondMedlyn(anet,ft,veg_esat,ceair,stomatal_intercept_btran,leaf_co2_ppress,can_press, gb_mol,gs_mol)

    ! Input
    real(r8), intent(in) :: anet                     ! net leaf photosynthesis (umol CO2/m**2/s)
    integer, intent(in)  :: ft                       ! plant functional type index
    real(r8), intent(in) :: veg_esat                 ! saturation vapor pressure at veg_tempk (Pa)
    real(r8), intent(in) :: can_press                ! Air pressure NEAR the surface of the leaf (Pa)
    real(r8), intent(in) :: gb_mol                   ! leaf boundary layer conductance (umol /m**2/s)
    real(r8), intent(in) :: ceair                    ! vapor pressure of air, constrained (Pa)
    real(r8), intent(in) :: stomatal_intercept_btran ! water-stressed minimum stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(in) :: leaf_co2_ppress          ! CO2 partial pressure at leaf surface (Pa)

    ! Output
    real(r8) :: gs_mol            ! leaf stomatal conductance (umol H2O/m**2/s)

    ! locals
    real(r8) :: vpd                  ! water vapor deficit in Medlyn stomatal model (KPa)
    real(r8) :: term                 ! intermediate term used to simplify equations
    real(r8) :: aquad,bquad,cquad    ! quadradic solve terms
    real(r8) :: r1,r2                ! quadradic solve roots

    ! Evaluate trival solution, if there is no positive net assimiolation
    ! the stomatal conductance is the intercept conductance
    if (anet <= nearzero) then
       gs_mol = stomatal_intercept_btran
       return
    end if

    ! stomatal conductance calculated from Medlyn et al. (2011), the numerical &
    ! implementation was adapted from the equations in CLM5.0 [kPa]

    vpd =  max((veg_esat - ceair), 50._r8) * kpa_per_pa       !addapted from CLM5. Put some constraint on VPD

    if(zero_bl_resist) then
       
       ! We assume zero resistance in the leaf boundary layer, and that humidity at
       ! the leaf surface is equal to humidity outside the boundary layer
       gs_mol = h2o_co2_stoma_diffuse_ratio*(1._r8 + lb_params%medlyn_slope(ft)/sqrt(vpd))*anet/leaf_co2_ppress*can_press + stomatal_intercept_btran

    else
       
       !when Medlyn stomatal conductance is being used, the unit is KPa. Ignoring the constraint will cause errors when model runs.
       term = h2o_co2_stoma_diffuse_ratio * anet / (leaf_co2_ppress / can_press)
       aquad = 1.0_r8
       bquad = -(2.0 * (stomatal_intercept_btran+ term) + (lb_params%medlyn_slope(ft) * term)**2 / &
            (gb_mol * vpd ))
       cquad = stomatal_intercept_btran*stomatal_intercept_btran + &
            (2.0*stomatal_intercept_btran + term * &
            (1.0 - lb_params%medlyn_slope(ft)* lb_params%medlyn_slope(ft) / vpd)) * term
       
       call QuadraticRoots(aquad, bquad, cquad, r1, r2)
       gs_mol = max(r1,r2)
    end if
      
  end subroutine StomatalCondMedlyn

  ! =======================================================================================

  subroutine StomatalCondBallBerry(a_gs,ft,veg_esat,ceair,stomatal_intercept_btran,leaf_co2_ppress,can_press, gb_mol, a_gs, gs_mol)
                

    ! Input
    integer, intent(in)  :: ft                       ! plant functional type index
    real(r8), intent(in) :: veg_esat                 ! saturation vapor pressure at veg_tempk (Pa)
    real(r8), intent(in) :: can_press                ! Air pressure NEAR the surface of the leaf (Pa)
    real(r8), intent(in) :: gb_mol                   ! leaf boundary layer conductance (umol /m**2/s)
    real(r8), intent(in) :: ceair                    ! vapor pressure of air, constrained (Pa)
    real(r8), intent(in) :: stomatal_intercept_btran ! water-stressed minimum stomatal conductance (umol H2O/m**2/s)
    real(r8), intent(in) :: leaf_co2_ppress          ! CO2 partial pressure at leaf surface (Pa)
    real(r8), intent(in) :: a_gs                     ! The assimilation (a) for calculating conductance (gs)
                                                     ! is either = to anet or agross
    
    ! Output
    real(r8) :: gs_mol            ! leaf stomatal conductance (umol H2O/m**2/s)

    ! locals
    real(r8) :: hs                   ! vapor pressure over saturation vapor pressure ratio
    real(r8) :: aquad,bquad,cquad    ! quadradic solve terms
    real(r8) :: r1,r2                ! quadradic solve roots

    if_pos_anet_cond: if (a_gs <= nearzero) then
       gs_mol = stomatal_intercept_btran
       return
    end if if_pos_anet_cond
                 
                    
    if(zero_bl_resist) then
                    
       hs = (ceair/ veg_esat)  
       gs_mol = lb_params%bb_slope(ft)*a_gs*hs/leaf_co2_ppress*can_press + stomatal_intercept_btran
       
    else
       aquad = leaf_co2_ppress
       bquad = leaf_co2_ppress*(gb_mol - stomatal_intercept_btran) - lb_params%bb_slope(ft) * a_gs * can_press
       cquad = -gb_mol*(leaf_co2_ppress * stomatal_intercept_btran + &
            lb_params%bb_slope(ft) * a_gs * can_press * ceair/ veg_esat )
                    
       call QuadraticRoots(aquad, bquad, cquad, r1, r2)
       gs_mol = max(r1,r2)
    end if
                    
                    
    return
  end subroutine StomatalCondBallBerry
  
  ! =======================================================================================

  subroutine LeafLayerPhotosynthesis(f_sun,         &  ! in
                                     parsun,        &  ! in
                                     parsha,        &  ! in
                                     leaf_area_sun,     &  ! in
                                     leaf_area_sha,     &  ! in
                                     ft,                &  ! in
                                     vcmax,             &  ! in
                                     jmax,              &  ! in
                                     co2_rcurve_islope, &  ! in
                                     veg_tempk,         &  ! in
                                     veg_esat,          &  ! in
                                     can_press,         &  ! in
                                     can_co2_ppress,    &  ! in
                                     can_o2_ppress,     &  ! in
                                     btran,             &  ! in
                                     stomatal_intercept_btran,  &  ! in
                                     cf,                &  ! in
                                     gb_mol,            &  ! in
                                     ceair,             &  ! in
                                     mm_kco2,           &  ! in
                                     mm_ko2,            &  ! in
                                     co2_cpoint,        &  ! in
                                     lmr,               &  ! in
                                     leaf_psi,          &  ! in
                                     rb,                &  ! in
                                     psn_out,           &  ! out
                                     rstoma_out,        &  ! out
                                     anet_av_out,       &  ! out
                                     c13disc_z)            ! out


  ! ------------------------------------------------------------------------------------
  ! This subroutine calculates photosynthesis and stomatal conductance within each leaf
  ! sublayer.
  ! A note on naming conventions: As this subroutine is called for every
  ! leaf-sublayer, many of the arguments are specific to that "leaf sub layer"
  ! (LSL), those variables are given a dimension tag "_lsl"
  ! Other arguments or variables may be indicative of scales broader than the LSL.
  ! ------------------------------------------------------------------------------------

  use EDParamsMod       , only : theta_cj_c3, theta_cj_c4


  ! Arguments
  ! ------------------------------------------------------------------------------------
  real(r8), intent(in) :: f_sun_lsl         !
  real(r8), intent(in) :: parsun_lsl        ! Absorbed PAR in sunlist leaves
  real(r8), intent(in) :: parsha_lsl        ! Absorved PAR in shaded leaves
  real(r8), intent(in) :: leaf_area_sun     ! leaf area in sunlit leaves [m2]
  real(r8), intent(in) :: leaf_area_sha     ! leaf area in shaded leaves [m2]
  integer,  intent(in) :: ft                ! (plant) Functional Type Index
  real(r8), intent(in) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
  real(r8), intent(in) :: jmax              ! maximum electron transport rate (umol electrons/m**2/s)
  real(r8), intent(in) :: co2_rcurve_islope ! initial slope of CO2 response curve (C4 plants)
  real(r8), intent(in) :: veg_tempk         ! vegetation temperature
  real(r8), intent(in) :: veg_esat          ! saturation vapor pressure at veg_tempk (Pa)

  ! Important Note on the following gas pressures.  This photosynthesis scheme will iteratively
  ! solve for the co2 partial pressure at the leaf surface (ie in the stomata). The reference
  ! point for these input values are NOT within that boundary layer that separates the stomata from
  ! the canopy air space.  The reference point for these is on the outside of that boundary
  ! layer.  This routine, which operates at the leaf scale, makes no assumptions about what the
  ! scale of the refernce is, it could be lower atmosphere, it could be within the canopy
  ! but most likely it is the closest value one can get to the edge of the leaf's boundary
  ! layer.  We use the convention "can_" because a reference point of within the canopy
  ! ia a best reasonable scenario of where we can get that information from.

  real(r8), intent(in) :: can_press       ! Air pressure NEAR the surface of the leaf (Pa)
  real(r8), intent(in) :: can_co2_ppress  ! Partial pressure of CO2 NEAR the leaf surface (Pa)
  real(r8), intent(in) :: can_o2_ppress   ! Partial pressure of O2 NEAR the leaf surface (Pa)
  real(r8), intent(in) :: btran           ! transpiration wetness factor (0 to 1)
  real(r8), intent(in) :: stomatal_intercept_btran !water-stressed minimum stomatal conductance (umol H2O/m**2/s)
  real(r8), intent(in) :: cf              ! s m**2/umol -> s/m (ideal gas conversion) [umol/m3]
  real(r8), intent(in) :: gb_mol          ! leaf boundary layer conductance (umol /m**2/s)
  real(r8), intent(in) :: ceair           ! vapor pressure of air, constrained (Pa)
  real(r8), intent(in) :: mm_kco2         ! Michaelis-Menten constant for CO2 (Pa)
  real(r8), intent(in) :: mm_ko2          ! Michaelis-Menten constant for O2 (Pa)
  real(r8), intent(in) :: co2_cpoint      ! CO2 compensation point (Pa)
  real(r8), intent(in) :: lmr             ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
  real(r8), intent(in) :: leaf_psi        ! Leaf water potential [MPa]
  real(r8), intent(in) :: rb              ! Boundary Layer resistance of leaf [s/m]
  
  real(r8), intent(out) :: psn_out        ! carbon assimilated in this leaf layer umolC/m2/s
  real(r8), intent(out) :: rstoma_out     ! stomatal resistance (1/gs_lsl) (s/m)
  real(r8), intent(out) :: anet_av_out    ! net leaf photosynthesis (umol CO2/m**2/s)
  ! averaged over sun and shade leaves.
  real(r8), intent(out) :: c13disc_z      ! carbon 13 in newly assimilated carbon
  real(r8), intent(out) :: test_out       ! for testing

 

  
  ! Locals
  ! ------------------------------------------------------------------------
  integer :: c3c4_path_index    ! Index for which photosynthetic pathway
  ! is active.  C4 = 0,  C3 = 1
  integer :: sunsha             ! Index for differentiating sun and shade
  real(r8) :: gstoma            ! Stomatal Conductance of this leaf layer (m/s)
  real(r8) :: agross            ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
  real(r8) :: anet              ! net leaf photosynthesis (umol CO2/m**2/s)
  real(r8) :: a_gs              ! The assimilation (a) for calculating conductance (gs)
                                ! is either = to anet or agross
  real(r8) :: je                ! electron transport rate (umol electrons/m**2/s)
  real(r8) :: qabs              ! PAR absorbed by PS II (umol photons/m**2/s)
  real(r8) :: aquad,bquad,cquad ! terms for quadratic equations
  real(r8) :: r1,r2             ! roots of quadratic equation
  real(r8) :: co2_inter_c       ! intercellular leaf CO2 (Pa)
  real(r8) :: co2_inter_c_old   ! intercellular leaf CO2 (Pa) (previous iteration)
  logical  :: loop_continue     ! Loop control variable
  integer  :: niter             ! iteration loop index
  real(r8) :: gs_mol            ! leaf stomatal conductance (umol H2O/m**2/s)
  real(r8) :: gs                ! leaf stomatal conductance (m/s)
  real(r8) :: hs                ! fractional humidity at leaf surface (dimensionless)
  real(r8) :: ac                ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
  real(r8) :: aj                ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
  real(r8) :: ap                ! product-limited (C3) or CO2-limited
  ! (C4) gross photosynthesis (umol CO2/m**2/s)
  real(r8) :: ai                ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
  real(r8) :: leaf_co2_ppress   ! CO2 partial pressure at leaf surface (Pa)
  real(r8) :: init_co2_inter_c  ! First guess intercellular co2 specific to C path
  real(r8) :: term                 ! intermediate variable in Medlyn stomatal conductance model
  real(r8) :: vpd                  ! water vapor deficit in Medlyn stomatal model (KPa)
  real(r8) :: leaf_area_sun_lsl
  real(r8) :: leaf_area_sha_lsl
  real(r8) :: par


  ! Parameters
  ! ------------------------------------------------------------------------
  ! Fraction of light absorbed by non-photosynthetic pigments
  real(r8),parameter :: fnps = 0.15_r8

  ! term accounting that two photons are needed to fully transport a single 
  ! electron in photosystem 2
  real(r8), parameter :: photon_to_e = 0.5_r8
  
  ! For plants with no leaves, a miniscule amount of conductance
  ! can happen through the stems, at a partial rate of cuticular conductance
  real(r8),parameter :: stem_cuticle_loss_frac = 0.1_r8

  ! empirical curvature parameter for electron transport rate
  real(r8),parameter :: theta_psii = 0.7_r8

  ! First guess on ratio between intercellular co2 and the atmosphere
  ! an iterator converges on actual
  real(r8),parameter :: init_a2l_co2_c3 = 0.7_r8
  real(r8),parameter :: init_a2l_co2_c4 = 0.4_r8

  ! quantum efficiency, used only for C4 (mol CO2 / mol photons) (index 0)
  real(r8),parameter,dimension(0:1) :: quant_eff = [0.05_r8,0.0_r8]
  integer, parameter :: max_iters = 5

  ! empirical curvature parameter for ap photosynthesis co-limitation
  real(r8),parameter :: theta_ip = 0.999_r8

  ! Set this to true if you want to assume zero resistance between
  ! leaf surface and boundary layer, which then assumes leaf surface
  ! humidity is the same as humidity on the other side of the boundary
  ! layer, essentially bypassing the use
  ! of Fick's law while calculating stomatal conductance.
  
  
  
  
  associate( bb_slope  => lb_params%bb_slope      ,& ! slope of BB relationship, unitless
       medlyn_slope=> lb_params%medlyn_slope          , & ! Slope for Medlyn stomatal conductance model method, the unit is KPa^0.5
       stomatal_intercept=> lb_params%stomatal_intercept )  !Unstressed minimum stomatal conductance, the unit is umol/m**2/s

  ! photosynthetic pathway: 0. = c4, 1. = c3
  c3c4_path_index = nint(lb_params%photopath(ft))

  if (c3c4_path_index == c3_path_index) then
     init_co2_inter_c = init_a2l_co2_c3 * can_co2_ppress
  else
     init_co2_inter_c = init_a2l_co2_c4 * can_co2_ppress
  end if

  ! Part III: Photosynthesis and Conductance
  ! ----------------------------------------------------------------------------------

  if_daytime: if ( parsun_lsl <= 0._r8 ) then  ! night time
 
     anet_av_out = -lmr
     psn_out     = 0._r8

     ! The cuticular conductance already factored in maximum resistance as a bound
     ! no need to re-bound it

     rstoma_out = cf/stomatal_intercept_btran

     c13disc_z = 0.0_r8    !carbon 13 discrimination in night time carbon flux, note value of 1.0 is used in CLM

  else ! day time (a little bit more complicated ...)

     ! Is there leaf area? - (NV can be larger than 0 with only stem area if deciduous)
      !leaf_area_sun_lsl = laisun_lsl*canopy_area_lsl
      !leaf_area_sha_lsl = laisha_lsl*canopy_area_lsl
     
     if_leafarea: if ( leaf_area_sun + leaf_area_sha > min_la_to_solve ) then
      
        !Loop aroun shaded and unshaded leaves
        psn_out     = 0._r8    ! psn is accumulated across sun and shaded leaves.
        rstoma_out  = 0._r8    ! 1/rs is accumulated across sun and shaded leaves.
        anet_av_out = 0._r8
        gstoma  = 0._r8

        do  sunsha = 1,2
           ! Electron transport rate for C3 plants.
           ! Convert par from W/m2 to umol photons/m**2/s
           ! Convert from units of par absorbed per unit ground area to par
           ! absorbed per unit leaf area
            if (sunsha == 1) then ! sunlit
               qabs = parsun_lsl
            else
               qabs = parsha_lsl             
            end if
            
            qabs = qabs*photon_to_e*(1.0_r8 - fnps)

           !convert the absorbed par into absorbed par per m2 of leaf,
           ! so it is consistant with the vcmax and lmr numbers.
           aquad = theta_psii
           bquad = -(qabs + jmax)
           cquad = qabs * jmax
           call QuadraticRoots(aquad, bquad, cquad, r1, r2)
           je = min(r1,r2)
             
           ! Initialize intercellular co2
           co2_inter_c = init_co2_inter_c

           niter = 0
           loop_continue = .true.
           iter_loop: do while(loop_continue)
              ! Increment iteration counter. Stop if too many iterations
              niter = niter + 1

              ! Save old co2_inter_c
              co2_inter_c_old = co2_inter_c

              ! Photosynthesis limitation rate calculations
              if (c3c4_path_index == c3_path_index)then

                 ! C3: Rubisco-limited photosynthesis
                 ac = vcmax * max(co2_inter_c-co2_cpoint, 0._r8) / &
                      (co2_inter_c+mm_kco2 * (1._r8+can_o2_ppress / mm_ko2 ))
                      
                 ! C3: RuBP-limited photosynthesis
                 aj = je * max(co2_inter_c-co2_cpoint, 0._r8) / &
                      (4._r8*co2_inter_c+8._r8*co2_cpoint)
                      
                    ! Gross photosynthesis smoothing calculations. Co-limit ac and aj.
                    aquad = theta_cj_c3
                    bquad = -(ac + aj)
                    cquad = ac * aj
                    call QuadraticRoots(aquad, bquad, cquad, r1, r2)
                    agross = min(r1,r2)
                    
              else

                 ! C4: Rubisco-limited photosynthesis
                 ac = vcmax


                 ! C4: RuBP-limited photosynthesis
                 if(sunsha == 1)then !sunlit
                    !guard against /0's in the night.
                    if(leaf_area_sun > min_la_to_solve) then
                       !par = ConvertPar(leaf_area_sun_lsl, parsun_lsl)
                        par = parsun_lsl
                       aj = quant_eff(c3c4_path_index)*par
                    else
                       aj = 0._r8
                    end if
                 else
                   ! par = ConvertPar(leaf_area_sha_lsl, parsha_lsl)
                     par = parsha_lsl
                    aj = quant_eff(c3c4_path_index)*par
                 end if
               

                 ! C4: PEP carboxylase-limited (CO2-limited)
                 ap = co2_rcurve_islope * max(co2_inter_c, 0._r8) / can_press
                 
                 ! Gross photosynthesis smoothing calculations. First co-limit ac and aj. Then co-limit ap

                 aquad = theta_cj_c4
                 bquad = -(ac + aj)
                 cquad = ac * aj
                 call QuadraticRoots(aquad, bquad, cquad, r1, r2)
                 ai = min(r1,r2)
                 
                 aquad = theta_ip
                 bquad = -(ai + ap)
                 cquad = ai * ap
                 call QuadraticRoots(aquad, bquad, cquad, r1, r2)
                 agross = min(r1,r2)
                 


              end if
              
              ! Calculate anet, only exit iteration with negative anet when
              ! using anet in calculating gs this is version B  
              anet = agross  - lmr
              
              if ( stomatal_assim_model == gross_assim_model ) then
                 if ( stomatal_model == medlyn_model ) then
                    write (fates_log(),*) 'Gross Assimilation conductance is incompatible with the Medlyn model'
                    call endrun(msg=errMsg(sourcefile, __LINE__))
                 end if
                 a_gs = agross
              else
                 if (anet < 0._r8) then
                    loop_continue = .false.
                 end if
                 a_gs = anet
              end if


              leaf_co2_ppress = can_co2_ppress - h2o_co2_bl_diffuse_ratio/gb_mol * a_gs * can_press 		   

              ! This does not seem necessary. THere would have to be a massive resistance
              ! between the two, no?
              if(use_mincap_leafco2) leaf_co2_ppress = max(leaf_co2_ppress,1.e-06_r8)
              
              ! A note about the use of the quadratic equations for calculating stomatal conductance
              ! ------------------------------------------------------------------------------------
              ! These two following models calculate the conductance between the intercellular leaf
              ! space and the leaf surface, not the canopy air space.  Transport between the leaf
              ! surface and the canopy air space is governed by the leaf boundary layer conductance.
              ! However, we need to estimate the properties at the surface of the leaf to solve for
              ! the stomatal conductance. We do this by using Fick's law (gradient resistance
              ! approximation of diffusion) to estimate the flux of water vapor across the
              ! leaf boundary layer, and balancing that with the flux across the stomata. It
              ! results in the following equation for leaf surface humidity:
              !
              ! e_s = (e_i g_s + e_c g_b)/(g_b + g_s)
              !
              ! The leaf surface humidity (e_s) becomes an expression of canopy humidity (e_c),
              ! intercellular humidity (e_i, which is the saturation humidity at leaf temperature),
              ! boundary layer conductance (g_b) (these are all known) and stomatal conductance
              ! (g_s) (this is still unknown).  This expression is substituted into the stomatal
              ! conductance equation. The resulting form of these equations becomes a quadratic.
              !
              ! For a detailed explanation, see the FATES technical note, section
              ! "1.11 Stomatal Conductance"
              !
              ! ------------------------------------------------------------------------------------
              
              if ( stomatal_model == medlyn_model ) then
                 call StomatalCondMedlyn(anet,ft,veg_esat,ceair,stomatal_intercept_btran,leaf_co2_ppress,can_press, gb_mol,gs_mol)
              else
                 call StomatalCondBallBerry()
              end if

         

            
                 



            



              
              ! Derive new estimate for co2_inter_c
              co2_inter_c = can_co2_ppress - anet * can_press * &
                   (h2o_co2_bl_diffuse_ratio*gs_mol+h2o_co2_stoma_diffuse_ratio*gb_mol) / (gb_mol*gs_mol)
               if (sunsha == 1) test_out = co2_inter_c
                   
               

              ! Check for co2_inter_c convergence. Delta co2_inter_c/pair = mol/mol.
              ! Multiply by 10**6 to convert to umol/mol (ppm). Exit iteration if
              ! convergence criteria of +/- 1 x 10**-6 ppm is met OR if at least ten
              ! iterations (niter=10) are completed

              if ((abs(co2_inter_c-co2_inter_c_old)/can_press*1.e06_r8 <=  2.e-06_r8) &
                   .or. niter >= max_iters) then
                 loop_continue = .false.
                 !if (.not. loop_continue) print *, veg_tempk - 273.15
              end if
           end do iter_loop
           

           ! End of co2_inter_c iteration.  Check for an < 0, in which case gs_mol = bbb
           ! And Final estimates for leaf_co2_ppress and co2_inter_c 
           ! (needed for early exit of co2_inter_c iteration when an < 0)	 
           if (anet < 0._r8) then
              gs_mol = stomatal_intercept_btran
           end if

           ! Final estimates for leaf_co2_ppress and co2_inter_c
           leaf_co2_ppress = can_co2_ppress - h2o_co2_bl_diffuse_ratio/gb_mol * anet * can_press
           leaf_co2_ppress = max(leaf_co2_ppress,1.e-06_r8)
           co2_inter_c = can_co2_ppress - anet * can_press * &
                (h2o_co2_bl_diffuse_ratio*gs_mol+h2o_co2_stoma_diffuse_ratio*gb_mol) / (gb_mol*gs_mol)

           ! Convert gs_mol (umol /m**2/s) to gs (m/s) and then to rs (s/m)
           gs = gs_mol / cf

           ! estimate carbon 13 discrimination in leaf level carbon
           ! flux Liang WEI and Hang ZHOU 2018, based on
           ! Ubierna and Farquhar, 2014 doi:10.1111/pce.12346, using the simplified model:
           ! $\Delta ^{13} C = \alpha_s + (b - \alpha_s) \cdot \frac{C_i}{C_a}$
           ! just hard code b and \alpha_s for now, might move to parameter set in future
           ! b = 27.0 alpha_s = 4.4
           ! TODO, not considering C4 or CAM right now, may need to address this
           ! note co2_inter_c is intracelluar CO2, not intercelluar
           c13disc_z = 4.4_r8 + (27.0_r8 - 4.4_r8) * &
                min (can_co2_ppress, max (co2_inter_c, 0._r8)) / can_co2_ppress

           ! Accumulate total photosynthesis umol/m2 ground/s-1.
           ! weight per unit sun and sha leaves.
           if(sunsha == 1)then !sunlit
              psn_out     = psn_out + agross * f_sun_lsl
              anet_av_out = anet_av_out + anet * f_sun_lsl
              gstoma  = gstoma + 1._r8/(min(1._r8/gs, rsmax0)) * f_sun_lsl
           else
              psn_out = psn_out + agross * (1.0_r8-f_sun_lsl)
              anet_av_out = anet_av_out + anet * (1.0_r8-f_sun_lsl)
              gstoma  = gstoma + &
                   1._r8/(min(1._r8/gs, rsmax0)) * (1.0_r8-f_sun_lsl)
           end if

           ! Make sure iterative solution is correct
           if (gs_mol < 0._r8) then
              write (fates_log(),*)'Negative stomatal conductance:'
              write (fates_log(),*)'gs_mol= ',gs_mol
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end if

        enddo !sunsha loop

        ! Stomatal resistance of the leaf-layer
        if ( (hlm_use_planthydro.eq.itrue .and. EDPftvarcon_inst%hydr_k_lwp(ft)>nearzero) ) then
           rstoma_out = LeafHumidityStomaResis(leaf_psi, veg_tempk, ceair, can_press, veg_esat, &
                                               rb, gstoma, ft)
        else
           rstoma_out = 1._r8/gstoma
        end if
           
        
     else

        ! No leaf area. This layer is present only because of stems.
        ! Net assimilation is zero, not negative because there are
        ! no leaves to even respire
        ! (leaves are off, or have reduced to 0)

        psn_out     = 0._r8
        anet_av_out = 0._r8

        rstoma_out  = min(rsmax0,cf/(stem_cuticle_loss_frac*stomatal_intercept(ft)))
        c13disc_z = 0.0_r8

     end if if_leafarea !is there leaf area?


   end if if_daytime    ! night or day


   end associate
   return
  end subroutine LeafLayerPhotosynthesis

  ! =======================================================================================

  function LeafHumidityStomaResis(leaf_psi, veg_tempk, ceair, can_press, veg_esat, &
       rb, gstoma, ft) result(rstoma_out)

    ! -------------------------------------------------------------------------------------
    ! This calculates inner leaf humidity as a function of mesophyll water potential 
    ! Adopted from  Vesala et al., 2017 https://www.frontiersin.org/articles/10.3389/fpls.2017.00054/full
    !
    ! Equation 1 in Vesala et al:
    ! lwp_star = wi/w0 = exp( k_lwp*leaf_psi*molar_mass_water/(rgas_J_k_mol * veg_tempk) )
    !
    ! Terms:
    ! leaf_psi: leaf water potential [MPa]
    ! k_lwp: inner leaf humidity scaling coefficient [-]
    ! rgas_J_K_mol: universal gas constant, [J/K/mol], 8.3144598
    ! molar_mass_water, molar mass of water, [g/mol]: 18.0
    !
    ! Unit conversions:
    ! 1 Pa = 1 N/m2 = 1 J/m3
    ! density of liquid water [kg/m3] = 1000
    ! 
    ! units of equation 1:  exp( [MPa]*[g/mol]/( [J/K/mol] * [K] ) )
    !                            [MJ/m3]*[g/mol]*[m3/kg]*[kg/g]*[J/MJ]  / ([J/mol])
    ! dimensionless:             [J/g]*[g/mol]/([J/mol])
    !
    ! Note: unit conversions drop out b/c [m3/kg]*[kg/g]*[J/MJ] = 1e-3*1.e-3*1e6 = 1.0
    !
    ! Junyan Ding 2021
    ! -------------------------------------------------------------------------------------

    ! Arguments
    real(r8) :: leaf_psi   ! Leaf water potential [MPa]
    real(r8) :: veg_tempk  ! Leaf temperature     [K]
    real(r8) :: ceair      ! vapor pressure of air, constrained [Pa]
    real(r8) :: can_press  ! Atmospheric pressure of canopy [Pa]
    real(r8) :: veg_esat   ! Saturated vapor pressure at veg surf [Pa]
    real(r8) :: rb         ! Leaf Boundary layer resistance [s/m]
    real(r8) :: gstoma     ! Stomatal Conductance of this leaf layer [m/s]
    integer  :: ft         ! Plant Functional Type
    real(r8) :: rstoma_out ! Total Stomatal resistance (stoma and BL) [s/m]

    ! Locals
    real(r8) :: k_lwp      ! Scaling coefficient for the ratio of leaf xylem
    ! water potential to mesophyll water potential
    real(r8) :: qs         ! Specific humidity [g/kg]
    real(r8) :: qsat       ! Saturation specific humidity  [g/kg]
    real(r8) :: qsat_adj   ! Adjusted saturation specific humidity  [g/kg]
    real(r8) :: lwp_star   ! leaf water potential scaling coefficient
    ! for inner leaf humidity, 0 means total dehydroted
    ! leaf, 1 means total saturated leaf

    ! Note: to disable this control, set k_lwp to zero, LWP_star will be 1
    k_lwp = EDPftvarcon_inst%hydr_k_lwp(ft)
    if (leaf_psi<0._r8) then
       lwp_star = exp(k_lwp*leaf_psi*molar_mass_water/(rgas_J_K_mol *veg_tempk))
    else 
       lwp_star = 1._r8
    end if

    ! compute specific humidity from vapor pressure
    ! q = molar_mass_ratio_vapdry*e/(can_press - (1-molar_mass_ratio_vapdry)*e) 
    ! source https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
    ! now adjust inner leaf humidity by LWP_star

    qs = molar_mass_ratio_vapdry * ceair / (can_press - (1._r8-molar_mass_ratio_vapdry) * ceair)
    qsat = molar_mass_ratio_vapdry * veg_esat / (can_press - (1._r8-molar_mass_ratio_vapdry) * veg_esat)
    qsat_adj = qsat*lwp_star

    ! Adjusting gs (compute a virtual gs) that will be passed to host model

    if ( qsat_adj < qs ) then

       ! if inner leaf vapor pressure is less then or equal to that at leaf surface
       ! then set stomata resistance to be very large to stop the transpiration or back flow of vapor
       rstoma_out = rsmax0

    else

       rstoma_out = (qsat-qs)*( 1/gstoma + rb)/(qsat_adj - qs)-rb

    end if

    if (rstoma_out < nearzero ) then
       write (fates_log(),*) 'qsat:', qsat, 'qs:', qs
       write (fates_log(),*) 'LWP :', leaf_psi
       write (fates_log(),*) 'ceair:', ceair, 'veg_esat:', veg_esat            
       write (fates_log(),*) 'rstoma_out:', rstoma_out, 'rb:', rb  
       write (fates_log(),*) 'LWP_star', lwp_star 
       call endrun(msg=errMsg(sourcefile, __LINE__))                  
    end if

  end function LeafHumidityStomaResis


  ! =====================================================================================

  subroutine ScaleLeafLayerFluxToCohort(nv,          & ! in   currentCohort%nv
       psn_llz,     & ! in   %psn_z(1:currentCohort%nv,ft,cl)
       lmr_llz,     & ! in   lmr_z(1:currentCohort%nv,ft,cl)
       rs_llz,      & ! in   rs_z(1:currentCohort%nv,ft,cl)
       elai_llz,    & ! in   %elai_profile(cl,ft,1:currentCohort%nv)
       c13disc_llz, & ! in   c13disc_z(cl, ft, 1:currentCohort%nv)
       c_area,      & ! in   currentCohort%c_area
       nplant,      & ! in   currentCohort%n
       rb,          & ! in   bc_in(s)%rb_pa(ifp)
       maintresp_reduction_factor, & ! in
       g_sb_laweight, & ! out  currentCohort%g_sb_laweight [m/s] [m2-leaf]
       gpp,         &   ! out  currentCohort%gpp_tstep
       rdark,       &   ! out  currentCohort%rdark
       c13disc_clm, &   ! out currentCohort%c13disc_clm
       cohort_eleaf_area ) ! out [m2]

    ! ------------------------------------------------------------------------------------
    ! This subroutine effectively integrates leaf carbon fluxes over the
    ! leaf layers to give cohort totals.
    ! Some arguments have the suffix "_llz".  This indicates that the vector
    ! is stratefied in the leaf-layer (ll)  dimension, and is a portion of the calling
    ! array which has the "_z" tag, thus "llz".
    ! ------------------------------------------------------------------------------------

    use FatesConstantsMod, only : umolC_to_kgC

    ! Arguments
    integer, intent(in)  :: nv               ! number of active leaf layers
    real(r8), intent(in) :: psn_llz(nv)      ! layer photosynthesis rate (GPP) [umolC/m2leaf/s]
    real(r8), intent(in) :: lmr_llz(nv)      ! layer dark respiration rate [umolC/m2leaf/s]
    real(r8), intent(in) :: rs_llz(nv)       ! leaf layer stomatal resistance [s/m]
    real(r8), intent(in) :: elai_llz(nv)     ! exposed LAI per layer [m2 leaf/ m2 pft footprint]
    real(r8), intent(in) :: c13disc_llz(nv)  ! leaf layer c13 discrimination, weighted mean
    real(r8), intent(in) :: c_area           ! crown area m2/m2
    real(r8), intent(in) :: nplant           ! indiv/m2
    real(r8), intent(in) :: rb               ! leaf boundary layer resistance (s/m)
    real(r8), intent(in) :: maintresp_reduction_factor  ! factor by which to reduce maintenance respiration
    real(r8), intent(out) :: g_sb_laweight     ! Combined conductance (stomatal + boundary layer) for the cohort
    ! weighted by leaf area [m/s]*[m2]
    real(r8), intent(out) :: gpp               ! GPP (kgC/indiv/s)
    real(r8), intent(out) :: rdark             ! Dark Leaf Respiration (kgC/indiv/s)
    real(r8), intent(out) :: cohort_eleaf_area ! Effective leaf area of the cohort [m2]
    real(r8), intent(out) :: c13disc_clm       ! unpacked Cohort level c13 discrimination
    real(r8)              :: sum_weight        ! sum of weight for unpacking d13c flux (c13disc_z) from
    ! (canopy_layer, pft, leaf_layer) matrix to cohort (c13disc_clm)

    ! GPP IN THIS SUBROUTINE IS A RATE. THE CALLING ARGUMENT IS GPP_TSTEP. AFTER THIS
    ! CALL THE RATE WILL BE MULTIPLIED BY THE INTERVAL TO GIVE THE INTEGRATED QUANT.

    ! Locals
    integer  :: il                       ! leaf layer index
    real(r8) :: cohort_layer_eleaf_area  ! the effective leaf area of the cohort's current layer [m2]

    cohort_eleaf_area = 0.0_r8
    g_sb_laweight     = 0.0_r8
    gpp               = 0.0_r8
    rdark             = 0.0_r8

    do il = 1, nv        ! Loop over the leaf layers this cohort participates in


       ! Cohort's total effective leaf area in this layer [m2]
       ! leaf area index of the layer [m2/m2 ground] * [m2 ground]
       ! elai_llz is the LAI for the whole PFT. Multiplying this by the ground
       ! area this cohort contributes, give the cohort's portion of the leaf
       ! area in this layer
       cohort_layer_eleaf_area = elai_llz(il) * c_area

       ! Increment the cohort's total effective leaf area [m2]
       cohort_eleaf_area       = cohort_eleaf_area + cohort_layer_eleaf_area

       ! Leaf conductance (stomatal and boundary layer)
       ! This should be the weighted average over the leaf surfaces.
       ! Since this is relevant to the stomata, its weighting should be based
       ! on total leaf area, and not really footprint area
       ! [m/s] * [m2 cohort's leaf layer]
       g_sb_laweight = g_sb_laweight + 1.0_r8/(rs_llz(il)+rb) * cohort_layer_eleaf_area

       ! GPP    [umolC/m2leaf/s] * [m2 leaf ] -> [umolC/s]
       gpp = gpp + psn_llz(il) * cohort_layer_eleaf_area

       ! Dark respiration
       ! [umolC/m2leaf/s] * [m2 leaf] 
       rdark = rdark + lmr_llz(il) * cohort_layer_eleaf_area

    end do



    if (nv > 1) then
       ! cohort%c13disc_clm as weighted mean of d13c flux at all related leave layers
       sum_weight = sum(psn_llz(1:nv-1) * elai_llz(1:nv-1))
       if (sum_weight .eq. 0.0_r8) then
          c13disc_clm = 0.0
       else
          c13disc_clm = sum(c13disc_llz(1:nv-1) * psn_llz(1:nv-1) * elai_llz(1:nv-1)) / sum_weight
       end if

    end if


    ! -----------------------------------------------------------------------------------
    ! We DO NOT normalize g_sb_laweight.
    ! The units that we are passing back are [m/s] * [m2 effective leaf]
    ! We will add these up over the whole patch, and then normalized
    ! by the patch's total leaf area in the calling routine
    ! -----------------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------
    ! Convert dark respiration and GPP from [umol/s] to [kgC/plant/s]
    ! Also, apply the maintenance respiration reduction factor
    ! -----------------------------------------------------------------------------------

    rdark     = rdark * umolC_to_kgC * maintresp_reduction_factor / nplant
    gpp       = gpp * umolC_to_kgC / nplant

    if ( debug ) then
       write(fates_log(),*) 'EDPhoto 816 ', gpp
       write(fates_log(),*) 'EDPhoto 817 ', psn_llz(1:nv)
       write(fates_log(),*) 'EDPhoto 820 ', nv
       write(fates_log(),*) 'EDPhoto 821 ', elai_llz(1:nv)
       write(fates_log(),*) 'EDPhoto 843 ', rdark
       write(fates_log(),*) 'EDPhoto 873 ', nv
       write(fates_log(),*) 'EDPhoto 874 ', cohort_eleaf_area
    endif

    return
  end subroutine ScaleLeafLayerFluxToCohort

  ! =====================================================================================

  function ft1_f(tl, ha) result(ans)
    !
    !!DESCRIPTION:
    ! photosynthesis temperature response
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    !!USES
    use FatesConstantsMod, only : rgas => rgas_J_K_kmol

    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temperature function (K)
    real(r8), intent(in) :: ha  ! activation energy in photosynthesis temperature function (J/mol)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = exp( ha / (rgas*1.e-3_r8*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )

    return
  end function ft1_f

  ! =====================================================================================

  function fth_f(tl,hd,se,scaleFactor) result(ans)
    !
    !!DESCRIPTION:
    !photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    use FatesConstantsMod, only : rgas => rgas_J_K_kmol

    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temp function (K)
    real(r8), intent(in) :: hd  ! deactivation energy in photosynthesis temp function (J/mol)
    real(r8), intent(in) :: se  ! entropy term in photosynthesis temp function (J/mol/K)
    real(r8), intent(in) :: scaleFactor  ! scaling factor for high temp inhibition (25 C = 1.0)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = scaleFactor / ( 1._r8 + exp( (-hd+se*tl) / (rgas*1.e-3_r8*tl) ) )

    return
  end function fth_f

  ! =====================================================================================

  function fth25_f(hd,se)result(ans)
    !
    !!DESCRIPTION:
    ! scaling factor for photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY:
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    !!USES

    use FatesConstantsMod, only : rgas => rgas_J_K_kmol

    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: hd    ! deactivation energy in photosynthesis temp function (J/mol)
    real(r8), intent(in) :: se    ! entropy term in photosynthesis temp function (J/mol/K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = 1._r8 + exp( (-hd+se*(tfrz+25._r8)) / (rgas*1.e-3_r8*(tfrz+25._r8)) )

    return
  end function fth25_f

  ! =====================================================================================

  subroutine UpdateCanopyNCanNRadPresent(currentPatch)

    ! ---------------------------------------------------------------------------------
    ! This subroutine calculates two patch level quanities:
    ! currentPatch%ncan   and
    ! currentPatch%canopy_mask
    !
    ! currentPatch%ncan(:,:) is a two dimensional array that indicates
    ! the total number of leaf layers (including those that are not exposed to light)
    ! in each canopy layer and for each functional type.
    !
    ! currentPatch%nrad(:,:) is a two dimensional array that indicates
    ! the total number of EXPOSED leaf layers, but for all intents and purposes
    ! in the photosynthesis routine, this appears to be the same as %ncan...
    !
    ! currentPatch%canopy_mask(:,:) has the same dimensions, is binary, and
    ! indicates whether or not leaf layers are present (by evaluating the canopy area
    ! profile).
    ! ---------------------------------------------------------------------------------

    ! Arguments
    type(fates_patch_type), target :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort

    ! Locals
    integer :: cl  ! Canopy Layer Index
    integer :: ft  ! Function Type Index
    integer :: iv  ! index of the exposed leaf layer for each canopy layer and pft

    ! Loop through the cohorts in this patch, associate each cohort with a layer and PFT
    ! and use the cohort's memory of how many layer's it takes up to assign the maximum
    ! of the layer/pft index it is in
    ! ---------------------------------------------------------------------------------

    currentPatch%ncan(:,:) = 0
    ! redo the canopy structure algorithm to get round a
    ! bug that is happening for site 125, FT13.
    currentCohort => currentPatch%tallest
    do while(associated(currentCohort))

       currentPatch%ncan(currentCohort%canopy_layer,currentCohort%pft) = &
            max(currentPatch%ncan(currentCohort%canopy_layer,currentCohort%pft), &
            currentCohort%NV)

       currentCohort => currentCohort%shorter

    enddo !cohort

    ! NRAD = NCAN ...
    currentPatch%nrad = currentPatch%ncan

    ! Now loop through and identify which layer and pft combo has scattering elements
    do cl = 1,nclmax
       do ft = 1,numpft
          currentPatch%canopy_mask(cl,ft) = 0
          do iv = 1, currentPatch%nrad(cl,ft);
             if(currentPatch%canopy_area_profile(cl,ft,iv) > 0._r8)then
                currentPatch%canopy_mask(cl,ft) = 1
             end if
          end do !iv
       enddo !ft
    enddo !cl

    return
  end subroutine UpdateCanopyNCanNRadPresent

  ! ====================================================================================

  subroutine GetCanopyGasParameters(can_press, &
       can_o2_partialpress, &
       veg_tempk, &
       air_tempk, &
       air_vpress, &
       veg_esat,   &
       rb,        &
       mm_kco2,   &
       mm_ko2,    &
       co2_cpoint, &
       cf,         &
       gb_mol, &
       ceair)

    ! ---------------------------------------------------------------------------------
    ! This subroutine calculates the specific Michaelis Menten Parameters (pa) for CO2
    ! and O2, as well as the CO2 compentation point.
    ! ---------------------------------------------------------------------------------

    use FatesConstantsMod, only: umol_per_mol
    use FatesConstantsMod, only: mmol_per_mol
    use FatesConstantsMod, only: umol_per_kmol
    use FatesConstantsMod, only : rgas => rgas_J_K_kmol

    ! Arguments
    real(r8), intent(in) :: can_press           ! Air pressure within the canopy (Pa)
    real(r8), intent(in) :: can_o2_partialpress ! Partial press of o2 in the canopy (Pa)
    real(r8), intent(in) :: veg_tempk           ! The temperature of the vegetation (K)
    real(r8), intent(in) :: air_tempk           ! Temperature of canopy air (K)
    real(r8), intent(in) :: air_vpress          ! Vapor pressure of canopy air (Pa)
    real(r8), intent(in) :: veg_esat            ! Saturated vapor pressure at veg surf (Pa)
    real(r8), intent(in) :: rb                  ! Leaf Boundary layer resistance (s/m)

    real(r8), intent(out) :: mm_kco2       ! Michaelis-Menten constant for CO2 (Pa)
    real(r8), intent(out) :: mm_ko2        !  Michaelis-Menten constant for O2 (Pa)
    real(r8), intent(out) :: co2_cpoint    !  CO2 compensation point (Pa)
    real(r8), intent(out) :: cf            ! conversion factor between molar form and velocity form
    ! of conductance and resistance: [umol/m3]
    real(r8), intent(out) :: gb_mol        ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(r8), intent(out) :: ceair         ! vapor pressure of air, constrained (Pa)

    ! Locals
    real(r8) :: kc25                ! Michaelis-Menten constant for CO2 at 25C (Pa)
    real(r8) :: ko25                ! Michaelis-Menten constant for O2 at 25C (Pa)
    real(r8) :: sco                 ! relative specificity of rubisco
    real(r8) :: cp25                ! CO2 compensation point at 25C (Pa)

    ! ---------------------------------------------------------------------------------
    ! Intensive values (per mol of air)
    ! kc, ko, currentPatch, from: Bernacchi et al (2001)
    ! Plant, Cell and Environment 24:253-259
    ! ---------------------------------------------------------------------------------

    real(r8), parameter :: mm_kc25_umol_per_mol       = 404.9_r8
    real(r8), parameter :: mm_ko25_mmol_per_mol       = 278.4_r8
    real(r8), parameter :: co2_cpoint_umol_per_mol    = 42.75_r8

    ! Activation energy, from:
    ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
    ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430

    real(r8), parameter :: kcha    = 79430._r8  ! activation energy for kc (J/mol)
    real(r8), parameter :: koha    = 36380._r8  ! activation energy for ko (J/mol)
    real(r8), parameter :: cpha    = 37830._r8  ! activation energy for cp (J/mol)


    ! Derive sco from currentPatch and O2 using present-day O2 (0.209 mol/mol) and re-calculate
    ! currentPatch to account for variation in O2 using currentPatch = 0.5 O2 / sco

    ! FIXME (RGK 11-30-2016 THere are more constants here, but I don't have enough information
    ! about what they are or do, so I can't give them more descriptive names. Someone please
    ! fill this in when possible)

    kc25 = ( mm_kc25_umol_per_mol / umol_per_mol ) * can_press
    ko25 = ( mm_ko25_mmol_per_mol / mmol_per_mol ) * can_press
    sco  = 0.5_r8 * 0.209_r8 / (co2_cpoint_umol_per_mol / umol_per_mol )
    cp25 = 0.5_r8 * can_o2_partialpress / sco

    if( veg_tempk.gt.150_r8 .and. veg_tempk.lt.350_r8 )then
       mm_kco2       = kc25 * ft1_f(veg_tempk, kcha)
       mm_ko2         = ko25 * ft1_f(veg_tempk, koha)
       co2_cpoint     = cp25 * ft1_f(veg_tempk, cpha)
    else
       mm_kco2    = 1.0_r8
       mm_ko2     = 1.0_r8
       co2_cpoint = 1.0_r8
    end if

    ! ---------------------------------------------------------------------------------
    !
    ! cf is the conversion factor between molar form and velocity form
    ! of conductance and resistance: [umol/m3]
    !
    ! i.e.
    ! [m/s] * [umol/m3] -> [umol/m2/s]
    !
    ! Breakdown of the conversion factor: [ umol / m3 ]
    !
    ! Rgas [J /K /kmol]
    ! Air Potential Temperature [ K ]
    ! Canopy Pressure      [ Pa ]
    ! conversion: umol/kmol =  1e9
    !
    ! [ Pa * K * kmol umol/kmol  /  J K ] = [ Pa * umol / J ]
    ! since: 1 Pa = 1 N / m2
    ! [ Pa * umol / J ] = [ N * umol / J m2 ]
    ! since: 1 J = 1 N * m
    ! [ N * umol / J m2 ] = [ N * umol / N m3 ]
    ! [ umol / m3 ]
    !
    ! --------------------------------------------------------------------------------

    cf = can_press/(rgas * air_tempk )*umol_per_kmol
    gb_mol = (1._r8/ rb) * cf

    ! Constrain eair >= 0.05*esat_tv so that solution does not blow up. This ensures
    ! that hs does not go to zero. Also eair <= veg_esat so that hs <= 1
    ceair = min( max(air_vpress, 0.05_r8*veg_esat ),veg_esat )



    return
  end subroutine GetCanopyGasParameters

  ! ====================================================================================

  subroutine LeafLayerMaintenanceRespiration_Ryan_1991(lnc_top, &
       nscaler,       &
       ft,            &
       veg_tempk,     &
       lmr)

    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod, only : umolC_to_kgC
    use FatesConstantsMod, only : g_per_kg
    use EDPftvarcon      , only : EDPftvarcon_inst

    ! -----------------------------------------------------------------------
    ! Base maintenance respiration rate for plant tissues maintresp_leaf_ryan1991_baserate
    ! M. Ryan, 1991. Effects of climate change on plant respiration.
    ! Ecological Applications, 1(2), 157-167.
    ! Original expression is br = 0.0106 molC/(molN h)
    ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
    ! Which is the default value of maintresp_nonleaf_baserate

    ! Arguments
    real(r8), intent(in)  :: lnc_top      ! Leaf nitrogen content per unit area at canopy top [gN/m2]
    real(r8), intent(in)  :: nscaler      ! Scale for leaf nitrogen profile
    integer,  intent(in)  :: ft           ! (plant) Functional Type Index
    real(r8), intent(in)  :: veg_tempk    ! vegetation temperature
    real(r8), intent(out) :: lmr          ! Leaf Maintenance Respiration  (umol CO2/m**2/s)

    ! Locals
    real(r8) :: lmr25   ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25top  ! canopy top leaf maint resp rate at 25C for this pft (umol CO2/m**2/s)
    integer :: c3c4_path_index    ! Index for which photosynthetic pathway

    ! Parameter
    real(r8), parameter :: lmrha = 46390._r8    ! activation energy for lmr (J/mol)
    real(r8), parameter :: lmrhd = 150650._r8   ! deactivation energy for lmr (J/mol)
    real(r8), parameter :: lmrse = 490._r8      ! entropy term for lmr (J/mol/K)
    real(r8), parameter :: lmrc = 1.15912391_r8 ! scaling factor for high
    ! temperature inhibition (25 C = 1.0)

    lmr25top = EDPftvarcon_inst%maintresp_leaf_ryan1991_baserate(ft) * (1.5_r8 ** ((25._r8 - 20._r8)/10._r8))
    lmr25top = lmr25top * lnc_top / (umolC_to_kgC * g_per_kg)


    ! Part I: Leaf Maintenance respiration: umol CO2 / m**2 [leaf] / s
    ! ----------------------------------------------------------------------------------
    lmr25 = lmr25top * nscaler

    ! photosynthetic pathway: 0. = c4, 1. = c3
    c3c4_path_index = nint(lb_params%photopath(ft))

    if (c3c4_path_index == c3_path_index) then
       ! temperature sensitivity of C3 plants
       lmr = lmr25 * ft1_f(veg_tempk, lmrha) * &
            fth_f(veg_tempk, lmrhd, lmrse, lmrc)
    else
       ! temperature sensitivity of C4 plants
       lmr = lmr25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
       lmr = lmr / (1._r8 + exp( 1.3_r8*(veg_tempk-(tfrz+55._r8)) ))
    endif

    ! Any hydrodynamic limitations could go here, currently none
    ! lmr = lmr * (nothing)

  end subroutine LeafLayerMaintenanceRespiration_Ryan_1991

  ! ====================================================================================   

  subroutine LeafLayerMaintenanceRespiration_Atkin_etal_2017(lnc_top, &
       cumulative_lai, &
       vcmax25top,     &
       ft,             &
       veg_tempk,      &
       tgrowth,        &
       lmr)

    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod, only : umolC_to_kgC
    use FatesConstantsMod, only : g_per_kg
    use FatesConstantsMod, only : lmr_b
    use FatesConstantsMod, only : lmr_c
    use FatesConstantsMod, only : lmr_TrefC
    use FatesConstantsMod, only : lmr_r_1
    use FatesConstantsMod, only : lmr_r_2
    use EDPftvarcon      , only : EDPftvarcon_inst

    ! Arguments
    real(r8), intent(in)  :: lnc_top          ! Leaf nitrogen content per unit area at canopy top [gN/m2]
    integer,  intent(in)  :: ft               ! (plant) Functional Type Index
    real(r8), intent(in)  :: vcmax25top       ! top of canopy vcmax
    real(r8), intent(in)  :: cumulative_lai   ! cumulative lai above the current leaf layer
    real(r8), intent(in)  :: veg_tempk        ! vegetation temperature  (degrees K)
    real(r8), intent(in)  :: tgrowth          ! lagged vegetation temperature averaged over acclimation timescale (degrees K)
    real(r8), intent(out) :: lmr              ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
    
    ! Locals
    real(r8) :: lmr25   ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: r_0     ! base respiration rate, PFT-dependent (umol CO2/m**2/s)
    real(r8) :: r_t_ref ! acclimated ref respiration rate (umol CO2/m**2/s)
    real(r8) :: lmr25top  ! canopy top leaf maint resp rate at 25C for this pft (umol CO2/m**2/s)

    real(r8) :: rdark_scaler ! negative exponential scaling of rdark
    real(r8) :: kn           ! decay coefficient
   
    ! parameter values of r_0 as listed in Atkin et al 2017: (umol CO2/m**2/s) 
    ! Broad-leaved trees  1.7560
    ! Needle-leaf trees   1.4995
    ! Shrubs              2.0749
    ! C3 herbs/grasses    2.1956
    ! In the absence of better information, we use the same value for C4 grasses as C3 grasses.

    ! r_0 currently put into the EDPftvarcon_inst%dev_arbitrary_pft
    ! all figs in Atkin et al 2017 stop at zero Celsius so we will assume acclimation is fixed below that
    r_0 = EDPftvarcon_inst%maintresp_leaf_atkin2017_baserate(ft)

    ! This code uses the relationship between leaf N and respiration from Atkin et al 
    ! for the top of the canopy, but then scales through the canopy based on a rdark_scaler.
    ! To assume proportionality with N through the canopy following Lloyd et al. 2010, use the
    ! default parameter value of 2.43, which results in the scaling of photosynthesis and respiration
    ! being proportional through the canopy. To have a steeper decrease in respiration than photosynthesis
    ! this number can be smaller. There is some observational evidence for this being the case
    ! in Lamour et al. 2023. 

    kn = decay_coeff_vcmax(vcmax25top, &
                           EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff1(ft), &
                           EDPftvarcon_inst%maintresp_leaf_vert_scaler_coeff2(ft))

    rdark_scaler = exp(-kn * cumulative_lai)
    
    r_t_ref = max(0._r8, rdark_scaler * (r_0 + lmr_r_1 * lnc_top + lmr_r_2 * max(0._r8, (tgrowth - tfrz) )) )

    if (r_t_ref .eq. 0._r8) then
       warn_msg = 'Rdark is negative at this temperature and is capped at 0. tgrowth (C): '//trim(N2S(tgrowth-tfrz))//' pft: '//trim(I2S(ft))
       call FatesWarn(warn_msg,index=4)            
    end if

    lmr = r_t_ref * exp(lmr_b * (veg_tempk - tfrz - lmr_TrefC) + lmr_c * &
         ((veg_tempk-tfrz)**2 - lmr_TrefC**2))

  end subroutine LeafLayerMaintenanceRespiration_Atkin_etal_2017

  ! ====================================================================================

  subroutine LeafLayerBiophysicalRates( parsun_per_la, &
       ft,            &
       vcmax25top_ft, &
       jmax25top_ft, &
       co2_rcurve_islope25top_ft, &
       nscaler,    &
       veg_tempk,      &
       dayl_factor, &
       t_growth,   &
       t_home,     &
       btran, &
       vcmax, &
       jmax, &
       co2_rcurve_islope )

    ! ---------------------------------------------------------------------------------
    ! This subroutine calculates the localized rates of several key photosynthesis
    ! rates.  By localized, we mean specific to the plant type and leaf layer,
    ! which factors in leaf physiology, as well as environmental effects.
    ! This procedure should be called prior to iterative solvers, and should
    ! have pre-calculated the reference rates for the pfts before this.
    !
    ! The output biophysical rates are:
    ! vcmax: maximum rate of carboxilation,
    ! jmax: maximum electron transport rate,
    ! co2_rcurve_islope: initial slope of CO2 response curve (C4 plants)
    ! ---------------------------------------------------------------------------------

    use EDPftvarcon         , only : EDPftvarcon_inst

    ! Arguments
    ! ------------------------------------------------------------------------------

    real(r8), intent(in) :: parsun_per_la   ! PAR absorbed per sunlit leaves for this layer
    integer,  intent(in) :: ft              ! (plant) Functional Type Index
    real(r8), intent(in) :: nscaler         ! Scale for leaf nitrogen profile
    real(r8), intent(in) :: vcmax25top_ft   ! canopy top maximum rate of carboxylation at 25C
    ! for this pft (umol CO2/m**2/s)
    real(r8), intent(in) :: jmax25top_ft    ! canopy top maximum electron transport rate at 25C
    ! for this pft (umol electrons/m**2/s)
    real(r8), intent(in) :: co2_rcurve_islope25top_ft ! initial slope of CO2 response curve
    ! (C4 plants) at 25C, canopy top, this pft
    real(r8), intent(in) :: veg_tempk           ! vegetation temperature
    real(r8), intent(in) :: dayl_factor         ! daylength scaling factor (0-1)
    real(r8), intent(in) :: t_growth            ! T_growth (short-term running mean temperature) (K)
    real(r8), intent(in) :: t_home              ! T_home (long-term running mean temperature) (K)
    real(r8), intent(in) :: btran           ! transpiration wetness factor (0 to 1)

    real(r8), intent(out) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8), intent(out) :: jmax              ! maximum electron transport rate
    ! (umol electrons/m**2/s)
    real(r8), intent(out) :: co2_rcurve_islope ! initial slope of CO2 response curve (C4 plants)

    ! Locals
    ! -------------------------------------------------------------------------------
    real(r8) :: vcmax25             ! leaf layer: maximum rate of carboxylation at 25C
    ! (umol CO2/m**2/s)
    real(r8) :: jmax25              ! leaf layer: maximum electron transport rate at 25C
    ! (umol electrons/m**2/s)
    real(r8) :: co2_rcurve_islope25 ! leaf layer: Initial slope of CO2 response curve
    ! (C4 plants) at 25C
    integer :: c3c4_path_index      ! Index for which photosynthetic pathway
    real(r8) :: dayl_factor_local   ! Local version of daylength factor

    ! Parameters
    ! ---------------------------------------------------------------------------------
    real(r8) :: vcmaxha        ! activation energy for vcmax (J/mol)
    real(r8) :: jmaxha         ! activation energy for jmax (J/mol)
    real(r8) :: vcmaxhd        ! deactivation energy for vcmax (J/mol)
    real(r8) :: jmaxhd         ! deactivation energy for jmax (J/mol)
    real(r8) :: vcmaxse        ! entropy term for vcmax (J/mol/K)
    real(r8) :: jmaxse         ! entropy term for jmax (J/mol/K)
    real(r8) :: t_growth_celsius ! average growing temperature
    real(r8) :: t_home_celsius   ! average home temperature
    real(r8) :: jvr            ! ratio of Jmax25 / Vcmax25
    real(r8) :: vcmaxc         ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(r8) :: jmaxc          ! scaling factor for high temperature inhibition (25 C = 1.0)

    select case(photo_tempsens_model)
    case (photosynth_acclim_model_none) !No temperature acclimation
       vcmaxha = EDPftvarcon_inst%vcmaxha(FT)
       jmaxha  = EDPftvarcon_inst%jmaxha(FT)
       vcmaxhd = EDPftvarcon_inst%vcmaxhd(FT)
       jmaxhd  = EDPftvarcon_inst%jmaxhd(FT)
       vcmaxse = EDPftvarcon_inst%vcmaxse(FT)
       jmaxse  = EDPftvarcon_inst%jmaxse(FT)
    case (photosynth_acclim_model_kumarathunge_etal_2019) !Kumarathunge et al. temperature acclimation, Thome=30-year running mean
       t_growth_celsius = t_growth-tfrz
       t_home_celsius = t_home-tfrz
       vcmaxha = (42.6_r8 + (1.14_r8*t_growth_celsius))*1e3_r8 !J/mol
       jmaxha = 40.71_r8*1e3_r8 !J/mol
       vcmaxhd = 200._r8*1e3_r8 !J/mol
       jmaxhd = 200._r8*1e3_r8 !J/mol
       vcmaxse = (645.13_r8 - (0.38_r8*t_growth_celsius))
       jmaxse = 658.77_r8 - (0.84_r8*t_home_celsius) - 0.52_r8*(t_growth_celsius-t_home_celsius)
       jvr = 2.56_r8 - (0.0375_r8*t_home_celsius)-(0.0202_r8*(t_growth_celsius-t_home_celsius))
    case default
       write (fates_log(),*)'error, incorrect leaf photosynthesis temperature acclimation model specified'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end select

    vcmaxc = fth25_f(vcmaxhd, vcmaxse)
    jmaxc  = fth25_f(jmaxhd, jmaxse)

    if(parsun_per_la <= 0._r8) then
       vcmax             = 0._r8
       jmax              = 0._r8
       co2_rcurve_islope = 0._r8
    else                                     ! day time

       ! update the daylength factor local variable if the switch is on
       if ( dayl_switch == itrue ) then
          dayl_factor_local = dayl_factor
       else
          dayl_factor_local = 1.0_r8
       endif

       ! Vcmax25top was already calculated to derive the nscaler function
       vcmax25 = vcmax25top_ft * nscaler * dayl_factor_local
       select case(photo_tempsens_model)
       case (photosynth_acclim_model_none)
          jmax25  = jmax25top_ft * nscaler * dayl_factor_local
       case (photosynth_acclim_model_kumarathunge_etal_2019) 
          jmax25 = vcmax25*jvr
       case default
          write (fates_log(),*)'error, incorrect leaf photosynthesis temperature acclimation model specified'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select

       co2_rcurve_islope25 = co2_rcurve_islope25top_ft * nscaler

       ! Adjust for temperature
       ! photosynthetic pathway: 0. = c4, 1. = c3
       c3c4_path_index = nint(lb_params%photopath(ft))

       if (c3c4_path_index == c3_path_index) then
          vcmax = vcmax25 * ft1_f(veg_tempk, vcmaxha) * fth_f(veg_tempk, vcmaxhd, vcmaxse, vcmaxc)
       else
          vcmax = vcmax25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
          vcmax = vcmax / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-veg_tempk ) ))
          vcmax = vcmax / (1._r8 + exp( 0.3_r8*(veg_tempk-(tfrz+40._r8)) ))
       end if

       jmax  = jmax25 * ft1_f(veg_tempk, jmaxha) * fth_f(veg_tempk, jmaxhd, jmaxse, jmaxc)

       !q10 response of product limited psn.
       co2_rcurve_islope = co2_rcurve_islope25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
    end if

    ! Adjust for water limitations
    vcmax = vcmax * btran

    return

  end subroutine LeafLayerBiophysicalRates

  subroutine lowstorage_maintresp_reduction(frac, pft, maintresp_reduction_factor)

    ! This subroutine reduces maintenance respiration rates when storage pool is low.  The premise
    ! of this is that mortality of plants increases when storage is low because they are not able
    ! to repair tissues, generate defense compounds, etc.  This reduction is reflected in a reduced
    ! maintenance demand.  The output of this function takes the form of a curve between 0 and 1,
    ! and the curvature of the function is determined by a parameter.

    ! Uses
    use EDPftvarcon         , only : EDPftvarcon_inst

    ! Arguments
    ! ------------------------------------------------------------------------------
    real(r8), intent(in) :: frac      ! ratio of storage to target leaf biomass
    integer,  intent(in) :: pft       ! what pft is this cohort?
    real(r8), intent(out) :: maintresp_reduction_factor  ! the factor by which to reduce maintenance respiration

    ! --------------------------------------------------------------------------------
    ! Parameters are at the PFT level:
    ! fates_maintresp_reduction_curvature controls the curvature of this.
    ! If this parameter is zero, then there is no reduction until the plant dies at storage = 0.
    ! If this parameter is one, then there is a linear reduction in respiration below the storage point.
    ! Intermediate values will give some (concave-downwards) curvature.
    !
    ! maintresp_reduction_intercept controls the maximum amount of throttling.
    ! zero means no throttling at any point, so it turns this mechanism off completely and so
    ! allows an entire cohort to die via negative carbon-induced termination mortality.
    ! one means complete throttling, so no maintenance respiration at all, when out of carbon.
    ! ---------------------------------------------------------------------------------

    if( frac .lt. 1._r8 )then
       if ( abs(EDPftvarcon_inst%maintresp_reduction_curvature(pft)-1._r8) > nearzero ) then
          maintresp_reduction_factor = (1._r8 - EDPftvarcon_inst%maintresp_reduction_intercept(pft)) + &
               EDPftvarcon_inst%maintresp_reduction_intercept(pft) * &
               (1._r8 - EDPftvarcon_inst%maintresp_reduction_curvature(pft)**frac) &
               / (1._r8-EDPftvarcon_inst%maintresp_reduction_curvature(pft))
       else  ! avoid nan answer for linear case
          maintresp_reduction_factor = (1._r8 - EDPftvarcon_inst%maintresp_reduction_intercept(pft)) + &
               EDPftvarcon_inst%maintresp_reduction_intercept(pft) * frac
       endif

    else
       maintresp_reduction_factor = 1._r8
    endif


  end subroutine lowstorage_maintresp_reduction
  
   ! =====================================================================================

   real(r8) function ConvertPar(leaf_area, par_wm2) result(par_umolm2s)
      !
      ! DESCRIPTION:
      ! Convert par from W/m2 to umol photons/m2/s
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: leaf_area ! leaf area index [m2]
      real(r8), intent(in) :: par_wm2   ! absorbed PAR [W/m2]

      if (par_wm2 > nearzero .and. leaf_area > min_la_to_solve) then
         par_umolm2s = par_wm2/leaf_area*wm2_to_umolm2s
      else                 
         ! The radiative transfer schemes are imperfect
         ! they can sometimes generate negative values here if par or leaf area is 0.0
         par_umolm2s = 0.0_r8
      end if

   end function ConvertPar

end module FATESPlantRespPhotosynthMod
