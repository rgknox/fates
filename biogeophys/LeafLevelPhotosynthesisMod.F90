module LeafLevelPhotosynthesisMod

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
  
  implicit none
  private
  
  procedure, public :: LeafLayerPhotosnthesis


  type,public :: psn_params_type

     ! Parameters dimensioned by PFT
     real(r8), allocatable :: bb_slope(:)           ! slope of BB relationship, unitless
     real(r8), allocatable :: medlyn_slope(:)       ! slope for Medlyn stomatal conductance model method, the unit is KPa^0.5
     real(r8), allocatable :: stomatal_intercept(:) ! Unstressed minimum stomatal conductance, the unit is umol/m**2/s
     integer,  allocatable :: psn_path_index(:)     ! Index specifying if the plant of interest is C3/C4

     real(r8), allocatable :: vcmaxha(:)   ! activation energy for vcmax (J/mol)
     real(r8), allocatable :: jmaxha(:)    ! activation energy for jmax (J/mol)
     real(r8), allocatable :: tpuha(:)     ! activation energy for tpu (J/mol)
     real(r8), allocatable :: vcmaxhd(:)   ! deactivation energy for vcmax (J/mol)
     real(r8), allocatable :: jmaxhd(:)    ! deactivation energy for jmax (J/mol)
     real(r8), allocatable :: tpuhd(:)     ! deactivation energy for tpu (J/mol)
     real(r8), allocatable :: vcmaxse(:)   ! entropy term for vcmax (J/mol/K)
     real(r8), allocatable :: jmaxse(:)    ! entropy term for jmax (J/mol/K)
     real(r8), allocatable :: tpuse(:)     ! entropy term for tpu (J/mol/K)

     real(r8), allocatable :: vcmaxc(:)    ! scaling factor for high temperature inhibition (25 C = 1.0)
     real(r8), allocatable :: jmaxc(:)     ! scaling factor for high temperature inhibition (25 C = 1.0)
     real(r8), allocatable :: tpuc(:)      ! scaling factor for high temperature inhibition (25 C = 1.0)
     
     ! Scalar parameters
     integer               :: stoma_model           ! Switch: 1 for Ball-Berry, 2 for Medlyn
     
  end type psn_params_type


  type(psn_params_type), public :: psn_params

  
  integer, parameter :: bb_stoma_model = 1
  integer, parameter :: medlyn_stoma_model = 2

  integer, parameter :: c3_path = 1
  integer, parameter :: c4_path = 2

  
  ! maximum stomatal resistance [s/m] (used across several procedures)
  real(r8),public, parameter :: rsmax0 =  2.e8_r8                
   
contains


  subroutine InitHighTempInhibition()


    integer :: ft

    do ft = 1, size(psm_params%vcmaxhd)
       psn_params%vcmaxc(ft) = fth25_f(psn_params%vcmaxhd(ft), psn_params%vcmaxse(ft))
       psn_params%jmaxc(ft)  = fth25_f(psn_params%jmaxhd(ft), psn_params%jmaxse(ft))
       psn_params%tpuc(ft)   = fth25_f(psn_params%tpuhd(ft), psn_params%tpuse(ft))
    end do
       
    return
  end subroutine InitHighTempInhibition


  

  
  subroutine LeafLayerPhotosynthesis(f_sun_lsl,         &  ! in
       parsun_leaf,       &  ! in
       parsha_leaf,       &  ! in
       ft,                &  ! in
       vcmax,             &  ! in
       jmax,              &  ! in
       tpu,               &  ! in
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


    ! Arguments
    ! ------------------------------------------------------------------------------------
    real(r8), intent(in) :: f_sun_lsl         ! 
    real(r8), intent(in) :: parsun_leaf       ! Absorbed PAR in sunlist leaves [W/m2 leaf]
    real(r8), intent(in) :: parsha_leaf       ! Absorved PAR in shaded leaves  [W/m2 leaf]
    integer,  intent(in) :: ft                ! (plant) Functional Type Index
    real(r8), intent(in) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
    real(r8), intent(in) :: jmax              ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8), intent(in) :: tpu               ! triose phosphate utilization rate (umol CO2/m**2/s)
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
    real(r8), intent(in) :: stoma_slope     ! either BB or Medlyn stomatal slope

    real(r8), intent(out) :: psn_out        ! carbon assimilated in this leaf layer umolC/m2/s
    real(r8), intent(out) :: rstoma_out     ! stomatal resistance (1/gs_lsl) (s/m)
    real(r8), intent(out) :: anet_av_out    ! net leaf photosynthesis (umol CO2/m**2/s) 
                                            ! averaged over sun and shade leaves.
    real(r8), intent(out) :: c13disc_z      ! carbon 13 in newly assimilated carbon 				     

    ! Locals
    integer :: sunsha             ! Index for differentiating sun and shade
    real(r8) :: gstoma            ! Stomatal Conductance of this leaf layer (m/s)
    real(r8) :: agross            ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
    real(r8) :: anet              ! net leaf photosynthesis (umol CO2/m**2/s)
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
    real(r8) :: gs_mol_err        ! gs_mol for error check
    real(r8) :: ac                ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    real(r8) :: aj                ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    real(r8) :: ap                ! product-limited (C3) or CO2-limited  
                                  ! (C4) gross photosynthesis (umol CO2/m**2/s)
    real(r8) :: ai                ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
    real(r8) :: leaf_co2_ppress   ! CO2 partial pressure at leaf surface (Pa)
    real(r8) :: init_co2_inter_c  ! First guess intercellular co2 specific to C path
    real(r8) :: term              ! intermediate variable in Medlyn stomatal conductance model
    real(r8) :: vpd               ! water vapor deficit in Medlyn stomatal model (KPa)


    ! Parameters
    ! ------------------------------------------------------------------------
    ! Fraction of light absorbed by non-photosynthetic pigments
    real(r8),parameter :: fnps = 0.15_r8       

    ! For plants with no leaves, a miniscule amount of conductance
    ! can happen through the stems, at a partial rate of cuticular conductance
    real(r8),parameter :: stem_cuticle_loss_frac = 0.1_r8

    ! empirical curvature parameter for electron transport rate
    real(r8),parameter :: theta_psii = 0.7_r8   

    ! quantum efficiency, used only for C4 (mol CO2 / mol photons) (index 2)
    real(r8),parameter,dimension(2) :: quant_eff = [0.0_r8,0.05_r8]

    ! empirical curvature parameter for ac, aj photosynthesis co-limitation. 
    ! Changed theta_cj and theta_ip to 0.999 to effectively remove smoothing logic 
    ! following Anthony Walker's findings from MAAT. 
    real(r8),parameter,dimension(2) :: theta_cj  = [0.999_r8,0.999_r8]

    ! First guess on ratio between intercellular co2 and the atmosphere
    ! an iterator converges on actual
    real(r8),parameter,dimension(2) :: init_a2l_co2 = [0.7_r8,0.4_r8]


    ! empirical curvature parameter for ap photosynthesis co-limitation
    real(r8),parameter :: theta_ip = 0.999_r8



    init_co2_inter_c = init_a2l_co2(psn_params%psn_path_index(ft)) * can_co2_ppress



    if ( (parsun_leaf+parsha_leaf) < rsnbl_math_prec ) then  ! night time

       ! Trivial case: night-time, no solar radiation

       anet_av_out = -lmr
       psn_out     = 0._r8

       ! The cuticular conductance already factored in maximum resistance as a bound
       ! no need to re-bound it

       rstoma_out = cf/stomatal_intercept_btran
       c13disc_z = 0.0_r8    !carbon 13 discrimination in night time carbon flux, note value of 1.0 is used in CLM

    else

       ! Non-trivial case, there is both some light absorbed
       ! and some leaf tissue to absorb it.

       !Loop aroun shaded and unshaded leaves          
       psn_out     = 0._r8    ! psn is accumulated across sun and shaded leaves. 
       rstoma_out  = 0._r8    ! 1/rs is accumulated across sun and shaded leaves. 
       anet_av_out = 0._r8
       gstoma      = 0._r8

       sunsha_loop: do  sunsha = 1,2      

          ! Electron transport rate for C3 plants. 
          ! Convert par from W/m2 to umol photons/m**2/s using the factor 4.6
          ! Convert from units of par absorbed per unit ground area to par 
          ! absorbed per unit leaf area. 

          if(sunsha == 1)then !sunlit
             qabs = parsun_leaf
          else
             qabs = parsha_leaf
          end if

          qabs = qabs * 0.5_r8 * (1._r8 - fnps) *  4.6_r8 

          !convert the absorbed par into absorbed par per m2 of leaf, 
          ! so it is consistant with the vcmax and lmr numbers. 
          aquad = theta_psii
          bquad = -(qabs + jmax)
          cquad = qabs * jmax
          call quadratic_f (aquad, bquad, cquad, r1, r2)
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
             if (psn_params%psn_path_index(ft) == c3_path)then    

                ! C3: Rubisco-limited photosynthesis
                ac = vcmax * max(co2_inter_c-co2_cpoint, 0._r8) / &
                     (co2_inter_c+mm_kco2 * (1._r8+can_o2_ppress / mm_ko2 ))

                ! C3: RuBP-limited photosynthesis
                aj = je * max(co2_inter_c-co2_cpoint, 0._r8) / &
                     (4._r8*co2_inter_c+8._r8*co2_cpoint)

                ! C3: Product-limited photosynthesis 
                ap = 3._r8 * tpu

             else

                ! C4: Rubisco-limited photosynthesis
                ac = vcmax

                ! C4: RuBP-limited photosynthesis
                if(sunsha == 1)then !sunlit

                   aj = quant_eff(psn_params%psn_path_index(ft)) * parsun_leaf * 4.6_r8

                else

                   aj = quant_eff(psn_params%psn_path_index(ft)) * parsha_leaf * 4.6_r8

                end if

                ! C4: PEP carboxylase-limited (CO2-limited)
                ap = co2_rcurve_islope * max(co2_inter_c, 0._r8) / can_press  

             end if

             ! Gross photosynthesis smoothing calculations. First co-limit ac and aj.
             ! Then co-limit ap
             
             aquad = theta_cj(psn_params%psn_path_index(ft))
             bquad = -(ac + aj)
             cquad = ac * aj
             call quadratic_f (aquad, bquad, cquad, r1, r2)
             ai = min(r1,r2)

             aquad = theta_ip
             bquad = -(ai + ap)
             cquad = ai * ap
             call quadratic_f (aquad, bquad, cquad, r1, r2)
             agross = min(r1,r2)

             ! Net carbon assimilation. Exit iteration if an < 0
             anet = agross  - lmr
             if (anet < 0._r8) then
                loop_continue = .false.
             end if

             ! Quadratic gs_mol calculation with an known. Valid for an >= 0.
             ! With an <= 0, then gs_mol = stomatal_intercept_btran                 
             leaf_co2_ppress = can_co2_ppress - h2o_co2_bl_diffuse_ratio/gb_mol * anet * can_press 
             leaf_co2_ppress = max(leaf_co2_ppress,1.e-06_r8)
             if ( psn_params%stoma_model == medlyn_stoma_model ) then

                ! stomatal conductance calculated from Medlyn et al. (2011), the numerical &
                ! implementation was adapted from the equations in CLM5.0
                ! addapted from CLM5. Put some constraint on VPD
                vpd =  max((veg_esat - ceair), 50._r8) * 0.001_r8          

                ! when Medlyn stomatal conductance is being used,
                ! the unit is KPa. Ignoring the constraint will cause errors when model runs.          
                term = h2o_co2_stoma_diffuse_ratio * anet / (leaf_co2_ppress / can_press)
                aquad = 1.0_r8
                bquad = -(2.0 * (stomatal_intercept_btran + term) + (psn_params%medlyn_slope(ft) * term)**2. / &
                     (gb_mol * vpd ))
                cquad = stomatal_intercept_btran*stomatal_intercept_btran + &
                     (2.0*stomatal_intercept_btran + term * &
                     (1.0 - psn_params%medlyn_slope(ft)**2. / vpd)) * term

                call quadratic_f (aquad, bquad, cquad, r1, r2)
                gs_mol = max(r1,r2)

             else if ( psn_params%stoma_model == bb_stoma_model ) then         !stomatal conductance calculated from Ball et al. (1987)
                aquad = leaf_co2_ppress
                bquad = leaf_co2_ppress*(gb_mol - stomatal_intercept_btran) - psn_params%bb_slope(ft) * anet * can_press
                cquad = -gb_mol*(leaf_co2_ppress*stomatal_intercept_btran + &
                     psn_params%bb_slope(ft)*anet*can_press * ceair/ veg_esat )

                call quadratic_f (aquad, bquad, cquad, r1, r2)
                gs_mol = max(r1,r2)
             end if
             ! Derive new estimate for co2_inter_c
             co2_inter_c = can_co2_ppress - anet * can_press * &
                  (h2o_co2_bl_diffuse_ratio*gs_mol+h2o_co2_stoma_diffuse_ratio*gb_mol) / (gb_mol*gs_mol)

             ! Check for co2_inter_c convergence. Delta co2_inter_c/pair = mol/mol. 
             ! Multiply by 10**6 to convert to umol/mol (ppm). Exit iteration if 
             ! convergence criteria of +/- 1 x 10**-6 ppm is met OR if at least ten 
             ! iterations (niter=10) are completed

             if ((abs(co2_inter_c-co2_inter_c_old)/can_press*1.e06_r8 <=  2.e-06_r8) &
                  .or. niter == 5) then
                loop_continue = .false.
             end if

          end do iter_loop

          ! End of co2_inter_c iteration.  Check for an < 0, in which case
          ! gs_mol =stomatal_intercept_btran 
          if (anet < 0._r8) then
             gs_mol = stomatal_intercept_btran
          end if

          ! Final estimates for leaf_co2_ppress and co2_inter_c 
          ! (needed for early exit of co2_inter_c iteration when an < 0)
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

          ! Compare with Medlyn model: gs_mol = 1.6*(1+m/sqrt(vpd)) * an/leaf_co2_ppress*p + b
          if ( psn_params%stoma_model == medlyn_stoma_model ) then
             gs_mol_err = h2o_co2_stoma_diffuse_ratio * &
                  (1 + psn_params%medlyn_slope(ft)/sqrt(vpd)) * &
                  max(anet,0._r8)/leaf_co2_ppress*can_press + stomatal_intercept_btran
             ! Compare with Ball-Berry model: gs_mol = m * an * hs/leaf_co2_ppress*p + b             
          else if ( psn_params%stoma_model == bb_stoma_model ) then 
             hs = (gb_mol*ceair + gs_mol* veg_esat ) / ((gb_mol+gs_mol)*veg_esat )
             gs_mol_err = psn_params%bb_slope(ft)*max(anet, 0._r8)*hs/leaf_co2_ppress*can_press + stomatal_intercept_btran
          end if

          if (abs(gs_mol-gs_mol_err) > 1.e-01_r8) then
             write (fates_log(),*) 'Stomatal model error check - stomatal conductance error:'
             write (fates_log(),*) gs_mol, gs_mol_err
          end if

       enddo sunsha_loop

       ! This is the stomatal resistance of the leaf layer
       rstoma_out = 1._r8/gstoma


    end if

    return
  end subroutine LeafLayerPhotosynthesis

  ! =====================================================================================
  
  subroutine LeafLayerMaintenanceRespiration(lnc,   &
                                             nscaler,   &
                                             ft,        &
                                             veg_tempk,     &
                                             lmr)

    
 
      ! Arguments
      real(r8), intent(in)  :: lnc          ! leaf nitrogen concentration of this layer
      real(r8), intent(in)  :: lmr25        ! canopy top leaf maint resp rate at 25C 
                                            ! for this pft (umol CO2/m**2/s)
      integer,  intent(in)  :: ft           ! (plant) Functional Type Index
      real(r8), intent(in)  :: nscaler      ! Scale for leaf nitrogen profile
      real(r8), intent(in)  :: veg_tempk    ! vegetation temperature
      real(r8), intent(out) :: lmr          ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
      
      ! Locals
      real(r8) :: lmr25   ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
      
      ! Parameter
      real(r8), parameter :: lmrha = 46390._r8    ! activation energy for lmr (J/mol)
      real(r8), parameter :: lmrhd = 150650._r8   ! deactivation energy for lmr (J/mol)
      real(r8), parameter :: lmrse = 490._r8      ! entropy term for lmr (J/mol/K)
      real(r8), parameter :: lmrc = 1.15912391_r8 ! scaling factor for high 
                                                  ! temperature inhibition (25 C = 1.0)

      ! REF FOR THIS EQUATION!!?
      
      lmr25 = 2.525e-6_r8 * (1.5_r8 ** ((25._r8 - 20._r8)/10._r8))
      lmr25 = lmr25 * lnc / (umolC_to_kgC * g_per_kg)


      ! Leaf Maintenance respiration: umol CO2 / m**2 [leaf] / s
      ! ----------------------------------------------------------------------------------
      
      if ( psn_params%psn_path_index(ft) == c3_path )then
         lmr = lmr25 * ft1_f(veg_tempk, lmrha) * &
               fth_f(veg_tempk, lmrhd, lmrse, lmrc)
      else
         lmr = lmr25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
         lmr = lmr / (1._r8 + exp( 1.3_r8*(veg_tempk-(tfrz+55._r8)) ))
      end if
      
      ! Any hydrodynamic limitations could go here, currently none
      ! lmr = lmr * (nothing)
      
    end subroutine LeafLayerMaintenanceRespiration
    
  ! =====================================================================================
  

   subroutine LeafLayerBiophysicalRates( parsun_lsl, &
                                         ft,            &
                                         vcmax25top_ft, &
                                         jmax25top_ft, &
                                         tpu25top_ft, &
                                         co2_rcurve_islope25top_ft, &
                                         nscaler,    &
                                         veg_tempk,      &
                                         btran, &
                                         vcmax, &
                                         jmax, &
                                         tpu, &
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
      ! tpu: triose phosphate utilization rate and
      ! co2_rcurve_islope: initial slope of CO2 response curve (C4 plants)
      ! ---------------------------------------------------------------------------------


      ! Arguments
      ! ------------------------------------------------------------------------------
      
      real(r8), intent(in) :: parsun_lsl      ! PAR absorbed in sunlit leaves for this layer
      integer,  intent(in) :: ft              ! (plant) Functional Type Index
      real(r8), intent(in) :: nscaler         ! Scale for leaf nitrogen profile
      real(r8), intent(in) :: vcmax25top_ft   ! canopy top maximum rate of carboxylation at 25C 
                                              ! for this pft (umol CO2/m**2/s)
      real(r8), intent(in) :: jmax25top_ft    ! canopy top maximum electron transport rate at 25C 
                                              ! for this pft (umol electrons/m**2/s)
      real(r8), intent(in) :: tpu25top_ft     ! canopy top triose phosphate utilization rate at 25C 
                                              ! for this pft (umol CO2/m**2/s)
      real(r8), intent(in) :: co2_rcurve_islope25top_ft ! initial slope of CO2 response curve
                                              ! (C4 plants) at 25C, canopy top, this pft
      real(r8), intent(in) :: veg_tempk       ! vegetation temperature
      real(r8), intent(in) :: btran           ! transpiration wetness factor (0 to 1) 
                                                    
      real(r8), intent(out) :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
      real(r8), intent(out) :: jmax              ! maximum electron transport rate 
                                                 ! (umol electrons/m**2/s)
      real(r8), intent(out) :: tpu               ! triose phosphate utilization rate 
                                                 ! (umol CO2/m**2/s)
      real(r8), intent(out) :: co2_rcurve_islope ! initial slope of CO2 response curve (C4 plants)
      
      ! Locals
      ! -------------------------------------------------------------------------------
      real(r8) :: vcmax25             ! leaf layer: maximum rate of carboxylation at 25C
                                      ! (umol CO2/m**2/s)
      real(r8) :: jmax25              ! leaf layer: maximum electron transport rate at 25C 
                                      ! (umol electrons/m**2/s)
      real(r8) :: tpu25               ! leaf layer: triose phosphate utilization rate at 25C 
                                      ! (umol CO2/m**2/s)
      real(r8) :: co2_rcurve_islope25 ! leaf layer: Initial slope of CO2 response curve 
                                      ! (C4 plants) at 25C
      
      ! Parameters
      ! ---------------------------------------------------------------------------------
      
      

      if ( parsun_lsl <= 0._r8) then           ! night time
         vcmax             = 0._r8
         jmax              = 0._r8
         tpu               = 0._r8
         co2_rcurve_islope = 0._r8
      else                                     ! day time

         ! Scale by attenuation of nitrogen in canopy
         vcmax25 = vcmax25top_ft * nscaler
         jmax25  = jmax25top_ft * nscaler
         tpu25   = tpu25top_ft * nscaler
         co2_rcurve_islope25 = co2_rcurve_islope25top_ft * nscaler
         
         ! Adjust for temperature
         vcmax = vcmax25 * ft1_f(veg_tempk, psn_params%vcmaxha) * &
              fth_f(veg_tempk, psn_params%vcmaxhd(ft), psn_params%vcmaxse(ft), psn_params%vcmaxc(ft))
         jmax  = jmax25 * ft1_f(veg_tempk, psn_params%jmaxha) * &
              fth_f(veg_tempk, psn_params%jmaxhd(ft), psn_params%jmaxse(ft), psn_params%jmaxc(ft))
         tpu   = tpu25 * ft1_f(veg_tempk, psn_params%tpuha) * &
              fth_f(veg_tempk, psn_params%tpuhd(ft), psn_params%tpuse,(ft) psn_params%tpuc(ft))
         
         if ( psn_params%psn_path_index(ft)  /= c3_path ) then
            vcmax = vcmax25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
            vcmax = vcmax / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-veg_tempk ) ))
            vcmax = vcmax / (1._r8 + exp( 0.3_r8*(veg_tempk-(tfrz+40._r8)) ))
         end if
         !q10 response of product limited psn. 
         co2_rcurve_islope = co2_rcurve_islope25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8) 
      end if
      
      ! Adjust for water limitations 
      vcmax = vcmax * btran
      
      return
    end subroutine LeafLayerBiophysicalRates

    ! ====================================================================================
    
  subroutine quadratic_f (a, b, c, r1, r2)

    !------------------------------------------------------------------------------
    ! !REVISION HISTORY:
    ! 4/5/10: Adapted from /home/bonan/ecm/psn/An_gs_iterative.f90 by Keith Oleson
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: a,b,c       ! Terms for quadratic equation
    real(r8), intent(out) :: r1,r2       ! Roots of quadratic equation
    !
    ! !LOCAL VARIABLES:
    real(r8) :: q                        ! Temporary term for quadratic solution
    !------------------------------------------------------------------------------

    if (a == 0._r8) then
       write (fates_log(),*) 'Quadratic solution error: a = ',a
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if (b >= 0._r8) then
       q = -0.5_r8 * (b + sqrt(b*b - 4._r8*a*c))
    else
       q = -0.5_r8 * (b - sqrt(b*b - 4._r8*a*c))
    end if

    r1 = q / a
    if (q /= 0._r8) then
       r2 = c / q
    else
       r2 = 1.e36_r8
    end if

  end subroutine quadratic_f

  ! ====================================================================================

  subroutine quadratic_fast (a, b, c, r1, r2)

    !-----------------------------------------------------------------------------
    ! Solution from Press et al (1986) Numerical Recipes: The Art of Scientific
    ! Computing (Cambridge University Press, Cambridge), pp. 145.
    !
    ! !REVISION HISTORY:
    ! 4/5/10: Adapted from /home/bonan/ecm/psn/An_gs_iterative.f90 by Keith Oleson
    ! 7/23/16: Copied over from CLM by Ryan Knox
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: a,b,c       ! Terms for quadratic equation
    real(r8), intent(out) :: r1,r2       ! Roots of quadratic equation
    !
    ! !LOCAL VARIABLES:
    real(r8) :: q                        ! Temporary term for quadratic solution
    real(r8) :: term
    real(r8) :: onesign
    !------------------------------------------------------------------------------

    !  if (a == 0._r8) then
    !     write (fates_log(),*) 'Quadratic solution error: a = ',a
    !     call endrun(msg=errMsg(sourcefile, __LINE__))
    !  end if

    onesign = 1.0_r8
    q = -0.5_r8 * (b + sign(onesign,b)*sqrt(b*b - 4._r8*a*c))

    
    !if (b >= 0._r8) then
       
    !   q = -0.5_r8 * (b + sqrt(b*b - 4._r8*a*c))
    !else
    !   q = -0.5_r8 * (b - sqrt(b*b - 4._r8*a*c))
    !end if

    r1 = q / a
    !  if (q /= 0._r8) then
    r2 = c / q
    !  else
    !     r2 = 1.e36_r8
    !  end if

  end subroutine quadratic_fast

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
     use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    
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
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
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
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm

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

end module LeafLevelPhotosynthesisMod
