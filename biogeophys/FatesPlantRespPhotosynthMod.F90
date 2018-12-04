module FATESPlantRespPhotosynthMod
   
   !-------------------------------------------------------------------------------------
   ! !DESCRIPTION:
   ! Calculates the plant respiration and photosynthetic fluxes for the FATES model
   ! This code is similar to and was originally based off of the 'photosynthesis' 
   ! subroutine in the CLM model.
   !
   ! Parameter for activation and deactivation energies were taken from:
   ! Activation energy, from:
   ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
   ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
   ! except TPU from: Harley et al (1992) Plant, Cell and Environment 15:271-282
   ! High temperature deactivation, from:
   ! Leuning (2002) Plant, Cell and Environment 25:1205-1210
   ! The factor "c" scales the deactivation to a value of 1.0 at 25C
   ! Photosynthesis and stomatal conductance parameters, from:
   ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
   ! ------------------------------------------------------------------------------------
   
   ! !USES:

   use FatesGlobals,      only : endrun => fates_endrun
   use FatesGlobals,      only : fates_log
   use FatesConstantsMod, only : r8 => fates_r8
   use FatesConstantsMod, only : itrue
   use FatesConstantsMod, only : mol_per_umol
   use FatesConstantsMod, only : umol_per_mol
   use FatesConstantsMod, only : kpa_per_pa
   use FatesInterfaceMod, only : hlm_use_planthydro
   use FatesInterfaceMod, only : hlm_parteh_mode
   use FatesInterfaceMod, only : numpft
   use FatesConstantsMod, only : nearzero
   use EDPftvarcon,       only : EDPftvarcon_inst 
   use EDTypesMod,        only : maxpft
   use EDTypesMod,        only : nlevleaf
   use EDTypesMod,        only : nclmax

   use PRTGenericMod,     only : prt_carbon_allom_hyp
   use PRTGenericMod,     only : prt_cnp_flex_allom_hyp 
   use PRTGenericMod,     only : all_carbon_elements
   use PRTGenericMod,     only : nitrogen_element
   use PRTGenericMod,     only : leaf_organ
   use PRTGenericMod,     only : fnrt_organ
   use PRTGenericMod,     only : sapw_organ
   use PRTGenericMod,     only : store_organ
   use PRTGenericMod,     only : repro_organ
   use PRTGenericMod,     only : struct_organ

   ! CIME Globals
   use shr_log_mod , only      : errMsg => shr_log_errMsg

   implicit none
   private
   
   public :: FatesPlantRespPhotosynthDrive ! Called by the HLM-Fates interface
   
   character(len=*), parameter, private :: sourcefile = &
         __FILE__
   !-------------------------------------------------------------------------------------
   
   ! maximum stomatal resistance [s/m] (used across several procedures)
   real(r8),parameter :: rsmax0 =  2.e4_r8                    
   
   logical   ::  debug = .false.

   ! quadratic result definition
   integer, parameter :: upper_root = 1
   integer, parameter :: lower_root = 2


   ! ------------------------------------------------------------------------
   ! Fraction of light absorbed by non-photosynthetic pigments
   real(r8),parameter :: fnps = 0.15_r8       
   
   ! empirical curvature parameter for electron transport rate
   real(r8),parameter :: theta_psii = 0.7_r8   
   
   ! First guess on ratio between intercellular co2 and the atmosphere
   ! an iterator converges on actual
   real(r8),parameter :: init_a2l_co2_c3 = 0.7_r8
   real(r8),parameter :: init_a2l_co2_c4 = 0.4_r8

   ! quantum efficiency, used only for C4 (mol CO2 / mol photons) (index 0)
   real(r8),parameter,dimension(0:1) :: quant_eff = [0.05_r8,0.0_r8]

   ! empirical curvature parameter for ac, aj photosynthesis co-limitation
   real(r8),parameter,dimension(0:1) :: theta_cj  = [0.80_r8,0.98_r8]

   ! empirical curvature parameter for ap photosynthesis co-limitation
   real(r8),parameter :: theta_ip = 0.95_r8

   ! Ratio of H2O/CO2 gas diffusion in stomatal airspace (approximate)
   real(r8),parameter :: h2o_co2_stoma_diffuse_ratio = 1.6_r8
   
   ! Ratio of H2O/CO2 gass diffusion in the leaf boundary layer (approximate) 
   real(r8),parameter :: h2o_co2_bl_diffuse_ratio = 1.4_r8

   ! THIS SHOULD BE A UNIFIED un_initialized in fates constants
   real(r8),parameter :: un_initialized = -9.9e32_r8


contains
  
  !--------------------------------------------------------------------------------------
  
  subroutine FatesPlantRespPhotosynthDrive (nsites, sites,bc_in,bc_out,dtime)

    ! -----------------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Leaf photosynthesis and stomatal conductance calculation as described by
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
    ! a multi-layer canopy
    ! -----------------------------------------------------------------------------------


    ! !USES:

    

    use FatesSynchronizedParamsMod , only : FatesSynchronizedParamsInst
    use EDTypesMod        , only : ed_patch_type
    use EDTypesMod        , only : ed_cohort_type
    use EDTypesMod        , only : ed_site_type
    use EDTypesMod        , only : maxpft
    use EDTypesMod        , only : dinc_ed
    use FatesInterfaceMod , only : bc_in_type
    use FatesInterfaceMod , only : bc_out_type
    use EDCanopyStructureMod, only : calc_areaindex
    use FatesConstantsMod, only : umolC_to_kgC
    use FatesConstantsMod, only : g_per_kg
    use FatesConstantsMod, only : umol_per_mmol
    use FatesConstantsMod, only : rgas => rgas_J_K_kmol
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesParameterDerivedMod, only : param_derived
    use EDParamsMod, only : ED_val_bbopt_c3, ED_val_bbopt_c4, ED_val_base_mr_20
    use FatesAllometryMod, only : bleaf
    use FatesAllometryMod, only : storage_fraction_of_target
    use FatesAllometryMod, only : set_root_fraction
    use FatesAllometryMod, only : i_hydro_rootprof_context
    use FatesAllometryMod, only : decay_coeff_kn

    ! ARGUMENTS:
    ! -----------------------------------------------------------------------------------
    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)
    real(r8),intent(in)                     :: dtime


    ! LOCAL VARIABLES:
    ! -----------------------------------------------------------------------------------
    type (ed_patch_type) , pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort

    ! -----------------------------------------------------------------------------------
    ! These three arrays hold leaf-level biophysical rates that are calculated
    ! in one loop and then sent to the cohorts in another loop.  If hydraulics are
    ! on, we calculate a unique solution for each level-cohort-layer combination.
    ! If we are not using hydraulics, we calculate a unique solution for each
    ! level-pft-layer combination.  Thus the following three arrays are statically 
    ! allocated for the maximum space of the two cases (numCohortsPerPatch)
    ! The "_z" suffix indicates these variables are discretized at the "leaf_layer"
    ! scale.
    ! Note: For these temporary arrays, we have the leaf layer dimension first
    ! and the canopy layer last. This order is chosen for efficiency. The arrays
    ! such as leaf area that are bound to the patch structure DO NOT follow this order
    ! as they are used in many other parts of the code with different looping, we
    ! are not modifying its order now.
    ! -----------------------------------------------------------------------------------

    ! leaf maintenance (dark) respiration [umol CO2/m**2/s]
    real(r8) :: lmr_z(nlevleaf,maxpft,nclmax)

    ! stomatal resistance [s/m]
    real(r8) :: rs_z(nlevleaf,maxpft,nclmax)    

    ! net leaf photosynthesis averaged over sun and shade leaves. [umol CO2/m**2/s]
    real(r8) :: anet_av_z(nlevleaf,maxpft,nclmax)  
    
    ! Mask used to determine which leaf-layer biophysical rates have been
    ! used already
    logical :: rate_mask_z(nlevleaf,maxpft,nclmax)

    real(r8) :: vcmax_z            ! leaf layer maximum rate of carboxylation 
                                   ! (umol co2/m**2/s)
    real(r8) :: jmax_z             ! leaf layer maximum electron transport rate 
                                   ! (umol electrons/m**2/s)
    real(r8) :: tpu_z              ! leaf layer triose phosphate utilization rate 
                                   ! (umol CO2/m**2/s)
    real(r8) :: kp_z               ! leaf layer initial slope of CO2 response 
                                   ! curve (C4 plants)
   
    real(r8) :: mm_kco2            ! Michaelis-Menten constant for CO2 (Pa)
    real(r8) :: mm_ko2             ! Michaelis-Menten constant for O2 (Pa)
    real(r8) :: mm_km              ! Michaelis-menten combo term applied in Farquar 1980
    real(r8) :: co2_cpoint         ! CO2 compensation point (Pa)
    real(r8) :: btran_eff          ! effective transpiration wetness factor (0 to 1) 
    real(r8) :: bbb                ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(r8) :: kn(maxpft)         ! leaf nitrogen decay coefficient
    real(r8) :: cf                 ! s m**2/umol -> s/m (ideal gas conversion) [umol/m3]
    real(r8) :: gb_mol             ! leaf boundary layer conductance (molar form: [umol /m**2/s])
    real(r8) :: ceair              ! vapor pressure of air, constrained (Pa)
    real(r8) :: nscaler            ! leaf nitrogen scaling coefficient
    real(r8) :: leaf_frac          ! ratio of to leaf biomass to total alive biomass
    real(r8) :: tcsoi              ! Temperature response function for root respiration. 
    real(r8) :: tcwood             ! Temperature response function for wood
    real(r8) :: vpd_kpa            ! Vapor Pressure Deficit (kPa)
    real(r8) :: elai               ! exposed LAI (patch scale)
    real(r8) :: live_stem_n        ! Live stem (above-ground sapwood) 
                                   ! nitrogen content (kgN/plant)
    real(r8) :: live_croot_n       ! Live coarse root (below-ground sapwood) 
                                   ! nitrogen content (kgN/plant)
    real(r8) :: sapw_c             ! Sapwood carbon (kgC/plant)
    real(r8) :: fnrt_c             ! Fine root carbon (kgC/plant)
    real(r8) :: fnrt_n             ! Fine root nitrogen content (kgN/plant)
    real(r8) :: leaf_c             ! Leaf carbon (kgC/plant)
    real(r8) :: leaf_n             ! leaf nitrogen content (kgN/plant)
    real(r8) :: g_sb_leaves        ! Mean combined (stomata+boundary layer) leaf conductance [m/s]
                                   ! over all of the patch's leaves.  The "sb" refers to the combined
                                   ! "s"tomatal and "b"oundary layer.
                                   ! This quantity is relevant on leaf surfaces. It does not
                                   ! have units of /m2 leaf per say, but is implicitly on leaf surfaces
    real(r8) :: r_sb_leaves        ! Mean leaf resistance over all the patch's leaves [s/m]
                                   ! This is the direct reciprocal of g_sb_leaves
    real(r8) :: r_stomata          ! Mean stomatal resistance across all leaves in the patch [s/m]


    real(r8) :: maintresp_reduction_factor  ! factor by which to reduce maintenance 
                                            ! respiration when storage pools are low
    real(r8) :: b_leaf             ! leaf biomass kgC
    real(r8) :: frac               ! storage pool as a fraction of target leaf biomass
    real(r8) :: check_elai         ! This is a check on the effective LAI that is calculated
                                   ! over each cohort x layer.
    real(r8) :: cohort_eleaf_area  ! This is the effective leaf area [m2] reported by each cohort
    real(r8) :: lnc_top            ! Leaf nitrogen content per unit area at canopy top [gN/m2]
    real(r8) :: lmr25top           ! canopy top leaf maint resp rate at 25C 
                                   ! for this plant or pft (umol CO2/m**2/s)
    real(r8) :: leaf_inc           ! LAI-only portion of the vegetation increment of dinc_ed
    real(r8) :: lai_canopy_above   ! the LAI in the canopy layers above the layer of interest
    real(r8) :: lai_layers_above   ! the LAI in the leaf layers, within the current canopy, 
                                   ! above the leaf layer of interest
    real(r8) :: lai_current        ! the LAI in the current leaf layer
    real(r8) :: cumulative_lai     ! the cumulative LAI, top down, to the leaf layer of interest


    
    ! -----------------------------------------------------------------------------------
    ! Keeping these two definitions in case they need to be added later
    !
    ! -----------------------------------------------------------------------------------
    !real(r8) :: psncanopy_pa  ! patch sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    !real(r8) :: lmrcanopy_pa  ! patch sunlit leaf maintenance respiration rate (umol CO2/m**2/s) 

    integer  :: cl,s,iv,j,ps,ft,ifp ! indices
    integer  :: nv                  ! number of leaf layers
    integer  :: NCL_p               ! number of canopy layers in patch

    ! Choose between the to photosynthesis solvers 
    integer, parameter :: iterative_quad    = 1
    integer, parameter :: semianalytic_quad = 2

    integer, parameter :: photo_solver = semianalytic_quad


    ! Parameters
    ! -----------------------------------------------------------------------
    ! Base maintenance respiration rate for plant tissues base_mr_20
    ! M. Ryan, 1991. Effects of climate change on plant respiration.
    ! Ecological Applications, 1(2), 157-167.
    ! Original expression is br = 0.0106 molC/(molN h)
    ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
    !
    ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
    ! (gC/gN/s)
    ! ------------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------
    ! Photosynthesis and stomatal conductance parameters, from:
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
    ! -----------------------------------------------------------------------------------

    ! Ball-Berry minimum leaf conductance, unstressed (umol H2O/m**2/s)
    ! For C3 and C4 plants
    ! -----------------------------------------------------------------------------------
    real(r8), dimension(0:1) :: bbbopt 

    associate(  &
         c3psn     => EDPftvarcon_inst%c3psn  , &
         slatop    => EDPftvarcon_inst%slatop , & ! specific leaf area at top of canopy, 
                                                  ! projected area basis [m^2/gC]
         woody     => EDPftvarcon_inst%woody  , & ! Is vegetation woody or not? 
         q10       => FatesSynchronizedParamsInst%Q10 )

      bbbopt(0) = ED_val_bbopt_c4
      bbbopt(1) = ED_val_bbopt_c3
      
      do s = 1,nsites

         ! Multi-layer parameters scaled by leaf nitrogen profile.
         ! Loop through each canopy layer to calculate nitrogen profile using
         ! cumulative lai at the midpoint of the layer
         
         ifp = 0
         currentpatch => sites(s)%oldest_patch
         do while (associated(currentpatch))  

            ifp   = ifp+1
            NCL_p = currentPatch%NCL_p
            
            ! Part I. Zero output boundary conditions
            ! ---------------------------------------------------------------------------
            bc_out(s)%rssun_pa(ifp)     = 0._r8
            bc_out(s)%rssha_pa(ifp)     = 0._r8

            g_sb_leaves = 0._r8
            check_elai  = 0._r8

            ! Part II. Filter out patches 
            ! Patch level filter flag for photosynthesis calculations
            ! has a short memory, flags:
            ! 1 = patch has not been called
            ! 2 = patch is currently marked for photosynthesis
            ! 3 = patch has been called for photosynthesis already
            ! ---------------------------------------------------------------------------
            if(bc_in(s)%filter_photo_pa(ifp)==2)then


               ! Part III. Calculate the number of sublayers for each pft and layer.  
               ! And then identify which layer/pft combinations have things in them.  
               ! Output:
               ! currentPatch%ncan(:,:)
               ! currentPatch%canopy_mask(:,:)
               call UpdateCanopyNCanNRadPresent(currentPatch)


               ! Part IV.  Identify some environmentally derived parameters:
               !           These quantities are biologically irrelevant
               !  Michaelis-Menten constant for CO2 (Pa)
               !  Michaelis-Menten constant for O2 (Pa)
               !  CO2 compensation point (Pa)
               !  leaf boundary layer conductance of h20
               !  constrained vapor pressure
               call GetCanopyGasParameters(bc_in(s)%forc_pbot,       & ! in
                                           bc_in(s)%oair_pa(ifp),    & ! in
                                           bc_in(s)%t_veg_pa(ifp),   & ! in
                                           bc_in(s)%tgcm_pa(ifp),    & ! in
                                           bc_in(s)%eair_pa(ifp),    & ! in
                                           bc_in(s)%esat_tv_pa(ifp), & ! in
                                           bc_in(s)%rb_pa(ifp),      & ! in
                                           mm_kco2,                  & ! out              
                                           mm_ko2,                   & ! out
                                           mm_km,                    & ! out
                                           co2_cpoint,               & ! out
                                           cf,                       & ! out
                                           gb_mol,                   & ! out
                                           ceair)                      ! out

               ! Part V.  Pre-process some variables that are PFT dependent
               ! but not environmentally dependent
               ! ------------------------------------------------------------------------

               do ft = 1,numpft

                  ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
                  ! kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al 
                  ! (2010) Biogeosciences, 7, 1833-1859
                  
                  kn(ft) = decay_coeff_kn(ft)
                  

                  ! This is probably unnecessary and already calculated
                  ! ALSO, THIS ROOTING PROFILE IS USED TO CALCULATE RESPIRATION
                  ! YET IT USES THE PROFILE THAT IS CONSISTENT WITH WATER UPTAKE
                  ! AND NOT THE PROFILE WE USE FOR DECOMPOSITION
                  ! SEEMS LIKE THE LATTER WOULD BE MORE APPROPRIATE, RIGHT? (RGK 05-2018)
                  call set_root_fraction(currentPatch%rootfr_ft(ft,1:bc_in(s)%nlevsoil), ft, &
                       bc_in(s)%zi_sisl,lowerb=lbound(bc_in(s)%zi_sisl,1), &
                       icontext = i_hydro_rootprof_context)
                  
               end do !ft 

               

               ! ------------------------------------------------------------------------
               ! Part VI: Loop over all leaf layers.
               ! The concept of leaf layers is a result of the radiative transfer scheme.
               ! A leaf layer has uniform radiation environment.  Leaf layers are a group
               ! of vegetation surfaces (stems and leaves) which inhabit the same 
               ! canopy-layer "CL", have the same functional type "ft" and within those
               ! two partitions are further partitioned into vertical layers where
               ! downwelling radiation attenuates in order.
               ! In this phase we loop over the leaf layers and calculate the
               ! photosynthesis and respiration of the layer (since all biophysical
               ! properties are homogeneous).  After this step, we can loop through
               ! our cohort list, associate each cohort with its list of leaf-layers
               ! and transfer these quantities to the cohort.
               ! With plant hydraulics, we must realize that photosynthesis and
               ! respiration will be different for leaves of each cohort in the leaf
               ! layers, as they will have there own hydraulic limitations.
               ! NOTE: Only need to flush mask on the number of used pfts, not the whole 
               ! scratch space.
               ! ------------------------------------------------------------------------
               rate_mask_z(:,1:numpft,:) = .false.

               if(currentPatch%countcohorts > 0.0)then   ! Ignore empty patches

                  currentCohort => currentPatch%tallest
                  do while (associated(currentCohort)) ! Cohort loop
                     
                     ! Identify the canopy layer (cl), functional type (ft)
                     ! and the leaf layer (IV) for this cohort
                     ft = currentCohort%pft
                     cl = currentCohort%canopy_layer
                     
                     call bleaf(currentCohort%dbh,currentCohort%pft,currentCohort%canopy_trim,b_leaf)
                     call storage_fraction_of_target(b_leaf, &
                           currentCohort%prt%GetState(store_organ, all_carbon_elements), &
                           frac)
                     call lowstorage_maintresp_reduction(frac,currentCohort%pft, &
                          maintresp_reduction_factor)

                     ! are there any leaves of this pft in this layer?
                     if(currentPatch%canopy_mask(cl,ft) == 1)then 
                        
                        ! Loop over leaf-layers
                        do iv = 1,currentCohort%nv
                           
                           ! ------------------------------------------------------------
                           ! If we are doing plant hydro-dynamics (or any run-type
                           ! where cohorts may generate different photosynthetic rates
                           ! of other cohorts in the same canopy-pft-layer combo), 
                           ! we re-calculate the leaf biophysical rates for the 
                           ! cohort-layer combo of interest.  
                           ! but in the vanilla case, we only re-calculate if it has
                           ! not been done yet.
                           ! ------------------------------------------------------------
                           
                           if ( .not.rate_mask_z(iv,ft,cl) .or. &
                                 (hlm_use_planthydro.eq.itrue) .or. &
                                 (hlm_parteh_mode .ne. prt_carbon_allom_hyp )   ) then
                              
                              if (hlm_use_planthydro.eq.itrue) then


                                 bbb       = max (bbbopt(nint(c3psn(ft)))*currentCohort%co_hydr%btran(1), 1._r8)
                                 btran_eff = currentCohort%co_hydr%btran(1) 
                                 
                                 ! dinc_ed is the total vegetation area index of each "leaf" layer
                                 ! we convert to the leaf only portion of the increment
                                 ! ------------------------------------------------------
                                 leaf_inc    = dinc_ed * &
                                               currentCohort%treelai/(currentCohort%treelai+currentCohort%treesai)
                                 
                                 ! Now calculate the cumulative top-down lai of the current layer's midpoint
                                 lai_canopy_above  = sum(currentPatch%canopy_layer_tlai(1:cl-1)) 
                                 lai_layers_above  = leaf_inc * (iv-1)
                                 lai_current       = min(leaf_inc, currentCohort%treelai - lai_layers_above)
                                 cumulative_lai    = lai_canopy_above + lai_layers_above + 0.5*lai_current 

                              else
                                 
                                 bbb   = max (bbbopt(nint(c3psn(ft)))*currentPatch%btran_ft(ft), 1._r8)
                                 btran_eff = currentPatch%btran_ft(ft)
                                 ! For consistency sake, we use total LAI here, and not exposed
                                 ! if the plant is under-snow, it will be effectively dormant for 
                                 ! the purposes of nscaler

                                 cumulative_lai = sum(currentPatch%canopy_layer_tlai(1:cl-1))  + &
                                                  sum(currentPatch%tlai_profile(cl,ft,1:iv-1)) + &
                                                  0.5*currentPatch%tlai_profile(cl,ft,iv)
                                           

                              end if
                           
                              
                              ! Scale for leaf nitrogen profile
                              nscaler = exp(-kn(ft) * cumulative_lai)
                              
                              ! Leaf maintenance respiration to match the base rate used in CN
                              ! but with the new temperature functions for C3 and C4 plants.

                              ! CN respiration has units:  g C / g N [leaf] / s. This needs to be
                              ! converted from g C / g N [leaf] / s to umol CO2 / m**2 [leaf] / s
                              
                              ! Then scale this value at the top of the canopy for canopy depth
                              ! Leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
                              select case(hlm_parteh_mode)
                              case (prt_carbon_allom_hyp)

                                 lnc_top  = EDPftvarcon_inst%prt_nitr_stoich_p1(ft,leaf_organ)/slatop(ft)
                                 
                              case (prt_cnp_flex_allom_hyp)

                                 leaf_c  = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
                                 leaf_n  = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
                                 lnc_top = leaf_n / (slatop(ft) * leaf_c )

                              end select

                              lmr25top = 2.525e-6_r8 * (1.5_r8 ** ((25._r8 - 20._r8)/10._r8))
                              lmr25top = lmr25top * lnc_top / (umolC_to_kgC * g_per_kg)
                              

                              ! Part VII: Calculate dark respiration (leaf maintenance) for this layer
                              call LeafLayerMaintenanceRespiration( lmr25top,                 &  ! in
                                                                    nscaler,                  &  ! in
                                                                    ft,                       &  ! in
                                                                    bc_in(s)%t_veg_pa(ifp),   &  ! in
                                                                    lmr_z(iv,ft,cl))             ! out
                              
                              ! Part VII: Calculate (1) maximum rate of carboxylation (vcmax), 
                              ! (2) maximum electron transport rate, (3) triose phosphate 
                              ! utilization rate and (4) the initial slope of CO2 response curve 
                              ! (C4 plants). Earlier we calculated their base rates as dictated
                              ! by their plant functional type and some simple scaling rules for
                              ! nitrogen limitation baesd on canopy position (not prognostic).
                              ! These rates are the specific rates used in the actual photosynthesis
                              ! calculations that take localized environmental effects (temperature)
                              ! into consideration.

                              call LeafLayerBiophysicalRates(currentPatch%ed_parsun_z(cl,ft,iv), &  ! in
                                                             ft,                                 &  ! in
                                                             EDPftvarcon_inst%vcmax25top(ft),    &  ! in
                                                             param_derived%jmax25top(ft),        &  ! in
                                                             param_derived%tpu25top(ft),         &  ! in
                                                             param_derived%kp25top(ft),          &  ! in
                                                             nscaler,                            &  ! in
                                                             bc_in(s)%t_veg_pa(ifp),             &  ! in
                                                             btran_eff,                          &  ! in
                                                             vcmax_z,                            &  ! out
                                                             jmax_z,                             &  ! out
                                                             tpu_z,                              &  ! out
                                                             kp_z )                                 ! out
                              
                              ! Part IX: This call calculates the actual photosynthesis for the 
                              ! leaf layer, as well as the stomatal resistance and the net assimilated carbon.

!                              select case(photo_solver)
                              
!                              case(iterative_quad)
                              call LeafLayerPhotosynthesis(currentPatch%f_sun(cl,ft,iv),    &  ! in
                                                        currentPatch%ed_parsun_z(cl,ft,iv), &  ! in
                                                        currentPatch%ed_parsha_z(cl,ft,iv), &  ! in
                                                        currentPatch%ed_laisun_z(cl,ft,iv), &  ! in
                                                        currentPatch%ed_laisha_z(cl,ft,iv), &  ! in
                                                currentPatch%canopy_area_profile(cl,ft,iv), &  ! in
                                                        ft,                                 &  ! in
                                                        vcmax_z,                            &  ! in
                                                        jmax_z,                             &  ! in
                                                        tpu_z,                              &  ! in
                                                        kp_z,                               &  ! in
                                                        bc_in(s)%t_veg_pa(ifp),             &  ! in
                                                        bc_in(s)%esat_tv_pa(ifp),           &  ! in
                                                        bc_in(s)%forc_pbot,                 &  ! in
                                                        bc_in(s)%cair_pa(ifp),              &  ! in
                                                        bc_in(s)%oair_pa(ifp),              &  ! in
                                                        btran_eff,                          &  ! in
                                                        bbb,                                &  ! in
                                                        cf,                                 &  ! in
                                                        gb_mol,                             &  ! in
                                                        ceair,                              &  ! in
                                                        mm_kco2,                            &  ! in
                                                        mm_ko2,                             &  ! in
                                                        co2_cpoint,                         &  ! in
                                                        lmr_z(iv,ft,cl),                    &  ! in
                                                        currentPatch%psn_z(cl,ft,iv),       &  ! out
                                                        rs_z(iv,ft,cl),                     &  ! out
                                                        anet_av_z(iv,ft,cl))                   ! out

!                              case(semianalytic_quad)

                              print*,"esat: ", bc_in(s)%esat_tv_pa(ifp)
                              print*,"ceair: ",ceair

                              ! The medlyn model likes vpd in kpa
                              vpd_kpa = max(nearzero,(bc_in(s)%esat_tv_pa(ifp) - ceair)) * kpa_per_pa

                              call SemiAnalyticalQuad(bc_in(s)%forc_pbot,                 &
                                                      bc_in(s)%cair_pa(ifp),              &  ! in
                                                      mm_km,                              &  ! in
                                                      vpd_kpa,                            &  ! vpd
                                                      gb_mol,                             &  ! in
                                                      cf,                                 &  ! in
                                                      vcmax_z,                            &  ! in
                                                      jmax_z,                             &  ! in
                                                      tpu_z,                              &  ! in
                                                      co2_cpoint,                         &  ! in
                                                      lmr_z(iv,ft,cl),                    &  ! in
                                                      currentPatch%ed_parsun_z(cl,ft,iv), &  ! in
                                                      currentPatch%ed_parsha_z(cl,ft,iv), &  ! in
                                                      currentPatch%ed_laisun_z(cl,ft,iv), &  ! in
                                                      currentPatch%ed_laisha_z(cl,ft,iv), &  ! in
                                              currentPatch%canopy_area_profile(cl,ft,iv), &  ! in
                                                      currentPatch%f_sun(cl,ft,iv),    &  ! in
                                                      ft,                                 &  ! in
                                                      anet_av_z(iv,ft,cl),                &  ! out
                                                      currentPatch%psn_z(cl,ft,iv),       &  ! out
                                                      rs_z(iv,ft,cl))

!                              case default
!                                 write(fates_log(),*) 'A photo solver was chosen that DNE'
!                                 call endrun(msg=errMsg(sourcefile, __LINE__))
!                              end select


                              rate_mask_z(iv,ft,cl) = .true.
                           end if
                        end do

                        ! Zero cohort flux accumulators.
                        currentCohort%npp_tstep  = 0.0_r8
                        currentCohort%resp_tstep = 0.0_r8
                        currentCohort%gpp_tstep  = 0.0_r8
                        currentCohort%rdark      = 0.0_r8
                        currentCohort%resp_m     = 0.0_r8
                        currentCohort%ts_net_uptake = 0.0_r8

                        ! ---------------------------------------------------------------
                        ! Part VII: Transfer leaf flux rates (like maintenance respiration,
                        ! carbon assimilation and conductance) that are defined by the 
                        ! leaf layer (which is area independent, ie /m2) onto each cohort
                        ! (where the rates become per cohort, ie /individual). Most likely
                        ! a sum over layers.
                        ! ---------------------------------------------------------------
                        nv = currentCohort%nv
                        call ScaleLeafLayerFluxToCohort(nv,                                    & !in
                                                        currentPatch%psn_z(cl,ft,1:nv),        & !in
                                                        lmr_z(1:nv,ft,cl),                     & !in
                                                        rs_z(1:nv,ft,cl),                      & !in
                                                        currentPatch%elai_profile(cl,ft,1:nv), & !in
                                                        currentCohort%c_area,                  & !in
                                                        currentCohort%n,                       & !in
                                                        bc_in(s)%rb_pa(ifp),                   & !in
                                                        maintresp_reduction_factor,            & !in
                                                        currentCohort%g_sb_laweight,           & !out
                                                        currentCohort%gpp_tstep,               & !out
                                                        currentCohort%rdark,                   & !out
                                                        cohort_eleaf_area)                       !out
                        
                        ! Net Uptake does not need to be scaled, just transfer directly
                        currentCohort%ts_net_uptake(1:nv) = anet_av_z(1:nv,ft,cl) * umolC_to_kgC

                     else
                        
                        ! In this case, the cohort had no leaves, 
                        ! so no productivity,conductance, transpiration uptake
                        ! or dark respiration
                        cohort_eleaf_area       = 0.0_r8
                        currentCohort%gpp_tstep = 0.0_r8 
                        currentCohort%rdark = 0.0_r8 
                        currentCohort%g_sb_laweight = 0.0_r8 
                        currentCohort%ts_net_uptake(:) = 0.0_r8
                        
                     end if  ! if(currentPatch%canopy_mask(cl,ft) == 1)then
                     

                     ! ------------------------------------------------------------------
                     ! Part VIII: Calculate maintenance respiration in the sapwood and
                     ! fine root pools.
                     ! ------------------------------------------------------------------
                
                     ! Calculate the amount of nitrogen in the above and below ground 
                     ! stem and root pools, used for maint resp
                     ! We are using the fine-root C:N ratio as an approximation for
                     ! the sapwood pools.
                     ! Units are in (kgN/plant)
                     ! ------------------------------------------------------------------

                     sapw_c   = currentCohort%prt%GetState(sapw_organ, all_carbon_elements)
                     fnrt_c   = currentCohort%prt%GetState(fnrt_organ, all_carbon_elements)

                     select case(hlm_parteh_mode)
                     case (prt_carbon_allom_hyp)

                        live_stem_n = EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) * &
                              sapw_c * EDPftvarcon_inst%prt_nitr_stoich_p1(ft,sapw_organ)
                        
                        live_croot_n = (1.0_r8-EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)) * &
                              sapw_c * EDPftvarcon_inst%prt_nitr_stoich_p1(ft,sapw_organ)

                        fnrt_n = fnrt_c * EDPftvarcon_inst%prt_nitr_stoich_p1(ft,fnrt_organ)

                     case(prt_cnp_flex_allom_hyp) 
                     
                        live_stem_n = EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) * &
                             currentCohort%prt%GetState(sapw_organ, nitrogen_element)

                        live_croot_n = (1.0_r8-EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)) * &
                             currentCohort%prt%GetState(sapw_organ, nitrogen_element)

                        fnrt_n = currentCohort%prt%GetState(fnrt_organ, nitrogen_element)

                     case default
                        

                     end select

                     !------------------------------------------------------------------------------
                     ! Calculate Whole Plant Respiration
                     ! (this doesn't really need to be in this iteration at all, surely?)
                     ! Response: (RGK 12-2016): I think the positioning of these calls is 
                     ! appropriate as of now.  Maintenance calculations in sapwood and roots
                     ! vary by cohort and with changing temperature at the minimum, and there are
                     ! no sub-pools chopping up those pools any finer that need to be dealt with.
                     !------------------------------------------------------------------------------
                     
                     ! Live stem MR (kgC/plant/s) (above ground sapwood)
                     ! ------------------------------------------------------------------
                     if (woody(ft) == 1) then
                        tcwood = q10**((bc_in(s)%t_veg_pa(ifp)-tfrz - 20.0_r8)/10.0_r8) 
                        ! kgC/s = kgN * kgC/kgN/s
                        currentCohort%livestem_mr  = live_stem_n * ED_val_base_mr_20 * tcwood * maintresp_reduction_factor
                     else
                        currentCohort%livestem_mr  = 0._r8
                     end if
                     
                     
                     ! Fine Root MR  (kgC/plant/s)
                     ! ------------------------------------------------------------------
                     currentCohort%froot_mr = 0._r8
                     do j = 1,bc_in(s)%nlevsoil
                        tcsoi  = q10**((bc_in(s)%t_soisno_sl(j)-tfrz - 20.0_r8)/10.0_r8)
                        currentCohort%froot_mr = currentCohort%froot_mr + &
                              fnrt_n * ED_val_base_mr_20 * tcsoi * currentPatch%rootfr_ft(ft,j) * maintresp_reduction_factor
                     enddo
                     
                     ! Coarse Root MR (kgC/plant/s) (below ground sapwood)
                     ! ------------------------------------------------------------------
                     if (woody(ft) == 1) then
                        currentCohort%livecroot_mr = 0._r8
                        do j = 1,bc_in(s)%nlevsoil
                           ! Soil temperature used to adjust base rate of MR
                           tcsoi  = q10**((bc_in(s)%t_soisno_sl(j)-tfrz - 20.0_r8)/10.0_r8)
                           currentCohort%livecroot_mr = currentCohort%livecroot_mr + &
                                 live_croot_n * ED_val_base_mr_20 * tcsoi * &
                                 currentPatch%rootfr_ft(ft,j) * maintresp_reduction_factor
                        enddo
                     else
                        currentCohort%livecroot_mr = 0._r8    
                     end if


                     ! ------------------------------------------------------------------
                     ! Part IX: Perform some unit conversions (rate to integrated) and
                     ! calcualate some fluxes that are sums and nets of the base fluxes
                     ! ------------------------------------------------------------------
                     
                     if ( debug ) write(fates_log(),*) 'EDPhoto 904 ', currentCohort%resp_m
                     if ( debug ) write(fates_log(),*) 'EDPhoto 905 ', currentCohort%rdark
                     if ( debug ) write(fates_log(),*) 'EDPhoto 906 ', currentCohort%livestem_mr
                     if ( debug ) write(fates_log(),*) 'EDPhoto 907 ', currentCohort%livecroot_mr
                     if ( debug ) write(fates_log(),*) 'EDPhoto 908 ', currentCohort%froot_mr
                        
                    

                     ! add on whole plant respiration values in kgC/indiv/s-1  
                     currentCohort%resp_m = currentCohort%livestem_mr + &
                                            currentCohort%livecroot_mr + &
                                            currentCohort%froot_mr

                     ! no drought response right now.. something like:
                     ! resp_m = resp_m * (1.0_r8 - currentPatch%btran_ft(currentCohort%pft) * &
                     !                    EDPftvarcon_inst%resp_drought_response(ft))   

                     currentCohort%resp_m = currentCohort%resp_m + currentCohort%rdark
                     
                     ! convert from kgC/indiv/s to kgC/indiv/timestep       
                     currentCohort%resp_m        = currentCohort%resp_m  * dtime
                     currentCohort%gpp_tstep     = currentCohort%gpp_tstep * dtime
                     currentCohort%ts_net_uptake = currentCohort%ts_net_uptake * dtime

                     if ( debug ) write(fates_log(),*) 'EDPhoto 911 ', currentCohort%gpp_tstep
                     if ( debug ) write(fates_log(),*) 'EDPhoto 912 ', currentCohort%resp_tstep
                     if ( debug ) write(fates_log(),*) 'EDPhoto 913 ', currentCohort%resp_m
                     
                     currentCohort%resp_g     = EDPftvarcon_inst%grperc(ft) * &
                                                (max(0._r8,currentCohort%gpp_tstep - &
                                                currentCohort%resp_m))
                     currentCohort%resp_tstep = currentCohort%resp_m + &
                                                currentCohort%resp_g ! kgC/indiv/ts
                     currentCohort%npp_tstep  = currentCohort%gpp_tstep - &
                                                currentCohort%resp_tstep  ! kgC/indiv/ts
                     
                     ! Accumulate the combined conductance (stomatal+leaf boundary layer)
                     ! Note that currentCohort%g_sb_laweight is weighted by the leaf area 
                     ! of each cohort and has units of [m/s] * [m2 leaf]

                     g_sb_leaves  = g_sb_leaves + currentCohort%g_sb_laweight
                     
                     ! Accumulate the total effective leaf area from all cohorts
                     ! in this patch. Normalize by canopy area outside the loop
                     check_elai   = check_elai  + cohort_eleaf_area
                     
                     currentCohort => currentCohort%shorter
                     
                  enddo  ! end cohort loop.   
               end if !count_cohorts is more than zero.
               
               check_elai = check_elai / currentPatch%total_canopy_area
               elai       = calc_areaindex(currentPatch,'elai')

               ! Normalize canopy total conductance by the effective LAI
               ! The value here was integrated over each cohort x leaf layer
               ! and was weighted by m2 of effective leaf area for each layer
               
               if(check_elai>tiny(check_elai)) then
                  
                  ! Normalize the leaf-area weighted canopy conductance
                  ! The denominator is the total effective leaf area in the canopy,
                  ! units of [m/s]*[m2] / [m2] = [m/s]
                  g_sb_leaves = g_sb_leaves / (elai*currentPatch%total_canopy_area)
                  
                  if( g_sb_leaves > (1._r8/rsmax0) ) then 
                     
                     ! Combined mean leaf resistance is the inverse of mean leaf conductance
                     r_sb_leaves  = 1.0_r8/g_sb_leaves
                     
                     if (r_sb_leaves<bc_in(s)%rb_pa(ifp)) then
                        write(fates_log(),*) 'Combined canopy resistance was somehow smaller than'
                        write(fates_log(),*) 'its boundary layer resistance component'
                        write(fates_log(),*) 'r_sb_leaves [s/m]: ',r_sb_leaves
                        write(fates_log(),*) 'bc_in(s)%rb_pa(ifp) [s/m]: ',bc_in(s)%rb_pa(ifp)
                        call endrun(msg=errMsg(sourcefile, __LINE__))
                     end if
                     
                     ! Mean leaf stomatal resistance for all patch leaves
                     r_stomata = (r_sb_leaves - bc_in(s)%rb_pa(ifp))
                                          
                  else
                     
                     ! Here we prevent super high resistances
                     ! and use a nominal value when conductance is low
                     r_stomata = rsmax0
                     
                  end if
                  
                  ! This will be multiplied by scaled by effective LAI in the host model
                  ! when it comes time to calculate a flux rate per unit ground
                  bc_out(s)%rssun_pa(ifp) = r_stomata
                  bc_out(s)%rssha_pa(ifp) = r_stomata
                  
                  ! This value is used for diagnostics, the molar form of conductance
                  ! is what is used in the field usually, so we track that form
                  currentPatch%c_stomata  = cf / r_stomata
                  
               else
                  
                  ! But this will prevent it from using an unintialized value
                  bc_out(s)%rssun_pa(ifp) = rsmax0
                  bc_out(s)%rssha_pa(ifp) = rsmax0

                  ! This value is used for diagnostics, the molar form of conductance
                  ! is what is used in the field usually, so we track that form
                  currentPatch%c_stomata  = cf / rsmax0
                  
               end if
               
               ! This value is used for diagnostics, the molar form of conductance
               ! is what is used in the field usually, so we track that form
               currentPatch%c_lblayer = cf / bc_in(s)%rb_pa(ifp)
               
            end if
            
            currentPatch => currentPatch%younger
            
         end do
         
      end do !site loop
      
    end associate
  end subroutine FatesPlantRespPhotosynthDrive
  
  ! =======================================================================================
  
  subroutine LeafLayerPhotosynthesis(f_sun_lsl,         &  ! in
                                     parsun_lsl,        &  ! in
                                     parsha_lsl,        &  ! in
                                     laisun_lsl,        &  ! in
                                     laisha_lsl,        &  ! in
                                     canopy_area_lsl,   &  ! in
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
                                     bbb,               &  ! in
                                     cf,                &  ! in
                                     gb_mol,            &  ! in
                                     ceair,             &  ! in
                                     mm_kco2,           &  ! in
                                     mm_ko2,            &  ! in
                                     co2_cpoint,        &  ! in
                                     lmr,               &  ! in
                                     psn_out,           &  ! out
                                     rstoma_out,        &  ! out
                                     anet_av_out)          ! out

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
    real(r8), intent(in) :: parsun_lsl        ! Absorbed PAR in sunlist leaves
    real(r8), intent(in) :: parsha_lsl        ! Absorved PAR in shaded leaves
    real(r8), intent(in) :: laisun_lsl        ! LAI in sunlit leaves
    real(r8), intent(in) :: laisha_lsl        ! LAI in shaded leaves
    real(r8), intent(in) :: canopy_area_lsl   ! 
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
   real(r8), intent(in) :: bbb             ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
   real(r8), intent(in) :: cf              ! s m**2/umol -> s/m (ideal gas conversion) [umol/m3]
   real(r8), intent(in) :: gb_mol          ! leaf boundary layer conductance (umol /m**2/s)
   real(r8), intent(in) :: ceair           ! vapor pressure of air, constrained (Pa)
   real(r8), intent(in) :: mm_kco2         ! Michaelis-Menten constant for CO2 (Pa)
   real(r8), intent(in) :: mm_ko2          ! Michaelis-Menten constant for O2 (Pa)
   real(r8), intent(in) :: co2_cpoint      ! CO2 compensation point (Pa)
   real(r8), intent(in) :: lmr             ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
   
   real(r8), intent(out) :: psn_out        ! carbon assimilated in this leaf layer umolC/m2/s
   real(r8), intent(out) :: rstoma_out     ! stomatal resistance (1/gs_lsl) (s/m)
   real(r8), intent(out) :: anet_av_out    ! net leaf photosynthesis (umol CO2/m**2/s) 
                                           ! averaged over sun and shade leaves.  

   ! Locals
   ! ------------------------------------------------------------------------
   integer :: c3c4_path_index    ! Index for which photosynthetic pathway 
                                 ! is active.  C4 = 0,  C3 = 1
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
   

   associate( bb_slope  => EDPftvarcon_inst%BB_slope)    ! slope of BB relationship

     ! photosynthetic pathway: 0. = c4, 1. = c3
     c3c4_path_index = nint(EDPftvarcon_inst%c3psn(ft))
     
     if (c3c4_path_index == 1) then
        init_co2_inter_c = init_a2l_co2_c3 * can_co2_ppress
     else
        init_co2_inter_c = init_a2l_co2_c4 * can_co2_ppress
     end if

     ! Part III: Photosynthesis and Conductance
     ! ----------------------------------------------------------------------------------
     
     if ( parsun_lsl <= 0._r8 ) then  ! night time

        anet_av_out = -lmr
        psn_out     = 0._r8
        rstoma_out  = min(rsmax0, 1._r8/bbb * cf)
        
     else ! day time (a little bit more complicated ...)
        
        !is there leaf area? - (NV can be larger than 0 with only stem area if deciduous)
        if ( laisun_lsl + laisha_lsl > 0._r8 ) then 

           !Loop aroun shaded and unshaded leaves          
           psn_out     = 0._r8    ! psn is accumulated across sun and shaded leaves. 
           rstoma_out  = 0._r8    ! 1/rs is accumulated across sun and shaded leaves. 
           anet_av_out = 0._r8
           gstoma      = 0._r8
           
           do  sunsha = 1,2      
              ! Electron transport rate for C3 plants. 
              ! Convert par from W/m2 to umol photons/m**2/s using the factor 4.6
              ! Convert from units of par absorbed per unit ground area to par 
              ! absorbed per unit leaf area. 
              
              if(sunsha == 1)then !sunlit
                 if(( laisun_lsl * canopy_area_lsl) > 0.0000000001_r8)then

                    qabs = parsun_lsl / (laisun_lsl * canopy_area_lsl )   

                    print*,"photo0", sunsha, qabs,parsun_lsl,laisun_lsl,canopy_area_lsl
                    qabs = qabs * 0.5_r8 * (1._r8 - fnps) *  4.6_r8 
                    
                 else
                    qabs = 0.0_r8
                    print*,"photo0", sunsha, qabs,parsun_lsl,laisun_lsl,canopy_area_lsl
                 end if
              else

                 qabs = parsha_lsl / (laisha_lsl * canopy_area_lsl)  
                 print*,"photo0", sunsha, qabs,parsun_lsl,laisun_lsl,canopy_area_lsl
                 qabs = qabs * 0.5_r8 * (1._r8 - fnps) *  4.6_r8 

              end if

              !convert the absorbed par into absorbed par per m2 of leaf, 
              ! so it is consistant with the vcmax and lmr numbers. 
              aquad = theta_psii
              bquad = -(qabs + jmax)
              cquad = qabs * jmax
              call quadratic_f (aquad, bquad, cquad, r1, r2)
              je = min(r1,r2)

             
              print*,"je: ",je,jmax

              ! Initialize intercellular co2
              co2_inter_c = init_co2_inter_c

              niter = 0
              loop_continue = .true.
              do while(loop_continue)                 
                 ! Increment iteration counter. Stop if too many iterations
                 niter = niter + 1
                 
                 ! Save old co2_inter_c
                 co2_inter_c_old = co2_inter_c
                 
                 ! Photosynthesis limitation rate calculations 
                 if (c3c4_path_index == 1)then    

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
                       !guard against /0's in the night.
                       if((laisun_lsl * canopy_area_lsl) > 0.0000000001_r8) then   
                          aj = quant_eff(c3c4_path_index) * parsun_lsl * 4.6_r8
                          !convert from per cohort to per m2 of leaf)
                          aj = aj / (laisun_lsl * canopy_area_lsl)
                       else
                          aj = 0._r8
                       end if
                    else
                       aj = quant_eff(c3c4_path_index) * parsha_lsl * 4.6_r8
                       aj = aj / (laisha_lsl * canopy_area_lsl)
                    end if

                    ! C4: PEP carboxylase-limited (CO2-limited)
                    ap = co2_rcurve_islope * max(co2_inter_c, 0._r8) / can_press  
                    
                 end if

                 ! Gross photosynthesis smoothing calculations. First co-limit ac and aj. Then co-limit ap
                 aquad = theta_cj(c3c4_path_index)
                 bquad = -(ac + aj)
                 cquad = ac * aj
                 call quadratic_f (aquad, bquad, cquad, r1, r2)
                 ai = min(r1,r2)

                 aquad = theta_ip
                 bquad = -(ai + ap)
                 cquad = ai * ap
                 call quadratic_f (aquad, bquad, cquad, r1, r2)
                 agross = min(r1,r2)

!                 agross = SmoothLimitCollatz1991(ac,aj,ap,c3c4_path_index)

                 ! Net carbon assimilation. Exit iteration if an < 0
                 anet = agross  - lmr
                 if (anet < 0._r8) then
                    loop_continue = .false.
                 end if

                 print*,"cc = ",co2_inter_c
                 print*,"anet(1) = ",anet

                 ! Quadratic gs_mol calculation with an known. Valid for an >= 0.
                 ! With an <= 0, then gs_mol = bbb
                 
                 leaf_co2_ppress = can_co2_ppress- 1.4_r8/gb_mol * anet * can_press 
                 leaf_co2_ppress = max(leaf_co2_ppress,1.e-06_r8)
                 aquad = leaf_co2_ppress
                 bquad = leaf_co2_ppress*(gb_mol - bbb) - bb_slope(ft) * anet * can_press
                 cquad = -gb_mol*(leaf_co2_ppress*bbb + &
                                  bb_slope(ft)*anet*can_press * ceair/ veg_esat )

                 call quadratic_f (aquad, bquad, cquad, r1, r2)
                 gs_mol = max(r1,r2)
                 
                 ! Derive new estimate for co2_inter_c
                 co2_inter_c = can_co2_ppress - anet * can_press * &
                       (h2o_co2_bl_diffuse_ratio*gs_mol + h2o_co2_stoma_diffuse_ratio*gb_mol) / (gb_mol*gs_mol)

                 ! Check for co2_inter_c convergence. Delta co2_inter_c/pair = mol/mol. 
                 ! Multiply by 10**6 to convert to umol/mol (ppm). Exit iteration if 
                 ! convergence criteria of +/- 1 x 10**-6 ppm is met OR if at least ten 
                 ! iterations (niter=10) are completed
                 
                 if ((abs(co2_inter_c-co2_inter_c_old)/can_press*1.e06_r8 <=  2.e-06_r8) &
                       .or. niter == 5) then
                    loop_continue = .false.
                 end if
              end do !iteration loop
              
              ! End of co2_inter_c iteration.  Check for an < 0, in which case gs_mol = bbb
              if (anet < 0._r8) then
                 gs_mol = bbb
              end if
              
              ! Final estimates for leaf_co2_ppress and co2_inter_c 
              ! (needed for early exit of co2_inter_c iteration when an < 0)
              leaf_co2_ppress = can_co2_ppress - 1.4_r8/gb_mol * anet * can_press
              leaf_co2_ppress = max(leaf_co2_ppress,1.e-06_r8)
              co2_inter_c = can_co2_ppress - anet * can_press * &
                            (h2o_co2_bl_diffuse_ratio*gs_mol + h2o_co2_stoma_diffuse_ratio*gb_mol) / (gb_mol*gs_mol)
              
              ! Convert gs_mol (umol /m**2/s) to gs (m/s) and then to rs (s/m)
              gs = gs_mol / cf
              
!              if ( debug ) write(fates_log(),*) 'EDPhoto 737 ', psn_out
!              if ( debug ) write(fates_log(),*) 'EDPhoto 738 ', agross
!              if ( debug ) write(fates_log(),*) 'EDPhoto 739 ', f_sun_lsl

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

!              if ( debug ) write(fates_log(),*) 'EDPhoto 758 ', psn_out
!              if ( debug ) write(fates_log(),*) 'EDPhoto 759 ', agross
!              if ( debug ) write(fates_log(),*) 'EDPhoto 760 ', f_sun_lsl
              
              ! Make sure iterative solution is correct
              if (gs_mol < 0._r8) then
                 write (fates_log(),*)'Negative stomatal conductance:'
                 write (fates_log(),*)'gs_mol= ',gs_mol
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if
              
              ! Compare with Ball-Berry model: gs_mol = m * an * hs/leaf_co2_ppress p + b
              hs = (gb_mol*ceair + gs_mol* veg_esat ) / ((gb_mol+gs_mol)*veg_esat )
              gs_mol_err = bb_slope(ft)*max(anet, 0._r8)*hs/leaf_co2_ppress*can_press + bbb
              
              if (abs(gs_mol-gs_mol_err) > 1.e-01_r8) then
                 write (fates_log(),*) 'CF: Ball-Berry error check - stomatal conductance error:'
                 write (fates_log(),*) gs_mol, gs_mol_err
              end if
              
           enddo !sunsha loop

           ! This is the stomatal resistance of the leaf layer
           rstoma_out = 1._r8/gstoma
           
        else
           !No leaf area. This layer is present only because of stems. 
           ! (leaves are off, or have reduced to 0)
           psn_out = 0._r8
           rstoma_out = min(rsmax0, 1._r8/bbb * cf)
           
        end if !is there leaf area? 
        
        
     end if    ! night or day 
   end associate
   return
  end subroutine LeafLayerPhotosynthesis
 
  ! =====================================================================================

  subroutine ScaleLeafLayerFluxToCohort(nv,          & ! in   currentCohort%nv
                                        psn_llz,     & ! in   %psn_z(1:currentCohort%nv,ft,cl)
                                        lmr_llz,     & ! in   lmr_z(1:currentCohort%nv,ft,cl)
                                        rs_llz,      & ! in   rs_z(1:currentCohort%nv,ft,cl)
                                        elai_llz,    & ! in   %elai_profile(cl,ft,1:currentCohort%nv)
                                        c_area,      & ! in   currentCohort%c_area
                                        nplant,      & ! in   currentCohort%n
                                        rb,          & ! in   bc_in(s)%rb_pa(ifp)
                                        maintresp_reduction_factor, & ! in 
                                        g_sb_laweight, & ! out  currentCohort%g_sb_laweight [m/s] [m2-leaf]
                                        gpp,         &   ! out  currentCohort%gpp_tstep
                                        rdark,       &   ! out  currentCohort%rdark
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
    real(r8), intent(in) :: c_area           ! crown area m2/m2
    real(r8), intent(in) :: nplant           ! indiv/m2
    real(r8), intent(in) :: rb               ! leaf boundary layer resistance (s/m)
    real(r8), intent(in) :: maintresp_reduction_factor  ! factor by which to reduce maintenance respiration
    real(r8), intent(out) :: g_sb_laweight      ! Combined conductance (stomatal + boundary layer) for the cohort 
                                        ! weighted by leaf area [m/s]*[m2]
    real(r8), intent(out) :: gpp        ! GPP (kgC/indiv/s)
    real(r8), intent(out) :: rdark      ! Dark Leaf Respiration (kgC/indiv/s)
    real(r8), intent(out) :: cohort_eleaf_area  ! Effective leaf area of the cohort [m2]
    
    ! GPP IN THIS SUBROUTINE IS A RATE. THE CALLING ARGUMENT IS GPP_TSTEP. AFTER THIS
    ! CALL THE RATE WILL BE MULTIPLIED BY THE INTERVAL TO GIVE THE INTEGRATED QUANT.
    
    ! Locals
    integer  :: il                       ! leaf layer index
    real(r8) :: cohort_layer_eleaf_area  ! the effective leaf area of the cohort's current layer [m2]
    
    cohort_eleaf_area = 0.0_r8
    g_sb_laweight             = 0.0_r8
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
       
       ! GPP    [umolC/m2leaf/s] * [m2 leaf ] -> [umolC/s]   (This is cohort group sum)
       gpp = gpp + psn_llz(il) * cohort_layer_eleaf_area
       
       ! Dark respiration
       ! [umolC/m2leaf/s] * [m2 leaf]    (This is the cohort group sum)
       rdark = rdark + lmr_llz(il) * cohort_layer_eleaf_area
       
    end do

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
  
  subroutine quadratic_f (a, b, c, r1, r2)
     !
     ! !DESCRIPTION:
     !==============================================================================!
     !----------------- Solve quadratic equation for its two roots -----------------!
     !==============================================================================!
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
     !
     ! !DESCRIPTION:
     !==============================================================================!
     !----------------- Solve quadratic equation for its two roots -----------------!
     ! THIS METHOD SIMPLY REMOVES THE DIV0 CHECK AND ERROR REPORTING                !
     !==============================================================================!
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
     !------------------------------------------------------------------------------
    
   !  if (a == 0._r8) then
   !     write (fates_log(),*) 'Quadratic solution error: a = ',a
   !     call endrun(msg=errMsg(sourcefile, __LINE__))
   !  end if
   
     if (b >= 0._r8) then
        q = -0.5_r8 * (b + sqrt(b*b - 4._r8*a*c))
     else
        q = -0.5_r8 * (b - sqrt(b*b - 4._r8*a*c))
     end if
   
     r1 = q / a
   !  if (q /= 0._r8) then
     r2 = c / q
   !  else
   !     r2 = 1.e36_r8
   !  end if
     
   end subroutine quadratic_fast


   ! ====================================================================================

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
      
      
      use EDTypesMod , only : ed_patch_type
      use EDTypesMod , only : ed_cohort_type

      ! Arguments
      type(ed_patch_type), target :: currentPatch
      type(ed_cohort_type), pointer :: currentCohort
      
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
                                     can_o2_ppress, &
                                     veg_tempk, &
                                     air_tempk, &
                                     air_vpress, &
                                     veg_esat,   &
                                     rb,        &
                                     mm_kco2,   &
                                     mm_ko2,    &
                                     mm_km,      &
                                     co2_cpoint, &
                                     cf,         &
                                     gb_mol, &
                                     ceair)
      
      ! ---------------------------------------------------------------------------------
      ! This subroutine calculates the specific Michaelis Menten Parameters (pa) for CO2
      ! and O2, as well as the CO2 compentation point.
      ! ---------------------------------------------------------------------------------
      
      use FatesConstantsMod, only: mmol_per_mol
      use FatesConstantsMod, only: umol_per_kmol
      use FatesConstantsMod, only : rgas => rgas_J_K_kmol

      ! Arguments
      real(r8), intent(in) :: can_press           ! Air pressure within the canopy (Pa)
      real(r8), intent(in) :: can_o2_ppress      ! Partial press of o2 in the canopy (Pa)
      real(r8), intent(in) :: veg_tempk           ! The temperature of the vegetation (K)
      real(r8), intent(in) :: air_tempk           ! Temperature of canopy air (K)
      real(r8), intent(in) :: air_vpress          ! Vapor pressure of canopy air (Pa)
      real(r8), intent(in) :: veg_esat            ! Saturated vapor pressure at veg surf (Pa)
      real(r8), intent(in) :: rb                  ! Leaf Boundary layer resistance (s/m)

      real(r8), intent(out) :: mm_kco2       ! Michaelis-Menten constant for CO2 (Pa)
      real(r8), intent(out) :: mm_ko2        ! Michaelis-Menten constant for O2 (Pa)
      real(r8), intent(out) :: mm_km         ! Michaelis-Menten combo term 
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
      ! except TPU from: Harley et al (1992) Plant, Cell and Environment 15:271-282

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
      cp25 = 0.5_r8 * can_o2_ppress / sco
      
      if( veg_tempk.gt.150_r8 .and. veg_tempk.lt.350_r8 )then
         mm_kco2       = kc25 * ft1_f(veg_tempk, kcha)
         mm_ko2         = ko25 * ft1_f(veg_tempk, koha)
         co2_cpoint     = cp25 * ft1_f(veg_tempk, cpha)
      else
         mm_kco2    = 1.0_r8
         mm_ko2     = 1.0_r8
         co2_cpoint = 1.0_r8
      end if
      
      ! Combined term from Michaeles Menten used in Farquar 1980 Gross
      ! photosynthesis calculation
      
      mm_km = mm_kco2 * (1._r8 + can_o2_ppress / mm_ko2 )

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
   
   subroutine LeafLayerMaintenanceRespiration(lmr25top_ft, &
                                              nscaler,   &
                                              ft,        &
                                              veg_tempk,     &
                                              lmr)

      use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
 
      ! Arguments
      real(r8), intent(in)  :: lmr25top_ft  ! canopy top leaf maint resp rate at 25C 
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


 


      ! Part I: Leaf Maintenance respiration: umol CO2 / m**2 [leaf] / s
      ! ----------------------------------------------------------------------------------
      lmr25 = lmr25top_ft * nscaler 
      
      if ( nint(EDpftvarcon_inst%c3psn(ft)) == 1)then
         lmr = lmr25 * ft1_f(veg_tempk, lmrha) * &
               fth_f(veg_tempk, lmrhd, lmrse, lmrc)
      else
         lmr = lmr25 * 2._r8**((veg_tempk-(tfrz+25._r8))/10._r8)
         lmr = lmr / (1._r8 + exp( 1.3_r8*(veg_tempk-(tfrz+55._r8)) ))
      end if
      
      ! Any hydrodynamic limitations could go here, currently none
      ! lmr = lmr * (nothing)
      
    end subroutine LeafLayerMaintenanceRespiration
     
   ! ====================================================================================

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

      use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm

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
      real(r8), intent(in) :: co2_rcurve_islope25top_ft      ! initial slope of CO2 response curve
                                              ! (C4 plants) at 25C, canopy top, this pft
      real(r8), intent(in) :: veg_tempk           ! vegetation temperature
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
      real(r8) :: vcmaxha        ! activation energy for vcmax (J/mol)
      real(r8) :: jmaxha         ! activation energy for jmax (J/mol)
      real(r8) :: tpuha          ! activation energy for tpu (J/mol)
      real(r8) :: vcmaxhd        ! deactivation energy for vcmax (J/mol)
      real(r8) :: jmaxhd         ! deactivation energy for jmax (J/mol)
      real(r8) :: tpuhd          ! deactivation energy for tpu (J/mol)
      real(r8) :: vcmaxse        ! entropy term for vcmax (J/mol/K)
      real(r8) :: jmaxse         ! entropy term for jmax (J/mol/K)
      real(r8) :: tpuse          ! entropy term for tpu (J/mol/K)
      real(r8) :: vcmaxc         ! scaling factor for high temperature inhibition (25 C = 1.0)
      real(r8) :: jmaxc          ! scaling factor for high temperature inhibition (25 C = 1.0)
      real(r8) :: tpuc           ! scaling factor for high temperature inhibition (25 C = 1.0)

      vcmaxha = EDPftvarcon_inst%vcmaxha(FT)
      jmaxha  = EDPftvarcon_inst%jmaxha(FT)
      tpuha   = EDPftvarcon_inst%tpuha(FT)
      
      vcmaxhd = EDPftvarcon_inst%vcmaxhd(FT)
      jmaxhd  = EDPftvarcon_inst%jmaxhd(FT)
      tpuhd   = EDPftvarcon_inst%tpuhd(FT)
      
      vcmaxse = EDPftvarcon_inst%vcmaxse(FT)
      jmaxse  = EDPftvarcon_inst%jmaxse(FT)
      tpuse   = EDPftvarcon_inst%tpuse(FT)

      vcmaxc = fth25_f(vcmaxhd, vcmaxse)
      jmaxc  = fth25_f(jmaxhd, jmaxse)
      tpuc   = fth25_f(tpuhd, tpuse)

      if ( parsun_lsl <= 0._r8) then           ! night time
         vcmax             = 0._r8
         jmax              = 0._r8
         tpu               = 0._r8
         co2_rcurve_islope = 0._r8
      else                                     ! day time
         vcmax25 = vcmax25top_ft * nscaler
         jmax25  = jmax25top_ft * nscaler
         tpu25   = tpu25top_ft * nscaler
         co2_rcurve_islope25 = co2_rcurve_islope25top_ft * nscaler
         
         ! Adjust for temperature
         vcmax = vcmax25 * ft1_f(veg_tempk, vcmaxha) * fth_f(veg_tempk, vcmaxhd, vcmaxse, vcmaxc)
         jmax  = jmax25 * ft1_f(veg_tempk, jmaxha) * fth_f(veg_tempk, jmaxhd, jmaxse, jmaxc)
         tpu   = tpu25 * ft1_f(veg_tempk, tpuha) * fth_f(veg_tempk, tpuhd, tpuse, tpuc)
         
         if (nint(EDPftvarcon_inst%c3psn(ft))  /=  1) then
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

    subroutine lowstorage_maintresp_reduction(frac, pft, maintresp_reduction_factor)

      ! This subroutine reduces maintenance respiration rates when storage pool is low.  The premise
      ! of this is that mortality of plants increases when storage is low because they are not able
      ! to repair tissues, generate defense compounds, etc.  This reduction is reflected in a reduced
      ! maintenance demand.  The output of this function takes the form of a curve between 0 and 1, 
      ! and the curvature of the function is determined by a parameter.

      ! Uses

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
          if ( EDPftvarcon_inst%maintresp_reduction_curvature(pft) .ne. 1._r8 ) then
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

    
    ! ===================================================================================

    ! Variable Definitions from MAAT (state)
    !
    ! environmental state
    !  oi = numeric(1),                 # atmospheric & internal O2  (kPa)
    !  ca = numeric(1),                 # atmospheric CO2            ( Pa)
    !  cb = numeric(1),                 # boundary layer CO2         ( Pa)
    !  ci = numeric(1),                 # leaf internal CO2          ( Pa) 
    !  cc = numeric(1),                 # chloroplast CO2            ( Pa)
    !  leaf_temp = numeric(1),          # leaf temperature           (oC)
      
    !  #leaf state - calculated by canopy object so need initialisation
    !  leafN_area = 2,                  # leaf N per unit area       (g N m-2)
    !  fwdw_ratio = 5,                  # fresh weight dry weight ratio, used for Sphagnum conductance term 
      
    !  #calculated state
    !  J   = numeric(1),                # electron transport rate                            (umol electrons m-2 s-1) 
    !  Acg = numeric(1),                # Carboxylaton limited rate of net asssimilation     (umol m-2 s-1)
    !  Ajg = numeric(1),                # light limited rate of carboxylation                (umol m-2 s-1)
    !  Apg = numeric(1),                # TPU limited rate of carboxylation                  (umol m-2 s-1)
    !  A   = numeric(1),                # actual rate of carboxylation                       (umol m-2 s-1)
    !  rd  = numeric(1),                # actual rate of respiration                         (umol m-2 s-1)
    !  lim = numeric(1),                # flag indicationg limitation state of assimilation, wc = wc limited, wj = wj limited, wp = wp limited
    !  # diagnostic state
    !  A_ana_rbzero   = numeric(1),     # rate of carboxylation assuming zero boundary layer resistance to CO2 diffusion (umol m-2 s-1)
    !  A_ana_rbg0zero = numeric(1),     # rate of carboxylation assuming zero boundary layer resistance & zero minimum stomatal conductance to CO2 diffusion (umol m-2 s-1)
    !  aguess    = numeric(4),          # value of three guesses and solution from semi-analytical solver  (umol m-2 s-1)
    !  faguess   = numeric(4),          # value of solver function on three guesses and solution from semi-analytical solver  (umol m-2 s-1)
    !  iter      = numeric(1),          # number of iterations to solve solver function in Brent uniroot
    !  estimprec = numeric(1),          # estimated precision of solve from uniroot solver function 
    !  assim     = numeric(2),          # roots of fitted quadratic in semi-analytical solver
    !  fA_ana_final   = numeric(2),     # value of solver function for final gusee from semi-analytical solver  (umol m-2 s-1)
    !  A_noR          = numeric(1),     # rate of carboxylation assuming zero resistance to CO2 diffusion (umol m-2 s-1)
    !  transition     = numeric(1)      # cc at the transition point where wc = wj                        (Pa)


    function quad_sol(a,b,c,root_type_in) result(root)
       
       real(r8), intent(in) :: a
       real(r8), intent(in) :: b
       real(r8), intent(in) :: c
       integer,  optional, intent(in) :: root_type_in
       
       integer  :: root_type
       real(r8) :: q
       real(r8) :: root1
       real(r8) :: root2
       real(r8) :: root

       if( present(root_type_in)) then
          root_type = root_type_in
       else
          root_type = lower_root
       end if

       q     = -0.5 * ( b + b/abs(b) * sqrt(b*b - 4*a*c) )

       root1 = q / a
       if (q /= 0._r8) then
          root2 = c / q
       else
          root2 = 1.e36_r8
       end if
       
       if(root_type.eq.upper_root) then
          root = max(root1,root2)
       elseif(root_type.eq.lower_root) then
          root = min(root1,root2)
       end if
       
     end function quad_sol

    ! ===================================================================================

     function SmoothLimitCollatz1991(a_cbox, a_rubp, a_prod, c3c4_path_index ) result(a_gross)

       ! --------------------------------------------------------------------------------
       ! Fine a smoothed solution of all three possible limiting states
       ! via Collatz et al. 1991
       ! --------------------------------------------------------------------------------

       real(r8), intent(in) :: a_cbox          ! Rubisco-limited (carboxylation) gross 
                                               ! photosynthesis (assimilate) (umol m-2 s-1)
       real(r8), intent(in) :: a_rubp          ! RuBP-limited gross photosynthesis    
                                               ! (umol m-2 s-1)
       real(r8), intent(in) :: a_prod          ! product-limited (C3) or CO2-limited  
                                               ! (umol m-2 s-1)
                                               ! (C4) gross photosynthesis
       integer,  intent(in) :: c3c4_path_index ! C3 = 1
                                               ! C4 = 0
       ! Output
       real(r8)             :: a_gross  ! Gross co-limited photosynthesis (umol m-2 s-1)

       real(r8) :: a, b, c    ! quadratic terms

       ! Find smooth solution for just light and carboxylation limitations first...
       
       a  = theta_cj(c3c4_path_index)
       b  = -(a_cbox + a_rubp)
       c  = a_cbox * a_rubp
       
       a_gross = quad_sol(a,b,c,lower_root)

       ! If product limited assimilation is available, combine that with the other two
       
       if( a_prod > 0._r8 ) then
          
          a =  theta_ip
          b =  -(a_prod + a_gross)
          c =  a_prod * a_gross
          
          a_gross = quad_sol(a,b,c,lower_root)
          
       end if
       
     end function SmoothLimitCollatz1991

    ! ===================================================================================
    
    function RStomaH2OMedlyn2011Fe(vpd) result(r_stomata_fe)
       
      ! f(e) component of rs from Medlyn 2011  
      
      real(r8), intent(in) :: vpd                ! Vapor pressure deficit (kPa)
      real(r8)             :: r_stomata_fe       ! Stomatal Conductance (molm-2 s-1)
      real(r8), parameter  :: g1_medlyn = 6._r8  ! Medlyn 2011 gs slope (kPa^0.5)
       
      print*,"VPD:",vpd

      r_stomata_fe = ( 1._r8 + g1_medlyn / sqrt(vpd) )

    end function RStomaH2OMedlyn2011Fe

    ! ===================================================================================

    function RStomaH2OMedlyn2011(vpd, can_co2_ppress, can_press, a_net) result(r_stomata)

      ! This function calculates the stomatal conductance to water vapor
      ! following Medlyn et al. 2011

      real(r8) :: vpd            ! vapor pressure deficit [kPa]
      real(r8) :: can_co2_ppress ! canopy co2 partial pressure [Pa]
      real(r8) :: can_press      ! atmospheric pressure in the canopy [Pa]
      real(r8) :: a_net          ! net assimilation rate [umol m-2 s-1]
      real(r8) :: r_stomata      ! Stomatal resistance to h2o [m2 s mol-1 h2o]
      real(r8), parameter :: g_stomata_min = 0.01 ! Medlyn 2011 minimum stomatal cond
                                                  ! [mol h2o m-2 s-1]
      
      real(r8) :: mol_frac       ! molar fraction of co2/air [umol/mol]
                                 ! this is equal to the partial pressure / total pressure
      
      mol_frac = can_co2_ppress / can_press * umol_per_mol
      

      r_stomata = 1.0_r8 / ( g_stomata_min + &
           RStomaH2OMedlyn2011Fe(vpd) * a_net / mol_frac)

      return
    end function RStomaH2OMedlyn2011

    ! ===================================================================================

    function FicksLawDiffusion(c_co2_in,a_net,press,resis) result(c_co2)

      ! This function is used to calculate CO2 partial pressure (concentration)
      ! of the chloroplast, intercellular space, or the boundary layer, based on 
      ! Ficks Law diffusion given a concentration on one side of a resistance pathway.

      real(r8) :: c_co2_in     ! known co2 concentration on other
                               ! side of boundary/resistor [Pa]
      real(r8) :: a_net        ! net assimilation        [umol m-2 s-1]
      real(r8) :: press        ! air pressure in the 
                               ! vicinity of the leaf      [Pa]
      real(r8) :: resis        ! resistance of the boundary
                               ! of interest               [m2 s mol-1]
                               ! could be either: 1 boundary layer resistance
                               !                  2 stomatal resistance
                               !       or         3 their series
      real(r8) :: c_co2        ! co2 concentration at the point of interest


      ! [Pa CO2] = [Pa CO2] - [umol CO2 m-2 s-1] * 
      !            [m2 s mol-1 air] * [Pa air] * [mol / umol]

      c_co2 = c_co2_in - a_net * resis * press * mol_per_umol

    end function FicksLawDiffusion

    ! ===================================================================================
    
    function AGrossFarquhar1980(vcmax,cplast_co2,mm_km) result(a_cbox)
      
      real(r8) :: vcmax  ! Maximum Carboxylation rate of RuBisCO [umol m-2 s-1]
                         ! functional on temp, and water stress
      
      real(r8) :: cplast_co2 ! Chloroplast CO2 partial pressure [Pa]

      real(r8) :: mm_km      ! Michaeles-Menten combined term that incorporates
                             ! the MM coefficients, the co2 compensation point
                             ! and the o2 partial pressure

      real(r8) :: a_cbox     ! Gross carboxylation limited assimilation 
                             ! [umol m-2 s-1 pa-1]
                             ! Normalized by chloroplastic co2 concentration

      a_cbox = vcmax / (cplast_co2 + mm_km)
      
    end function AGrossFarquhar1980

    ! ===================================================================================

    function AGrossRuBPLim(je,cplast_co2,co2_cpoint) result(a_rubp)
      
      !  Generic method calculate light limited photosynthetic rate 
      !  currently no other formulation exists (other than slighly 
      !  different parameters in the denominator)
           
      real(r8) :: je         ! electron transport rate [umol electrons m-2 s-1]
      real(r8) :: cplast_co2 ! chloroplast co2 concentration [Pa]
      real(r8) :: co2_cpoint ! CO2 compensation point [Pa]
           
      real(r8) :: a_rubp     ! Gross RuBp limited assimilation [umol m-2 s-1 pa-1]
                             ! Normalized by chloroplastic co2 concentration
      
      real(r8), parameter :: e_per_co2_fix = 4.0_r8  ! number of fully transported
                                                     ! e needed to fix 1 co2 molecule
      
      a_rubp = je / (4._r8*(cplast_co2 + 2._r8 * co2_cpoint))
      
      
    end function AGrossRuBPLim

    ! ===================================================================================
    
    function ERateLimFarquharWong1984(par,jmax) result(je)
      
      ! calculates the electron transport rate given Jmax
      ! I - irradiance & alpha - electrons transported per photon
      
      ! in von Caemmerer 2000    - qabs = par*a*(1-f)/2
      ! where a is leaf light absorptance, 1-f corrects for spectral 
      ! quality and light not absorbed by photosystems, and /2 accounts 
      ! for the 2 photons needed to fully tranport 1 electron 
      ! this is presumably the same as I * alpha
      ! if so von C 2000's alpha is 0.36125
      
      ! Convert par from W/m2 to umol photons/m**2/s using the factor 4.6
      
      real(r8) :: par                  ! Absorbed photosynthetically active 
                                       ! radiation by the leaf tissues [W/m2]
      real(r8) :: jmax                 !
      real(r8) :: qabs                 ! Absorbed photosynthetically active radiation
                                       ! by photosystem 2 [umol photons/m**2/s]
      real(r8) :: aquad, bquad, cquad  ! quadratic smooth terms
      real(r8) :: je                   ! electron transport rate 
    
      qabs = par * 0.5_r8 * (1._r8 - fnps) *  4.6_r8 
      
      aquad   = theta_psii
      bquad   = -(qabs + jmax)
      cquad   = qabs * jmax
      
      je = quad_sol(aquad,bquad,cquad)
      
    end function ERateLimFarquharWong1984

    ! ===================================================================================

    function AGrossTPU(tpu,cplast_co2,co2_cpoint ) result(a_prod)

      real(r8) :: cplast_co2 ! chloroplast co2 concentration [Pa]
      real(r8) :: co2_cpoint ! CO2 compensation point [Pa]
      real(r8) :: tpu        ! Triose phosphate utilization rate 
                             ! (umol CO2/m**2/s)

      real(r8) :: a_prod     ! Gross product limited assimilation
                             ! [umol m-2 s-1 pa-1]

      ! SHOULD THIS BE MULT BY 3?


      if (cplast_co2 < co2_cpoint) then
         a_prod = un_initialized
      else
         a_prod = tpu / ( cplast_co2-co2_cpoint)
      end if
    
    end function AGrossTPU


    ! ==============================================================================
    
    subroutine GetResidualFromAnet(a_net_in, r_blayer_co2, vpd, can_press, &
                                   can_co2_ppress, vcmax, je, tpu, mm_km,  &
                                   co2_cpoint, c3c4_path_index, lmr, a_net_resid, cplast_co2)
      
      ! combines A, rs, ri, ci & cc eqs to a single f(A), 
      ! combines all rate limiting processes
      !  -- for use with uniroot solver
      !  -- A is pased to this equation by the uniroot solver and is solved to find the root of this equation
      
      ! calculate cc from ca, rb, rs, and ri
      ! total resistance of a set of resistors in series is simply their sum 
      ! assumes boundary layer and stomatal resistance terms are in h2o units
      ! assumes mesophyll resistance is in co2 units
      
      real(r8),intent(in)  :: a_net_in
      real(r8),intent(in)  :: r_blayer_co2
      real(r8),intent(in)  :: vpd                ! Vapor Pressure Deficit [kpa]
      real(r8),intent(in)  :: can_press
      real(r8),intent(in)  :: can_co2_ppress
      real(r8),intent(in)  :: vcmax
      real(r8),intent(in)  :: je
      real(r8),intent(in)  :: tpu
      real(r8),intent(in)  :: mm_km
      real(r8),intent(in)  :: co2_cpoint
      integer,intent(in)   :: c3c4_path_index
      real(r8),intent(in)  :: lmr
      
      real(r8),intent(out) :: a_net_resid
      real(r8),intent(out) :: cplast_co2
      
      real(r8) :: r_stomata_co2
      real(r8) :: a_net_out  
      
      r_stomata_co2 = h2o_co2_stoma_diffuse_ratio * &
           RStomaH2OMedlyn2011(vpd, can_co2_ppress, can_press, a_net_in)
      
      ! Use ficks law to estimate the chloroplast CO2 concentration
      
      cplast_co2 = FicksLawDiffusion(can_co2_ppress,a_net_in, &
                                     can_press,r_stomata_co2+r_blayer_co2)


      ! Re-calculate co-limited CO2 assimilation with the updated co2 concentration
      call CoLimitedAssimilation(cplast_co2,vcmax,je,tpu,mm_km,co2_cpoint,lmr,c3c4_path_index,a_net_out)


      ! The residual is the difference, the re-calculated rate minus the input rate
      a_net_resid = a_net_out - a_net_in


      return
    end subroutine GetResidualFromAnet

    ! ===================================================================================
    
    subroutine CoLimitedAssimilation(cplast_co2,vcmax,je,tpu,mm_km, &
                                     co2_cpoint,lmr,c3c4_path_index,a_net_out)

      ! ---------------------------------------------------------------------------------
      ! DESCRIPTION
      ! ---------------------------------------------------------------------------------
      
      real(r8),intent(in) :: cplast_co2
      real(r8),intent(in) :: vcmax
      real(r8),intent(in) :: je
      
      real(r8),intent(in) :: tpu             ! Triose phosphate utilization rate 
                                             ! (umol CO2/m**2/s)
      real(r8),intent(in) :: mm_km           ! Michaelis-Menten O2 dependent
                                             ! combined term.
      real(r8),intent(in) :: co2_cpoint      ! CO2 compensation point (Pa)
      real(r8),intent(in) :: lmr             ! Leaf Maintenance Respiration (umol CO2 m-2 s-1)
      integer,intent(in)  :: c3c4_path_index ! C3 = 1, C4 = 0

      real(r8),intent(out) :: a_net_out      ! Net carbon assimilation (umol CO2 m-2 s-1)

      ! Locals
      real(r8) :: a_cbox         ! Rubisco-limited (carboxylation) gross 
                                 ! photosynthesis (assimilate) (umol m-2 s-1 pa-1)
      real(r8) :: a_rubp         ! RuBP-limited gross photosynthesis    
                                 ! (umol m-2 s-1 pa-1)
      real(r8) :: a_prod         ! product-limited (C3) or CO2-limited  
                                 ! (umol m-2 s-1 pa-1)
                                 ! (C4) gross photosynthesis
      real(r8) :: a_gross_min    ! Co-limited gross assimilation
                                 ! normalized by co2 

         
      ! Rubisco-limited (carboxylation) gross photosynthesis
      
      a_cbox = AGrossFarquhar1980(vcmax,cplast_co2,mm_km)
      
      ! RuBP-limited gross photosynthesis
      
      a_rubp = AGrossRuBPLim(je,cplast_co2,co2_cpoint) 

      ! product-limited (C3) or CO2-limited 
      
      a_prod = AGrossTPU(tpu,cplast_co2,co2_cpoint)
      
      ! Use quadratic smoothing to find a co-limited product
      a_gross_min = SmoothLimitCollatz1991(a_cbox, a_rubp, a_prod, c3c4_path_index )

      ! Convert to net assimilation [umol co2 m-2 s-1]
      a_net_out = a_gross_min * (cplast_co2 - co2_cpoint) - lmr


      return
    end subroutine CoLimitedAssimilation


    ! ===================================================================================

    subroutine AnalyticalZeroG(can_co2_ppress, vpd, vcmax, &
                               je, tpu, mm_km, co2_cpoint, lmr, c3c4_path_index, a_net )

       ! --------------------------------------------------------------------------------
       ! - finds the analytical solution assuming rb and rs are zero
       ! 
       ! Anthony Walker, 2018  (converted over from MAAT by Ryan Knox)
       !
       ! Original function: f_A_r_leaf_analytical_quad
       !
       ! For more on MAAT: https://github.com/walkeranthonyp/MAAT
       ! 
       ! ---------------------------------------------------------------------------------

       real(r8), intent(in) :: can_co2_ppress   ! co2 concentration of canopy air space (Pa)
       real(r8), intent(in) :: vpd             ! Vapor pressure deficit (kPa)
       real(r8), intent(in) :: vcmax           ! Temperature corrected vcmax
       real(r8), intent(in) :: je
       real(r8), intent(in) :: tpu
       real(r8), intent(in) :: mm_km
       real(r8), intent(in) :: co2_cpoint
       real(r8), intent(in) :: lmr
       integer,intent(in)   :: c3c4_path_index ! C3 = 1
                                               ! C4 = 0

       real(r8), intent(out) :: a_net          ! Net assimilation (umol m-2 s-1)


       ! Locals
       real(r8) :: r_stomata_fe 
       real(r8) :: cplast_co2


       ! Approximation of chloroplastic CO2
       ! --------------------------------------------------------------------------------

       print*,"vpd: ",vpd
       
       r_stomata_fe = RStomaH2OMedlyn2011Fe(vpd)

       print*,"rs:",r_stomata_fe
       
       if (r_stomata_fe < h2o_co2_stoma_diffuse_ratio) then
          write(fates_log(),*) 'Stomatal conductance environmental term is bogus'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       cplast_co2 = can_co2_ppress * (1._r8 - (h2o_co2_stoma_diffuse_ratio/r_stomata_fe))

       print*,"cc = ",cplast_co2
       
       call CoLimitedAssimilation(cplast_co2,vcmax,je,tpu,mm_km, &
                                  co2_cpoint,lmr,c3c4_path_index,a_net)


       ! WHY CALCULATE RS?
       ! .$state_pars$rs <- .$state$ca / (fe * Anet * .$env$atm_press*1e-6) 

       return
    end subroutine AnalyticalZeroG

    ! ===================================================================================

    subroutine SemiAnalyticalQuad(can_press,       &
                                  can_co2_ppress,  &
                                  mm_km,           & 
                                  vpd,             & 
                                  g_blayer_h2o,    &   ! aka gb_mol
                                  cf,              &
                                  vcmax,           & 
                                  jmax,            &
                                  tpu,             &
                                  co2_cpoint,      &
                                  lmr,             &
                                  parsun_lsl,      &
                                  parsha_lsl,      &
                                  laisun_lsl,      &
                                  laisha_lsl,      &
                                  canopy_area_lsl, &
                                  f_sun_lsl,       &
                                  ipft,            &
                                  a_net_out,       &
                                  psn_out,         &
                                  r_stomata_out)

      ! --------------------------------------------------------------------------------
      !  IF THIS PRODUCES A NEGATIVE A< THIS ROUTINE NEEDS TO BE CALLED
      !  WITH A REDUCED CONDUCTANCE.
      !
      ! - finds the analytical solution assuming rb and ri are zero to use as first 
      !   guess (a1)
      ! - make a second guess (a2) by calculating Cc using a1 and actual values of
      !   rb and ri   
      ! - make a third guess (a3) by taking the mean of a1 and a2  
      ! - calculate residual function value for these three guesses (fa1, fa2, fa3) 
      ! - check fa1 and fa2 span the root
      ! - if not, sequentially add 0.01 to the third guess (a3) until fa1 and fa3 
      !   span the root
      ! - fit a quadratic through the three sets of co-ordinates to find the root 
      ! 
      ! Anthony Walker, 2018  (converted over from MAAT w/ Ryan Knox)
      !
      ! Original function: f_A_r_leaf_semiana_quad <- function(.) {}
      !
      ! For more on MAAT: https://github.com/walkeranthonyp/MAAT
      ! 
      ! Carbon Assimilation is the solved quantity
      !
      !
      ! TO-DO: IS THERE ANY PRE-PROCESSING THAT NEEDS TO HAPPEN TO THE 
      ! FINAL CONDUCTANCE AFTER A_NET IS DETERMINED?
      !
      ! TO-DO: CHECK TPU...
      ! ---------------------------------------------------------------------------------
      
      ! Input Arguments
      real(r8),intent(in) :: can_press       ! Canopy atmospheric Pressure [Pa]
      real(r8),intent(in) :: can_co2_ppress  ! Canopy co2 partial pressure [Pa]
      real(r8),intent(in) :: mm_km           ! Michaelis-Menten O2 dependent
                                             ! combined term.
      real(r8),intent(in) :: vpd             ! Vapor pressure deficit [Pa]
      real(r8),intent(in) :: g_blayer_h2o    ! Leaf boundary layer conductance of 
                                             ! H2O  [mol h2o m-2 s-1]
      real(r8),intent(in) :: cf              ! ideal gas conversion factor between
                                             ! molar and velocity form conductance
      real(r8),intent(in) :: vcmax           !
      real(r8),intent(in) :: jmax            !
      real(r8),intent(in) :: tpu             ! Triose phosphate utilization rate 
                                             ! (umol CO2/m**2/s)
      real(r8),intent(in) :: co2_cpoint      ! CO2 compensation point (Pa)
      real(r8),intent(in) :: lmr             ! Leaf Maintenance Respiration (umol CO2 m-2 s-1)
      real(r8),intent(in) :: parsun_lsl      !
      real(r8),intent(in) :: parsha_lsl      !
      real(r8),intent(in) :: laisun_lsl      !
      real(r8),intent(in) :: laisha_lsl      !
      real(r8),intent(in) :: canopy_area_lsl !
      real(r8),intent(in) :: f_sun_lsl       ! fraction of leaves exposed to direct sun
      integer,intent(in)  :: ipft            ! PFT Index

      ! Output Arguments
      real(r8),intent(out) :: a_net_out       ! Net carbon assimilation (umol CO2 m-2 s-1)
      real(r8)             :: a_net_resid_out ! Weighted sun/shaded residual of a_net_out's
      real(r8),intent(out) :: psn_out         ! carbon assimilated in this leaf layer umolC/m2/s
      real(r8),intent(out) :: r_stomata_out   ! stomatal resistance (1/gs_lsl) (s/m)

      ! Locals
      integer  :: c3c4_path_index            ! C3 = 1, C4  = 0
      integer  :: sunsha                     ! counting index for sunlit/shaded calcs
      real(r8) :: qabs                       ! absorbed PAR by leaf tissues [W/m2]
      real(r8) :: je                         ! electron transport rate (umol electrons/m**2/s)
      real(r8) :: r_stomata_h2o              ! stomatal resistance to h2o [m2 s mol-1 h2o]
      real(r8) :: r_stomata_co2              ! stomatal resistance to co2 [m2 s mol-1 co2]
      real(r8) :: g_stomata_h2o              ! weighted average of final stomatal conductance
                                             ! of across sunlit/shaded leaves [m/s]
      real(r8) :: r_blayer_co2               ! boundary layer resistance to co2 [mol co2 m-2 s-1]
      real(r8) :: cplast_co2                 ! CO2 concentration in the chloroplast [Pa]
      real(r8) :: a_net(3)                   ! Net photosynthesis (assimilation) for the
                                             ! 3 intermediate steps (umol m-2 s-1)
      real(r8) :: a_net_resid(3)             ! Residual of each 3 net photosynthesis
                                             ! calculations (umol m-2 s-1)
      real(r8) :: a_net_final                ! final net photosynthesis of each sun/shade
      real(r8) :: a_net_temp                 ! temp a_net for final residual 
      real(r8) :: r1,r2                      ! temporary, roots from quadratic
      real(r8) :: ax,bx,cx                   ! intermediate smoothing terms 
                                             ! via Muller's method
      
      ! photosynthetic pathway: 0. = c4, 1. = c3
      c3c4_path_index = nint(EDPftvarcon_inst%c3psn(ipft))


      ! Solve for trivial case ... no light!
      ! ---------------------------------------------------------------------------------
      if ( parsun_lsl <= nearzero ) then  ! night time
         a_net_out       = -lmr
         psn_out         = 0._r8
         a_net_resid_out = un_initialized
         r_stomata_out   = min(rsmax0, &
                               RStomaH2OMedlyn2011(vpd, can_co2_ppress, can_press, 0._r8)   * cf)
         return
      end if
      
      ! Solve for trivial case ... no leaves!
      ! ---------------------------------------------------------------------------------
      if ( (laisun_lsl + laisha_lsl) <= nearzero ) then 
         a_net_out       = 0._r8
         a_net_resid_out = un_initialized
         psn_out         = 0._r8
         ! This stomatal resistance is really just
         ! an initialization. It will not effect
         ! canopy conductance since it will be multiplied
         ! by zero leaf area.
         r_stomata_out  = min(rsmax0, &
              RStomaH2OMedlyn2011(vpd, can_co2_ppress, can_press, 0._r8)  * cf)
         return
      end if
      
      ! ---------------------------------------------------------------------------------
      ! Non-trival case.  Leaves and light exist.
      ! Photosynthesis calculations are performed on both the
      ! sunlit and shaded portions of this leaf layer
      ! ---------------------------------------------------------------------------------

      ! Convert leaf boundary layer resistance units to
      ! diffusion of CO2 from H2O
      r_blayer_co2  = h2o_co2_bl_diffuse_ratio * 1._r8/g_blayer_h2o


      ! These values are accumulated over the sunlit and shaded portions. Initialize 0
      psn_out         = 0._r8
      a_net_out       = 0._r8
      g_stomata_h2o   = 0._r8
      a_net_resid_out = 0._r8
      
      do  sunsha = 1,2      

         ! Initialize the net photosynthesis intermediates
         a_net(1:3)       = un_initialized
         a_net_resid(1:3) = un_initialized
         

         ! First Step
         ! Convert from units of par absorbed per unit ground area to par 
         ! absorbed per unit leaf area. 
         
         if(sunsha == 1)then !sunlit
            if(( laisun_lsl * canopy_area_lsl) > 0.0000000001_r8)then
               qabs = parsun_lsl / (laisun_lsl * canopy_area_lsl )   
            else
               qabs = 0._r8
            end if
         else
            if ((laisha_lsl * canopy_area_lsl) > 0.0000000001_r8) then
               qabs = parsha_lsl / (laisha_lsl * canopy_area_lsl)  
            else
               qabs = 0._r8
            end if
         end if

         if(qabs<nearzero) then
             psn_out         = psn_out  + 0._r8
             a_net_out       = a_net_out - lmr*f_sun_lsl
             a_net_resid_out = a_net_resid_out + 0._r8
             g_stomata_h2o   = g_stomata_h2o + f_sun_lsl * &
                  (1._r8 / (RStomaH2OMedlyn2011(vpd, can_co2_ppress, can_press,0._r8) * cf))
            cycle    ! Go to 
         end if


         ! Update the electron transport rate (light limited) for C3 plants
         ! Determine electron transport rate, this does not need to be iterated, as
         ! it is not dependent on the chloroplast CO2 concentration or stomatal
         ! resistances.
         
         je = ERateLimFarquharWong1984(qabs,jmax)

         print*,"photo1", sunsha, qabs,parsun_lsl,laisun_lsl,canopy_area_lsl
         print*,"je: ",je,jmax

         ! ---------------------------------------------------------------------------------
         ! Get first guess (or upper bound for net photosynthesis)
         ! ---------------------------------------------------------------------------------

         ! find the analytical solution assuming rb and ri are zero to use as first guess (a1)
         ! rb is boundary layer resistance
         ! ri is internal resistance 
         ! umol Co2 m-2 s-1 


         call AnalyticalZeroG(can_co2_ppress, vpd, vcmax, &
                              je, tpu, mm_km, co2_cpoint, lmr, c3c4_path_index, a_net(1) )

      
         
         print*,"anet(1): ",a_net(1)

         if( a_net(1) < 0._r8 ) then
            psn_out         = psn_out  + a_net(1)
            a_net_out       = a_net_out + a_net(1)*f_sun_lsl
            a_net_resid_out = a_net_resid_out + un_initialized
            g_stomata_h2o   = g_stomata_h2o + f_sun_lsl * &
                 (1._r8 / (RStomaH2OMedlyn2011(vpd, can_co2_ppress, can_press, 0._r8) * cf))
            cycle
         end if
         

         ! Determine the convergence residual of this first guess (residual umol CO2 m-2 s-1 )

         call GetResidualFromAnet(a_net(1), r_blayer_co2, vpd, can_press, &
                                  can_co2_ppress, vcmax, je, tpu, mm_km, co2_cpoint, &
                                  c3c4_path_index, lmr, a_net_resid(1), cplast_co2)

         
         if(abs(a_net_resid(1)) < 1.e-8_r8 ) then
             psn_out         = psn_out  + a_net(1)
             a_net_out       = a_net_out + a_net(1)*f_sun_lsl
             a_net_resid_out = a_net_resid_out + a_net_resid(1)
             g_stomata_h2o   = g_stomata_h2o + f_sun_lsl * &
                  (1._r8 / (RStomaH2OMedlyn2011(vpd, can_co2_ppress, can_press, a_net(1)) * cf))
            cycle
         end if
      
         ! ---------------------------------------------------------------------------------
         ! Get second guess (or lower bound that brackets solution
         ! ---------------------------------------------------------------------------------
         
      
         ! Calculate the stomatal resistance of co2
         ! Stomatal conductance submodels report in resistance to h2o, so
         ! convert from resistance of h2o to co2

         a_net(2) = a_net(1) + a_net_resid(1)
         
!!         r_stomata_h2o =  RStomaH2OMedlyn2011(vpd, can_co2_ppress, can_press, a_net(1))
!!         r_stomata_co2 =  r_stomata_h2o * h2o_co2_stoma_diffuse_ratio 
             
         ! Use ficks law to estimate the chloroplast CO2 concentration
         
!!         cplast_co2 = FicksLawDiffusion(can_co2_ppress, &
!!                                        a_net(1),       &
!!                                        can_press,      &
!!                                        r_stomata_co2 + r_blayer_co2)
         

         ! Update our assimilation rate based on our new cplast_co2 [umol m-2 s-1]
         
!!         call CoLimitedAssimilation(cplast_co2,vcmax,je,tpu,mm_km, &
!!                                    co2_cpoint,lmr,c3c4_path_index,a_net(2))
         

         if( a_net(2) < 0._r8 ) then
            write(fates_log(),*) ' this was unexpected'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         
         ! Calculate the 2nd residual
         
         call GetResidualFromAnet(a_net(2), r_blayer_co2, vpd, can_press, &
                                  can_co2_ppress, vcmax, je, tpu, mm_km, co2_cpoint, &
                                  c3c4_path_index, lmr, a_net_resid(2), cplast_co2)

         ! if both a1 and a2 are below zero then 
         ! return -1 and leafsys routine will calculate A assuming gs = g0
         if( (a_net_resid(1)*a_net_resid(2)) > 0.0_r8 ) then
            write(fates_log(),*) 'The semi-analytical photosynthesis solver did not generate'
            write(fates_log(),*) ' reasonable starting bounds to the solution. Exiting'
            write(fates_log(),*) ' a_net[1] = ',a_net(1)
            write(fates_log(),*) ' a_net_resid[1] = ',a_net_resid(1)
            write(fates_log(),*) ' a_net[2] = ',a_net(2)
            write(fates_log(),*) ' a_net_resid[2] = ',a_net_resid(2)
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
      
         ! ---------------------------------------------------------------------------------
         ! Calculate the 3rd bracketed guess
         ! ---------------------------------------------------------------------------------
         
         ! Calculate the mean of net photosynthesis of the two edge cases
         a_net(3) = 0.5_r8 * (a_net(1) + a_net(2))

         call GetResidualFromAnet(a_net(3), r_blayer_co2, vpd, can_press, &
                                  can_co2_ppress, vcmax, je, tpu, mm_km, co2_cpoint, &
                                  c3c4_path_index, lmr, a_net_resid(3), cplast_co2)
      
         ! ---------------------------------------------------------------------------------
         ! fit a quadratic through the three sets of co-ordinates a la Lomas 
         ! (one step of Muller's method) 
         ! Muller, David E., "A Method for Solving Algebraic Equations Using an Automatic 
         ! Computer," Mathematical Tables and Other Aids to Computation, 10 (1956)    
         ! bx <- ((a1+a0)*(fa2-fa1)/(a2-a1)-(a2+a1)*(fa1-fa0)/(a1-a0))/(a0-a2)
         ! ax <- (fa2-fa1-bx*(a2-a1))/(a2**2-a1**2)
         ! cx <- fa1-bx*a1-ax*a1**2
         ! ---------------------------------------------------------------------------------
         
         bx = ( (a_net(2)+a_net(1))*(a_net_resid(3)-a_net_resid(2)) / (a_net(3)-a_net(2)) - &
              (a_net(3)+a_net(2))*(a_net_resid(2)-a_net_resid(1))/(a_net(2)-a_net(1)) ) / & 
              (a_net(1)-a_net(3))
         
         ax = ( a_net_resid(3) - a_net_resid(2) - bx * (a_net(3)-a_net(2))) / &
              (a_net(3)**2._r8 - a_net(2)**2._r8)
         cx = a_net_resid(2)-bx*a_net(2)-ax*a_net(2)**2._r8
      
         ! find the root of the quadratic that lies between guesses 
         call quadratic_f(ax,bx,cx,r1,r2)
         
         a_net_final = max(r1,r2)

         ! Calculate stomatal conductance baesd on final a_net
         r_stomata_co2 = h2o_co2_stoma_diffuse_ratio * &
              RStomaH2OMedlyn2011(vpd, can_co2_ppress, can_press, a_net_final)
      
         ! Use ficks law to estimate the chloroplast CO2 concentration
         cplast_co2 = FicksLawDiffusion(can_co2_ppress,a_net_final,can_press,r_stomata_co2+r_blayer_co2)

         
         ! Re-calculate co-limited CO2 assimilation one last time
         ! to estimate the residual
         call CoLimitedAssimilation(cplast_co2,vcmax,je,tpu,mm_km, &
                                    co2_cpoint,lmr,c3c4_path_index,a_net_temp)
         
         
         ! The residual is the difference, the re-calculated rate minus the input rate
         a_net_resid = a_net_temp - a_net_final

         
         ! Update averages of the sun/shaded weighted 
         ! photosynthes, net assimilation and conductance output variables
         
         if(sunsha == 1)then
            psn_out         = psn_out  + (a_net_final+lmr) * f_sun_lsl
            a_net_out       = a_net_out + a_net_final * f_sun_lsl
            a_net_resid_out = a_net_resid_out + (a_net_temp - a_net_final) * f_sun_lsl
            g_stomata_h2o   = g_stomata_h2o + f_sun_lsl * &
                 (1._r8 / (RStomaH2OMedlyn2011(vpd, can_co2_ppress, can_press, a_net_out) * cf))
         else
            psn_out         = psn_out + (a_net_final+lmr) * (1._r8-f_sun_lsl)
            a_net_out       = a_net_out + a_net_final * (1._r8-f_sun_lsl) 
            a_net_resid_out = a_net_resid_out + (a_net_temp - a_net_final) * (1._r8-f_sun_lsl)
            g_stomata_h2o   = g_stomata_h2o + (1._r8-f_sun_lsl) * &
                 (1._r8 / (RStomaH2OMedlyn2011(vpd, can_co2_ppress, can_press, a_net_out) * cf))
         end if
         
         
      end do
      
      
      ! The resistance passed back is in units [s/m].
      r_stomata_out = 1._r8 / g_stomata_h2o
      
      
      return
    end subroutine SemiAnalyticalQuad
    


 end module FATESPlantRespPhotosynthMod
