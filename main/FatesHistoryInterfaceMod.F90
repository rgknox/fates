module FatesHistoryInterfaceMod

  use FatesConstantsMod        , only : r8 => fates_r8
  use FatesConstantsMod        , only : fates_avg_flag_length
  use FatesConstantsMod        , only : fates_short_string_length
  use FatesConstantsMod        , only : fates_long_string_length
  use FatesConstantsMod        , only : itrue,ifalse
  use FatesConstantsMod        , only : calloc_abs_error
  use FatesConstantsMod        , only : mg_per_kg
  use FatesConstantsMod        , only : pi_const
  use FatesGlobals             , only : fates_log
  use FatesGlobals             , only : endrun => fates_endrun
  use EDTypesMod               , only : nclmax
  use EDTypesMod               , only : ican_upper
  use PRTGenericMod            , only : element_pos
  use PRTGenericMod            , only : num_elements
  use EDTypesMod               , only : site_fluxdiags_type
  use EDtypesMod               , only : ed_site_type
  use EDtypesMod               , only : ed_cohort_type
  use EDtypesMod               , only : ed_patch_type  
  use EDtypesMod               , only : AREA
  use EDtypesMod               , only : AREA_INV
  use EDTypesMod               , only : numWaterMem
  use EDTypesMod               , only : num_vegtemp_mem
  use EDTypesMod               , only : site_massbal_type
  use PRTGenericMod            , only : element_list
  use EDTypesMod               , only : N_DIST_TYPES
  use EDTypesMod               , only : dtype_ifall
  use EDTypesMod               , only : dtype_ifire
  use EDTypesMod               , only : dtype_ilog
  use FatesIODimensionsMod     , only : fates_io_dimension_type
  use FatesIOVariableKindMod   , only : fates_io_variable_kind_type
  use FatesHistoryVariableType , only : fates_history_variable_type
  use FatesInterfaceTypesMod        , only : hlm_hio_ignore_val
  use FatesInterfaceTypesMod        , only : hlm_use_planthydro
  use FatesInterfaceTypesMod        , only : hlm_use_ed_st3
  use FatesInterfaceTypesMod        , only : hlm_use_cohort_age_tracking
  use FatesInterfaceTypesMod        , only : numpft
  use FatesInterfaceTypesMod        , only : hlm_freq_day
  use FatesInterfaceTypesMod        , only : hlm_parteh_mode
  use EDParamsMod              , only : ED_val_comp_excln
  use EDParamsMod              , only : ED_val_phen_coldtemp
  use FatesInterfaceTypesMod        , only : nlevsclass, nlevage
  use FatesInterfaceTypesMod        , only : nlevheight
  use FatesInterfaceTypesMod        , only : bc_in_type
  use FatesInterfaceTypesMod        , only : hlm_model_day
  use FatesInterfaceTypesMod        , only : nlevcoage

  ! FIXME(bja, 2016-10) need to remove CLM dependancy 
  use EDPftvarcon              , only : EDPftvarcon_inst
  use PRTParametersMod         , only : prt_params
  
  ! CIME Globals
  use shr_log_mod              , only : errMsg => shr_log_errMsg
  use shr_infnan_mod           , only : isnan => shr_infnan_isnan
  use FatesConstantsMod        , only : g_per_kg
  use FatesConstantsMod        , only : ha_per_m2
  use FatesConstantsMod        , only : days_per_sec
  use FatesConstantsMod        , only : sec_per_day
  use FatesConstantsMod        , only : days_per_year
  use FatesConstantsMod        , only : years_per_day
  use FatesLitterMod           , only : litter_type
  use FatesConstantsMod        , only : secondaryforest

  use PRTGenericMod            , only : leaf_organ, fnrt_organ, sapw_organ
  use PRTGenericMod            , only : struct_organ, store_organ, repro_organ
  use PRTGenericMod            , only : all_carbon_elements
  use PRTGenericMod            , only : carbon12_element
  use PRTGenericMod            , only : nitrogen_element, phosphorus_element
  use PRTGenericMod            , only : prt_carbon_allom_hyp

  implicit none
  private          ! By default everything is private

  ! These variables hold the index of the history output structure so we don't
  ! have to constantly do name lookup when we want to populate the dataset
  ! These indices are set during "define_history_vars()" call to "set_history_var()"
  ! during the initialize phase.  Definitions are not provided, for an explanation of
  ! the variable go to its registry.  (IH_ signifies "index history")
  !
  ! Because of the complex sub-gridscale structure of FATES, in which multiple patches and cohorts
  ! exist within a gridcell, along with vertical gradients within and between canopy layers, as well
  ! as distinct classes such as PFTs or fuel size bins, there are multiple different dimensions in
  ! which it is possible to output history variables to better understand what's going on.
  !
  ! a key point is that, while the number of patches or cohorts can in principle be large, and 
  ! the age and size indices of a given patch or cohort can be finely resolved, we collapse these 
  ! continuously varying indices into bins of time-invariant width for the purposes of history 
  ! outputting.  This is because a given patch or cohort may not persist across a given interval
  ! of history averaging, so it is better to output all patches of cohorts whose index is within 
  ! a given interval along the size or age bin.
  !
  ! Another particularity of the issue of FATES shifting its subgrid structure frequently 
  ! and possibly having multiple (or zero) patches or cohorts within a given bin is that, if you
  ! want to output an average quantities across some dimension, such as a mean carbon flux across 
  ! patch area of a given age, in general it is better to output both the numerator and denominator
  ! of the averaging calculation separately, rather than the average itself, and then calculate 
  ! the average in post-processing. So, e.g. this means outputting both the patch area and the 
  ! product of the flux within each patch and the patch area as separate variables.  Doing this 
  ! allows conservation even when the weights are changing rapidly and simplifies the logic when
  ! the number of patches or cohorts may be anywhere from zero to a large number.
  !
  ! So what this means is that anything that is disaggregated at the patch area requires 
  ! outputting the patch age distribution (in units of patch area / site area) as the denominator
  ! of the average and then calculating the numerator of the average as XXX times the patch 
  ! area so (so in units of XXX * patch area / site area). For cohort-level quantities,
  ! this requires outputting the number density (in units of individuals per site area), etc.
  !
  ! For reference, some standardized abbreviations of the FATES dimensions are listed here:
  ! sz = size-class dimension
  ! cacls = cohort age-class dimension
  ! pft  = the pft dimension
  ! age  = the age bin dimension
  ! height = the height bin dimension
  ! cwdsc  = the coarse woody debris size class dimension
  ! 
  ! Since the netcdf interface can only handle variables with a certain number of dimensions,
  ! we have create some "multiplexed" dimensions that combine two or more dimensions into a
  ! single dimension.  Examples of these are the following:
  ! scpf = size class x PFT
  ! cacpf = cohort age class x PFT
  ! cnlf = canopy layer x leaf layer
  ! cnlfpft = canopy layer x leaf layer x PFT
  ! scag = size class bin x age bin
  ! scagpft = size class bin x age bin x PFT
  ! agepft  = age bin x PFT
 

  ! A recipe for adding a new history variable to this module:
  ! (1) decide what time frequency it makes sense to update the variable at, and what dimension(s)
  !     you want to output the variable on
  ! (2) add the ih_ integer variable in the immediately following section of the module.  
  !     use the suffix as outlined above for the dimension you are using.
  ! (3) define a corresponding hio_ variable by associating it to the ih_ variable 
  !     in the associate section of the subroutine that corresponds to the time-updating 
  !     frequency that you've chosen
  !     (i.e. if half-hourly, then work in subroutine update_history_prod; if daily, 
  !     then work in subroutine update_history_dyn)
  ! (4) within that subroutine, add the logic that passes the information from the 
  !     fates-native variable (possibly on a patch or cohort structure) to the history 
  !     hio_ variable that you've associated to.
  ! (5) add the variable name, metadata, units, dimension, updating frequency, the ih_ variable 
  !     index, etc via a call to the set_history_var method in the subroutine define_history_vars.
  !
  
  ! Indices to 1D Patch variables

  integer :: ih_storec
  integer :: ih_leafc
  integer :: ih_sapwc
  integer :: ih_fnrtc
  integer :: ih_reproc
  integer :: ih_totvegc

  integer :: ih_storen
  integer :: ih_leafn
  integer :: ih_sapwn
  integer :: ih_fnrtn
  integer :: ih_repron
  integer :: ih_totvegn

  integer :: ih_storep
  integer :: ih_leafp
  integer :: ih_sapwp
  integer :: ih_fnrtp
  integer :: ih_reprop
  integer :: ih_totvegp

  integer :: ih_nuptake
  integer :: ih_puptake
  integer :: ih_cefflux
  integer :: ih_nefflux
  integer :: ih_pefflux
  integer :: ih_nneedgrow
  integer :: ih_nneedmax
  integer :: ih_pneedgrow
  integer :: ih_pneedmax
  
  integer :: ih_trimming
  integer :: ih_area_plant
  integer :: ih_area_trees

  integer :: ih_cwd_elcwd

  integer :: ih_litter_in    ! carbon only
  integer :: ih_litter_out   ! carbon only
  integer :: ih_seed_bank    ! carbon only
  integer :: ih_seeds_in     ! carbon only

  integer :: ih_litter_in_elem
  integer :: ih_litter_out_elem
  integer :: ih_seed_bank_elem
  integer :: ih_seeds_in_local_elem
  integer :: ih_seeds_in_extern_elem
  integer :: ih_seed_decay_elem
  integer :: ih_seed_germ_elem

  integer :: ih_fines_ag_elem
  integer :: ih_fines_bg_elem
  integer :: ih_cwd_ag_elem
  integer :: ih_cwd_bg_elem
  integer :: ih_burn_flux_elem

  ! Size-class x PFT mass states
  
  integer :: ih_bstor_canopy_scpf
  integer :: ih_bstor_understory_scpf
  integer :: ih_bleaf_canopy_scpf
  integer :: ih_bleaf_understory_scpf



  integer :: ih_totvegn_scpf
  integer :: ih_leafn_scpf
  integer :: ih_fnrtn_scpf
  integer :: ih_storen_scpf
  integer :: ih_sapwn_scpf
  integer :: ih_repron_scpf
  integer :: ih_nuptake_scpf
  integer :: ih_nefflux_scpf
  integer :: ih_nneedgrow_scpf
  integer :: ih_nneedmax_scpf

  integer :: ih_totvegc_scpf
  integer :: ih_leafc_scpf
  integer :: ih_fnrtc_scpf
  integer :: ih_storec_scpf
  integer :: ih_sapwc_scpf
  integer :: ih_reproc_scpf
  integer :: ih_cefflux_scpf

  integer :: ih_totvegp_scpf
  integer :: ih_leafp_scpf
  integer :: ih_fnrtp_scpf
  integer :: ih_reprop_scpf
  integer :: ih_storep_scpf
  integer :: ih_sapwp_scpf
  integer :: ih_puptake_scpf
  integer :: ih_pefflux_scpf
  integer :: ih_pneedgrow_scpf
  integer :: ih_pneedmax_scpf

  integer :: ih_daily_temp
  integer :: ih_daily_rh
  integer :: ih_daily_prec
 
  integer :: ih_agb
  integer :: ih_npp
  integer :: ih_gpp
  integer :: ih_aresp
  integer :: ih_maint_resp
  integer :: ih_growth_resp
  integer :: ih_ar_canopy
  integer :: ih_gpp_canopy
  integer :: ih_ar_understory
  integer :: ih_gpp_understory
  integer :: ih_canopy_biomass
  integer :: ih_understory_biomass

  integer :: ih_primaryland_fusion_error
  integer :: ih_disturbance_rate_p2p
  integer :: ih_disturbance_rate_p2s
  integer :: ih_disturbance_rate_s2s
  integer :: ih_fire_disturbance_rate
  integer :: ih_logging_disturbance_rate
  integer :: ih_fall_disturbance_rate
  integer :: ih_potential_disturbance_rate
  integer :: ih_harvest_carbonflux

  ! Indices to site by size-class by age variables
  integer :: ih_nplant_scag
  integer :: ih_nplant_canopy_scag
  integer :: ih_nplant_understory_scag
  integer :: ih_ddbh_canopy_scag
  integer :: ih_ddbh_understory_scag
  integer :: ih_mortality_canopy_scag
  integer :: ih_mortality_understory_scag

  ! Indices to site by size-class by age by pft variables
  integer :: ih_nplant_scagpft

  ! Indices to site by patch age by pft variables
  integer :: ih_biomass_agepft
  integer :: ih_npp_agepft
  integer :: ih_scorch_height_agepft

  ! Indices to (site) variables

  integer :: ih_nep

  integer :: ih_c_stomata
  integer :: ih_c_lblayer

  integer :: ih_fire_c_to_atm


  integer :: ih_cbal_err_fates
  integer :: ih_err_fates

  integer :: ih_npatches
  integer :: ih_ncohorts
  integer :: ih_demotion_carbonflux
  integer :: ih_promotion_carbonflux
  integer :: ih_canopy_mortality_carbonflux
  integer :: ih_understory_mortality_carbonflux
  integer :: ih_canopy_spread
  integer :: ih_npp_leaf
  integer :: ih_npp_seed
  integer :: ih_npp_stem
  integer :: ih_npp_froot
  integer :: ih_npp_croot
  integer :: ih_npp_stor
  integer :: ih_leaf_mr
  integer :: ih_froot_mr
  integer :: ih_livestem_mr
  integer :: ih_livecroot_mr
  integer :: ih_fraction_secondary_forest
  integer :: ih_biomass_secondary_forest
  integer :: ih_woodproduct
  integer :: ih_h2oveg
  integer :: ih_h2oveg_dead
  integer :: ih_h2oveg_recruit
  integer :: ih_h2oveg_growturn_err
  integer :: ih_h2oveg_pheno_err
  integer :: ih_h2oveg_hydro_err
  
  integer :: ih_cstatus
  integer :: ih_dstatus
  integer :: ih_gdd
  integer :: ih_nchilldays
  integer :: ih_ncolddays
  integer :: ih_cleafoff
  integer :: ih_cleafon
  integer :: ih_dleafoff
  integer :: ih_dleafon
  integer :: ih_meanliqvol

  integer :: ih_nesterov_fire_danger
  integer :: ih_fire_nignitions
  integer :: ih_fire_fdi
  integer :: ih_fire_intensity_area_product
  integer :: ih_spitfire_ros
  integer :: ih_fire_ros_area_product
  integer :: ih_effect_wspeed
  integer :: ih_tfc_ros
  integer :: ih_tfc_ros_area_product
  integer :: ih_fire_intensity
  integer :: ih_fire_area
  integer :: ih_fire_fuel_bulkd
  integer :: ih_fire_fuel_eff_moist
  integer :: ih_fire_fuel_sav
  integer :: ih_fire_fuel_mef
  integer :: ih_sum_fuel
  integer :: ih_fragmentation_scaler

  integer :: ih_nplant_scpf
  integer :: ih_gpp_scpf
  integer :: ih_npp_totl_scpf
  integer :: ih_npp_leaf_scpf
  integer :: ih_npp_seed_scpf
  integer :: ih_npp_fnrt_scpf
  integer :: ih_npp_bgsw_scpf
  integer :: ih_npp_bgdw_scpf
  integer :: ih_npp_agsw_scpf
  integer :: ih_npp_agdw_scpf
  integer :: ih_npp_stor_scpf
  
 
  integer :: ih_mortality_canopy_scpf
  integer :: ih_mortality_understory_scpf
  integer :: ih_nplant_canopy_scpf
  integer :: ih_nplant_understory_scpf
  integer :: ih_ddbh_canopy_scpf
  integer :: ih_ddbh_understory_scpf
  integer :: ih_gpp_canopy_scpf
  integer :: ih_gpp_understory_scpf
  integer :: ih_ar_canopy_scpf
  integer :: ih_ar_understory_scpf

  integer :: ih_ddbh_scpf
  integer :: ih_growthflux_scpf
  integer :: ih_growthflux_fusion_scpf
  integer :: ih_ba_scpf
  integer :: ih_agb_scpf
  integer :: ih_m1_scpf
  integer :: ih_m2_scpf
  integer :: ih_m3_scpf
  integer :: ih_m4_scpf
  integer :: ih_m5_scpf
  integer :: ih_m6_scpf
  integer :: ih_m7_scpf  
  integer :: ih_m8_scpf
  integer :: ih_m9_scpf
  integer :: ih_m10_scpf
  integer :: ih_crownfiremort_scpf
  integer :: ih_cambialfiremort_scpf

  integer :: ih_m10_capf
  integer :: ih_nplant_capf

  integer :: ih_ar_scpf
  integer :: ih_ar_grow_scpf
  integer :: ih_ar_maint_scpf
  integer :: ih_ar_darkm_scpf
  integer :: ih_ar_agsapm_scpf
  integer :: ih_ar_crootm_scpf
  integer :: ih_ar_frootm_scpf
  
  integer :: ih_c13disc_scpf

  ! indices to (site x size [size class bins]) variables
  integer :: ih_ba_sz
  integer :: ih_nplant_sz
  integer :: ih_nplant_canopy_sz
  integer :: ih_nplant_understory_sz
  integer :: ih_lai_canopy_sz
  integer :: ih_lai_understory_sz
  integer :: ih_sai_canopy_sz
  integer :: ih_sai_understory_sz
  integer :: ih_mortality_canopy_sz
  integer :: ih_mortality_understory_sz
  integer :: ih_demotion_rate_sz
  integer :: ih_promotion_rate_sz
  integer :: ih_trimming_canopy_sz
  integer :: ih_trimming_understory_sz
  integer :: ih_crown_area_canopy_sz
  integer :: ih_crown_area_understory_sz
  integer :: ih_ddbh_canopy_sz
  integer :: ih_ddbh_understory_sz
  integer :: ih_agb_sz
  integer :: ih_biomass_sz

  ! mortality vars
  integer :: ih_m1_sz
  integer :: ih_m2_sz
  integer :: ih_m3_sz
  integer :: ih_m4_sz
  integer :: ih_m5_sz
  integer :: ih_m6_sz
  integer :: ih_m7_sz  
  integer :: ih_m8_sz
  integer :: ih_m9_sz
  integer :: ih_m10_sz

  integer :: ih_m10_cacls
  integer :: ih_nplant_cacls

  ! lots of non-default diagnostics for understanding canopy versus understory carbon balances
  integer :: ih_rdark_canopy_sz
  integer :: ih_livestem_mr_canopy_sz
  integer :: ih_livecroot_mr_canopy_sz
  integer :: ih_froot_mr_canopy_sz
  integer :: ih_resp_g_canopy_sz
  integer :: ih_resp_m_canopy_sz
  integer :: ih_leaf_md_canopy_sz
  integer :: ih_root_md_canopy_sz
  integer :: ih_carbon_balance_canopy_sz
  integer :: ih_bstore_md_canopy_sz
  integer :: ih_bdead_md_canopy_sz
  integer :: ih_bsw_md_canopy_sz
  integer :: ih_seed_prod_canopy_sz
  integer :: ih_npp_leaf_canopy_sz
  integer :: ih_npp_fnrt_canopy_sz
  integer :: ih_npp_sapw_canopy_sz
  integer :: ih_npp_dead_canopy_sz
  integer :: ih_npp_seed_canopy_sz
  integer :: ih_npp_stor_canopy_sz

  integer :: ih_rdark_understory_sz
  integer :: ih_livestem_mr_understory_sz
  integer :: ih_livecroot_mr_understory_sz
  integer :: ih_froot_mr_understory_sz
  integer :: ih_resp_g_understory_sz
  integer :: ih_resp_m_understory_sz
  integer :: ih_leaf_md_understory_sz
  integer :: ih_root_md_understory_sz
  integer :: ih_carbon_balance_understory_sz
  integer :: ih_bsw_md_understory_sz
  integer :: ih_bdead_md_understory_sz
  integer :: ih_bstore_md_understory_sz
  integer :: ih_seed_prod_understory_sz
  integer :: ih_npp_leaf_understory_sz
  integer :: ih_npp_fnrt_understory_sz
  integer :: ih_npp_sapw_understory_sz
  integer :: ih_npp_dead_understory_sz
  integer :: ih_npp_seed_understory_sz
  integer :: ih_npp_stor_understory_sz

  integer :: ih_yesterdaycanopylevel_canopy_sz
  integer :: ih_yesterdaycanopylevel_understory_sz

  ! indices to (site x pft) variables
  integer :: ih_biomass_pft
  integer :: ih_leafbiomass_pft
  integer :: ih_storebiomass_pft
  integer :: ih_nindivs_pft
  integer :: ih_recruitment_pft
  integer :: ih_mortality_pft
  integer :: ih_crownarea_pft
  integer :: ih_canopycrownarea_pft

  ! indices to (site x patch-age) variables
  integer :: ih_area_age
  integer :: ih_lai_age
  integer :: ih_canopy_area_age
  integer :: ih_gpp_age
  integer :: ih_npp_age
  integer :: ih_ncl_age
  integer :: ih_npatches_age
  integer :: ih_zstar_age
  integer :: ih_biomass_age
  integer :: ih_c_stomata_age
  integer :: ih_c_lblayer_age
  integer :: ih_agesince_anthrodist_age
  integer :: ih_secondaryforest_area_age
  integer :: ih_area_burnt_age
  ! integer :: ih_fire_rate_of_spread_front_age
  integer :: ih_fire_intensity_age
  integer :: ih_fire_sum_fuel_age

  ! indices to (site x height) variables
  integer :: ih_canopy_height_dist_height
  integer :: ih_leaf_height_dist_height

  ! Indices to hydraulics variables
  
  integer :: ih_errh2o_scpf
  integer :: ih_tran_scpf

!  integer :: ih_h2osoi_scagpft  ! hijacking the scagpft dimension instead of creating a new shsl dimension
  integer :: ih_sapflow_scpf
  integer :: ih_sapflow
  integer :: ih_iterh1_scpf          
  integer :: ih_iterh2_scpf           
  integer :: ih_supsub_scpf              
  integer :: ih_ath_scpf               
  integer :: ih_tth_scpf               
  integer :: ih_sth_scpf                     
  integer :: ih_lth_scpf                     
  integer :: ih_awp_scpf                     
  integer :: ih_twp_scpf  
  integer :: ih_swp_scpf                     
  integer :: ih_lwp_scpf  
  integer :: ih_aflc_scpf                     
  integer :: ih_tflc_scpf  
  integer :: ih_sflc_scpf                     
  integer :: ih_lflc_scpf                   
  integer :: ih_btran_scpf
  
  ! Hydro: Soil water states
  integer :: ih_rootwgt_soilvwc
  integer :: ih_rootwgt_soilvwcsat
  integer :: ih_rootwgt_soilmatpot

  ! Hydro: Soil water state by layer
  integer :: ih_soilmatpot_sl
  integer :: ih_soilvwc_sl
  integer :: ih_soilvwcsat_sl
  
  ! Hydro: Root water Uptake rates
  integer :: ih_rootuptake
  integer :: ih_rootuptake_sl
  integer :: ih_rootuptake0_scpf
  integer :: ih_rootuptake10_scpf
  integer :: ih_rootuptake50_scpf
  integer :: ih_rootuptake100_scpf

  
  ! indices to (site x fuel class) variables
  integer :: ih_litter_moisture_fuel
  integer :: ih_burnt_frac_litter_fuel
  integer :: ih_fuel_amount_fuel

  ! indices to (site x cwd size class) variables
  integer :: ih_cwd_ag_cwdsc
  integer :: ih_cwd_bg_cwdsc
  integer :: ih_cwd_ag_in_cwdsc
  integer :: ih_cwd_bg_in_cwdsc
  integer :: ih_cwd_ag_out_cwdsc
  integer :: ih_cwd_bg_out_cwdsc

  ! indices to (site x [canopy layer x leaf layer]) variables
  integer :: ih_parsun_z_cnlf
  integer :: ih_parsha_z_cnlf
  integer :: ih_laisun_z_cnlf
  integer :: ih_laisha_z_cnlf
  integer :: ih_fabd_sun_cnlf
  integer :: ih_fabd_sha_cnlf
  integer :: ih_fabi_sun_cnlf
  integer :: ih_fabi_sha_cnlf
  integer :: ih_ts_net_uptake_cnlf
  integer :: ih_crownarea_cnlf
  integer :: ih_parprof_dir_cnlf
  integer :: ih_parprof_dif_cnlf

  ! indices to (site x [canopy layer x leaf layer x pft]) variables
  integer :: ih_parsun_z_cnlfpft
  integer :: ih_parsha_z_cnlfpft
  integer :: ih_laisun_z_cnlfpft
  integer :: ih_laisha_z_cnlfpft
  integer :: ih_fabd_sun_cnlfpft
  integer :: ih_fabd_sha_cnlfpft
  integer :: ih_fabi_sun_cnlfpft
  integer :: ih_fabi_sha_cnlfpft
  integer :: ih_parprof_dir_cnlfpft
  integer :: ih_parprof_dif_cnlfpft

  ! indices to (site x canopy layer) variables
  integer :: ih_parsun_top_can
  integer :: ih_parsha_top_can
  integer :: ih_laisun_top_can
  integer :: ih_laisha_top_can
  integer :: ih_fabd_sun_top_can
  integer :: ih_fabd_sha_top_can
  integer :: ih_fabi_sun_top_can
  integer :: ih_fabi_sha_top_can
  integer :: ih_crownarea_can

  ! The number of variable dim/kind types we have defined (static)

  integer, parameter, public :: fates_history_num_dimensions = 50
  integer, parameter, public :: fates_history_num_dim_kinds = 50

  ! This structure is allocated by thread, and must be calculated after the FATES
  ! sites are allocated, and their mapping to the HLM is identified.  This structure
  ! is not combined with iovar_bounds, because that one is multi-instanced.  This
  ! structure is used more during the update phase, wherease _bounds is used
  ! more for things like flushing
  type, public :: iovar_map_type
     integer, allocatable :: site_index(:)   ! maps site indexes to the HIO site position
     integer, allocatable :: patch1_index(:) ! maps site index to the HIO
                                             ! (deprecated, but required for coupling)
  end type iovar_map_type


  type, public :: fates_history_interface_type
     
     ! Instance of the list of history output varialbes
     type(fates_history_variable_type), allocatable :: hvars(:)
     integer, private :: num_history_vars_
     
     ! Instanteat one registry of the different dimension/kinds (dk)
     ! All output variables will have a pointer to one of these dk's
     type(fates_io_variable_kind_type) :: dim_kinds(fates_history_num_dim_kinds)
     
     ! This is a structure that explains where FATES patch boundaries
     ! on each thread point to in the host IO array, this structure is
     ! allocated by number of threads. This could be dynamically
     ! allocated, but is unlikely to change...?
     type(fates_io_dimension_type) :: dim_bounds(fates_history_num_dimensions)
     
     type(iovar_map_type), pointer :: iovar_map(:)

   
     !! THESE WERE EXPLICITLY PRIVATE WHEN TYPE WAS PUBLIC
     integer, private :: patch_index_, column_index_, levgrnd_index_, levscpf_index_
     integer, private :: levscls_index_, levpft_index_, levage_index_
     integer, private :: levfuel_index_, levcwdsc_index_, levscag_index_
     integer, private :: levcan_index_, levcnlf_index_, levcnlfpft_index_
     integer, private :: levscagpft_index_, levagepft_index_
     integer, private :: levheight_index_
     integer, private :: levelem_index_, levelpft_index_
     integer, private :: levelcwd_index_, levelage_index_
     integer, private :: levcacls_index_, levcapf_index_

     
   contains
     
     procedure :: Init
     procedure :: SetThreadBoundsEach
     procedure :: initialize_history_vars
     procedure :: assemble_history_output_types
     
     procedure :: update_history_dyn
     procedure :: update_history_hifrq
     procedure :: update_history_hydraulics

     ! 'get' methods used by external callers to access private read only data

     procedure :: num_history_vars
     procedure :: patch_index
     procedure :: column_index
     procedure :: levgrnd_index
     procedure :: levscpf_index
     procedure :: levscls_index
     procedure :: levcapf_index
     procedure :: levcacls_index
     procedure :: levpft_index
     procedure :: levage_index
     procedure :: levfuel_index
     procedure :: levcwdsc_index
     procedure :: levcan_index
     procedure :: levcnlf_index
     procedure :: levcnlfpft_index
     procedure :: levscag_index
     procedure :: levscagpft_index
     procedure :: levagepft_index
     procedure :: levheight_index
     procedure :: levelem_index
     procedure :: levelpft_index
     procedure :: levelcwd_index
     procedure :: levelage_index

     ! private work functions
     procedure, private :: define_history_vars
     procedure, private :: set_history_var
     procedure, private :: init_dim_kinds_maps
     procedure, private :: set_dim_indices
     procedure, private :: flush_hvars

     procedure, private :: set_patch_index
     procedure, private :: set_column_index
     procedure, private :: set_levgrnd_index
     procedure, private :: set_levscpf_index
     procedure, private :: set_levcacls_index
     procedure, private :: set_levcapf_index
     procedure, private :: set_levscls_index
     procedure, private :: set_levpft_index
     procedure, private :: set_levage_index
     procedure, private :: set_levfuel_index
     procedure, private :: set_levcwdsc_index
     procedure, private :: set_levcan_index
     procedure, private :: set_levcnlf_index
     procedure, private :: set_levcnlfpft_index
     procedure, private :: set_levscag_index
     procedure, private :: set_levscagpft_index
     procedure, private :: set_levagepft_index
     procedure, private :: set_levheight_index
     
     procedure, private :: set_levelem_index
     procedure, private :: set_levelpft_index
     procedure, private :: set_levelcwd_index
     procedure, private :: set_levelage_index


  end type fates_history_interface_type
   
  character(len=*), parameter :: sourcefile = &
         __FILE__

contains

  ! ======================================================================
  
  subroutine Init(this, num_threads, fates_bounds)

    use FatesIODimensionsMod, only : patch, column, levgrnd, levscpf
    use FatesIODimensionsMod, only : levscls, levpft, levage
    use FatesIODimensionsMod, only : levcacls, levcapf
    use FatesIODimensionsMod, only : levfuel, levcwdsc, levscag
    use FatesIODimensionsMod, only : levscagpft, levagepft
    use FatesIODimensionsMod, only : levcan, levcnlf, levcnlfpft
    use FatesIODimensionsMod, only : fates_bounds_type
    use FatesIODimensionsMod, only : levheight
    use FatesIODimensionsMod, only : levelem, levelpft
    use FatesIODimensionsMod, only : levelcwd, levelage

    implicit none

    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: num_threads
    type(fates_bounds_type), intent(in) :: fates_bounds

    integer :: dim_count = 0

    dim_count = dim_count + 1
    call this%set_patch_index(dim_count)
    call this%dim_bounds(dim_count)%Init(patch, num_threads, &
         fates_bounds%patch_begin, fates_bounds%patch_end)

    dim_count = dim_count + 1
    call this%set_column_index(dim_count)
    call this%dim_bounds(dim_count)%Init(column, num_threads, &
         fates_bounds%column_begin, fates_bounds%column_end)

    dim_count = dim_count + 1
    call this%set_levgrnd_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levgrnd, num_threads, &
         fates_bounds%ground_begin, fates_bounds%ground_end)

    dim_count = dim_count + 1
    call this%set_levscpf_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levscpf, num_threads, &
         fates_bounds%sizepft_class_begin, fates_bounds%sizepft_class_end)

    dim_count = dim_count + 1
    call this%set_levscls_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levscls, num_threads, &
         fates_bounds%size_class_begin, fates_bounds%size_class_end)

    dim_count = dim_count + 1
    call this%set_levcacls_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcacls, num_threads, &
         fates_bounds%coage_class_begin, fates_bounds%coage_class_end)

    dim_count = dim_count + 1
    call this%set_levcapf_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcapf, num_threads, &
         fates_bounds%coagepf_class_begin, fates_bounds%coagepf_class_end)

    dim_count = dim_count + 1
    call this%set_levpft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levpft, num_threads, &
         fates_bounds%pft_class_begin, fates_bounds%pft_class_end)

    dim_count = dim_count + 1
    call this%set_levage_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levage, num_threads, &
         fates_bounds%age_class_begin, fates_bounds%age_class_end)

    dim_count = dim_count + 1
    call this%set_levfuel_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levfuel, num_threads, &
         fates_bounds%fuel_begin, fates_bounds%fuel_end)

    dim_count = dim_count + 1
    call this%set_levcwdsc_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcwdsc, num_threads, &
         fates_bounds%cwdsc_begin, fates_bounds%cwdsc_end)

    dim_count = dim_count + 1
    call this%set_levcan_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcan, num_threads, &
         fates_bounds%can_begin, fates_bounds%can_end)

    dim_count = dim_count + 1
    call this%set_levcnlf_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcnlf, num_threads, &
         fates_bounds%cnlf_begin, fates_bounds%cnlf_end)

    dim_count = dim_count + 1
    call this%set_levcnlfpft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcnlfpft, num_threads, &
         fates_bounds%cnlfpft_begin, fates_bounds%cnlfpft_end)

    dim_count = dim_count + 1
    call this%set_levscag_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levscag, num_threads, &
         fates_bounds%sizeage_class_begin, fates_bounds%sizeage_class_end)
    
    dim_count = dim_count + 1
    call this%set_levscagpft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levscagpft, num_threads, &
         fates_bounds%sizeagepft_class_begin, fates_bounds%sizeagepft_class_end)
    
    dim_count = dim_count + 1
    call this%set_levagepft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levagepft, num_threads, &
         fates_bounds%agepft_class_begin, fates_bounds%agepft_class_end)
    
    dim_count = dim_count + 1
    call this%set_levheight_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levheight, num_threads, &
         fates_bounds%height_begin, fates_bounds%height_end)

    dim_count = dim_count + 1
    call this%set_levelem_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levelem, num_threads, &
         fates_bounds%elem_begin, fates_bounds%elem_end)

    dim_count = dim_count + 1
    call this%set_levelpft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levelpft, num_threads, &
          fates_bounds%elpft_begin, fates_bounds%elpft_end)
    
    dim_count = dim_count + 1
    call this%set_levelcwd_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levelcwd, num_threads, &
         fates_bounds%elcwd_begin, fates_bounds%elcwd_end)

    dim_count = dim_count + 1
    call this%set_levelage_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levelage, num_threads, &
          fates_bounds%elage_begin, fates_bounds%elage_end)
    

    ! FIXME(bja, 2016-10) assert(dim_count == FatesHistorydimensionmod::num_dimension_types)

    ! Allocate the mapping between FATES indices and the IO indices
    allocate(this%iovar_map(num_threads))
    
  end subroutine Init

  ! ======================================================================
  subroutine SetThreadBoundsEach(this, thread_index, thread_bounds)

    use FatesIODimensionsMod, only : fates_bounds_type

    implicit none

    class(fates_history_interface_type), intent(inout) :: this

    integer, intent(in) :: thread_index
    type(fates_bounds_type), intent(in) :: thread_bounds

    integer :: index
    
    index = this%patch_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%patch_begin, thread_bounds%patch_end)

    index = this%column_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%column_begin, thread_bounds%column_end)

    index = this%levgrnd_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%ground_begin, thread_bounds%ground_end)

    index = this%levscpf_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%sizepft_class_begin, thread_bounds%sizepft_class_end)

    index = this%levscls_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%size_class_begin, thread_bounds%size_class_end)

    index = this%levcacls_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%coage_class_begin, thread_bounds%coage_class_end)

    index = this%levcapf_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%coagepf_class_begin, thread_bounds%coagepf_class_end)

    index = this%levpft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%pft_class_begin, thread_bounds%pft_class_end)
    
    index = this%levage_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%age_class_begin, thread_bounds%age_class_end)
    
    index = this%levfuel_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%fuel_begin, thread_bounds%fuel_end)
    
    index = this%levcwdsc_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cwdsc_begin, thread_bounds%cwdsc_end)
    
    index = this%levcan_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%can_begin, thread_bounds%can_end)
    
    index = this%levcnlf_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cnlf_begin, thread_bounds%cnlf_end)
    
    index = this%levcnlfpft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
          thread_bounds%cnlfpft_begin, thread_bounds%cnlfpft_end)
    
    index = this%levscag_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
          thread_bounds%sizeage_class_begin, thread_bounds%sizeage_class_end)
    
    index = this%levscagpft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
          thread_bounds%sizeagepft_class_begin, thread_bounds%sizeagepft_class_end)
    
    index = this%levagepft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
          thread_bounds%agepft_class_begin, thread_bounds%agepft_class_end)
    
    index = this%levheight_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
          thread_bounds%height_begin, thread_bounds%height_end)

    index = this%levelem_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%elem_begin, thread_bounds%elem_end)

    index = this%levelpft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%elpft_begin, thread_bounds%elpft_end)
    
    index = this%levelcwd_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%elcwd_begin, thread_bounds%elcwd_end)

    index = this%levelage_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%elage_begin, thread_bounds%elage_end)
    

    


    
  end subroutine SetThreadBoundsEach
  
  ! ===================================================================================
  subroutine assemble_history_output_types(this)

    use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
    use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_coage_r8, site_coage_pft_r8
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
    use FatesIOVariableKindMod, only : site_height_r8
    use FatesIOVariableKindMod, only : site_elem_r8, site_elpft_r8
    use FatesIOVariableKindMod, only : site_elcwd_r8, site_elage_r8

   implicit none

    class(fates_history_interface_type), intent(inout) :: this

    call this%init_dim_kinds_maps()

    call this%set_dim_indices(patch_r8, 1, this%patch_index())

    call this%set_dim_indices(site_r8, 1, this%column_index())

    call this%set_dim_indices(patch_ground_r8, 1, this%patch_index())
    call this%set_dim_indices(patch_ground_r8, 2, this%levgrnd_index())

    call this%set_dim_indices(site_ground_r8, 1, this%column_index())
    call this%set_dim_indices(site_ground_r8, 2, this%levgrnd_index())

    call this%set_dim_indices(patch_size_pft_r8, 1, this%patch_index())
    call this%set_dim_indices(patch_size_pft_r8, 2, this%levscpf_index())

    call this%set_dim_indices(site_size_pft_r8, 1, this%column_index())
    call this%set_dim_indices(site_size_pft_r8, 2, this%levscpf_index())

    call this%set_dim_indices(site_size_r8, 1, this%column_index())
    call this%set_dim_indices(site_size_r8, 2, this%levscls_index())

    call this%set_dim_indices(site_coage_r8, 1, this%column_index())
    call this%set_dim_indices(site_coage_r8, 2, this%levcacls_index())

    call this%set_dim_indices(site_coage_pft_r8, 1, this%column_index())
    call this%set_dim_indices(site_coage_pft_r8, 2, this%levcapf_index())

    call this%set_dim_indices(site_pft_r8, 1, this%column_index())
    call this%set_dim_indices(site_pft_r8, 2, this%levpft_index())

    call this%set_dim_indices(site_age_r8, 1, this%column_index())
    call this%set_dim_indices(site_age_r8, 2, this%levage_index())

    call this%set_dim_indices(site_fuel_r8, 1, this%column_index())
    call this%set_dim_indices(site_fuel_r8, 2, this%levfuel_index())

    call this%set_dim_indices(site_cwdsc_r8, 1, this%column_index())
    call this%set_dim_indices(site_cwdsc_r8, 2, this%levcwdsc_index())

    call this%set_dim_indices(site_can_r8, 1, this%column_index())
    call this%set_dim_indices(site_can_r8, 2, this%levcan_index())

    call this%set_dim_indices(site_cnlf_r8, 1, this%column_index())
    call this%set_dim_indices(site_cnlf_r8, 2, this%levcnlf_index())

    call this%set_dim_indices(site_cnlfpft_r8, 1, this%column_index())
    call this%set_dim_indices(site_cnlfpft_r8, 2, this%levcnlfpft_index())

    call this%set_dim_indices(site_scag_r8, 1, this%column_index())
    call this%set_dim_indices(site_scag_r8, 2, this%levscag_index())

    call this%set_dim_indices(site_scagpft_r8, 1, this%column_index())
    call this%set_dim_indices(site_scagpft_r8, 2, this%levscagpft_index())

    call this%set_dim_indices(site_agepft_r8, 1, this%column_index())
    call this%set_dim_indices(site_agepft_r8, 2, this%levagepft_index())

    call this%set_dim_indices(site_height_r8, 1, this%column_index())
    call this%set_dim_indices(site_height_r8, 2, this%levheight_index())

    call this%set_dim_indices(site_elem_r8, 1, this%column_index())
    call this%set_dim_indices(site_elem_r8, 2, this%levelem_index())
    
    call this%set_dim_indices(site_elpft_r8, 1, this%column_index())
    call this%set_dim_indices(site_elpft_r8, 2, this%levelpft_index())

    call this%set_dim_indices(site_elcwd_r8, 1, this%column_index())
    call this%set_dim_indices(site_elcwd_r8, 2, this%levelcwd_index())
    
    call this%set_dim_indices(site_elage_r8, 1, this%column_index())
    call this%set_dim_indices(site_elage_r8, 2, this%levelage_index())
    

  end subroutine assemble_history_output_types
  
  ! ===================================================================================
  
  subroutine set_dim_indices(this, dk_name, idim, dim_index)

    use FatesIOVariableKindMod , only : iotype_index

    implicit none

    ! arguments
    class(fates_history_interface_type), intent(inout) :: this
    character(len=*), intent(in)     :: dk_name
    integer, intent(in)              :: idim  ! dimension index
    integer, intent(in) :: dim_index


    ! local
    integer :: ityp

    ityp = iotype_index(trim(dk_name), fates_history_num_dim_kinds, this%dim_kinds)

    ! First check to see if the dimension is allocated
    if (this%dim_kinds(ityp)%ndims < idim) then
       write(fates_log(), *) 'Trying to define dimension size to a dim-type structure'
       write(fates_log(), *) 'but the dimension index does not exist'
       write(fates_log(), *) 'type: ',dk_name,' ndims: ',this%dim_kinds(ityp)%ndims,' input dim:',idim
       stop
       !end_run
    end if

    if (idim == 1) then
       this%dim_kinds(ityp)%dim1_index = dim_index
    else if (idim == 2) then
       this%dim_kinds(ityp)%dim2_index = dim_index
    end if

    ! With the map, we can set the dimension size
    this%dim_kinds(ityp)%dimsize(idim) = this%dim_bounds(dim_index)%upper_bound - &
         this%dim_bounds(dim_index)%lower_bound + 1

 end subroutine set_dim_indices
  
 ! =======================================================================
 subroutine set_patch_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%patch_index_ = index
 end subroutine set_patch_index

 integer function patch_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   patch_index = this%patch_index_
 end function patch_index

 ! =======================================================================
 subroutine set_column_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%column_index_ = index
 end subroutine set_column_index

 integer function column_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   column_index = this%column_index_
 end function column_index

 ! =======================================================================
 subroutine set_levgrnd_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levgrnd_index_ = index
 end subroutine set_levgrnd_index

 integer function levgrnd_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levgrnd_index = this%levgrnd_index_
 end function levgrnd_index

 ! =======================================================================
 subroutine set_levscpf_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levscpf_index_ = index
 end subroutine set_levscpf_index

 integer function levscpf_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levscpf_index = this%levscpf_index_
 end function levscpf_index

 ! =======================================================================
 subroutine set_levscls_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levscls_index_ = index
 end subroutine set_levscls_index

 integer function levscls_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levscls_index = this%levscls_index_
 end function levscls_index

!=========================================================================
 subroutine set_levcacls_index(this, index)
  implicit none
  class(fates_history_interface_type), intent(inout) :: this
  integer, intent(in) :: index
  this%levcacls_index_ = index
end subroutine set_levcacls_index

integer function levcacls_index(this)
  implicit none
  class(fates_history_interface_type), intent(in) :: this
  levcacls_index = this%levcacls_index_
end function levcacls_index

!=========================================================================
 subroutine set_levcapf_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levcapf_index_ = index
 end subroutine set_levcapf_index

integer function levcapf_index(this)
  implicit none
  class(fates_history_interface_type), intent(in) :: this
  levcapf_index = this%levcapf_index_
end function levcapf_index

 ! =======================================================================
 subroutine set_levpft_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levpft_index_ = index
 end subroutine set_levpft_index

 integer function levpft_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levpft_index = this%levpft_index_
 end function levpft_index

 ! =======================================================================
 subroutine set_levage_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levage_index_ = index
 end subroutine set_levage_index

 integer function levage_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levage_index = this%levage_index_
 end function levage_index

 ! =======================================================================
 subroutine set_levfuel_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levfuel_index_ = index
 end subroutine set_levfuel_index

 integer function levfuel_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levfuel_index = this%levfuel_index_
 end function levfuel_index

 ! =======================================================================
 subroutine set_levcwdsc_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levcwdsc_index_ = index
 end subroutine set_levcwdsc_index

 integer function levcwdsc_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levcwdsc_index = this%levcwdsc_index_
 end function levcwdsc_index

 ! =======================================================================
 subroutine set_levcan_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levcan_index_ = index
 end subroutine set_levcan_index

 integer function levcan_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levcan_index = this%levcan_index_
 end function levcan_index

 ! =======================================================================
 subroutine set_levcnlf_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levcnlf_index_ = index
 end subroutine set_levcnlf_index

 integer function levcnlf_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levcnlf_index = this%levcnlf_index_
 end function levcnlf_index

 ! =======================================================================
 subroutine set_levcnlfpft_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levcnlfpft_index_ = index
 end subroutine set_levcnlfpft_index

 integer function levcnlfpft_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levcnlfpft_index = this%levcnlfpft_index_
 end function levcnlfpft_index

 ! ======================================================================================
 subroutine set_levscag_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levscag_index_ = index
 end subroutine set_levscag_index

 integer function levscag_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levscag_index = this%levscag_index_
 end function levscag_index

 ! ======================================================================================
 subroutine set_levscagpft_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levscagpft_index_ = index
 end subroutine set_levscagpft_index

 integer function levscagpft_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levscagpft_index = this%levscagpft_index_
 end function levscagpft_index

 ! ======================================================================================
 subroutine set_levagepft_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levagepft_index_ = index
 end subroutine set_levagepft_index

 integer function levagepft_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levagepft_index = this%levagepft_index_
 end function levagepft_index

 ! ======================================================================================
 subroutine set_levheight_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levheight_index_ = index
 end subroutine set_levheight_index

 integer function levheight_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levheight_index = this%levheight_index_
 end function levheight_index

 ! ======================================================================================

 subroutine set_levelem_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levelem_index_ = index
 end subroutine set_levelem_index

 integer function levelem_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levelem_index = this%levelem_index_
  end function levelem_index

 ! ======================================================================================
       
 subroutine set_levelpft_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levelpft_index_ = index
 end subroutine set_levelpft_index

 integer function levelpft_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levelpft_index = this%levelpft_index_
 end function levelpft_index

 ! ======================================================================================

 subroutine set_levelcwd_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levelcwd_index_ = index
 end subroutine set_levelcwd_index

 integer function levelcwd_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levelcwd_index = this%levelcwd_index_
  end function levelcwd_index

 ! ======================================================================================

 subroutine set_levelage_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levelage_index_ = index
 end subroutine set_levelage_index

 integer function levelage_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levelage_index = this%levelage_index_
 end function levelage_index

 ! ======================================================================================

 subroutine flush_hvars(this,nc,upfreq_in)
 
   class(fates_history_interface_type) :: this
   integer,intent(in)                  :: nc
   integer,intent(in)                  :: upfreq_in
   integer                             :: ivar
   integer                             :: lb1,ub1,lb2,ub2

   do ivar=1,ubound(this%hvars,1)
      if (this%hvars(ivar)%upfreq == upfreq_in) then ! Only flush variables with update on dynamics step
         
         if(this%dim_kinds(this%hvars(ivar)%dim_kinds_index)%ndims == 1) then
            this%hvars(ivar)%r81d(this%iovar_map(nc)%site_index(:)) = this%hvars(ivar)%flushval
         else
            this%hvars(ivar)%r82d(this%iovar_map(nc)%site_index(:),:) = this%hvars(ivar)%flushval
         end if
         
      end if
   end do
   
end subroutine flush_hvars

  
  ! =====================================================================================
   
  subroutine set_history_var(this, vname, units, long, use_default, avgflag, vtype, &
       hlms, flushval, upfreq, ivar, initialize, index)

    use FatesUtilsMod, only     : check_hlm_list
    use FatesInterfaceTypesMod, only : hlm_name

    implicit none
    
    ! arguments
    class(fates_history_interface_type), intent(inout) :: this
    character(len=*), intent(in)  :: vname
    character(len=*), intent(in)  :: units
    character(len=*), intent(in)  :: long
    character(len=*), intent(in)  :: use_default
    character(len=*), intent(in)  :: avgflag
    character(len=*), intent(in)  :: vtype
    character(len=*), intent(in)  :: hlms
    real(r8), intent(in)          :: flushval ! IF THE TYPE IS AN INT WE WILL round with NINT
    integer, intent(in)           :: upfreq
    logical, intent(in) :: initialize
    integer, intent(inout)       :: ivar
    integer, intent(inout)       :: index  ! This is the index for the variable of
                                           ! interest that is associated with an
                                           ! explict name (for fast reference during update)
                                           ! A zero is passed back when the variable is
                                           ! not used
    

    ! locals
    integer :: ub1, lb1, ub2, lb2    ! Bounds for allocating the var
    integer :: ityp

    logical :: write_var

    write_var = check_hlm_list(trim(hlms), trim(hlm_name))
    if( write_var ) then
       ivar  = ivar+1
       index = ivar    
       
       if (initialize) then
          call this%hvars(ivar)%Init(vname, units, long, use_default, &
               vtype, avgflag, flushval, upfreq, &
               fates_history_num_dim_kinds, this%dim_kinds, this%dim_bounds)
       end if
    else
       index = 0
    end if
    
    return
  end subroutine set_history_var
  
  ! ====================================================================================
  
  subroutine init_dim_kinds_maps(this)
    
    ! ----------------------------------------------------------------------------------
    ! This subroutine simply initializes the structures that define the different
    ! array and type formats for different IO variables
    !
    ! PA_R8   : 1D patch scale 8-byte reals
    ! SI_R8   : 1D site scale 8-byte reals
    !
    ! The allocation on the structures is not dynamic and should only add up to the
    ! number of entries listed here.
    !
    ! ----------------------------------------------------------------------------------
    use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
    use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_coage_r8, site_coage_pft_r8
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
    use FatesIOVariableKindMod, only : site_height_r8
    use FatesIOVariableKindMod, only : site_elem_r8, site_elpft_r8
    use FatesIOVariableKindMod, only : site_elcwd_r8, site_elage_r8

    implicit none
    
    ! Arguments
    class(fates_history_interface_type), intent(inout) :: this
       

    integer :: index

    ! 1d Patch
    index = 1
    call this%dim_kinds(index)%Init(patch_r8, 1)

    ! 1d Site
    index = index + 1
    call this%dim_kinds(index)%Init(site_r8, 1)

    ! patch x ground
    index = index + 1
    call this%dim_kinds(index)%Init(patch_ground_r8, 2)

    ! patch x size-class/pft
    index = index + 1
    call this%dim_kinds(index)%Init(patch_size_pft_r8, 2)

    ! site x ground
    index = index + 1
    call this%dim_kinds(index)%Init(site_ground_r8, 2)

    ! site x size-class/pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_size_pft_r8, 2)

    ! site x size-class
    index = index + 1
    call this%dim_kinds(index)%Init(site_size_r8, 2)

    ! site x cohort age-class/pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_coage_pft_r8, 2)

    ! site x cohort age-class
    index = index + 1
    call this%dim_kinds(index)%Init(site_coage_r8, 2)

    ! site x pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_pft_r8, 2)

    ! site x patch-age class
    index = index + 1
    call this%dim_kinds(index)%Init(site_age_r8, 2)

    ! site x fuel size class
    index = index + 1
    call this%dim_kinds(index)%Init(site_fuel_r8, 2)

    ! site x cwd size class
    index = index + 1
    call this%dim_kinds(index)%Init(site_cwdsc_r8, 2)

    ! site x can class
    index = index + 1
    call this%dim_kinds(index)%Init(site_can_r8, 2)

    ! site x cnlf class
    index = index + 1
    call this%dim_kinds(index)%Init(site_cnlf_r8, 2)

    ! site x cnlfpft class
    index = index + 1
    call this%dim_kinds(index)%Init(site_cnlfpft_r8, 2)

    ! site x size-class x age class
    index = index + 1
    call this%dim_kinds(index)%Init(site_scag_r8, 2)

    ! site x size-class x age class x pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_scagpft_r8, 2)

    ! site x age class x pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_agepft_r8, 2)

    ! site x height
    index = index + 1
    call this%dim_kinds(index)%Init(site_height_r8, 2)

    ! site x elemenet
    index = index + 1
    call this%dim_kinds(index)%Init(site_elem_r8, 2)

    ! site x element x pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_elpft_r8, 2)
    
    ! site x element x cwd
    index = index + 1
    call this%dim_kinds(index)%Init(site_elcwd_r8, 2)

    ! site x element x age
    index = index + 1
    call this%dim_kinds(index)%Init(site_elage_r8, 2)


    ! FIXME(bja, 2016-10) assert(index == fates_history_num_dim_kinds)
  end subroutine init_dim_kinds_maps

 ! =======================================================================



  ! ====================================================================================
  
  subroutine update_history_dyn(this,nc,nsites,sites)
    
    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after Ecosystem Dynamics have been processed.
    ! ---------------------------------------------------------------------------------
    

    use EDtypesMod          , only : nfsc
    use FatesLitterMod      , only : ncwd
    use EDtypesMod          , only : ican_upper
    use EDtypesMod          , only : ican_ustory
    use FatesSizeAgeTypeIndicesMod, only : get_sizeage_class_index
    use FatesSizeAgeTypeIndicesMod, only : get_sizeagepft_class_index
    use FatesSizeAgeTypeIndicesMod, only : get_agepft_class_index
    use FatesSizeAgeTypeIndicesMod, only : get_age_class_index
    use FatesSizeAgeTypeIndicesMod, only : get_height_index
    use FatesSizeAgeTypeIndicesMod, only : sizetype_class_index
    use FatesSizeAgeTypeIndicesMod, only : coagetype_class_index
    use EDTypesMod        , only : nlevleaf
    use EDParamsMod,           only : ED_val_history_height_bin_edges

    ! Arguments
    class(fates_history_interface_type)             :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    
    ! Locals
    type(litter_type), pointer         :: litt_c   ! Pointer to the carbon12 litter pool
    type(litter_type), pointer         :: litt     ! Generic pointer to any litter pool
    type(site_fluxdiags_type), pointer :: flux_diags
    type(site_fluxdiags_type), pointer :: flux_diags_c
    type(site_massbal_type), pointer :: site_mass

    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ipa, ipa2 ! The local "I"ndex of "PA"tches 
    integer  :: io_soipa 
    integer  :: lb1,ub1,lb2,ub2  ! IO array bounds for the calling thread
    integer  :: ivar             ! index of IO variable object vector
    integer  :: ft               ! functional type index
    integer  :: cwd
    integer  :: elcwd, elpft            ! combined index of element and pft or cwd
    integer  :: i_scpf,i_pft,i_sz     ! iterators for scpf, pft, and size-class dims
    integer  :: i_cacls, i_capf      ! iterators for cohort age and cohort age x pft
    integer  :: i_cwd,i_fuel            ! iterators for cwd and fuel dims
    integer  :: iscag        ! size-class x age index
    integer  :: iscagpft     ! size-class x age x pft index
    integer  :: iagepft      ! age x pft index
    integer  :: ican, ileaf, cnlf_indx  ! iterators for leaf and canopy level
    integer  :: height_bin_max, height_bin_min   ! which height bin a given cohort's canopy is in
    integer  :: i_heightbin  ! iterator for height bins
    integer  :: el           ! Loop index for elements
    integer  :: model_day_int ! integer model day from reference 
    integer  :: ageclass_since_anthrodist  ! what is the equivalent age class for
                                           ! time-since-anthropogenic-disturbance of secondary forest

    
    real(r8) :: n_perm2     ! individuals per m2 for the whole column
    real(r8) :: dbh         ! diameter ("at breast height")
    real(r8) :: coage       ! cohort age 
    real(r8) :: npp_partition_error ! a check that the NPP partitions sum to carbon allocation
    real(r8) :: frac_canopy_in_bin  ! fraction of a leaf's canopy that is within a given height bin
    real(r8) :: binbottom,bintop    ! edges of height bins
    
    real(r8) :: gpp_cached ! variable used to cache gpp value in previous time step; for C13 discrimination

    ! The following are all carbon states, turnover and net allocation flux variables
    ! the organs of relevance should be self explanatory
    real(r8) :: sapw_m    ! Sapwood mass (elemental, c,n or p) [kg/plant]
    real(r8) :: struct_m  ! Structural mass ""
    real(r8) :: leaf_m    ! Leaf mass ""
    real(r8) :: fnrt_m    ! Fineroot mass ""
    real(r8) :: store_m   ! Storage mass ""
    real(r8) :: alive_m   ! Alive biomass (sap+leaf+fineroot+repro+storage) ""
    real(r8) :: total_m   ! Total vegetation mass
    real(r8) :: repro_m   ! Total reproductive mass (on plant) ""
    real(r8) :: sapw_m_turnover
    real(r8) :: store_m_turnover
    real(r8) :: leaf_m_turnover
    real(r8) :: fnrt_m_turnover
    real(r8) :: struct_m_turnover
    real(r8) :: sapw_m_net_alloc
    real(r8) :: store_m_net_alloc
    real(r8) :: leaf_m_net_alloc
    real(r8) :: fnrt_m_net_alloc
    real(r8) :: struct_m_net_alloc
    real(r8) :: repro_m_net_alloc
    real(r8) :: area_frac

    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort

    real(r8), parameter :: tiny = 1.e-5_r8      ! some small number
    real(r8), parameter :: reallytalltrees = 1000.   ! some large number (m)
    
    integer :: tmp

    associate( hio_npatches         => this%hvars(ih_npatches)%r81d, &
               hio_ncohorts         => this%hvars(ih_ncohorts)%r81d, &
               hio_trimming         => this%hvars(ih_trimming)%r81d, &
               hio_area_plant       => this%hvars(ih_area_plant)%r81d, &
               hio_area_trees  => this%hvars(ih_area_trees)%r81d, & 
               hio_canopy_spread    => this%hvars(ih_canopy_spread)%r81d, &
               hio_biomass_pft      => this%hvars(ih_biomass_pft)%r82d, &
               hio_leafbiomass_pft  => this%hvars(ih_leafbiomass_pft)%r82d, &
               hio_storebiomass_pft => this%hvars(ih_storebiomass_pft)%r82d, &
               hio_nindivs_pft      => this%hvars(ih_nindivs_pft)%r82d, &
               hio_recruitment_pft  => this%hvars(ih_recruitment_pft)%r82d, &
               hio_mortality_pft    => this%hvars(ih_mortality_pft)%r82d, &
               hio_crownarea_pft    => this%hvars(ih_crownarea_pft)%r82d, &
               hio_canopycrownarea_pft  => this%hvars(ih_canopycrownarea_pft)%r82d, &
               hio_nesterov_fire_danger => this%hvars(ih_nesterov_fire_danger)%r81d, &
               hio_fire_nignitions => this%hvars(ih_fire_nignitions)%r81d, &
               hio_fire_fdi => this%hvars(ih_fire_fdi)%r81d, &
               hio_spitfire_ros     => this%hvars(ih_spitfire_ros)%r81d, &
               hio_fire_ros_area_product=> this%hvars(ih_fire_ros_area_product)%r81d, &
               hio_tfc_ros          => this%hvars(ih_tfc_ros)%r81d, &
               hio_tfc_ros_area_product => this%hvars(ih_tfc_ros_area_product)%r81d, &
               hio_effect_wspeed    => this%hvars(ih_effect_wspeed)%r81d, &
               hio_fire_intensity   => this%hvars(ih_fire_intensity)%r81d, &
               hio_fire_intensity_area_product => this%hvars(ih_fire_intensity_area_product)%r81d, &
               hio_fire_area        => this%hvars(ih_fire_area)%r81d, &
               hio_fire_fuel_bulkd  => this%hvars(ih_fire_fuel_bulkd)%r81d, &
               hio_fire_fuel_eff_moist => this%hvars(ih_fire_fuel_eff_moist)%r81d, &
               hio_fire_fuel_sav    => this%hvars(ih_fire_fuel_sav)%r81d, &
               hio_fire_fuel_mef    => this%hvars(ih_fire_fuel_mef)%r81d, &
               hio_sum_fuel         => this%hvars(ih_sum_fuel)%r81d,  &
               hio_fragmentation_scaler  => this%hvars(ih_fragmentation_scaler)%r81d,  &
               hio_litter_in        => this%hvars(ih_litter_in)%r81d, &
               hio_litter_out       => this%hvars(ih_litter_out)%r81d, &
               hio_seed_bank        => this%hvars(ih_seed_bank)%r81d, &
               hio_seeds_in         => this%hvars(ih_seeds_in)%r81d, &
               hio_litter_in_elem      => this%hvars(ih_litter_in_elem)%r82d, &
               hio_litter_out_elem     => this%hvars(ih_litter_out_elem)%r82d, &
               hio_seed_bank_elem      => this%hvars(ih_seed_bank_elem)%r82d, &
               hio_seeds_in_local_elem => this%hvars(ih_seeds_in_local_elem)%r82d, &
               hio_seed_in_extern_elem => this%hvars(ih_seeds_in_extern_elem)%r82d, & 
               hio_seed_decay_elem     => this%hvars(ih_seed_decay_elem)%r82d, &
               hio_seed_germ_elem      => this%hvars(ih_seed_germ_elem)%r82d, &
               hio_agb              => this%hvars(ih_agb)%r81d, &
               hio_canopy_biomass   => this%hvars(ih_canopy_biomass)%r81d, &
               hio_understory_biomass   => this%hvars(ih_understory_biomass)%r81d, &
               hio_primaryland_fusion_error    => this%hvars(ih_primaryland_fusion_error)%r81d, &
               hio_disturbance_rate_p2p       => this%hvars(ih_disturbance_rate_p2p)%r81d, &
               hio_disturbance_rate_p2s       => this%hvars(ih_disturbance_rate_p2s)%r81d, &
               hio_disturbance_rate_s2s       => this%hvars(ih_disturbance_rate_s2s)%r81d, &
               hio_fire_disturbance_rate      => this%hvars(ih_fire_disturbance_rate)%r81d, &
               hio_logging_disturbance_rate   => this%hvars(ih_logging_disturbance_rate)%r81d, &
               hio_fall_disturbance_rate      => this%hvars(ih_fall_disturbance_rate)%r81d, &
               hio_potential_disturbance_rate => this%hvars(ih_potential_disturbance_rate)%r81d, &
               hio_harvest_carbonflux => this%hvars(ih_harvest_carbonflux)%r81d, &
               hio_gpp_scpf         => this%hvars(ih_gpp_scpf)%r82d, &
               hio_npp_totl_scpf    => this%hvars(ih_npp_totl_scpf)%r82d, &
               hio_npp_leaf_scpf    => this%hvars(ih_npp_leaf_scpf)%r82d, &
               hio_npp_seed_scpf    => this%hvars(ih_npp_seed_scpf)%r82d, &
               hio_npp_fnrt_scpf    => this%hvars(ih_npp_fnrt_scpf)%r82d, &
               hio_npp_bgsw_scpf    => this%hvars(ih_npp_bgsw_scpf)%r82d, &
               hio_npp_bgdw_scpf    => this%hvars(ih_npp_bgdw_scpf)%r82d, &
               hio_npp_agsw_scpf    => this%hvars(ih_npp_agsw_scpf)%r82d, &
               hio_npp_agdw_scpf    => this%hvars(ih_npp_agdw_scpf)%r82d, &
               hio_npp_stor_scpf    => this%hvars(ih_npp_stor_scpf)%r82d, &
               hio_npp_leaf         => this%hvars(ih_npp_leaf)%r81d, &
               hio_npp_seed         => this%hvars(ih_npp_seed)%r81d, &
               hio_npp_stem         => this%hvars(ih_npp_stem)%r81d, &
               hio_npp_froot        => this%hvars(ih_npp_froot)%r81d, &
               hio_npp_croot        => this%hvars(ih_npp_croot)%r81d, &
               hio_npp_stor         => this%hvars(ih_npp_stor)%r81d, &
               hio_bstor_canopy_scpf      => this%hvars(ih_bstor_canopy_scpf)%r82d, &
               hio_bstor_understory_scpf  => this%hvars(ih_bstor_understory_scpf)%r82d, &
               hio_bleaf_canopy_scpf      => this%hvars(ih_bleaf_canopy_scpf)%r82d, &
               hio_bleaf_understory_scpf  => this%hvars(ih_bleaf_understory_scpf)%r82d, &
               hio_mortality_canopy_scpf         => this%hvars(ih_mortality_canopy_scpf)%r82d, &
               hio_mortality_understory_scpf     => this%hvars(ih_mortality_understory_scpf)%r82d, &
               hio_nplant_canopy_scpf     => this%hvars(ih_nplant_canopy_scpf)%r82d, &
               hio_nplant_understory_scpf => this%hvars(ih_nplant_understory_scpf)%r82d, &
               hio_ddbh_canopy_scpf       => this%hvars(ih_ddbh_canopy_scpf)%r82d, &
               hio_ddbh_understory_scpf   => this%hvars(ih_ddbh_understory_scpf)%r82d, &
               hio_ddbh_canopy_sz       => this%hvars(ih_ddbh_canopy_sz)%r82d, &
               hio_ddbh_understory_sz   => this%hvars(ih_ddbh_understory_sz)%r82d, &
               hio_gpp_canopy_scpf        => this%hvars(ih_gpp_canopy_scpf)%r82d, &
               hio_gpp_understory_scpf    => this%hvars(ih_gpp_understory_scpf)%r82d, &
               hio_ar_canopy_scpf         => this%hvars(ih_ar_canopy_scpf)%r82d, &
               hio_ar_understory_scpf     => this%hvars(ih_ar_understory_scpf)%r82d, &
               hio_ddbh_scpf        => this%hvars(ih_ddbh_scpf)%r82d, &
               hio_growthflux_scpf        => this%hvars(ih_growthflux_scpf)%r82d, &
               hio_growthflux_fusion_scpf        => this%hvars(ih_growthflux_fusion_scpf)%r82d, &
               hio_ba_scpf          => this%hvars(ih_ba_scpf)%r82d, &
               hio_agb_scpf         => this%hvars(ih_agb_scpf)%r82d, &
               hio_nplant_scpf      => this%hvars(ih_nplant_scpf)%r82d, &
               hio_nplant_capf      => this%hvars(ih_nplant_capf)%r82d, &
               
               hio_m1_scpf          => this%hvars(ih_m1_scpf)%r82d, &
               hio_m2_scpf          => this%hvars(ih_m2_scpf)%r82d, &
               hio_m3_scpf          => this%hvars(ih_m3_scpf)%r82d, &
               hio_m4_scpf          => this%hvars(ih_m4_scpf)%r82d, &
               hio_m5_scpf          => this%hvars(ih_m5_scpf)%r82d, &
               hio_m6_scpf          => this%hvars(ih_m6_scpf)%r82d, &
               hio_m7_scpf          => this%hvars(ih_m7_scpf)%r82d, &                  
               hio_m8_scpf          => this%hvars(ih_m8_scpf)%r82d, &
               hio_m9_scpf          => this%hvars(ih_m9_scpf)%r82d, &
               hio_m10_scpf         => this%hvars(ih_m10_scpf)%r82d, &
               hio_m10_capf         => this%hvars(ih_m10_capf)%r82d, &
      
               hio_crownfiremort_scpf     => this%hvars(ih_crownfiremort_scpf)%r82d, &
               hio_cambialfiremort_scpf   => this%hvars(ih_cambialfiremort_scpf)%r82d, &

               hio_fire_c_to_atm  => this%hvars(ih_fire_c_to_atm)%r81d, &
               hio_burn_flux_elem    => this%hvars(ih_burn_flux_elem)%r82d, &

               hio_m1_sz          => this%hvars(ih_m1_sz)%r82d, &
               hio_m2_sz          => this%hvars(ih_m2_sz)%r82d, &
               hio_m3_sz          => this%hvars(ih_m3_sz)%r82d, &
               hio_m4_sz          => this%hvars(ih_m4_sz)%r82d, &
               hio_m5_sz          => this%hvars(ih_m5_sz)%r82d, &
               hio_m6_sz          => this%hvars(ih_m6_sz)%r82d, &
               hio_m7_sz          => this%hvars(ih_m7_sz)%r82d, &
               hio_m8_sz          => this%hvars(ih_m8_sz)%r82d, &
               hio_m9_sz          => this%hvars(ih_m9_sz)%r82d, &
               hio_m10_sz         => this%hvars(ih_m10_sz)%r82d, &
               hio_m10_cacls        => this%hvars(ih_m10_cacls)%r82d, &
              
               hio_c13disc_scpf     => this%hvars(ih_c13disc_scpf)%r82d, &                    

               hio_cwd_elcwd           => this%hvars(ih_cwd_elcwd)%r82d, &
               hio_cwd_ag_elem         => this%hvars(ih_cwd_ag_elem)%r82d, &
               hio_cwd_bg_elem         => this%hvars(ih_cwd_bg_elem)%r82d, &
               hio_fines_ag_elem       => this%hvars(ih_fines_bg_elem)%r82d, &
               hio_fines_bg_elem       => this%hvars(ih_fines_ag_elem)%r82d, &
               hio_ba_sz          => this%hvars(ih_ba_sz)%r82d, &
               hio_agb_sz          => this%hvars(ih_agb_sz)%r82d, &
               hio_biomass_sz          => this%hvars(ih_biomass_sz)%r82d, &
               hio_nplant_sz         => this%hvars(ih_nplant_sz)%r82d, &
               hio_nplant_cacls        => this%hvars(ih_nplant_cacls)%r82d, &
               hio_nplant_canopy_sz         => this%hvars(ih_nplant_canopy_sz)%r82d, &
               hio_nplant_understory_sz     => this%hvars(ih_nplant_understory_sz)%r82d, &
               hio_lai_canopy_sz         => this%hvars(ih_lai_canopy_sz)%r82d, &
               hio_lai_understory_sz     => this%hvars(ih_lai_understory_sz)%r82d, &
               hio_sai_canopy_sz         => this%hvars(ih_sai_canopy_sz)%r82d, &
               hio_sai_understory_sz     => this%hvars(ih_sai_understory_sz)%r82d, &
               hio_mortality_canopy_sz      => this%hvars(ih_mortality_canopy_sz)%r82d, &
               hio_mortality_understory_sz  => this%hvars(ih_mortality_understory_sz)%r82d, &
               hio_demotion_rate_sz         => this%hvars(ih_demotion_rate_sz)%r82d, &
               hio_demotion_carbonflux        => this%hvars(ih_demotion_carbonflux)%r81d, &
               hio_promotion_rate_sz        => this%hvars(ih_promotion_rate_sz)%r82d, &
               hio_trimming_canopy_sz         => this%hvars(ih_trimming_canopy_sz)%r82d, &
               hio_trimming_understory_sz     => this%hvars(ih_trimming_understory_sz)%r82d, &
               hio_crown_area_canopy_sz         => this%hvars(ih_crown_area_canopy_sz)%r82d, &
               hio_crown_area_understory_sz     => this%hvars(ih_crown_area_understory_sz)%r82d, &
               hio_promotion_carbonflux       => this%hvars(ih_promotion_carbonflux)%r81d, &
               hio_canopy_mortality_carbonflux     => this%hvars(ih_canopy_mortality_carbonflux)%r81d, &
               hio_understory_mortality_carbonflux => this%hvars(ih_understory_mortality_carbonflux)%r81d, &
               hio_leaf_md_canopy_sz           => this%hvars(ih_leaf_md_canopy_sz)%r82d, &
               hio_root_md_canopy_sz           => this%hvars(ih_root_md_canopy_sz)%r82d, &
               hio_carbon_balance_canopy_sz    => this%hvars(ih_carbon_balance_canopy_sz)%r82d, &
               hio_bsw_md_canopy_sz            => this%hvars(ih_bsw_md_canopy_sz)%r82d, &
               hio_bdead_md_canopy_sz          => this%hvars(ih_bdead_md_canopy_sz)%r82d, &
               hio_bstore_md_canopy_sz         => this%hvars(ih_bstore_md_canopy_sz)%r82d, &
               hio_seed_prod_canopy_sz         => this%hvars(ih_seed_prod_canopy_sz)%r82d, &
               hio_npp_leaf_canopy_sz          => this%hvars(ih_npp_leaf_canopy_sz)%r82d, &
               hio_npp_fnrt_canopy_sz         => this%hvars(ih_npp_fnrt_canopy_sz)%r82d, &
               hio_npp_sapw_canopy_sz           => this%hvars(ih_npp_sapw_canopy_sz)%r82d, &
               hio_npp_dead_canopy_sz         => this%hvars(ih_npp_dead_canopy_sz)%r82d, &
               hio_npp_seed_canopy_sz         => this%hvars(ih_npp_seed_canopy_sz)%r82d, &
               hio_npp_stor_canopy_sz         => this%hvars(ih_npp_stor_canopy_sz)%r82d, &
               hio_leaf_md_understory_sz       => this%hvars(ih_leaf_md_understory_sz)%r82d, &
               hio_root_md_understory_sz       => this%hvars(ih_root_md_understory_sz)%r82d, &
               hio_carbon_balance_understory_sz=> this%hvars(ih_carbon_balance_understory_sz)%r82d, &
               hio_bstore_md_understory_sz     => this%hvars(ih_bstore_md_understory_sz)%r82d, &
               hio_bsw_md_understory_sz        => this%hvars(ih_bsw_md_understory_sz)%r82d, &
               hio_bdead_md_understory_sz      => this%hvars(ih_bdead_md_understory_sz)%r82d, &
               hio_seed_prod_understory_sz     => this%hvars(ih_seed_prod_understory_sz)%r82d, &
               hio_npp_leaf_understory_sz      => this%hvars(ih_npp_leaf_understory_sz)%r82d, &
               hio_npp_fnrt_understory_sz     => this%hvars(ih_npp_fnrt_understory_sz)%r82d, &
               hio_npp_sapw_understory_sz       => this%hvars(ih_npp_sapw_understory_sz)%r82d, &
               hio_npp_dead_understory_sz     => this%hvars(ih_npp_dead_understory_sz)%r82d, &
               hio_npp_seed_understory_sz     => this%hvars(ih_npp_seed_understory_sz)%r82d, &
               hio_npp_stor_understory_sz     => this%hvars(ih_npp_stor_understory_sz)%r82d, &
               hio_nplant_scagpft                => this%hvars(ih_nplant_scagpft)%r82d, &
               hio_npp_agepft                    => this%hvars(ih_npp_agepft)%r82d, &
               hio_biomass_agepft                => this%hvars(ih_biomass_agepft)%r82d, &
               hio_scorch_height_agepft          => this%hvars(ih_scorch_height_agepft)%r82d, &
               hio_yesterdaycanopylevel_canopy_sz     => this%hvars(ih_yesterdaycanopylevel_canopy_sz)%r82d, &
               hio_yesterdaycanopylevel_understory_sz => this%hvars(ih_yesterdaycanopylevel_understory_sz)%r82d, &
               hio_area_age         => this%hvars(ih_area_age)%r82d, &
               hio_lai_age          => this%hvars(ih_lai_age)%r82d, &
               hio_canopy_area_age  => this%hvars(ih_canopy_area_age)%r82d, &
               hio_ncl_age          => this%hvars(ih_ncl_age)%r82d, &
               hio_npatches_age     => this%hvars(ih_npatches_age)%r82d, &
               hio_zstar_age        => this%hvars(ih_zstar_age)%r82d, &
               hio_biomass_age        => this%hvars(ih_biomass_age)%r82d, &
               hio_fraction_secondary_forest   => this%hvars(ih_fraction_secondary_forest)%r81d, &
               hio_biomass_secondary_forest    => this%hvars(ih_biomass_secondary_forest)%r81d, &
               hio_woodproduct                 => this%hvars(ih_woodproduct)%r81d, &
               hio_agesince_anthrodist_age     => this%hvars(ih_agesince_anthrodist_age)%r82d, &
               hio_secondaryforest_area_age    => this%hvars(ih_secondaryforest_area_age)%r82d, &
               hio_area_burnt_age              => this%hvars(ih_area_burnt_age)%r82d, &
               ! hio_fire_rate_of_spread_front_age  => this%hvars(ih_fire_rate_of_spread_front_age)%r82d, &
               hio_fire_intensity_age          => this%hvars(ih_fire_intensity_age)%r82d, &
               hio_fire_sum_fuel_age           => this%hvars(ih_fire_sum_fuel_age)%r82d, &
               hio_burnt_frac_litter_fuel      => this%hvars(ih_burnt_frac_litter_fuel)%r82d, &
               hio_fuel_amount_fuel            => this%hvars(ih_fuel_amount_fuel)%r82d, &
               hio_canopy_height_dist_height   => this%hvars(ih_canopy_height_dist_height)%r82d, &
               hio_leaf_height_dist_height     => this%hvars(ih_leaf_height_dist_height)%r82d, &
               hio_litter_moisture_fuel        => this%hvars(ih_litter_moisture_fuel)%r82d, &
               hio_cwd_ag_cwdsc                  => this%hvars(ih_cwd_ag_cwdsc)%r82d, &
               hio_cwd_bg_cwdsc                  => this%hvars(ih_cwd_bg_cwdsc)%r82d, &
               hio_cwd_ag_in_cwdsc               => this%hvars(ih_cwd_ag_in_cwdsc)%r82d, &
               hio_cwd_bg_in_cwdsc               => this%hvars(ih_cwd_bg_in_cwdsc)%r82d, &
               hio_cwd_ag_out_cwdsc              => this%hvars(ih_cwd_ag_out_cwdsc)%r82d, &
               hio_cwd_bg_out_cwdsc              => this%hvars(ih_cwd_bg_out_cwdsc)%r82d, &
               hio_crownarea_cnlf                => this%hvars(ih_crownarea_cnlf)%r82d, &
               hio_crownarea_can                 => this%hvars(ih_crownarea_can)%r82d, &
               hio_nplant_scag                   => this%hvars(ih_nplant_scag)%r82d, &
               hio_nplant_canopy_scag            => this%hvars(ih_nplant_canopy_scag)%r82d, &
               hio_nplant_understory_scag        => this%hvars(ih_nplant_understory_scag)%r82d, &
               hio_ddbh_canopy_scag              => this%hvars(ih_ddbh_canopy_scag)%r82d, &
               hio_ddbh_understory_scag          => this%hvars(ih_ddbh_understory_scag)%r82d, &
               hio_mortality_canopy_scag         => this%hvars(ih_mortality_canopy_scag)%r82d, &
               hio_mortality_understory_scag     => this%hvars(ih_mortality_understory_scag)%r82d, &
               hio_cstatus                  => this%hvars(ih_cstatus)%r81d, &
               hio_dstatus                  => this%hvars(ih_dstatus)%r81d, &
               hio_gdd                           => this%hvars(ih_gdd)%r81d, &
               hio_ncolddays                => this%hvars(ih_ncolddays)%r81d, &
               hio_nchilldays               => this%hvars(ih_nchilldays)%r81d, &
               hio_cleafoff                      => this%hvars(ih_cleafoff)%r81d, &
               hio_cleafon                       => this%hvars(ih_cleafon)%r81d, &
               hio_dleafoff                      => this%hvars(ih_dleafoff)%r81d, &
               hio_dleafon                       => this%hvars(ih_dleafoff)%r81d, &
               hio_meanliqvol                    => this%hvars(ih_meanliqvol)%r81d, &
               hio_cbal_err_fates                => this%hvars(ih_cbal_err_fates)%r81d, &
               hio_err_fates                     => this%hvars(ih_err_fates)%r82d )

               
      ! ---------------------------------------------------------------------------------
      ! Flush arrays to values defined by %flushval (see registry entry in
      ! subroutine define_history_vars()
      ! ---------------------------------------------------------------------------------
      call this%flush_hvars(nc,upfreq_in=1)


      ! If we don't have dynamics turned on, we just abort these diagnostics
      if (hlm_use_ed_st3.eq.itrue) return

      model_day_int = nint(hlm_model_day)

      ! ---------------------------------------------------------------------------------
      ! Loop through the FATES scale hierarchy and fill the history IO arrays
      ! ---------------------------------------------------------------------------------
      
      do s = 1,nsites
         
         io_si  = this%iovar_map(nc)%site_index(s)

         if(nsites /= ubound(this%iovar_map(nc)%site_index,dim=1) ) then
            print*,"NSITE DNE SIZE OF SITE_INDEX"
            stop
         end if
         

         ! Total carbon model error [kgC/day -> mgC/day]
         hio_cbal_err_fates(io_si) = &
               sites(s)%mass_balance(element_pos(carbon12_element))%err_fates * mg_per_kg

         ! Total carbon lost to atmosphere from burning (kgC/site/day -> gC/m2/s)
         hio_fire_c_to_atm(io_si) = &
              sites(s)%mass_balance(element_pos(carbon12_element))%burn_flux_to_atm * &
              g_per_kg * ha_per_m2 * days_per_sec

         ! Total model error [kg/day -> mg/day]  (all elements)
         do el = 1, num_elements

             hio_err_fates(io_si,el) = sites(s)%mass_balance(el)%err_fates * mg_per_kg

             ! Total element lost to atmosphere from burning (kg/site/day -> g/m2/s)
             hio_burn_flux_elem(io_si,el) = &
                  sites(s)%mass_balance(el)%burn_flux_to_atm * &
                  g_per_kg * ha_per_m2 * days_per_sec

         end do

         hio_canopy_spread(io_si)        = sites(s)%spread

         ! Update the site statuses (stati?)
         hio_cstatus(io_si)   = real(sites(s)%cstatus,r8)
         hio_dstatus(io_si)   = real(sites(s)%dstatus,r8)

         !count number of days for leaves off
         hio_nchilldays(io_si) = real(sites(s)%nchilldays,r8)
         hio_ncolddays(io_si)  = real(sites(s)%ncolddays,r8)

            
         hio_gdd(io_si)      = sites(s)%grow_deg_days
         hio_cleafoff(io_si) = real(model_day_int - sites(s)%cleafoffdate,r8)
         hio_cleafon(io_si)  = real(model_day_int - sites(s)%cleafondate,r8)
         hio_dleafoff(io_si) = real(model_day_int - sites(s)%dleafoffdate,r8)
         hio_dleafon(io_si)  = real(model_day_int - sites(s)%dleafondate,r8)

         if(model_day_int>numWaterMem)then
            hio_meanliqvol(io_si) = &
                 sum(sites(s)%water_memory(1:numWaterMem))/real(numWaterMem,r8)
         end if

         ! track total wood product accumulation at the site level
         hio_woodproduct(io_si)          = sites(s)%resources_management%trunk_product_site &
              * AREA_INV * g_per_kg
         
         ! site-level fire variables
         hio_nesterov_fire_danger(io_si) = sites(s)%acc_NI
         hio_fire_nignitions(io_si) = sites(s)%NF
         hio_fire_fdi(io_si) = sites(s)%FDI

         ! If hydraulics are turned on, track the error terms
         ! associated with dynamics

         if(hlm_use_planthydro.eq.itrue)then
            this%hvars(ih_h2oveg_dead)%r81d(io_si)         = sites(s)%si_hydr%h2oveg_dead
            this%hvars(ih_h2oveg_recruit)%r81d(io_si)      = sites(s)%si_hydr%h2oveg_recruit
            this%hvars(ih_h2oveg_growturn_err)%r81d(io_si) = sites(s)%si_hydr%h2oveg_growturn_err
            this%hvars(ih_h2oveg_pheno_err)%r81d(io_si)    = sites(s)%si_hydr%h2oveg_pheno_err
         end if

         ! error in primary lands from patch fusion
         hio_primaryland_fusion_error(io_si) = sites(s)%primary_land_patchfusion_error

         ! output site-level disturbance rates
         hio_disturbance_rate_p2p(io_si) = sum(sites(s)%disturbance_rates_primary_to_primary(1:N_DIST_TYPES))
         hio_disturbance_rate_p2s(io_si) = sum(sites(s)%disturbance_rates_primary_to_secondary(1:N_DIST_TYPES))
         hio_disturbance_rate_s2s(io_si) = sum(sites(s)%disturbance_rates_secondary_to_secondary(1:N_DIST_TYPES))

         hio_fire_disturbance_rate(io_si) = sites(s)%disturbance_rates_primary_to_primary(dtype_ifire) + &
              sites(s)%disturbance_rates_primary_to_secondary(dtype_ifire) + &
              sites(s)%disturbance_rates_secondary_to_secondary(dtype_ifire)

         hio_logging_disturbance_rate(io_si) = sites(s)%disturbance_rates_primary_to_primary(dtype_ilog) + &
              sites(s)%disturbance_rates_primary_to_secondary(dtype_ilog) + &
              sites(s)%disturbance_rates_secondary_to_secondary(dtype_ilog)

         hio_fall_disturbance_rate(io_si) = sites(s)%disturbance_rates_primary_to_primary(dtype_ifall) + &
              sites(s)%disturbance_rates_primary_to_secondary(dtype_ifall) + &
              sites(s)%disturbance_rates_secondary_to_secondary(dtype_ifall)

         hio_potential_disturbance_rate(io_si) = sum(sites(s)%potential_disturbance_rates(1:N_DIST_TYPES))

         hio_harvest_carbonflux(io_si) = sites(s)%harvest_carbon_flux

         ipa = 0
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            
            ! Increment the number of patches per site
            hio_npatches(io_si) = hio_npatches(io_si) + 1._r8

            cpatch%age_class  = get_age_class_index(cpatch%age)

            ! Increment the fractional area in each age class bin
            hio_area_age(io_si,cpatch%age_class) = hio_area_age(io_si,cpatch%age_class) &
                 + cpatch%area * AREA_INV

            ! Increment some patch-age-resolved diagnostics
            hio_lai_age(io_si,cpatch%age_class) = hio_lai_age(io_si,cpatch%age_class) &
                  + sum(cpatch%tlai_profile(:,:,:)) * cpatch%area

            hio_ncl_age(io_si,cpatch%age_class) = hio_ncl_age(io_si,cpatch%age_class) &
                  + cpatch%ncl_p * cpatch%area
            hio_npatches_age(io_si,cpatch%age_class) = hio_npatches_age(io_si,cpatch%age_class) + 1._r8
            if ( ED_val_comp_excln .lt. 0._r8 ) then ! only valid when "strict ppa" enabled
               hio_zstar_age(io_si,cpatch%age_class) = hio_zstar_age(io_si,cpatch%age_class) &
                    + cpatch%zstar * cpatch%area * AREA_INV
            endif

            ! some diagnostics on secondary forest area and its age distribution
            if ( cpatch%anthro_disturbance_label .eq. secondaryforest ) then
               hio_fraction_secondary_forest(io_si) = hio_fraction_secondary_forest(io_si) + &
                    cpatch%area * AREA_INV
               
               ageclass_since_anthrodist = get_age_class_index(cpatch%age_since_anthro_disturbance)
               
               hio_agesince_anthrodist_age(io_si,ageclass_since_anthrodist) = &
                    hio_agesince_anthrodist_age(io_si,ageclass_since_anthrodist)  &
                    + cpatch%area * AREA_INV

               hio_secondaryforest_area_age(io_si,cpatch%age_class) = &
                    hio_secondaryforest_area_age(io_si,cpatch%age_class)  &
                    + cpatch%area * AREA_INV
            endif
            
            !!! patch-age-resolved fire variables
            do i_pft = 1,numpft
               ! for scorch height, weight the value by patch area within any given age calss (in the event that there is
               ! more than one patch per age class.
               iagepft = cpatch%age_class + (i_pft-1) * nlevage
               hio_scorch_height_agepft(io_si,iagepft) = hio_scorch_height_agepft(io_si,iagepft) + &
                    cpatch%Scorch_ht(i_pft) * cpatch%area
            end do

            hio_area_burnt_age(io_si,cpatch%age_class) = hio_area_burnt_age(io_si,cpatch%age_class) + &
                 cpatch%frac_burnt * cpatch%area * AREA_INV

            ! hio_fire_rate_of_spread_front_age(io_si, cpatch%age_class) = hio_fire_rate_of_spread_age(io_si, cpatch%age_class) + &
            !      cpatch%ros_front * cpatch*frac_burnt * cpatch%area * AREA_INV

            hio_fire_intensity_age(io_si, cpatch%age_class) = hio_fire_intensity_age(io_si, cpatch%age_class) + &
                 cpatch%FI * cpatch%frac_burnt * cpatch%area * AREA_INV

            hio_fire_sum_fuel_age(io_si, cpatch%age_class) = hio_fire_sum_fuel_age(io_si, cpatch%age_class) + &
                 cpatch%sum_fuel * g_per_kg * cpatch%area * AREA_INV
             
            if(associated(cpatch%tallest))then
               hio_trimming(io_si) = hio_trimming(io_si) + cpatch%tallest%canopy_trim * cpatch%area * AREA_INV
            endif
            
            hio_area_plant(io_si) = hio_area_plant(io_si) + min(cpatch%total_canopy_area,cpatch%area) * AREA_INV

            hio_area_trees(io_si) = hio_area_trees(io_si) + min(cpatch%total_tree_area,cpatch%area) * AREA_INV
            
            ccohort => cpatch%shortest
            do while(associated(ccohort))
               
               ft = ccohort%pft

               call sizetype_class_index(ccohort%dbh, ccohort%pft, ccohort%size_class, ccohort%size_by_pft_class)
               call coagetype_class_index(ccohort%coage, ccohort%pft, &
                                          ccohort%coage_class, ccohort%coage_by_pft_class)
              
               ! Increment the number of cohorts per site
               hio_ncohorts(io_si) = hio_ncohorts(io_si) + 1._r8
               
               n_perm2   = ccohort%n * AREA_INV
                              
               hio_canopy_area_age(io_si,cpatch%age_class) = hio_canopy_area_age(io_si,cpatch%age_class) &
                    + ccohort%c_area * AREA_INV

               ! calculate leaf height distribution, assuming leaf area is evenly distributed thru crown depth
               height_bin_max = get_height_index(ccohort%hite)
               height_bin_min = get_height_index(ccohort%hite * (1._r8 - EDPftvarcon_inst%crown(ft)))
               do i_heightbin = height_bin_min, height_bin_max
                  binbottom = ED_val_history_height_bin_edges(i_heightbin)
                  if (i_heightbin .eq. nlevheight) then
                     bintop = reallytalltrees
                  else
                     bintop = ED_val_history_height_bin_edges(i_heightbin+1)
                  endif
                  ! what fraction of a cohort's crown is in this height bin?
                  frac_canopy_in_bin = (min(bintop,ccohort%hite) - &
                       max(binbottom,ccohort%hite * (1._r8 - EDPftvarcon_inst%crown(ft)))) / &
                       (ccohort%hite * EDPftvarcon_inst%crown(ft))
                  !
                  hio_leaf_height_dist_height(io_si,i_heightbin) = &
                       hio_leaf_height_dist_height(io_si,i_heightbin) + &
                       ccohort%c_area * AREA_INV * ccohort%treelai * frac_canopy_in_bin

                  ! if ( ( ccohort%c_area * AREA_INV * ccohort%treelai * frac_canopy_in_bin) .lt. 0._r8) then
                  !    write(fates_log(),*) ' negative hio_leaf_height_dist_height:'
                  !    write(fates_log(),*) '   c_area, treelai, frac_canopy_in_bin:', ccohort%c_area, ccohort%treelai, frac_canopy_in_bin
                  ! endif
               end do
               
               if (ccohort%canopy_layer .eq. 1) then
                  ! calculate the area of canopy that is within each height bin
                  hio_canopy_height_dist_height(io_si,height_bin_max) = &
                       hio_canopy_height_dist_height(io_si,height_bin_max) + ccohort%c_area * AREA_INV
               endif

               ! Update biomass components
               ! Mass pools [kgC]
               do el = 1, num_elements
            
                  sapw_m   = ccohort%prt%GetState(sapw_organ, element_list(el))
                  struct_m = ccohort%prt%GetState(struct_organ, element_list(el))
                  leaf_m   = ccohort%prt%GetState(leaf_organ, element_list(el))
                  fnrt_m   = ccohort%prt%GetState(fnrt_organ, element_list(el))
                  store_m  = ccohort%prt%GetState(store_organ, element_list(el))
                  repro_m  = ccohort%prt%GetState(repro_organ, element_list(el))
                  
                  alive_m  = leaf_m + fnrt_m + sapw_m
                  total_m  = alive_m + store_m + struct_m
                  
                  ! Plant multi-element states and fluxes
                  ! Zero states, and set the fluxes
                  if( element_list(el).eq.carbon12_element )then

                     this%hvars(ih_storec)%r81d(io_si)  = &
                          this%hvars(ih_storec)%r81d(io_si) + ccohort%n * store_m
                     this%hvars(ih_leafc)%r81d(io_si)   = &
                          this%hvars(ih_leafc)%r81d(io_si) + ccohort%n * leaf_m
                     this%hvars(ih_fnrtc)%r81d(io_si)   = &
                          this%hvars(ih_fnrtc)%r81d(io_si) + ccohort%n * fnrt_m
                     this%hvars(ih_reproc)%r81d(io_si)  = &
                          this%hvars(ih_reproc)%r81d(io_si)+ ccohort%n * repro_m
                     this%hvars(ih_sapwc)%r81d(io_si)   = &
                          this%hvars(ih_sapwc)%r81d(io_si)+ ccohort%n * sapw_m
                     this%hvars(ih_totvegc)%r81d(io_si) = &
                          this%hvars(ih_totvegc)%r81d(io_si)+ ccohort%n * total_m
                     
                     hio_agb(io_si)       = hio_agb(io_si) + n_perm2 * g_per_kg * &
                           ( leaf_m + (sapw_m + struct_m + store_m) * prt_params%allom_agb_frac(ccohort%pft) )
                     
                     
                     ! Update PFT partitioned biomass components
                     hio_leafbiomass_pft(io_si,ft) = hio_leafbiomass_pft(io_si,ft) + &
                           (ccohort%n * AREA_INV) * leaf_m     * g_per_kg
                     
                     hio_storebiomass_pft(io_si,ft) = hio_storebiomass_pft(io_si,ft) + &
                           (ccohort%n * AREA_INV) * store_m   * g_per_kg
                     
                     hio_nindivs_pft(io_si,ft) = hio_nindivs_pft(io_si,ft) + &
                           ccohort%n * AREA_INV
                     
                     hio_biomass_pft(io_si, ft) = hio_biomass_pft(io_si, ft) + &
                           (ccohort%n * AREA_INV) * total_m * g_per_kg

                     ! update total biomass per age bin
                     hio_biomass_age(io_si,cpatch%age_class) = hio_biomass_age(io_si,cpatch%age_class) &
                           + total_m * ccohort%n * AREA_INV
                     
                     ! track the total biomass on all secondary lands
                     if ( cpatch%anthro_disturbance_label .eq. secondaryforest ) then
                         hio_biomass_secondary_forest(io_si) = hio_biomass_secondary_forest(io_si) + &
                               total_m * ccohort%n * AREA_INV
                     endif
              
                  elseif(element_list(el).eq.nitrogen_element)then

                     this%hvars(ih_storen)%r81d(io_si)  = &
                          this%hvars(ih_storen)%r81d(io_si) + ccohort%n * store_m
                     this%hvars(ih_leafn)%r81d(io_si)   = &
                          this%hvars(ih_leafn)%r81d(io_si) + ccohort%n * leaf_m
                     this%hvars(ih_fnrtn)%r81d(io_si)   = &
                          this%hvars(ih_fnrtn)%r81d(io_si) + ccohort%n * fnrt_m
                     this%hvars(ih_repron)%r81d(io_si)  = &
                          this%hvars(ih_repron)%r81d(io_si) + ccohort%n * repro_m
                     this%hvars(ih_sapwn)%r81d(io_si)   = &
                          this%hvars(ih_sapwn)%r81d(io_si) + ccohort%n * sapw_m
                     this%hvars(ih_totvegn)%r81d(io_si) = &
                          this%hvars(ih_totvegn)%r81d(io_si) + ccohort%n * total_m

                     
                  elseif(element_list(el).eq.phosphorus_element) then
                     
                     this%hvars(ih_storep)%r81d(io_si)  = &
                          this%hvars(ih_storep)%r81d(io_si) + ccohort%n * store_m
                     this%hvars(ih_leafp)%r81d(io_si)   = &
                          this%hvars(ih_leafp)%r81d(io_si) + ccohort%n * leaf_m
                     this%hvars(ih_fnrtp)%r81d(io_si)   = &
                          this%hvars(ih_fnrtp)%r81d(io_si) + ccohort%n * fnrt_m
                     this%hvars(ih_reprop)%r81d(io_si)  = &
                          this%hvars(ih_reprop)%r81d(io_si) + ccohort%n * repro_m
                     this%hvars(ih_sapwp)%r81d(io_si)   = &
                          this%hvars(ih_sapwp)%r81d(io_si) + ccohort%n * sapw_m
                     this%hvars(ih_totvegp)%r81d(io_si) = &
                          this%hvars(ih_totvegp)%r81d(io_si)+ ccohort%n * total_m

                  end if
                     
               end do



               ! Update PFT crown area
               hio_crownarea_pft(io_si, ft) = hio_crownarea_pft(io_si, ft) + &
                    ccohort%c_area 

               if (ccohort%canopy_layer .eq. 1) then
                  ! Update PFT canopy crown area
                  hio_canopycrownarea_pft(io_si, ft) = hio_canopycrownarea_pft(io_si, ft) + &
                       ccohort%c_area 
               end if

               

               ! Site by Size-Class x PFT (SCPF) 
               ! ------------------------------------------------------------------------

               dbh = ccohort%dbh !-0.5*(1./365.25)*ccohort%ddbhdt

               ! Flux Variables (cohorts must had experienced a day before any of these values
               ! have any meaning, otherwise they are just inialization values
               if( .not.(ccohort%isnew) ) then

                  ! Turnover pools [kgC/day] * [day/yr] = [kgC/yr]
                  sapw_m_turnover   = ccohort%prt%GetTurnover(sapw_organ, carbon12_element) * days_per_year
                  store_m_turnover  = ccohort%prt%GetTurnover(store_organ, carbon12_element) * days_per_year
                  leaf_m_turnover   = ccohort%prt%GetTurnover(leaf_organ, carbon12_element) * days_per_year
                  fnrt_m_turnover   = ccohort%prt%GetTurnover(fnrt_organ, carbon12_element) * days_per_year
                  struct_m_turnover = ccohort%prt%GetTurnover(struct_organ, carbon12_element) * days_per_year
                  
                  ! Net change from allocation and transport [kgC/day] * [day/yr] = [kgC/yr]
                  sapw_m_net_alloc   = ccohort%prt%GetNetAlloc(sapw_organ, carbon12_element) * days_per_year
                  store_m_net_alloc  = ccohort%prt%GetNetAlloc(store_organ, carbon12_element) * days_per_year
                  leaf_m_net_alloc   = ccohort%prt%GetNetAlloc(leaf_organ, carbon12_element) * days_per_year
                  fnrt_m_net_alloc   = ccohort%prt%GetNetAlloc(fnrt_organ, carbon12_element) * days_per_year
                  struct_m_net_alloc = ccohort%prt%GetNetAlloc(struct_organ, carbon12_element) * days_per_year
                  repro_m_net_alloc  = ccohort%prt%GetNetAlloc(repro_organ, carbon12_element) * days_per_year

                  ! ecosystem-level, organ-partitioned NPP/allocation fluxes
                  hio_npp_leaf(io_si) = hio_npp_leaf(io_si) + leaf_m_net_alloc * n_perm2
                  hio_npp_seed(io_si) = hio_npp_seed(io_si) + repro_m_net_alloc * n_perm2
                  hio_npp_stem(io_si) = hio_npp_stem(io_si) + (sapw_m_net_alloc + struct_m_net_alloc) * n_perm2 * &
                       (prt_params%allom_agb_frac(ccohort%pft))
                  hio_npp_froot(io_si) = hio_npp_froot(io_si) + fnrt_m_net_alloc * n_perm2
                  hio_npp_croot(io_si) = hio_npp_croot(io_si) + (sapw_m_net_alloc + struct_m_net_alloc) * n_perm2 * &
                       (1._r8-prt_params%allom_agb_frac(ccohort%pft))
                  hio_npp_stor(io_si) = hio_npp_stor(io_si) + store_m_net_alloc * n_perm2
                  
                  associate( scpf => ccohort%size_by_pft_class, &

                       sz => ccohort%size_class, &
                       cacls => ccohort%coage_class, &
                       capf => ccohort%coage_by_pft_class)
     
			     
                    gpp_cached = hio_gpp_scpf(io_si,scpf)
      
                    hio_gpp_scpf(io_si,scpf)      = hio_gpp_scpf(io_si,scpf)      + &
                                                       n_perm2*ccohort%gpp_acc_hold  ! [kgC/m2/yr]
                    hio_npp_totl_scpf(io_si,scpf) = hio_npp_totl_scpf(io_si,scpf) + &
                                                       ccohort%npp_acc_hold *n_perm2
                    
                    
                    hio_npp_leaf_scpf(io_si,scpf) = hio_npp_leaf_scpf(io_si,scpf) + &
                                                       leaf_m_net_alloc*n_perm2
                    hio_npp_fnrt_scpf(io_si,scpf) = hio_npp_fnrt_scpf(io_si,scpf) + &
                                                       fnrt_m_net_alloc*n_perm2
                    hio_npp_bgsw_scpf(io_si,scpf) = hio_npp_bgsw_scpf(io_si,scpf) + &
                                                       sapw_m_net_alloc*n_perm2*           &
                                                       (1._r8-prt_params%allom_agb_frac(ccohort%pft))
                    hio_npp_agsw_scpf(io_si,scpf) = hio_npp_agsw_scpf(io_si,scpf) + &
                                                       sapw_m_net_alloc*n_perm2*           &
                                                       prt_params%allom_agb_frac(ccohort%pft)
                    hio_npp_bgdw_scpf(io_si,scpf) = hio_npp_bgdw_scpf(io_si,scpf) + &
                                                       struct_m_net_alloc*n_perm2*         &
                                                       (1._r8-prt_params%allom_agb_frac(ccohort%pft))
                    hio_npp_agdw_scpf(io_si,scpf) = hio_npp_agdw_scpf(io_si,scpf) + &
                                                       struct_m_net_alloc*n_perm2*         &
                                                       prt_params%allom_agb_frac(ccohort%pft)
                    hio_npp_seed_scpf(io_si,scpf) = hio_npp_seed_scpf(io_si,scpf) + &
                                                       repro_m_net_alloc*n_perm2
                    hio_npp_stor_scpf(io_si,scpf) = hio_npp_stor_scpf(io_si,scpf) + &
                                                       store_m_net_alloc*n_perm2

                    ! Woody State Variables (basal area growth increment)
                    if ( int(prt_params%woody(ft)) == itrue) then

                       ! basal area  [m2/ha]
                       hio_ba_scpf(io_si,scpf) = hio_ba_scpf(io_si,scpf) + &
                            0.25_r8*3.14159_r8*((dbh/100.0_r8)**2.0_r8)*ccohort%n

                       ! also by size class only
                       hio_ba_sz(io_si,sz) = hio_ba_sz(io_si,sz) + &
                            0.25_r8*3.14159_r8*((dbh/100.0_r8)**2.0_r8)*ccohort%n

                       ! growth increment
                       hio_ddbh_scpf(io_si,scpf) = hio_ddbh_scpf(io_si,scpf) + &
                            ccohort%ddbhdt*ccohort%n

                    end if

                    hio_m1_scpf(io_si,scpf) = hio_m1_scpf(io_si,scpf) + ccohort%bmort*ccohort%n
                    hio_m2_scpf(io_si,scpf) = hio_m2_scpf(io_si,scpf) + ccohort%hmort*ccohort%n
                    hio_m3_scpf(io_si,scpf) = hio_m3_scpf(io_si,scpf) + ccohort%cmort*ccohort%n
                    hio_m7_scpf(io_si,scpf) = hio_m7_scpf(io_si,scpf) + &
                         (ccohort%lmort_direct+ccohort%lmort_collateral+ccohort%lmort_infra) * ccohort%n
                    hio_m8_scpf(io_si,scpf) = hio_m8_scpf(io_si,scpf) + ccohort%frmort*ccohort%n
                    hio_m9_scpf(io_si,scpf) = hio_m9_scpf(io_si,scpf) + ccohort%smort*ccohort%n
                    
                    if (hlm_use_cohort_age_tracking .eq.itrue) then
                       hio_m10_scpf(io_si,scpf) = hio_m10_scpf(io_si,scpf) + ccohort%asmort*ccohort%n
                       hio_m10_capf(io_si,capf) = hio_m10_capf(io_si,capf) + ccohort%asmort*ccohort%n
                       hio_m10_sz(io_si,sz) = hio_m10_sz(io_si,sz) + ccohort%asmort*ccohort%n
                       hio_m10_cacls(io_si,cacls) = hio_m10_cacls(io_si,cacls)+ &
                            ccohort%asmort*ccohort%n
                    end if
                    
                    hio_m1_sz(io_si,sz) = hio_m1_sz(io_si,sz) + ccohort%bmort*ccohort%n
                    hio_m2_sz(io_si,sz) = hio_m2_sz(io_si,sz) + ccohort%hmort*ccohort%n
                    hio_m3_sz(io_si,sz) = hio_m3_sz(io_si,sz) + ccohort%cmort*ccohort%n
                    hio_m7_sz(io_si,sz) = hio_m7_sz(io_si,sz) + &
                         (ccohort%lmort_direct+ccohort%lmort_collateral+ccohort%lmort_infra) * ccohort%n
                    hio_m8_sz(io_si,sz) = hio_m8_sz(io_si,sz) + &
                         ccohort%frmort*ccohort%n
                    hio_m9_sz(io_si,sz) = hio_m9_sz(io_si,sz) + ccohort%smort*ccohort%n
                  
                  
                   
                    !C13 discrimination
                    if(gpp_cached + ccohort%gpp_acc_hold > 0.0_r8)then
                       hio_c13disc_scpf(io_si,scpf) = ((hio_c13disc_scpf(io_si,scpf) * gpp_cached) + &
                            (ccohort%c13disc_acc * ccohort%gpp_acc_hold)) / (gpp_cached + ccohort%gpp_acc_hold)
                    else
                       hio_c13disc_scpf(io_si,scpf) = 0.0_r8
                    endif

                    ! number density [/ha]
                    hio_nplant_scpf(io_si,scpf) = hio_nplant_scpf(io_si,scpf) + ccohort%n

                    ! number density along the cohort age dimension
                    if (hlm_use_cohort_age_tracking .eq.itrue) then
                       hio_nplant_capf(io_si,capf) = hio_nplant_capf(io_si,capf) + ccohort%n
                       hio_nplant_cacls(io_si,cacls) = hio_nplant_cacls(io_si,cacls)+ccohort%n
                    end if


                    ! Carbon only metrics
                    sapw_m   = ccohort%prt%GetState(sapw_organ, carbon12_element)
                    struct_m = ccohort%prt%GetState(struct_organ, carbon12_element)
                    leaf_m   = ccohort%prt%GetState(leaf_organ, carbon12_element)
                    fnrt_m   = ccohort%prt%GetState(fnrt_organ, carbon12_element)
                    store_m  = ccohort%prt%GetState(store_organ, carbon12_element)
                    repro_m  = ccohort%prt%GetState(repro_organ, carbon12_element)
                    alive_m  = leaf_m + fnrt_m + sapw_m
                    total_m  = alive_m + store_m + struct_m
                    
                    
                    ! number density by size and biomass
                    hio_agb_sz(io_si,sz) = hio_agb_sz(io_si,sz) + &
                          total_m * ccohort%n * prt_params%allom_agb_frac(ccohort%pft) * AREA_INV

                    hio_agb_scpf(io_si,scpf) = hio_agb_scpf(io_si,scpf) + &
                         total_m * ccohort%n * prt_params%allom_agb_frac(ccohort%pft) * AREA_INV


                    hio_biomass_sz(io_si,sz) = hio_biomass_sz(io_si,sz) + &
                          total_m * ccohort%n * AREA_INV

                    ! update size-class x patch-age related quantities

                    iscag = get_sizeage_class_index(ccohort%dbh,cpatch%age)
                    
                    hio_nplant_scag(io_si,iscag) = hio_nplant_scag(io_si,iscag) + ccohort%n

                    hio_nplant_sz(io_si,sz) = hio_nplant_sz(io_si,sz) + ccohort%n
                    
                  
                    ! update size, age, and PFT - indexed quantities

                    iscagpft = get_sizeagepft_class_index(ccohort%dbh,cpatch%age,ccohort%pft)
                    
                    hio_nplant_scagpft(io_si,iscagpft) = hio_nplant_scagpft(io_si,iscagpft) + ccohort%n

                    ! update age and PFT - indexed quantities

                    iagepft = get_agepft_class_index(cpatch%age,ccohort%pft)
                    
                    hio_npp_agepft(io_si,iagepft) = hio_npp_agepft(io_si,iagepft) + &
                         ccohort%n * ccohort%npp_acc_hold * AREA_INV

                    hio_biomass_agepft(io_si,iagepft) = hio_biomass_agepft(io_si,iagepft) + &
                          total_m * ccohort%n * AREA_INV

                    ! update SCPF/SZ- and canopy/subcanopy- partitioned quantities
                    if (ccohort%canopy_layer .eq. 1) then
                       hio_nplant_canopy_scag(io_si,iscag) = hio_nplant_canopy_scag(io_si,iscag) + ccohort%n
                       hio_mortality_canopy_scag(io_si,iscag) = hio_mortality_canopy_scag(io_si,iscag) + &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + & 
                            ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n
                       hio_ddbh_canopy_scag(io_si,iscag) = hio_ddbh_canopy_scag(io_si,iscag) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_bstor_canopy_scpf(io_si,scpf) = hio_bstor_canopy_scpf(io_si,scpf) + &
                             store_m * ccohort%n
                       hio_bleaf_canopy_scpf(io_si,scpf) = hio_bleaf_canopy_scpf(io_si,scpf) + &
                             leaf_m * ccohort%n

                       hio_canopy_biomass(io_si) = hio_canopy_biomass(io_si) + n_perm2 * total_m * g_per_kg

                       !hio_mortality_canopy_scpf(io_si,scpf) = hio_mortality_canopy_scpf(io_si,scpf)+ &
                       !    (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                       ! ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n

                       hio_mortality_canopy_scpf(io_si,scpf) = hio_mortality_canopy_scpf(io_si,scpf)+ &

                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%frmort + & 
                            ccohort%smort + ccohort%asmort) * ccohort%n + &
                            (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                            ccohort%n * sec_per_day * days_per_year

                       hio_nplant_canopy_scpf(io_si,scpf) = hio_nplant_canopy_scpf(io_si,scpf) + ccohort%n
                       hio_nplant_canopy_sz(io_si,sz) = hio_nplant_canopy_sz(io_si,sz) + ccohort%n
                       hio_lai_canopy_sz(io_si,sz) = hio_lai_canopy_sz(io_si,sz) + &
                                                            ccohort%treelai*ccohort%c_area * AREA_INV
                       hio_sai_canopy_sz(io_si,sz) = hio_sai_canopy_sz(io_si,sz) + &
                                                            ccohort%treesai*ccohort%c_area * AREA_INV
                       hio_trimming_canopy_sz(io_si,sz) = hio_trimming_canopy_sz(io_si,sz) + &
                            ccohort%n * ccohort%canopy_trim
                       hio_crown_area_canopy_sz(io_si,sz) = hio_crown_area_canopy_sz(io_si,sz) + &
                            ccohort%c_area
                       hio_gpp_canopy_scpf(io_si,scpf)      = hio_gpp_canopy_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%gpp_acc_hold
                       hio_ar_canopy_scpf(io_si,scpf)      = hio_ar_canopy_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%resp_acc_hold
                       ! growth increment
                       hio_ddbh_canopy_scpf(io_si,scpf) = hio_ddbh_canopy_scpf(io_si,scpf) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_ddbh_canopy_sz(io_si,sz) = hio_ddbh_canopy_sz(io_si,sz) + &
                            ccohort%ddbhdt*ccohort%n

                       ! sum of all mortality
                       hio_mortality_canopy_sz(io_si,sz) = hio_mortality_canopy_sz(io_si,sz) + &

                             (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                             ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n + &
                             (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                             ccohort%n * sec_per_day * days_per_year

                       hio_canopy_mortality_carbonflux(io_si) = hio_canopy_mortality_carbonflux(io_si) + &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + & 
                            ccohort%frmort + ccohort%smort + ccohort%asmort) * &
                            total_m * ccohort%n * g_per_kg * days_per_sec * years_per_day * ha_per_m2 + &
                            (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * total_m * &
                            ccohort%n * g_per_kg * ha_per_m2
                       

                       hio_carbon_balance_canopy_sz(io_si,sz) = hio_carbon_balance_canopy_sz(io_si,sz) + &
                             ccohort%n * ccohort%npp_acc_hold
                       
                      
                       hio_leaf_md_canopy_sz(io_si,sz) = hio_leaf_md_canopy_sz(io_si,sz) + &
                             leaf_m_turnover * ccohort%n
                       hio_root_md_canopy_sz(io_si,sz) = hio_root_md_canopy_sz(io_si,sz) + &
                             fnrt_m_turnover * ccohort%n
                       hio_bsw_md_canopy_sz(io_si,sz) = hio_bsw_md_canopy_sz(io_si,sz) + &
                             sapw_m_turnover * ccohort%n
                       hio_bstore_md_canopy_sz(io_si,sz) = hio_bstore_md_canopy_sz(io_si,sz) + &
                             store_m_turnover * ccohort%n
                       hio_bdead_md_canopy_sz(io_si,sz) = hio_bdead_md_canopy_sz(io_si,sz) + &
                             struct_m_turnover * ccohort%n
                       hio_seed_prod_canopy_sz(io_si,sz) = hio_seed_prod_canopy_sz(io_si,sz) + &
                             ccohort%seed_prod * ccohort%n

                       hio_npp_leaf_canopy_sz(io_si,sz) = hio_npp_leaf_canopy_sz(io_si,sz) + &
                             leaf_m_net_alloc * ccohort%n
                       hio_npp_fnrt_canopy_sz(io_si,sz) = hio_npp_fnrt_canopy_sz(io_si,sz) + &
                             fnrt_m_net_alloc * ccohort%n
                       hio_npp_sapw_canopy_sz(io_si,sz) = hio_npp_sapw_canopy_sz(io_si,sz) + &
                             sapw_m_net_alloc * ccohort%n
                       hio_npp_dead_canopy_sz(io_si,sz) = hio_npp_dead_canopy_sz(io_si,sz) + &
                             struct_m_net_alloc * ccohort%n
                       hio_npp_seed_canopy_sz(io_si,sz) = hio_npp_seed_canopy_sz(io_si,sz) + &
                             repro_m_net_alloc * ccohort%n
                       hio_npp_stor_canopy_sz(io_si,sz) = hio_npp_stor_canopy_sz(io_si,sz) + &
                             store_m_net_alloc * ccohort%n
                       
                       hio_yesterdaycanopylevel_canopy_sz(io_si,sz) = &
                            hio_yesterdaycanopylevel_canopy_sz(io_si,sz) + &
                            ccohort%canopy_layer_yesterday * ccohort%n
                    else
                       hio_nplant_understory_scag(io_si,iscag) = hio_nplant_understory_scag(io_si,iscag) + ccohort%n
                       hio_mortality_understory_scag(io_si,iscag) = hio_mortality_understory_scag(io_si,iscag) + &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                            ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n
                       hio_ddbh_understory_scag(io_si,iscag) = hio_ddbh_understory_scag(io_si,iscag) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_bstor_understory_scpf(io_si,scpf) = hio_bstor_understory_scpf(io_si,scpf) + &
                             store_m * ccohort%n
                       hio_bleaf_understory_scpf(io_si,scpf) = hio_bleaf_understory_scpf(io_si,scpf) + &
                             leaf_m  * ccohort%n
                       hio_understory_biomass(io_si) = hio_understory_biomass(io_si) + &
                             n_perm2 * total_m * g_per_kg

                       !hio_mortality_understory_scpf(io_si,scpf) = hio_mortality_understory_scpf(io_si,scpf)+ &
                        !    (ccohort%bmort + ccohort%hmort + ccohort%cmort + 
                       !      ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n

                       hio_mortality_understory_scpf(io_si,scpf) = hio_mortality_understory_scpf(io_si,scpf)+ &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + &
                            ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n + &
                            (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                            ccohort%n * sec_per_day * days_per_year

                       hio_nplant_understory_scpf(io_si,scpf) = hio_nplant_understory_scpf(io_si,scpf) + ccohort%n
                       hio_nplant_understory_sz(io_si,sz) = hio_nplant_understory_sz(io_si,sz) + ccohort%n
                       hio_lai_understory_sz(io_si,sz) = hio_lai_understory_sz(io_si,sz) + &
                                                                ccohort%treelai*ccohort%c_area  * AREA_INV
                       hio_sai_understory_sz(io_si,sz) = hio_sai_understory_sz(io_si,sz) + &
                                                                ccohort%treelai*ccohort%c_area  * AREA_INV
                       hio_trimming_understory_sz(io_si,sz) = hio_trimming_understory_sz(io_si,sz) + &
                            ccohort%n * ccohort%canopy_trim
                       hio_crown_area_understory_sz(io_si,sz) = hio_crown_area_understory_sz(io_si,sz) + &
                            ccohort%c_area
                       hio_gpp_understory_scpf(io_si,scpf)      = hio_gpp_understory_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%gpp_acc_hold
                       hio_ar_understory_scpf(io_si,scpf)      = hio_ar_understory_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%resp_acc_hold

                       ! growth increment
                       hio_ddbh_understory_scpf(io_si,scpf) = hio_ddbh_understory_scpf(io_si,scpf) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_ddbh_understory_sz(io_si,sz) = hio_ddbh_understory_sz(io_si,sz) + &
                            ccohort%ddbhdt*ccohort%n

                       ! sum of all mortality
                       hio_mortality_understory_sz(io_si,sz) = hio_mortality_understory_sz(io_si,sz) + &

                             (ccohort%bmort + ccohort%hmort + ccohort%cmort + & 
                             ccohort%frmort + ccohort%smort + ccohort%asmort) * ccohort%n + &
                             (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * &
                             ccohort%n * sec_per_day * days_per_year
                       
                       hio_understory_mortality_carbonflux(io_si) = hio_understory_mortality_carbonflux(io_si) + &
                             (ccohort%bmort + ccohort%hmort + ccohort%cmort + & 
                             ccohort%frmort + ccohort%smort + ccohort%asmort) * &
                             total_m * ccohort%n * g_per_kg * days_per_sec * years_per_day * ha_per_m2 + &
                             (ccohort%lmort_direct + ccohort%lmort_collateral + ccohort%lmort_infra) * total_m * &
                             ccohort%n * g_per_kg * ha_per_m2

                       hio_carbon_balance_understory_sz(io_si,sz) = hio_carbon_balance_understory_sz(io_si,sz) + &
                             ccohort%npp_acc_hold * ccohort%n

                       hio_leaf_md_understory_sz(io_si,sz) = hio_leaf_md_understory_sz(io_si,sz) + &
                            leaf_m_turnover * ccohort%n
                       hio_root_md_understory_sz(io_si,sz) = hio_root_md_understory_sz(io_si,sz) + &
                            fnrt_m_turnover * ccohort%n
                       hio_bsw_md_understory_sz(io_si,sz) = hio_bsw_md_understory_sz(io_si,sz) + &
                             sapw_m_turnover * ccohort%n
                       hio_bstore_md_understory_sz(io_si,sz) = hio_bstore_md_understory_sz(io_si,sz) + &
                             store_m_turnover * ccohort%n
                       hio_bdead_md_understory_sz(io_si,sz) = hio_bdead_md_understory_sz(io_si,sz) + &
                             struct_m_turnover * ccohort%n
                       hio_seed_prod_understory_sz(io_si,sz) = hio_seed_prod_understory_sz(io_si,sz) + &
                            ccohort%seed_prod * ccohort%n

                       hio_npp_leaf_understory_sz(io_si,sz) = hio_npp_leaf_understory_sz(io_si,sz) + &
                             leaf_m_net_alloc * ccohort%n
                       hio_npp_fnrt_understory_sz(io_si,sz) = hio_npp_fnrt_understory_sz(io_si,sz) + &
                             fnrt_m_net_alloc * ccohort%n
                       hio_npp_sapw_understory_sz(io_si,sz) = hio_npp_sapw_understory_sz(io_si,sz) + &
                             sapw_m_net_alloc * ccohort%n
                       hio_npp_dead_understory_sz(io_si,sz) = hio_npp_dead_understory_sz(io_si,sz) + &
                             struct_m_net_alloc * ccohort%n
                       hio_npp_seed_understory_sz(io_si,sz) = hio_npp_seed_understory_sz(io_si,sz) + &
                             repro_m_net_alloc * ccohort%n
                       hio_npp_stor_understory_sz(io_si,sz) = hio_npp_stor_understory_sz(io_si,sz) + &
                             store_m_net_alloc * ccohort%n
                       
                       hio_yesterdaycanopylevel_understory_sz(io_si,sz) = &
                            hio_yesterdaycanopylevel_understory_sz(io_si,sz) + &
                            ccohort%canopy_layer_yesterday * ccohort%n
                    endif
                    !
                    !
                    ccohort%canopy_layer_yesterday = real(ccohort%canopy_layer, r8)
                    !
                    ! growth flux of individuals into a given bin
                    ! track the actual growth here, the virtual growth from fusion lower down
                    if ( (sz - ccohort%size_class_lasttimestep ) .gt. 0) then
                       do i_sz = ccohort%size_class_lasttimestep + 1, sz
                          i_scpf = (ccohort%pft-1)*nlevsclass+i_sz
                          hio_growthflux_scpf(io_si,i_scpf) = hio_growthflux_scpf(io_si,i_scpf) + &
                               ccohort%n * days_per_year
                       end do
                    end if
                    ccohort%size_class_lasttimestep = sz


                    !
                  end associate
               else  ! i.e. cohort%isnew
                  !
                  ! if cohort is new, track its growth flux into the first size bin
                  i_scpf = (ccohort%pft-1)*nlevsclass+1
                  hio_growthflux_scpf(io_si,i_scpf) = hio_growthflux_scpf(io_si,i_scpf) + ccohort%n * days_per_year
                  ccohort%size_class_lasttimestep = 1
                                   
               end if

               ! resolve some canopy area profiles, both total and of occupied leaves
               ican = ccohort%canopy_layer
               !
               hio_crownarea_can(io_si, ican) = hio_crownarea_can(io_si, ican) + ccohort%c_area / AREA
               !
               do ileaf=1,ccohort%nv
                  cnlf_indx = ileaf + (ican-1) * nlevleaf
                  hio_crownarea_cnlf(io_si, cnlf_indx) = hio_crownarea_cnlf(io_si, cnlf_indx) + &
                       ccohort%c_area / AREA
               end do
               
               ccohort => ccohort%taller
            enddo ! cohort loop
            
            ! Patch specific variables that are already calculated
            ! These things are all duplicated. Should they all be converted to LL or array structures RF? 
            ! define scalar to counteract the patch albedo scaling logic for conserved quantities
                        
            ! Update Fire Variables
            hio_spitfire_ros(io_si)         = hio_spitfire_ros(io_si) + cpatch%ROS_front * cpatch%area * AREA_INV
            hio_fire_ros_area_product(io_si)= hio_fire_ros_area_product(io_si) + &
                 cpatch%frac_burnt * cpatch%ROS_front * cpatch%area * AREA_INV
            hio_effect_wspeed(io_si)        = hio_effect_wspeed(io_si) + cpatch%effect_wspeed * cpatch%area * AREA_INV
            hio_tfc_ros(io_si)              = hio_tfc_ros(io_si) + cpatch%TFC_ROS * cpatch%area * AREA_INV
            hio_tfc_ros_area_product(io_si) = hio_tfc_ros_area_product(io_si) + &
                 cpatch%frac_burnt * cpatch%TFC_ROS * cpatch%area * AREA_INV
            hio_fire_intensity(io_si)       = hio_fire_intensity(io_si) + cpatch%FI * cpatch%area * AREA_INV
            hio_fire_area(io_si)            = hio_fire_area(io_si) + cpatch%frac_burnt * cpatch%area * AREA_INV
            hio_fire_fuel_bulkd(io_si)      = hio_fire_fuel_bulkd(io_si) + cpatch%fuel_bulkd * cpatch%area * AREA_INV
            hio_fire_fuel_eff_moist(io_si)  = hio_fire_fuel_eff_moist(io_si) + cpatch%fuel_eff_moist * cpatch%area * AREA_INV
            hio_fire_fuel_sav(io_si)        = hio_fire_fuel_sav(io_si) + cpatch%fuel_sav * cpatch%area * AREA_INV
            hio_fire_fuel_mef(io_si)        = hio_fire_fuel_mef(io_si) + cpatch%fuel_mef * cpatch%area * AREA_INV
            hio_sum_fuel(io_si)             = hio_sum_fuel(io_si) + cpatch%sum_fuel * g_per_kg * cpatch%area * AREA_INV
            hio_fragmentation_scaler(io_si) = hio_fragmentation_scaler(io_si) + cpatch%fragmentation_scaler * cpatch%area * AREA_INV
            
            do i_fuel = 1,nfsc
               hio_litter_moisture_fuel(io_si, i_fuel) = hio_litter_moisture_fuel(io_si, i_fuel) + &
                    cpatch%litter_moisture(i_fuel) * cpatch%area * AREA_INV

               hio_fuel_amount_fuel(io_si, i_fuel) = hio_fuel_amount_fuel(io_si, i_fuel) + &
                    cpatch%fuel_frac(i_fuel) * cpatch%sum_fuel * cpatch%area * AREA_INV

               hio_burnt_frac_litter_fuel(io_si, i_fuel) = hio_burnt_frac_litter_fuel(io_si, i_fuel) + &
                    cpatch%burnt_frac_litter(i_fuel) * cpatch%frac_burnt * cpatch%area * AREA_INV
            end do


            hio_fire_intensity_area_product(io_si) = hio_fire_intensity_area_product(io_si) + &
                 cpatch%FI * cpatch%frac_burnt * cpatch%area * AREA_INV

            ! Update Litter Flux Variables

            litt_c       => cpatch%litter(element_pos(carbon12_element))
            flux_diags_c => sites(s)%flux_diags(element_pos(carbon12_element))
                         
            do i_cwd = 1, ncwd

                hio_cwd_ag_cwdsc(io_si, i_cwd) = hio_cwd_ag_cwdsc(io_si, i_cwd) + &
                      litt_c%ag_cwd(i_cwd)*cpatch%area * AREA_INV * g_per_kg
                hio_cwd_bg_cwdsc(io_si, i_cwd) = hio_cwd_bg_cwdsc(io_si, i_cwd) + &
                      sum(litt_c%bg_cwd(i_cwd,:)) * cpatch%area * AREA_INV * g_per_kg
                
                hio_cwd_ag_out_cwdsc(io_si, i_cwd) = hio_cwd_ag_out_cwdsc(io_si, i_cwd) + &
                      litt_c%ag_cwd_frag(i_cwd)*cpatch%area * AREA_INV * g_per_kg
                
                hio_cwd_bg_out_cwdsc(io_si, i_cwd) = hio_cwd_bg_out_cwdsc(io_si, i_cwd) + &
                      sum(litt_c%bg_cwd_frag(i_cwd,:)) * cpatch%area * AREA_INV * g_per_kg

            end do

            ipa = ipa + 1
            cpatch => cpatch%younger
         end do !patch loop

         ! divide so-far-just-summed but to-be-averaged patch-age-class variables by patch-age-class area to get mean values
         do ipa2 = 1, nlevage
            if (hio_area_age(io_si, ipa2) .gt. tiny) then
               hio_lai_age(io_si, ipa2) = hio_lai_age(io_si, ipa2) / (hio_area_age(io_si, ipa2)*AREA)
               hio_ncl_age(io_si, ipa2) = hio_ncl_age(io_si, ipa2) / (hio_area_age(io_si, ipa2)*AREA)
               do i_pft = 1, numpft
                  iagepft = ipa2 + (i_pft-1) * nlevage
                  hio_scorch_height_agepft(io_si, iagepft) = &
                       hio_scorch_height_agepft(io_si, iagepft) / (hio_area_age(io_si, ipa2)*AREA)
               enddo
            else
               hio_lai_age(io_si, ipa2) = 0._r8
               hio_ncl_age(io_si, ipa2) = 0._r8
            endif
         end do

         ! pass the cohort termination mortality as a flux to the history, and then reset the termination mortality buffer
         ! note there are various ways of reporting the total mortality, so pass to these as well
         do i_pft = 1, numpft
            do i_sz = 1,nlevsclass
               i_scpf = (i_pft-1)*nlevsclass + i_sz
               !
               ! termination mortality. sum of canopy and understory indices
               hio_m6_scpf(io_si,i_scpf) = (sites(s)%term_nindivs_canopy(i_sz,i_pft) + &
                                               sites(s)%term_nindivs_ustory(i_sz,i_pft)) * days_per_year

               hio_m6_sz(io_si,i_sz) = hio_m6_sz(io_si,i_sz) + &
                     (sites(s)%term_nindivs_canopy(i_sz,i_pft) + &
                      sites(s)%term_nindivs_ustory(i_sz,i_pft)) * days_per_year
                     

               !
               ! add termination mortality to canopy and understory mortality
               hio_mortality_canopy_sz(io_si,i_sz) = hio_mortality_canopy_sz(io_si,i_sz) + &
                    sites(s)%term_nindivs_canopy(i_sz,i_pft) * days_per_year

               hio_mortality_understory_sz(io_si,i_sz) = hio_mortality_understory_sz(io_si,i_sz) + &
                    sites(s)%term_nindivs_ustory(i_sz,i_pft) * days_per_year

               hio_mortality_canopy_scpf(io_si,i_scpf) = hio_mortality_canopy_scpf(io_si,i_scpf) + &
                     sites(s)%term_nindivs_canopy(i_sz,i_pft) * days_per_year

               hio_mortality_understory_scpf(io_si,i_scpf) = hio_mortality_understory_scpf(io_si,i_scpf) + &
                     sites(s)%term_nindivs_ustory(i_sz,i_pft) * days_per_year

               !
               ! imort on its own
               hio_m4_scpf(io_si,i_scpf) = sites(s)%imort_rate(i_sz, i_pft)
               hio_m4_sz(io_si,i_sz) = hio_m4_sz(io_si,i_sz) + sites(s)%imort_rate(i_sz, i_pft)
               !
               ! add imort to other mortality terms. consider imort as understory mortality even if it happens in 
               ! cohorts that may have been promoted as part of the patch creation, and use the pre-calculated site-level 
               ! values to avoid biasing the results by the dramatically-reduced number densities in cohorts that are subject to imort
               hio_mortality_understory_scpf(io_si,i_scpf) = hio_mortality_understory_scpf(io_si,i_scpf) + &
                    sites(s)%imort_rate(i_sz, i_pft)
               hio_mortality_understory_sz(io_si,i_sz) = hio_mortality_understory_sz(io_si,i_sz) + &
                    sites(s)%imort_rate(i_sz, i_pft)
               !
               iscag = i_sz ! since imort is by definition something that only happens in newly disturbed patches, treat as such
               hio_mortality_understory_scag(io_si,iscag) = hio_mortality_understory_scag(io_si,iscag) + &
                    sites(s)%imort_rate(i_sz, i_pft)

               ! fire mortality from the site-level diagnostic rates
               hio_m5_scpf(io_si,i_scpf) = sites(s)%fmort_rate_canopy(i_sz, i_pft) + &
                     sites(s)%fmort_rate_ustory(i_sz, i_pft)
               hio_m5_sz(io_si,i_sz) = hio_m5_sz(io_si,i_sz) + &
                     sites(s)%fmort_rate_canopy(i_sz, i_pft) +  sites(s)%fmort_rate_ustory(i_sz, i_pft)
               !
               hio_crownfiremort_scpf(io_si,i_scpf) = sites(s)%fmort_rate_crown(i_sz, i_pft)
               hio_cambialfiremort_scpf(io_si,i_scpf) = sites(s)%fmort_rate_cambial(i_sz, i_pft)
               !
               ! fire components of overall canopy and understory mortality
               hio_mortality_canopy_scpf(io_si,i_scpf) = hio_mortality_canopy_scpf(io_si,i_scpf) + &
                    sites(s)%fmort_rate_canopy(i_sz, i_pft)
               hio_mortality_canopy_sz(io_si,i_sz) = hio_mortality_canopy_sz(io_si,i_sz) + &
                    sites(s)%fmort_rate_canopy(i_sz, i_pft)

               ! the fire mortality rates for each layer are total dead, since the usable
               ! output will then normalize by the counts, we are allowed to sum over layers
               hio_mortality_understory_scpf(io_si,i_scpf) = hio_mortality_understory_scpf(io_si,i_scpf) + &
                     sites(s)%fmort_rate_ustory(i_sz, i_pft)

               hio_mortality_understory_sz(io_si,i_sz) = hio_mortality_understory_sz(io_si,i_sz) + &
                     sites(s)%fmort_rate_ustory(i_sz, i_pft)

               !
               ! carbon flux associated with mortality of trees dying by fire
               hio_canopy_mortality_carbonflux(io_si) = hio_canopy_mortality_carbonflux(io_si) + &
                     sites(s)%fmort_carbonflux_canopy
               
               hio_understory_mortality_carbonflux(io_si) = hio_understory_mortality_carbonflux(io_si) + &
                     sites(s)%fmort_carbonflux_ustory
               
               !
               ! for scag variables, also treat as happening in the newly-disurbed patch

               hio_mortality_canopy_scag(io_si,iscag) = hio_mortality_canopy_scag(io_si,iscag) + &
                    sites(s)%fmort_rate_canopy(i_sz, i_pft)
               hio_mortality_understory_scag(io_si,iscag) = hio_mortality_understory_scag(io_si,iscag) + &
                    sites(s)%fmort_rate_ustory(i_sz, i_pft)

               ! while in this loop, pass the fusion-induced growth rate flux to history
               hio_growthflux_fusion_scpf(io_si,i_scpf) = hio_growthflux_fusion_scpf(io_si,i_scpf) + &
                    sites(s)%growthflux_fusion(i_sz, i_pft) * days_per_year
            end do
         end do
         !
         
         ! treat carbon flux from imort the same way
         hio_understory_mortality_carbonflux(io_si) = hio_understory_mortality_carbonflux(io_si) + &
              sites(s)%imort_carbonflux
         !
         sites(s)%term_nindivs_canopy(:,:) = 0._r8
         sites(s)%term_nindivs_ustory(:,:) = 0._r8
         sites(s)%imort_carbonflux = 0._r8
         sites(s)%imort_rate(:,:) = 0._r8
         sites(s)%fmort_rate_canopy(:,:) = 0._r8
         sites(s)%fmort_rate_ustory(:,:) = 0._r8
         sites(s)%fmort_carbonflux_canopy = 0._r8
         sites(s)%fmort_carbonflux_ustory = 0._r8
         sites(s)%fmort_rate_cambial(:,:) = 0._r8
         sites(s)%fmort_rate_crown(:,:) = 0._r8
         sites(s)%growthflux_fusion(:,:) = 0._r8

         ! pass the recruitment rate as a flux to the history, and then reset the recruitment buffer
         do i_pft = 1, numpft
            hio_recruitment_pft(io_si,i_pft) = sites(s)%recruitment_rate(i_pft) * days_per_year
         end do
         sites(s)%recruitment_rate(:) = 0._r8

         ! summarize all of the mortality fluxes by PFT
         do i_pft = 1, numpft
            do i_sz = 1,nlevsclass
               i_scpf = (i_pft-1)*nlevsclass + i_sz

               hio_mortality_pft(io_si,i_pft) = hio_mortality_pft(io_si,i_pft) + &
                    hio_m1_scpf(io_si,i_scpf) + &
                    hio_m2_scpf(io_si,i_scpf) + &
                    hio_m3_scpf(io_si,i_scpf) + &
                    hio_m4_scpf(io_si,i_scpf) + &
                    hio_m5_scpf(io_si,i_scpf) + &
                    hio_m6_scpf(io_si,i_scpf) + &
		    hio_m7_scpf(io_si,i_scpf) + &
                    hio_m8_scpf(io_si,i_scpf) + &
                    hio_m9_scpf(io_si,i_scpf) + &
                    hio_m10_scpf(io_si,i_scpf) 

            end do
         end do
         
         ! ------------------------------------------------------------------------------
         ! Some carbon only litter diagnostics (legacy)
         ! ------------------------------------------------------------------------------

         flux_diags => sites(s)%flux_diags(element_pos(carbon12_element))

         hio_litter_in(io_si) = (sum(flux_diags%cwd_ag_input(:)) + &
              sum(flux_diags%cwd_bg_input(:)) + &
              sum(flux_diags%leaf_litter_input(:)) + &
              sum(flux_diags%root_litter_input(:))) * &
              g_per_kg * AREA_INV * days_per_sec

         hio_litter_out(io_si) = 0._r8
         hio_seed_bank(io_si)  = 0._r8
         hio_seeds_in(io_si)   = 0._r8

         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            
            litt => cpatch%litter(element_pos(carbon12_element))
            
            area_frac = cpatch%area * AREA_INV
            
            ! Sum up all output fluxes (fragmentation) kgC/m2/day -> gC/m2/s
            hio_litter_out(io_si) = hio_litter_out(io_si) + &
                 (sum(litt%leaf_fines_frag(:)) + &
                 sum(litt%root_fines_frag(:,:)) + &
                 sum(litt%ag_cwd_frag(:)) + &
                 sum(litt%bg_cwd_frag(:,:))) * &
                 area_frac * g_per_kg * days_per_sec

            ! Sum up total seed bank (germinated and ungerminated)
            hio_seed_bank(io_si) = hio_seed_bank(io_si) + &
                 (sum(litt%seed(:))+sum(litt%seed_germ(:))) * &
                 area_frac * g_per_kg * days_per_sec

            ! Sum up the input flux into the seed bank (local and external)
            hio_seeds_in(io_si) = hio_seeds_in(io_si) + &
                 (sum(litt%seed_in_local(:)) + sum(litt%seed_in_extern(:))) * &
                 area_frac * g_per_kg * days_per_sec
            
            cpatch => cpatch%younger
         end do
         

         ! ------------------------------------------------------------------------------
         ! Diagnostics discretized by element type
         ! ------------------------------------------------------------------------------

         hio_cwd_elcwd(io_si,:)   = 0._r8

         do el = 1, num_elements
            
            flux_diags => sites(s)%flux_diags(el)
            
            ! Sum up all input litter fluxes (above below, fines, cwd) [kg/ha/day]
            hio_litter_in_elem(io_si, el) =  & 
                 sum(flux_diags%cwd_ag_input(:)) + & 
                 sum(flux_diags%cwd_bg_input(:)) + &
                 sum(flux_diags%leaf_litter_input(:)) + &
                 sum(flux_diags%root_litter_input(:))

            hio_cwd_ag_elem(io_si,el)         = 0._r8
            hio_cwd_bg_elem(io_si,el)         = 0._r8
            hio_fines_ag_elem(io_si,el)       = 0._r8
            hio_fines_bg_elem(io_si,el)       = 0._r8

            hio_seed_bank_elem(io_si,el)      = 0._r8
            hio_seed_germ_elem(io_si,el)      = 0._r8
            hio_seed_decay_elem(io_si,el)     = 0._r8
            hio_seeds_in_local_elem(io_si,el) = 0._r8
            hio_seed_in_extern_elem(io_si,el) = 0._r8
            hio_litter_out_elem(io_si,el)     = 0._r8
            
            ! Plant multi-element states and fluxes
            ! Zero states, and set the fluxes
            if(element_list(el).eq.carbon12_element)then
               this%hvars(ih_totvegc_scpf)%r82d(io_si,:) = 0._r8
               this%hvars(ih_leafc_scpf)%r82d(io_si,:)   = 0._r8
               this%hvars(ih_fnrtc_scpf)%r82d(io_si,:)   = 0._r8
               this%hvars(ih_sapwc_scpf)%r82d(io_si,:)   = 0._r8
               this%hvars(ih_storec_scpf)%r82d(io_si,:)  = 0._r8
               this%hvars(ih_reproc_scpf)%r82d(io_si,:)  = 0._r8

               this%hvars(ih_cefflux_scpf)%r82d(io_si,:) = &
                    sites(s)%flux_diags(el)%nutrient_efflux_scpf(:)
               
               this%hvars(ih_cefflux)%r81d(io_si) = & 
                    sum(sites(s)%flux_diags(el)%nutrient_efflux_scpf(:),dim=1)
               
            elseif(element_list(el).eq.nitrogen_element)then
               
               this%hvars(ih_totvegn_scpf)%r82d(io_si,:) = 0._r8
               this%hvars(ih_leafn_scpf)%r82d(io_si,:)   = 0._r8
               this%hvars(ih_fnrtn_scpf)%r82d(io_si,:)   = 0._r8
               this%hvars(ih_sapwn_scpf)%r82d(io_si,:)   = 0._r8
               this%hvars(ih_storen_scpf)%r82d(io_si,:)  = 0._r8
               this%hvars(ih_repron_scpf)%r82d(io_si,:)  = 0._r8

               this%hvars(ih_nuptake_scpf)%r82d(io_si,:) = &
                    sites(s)%flux_diags(el)%nutrient_uptake_scpf(:)

               this%hvars(ih_nefflux_scpf)%r82d(io_si,:) = &
                    sites(s)%flux_diags(el)%nutrient_efflux_scpf(:)

               this%hvars(ih_nneedgrow_scpf)%r82d(io_si,:) = &
                    sites(s)%flux_diags(el)%nutrient_needgrow_scpf(:)

               this%hvars(ih_nneedmax_scpf)%r82d(io_si,:) = &
                    sites(s)%flux_diags(el)%nutrient_needmax_scpf(:)
               
               this%hvars(ih_nneedgrow)%r81d(io_si) = &
                    sum(sites(s)%flux_diags(el)%nutrient_needgrow_scpf(:),dim=1)

               this%hvars(ih_nneedmax)%r81d(io_si) = &
                    sum(sites(s)%flux_diags(el)%nutrient_needmax_scpf(:),dim=1)

               this%hvars(ih_nuptake)%r81d(io_si) = & 
                    sum(sites(s)%flux_diags(el)%nutrient_uptake_scpf(:),dim=1)

               this%hvars(ih_nefflux)%r81d(io_si) = & 
                    sum(sites(s)%flux_diags(el)%nutrient_efflux_scpf(:),dim=1)
               
               
            elseif(element_list(el).eq.phosphorus_element)then
               this%hvars(ih_totvegp_scpf)%r82d(io_si,:) = 0._r8
               this%hvars(ih_leafp_scpf)%r82d(io_si,:)   = 0._r8
               this%hvars(ih_fnrtp_scpf)%r82d(io_si,:)   = 0._r8
               this%hvars(ih_sapwp_scpf)%r82d(io_si,:)   = 0._r8
               this%hvars(ih_storep_scpf)%r82d(io_si,:)  = 0._r8
               this%hvars(ih_reprop_scpf)%r82d(io_si,:)  = 0._r8

               this%hvars(ih_puptake_scpf)%r82d(io_si,:) = &
                    sites(s)%flux_diags(el)%nutrient_uptake_scpf(:)

               this%hvars(ih_pefflux_scpf)%r82d(io_si,:) = &
                    sites(s)%flux_diags(el)%nutrient_efflux_scpf(:)
               
               this%hvars(ih_pneedgrow_scpf)%r82d(io_si,:) = &
                    sites(s)%flux_diags(el)%nutrient_needgrow_scpf(:)

               this%hvars(ih_pneedmax_scpf)%r82d(io_si,:) = &
                    sites(s)%flux_diags(el)%nutrient_needmax_scpf(:)

               this%hvars(ih_pneedgrow)%r81d(io_si) = &
                    sum(sites(s)%flux_diags(el)%nutrient_needgrow_scpf(:),dim=1)

               this%hvars(ih_pneedmax)%r81d(io_si) = &
                    sum(sites(s)%flux_diags(el)%nutrient_needmax_scpf(:),dim=1)

               this%hvars(ih_puptake)%r81d(io_si) = & 
                    sum(sites(s)%flux_diags(el)%nutrient_uptake_scpf(:),dim=1)
               
               this%hvars(ih_pefflux)%r81d(io_si) = & 
                    sum(sites(s)%flux_diags(el)%nutrient_efflux_scpf(:),dim=1)
               
            end if


            cpatch => sites(s)%oldest_patch
            do while(associated(cpatch))

               litt => cpatch%litter(el)

               area_frac = cpatch%area * AREA_INV

               ! Sum up all output fluxes (fragmentation)
               hio_litter_out_elem(io_si,el) = hio_litter_out_elem(io_si,el) + &
                    (sum(litt%leaf_fines_frag(:)) + &
                     sum(litt%root_fines_frag(:,:)) + &
                     sum(litt%ag_cwd_frag(:)) + & 
                     sum(litt%bg_cwd_frag(:,:))) * cpatch%area

               hio_seed_bank_elem(io_si,el) = hio_seed_bank_elem(io_si,el) + & 
                    sum(litt%seed(:)) * cpatch%area

               hio_seed_germ_elem(io_si,el) = hio_seed_germ_elem(io_si,el) + &
                    sum(litt%seed_germ(:)) *  cpatch%area
                    
               hio_seed_decay_elem(io_si,el) = hio_seed_decay_elem(io_si,el) + & 
                    sum(litt%seed_decay(:)) * cpatch%area

               hio_seeds_in_local_elem(io_si,el) = hio_seeds_in_local_elem(io_si,el) + & 
                    sum(litt%seed_in_local(:)) *  cpatch%area

               hio_seed_in_extern_elem(io_si,el) = hio_seed_in_extern_elem(io_si,el) + & 
                    sum(litt%seed_in_extern(:)) * cpatch%area

               ! Litter State Variables
               hio_cwd_ag_elem(io_si,el) = hio_cwd_ag_elem(io_si,el) + &
                     sum(litt%ag_cwd(:)) * cpatch%area
               
               hio_cwd_bg_elem(io_si,el) = hio_cwd_bg_elem(io_si,el) + &
                     sum(litt%bg_cwd(:,:)) * cpatch%area
               
               hio_fines_ag_elem(io_si,el) = hio_fines_ag_elem(io_si,el) + & 
                     sum(litt%leaf_fines(:)) * cpatch%area
               
               hio_fines_bg_elem(io_si,el) = hio_fines_bg_elem(io_si,el) + &
                     sum(litt%root_fines(:,:)) * cpatch%area

               do cwd=1,ncwd
                   elcwd = (el-1)*ncwd+cwd
                   hio_cwd_elcwd(io_si,elcwd) = hio_cwd_elcwd(io_si,elcwd) + & 
                         (litt%ag_cwd(cwd) + sum(litt%bg_cwd(cwd,:))) * cpatch%area

               end do

               ! Load Mass States
               ccohort => cpatch%tallest
               do while(associated(ccohort))

                  sapw_m   = ccohort%prt%GetState(sapw_organ, element_list(el))
                  struct_m = ccohort%prt%GetState(struct_organ, element_list(el))
                  leaf_m   = ccohort%prt%GetState(leaf_organ, element_list(el))
                  fnrt_m   = ccohort%prt%GetState(fnrt_organ, element_list(el))
                  store_m  = ccohort%prt%GetState(store_organ, element_list(el))
                  repro_m  = ccohort%prt%GetState(repro_organ, element_list(el))
                  total_m  = sapw_m+struct_m+leaf_m+fnrt_m+store_m+repro_m
                  
                  i_scpf = ccohort%size_by_pft_class

                  if(element_list(el).eq.carbon12_element)then
                     this%hvars(ih_totvegc_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_totvegc_scpf)%r82d(io_si,i_scpf) + total_m * ccohort%n
                     this%hvars(ih_leafc_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_leafc_scpf)%r82d(io_si,i_scpf) + leaf_m * ccohort%n
                     this%hvars(ih_fnrtc_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_fnrtc_scpf)%r82d(io_si,i_scpf) + fnrt_m * ccohort%n
                     this%hvars(ih_sapwc_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_sapwc_scpf)%r82d(io_si,i_scpf) + sapw_m * ccohort%n
                     this%hvars(ih_storec_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_storec_scpf)%r82d(io_si,i_scpf) + store_m * ccohort%n
                     this%hvars(ih_reproc_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_reproc_scpf)%r82d(io_si,i_scpf) + repro_m * ccohort%n
                  elseif(element_list(el).eq.nitrogen_element)then
                     this%hvars(ih_totvegn_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_totvegn_scpf)%r82d(io_si,i_scpf) + total_m * ccohort%n
                     this%hvars(ih_leafn_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_leafn_scpf)%r82d(io_si,i_scpf) + leaf_m * ccohort%n
                     this%hvars(ih_fnrtn_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_fnrtn_scpf)%r82d(io_si,i_scpf) + fnrt_m * ccohort%n
                     this%hvars(ih_sapwn_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_sapwn_scpf)%r82d(io_si,i_scpf) + sapw_m * ccohort%n
                     this%hvars(ih_storen_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_storen_scpf)%r82d(io_si,i_scpf) + store_m * ccohort%n
                     this%hvars(ih_repron_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_repron_scpf)%r82d(io_si,i_scpf) + repro_m * ccohort%n
                  elseif(element_list(el).eq.phosphorus_element)then
                     this%hvars(ih_totvegp_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_totvegp_scpf)%r82d(io_si,i_scpf) + total_m * ccohort%n
                     this%hvars(ih_leafp_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_leafp_scpf)%r82d(io_si,i_scpf) + leaf_m * ccohort%n
                     this%hvars(ih_fnrtp_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_fnrtp_scpf)%r82d(io_si,i_scpf) + fnrt_m * ccohort%n
                     this%hvars(ih_sapwp_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_sapwp_scpf)%r82d(io_si,i_scpf) + sapw_m * ccohort%n
                     this%hvars(ih_storep_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_storep_scpf)%r82d(io_si,i_scpf) + store_m * ccohort%n
                     this%hvars(ih_reprop_scpf)%r82d(io_si,i_scpf) = & 
                          this%hvars(ih_reprop_scpf)%r82d(io_si,i_scpf) + repro_m * ccohort%n
                  end if
                  
                  ccohort => ccohort%shorter
               end do
                    
               cpatch => cpatch%younger
            end do

         end do
         



         
         ! pass demotion rates and associated carbon fluxes to history
         do i_sz = 1,nlevsclass
            hio_demotion_rate_sz(io_si,i_sz) = sites(s)%demotion_rate(i_sz) * days_per_year
            hio_promotion_rate_sz(io_si,i_sz) = sites(s)%promotion_rate(i_sz) * days_per_year
         end do
         !
         ! convert kg C / ha / day to gc / m2 / sec
         hio_demotion_carbonflux(io_si) = sites(s)%demotion_carbonflux * g_per_kg * ha_per_m2 * days_per_sec
         hio_promotion_carbonflux(io_si) = sites(s)%promotion_carbonflux * g_per_kg * ha_per_m2 * days_per_sec
         !
         ! mortality-associated carbon fluxes
         
         hio_canopy_mortality_carbonflux(io_si) = hio_canopy_mortality_carbonflux(io_si) + &
               sites(s)%term_carbonflux_canopy * g_per_kg * days_per_sec * ha_per_m2
         
         hio_understory_mortality_carbonflux(io_si) = hio_understory_mortality_carbonflux(io_si) + &
               sites(s)%term_carbonflux_ustory * g_per_kg * days_per_sec * ha_per_m2

         ! and zero the site-level termination carbon flux variable
         sites(s)%term_carbonflux_canopy = 0._r8
         sites(s)%term_carbonflux_ustory = 0._r8
         !

         ! add the site-level disturbance-associated cwd and litter input fluxes to thir respective flux fields

         do i_cwd = 1, ncwd
             hio_cwd_ag_in_cwdsc(io_si, i_cwd) = hio_cwd_ag_in_cwdsc(io_si, i_cwd) + &
                   flux_diags_c%cwd_ag_input(i_cwd) * g_per_kg
             
             hio_cwd_bg_in_cwdsc(io_si, i_cwd) = hio_cwd_bg_in_cwdsc(io_si, i_cwd) + &
                   flux_diags_c%cwd_bg_input(i_cwd) * g_per_kg

         end do

         ! and reset the disturbance-related field buffers

         do el = 1, num_elements
             call sites(s)%flux_diags(el)%ZeroFluxDiags()
         end do

      enddo ! site loop
      
    end associate

    return
  end subroutine update_history_dyn
 
  subroutine update_history_hifrq(this,nc,nsites,sites,bc_in,dt_tstep)

    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after rapid timescale productivity calculations (gpp and respiration).
    ! ---------------------------------------------------------------------------------
    
    use EDTypesMod          , only : nclmax, nlevleaf
    !
    ! Arguments
    class(fates_history_interface_type)                 :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(nsites)
    real(r8)                , intent(in)            :: dt_tstep
    
    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ipa      ! The local "I"ndex of "PA"tches 
    integer  :: lb1,ub1,lb2,ub2  ! IO array bounds for the calling thread
    integer  :: ivar             ! index of IO variable object vector
    integer  :: ft               ! functional type index
    real(r8) :: n_density   ! individual of cohort per m2.
    real(r8) :: resp_g      ! growth respiration per timestep [kgC/indiv/step]
    real(r8) :: npp         ! npp for this time-step (adjusted for g resp) [kgC/indiv/step]
    real(r8) :: aresp       ! autotrophic respiration (adjusted for g resp) [kgC/indiv/step]
    real(r8) :: n_perm2     ! individuals per m2 for the whole column
    real(r8) :: patch_area_by_age(nlevage) ! patch area in each bin for normalizing purposes
    real(r8) :: canopy_area_by_age(nlevage) ! canopy area in each bin for normalizing purposes
    real(r8), parameter :: tiny = 1.e-5_r8      ! some small number
    integer  :: ipa2     ! patch incrementer
    integer :: cnlfpft_indx, cnlf_indx, ipft, ican, ileaf ! more iterators and indices
    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort
    real(r8) :: per_dt_tstep          ! Time step in frequency units (/s)

    associate( hio_gpp         => this%hvars(ih_gpp)%r81d, &
               hio_npp         => this%hvars(ih_npp)%r81d, &
               hio_aresp       => this%hvars(ih_aresp)%r81d, &
               hio_maint_resp  => this%hvars(ih_maint_resp)%r81d, &
               hio_growth_resp => this%hvars(ih_growth_resp)%r81d, &
               hio_c_stomata   => this%hvars(ih_c_stomata)%r81d, &
               hio_c_lblayer   => this%hvars(ih_c_lblayer)%r81d, &
               hio_nep         => this%hvars(ih_nep)%r81d, & 
               hio_ar_scpf     => this%hvars(ih_ar_scpf)%r82d, &
               hio_ar_grow_scpf   => this%hvars(ih_ar_grow_scpf)%r82d, &
               hio_ar_maint_scpf  => this%hvars(ih_ar_maint_scpf)%r82d, &
               hio_ar_agsapm_scpf => this%hvars(ih_ar_agsapm_scpf)%r82d, &
               hio_ar_darkm_scpf  => this%hvars(ih_ar_darkm_scpf)%r82d, &
               hio_ar_crootm_scpf => this%hvars(ih_ar_crootm_scpf)%r82d, &
               hio_ar_frootm_scpf => this%hvars(ih_ar_frootm_scpf)%r82d, &
               hio_gpp_canopy     => this%hvars(ih_gpp_canopy)%r81d, &
               hio_ar_canopy      => this%hvars(ih_ar_canopy)%r81d, &
               hio_gpp_understory => this%hvars(ih_gpp_understory)%r81d, &
               hio_ar_understory  => this%hvars(ih_ar_understory)%r81d, &
               hio_rdark_canopy_sz             => this%hvars(ih_rdark_canopy_sz)%r82d, &
               hio_livestem_mr_canopy_sz       => this%hvars(ih_livestem_mr_canopy_sz)%r82d, &
               hio_livecroot_mr_canopy_sz      => this%hvars(ih_livecroot_mr_canopy_sz)%r82d, &
               hio_froot_mr_canopy_sz          => this%hvars(ih_froot_mr_canopy_sz)%r82d, &
               hio_resp_g_canopy_sz            => this%hvars(ih_resp_g_canopy_sz)%r82d, &
               hio_resp_m_canopy_sz            => this%hvars(ih_resp_m_canopy_sz)%r82d, &
               hio_rdark_understory_sz         => this%hvars(ih_rdark_understory_sz)%r82d, &
               hio_livestem_mr_understory_sz   => this%hvars(ih_livestem_mr_understory_sz)%r82d, &
               hio_livecroot_mr_understory_sz  => this%hvars(ih_livecroot_mr_understory_sz)%r82d, &
               hio_froot_mr_understory_sz      => this%hvars(ih_froot_mr_understory_sz)%r82d, &
               hio_resp_g_understory_sz        => this%hvars(ih_resp_g_understory_sz)%r82d, &
               hio_resp_m_understory_sz        => this%hvars(ih_resp_m_understory_sz)%r82d, &
               hio_leaf_mr         => this%hvars(ih_leaf_mr)%r81d, &
               hio_froot_mr        => this%hvars(ih_froot_mr)%r81d, &
               hio_livecroot_mr    => this%hvars(ih_livecroot_mr)%r81d, &
               hio_livestem_mr     => this%hvars(ih_livestem_mr)%r81d, &
               hio_gpp_age         => this%hvars(ih_gpp_age)%r82d, &
               hio_npp_age         => this%hvars(ih_npp_age)%r82d, &
               hio_c_stomata_age   => this%hvars(ih_c_stomata_age)%r82d, &
               hio_c_lblayer_age   => this%hvars(ih_c_lblayer_age)%r82d, &
               hio_parsun_z_cnlf     => this%hvars(ih_parsun_z_cnlf)%r82d, &
               hio_parsha_z_cnlf     => this%hvars(ih_parsha_z_cnlf)%r82d, &
               hio_ts_net_uptake_cnlf => this%hvars(ih_ts_net_uptake_cnlf)%r82d, &
               hio_parsun_z_cnlfpft  => this%hvars(ih_parsun_z_cnlfpft)%r82d, &
               hio_parsha_z_cnlfpft  => this%hvars(ih_parsha_z_cnlfpft)%r82d, &
               hio_laisun_z_cnlf     => this%hvars(ih_laisun_z_cnlf)%r82d, &
               hio_laisha_z_cnlf     => this%hvars(ih_laisha_z_cnlf)%r82d, &
               hio_laisun_z_cnlfpft  => this%hvars(ih_laisun_z_cnlfpft)%r82d, &
               hio_laisha_z_cnlfpft  => this%hvars(ih_laisha_z_cnlfpft)%r82d, &
               hio_laisun_top_can     => this%hvars(ih_laisun_top_can)%r82d, &
               hio_laisha_top_can     => this%hvars(ih_laisha_top_can)%r82d, &
               hio_fabd_sun_cnlfpft  => this%hvars(ih_fabd_sun_cnlfpft)%r82d, &
               hio_fabd_sha_cnlfpft  => this%hvars(ih_fabd_sha_cnlfpft)%r82d, &
               hio_fabi_sun_cnlfpft  => this%hvars(ih_fabi_sun_cnlfpft)%r82d, &
               hio_fabi_sha_cnlfpft  => this%hvars(ih_fabi_sha_cnlfpft)%r82d, &
               hio_fabd_sun_cnlf  => this%hvars(ih_fabd_sun_cnlf)%r82d, &
               hio_fabd_sha_cnlf  => this%hvars(ih_fabd_sha_cnlf)%r82d, &
               hio_fabi_sun_cnlf  => this%hvars(ih_fabi_sun_cnlf)%r82d, &
               hio_fabi_sha_cnlf  => this%hvars(ih_fabi_sha_cnlf)%r82d, &
               hio_parprof_dir_cnlf  => this%hvars(ih_parprof_dir_cnlf)%r82d, &
               hio_parprof_dif_cnlf  => this%hvars(ih_parprof_dif_cnlf)%r82d, &
               hio_parprof_dir_cnlfpft  => this%hvars(ih_parprof_dir_cnlfpft)%r82d, &
               hio_parprof_dif_cnlfpft  => this%hvars(ih_parprof_dif_cnlfpft)%r82d, &
               hio_fabd_sun_top_can  => this%hvars(ih_fabd_sun_top_can)%r82d, &
               hio_fabd_sha_top_can  => this%hvars(ih_fabd_sha_top_can)%r82d, &
               hio_fabi_sun_top_can  => this%hvars(ih_fabi_sun_top_can)%r82d, &
               hio_fabi_sha_top_can  => this%hvars(ih_fabi_sha_top_can)%r82d, &
               hio_parsun_top_can     => this%hvars(ih_parsun_top_can)%r82d, &
               hio_parsha_top_can     => this%hvars(ih_parsha_top_can)%r82d &
               )


      ! Flush the relevant history variables 
      call this%flush_hvars(nc,upfreq_in=2)

      per_dt_tstep = 1.0_r8/dt_tstep

      do s = 1,nsites
         
         io_si  = this%iovar_map(nc)%site_index(s)
         hio_nep(io_si) = -bc_in(s)%tot_het_resp ! (gC/m2/s)
         
         ipa = 0
         cpatch => sites(s)%oldest_patch

         patch_area_by_age(1:nlevage) = 0._r8
         canopy_area_by_age(1:nlevage) = 0._r8

         do while(associated(cpatch))
            
            patch_area_by_age(cpatch%age_class)  = &
                 patch_area_by_age(cpatch%age_class) + cpatch%area

            canopy_area_by_age(cpatch%age_class) = &
                 canopy_area_by_age(cpatch%age_class) + cpatch%total_canopy_area

            ! Canopy resitance terms
            hio_c_stomata_age(io_si,cpatch%age_class) = &
                 hio_c_stomata_age(io_si,cpatch%age_class) + &
                 cpatch%c_stomata * cpatch%total_canopy_area
            
            hio_c_lblayer_age(io_si,cpatch%age_class) = &
                 hio_c_lblayer_age(io_si,cpatch%age_class) + &
                 cpatch%c_lblayer * cpatch%total_canopy_area
            
            hio_c_stomata(io_si) = hio_c_stomata(io_si) + &
                 cpatch%c_stomata * cpatch%total_canopy_area
            
            hio_c_lblayer(io_si) = hio_c_lblayer(io_si) + &
                 cpatch%c_lblayer * cpatch%total_canopy_area

            ccohort => cpatch%shortest
            do while(associated(ccohort))
               
               n_perm2   = ccohort%n * AREA_INV
               
               if ( .not. ccohort%isnew ) then

                   npp    = ccohort%npp_tstep
                   resp_g = ccohort%resp_g_tstep
                   aresp  = ccohort%resp_tstep
                                    
                  ! Calculate index for the scpf class
                  associate( scpf => ccohort%size_by_pft_class, &
                             sz => ccohort%size_class )
                    
                  ! scale up cohort fluxes to the site level
                  hio_npp(io_si) = hio_npp(io_si) + &
                        npp * g_per_kg * n_perm2 * per_dt_tstep
                  
                  hio_gpp(io_si) = hio_gpp(io_si) + &
                        ccohort%gpp_tstep * g_per_kg * n_perm2 * per_dt_tstep
                  hio_aresp(io_si) = hio_aresp(io_si) + &
                        aresp * g_per_kg * n_perm2 * per_dt_tstep
                  hio_growth_resp(io_si) = hio_growth_resp(io_si) + &
                        resp_g * g_per_kg * n_perm2 * per_dt_tstep
                  hio_maint_resp(io_si) = hio_maint_resp(io_si) + &
                        ccohort%resp_m * g_per_kg * n_perm2 * per_dt_tstep

                  ! Add up the total Net Ecosystem Production
                  ! for this timestep.  [gC/m2/s]
                  hio_nep(io_si) = hio_nep(io_si) + &
                       npp * g_per_kg * n_perm2 * per_dt_tstep

                  ! aggregate MR fluxes to the site level
                  hio_leaf_mr(io_si) = hio_leaf_mr(io_si) + ccohort%rdark &
                       * n_perm2 *  sec_per_day * days_per_year
                  hio_froot_mr(io_si) = hio_froot_mr(io_si) + ccohort%froot_mr &
                       * n_perm2 *  sec_per_day * days_per_year
                  hio_livecroot_mr(io_si) = hio_livecroot_mr(io_si) + ccohort%livecroot_mr &
                       * n_perm2 *  sec_per_day * days_per_year
                  hio_livestem_mr(io_si) = hio_livestem_mr(io_si) + ccohort%livestem_mr &
                       * n_perm2 *  sec_per_day * days_per_year

                  ! Total AR (kgC/m2/yr) = (kgC/plant/step) / (s/step) * (plant/m2) * (s/yr)
                  hio_ar_scpf(io_si,scpf)    =   hio_ar_scpf(io_si,scpf) + &
                        (ccohort%resp_tstep/dt_tstep) * n_perm2 * sec_per_day * days_per_year

                  ! Growth AR (kgC/m2/yr)
                  hio_ar_grow_scpf(io_si,scpf) = hio_ar_grow_scpf(io_si,scpf) + &
                        (resp_g/dt_tstep) * n_perm2 * sec_per_day * days_per_year

                  ! Maint AR (kgC/m2/yr)
                  hio_ar_maint_scpf(io_si,scpf) = hio_ar_maint_scpf(io_si,scpf) + &
                        (ccohort%resp_m/dt_tstep) * n_perm2 * sec_per_day * days_per_year
                  
                  ! Maintenance AR partition variables are stored as rates (kgC/plant/s)
                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_agsapm_scpf(io_si,scpf) = hio_ar_agsapm_scpf(io_si,scpf) + &
                        ccohort%livestem_mr * n_perm2 * sec_per_day * days_per_year

                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_darkm_scpf(io_si,scpf) = hio_ar_darkm_scpf(io_si,scpf) + &
                        ccohort%rdark * n_perm2 *  sec_per_day * days_per_year

                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_crootm_scpf(io_si,scpf) = hio_ar_crootm_scpf(io_si,scpf) + &
                        ccohort%livecroot_mr * n_perm2 * sec_per_day * days_per_year

                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_frootm_scpf(io_si,scpf) = hio_ar_frootm_scpf(io_si,scpf) + &
                        ccohort%froot_mr * n_perm2  * sec_per_day * days_per_year


                  ! accumulate fluxes per patch age bin
                  hio_gpp_age(io_si,cpatch%age_class) = hio_gpp_age(io_si,cpatch%age_class) &
                       + ccohort%gpp_tstep * ccohort%n * g_per_kg * per_dt_tstep
                  hio_npp_age(io_si,cpatch%age_class) = hio_npp_age(io_si,cpatch%age_class) &
                       + npp * ccohort%n * g_per_kg * per_dt_tstep

                  ! accumulate fluxes on canopy- and understory- separated fluxes
                  if (ccohort%canopy_layer .eq. 1) then
                     !
                     ! bulk fluxes are in gC / m2 / s
                     hio_gpp_canopy(io_si) = hio_gpp_canopy(io_si) + &
                          ccohort%gpp_tstep * g_per_kg * n_perm2 * per_dt_tstep                     
                     hio_ar_canopy(io_si) = hio_ar_canopy(io_si) + &
                          aresp * g_per_kg * n_perm2 * per_dt_tstep                     

                     !
                     ! size-resolved respiration fluxes are in kg C / ha / yr
                     hio_rdark_canopy_sz(io_si,sz) = hio_rdark_canopy_sz(io_si,sz) + &
                          ccohort%rdark  * ccohort%n * sec_per_day * days_per_year
                     hio_livestem_mr_canopy_sz(io_si,sz) = hio_livestem_mr_canopy_sz(io_si,sz) + &
                          ccohort%livestem_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_livecroot_mr_canopy_sz(io_si,sz) = hio_livecroot_mr_canopy_sz(io_si,sz) + &
                          ccohort%livecroot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_froot_mr_canopy_sz(io_si,sz) = hio_froot_mr_canopy_sz(io_si,sz) + &
                          ccohort%froot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_resp_g_canopy_sz(io_si,sz) = hio_resp_g_canopy_sz(io_si,sz) + &
                          resp_g  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                     hio_resp_m_canopy_sz(io_si,sz) = hio_resp_m_canopy_sz(io_si,sz) + &
                          ccohort%resp_m  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                  else
                     !
                     ! bulk fluxes are in gC / m2 / s
                     hio_gpp_understory(io_si) = hio_gpp_understory(io_si) + &
                          ccohort%gpp_tstep * g_per_kg * n_perm2 * per_dt_tstep                     
                     hio_ar_understory(io_si) = hio_ar_understory(io_si) + &
                          aresp * g_per_kg * n_perm2 * per_dt_tstep                     

                     !
                     ! size-resolved respiration fluxes are in kg C / ha / yr
                     hio_rdark_understory_sz(io_si,sz) = hio_rdark_understory_sz(io_si,sz) + &
                          ccohort%rdark  * ccohort%n * sec_per_day * days_per_year
                     hio_livestem_mr_understory_sz(io_si,sz) = hio_livestem_mr_understory_sz(io_si,sz) + &
                          ccohort%livestem_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_livecroot_mr_understory_sz(io_si,sz) = hio_livecroot_mr_understory_sz(io_si,sz) + &
                          ccohort%livecroot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_froot_mr_understory_sz(io_si,sz) = hio_froot_mr_understory_sz(io_si,sz) + &
                          ccohort%froot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_resp_g_understory_sz(io_si,sz) = hio_resp_g_understory_sz(io_si,sz) + &
                          resp_g  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                     hio_resp_m_understory_sz(io_si,sz) = hio_resp_m_understory_sz(io_si,sz) + &
                          ccohort%resp_m  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                  endif
                end associate
               endif

               !!! canopy leaf carbon balance
               ican = ccohort%canopy_layer
               do ileaf=1,ccohort%nv
                  cnlf_indx = ileaf + (ican-1) * nlevleaf
                  hio_ts_net_uptake_cnlf(io_si, cnlf_indx) = hio_ts_net_uptake_cnlf(io_si, cnlf_indx) + &
                       ccohort%ts_net_uptake(ileaf) * g_per_kg * per_dt_tstep * ccohort%c_area / AREA
               end do

               ccohort => ccohort%taller
            enddo ! cohort loop

            ! summarize radiation profiles through the canopy
            do ipft=1,numpft
               do ican=1,nclmax         !  cpatch%ncl_p  ?
                  do ileaf=1,nlevleaf   !  cpatch%ncan(ican,ipft) ?
                     ! calculate where we are on multiplexed dimensions
                     cnlfpft_indx = ileaf + (ican-1) * nlevleaf + (ipft-1) * nlevleaf * nclmax 
                     cnlf_indx = ileaf + (ican-1) * nlevleaf
                     !
                     ! first do all the canopy x leaf x pft calculations
                     hio_parsun_z_cnlfpft(io_si,cnlfpft_indx) = hio_parsun_z_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_parsun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_parsha_z_cnlfpft(io_si,cnlfpft_indx) = hio_parsha_z_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_parsha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_laisun_z_cnlfpft(io_si,cnlfpft_indx) = hio_laisun_z_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_laisun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_laisha_z_cnlfpft(io_si,cnlfpft_indx) = hio_laisha_z_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_laisha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_fabd_sun_cnlfpft(io_si,cnlfpft_indx) = hio_fabd_sun_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabd_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabd_sha_cnlfpft(io_si,cnlfpft_indx) = hio_fabd_sha_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabd_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sun_cnlfpft(io_si,cnlfpft_indx) = hio_fabi_sun_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabi_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sha_cnlfpft(io_si,cnlfpft_indx) = hio_fabi_sha_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabi_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_parprof_dir_cnlfpft(io_si,cnlfpft_indx) = hio_parprof_dir_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%parprof_pft_dir_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_parprof_dif_cnlfpft(io_si,cnlfpft_indx) = hio_parprof_dif_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%parprof_pft_dif_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     ! summarize across all PFTs
                     hio_parsun_z_cnlf(io_si,cnlf_indx) = hio_parsun_z_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_parsun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_parsha_z_cnlf(io_si,cnlf_indx) = hio_parsha_z_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_parsha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_laisun_z_cnlf(io_si,cnlf_indx) = hio_laisun_z_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_laisun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_laisha_z_cnlf(io_si,cnlf_indx) = hio_laisha_z_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_laisha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_fabd_sun_cnlf(io_si,cnlf_indx) = hio_fabd_sun_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabd_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabd_sha_cnlf(io_si,cnlf_indx) = hio_fabd_sha_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabd_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sun_cnlf(io_si,cnlf_indx) = hio_fabi_sun_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabi_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sha_cnlf(io_si,cnlf_indx) = hio_fabi_sha_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabi_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV

                  end do
                  !
                  ! summarize just the top leaf level across all PFTs, for each canopy level
                  hio_parsun_top_can(io_si,ican) = hio_parsun_top_can(io_si,ican) + &
                       cpatch%ed_parsun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_parsha_top_can(io_si,ican) = hio_parsha_top_can(io_si,ican) + &
                       cpatch%ed_parsha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  !
                  hio_laisun_top_can(io_si,ican) = hio_laisun_top_can(io_si,ican) + &
                       cpatch%ed_laisun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_laisha_top_can(io_si,ican) = hio_laisha_top_can(io_si,ican) + &
                       cpatch%ed_laisha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  !
                  hio_fabd_sun_top_can(io_si,ican) = hio_fabd_sun_top_can(io_si,ican) + &
                       cpatch%fabd_sun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_fabd_sha_top_can(io_si,ican) = hio_fabd_sha_top_can(io_si,ican) + &
                       cpatch%fabd_sha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_fabi_sun_top_can(io_si,ican) = hio_fabi_sun_top_can(io_si,ican) + &
                       cpatch%fabi_sun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_fabi_sha_top_can(io_si,ican) = hio_fabi_sha_top_can(io_si,ican) + &
                       cpatch%fabi_sha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  !
               end do
           end do

            ! PFT-mean radiation profiles
            do ican=1,nclmax
               do ileaf=1,nlevleaf
                  ! calculate where we are on multiplexed dimensions
                  cnlf_indx = ileaf + (ican-1) * nlevleaf
                  !
                  hio_parprof_dir_cnlf(io_si,cnlf_indx) = hio_parprof_dir_cnlf(io_si,cnlf_indx) + &
                       cpatch%parprof_dir_z(ican,ileaf) * cpatch%area * AREA_INV
                  hio_parprof_dif_cnlf(io_si,cnlf_indx) = hio_parprof_dif_cnlf(io_si,cnlf_indx) + &
                       cpatch%parprof_dif_z(ican,ileaf) * cpatch%area * AREA_INV
               end do
            end do

            ipa = ipa + 1
            cpatch => cpatch%younger
         end do !patch loop

         do ipa2 = 1, nlevage
            if (patch_area_by_age(ipa2) .gt. tiny) then
               hio_gpp_age(io_si, ipa2) = hio_gpp_age(io_si, ipa2) / (patch_area_by_age(ipa2))
               hio_npp_age(io_si, ipa2) = hio_npp_age(io_si, ipa2) / (patch_area_by_age(ipa2))
            else
               hio_gpp_age(io_si, ipa2) = 0._r8
               hio_npp_age(io_si, ipa2) = 0._r8
            endif

            ! Normalize resistance diagnostics
            if (canopy_area_by_age(ipa2) .gt. tiny) then
               hio_c_stomata_age(io_si,ipa2) = &
                    hio_c_stomata_age(io_si,ipa2) / canopy_area_by_age(ipa2)

               hio_c_lblayer_age(io_si,ipa2) = &
                    hio_c_lblayer_age(io_si,ipa2) / canopy_area_by_age(ipa2)
            else
               hio_c_stomata_age(io_si,ipa2) = 0._r8
               hio_c_lblayer_age(io_si,ipa2) = 0._r8
            end if
            
         end do
         
         ! Normalize resistance diagnostics
         if ( sum(canopy_area_by_age(1:nlevage)) .gt. tiny) then
            hio_c_stomata(io_si) = hio_c_stomata(io_si) / sum(canopy_area_by_age(1:nlevage))
            hio_c_lblayer(io_si) = hio_c_lblayer(io_si) / sum(canopy_area_by_age(1:nlevage))
         else
            hio_c_stomata(io_si) = 0._r8
            hio_c_lblayer(io_si) = 0._r8
         end if

     enddo ! site loop

   end associate
 
end subroutine update_history_hifrq

  ! =====================================================================================

  subroutine update_history_hydraulics(this,nc,nsites,sites,bc_in,dt_tstep)

    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after rapid timescale productivity calculations (gpp and respiration).
    ! ---------------------------------------------------------------------------------
    
    use FatesHydraulicsMemMod, only : ed_cohort_hydr_type, nshell
    use FatesHydraulicsMemMod, only : ed_site_hydr_type
    use EDTypesMod           , only : maxpft

    
    ! Arguments
    class(fates_history_interface_type)             :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(nsites)
    real(r8)                , intent(in)            :: dt_tstep
    
    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ipa      ! The local "I"ndex of "PA"tches 
    integer  :: ft               ! functional type index
!    integer  :: io_shsl  ! The combined "SH"ell "S"oil "L"ayer index in the IO array
    real(r8), parameter :: tiny = 1.e-5_r8      ! some small number
    real(r8) :: ncohort_scpf(nlevsclass*maxpft)  ! Bins to count up cohorts counts used in weighting
    ! should be "hio_nplant_scpf"
    real(r8) :: nplant_scpf(nlevsclass*maxpft)  ! Bins to count up cohorts counts used in weighting
                                                   ! should be "hio_nplant_scpf"
    real(r8) :: number_fraction
    real(r8) :: number_fraction_rate
    real(r8) :: mean_aroot
    integer  :: ipa2     ! patch incrementer
    integer  :: ix       ! histogram x (count) bin index
    integer  :: iscpf    ! index of the scpf group
    integer  :: ipft     ! index of the pft loop
    integer  :: isz    ! index of the size-class loop
    integer  :: k        ! rhizosphere shell index
    integer  :: jsoil    ! soil layer index
    integer  :: jrhiz    ! rhizosphere layer index
    integer  :: jr1, jr2 ! Rhizosphere top and bottom layers
    integer  :: nlevrhiz ! number of rhizosphere layers
    real(r8) :: mean_soil_vwc    ! mean soil volumetric water content [m3/m3]
    real(r8) :: mean_soil_vwcsat ! mean soil saturated volumetric water content [m3/m3]
    real(r8) :: mean_soil_matpot ! mean soil water potential [MPa]
    real(r8) :: layer_areaweight ! root area weighting factor for each soil layer
    real(r8) :: areaweight       ! root area weighting factor for column
    real(r8) :: vwc              ! volumetric water content of layer [m3/m3] = theta
    real(r8) :: vwc_sat          ! saturated water content of layer [m3/m3]
    real(r8) :: psi              ! matric potential of soil layer
    character(2) :: fmt_char
    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr
    type(ed_site_hydr_type), pointer :: site_hydr

    real(r8), parameter :: daysecs = 86400.0_r8 ! What modeler doesn't recognize 86400?
    real(r8), parameter :: yeardays = 365.0_r8  ! Should this be 365.25?

    integer,  parameter :: iterh2_nhist = 50
    real(r8), parameter :: iterh2_dx    = 1._r8
    real(r8)            ::  iterh2_histx(iterh2_nhist)
    real(r8)            ::  iterh2_histy(iterh2_nhist)
    
    logical, parameter :: print_iterations = .false.
    
    if(hlm_use_planthydro.eq.ifalse) return

    
    associate( hio_errh2o_scpf  => this%hvars(ih_errh2o_scpf)%r82d, &
          hio_tran_scpf         => this%hvars(ih_tran_scpf)%r82d, &
          hio_sapflow_scpf      => this%hvars(ih_sapflow_scpf)%r82d, &
          hio_sapflow        => this%hvars(ih_sapflow)%r81d, & 
          hio_iterh1_scpf       => this%hvars(ih_iterh1_scpf)%r82d, &          
          hio_iterh2_scpf       => this%hvars(ih_iterh2_scpf)%r82d, &           
          hio_ath_scpf          => this%hvars(ih_ath_scpf)%r82d, &               
          hio_tth_scpf          => this%hvars(ih_tth_scpf)%r82d, &               
          hio_sth_scpf          => this%hvars(ih_sth_scpf)%r82d, &                     
          hio_lth_scpf          => this%hvars(ih_lth_scpf)%r82d, &                     
          hio_awp_scpf          => this%hvars(ih_awp_scpf)%r82d, &                     
          hio_twp_scpf          => this%hvars(ih_twp_scpf)%r82d, &  
          hio_swp_scpf          => this%hvars(ih_swp_scpf)%r82d, &                     
          hio_lwp_scpf          => this%hvars(ih_lwp_scpf)%r82d, &  
          hio_aflc_scpf          => this%hvars(ih_aflc_scpf)%r82d, &                     
          hio_tflc_scpf          => this%hvars(ih_tflc_scpf)%r82d, &  
          hio_sflc_scpf          => this%hvars(ih_sflc_scpf)%r82d, &                     
          hio_lflc_scpf          => this%hvars(ih_lflc_scpf)%r82d, &                   
          hio_btran_scpf        => this%hvars(ih_btran_scpf)%r82d, &
          hio_h2oveg         => this%hvars(ih_h2oveg)%r81d, &
          hio_nplant_scpf    => this%hvars(ih_nplant_scpf)%r82d, &
          hio_nplant_capf    => this%hvars(ih_nplant_capf)%r82d, &
          hio_h2oveg_hydro_err   => this%hvars(ih_h2oveg_hydro_err)%r81d, &
          hio_rootwgt_soilvwc    => this%hvars(ih_rootwgt_soilvwc)%r81d, &
          hio_rootwgt_soilvwcsat => this%hvars(ih_rootwgt_soilvwcsat)%r81d, & 
          hio_rootwgt_soilmatpot => this%hvars(ih_rootwgt_soilmatpot)%r81d, &
          hio_soilmatpot_sl         => this%hvars(ih_soilmatpot_sl)%r82d, &
          hio_soilvwc_sl            => this%hvars(ih_soilvwc_sl)%r82d, &
          hio_soilvwcsat_sl         => this%hvars(ih_soilvwcsat_sl)%r82d, &
          hio_rootuptake         => this%hvars(ih_rootuptake)%r81d, &
          hio_rootuptake_sl         => this%hvars(ih_rootuptake_sl)%r82d, &
          hio_rootuptake0_scpf      => this%hvars(ih_rootuptake0_scpf)%r82d, &
          hio_rootuptake10_scpf     => this%hvars(ih_rootuptake10_scpf)%r82d, &
          hio_rootuptake50_scpf     => this%hvars(ih_rootuptake50_scpf)%r82d, &
          hio_rootuptake100_scpf    => this%hvars(ih_rootuptake100_scpf)%r82d )

      ! Flush the relevant history variables 
      call this%flush_hvars(nc,upfreq_in=4)

      if(print_iterations) then
          do iscpf = 1,iterh2_nhist
              iterh2_histx(iscpf) = iterh2_dx*real(iscpf-1,r8)
          end do
      end if
      do s = 1,nsites

         site_hydr => sites(s)%si_hydr
         nlevrhiz = site_hydr%nlevrhiz
         jr1 = site_hydr%i_rhiz_t
         jr2 = site_hydr%i_rhiz_b

         io_si  = this%iovar_map(nc)%site_index(s)
         
         hio_h2oveg(io_si)              = site_hydr%h2oveg
         hio_h2oveg_hydro_err(io_si)    = site_hydr%h2oveg_hydro_err

        
         
         ! Get column means of some soil diagnostics, these are weighted
         ! by the amount of fine-root surface area in each layer
         ! --------------------------------------------------------------------
         
         mean_soil_vwc    = 0._r8
         mean_soil_matpot = 0._r8
         mean_soil_vwcsat = 0._r8
         areaweight       = 0._r8
         
         do jrhiz=1,nlevrhiz
            
            jsoil = jrhiz + jr1-1
            vwc     = bc_in(s)%h2o_liqvol_sl(jsoil)
            psi     = site_hydr%wrf_soil(jrhiz)%p%psi_from_th(vwc)
            vwc_sat = bc_in(s)%watsat_sl(jsoil)
            layer_areaweight = site_hydr%l_aroot_layer(jrhiz)*pi_const*site_hydr%rs1(jrhiz)**2.0
            mean_soil_vwc    = mean_soil_vwc + vwc*layer_areaweight
            mean_soil_vwcsat = mean_soil_vwcsat + vwc_sat*layer_areaweight
            mean_soil_matpot = mean_soil_matpot + psi*layer_areaweight
            areaweight       = areaweight + layer_areaweight

            hio_soilmatpot_sl(io_si,jsoil) = psi
            hio_soilvwc_sl(io_si,jsoil)    = vwc
            hio_soilvwcsat_sl(io_si,jsoil) = vwc_sat
            
         end do
         
         hio_rootwgt_soilvwc(io_si)    = mean_soil_vwc/areaweight
         hio_rootwgt_soilvwcsat(io_si) = mean_soil_vwcsat/areaweight
         hio_rootwgt_soilmatpot(io_si) = mean_soil_matpot/areaweight
         
         hio_rootuptake(io_si) = sum(site_hydr%rootuptake_sl,dim=1)
         hio_rootuptake_sl(io_si,:) = 0._r8
         hio_rootuptake_sl(io_si,jr1:jr2) = site_hydr%rootuptake_sl(1:nlevrhiz)
         hio_rootuptake(io_si) = sum(site_hydr%sapflow_scpf)

         ! Normalization counters
         nplant_scpf(:) = 0._r8
         ncohort_scpf(:) = 0._r8
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            ccohort => cpatch%shortest
            do while(associated(ccohort))
               if ( .not. ccohort%isnew ) then
                  ! Calculate index for the scpf class
                  iscpf = ccohort%size_by_pft_class
                  nplant_scpf(iscpf) = nplant_scpf(iscpf) + ccohort%n
                  ncohort_scpf(iscpf) = ncohort_scpf(iscpf) + 1._r8
               end if
               ccohort => ccohort%taller
            enddo ! cohort loop
            cpatch => cpatch%younger
         end do !patch loop


         ! Generate a histogram of the the iteration counts
         if(print_iterations) then
             iterh2_histy(:) = 0._r8
             cpatch => sites(s)%oldest_patch
             do while(associated(cpatch))
                 ccohort => cpatch%shortest
                 do while(associated(ccohort))
                     ccohort_hydr => ccohort%co_hydr
                     ix = count((iterh2_histx(:)-0.00001_r8) < ccohort_hydr%iterh2 )
                     iterh2_histy(ix) = iterh2_histy(ix) + 1
                     ccohort => ccohort%taller
                 enddo ! cohort loop
                 cpatch => cpatch%younger
             end do !patch loop                     
         end if

         
         do ipft = 1, numpft
            do isz = 1,nlevsclass
               iscpf = (ipft-1)*nlevsclass + isz
               hio_sapflow_scpf(io_si,iscpf)       = site_hydr%sapflow_scpf(isz, ipft)
               hio_rootuptake0_scpf(io_si,iscpf)   = site_hydr%rootuptake0_scpf(isz,ipft)
               hio_rootuptake10_scpf(io_si,iscpf)  = site_hydr%rootuptake10_scpf(isz,ipft)
               hio_rootuptake50_scpf(io_si,iscpf)  = site_hydr%rootuptake50_scpf(isz,ipft)
               hio_rootuptake100_scpf(io_si,iscpf) = site_hydr%rootuptake100_scpf(isz,ipft)
               hio_iterh1_scpf(io_si,iscpf) = 0._r8
               hio_iterh2_scpf(io_si,iscpf) = 0._r8
            end do
         end do

         ipa = 0
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            
            ccohort => cpatch%shortest
            do while(associated(ccohort))

               ccohort_hydr => ccohort%co_hydr
               
               if ( .not. ccohort%isnew ) then

                  ! Calculate index for the scpf class
                  iscpf = ccohort%size_by_pft_class
                  
                  ! scale up cohort fluxes to their sites
                  number_fraction_rate = (ccohort%n / nplant_scpf(iscpf))/dt_tstep
                  
                  ! scale cohorts to mean quantity
                  number_fraction = (ccohort%n / nplant_scpf(iscpf))
                  
                  hio_errh2o_scpf(io_si,iscpf) = hio_errh2o_scpf(io_si,iscpf) + &
                        ccohort_hydr%errh2o * number_fraction_rate ! [kg/indiv/s]
                  
                  hio_tran_scpf(io_si,iscpf) = hio_tran_scpf(io_si,iscpf) + &
                        (ccohort_hydr%qtop) * number_fraction_rate ! [kg/indiv/s]
                  
                  hio_iterh1_scpf(io_si,iscpf)          = hio_iterh1_scpf(io_si,iscpf) + &
                        ccohort_hydr%iterh1/ncohort_scpf(iscpf)
                  
                  hio_iterh2_scpf(io_si,iscpf)          = hio_iterh2_scpf(io_si,iscpf) + &
                        ccohort_hydr%iterh2/ncohort_scpf(iscpf)
                  
                  mean_aroot = sum(ccohort_hydr%th_aroot(:)*ccohort_hydr%v_aroot_layer(:)) / &
                       sum(ccohort_hydr%v_aroot_layer(:))
                  
                  hio_ath_scpf(io_si,iscpf)             = hio_ath_scpf(io_si,iscpf) + &
                       mean_aroot * number_fraction      ! [m3 m-3]
                  
                  hio_tth_scpf(io_si,iscpf)             = hio_tth_scpf(io_si,iscpf) + &
                        ccohort_hydr%th_troot  * number_fraction         ! [m3 m-3]
                  
                  hio_sth_scpf(io_si,iscpf)             = hio_sth_scpf(io_si,iscpf) + &
                        ccohort_hydr%th_ag(2)  * number_fraction        ! [m3 m-3]
                  
                  hio_lth_scpf(io_si,iscpf)             =  hio_lth_scpf(io_si,iscpf) + &
                        ccohort_hydr%th_ag(1)  * number_fraction        ! [m3 m-3]

                  mean_aroot = sum(ccohort_hydr%psi_aroot(:)*ccohort_hydr%v_aroot_layer(:)) / &
                       sum(ccohort_hydr%v_aroot_layer(:))
                  
                  hio_awp_scpf(io_si,iscpf)             = hio_awp_scpf(io_si,iscpf) + &
                       mean_aroot * number_fraction     ! [MPa]
                  
                  hio_twp_scpf(io_si,iscpf)             = hio_twp_scpf(io_si,iscpf) + &
                        ccohort_hydr%psi_troot  * number_fraction       ! [MPa]
                  
                  hio_swp_scpf(io_si,iscpf)             = hio_swp_scpf(io_si,iscpf) + &
                        ccohort_hydr%psi_ag(2)  * number_fraction       ! [MPa]
                  
                  hio_lwp_scpf(io_si,iscpf)             = hio_lwp_scpf(io_si,iscpf) + &
                       ccohort_hydr%psi_ag(1)  * number_fraction       ! [MPa]

                  mean_aroot = sum(ccohort_hydr%ftc_aroot(:)*ccohort_hydr%v_aroot_layer(:)) / &
                       sum(ccohort_hydr%v_aroot_layer(:))
                  hio_aflc_scpf(io_si,iscpf)             = hio_aflc_scpf(io_si,iscpf) + &
                        mean_aroot   * number_fraction     
                  
                  hio_tflc_scpf(io_si,iscpf)             = hio_tflc_scpf(io_si,iscpf) + &
                        ccohort_hydr%ftc_troot  * number_fraction     
                  
                  hio_sflc_scpf(io_si,iscpf)             = hio_sflc_scpf(io_si,iscpf) + &
                       ccohort_hydr%ftc_ag(2)  * number_fraction       
                  
                  hio_lflc_scpf(io_si,iscpf)             = hio_lflc_scpf(io_si,iscpf) + &
                        ccohort_hydr%ftc_ag(1)  * number_fraction   
                  
                  hio_btran_scpf(io_si,iscpf)           = hio_btran_scpf(io_si,iscpf) + &
                        ccohort_hydr%btran  * number_fraction        ! [-]
                  
               endif

               ccohort => ccohort%taller
            enddo ! cohort loop
            ipa = ipa + 1
            cpatch => cpatch%younger
         end do !patch loop

         if(hlm_use_ed_st3.eq.ifalse) then
            do iscpf=1,nlevsclass*numpft
               if( abs(hio_nplant_scpf(io_si, iscpf)-nplant_scpf(iscpf)) > 1.0E-8_r8 ) then
                  write(fates_log(),*) 'numpft:',numpft
                  write(fates_log(),*) 'nlevsclass:',nlevsclass
                  write(fates_log(),*) 'scpf:',iscpf
                  write(fates_log(),*) 'io_si:',io_si
                  write(fates_log(),*) 'hio_nplant_scpf:',hio_nplant_scpf(io_si, iscpf)
                  write(fates_log(),*) 'nplant_scpf:',nplant_scpf(iscpf)
                  write(fates_log(),*) 'nplant check on hio_nplant_scpf fails during hydraulics history updates'
                  call endrun(msg=errMsg(sourcefile, __LINE__))
               end if
            end do
         end if

         if(print_iterations) then
!             print*,' Mean solves: ',sum(hio_iterh2_scpf(io_si,:))/real(count(ncohort_scpf(:)>0._r8),r8), &
!                   ' Mean failures: ',sum(hio_iterh1_scpf(io_si,:))/real(count(ncohort_scpf(:)>0._r8),r8)
             write(fmt_char,'(I2)') iterh2_nhist
             write(fates_log(),fmt='(A,'//fmt_char//'I5)') 'Solves: ',int(iterh2_histy(:))
             !write(*,*) 'Histogram: ',int(iterh2_histy(:))
         end if


         
      enddo ! site loop

    end associate
 
 end subroutine update_history_hydraulics

  ! ====================================================================================
  integer function num_history_vars(this)

    implicit none

    class(fates_history_interface_type), intent(in) :: this

    num_history_vars = this%num_history_vars_
    
  end function num_history_vars
  
  ! ====================================================================================
  
  subroutine initialize_history_vars(this)

    implicit none

    class(fates_history_interface_type), intent(inout) :: this

   ! Determine how many of the history IO variables registered in FATES
   ! are going to be allocated
   call this%define_history_vars(initialize_variables=.false.)

   ! Allocate the list of history output variable objects
   allocate(this%hvars(this%num_history_vars()))
   
   ! construct the object that defines all of the IO variables
   call this%define_history_vars(initialize_variables=.true.)
   
 end subroutine initialize_history_vars
  
  ! ====================================================================================
  
  subroutine define_history_vars(this, initialize_variables)
    
    ! ---------------------------------------------------------------------------------
    ! 
    !                    REGISTRY OF HISTORY OUTPUT VARIABLES
    !
    ! This subroutine is called in two contexts, either in count mode or inialize mode
    ! In count mode, we just walk through the list of registerred variables, compare
    ! if the variable of interest list the current host model and add it to the count
    ! if true.  This count is used just to allocate the variable space.  After this
    ! has been done, we go through the list a second time populating a memory structure.
    ! This phase is the "initialize" phase.  These two phases are differntiated by the
    ! string "callstep", which should be either "count" or "initialize".
    !
    ! Note 1 there are different ways you can flush or initialize the output fields.
    ! If you flush to a native type, (such as zero), the entire slab which covers
    ! indices which may not be relevant to FATES, are flushed to this value.  So
    ! in that case, lakes and crops that are not controlled by FATES will zero'd
    ! and when values are scaled up to the land-grid, the zero's for non FATES will
    ! be included.  This is good and correct if nothing is there.  
    !
    ! But, what if crops exist in the host model and occupy a fraction of the land-surface
    ! shared with natural vegetation? In that case, you want to flush your arrays
    ! with a value that the HLM treats as "do not average"
    ! 
    ! If your HLM makes use of, and you want, INTEGER OUTPUT, pass the flushval as
    ! a real.  The applied flush value will use the NINT() intrinsic function
    ! ---------------------------------------------------------------------------------

    use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
    use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8    
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_coage_pft_r8, site_coage_r8
    use FatesIOVariableKindMod, only : site_height_r8
    use FatesInterfaceTypesMod     , only : hlm_use_planthydro
    
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
    use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8
    use FatesIOVariableKindMod, only : site_elem_r8, site_elpft_r8
    use FatesIOVariableKindMod, only : site_elcwd_r8, site_elage_r8


    implicit none
    
    class(fates_history_interface_type), intent(inout) :: this
    logical, intent(in) :: initialize_variables  ! are we 'count'ing or 'initializ'ing?

    integer :: ivar
    character(len=10) :: tempstring 

    ! A maximum number of 32 characters is allowed for names
    ! We prepend EVERYTHING in FATES with a "FATES_".  This
    ! leaves us with 26 characters to use.
    ! 
    ! FATES_MORTALITY_UNDERSTORY_BY_SCPF
    !       12345678901234567890123456

    
    ivar=0
    
    ! Site level counting variables
    call this%set_history_var(vname='NUMPATCHES', units='none',                &
         long='Total number of ED patches per site', use_default='active',      &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_npatches)

    call this%set_history_var(vname='NUMCOHORTS', units='none',                &
         long='Total number of ED cohorts per site', use_default='active',      &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_ncohorts)
    
    ! Patch variables
    call this%set_history_var(vname='TRIMMING', units='none',                   &
         long='Degree to which canopy expansion is limited by leaf economics',  & 
         use_default='active', &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_trimming)
    
    call this%set_history_var(vname='AREA_PLANT', units='m2',                   &
         long='area occupied by all plants', use_default='active',              &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_area_plant)
    
    call this%set_history_var(vname='AREA_TREES', units='m2',                   &
         long='area occupied by woody plants', use_default='active',            &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_area_trees)

    call this%set_history_var(vname='COLD_STATUS', units='0,1,2', &
          long='Site level cold status, 0=not cold-dec, 1=too cold for leaves, 2=not-too cold',  &
          use_default='active',                                                  &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1, &
          ivar=ivar, initialize=initialize_variables, index = ih_cstatus )

    call this%set_history_var(vname='DROUGHT_STATUS', units='0,1,2,3', &
          long='Site level drought status, <2 too dry for leaves, >=2 not-too dry', &
          use_default='active',                                                  &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1, &
          ivar=ivar, initialize=initialize_variables, index = ih_dstatus)

    call this%set_history_var(vname='GDD', units='degC',  &
         long='site level growing degree days',                &
         use_default='active',                                                 &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_gdd)
    
    call this%set_history_var(vname='NUMCHILLDAYS', units = 'days', &
         long='site level number of chill days', &
         use_default='active',                                                 &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_nchilldays)

    call this%set_history_var(vname='NUMCOLDDAYS', units = 'days', &
         long='site level number of cold days', &
         use_default='active',                                                 &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_ncolddays)

    call this%set_history_var(vname='DAYSINCE_COLDLEAFOFF', units='days', &
         long='site level days elapsed since cold leaf drop', &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_cleafoff)

    call this%set_history_var(vname='DAYSINCE_COLDLEAFON', units='days', &
         long='site level days elapsed since cold leaf flush', &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_cleafon) 

    call this%set_history_var(vname='DAYSINCE_DROUGHTLEAFOFF', units='days', &
         long='site level days elapsed since drought leaf drop', &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_dleafoff)
    
    call this%set_history_var(vname='DAYSINCE_DROUGHTLEAFON', units='days', &
         long='site level days elapsed since drought leaf flush', &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_dleafon)

    call this%set_history_var(vname='MEANLIQVOL_DROUGHTPHEN', units='m3/m3', &
         long='site level mean liquid water volume for drought phen', &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_meanliqvol)

    call this%set_history_var(vname='CANOPY_SPREAD', units='0-1',               &
         long='Scaling factor between tree basal area and canopy area',         &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_canopy_spread)

    
    call this%set_history_var(vname='BIOMASS_BY_PFT', units='gC/m2',                   &
         long='total PFT level biomass', use_default='active',                     &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_biomass_pft )

    call this%set_history_var(vname='LEAFBIOMASS_BY_PFT', units='gC/m2',              &
         long='total PFT level leaf biomass', use_default='active',                &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_leafbiomass_pft )

    call this%set_history_var(vname='STOREBIOMASS_BY_PFT',  units='gC/m2',            &
         long='total PFT level stored biomass', use_default='active',              &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_storebiomass_pft )

    call this%set_history_var(vname='CROWNAREA_BY_PFT',  units='m2/ha',            &
         long='total PFT level crown area', use_default='inactive',              &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_crownarea_pft )
    
    call this%set_history_var(vname='CANOPYCROWNAREA_BY_PFT',  units='m2/ha',            &
         long='total PFT-level canopy-layer crown area', use_default='inactive',     &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_canopycrownarea_pft )
    
    call this%set_history_var(vname='NINDIVS_BY_PFT',  units='indiv / m2',            &
         long='total PFT level number of individuals', use_default='active',       &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_nindivs_pft )

    call this%set_history_var(vname='RECRUITMENT_BY_PFT',  units='indiv/ha/yr',            &
         long='Rate of recruitment by PFT', use_default='active',       &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_recruitment_pft )

    call this%set_history_var(vname='MORTALITY_BY_PFT',  units='indiv/ha/yr',            &
         long='Rate of total mortality by PFT', use_default='active',       &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_mortality_pft )


    
    ! patch age class variables
    call this%set_history_var(vname='PATCHAREA_BY_AGE', units='m2/m2',             &
         long='patch area by age bin', use_default='active',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_area_age )

    call this%set_history_var(vname='LAI_BY_AGE', units='m2/m2',                   &
         long='leaf area index by age bin', use_default='active',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_lai_age )

    call this%set_history_var(vname='CANOPYAREA_BY_AGE', units='m2/m2',             &
         long='canopy area by age bin', use_default='active',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_canopy_area_age )
    
    call this%set_history_var(vname='NCL_BY_AGE', units='--',                   &
         long='number of canopy levels by age bin', use_default='inactive',             &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_ncl_age )

    call this%set_history_var(vname='NPATCH_BY_AGE', units='--',                   &
         long='number of patches by age bin', use_default='inactive',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_npatches_age )

    call this%set_history_var(vname='BIOMASS_BY_AGE', units='kgC/m2',                   &
         long='Total Biomass within a given patch age bin', &
         use_default='inactive',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_biomass_age )
 
    if ( ED_val_comp_excln .lt. 0._r8 ) then ! only valid when "strict ppa" enabled
       tempstring = 'active'
    else
       tempstring = 'inactive'
    endif
    call this%set_history_var(vname='ZSTAR_BY_AGE', units='m',                   &
         long='product of zstar and patch area by age bin (divide by PATCH_AREA_BY_AGE to get mean zstar)', &
         use_default=trim(tempstring),                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_zstar_age )

    
    call this%set_history_var(vname='CANOPY_HEIGHT_DIST', units='m2/m2',                   &
         long='canopy height distribution', use_default='active',                     &
         avgflag='A', vtype=site_height_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_canopy_height_dist_height )

    call this%set_history_var(vname='LEAF_HEIGHT_DIST', units='m2/m2',                   &
         long='leaf height distribution', use_default='active',                     &
         avgflag='A', vtype=site_height_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_leaf_height_dist_height )

   

    ! Secondary forest area and age diagnostics

    call this%set_history_var(vname='SECONDARY_FOREST_FRACTION', units='m2/m2', &
         long='Secondary forest fraction', use_default='inactive', &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_fraction_secondary_forest )

    call this%set_history_var(vname='WOOD_PRODUCT', units='gC/m2', &
         long='Total wood product from logging', use_default='inactive', &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_woodproduct )

    call this%set_history_var(vname='SECONDARY_FOREST_BIOMASS', units='kgC/m2', &
         long='Biomass on secondary lands (per total site area, mult by SECONDARY_FOREST_FRACTION to get per secondary forest area)',&
         use_default='inactive', &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_biomass_secondary_forest )

    call this%set_history_var(vname='SECONDARY_AREA_AGE_ANTHRO', units='m2/m2', &
         long='Secondary forest patch area age distribution since anthropgenic disturbance', &
         use_default='inactive', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_agesince_anthrodist_age )

    call this%set_history_var(vname='SECONDARY_AREA_AGE_ANY', units='m2/m2', &
         long='Secondary forest patch area age distribution since any kind of disturbance', &
         use_default='inactive', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_secondaryforest_area_age )


    ! Fire Variables

    call this%set_history_var(vname='FIRE_NESTEROV_INDEX', units='none',       &
         long='nesterov_fire_danger index', use_default='active',               &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_nesterov_fire_danger)

    call this%set_history_var(vname='FIRE_IGNITIONS', units='number/km2/day',       &
         long='number of ignitions', use_default='active',               &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_nignitions)

    call this%set_history_var(vname='FIRE_FDI', units='none',       &
         long='probability that an ignition will lead to a fire', use_default='active',               &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fdi)

    call this%set_history_var(vname='FIRE_ROS', units='m/min',                 &
         long='fire rate of spread m/min', use_default='active',                &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_spitfire_ros)

    call this%set_history_var(vname='FIRE_ROS_AREA_PRODUCT', units='m/min',                 &
         long='product of fire rate of spread (m/min) and burned area (fraction)--divide by FIRE_AREA to get burned-area-weighted-mean ROS', use_default='active',                &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_ros_area_product)

    call this%set_history_var(vname='EFFECT_WSPEED', units='none',             &
         long ='effective windspeed for fire spread', use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_effect_wspeed )

    call this%set_history_var(vname='FIRE_TFC_ROS', units='kgC/m2',              &
         long ='total fuel consumed', use_default='active',                     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_tfc_ros )

    call this%set_history_var(vname='FIRE_TFC_ROS_AREA_PRODUCT', units='kgC/m2',              &
         long ='product of total fuel consumed and burned area--divide by FIRE_AREA to get burned-area-weighted-mean TFC', use_default='active',                     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_tfc_ros_area_product )

    call this%set_history_var(vname='FIRE_INTENSITY', units='kJ/m/s',          &
         long='spitfire fire intensity: kJ/m/s', use_default='active',          &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_intensity )

    call this%set_history_var(vname='FIRE_INTENSITY_AREA_PROD', units='kJ/m/s',          &
         long='spitfire product of fire intensity and burned area (divide by FIRE_AREA to get area-weighted mean intensity)', &
         use_default='active',          &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_intensity_area_product )

    call this%set_history_var(vname='FIRE_AREA', units='fraction',             &
         long='spitfire fire area burn fraction', use_default='active',                    &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_area )

    call this%set_history_var(vname='FIRE_FUEL_MEF', units='m',                &
         long='spitfire fuel moisture',  use_default='active',                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_mef )

    call this%set_history_var(vname='FIRE_FUEL_BULKD', units='kg biomass/m3',              &
         long='spitfire fuel bulk density',  use_default='active',              &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_bulkd )

    call this%set_history_var(vname='FIRE_FUEL_EFF_MOIST', units='m',          &
         long='spitfire fuel moisture', use_default='active',                   &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_eff_moist )

    call this%set_history_var(vname='FIRE_FUEL_SAV', units='per m',                &
         long='spitfire fuel surface/volume ',  use_default='active',           &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_sav )

    call this%set_history_var(vname='SUM_FUEL', units='gC m-2',                &
         long='total ground fuel related to ros (omits 1000hr fuels)',          & 
         use_default='active',                                                  & 
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_sum_fuel )

    call this%set_history_var(vname='FRAGMENTATION_SCALER', units='unitless (0-1)',                &
         long='factor by which litter/cwd fragmentation proceeds relative to max rate',          & 
         use_default='active',                                                  & 
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fragmentation_scaler )

    call this%set_history_var(vname='FUEL_MOISTURE_BY_NFSC', units='-',                &
         long='spitfire size-resolved fuel moisture', use_default='active',       &
         avgflag='A', vtype=site_fuel_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_moisture_fuel )

    call this%set_history_var(vname='FUEL_AMOUNT_BY_NFSC', units='kg C / m2',                &
         long='spitfire size-resolved fuel quantity', use_default='active',       &
         avgflag='A', vtype=site_fuel_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fuel_amount_fuel )

    call this%set_history_var(vname='AREA_BURNT_BY_AGE', units='m2/m2', &
         long='spitfire area burnt by patch age (divide by patch_area_by_age to get burnt fraction by age)', &
         use_default='active', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_area_burnt_age )

    call this%set_history_var(vname='FIRE_INTENSITY_BY_AGE', units='kJ/m/2', &
         long='product of fire intensity and burned area, resolved by patch age (so divide by AREA_BURNT_BY_PATCH_AGE to get burned-area-weighted-average intensity', &
         use_default='active', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_intensity_age )

    call this%set_history_var(vname='SUM_FUEL_BY_AGE', units='gC / m2 of site area', &
         long='spitfire ground fuel related to ros (omits 1000hr fuels) within each patch age bin (divide by patch_area_by_age to get fuel per unit area of that-age patch)', &
         use_default='active', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_sum_fuel_age )

    call this%set_history_var(vname='BURNT_FRACAREAPROD_BY_NFSC', units='fraction', &
         long='product of fraction of fuel burnt and burned area (divide by FIRE_AREA to get burned-area-weighted mean fraction fuel burnt)', &
         use_default='active', &
         avgflag='A', vtype=site_fuel_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_burnt_frac_litter_fuel )


    ! Litter Variables

    call this%set_history_var(vname='LITTER_IN', units='gC m-2 s-1',           &
         long='FATES litter flux in',  use_default='active',                   &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_in )

    call this%set_history_var(vname='LITTER_OUT', units='gC m-2 s-1',          &
         long='FATES litter flux out',  use_default='active',                  & 
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_out )

    call this%set_history_var(vname='SEED_BANK', units='gC m-2',               &
         long='Total Seed Mass of all PFTs',  use_default='active',             &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seed_bank )

    call this%set_history_var(vname='SEEDS_IN', units='gC m-2 s-1',            &
         long='Seed Production Rate',  use_default='active',                    &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seeds_in )

    call this%set_history_var(vname='LITTER_IN_BY_EL', units='kg ha-1 d-1',         &
         long='FATES litter flux in',  use_default='active',                      &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_in_elem )

    call this%set_history_var(vname='LITTER_OUT_BY_EL', units='kg ha-1 d-1',         &
         long='FATES litter flux out (fragmentation only)',  use_default='active',                      & 
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_out_elem )

    call this%set_history_var(vname='SEED_BANK_BY_EL', units='kg ha-1',             &
         long='Total Seed Mass of all PFTs',  use_default='active',               &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seed_bank_elem )

    call this%set_history_var(vname='SEEDS_IN_LOCAL_BY_EL', units='kg ha-1 d-1',     &
         long='Within Site Seed Production Rate',  use_default='active',           &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seeds_in_local_elem )

    call this%set_history_var(vname='SEEDS_IN_EXTERN_BY_EL', units='kg ha-1 d-1',     &
         long='External Seed Influx Rate',  use_default='active',                   &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seeds_in_extern_elem )

    call this%set_history_var(vname='SEED_GERM_BY_EL', units='kg ha-1 d-1',          &
         long='Seed mass converted into new cohorts', use_default='active',        &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seed_germ_elem )

    call this%set_history_var(vname='SEED_DECAY_BY_EL', units='kg ha-1 d-1',           &
         long='Seed mass decay', use_default='active',                          &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seed_decay_elem )

    
    ! SITE LEVEL CARBON STATE VARIABLES
    call this%set_history_var(vname='STOREC', units='kgC ha-1',                      &
         long='Total carbon in live plant storage', use_default='active',          &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_storec )
    
    call this%set_history_var(vname='TOTVEGC', units='kgC ha-1',                     &
         long='Total carbon in live plants', use_default='active',                 &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_totvegc )
    
    call this%set_history_var(vname='SAPWC', units='kgC ha-1',                       &
         long='Total carbon in live plant sapwood', use_default='active',          &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_sapwc )
    
    call this%set_history_var(vname='LEAFC', units='kgC ha-1',                       &
         long='Total carbon in live plant leaves', use_default='active',           &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_leafc )
    
    call this%set_history_var(vname='FNRTC', units='kgC ha-1',                       &
         long='Total carbon in live plant fine-roots', use_default='active',       &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fnrtc )

    call this%set_history_var(vname='REPROC', units='kgC ha-1',                          &
         long='Total carbon in live plant reproductive tissues', use_default='active', &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,       &
         ivar=ivar, initialize=initialize_variables, index = ih_reproc )

    call this%set_history_var(vname='CEFFLUX', units='kgC/ha/day', &
         long='carbon efflux, root to soil', use_default='active', &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cefflux )
    

    nitrogen_active_if: if(any(element_list(:)==nitrogen_element)) then
       call this%set_history_var(vname='STOREN', units='kgN ha-1',                      &
            long='Total nitrogen in live plant storage', use_default='active',          &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
            ivar=ivar, initialize=initialize_variables, index = ih_storen )
       
       call this%set_history_var(vname='TOTVEGN', units='kgN ha-1',                     &
            long='Total nitrogen in live plants', use_default='active',                 &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
            ivar=ivar, initialize=initialize_variables, index = ih_totvegn )
       
       call this%set_history_var(vname='SAPWN', units='kgN ha-1',                       &
            long='Total nitrogen in live plant sapwood', use_default='active',          &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
            ivar=ivar, initialize=initialize_variables, index = ih_sapwn )
       
       call this%set_history_var(vname='LEAFN', units='kgN ha-1',                       &
            long='Total nitrogen in live plant leaves', use_default='active',           &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
            ivar=ivar, initialize=initialize_variables, index = ih_leafn )
       
       call this%set_history_var(vname='FNRTN', units='kgN ha-1',                       &
            long='Total nitrogen in live plant fine-roots', use_default='active',       &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
            ivar=ivar, initialize=initialize_variables, index = ih_fnrtn )
       
       call this%set_history_var(vname='REPRON', units='kgN ha-1',                          &
            long='Total nitrogen in live plant reproductive tissues', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,       &
            ivar=ivar, initialize=initialize_variables, index = ih_repron )

       call this%set_history_var(vname='NUPTAKE', units='kgN d-1 ha-1',                          &
            long='Total nitrogen uptake by plants per sq meter per day', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,       &
            ivar=ivar, initialize=initialize_variables, index = ih_nuptake )

       call this%set_history_var(vname='NEFFLUX', units='kgN d-1 ha-1',                          &
            long='Nitrogen effluxed from plant (unused)', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,       &
            ivar=ivar, initialize=initialize_variables, index = ih_nefflux )

       call this%set_history_var(vname='NNEED_GROW', units='kgN d-1 ha-1',                          &
            long='(Approx) plant nitrogen needed to satisfy growth', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,       &
            ivar=ivar, initialize=initialize_variables, index = ih_nneedgrow )

       call this%set_history_var(vname='NNEED_MAX', units='kgN d-1 ha-1',                          &
            long='(Approx) plant nitrogen needed to reach maximum capacity', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,       &
            ivar=ivar, initialize=initialize_variables, index = ih_nneedmax )
       
    end if nitrogen_active_if

    
    phosphorus_active_if: if(any(element_list(:)==phosphorus_element)) then
       call this%set_history_var(vname='STOREP', units='kgP ha-1',                      &
            long='Total phosphorus in live plant storage', use_default='active',          &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
            ivar=ivar, initialize=initialize_variables, index = ih_storep )
       
       call this%set_history_var(vname='TOTVEGP', units='kgP ha-1',                     &
            long='Total phosphorus in live plants', use_default='active',                 &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
            ivar=ivar, initialize=initialize_variables, index = ih_totvegp )
       
       call this%set_history_var(vname='SAPWP', units='kgP ha-1',                       &
            long='Total phosphorus in live plant sapwood', use_default='active',          &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
            ivar=ivar, initialize=initialize_variables, index = ih_sapwp )
       
       call this%set_history_var(vname='LEAFP', units='kgP ha-1',                       &
            long='Total phosphorus in live plant leaves', use_default='active',           &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
            ivar=ivar, initialize=initialize_variables, index = ih_leafp )
       
       call this%set_history_var(vname='FNRTP', units='kgP ha-1',                       &
            long='Total phosphorus in live plant fine-roots', use_default='active',       &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
            ivar=ivar, initialize=initialize_variables, index = ih_fnrtp )
       
       call this%set_history_var(vname='REPROP', units='kgP ha-1',                          &
            long='Total phosphorus in live plant reproductive tissues', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,       &
            ivar=ivar, initialize=initialize_variables, index = ih_reprop )

       call this%set_history_var(vname='PUPTAKE', units='kgP ha-1 d-1',                          &
            long='Total phosphorus uptake by plants per sq meter per day', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,       &
            ivar=ivar, initialize=initialize_variables, index = ih_puptake )

       call this%set_history_var(vname='PEFFLUX', units='kgP ha-1 d-1',                          &
            long='Phosphorus effluxed from plant (unused)', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,       &
            ivar=ivar, initialize=initialize_variables, index = ih_pefflux )
       
       call this%set_history_var(vname='PNEED_GROW', units='kgP ha-1 d-1',                          &
            long='Plant phosphorus needed to satisfy growth', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,       &
            ivar=ivar, initialize=initialize_variables, index = ih_pneedgrow )

       call this%set_history_var(vname='PNEED_MAX', units='kgP ha-1 d-1', &
            long='Plant phosphorus needed to reach maximum capacity', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,       &
            ivar=ivar, initialize=initialize_variables, index = ih_pneedmax )
       
       
    end if phosphorus_active_if

    
    call this%set_history_var(vname='AGB', units='gC m-2',                  &
         long='Aboveground biomass',  use_default='active',                           &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_agb )

    call this%set_history_var(vname='TOTVEGC_CANOPY', units='gC m-2',                   &
         long='Biomass of canopy plants',  use_default='active',                            &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_canopy_biomass )

    call this%set_history_var(vname='TOTVEGC_USTORY', units='gC m-2',                   &
         long='Biomass of understory plants',  use_default='active',                            &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_understory_biomass )

    ! disturbance rates
    call this%set_history_var(vname='PRIMARY_PATCHFUSION_ERROR', units='m2 m-2 d-1',                   &
         long='Error in total primary lands associated with patch fusion',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_primaryland_fusion_error )

    call this%set_history_var(vname='DISTURBANCE_RATE_P2P', units='m2 m-2 d-1',                   &
         long='Disturbance rate from primary to primary lands',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_disturbance_rate_p2p )

    call this%set_history_var(vname='DISTURBANCE_RATE_P2S', units='m2 m-2 d-1',                   &
         long='Disturbance rate from primary to secondary lands',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_disturbance_rate_p2s )

    call this%set_history_var(vname='DISTURBANCE_RATE_S2S', units='m2 m-2 d-1',                   &
         long='Disturbance rate from secondary to secondary lands',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_disturbance_rate_s2s )

    call this%set_history_var(vname='DISTURBANCE_RATE_FIRE', units='m2 m-2 d-1',                   &
         long='Disturbance rate from fire',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_disturbance_rate )

    call this%set_history_var(vname='DISTURBANCE_RATE_LOGGING', units='m2 m-2 d-1',                   &
         long='Disturbance rate from logging',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_logging_disturbance_rate )

    call this%set_history_var(vname='DISTURBANCE_RATE_TREEFALL', units='m2 m-2 d-1',                   &
         long='Disturbance rate from treefall',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fall_disturbance_rate )

    call this%set_history_var(vname='DISTURBANCE_RATE_POTENTIAL', units='m2 m-2 d-1',                   &
         long='Potential (i.e., including unresolved) disturbance rate',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_potential_disturbance_rate )

    call this%set_history_var(vname='HARVEST_CARBON_FLUX', units='kg C m-2 d-1',                   &
         long='Harvest carbon flux',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_harvest_carbonflux )

    ! Canopy Resistance 

    call this%set_history_var(vname='COND_STOMATA', units='umol m-2 s-1',                   &
         long='mean stomatal conductance', use_default='active',                   &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_c_stomata )

    call this%set_history_var(vname='COND_LBLAYER', units='umol m-2 s-1',                   &
         long='mean leaf boundary layer conductance', use_default='active',                   &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_c_lblayer )


    ! Ecosystem Carbon Fluxes (updated rapidly, upfreq=2)

    call this%set_history_var(vname='NPP', units='gC/m^2/s',                &
         long='net primary production',  use_default='active',      &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_npp )

    call this%set_history_var(vname='GPP', units='gC/m^2/s',                   &
         long='gross primary production',  use_default='active',                &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp )

    call this%set_history_var(vname='AR', units='gC/m^2/s',                 &
         long='autotrophic respiration', use_default='active',                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_aresp )

    call this%set_history_var(vname='GROWTH_RESP', units='gC/m^2/s',           &
         long='growth respiration', use_default='active',                       &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_growth_resp )

    call this%set_history_var(vname='MAINT_RESP', units='gC/m^2/s',            &
         long='maintenance respiration', use_default='active',                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_maint_resp )

    ! Canopy resistance 

    call this%set_history_var(vname='COND_STOMATA_BY_AGE', units='umol m-2 s-1',                   &
         long='mean stomatal conductance - by patch age', use_default='inactive', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_c_stomata_age )

    call this%set_history_var(vname='COND_LBLAYER_BY_AGE', units='umol m-2 s-1',                   &
         long='mean leaf boundary layer conductance - by patch age', use_default='inactive', &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_c_lblayer_age )

    ! fast fluxes by age bin
    call this%set_history_var(vname='NPP_BY_AGE', units='gC/m^2/s',                   &
         long='net primary productivity by age bin', use_default='inactive',           &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2, &
         ivar=ivar, initialize=initialize_variables, index = ih_npp_age )

    call this%set_history_var(vname='GPP_BY_AGE', units='gC/m^2/s',                   &
         long='gross primary productivity by age bin', use_default='inactive',         &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2, &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp_age )

    ! fast fluxes separated canopy/understory
    call this%set_history_var(vname='GPP_CANOPY', units='gC/m^2/s',                   &
         long='gross primary production of canopy plants',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp_canopy )

    call this%set_history_var(vname='AR_CANOPY', units='gC/m^2/s',                 &
         long='autotrophic respiration of canopy plants', use_default='active',       &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_ar_canopy )

    call this%set_history_var(vname='GPP_USTORY', units='gC/m^2/s',                   &
         long='gross primary production of understory plants',  use_default='active',     &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp_understory )

    call this%set_history_var(vname='AR_USTORY', units='gC/m^2/s',                 &
         long='autotrophic respiration of understory plants', use_default='active',       &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_ar_understory )


    ! fast radiative fluxes resolved through the canopy
    call this%set_history_var(vname='PARSUN_BY_CNLF', units='W/m2',                 &
         long='PAR absorbed in the sun by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsun_z_cnlf )

    call this%set_history_var(vname='PARSHA_BY_CNLF', units='W/m2',                 &
         long='PAR absorbed in the shade by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsha_z_cnlf )

    call this%set_history_var(vname='PARSUN_BY_CNLFPFT', units='W/m2',                 &
         long='PAR absorbed in the sun by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsun_z_cnlfpft )

    call this%set_history_var(vname='PARSHA_BY_CNLFPFT', units='W/m2',                 &
         long='PAR absorbed in the shade by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsha_z_cnlfpft )

    call this%set_history_var(vname='PARSUN_CAN', units='W/m2',                 &
         long='PAR absorbed in the sun by top leaf layer in each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsun_top_can )

    call this%set_history_var(vname='PARSHA_CAN', units='W/m2',                 &
         long='PAR absorbed in the shade by top leaf layer in each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsha_top_can )

    call this%set_history_var(vname='LAISUN_BY_CNLF', units='m2/m2',                 &
         long='LAI in the sun by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisun_z_cnlf )

    call this%set_history_var(vname='LAISHA_BY_CNLF', units='m2/m2',                 &
         long='LAI in the shade by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisha_z_cnlf )

    call this%set_history_var(vname='LAISUN_BY_CNLFPFT', units='m2/m2',                 &
         long='LAI in the sun by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisun_z_cnlfpft )

    call this%set_history_var(vname='LAISHA_BY_CNLFPFT', units='m2/m2',                 &
         long='LAI in the shade by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisha_z_cnlfpft )

    call this%set_history_var(vname='LAISUN_TOP_CAN', units='m2/m2',                 &
         long='LAI in the sun by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisun_top_can )

    call this%set_history_var(vname='LAISHA_TOP_CAN', units='m2/m2',                 &
         long='LAI in the shade by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisha_top_can )

    call this%set_history_var(vname='FABD_SUN_BY_CNLFPFT', units='fraction',                 &
         long='sun fraction of direct light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sun_cnlfpft )

    call this%set_history_var(vname='FABD_SHA_BY_CNLFPFT', units='fraction',                 &
         long='shade fraction of direct light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sha_cnlfpft )

    call this%set_history_var(vname='FABI_SUN_BY_CNLFPFT', units='fraction',                 &
         long='sun fraction of indirect light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sun_cnlfpft )

    call this%set_history_var(vname='FABI_SHA_BY_CNLFPFT', units='fraction',                 &
         long='shade fraction of indirect light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sha_cnlfpft )

    call this%set_history_var(vname='FABD_SUN_BY_CNLF', units='fraction',                 &
         long='sun fraction of direct light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sun_cnlf )

    call this%set_history_var(vname='FABD_SHA_BY_CNLF', units='fraction',                 &
         long='shade fraction of direct light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sha_cnlf )

    call this%set_history_var(vname='FABI_SUN_CNLF', units='fraction',                 &
         long='sun fraction of indirect light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sun_cnlf )

    call this%set_history_var(vname='FABI_SHA_BY_CNLF', units='fraction',                 &
         long='shade fraction of indirect light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sha_cnlf )

    call this%set_history_var(vname='PARPROF_DIR_BY_CNLFPFT', units='W/m2',                 &
         long='Radiative profile of direct PAR through each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parprof_dir_cnlfpft )

    call this%set_history_var(vname='PARPROF_DIF_BY_CNLFPFT', units='W/m2',                 &
         long='Radiative profile of diffuse PAR through each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parprof_dif_cnlfpft )

    call this%set_history_var(vname='PARPROF_DIR_BY_CNLF', units='W/m2',                 &
         long='Radiative profile of direct PAR through each canopy and leaf layer (averaged across PFTs)', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parprof_dir_cnlf )

    call this%set_history_var(vname='PARPROF_DIF_BY_CNLF', units='W/m2',                 &
         long='Radiative profile of diffuse PAR through each canopy and leaf layer (averaged across PFTs)', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parprof_dif_cnlf )

    call this%set_history_var(vname='FABD_SUN_TOPLF_BY_CN', units='fraction',                 &
         long='sun fraction of direct light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sun_top_can )

    call this%set_history_var(vname='FABD_SHA_TOPLF_BY_CN', units='fraction',                 &
         long='shade fraction of direct light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sha_top_can )

    call this%set_history_var(vname='FABI_SUN_TOPLF_BY_CN', units='fraction',                 &
         long='sun fraction of indirect light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sun_top_can )

    call this%set_history_var(vname='FABI_SHA_TOPLF_BY_CN', units='fraction',                 &
         long='shade fraction of indirect light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sha_top_can )

    !!! canopy-resolved fluxes and structure
    call this%set_history_var(vname='NET_C_UPTAKE_BY_CNLF', units='gC/m2/s',                 &
         long='net carbon uptake by each canopy and leaf layer per unit ground area (i.e. divide by CROWNAREA_CNLF to make per leaf area)', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_ts_net_uptake_cnlf )

    call this%set_history_var(vname='CROWNAREA_BY_CNLF', units='m2/m2',                 &
         long='total crown area that is occupied by leaves in each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_crownarea_cnlf )

    call this%set_history_var(vname='CROWNAREA_BY_CN', units='m2/m2',                 &
         long='total crown area in each canopy layer', &
         use_default='active',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ELM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_crownarea_can )

    ! slow carbon fluxes associated with mortality from or transfer betweeen canopy and understory
    call this%set_history_var(vname='DEMOTION_CFLUX', units = 'gC/m2/s',               &
          long='demotion-associated biomass carbon flux from canopy to understory', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_demotion_carbonflux )

    call this%set_history_var(vname='PROMOTION_CFLUX', units = 'gC/m2/s',               &
          long='promotion-associated biomass carbon flux from understory to canopy', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_promotion_carbonflux )

    call this%set_history_var(vname='MORTALITY_CFLUX_CANOPY', units = 'gC/m2/s',               &
          long='flux of biomass carbon from live to dead pools from mortality of canopy plants', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_canopy_mortality_carbonflux )

    call this%set_history_var(vname='MORTALITY_CFLUX_USTORY', units = 'gC/m2/s',               &
          long='flux of biomass carbon from live to dead pools from mortality of understory plants',use_default='active',&
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_understory_mortality_carbonflux )

    ! size class by age dimensioned variables
    call this%set_history_var(vname='NPLANT_BY_SCAG',units = 'plants/ha',               &
          long='number of plants per hectare in each size x age class', use_default='active',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_scag )

    call this%set_history_var(vname='NPLANT_CANOPY_BY_SCAG',units = 'plants/ha',               &
          long='number of plants per hectare in canopy in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_canopy_scag )

    call this%set_history_var(vname='NPLANT_USTORY_BY_SCAG',units = 'plants/ha',               &
          long='number of plants per hectare in understory in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_understory_scag )

    call this%set_history_var(vname='DDBH_CANOPY_BY_SCAG',units = 'cm/yr/ha',               &
          long='growth rate of canopy plantsnumber of plants per hectare in canopy in each size x age class', &
          use_default='inactive', avgflag='A', vtype=site_scag_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_canopy_scag )

    call this%set_history_var(vname='DDBH_USTORY_BY_SCAG',units = 'cm/yr/ha',               &
          long='growth rate of understory plants in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_understory_scag )

    call this%set_history_var(vname='MORTALITY_CANOPY_BY_SCAG',units = 'plants/ha/yr',               &
          long='mortality rate of canopy plants in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_canopy_scag )

    call this%set_history_var(vname='MORTALITY_STORY_BY_SCAG',units = 'plants/ha/yr',               &
          long='mortality rate of understory plantsin each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_understory_scag )

    
    
    ! size x age x pft dimensioned
    call this%set_history_var(vname='NPLANT_BY_SCAGPFT',units = 'plants/ha',               &
          long='number of plants per hectare in each size x age x pft class', use_default='inactive',   &
          avgflag='A', vtype=site_scagpft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_scagpft )

    ! age x pft dimensioned
    call this%set_history_var(vname='NPP_BY_AGEPFT',units = 'kgC/m2/yr',               &
          long='NPP per PFT in each age bin', use_default='inactive',   &
          avgflag='A', vtype=site_agepft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_agepft )

    call this%set_history_var(vname='BIOMASS_BY_AGEPFT',units = 'kg C / m2',               &
          long='biomass per PFT in each age bin', use_default='inactive',   &
          avgflag='A', vtype=site_agepft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_biomass_agepft )

    call this%set_history_var(vname='SCORCH_HEIGHT_BY_AGEPFT',units = 'm',               &
          long='SPITFIRE Flame Scorch Height (calculated per PFT in each patch age bin)', &
          use_default='active',   &
          avgflag='A', vtype=site_agepft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_scorch_height_agepft )


    ! Carbon Flux (grid dimension x scpf) (THESE ARE DEFAULT INACTIVE!!!
    !                                     (BECAUSE THEY TAKE UP SPACE!!!
    ! ===================================================================================

    call this%set_history_var(vname='GPP_BY_SCPF', units='kgC/m2/yr',            &
          long='gross primary production by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_gpp_scpf )

    call this%set_history_var(vname='GPP_CANOPY_BY_SCPF', units='kgC/m2/yr',            &
          long='gross primary production of canopy plants by pft/size ', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_gpp_canopy_scpf )

    call this%set_history_var(vname='AR_CANOPY_BY_SCPF', units='kgC/m2/yr',            &
          long='autotrophic respiration of canopy plants by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ar_canopy_scpf )

    call this%set_history_var(vname='GPP_USTORY_BY_SCPF', units='kgC/m2/yr',            &
          long='gross primary production of understory plants by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_gpp_understory_scpf )

    call this%set_history_var(vname='AR_USTORY_BY_SCPF', units='kgC/m2/yr',            &
          long='autotrophic respiration of understory plants by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ar_understory_scpf )

    call this%set_history_var(vname='NPP_BY_SCPF', units='kgC/m2/yr',            &
          long='total net primary production by pft/size', use_default='inactive',       &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_totl_scpf )

    call this%set_history_var(vname='NPP_LEAF_BY_SCPF', units='kgC/m2/yr',       &
          long='NPP flux into leaves by pft/size', use_default='inactive',               &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_leaf_scpf )

   call this%set_history_var(vname='NPP_SEED_BY_SCPF', units='kgC/m2/yr',       &
         long='NPP flux into seeds by pft/size', use_default='inactive',                &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_seed_scpf )

   call this%set_history_var(vname='NPP_FNRT_BY_SCPF', units='kgC/m2/yr',       &
         long='NPP flux into fine roots by pft/size', use_default='inactive',           &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_fnrt_scpf )

   call this%set_history_var(vname='NPP_BGSW_BY_SCPF', units='kgC/m2/yr',       &
         long='NPP flux into below-ground sapwood by pft/size', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bgsw_scpf )

   call this%set_history_var(vname='NPP_BGDW_BY_SCPF', units='kgC/m2/yr',       &
         long='NPP flux into below-ground deadwood by pft/size', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bgdw_scpf )

   call this%set_history_var(vname='NPP_AGSW_BY_SCPF', units='kgC/m2/yr',       &
         long='NPP flux into above-ground sapwood by pft/size', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_agsw_scpf )

   call this%set_history_var(vname = 'NPP_AGDW_BY_SCPF', units='kgC/m2/yr',    &
         long='NPP flux into above-ground deadwood by pft/size', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_agdw_scpf )

   call this%set_history_var(vname = 'NPP_STOR_BY_SCPF', units='kgC/m2/yr',    &
         long='NPP flux into storage by pft/size', use_default='inactive',              &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_stor_scpf )

    call this%set_history_var(vname='DDBH_BY_SCPF', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='inactive',          &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,   &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_scpf )

    call this%set_history_var(vname='GROWTHFLUX_BY_SCPF', units = 'n/yr/ha',         &
          long='flux of individuals into a given size class bin via growth and recruitment',use_default='inactive',          &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,   &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_growthflux_scpf )

    call this%set_history_var(vname='GROWTHFLUX_FUSION_BY_SCPF', units = 'n/yr/ha',         &
          long='flux of individuals into a given size class bin via fusion',use_default='inactive',          &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,   &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_growthflux_fusion_scpf )

    call this%set_history_var(vname='DDBH_CANOPY_BY_SCPF', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_canopy_scpf )

    call this%set_history_var(vname='DDBH_USTORY_BY_SCPF', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_understory_scpf )

    call this%set_history_var(vname='BA_BY_SCPF', units = 'm2/ha',               &
          long='basal area by pft/size', use_default='inactive',   &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ba_scpf )

    call this%set_history_var(vname='AGB_BY_SCPF', units = 'kgC/m2', &
         long='Aboveground biomass by pft/size', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8, &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_agb_scpf ) 

    call this%set_history_var(vname='NPLANT_BY_SCPF', units = 'N/ha',         &
          long='stem number density by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_scpf )

    call this%set_history_var(vname='NPLANT_BY_CAPF', units = 'N/ha',       &
         long='stem number density by pft/coage', use_default='inactive', &
         avgflag='A', vtype=site_coage_pft_r8, hlms='CLM:ELM',flushval=0.0_r8,     &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_capf )

    call this%set_history_var(vname='M1_BY_SCPF', units = 'N/ha/yr',          &
          long='background mortality by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m1_scpf )
    
    call this%set_history_var(vname='M2_BY_SCPF', units = 'N/ha/yr',          &
          long='hydraulic mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m2_scpf )

    call this%set_history_var(vname='M3_BY_SCPF', units = 'N/ha/yr',          &
          long='carbon starvation mortality by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m3_scpf )

    call this%set_history_var(vname='M4_BY_SCPF', units = 'N/ha/yr',          &
          long='impact mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m4_scpf )

    call this%set_history_var(vname='M5_BY_SCPF', units = 'N/ha/yr',          &
          long='fire mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m5_scpf )

    call this%set_history_var(vname='CROWNFIREMORT_BY_SCPF', units = 'N/ha/yr',          &
          long='crown fire mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_crownfiremort_scpf )

    call this%set_history_var(vname='CAMBIALFIREMORT_BY_SCPF', units = 'N/ha/yr',          &
          long='cambial fire mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cambialfiremort_scpf )

    call this%set_history_var(vname='M6_BY_SCPF', units = 'N/ha/yr',          &
          long='termination mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m6_scpf )

    call this%set_history_var(vname='M7_BY_SCPF', units = 'N/ha/event',               &
          long='logging mortality by pft/size',use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m7_scpf )

    call this%set_history_var(vname='M8_BY_SCPF', units = 'N/ha/yr',          &
          long='freezing mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m8_scpf )

    call this%set_history_var(vname='M9_BY_SCPF', units = 'N/ha/yr',          &
          long='senescence mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m9_scpf )

    call this%set_history_var(vname='M10_BY_SCPF', units = 'N/ha/yr',         &
         long='age senescence mortality by pft/size',use_default='inactive', &
         avgflag='A', vtype =site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,     &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m10_scpf )
    
    call this%set_history_var(vname='M10_BY_CAPF',units='N/ha/yr',         &
         long='age senescence mortality by pft/cohort age',use_default='inactive', &
         avgflag='A', vtype =site_coage_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,         &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index =ih_m10_capf )

    call this%set_history_var(vname='MORTALITY_CANOPY_BY_SCPF', units = 'N/ha/yr',          &
          long='total mortality of canopy plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_canopy_scpf )

    call this%set_history_var(vname='C13DISC_BY_SCPF', units = 'per mil',               &
         long='C13 discrimination by pft/size',use_default='inactive',           &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_c13disc_scpf ) 

    call this%set_history_var(vname='BSTOR_CANOPY_BY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in storage pools of canopy plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bstor_canopy_scpf )

    call this%set_history_var(vname='BLEAF_CANOPY_BY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in leaf of canopy plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bleaf_canopy_scpf )

    call this%set_history_var(vname='NPLANT_CANOPY_BY_SCPF', units = 'N/ha',         &
          long='stem number of canopy plants density by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_canopy_scpf )

    call this%set_history_var(vname='MORTALITY_USTORY_BY_SCPF', units = 'N/ha/yr',          &
          long='total mortality of understory plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_understory_scpf )

    call this%set_history_var(vname='BSTOR_USTORY_BY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in storage pools of understory plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bstor_understory_scpf )

    call this%set_history_var(vname='BLEAF_USTORY_BY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in leaf of understory plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bleaf_understory_scpf )

    call this%set_history_var(vname='NPLANT_USTORY_BY_SCPF', units = 'N/ha',         &
          long='stem number of understory plants density by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_understory_scpf )

    call this%set_history_var(vname='CWD_AG_BY_CWDSC', units='gC/m^2', &
          long='size-resolved AG CWD stocks', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_ag_cwdsc )

    call this%set_history_var(vname='CWD_BG_BY_CWDSC', units='gC/m^2', &
          long='size-resolved BG CWD stocks', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_bg_cwdsc )

    call this%set_history_var(vname='CWD_AG_IN_BY_CWDSC', units='gC/m^2/y', &
          long='size-resolved AG CWD input', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_ag_in_cwdsc )

    call this%set_history_var(vname='CWD_BG_IN_BY_CWDSC', units='gC/m^2/y', &
          long='size-resolved BG CWD input', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_bg_in_cwdsc )

    call this%set_history_var(vname='CWD_AG_OUT_BY_CWDSC', units='gC/m^2/y', &
          long='size-resolved AG CWD output', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_ag_out_cwdsc )

    call this%set_history_var(vname='CWD_BG_OUT_BY_CWDSC', units='gC/m^2/y', &
          long='size-resolved BG CWD output', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_bg_out_cwdsc )

    ! Size structured diagnostics that require rapid updates (upfreq=2)

    call this%set_history_var(vname='AR_BY_SCPF',units = 'kgC/m2/yr',          &
          long='total autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_scpf )
    
    call this%set_history_var(vname='AR_GROW_BY_SCPF',units = 'kgC/m2/yr',          &
          long='growth autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_grow_scpf )

    call this%set_history_var(vname='AR_MAINT_BY_SCPF',units = 'kgC/m2/yr',          &
          long='maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_maint_scpf )

    call this%set_history_var(vname='AR_DARKM_BY_SCPF',units = 'kgC/m2/yr',          &
          long='dark portion of maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8,hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_darkm_scpf )

    call this%set_history_var(vname='AR_AGSAPM_BY_SCPF',units = 'kgC/m2/yr',          &
          long='above-ground sapwood maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8,hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_agsapm_scpf )
    
    call this%set_history_var(vname='AR_CROOTM_BY_SCPF',units = 'kgC/m2/yr',          &
          long='below-ground sapwood maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8,hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_crootm_scpf )

    call this%set_history_var(vname='AR_FROOTM_BY_SCPF',units = 'kgC/m2/yr',          &
          long='fine root maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_frootm_scpf )

    ! size-class only variables

    call this%set_history_var(vname='DDBH_CANOPY_BY_SZ', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_canopy_sz )

    call this%set_history_var(vname='DDBH_USTORY_BY_SZ', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_understory_sz )

    call this%set_history_var(vname='CANLEV_YEST_CANOPY_BY_SZ', units = 'indiv/ha',               &
          long='Yesterdays canopy level for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_yesterdaycanopylevel_canopy_sz )

    call this%set_history_var(vname='CANLEV_YEST_USTORY_BY_SZ', units = 'indiv/ha',               &
          long='Yesterdays canopy level for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_yesterdaycanopylevel_understory_sz )

    call this%set_history_var(vname='BA_BY_SZ', units = 'm2/ha',               &
          long='basal area by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ba_sz )

    call this%set_history_var(vname='AGB_BY_SZ', units = 'kgC/m2',               &
          long='Aboveground biomass by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_agb_sz )

    call this%set_history_var(vname='BIOMASS_BY_SZ', units = 'kgC/m2',               &
          long='Total biomass by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_biomass_sz )

    call this%set_history_var(vname='DEMOTION_RATE_BY_SZ', units = 'indiv/ha/yr',               &
          long='demotion rate from canopy to understory by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_demotion_rate_sz )

    call this%set_history_var(vname='PROMOTION_RATE_BY_SZ', units = 'indiv/ha/yr',               &
          long='promotion rate from understory to canopy by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_promotion_rate_sz )

    call this%set_history_var(vname='NPLANT_CANOPY_BY_SZ', units = 'indiv/ha',               &
          long='number of canopy plants by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_canopy_sz )
  
    call this%set_history_var(vname='LAI_CANOPY_BY_SZ', units = 'm2/m2',               &
          long='Leaf are index (LAI) by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_lai_canopy_sz )

    call this%set_history_var(vname='SAI_CANOPY_BY_SZ', units = 'm2/m2',               &
          long='stem area index(SAI) by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_sai_canopy_sz )

    call this%set_history_var(vname='MORTALITY_CANOPY_BY_SZ', units = 'indiv/ha/yr',               &
          long='total mortality of canopy trees by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_canopy_sz )

    call this%set_history_var(vname='NPLANT_USTORY_BY_SZ', units = 'indiv/ha',               &
          long='number of understory plants by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_understory_sz )

    call this%set_history_var(vname='LAI_USTORY_BY_SZ', units = 'indiv/ha',               &
          long='number of understory plants by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_lai_understory_sz )

    call this%set_history_var(vname='SAI_USTORY_BY_SZ', units = 'indiv/ha',               &
          long='number of understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_sai_understory_sz )

    call this%set_history_var(vname='NPLANT_BY_SZ', units = 'indiv/ha',               &
          long='number of plants by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_sz )

    call this%set_history_var(vname='NPLANT_BY_CACLS', units = 'indiv/ha',          &
         long='number of plants by coage class', use_default='active',   &
         avgflag='A', vtype=site_coage_r8, hlms='CLM:ELM', flushval=0.0_r8,     &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_cacls )

    call this%set_history_var(vname='M1_BY_SZ', units = 'N/ha/yr',          &
          long='background mortality by size', use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m1_sz )
    
    call this%set_history_var(vname='M2_BY_SZ', units = 'N/ha/yr',          &
          long='hydraulic mortality by size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m2_sz )

    call this%set_history_var(vname='M3_BY_SZ', units = 'N/ha/yr',          &
          long='carbon starvation mortality by size', use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m3_sz )

    call this%set_history_var(vname='M4_BY_SZ', units = 'N/ha/yr',          &
          long='impact mortality by size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m4_sz )

    call this%set_history_var(vname='M5_BY_SZ', units = 'N/ha/yr',          &
          long='fire mortality by size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m5_sz )

    call this%set_history_var(vname='M6_BY_SZ', units = 'N/ha/yr',          &
          long='termination mortality by size',use_default='active', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m6_sz )

    call this%set_history_var(vname='M7_BY_SZ', units = 'N/ha/event',               &
          long='logging mortality by size',use_default='active',           &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m7_sz )

    call this%set_history_var(vname='M8_BY_SZ', units = 'N/ha/event',               &
          long='freezing mortality by size',use_default='active',           &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m8_sz )

    call this%set_history_var(vname='M9_BY_SZ', units = 'N/ha/yr',              &
          long='senescence mortality by size',use_default='active',         &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m9_sz )

    call this%set_history_var(vname='M10_BY_SZ', units = 'N/ha/yr',              &
          long='age senescence mortality by size',use_default='active',         &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,     &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m10_sz ) 

    call this%set_history_var(vname='M10_BY_CACLS', units = 'N/ha/yr',             &
          long='age senescence mortality by cohort age',use_default='active',      &
          avgflag='A', vtype=site_coage_r8, hlms='CLM:ELM', flushval=0.0_r8,     &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m10_cacls )

    call this%set_history_var(vname='CBALANCE_CANOPY_BY_SZ', units = 'kg C / ha / yr', &
          long='carbon balance for canopy plants by size class', use_default='inactive',    &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_carbon_balance_canopy_sz )

    call this%set_history_var(vname='CBALANCE_USTORY_BY_SZ', units = 'kg C / ha / yr', &
          long='carbon balance for understory plants by size class', use_default='inactive',    &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_carbon_balance_understory_sz )
    
    call this%set_history_var(vname='MORTALITY_USTORY_BY_SZ', units = 'indiv/ha/yr',               &
          long='total mortality of understory trees by size class', use_default='active',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_understory_sz )

    call this%set_history_var(vname='TRIMMING_CANOPY_BY_SZ', units = 'indiv/ha',               &
          long='trimming term of canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_trimming_canopy_sz )

    call this%set_history_var(vname='TRIMMING_USTORY_BY_SZ', units = 'indiv/ha',               &
          long='trimming term of understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_trimming_understory_sz )

    call this%set_history_var(vname='CROWN_AREA_CANOPY_BY_SZ', units = 'm2/ha',               &
          long='total crown area of canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_crown_area_canopy_sz )

    call this%set_history_var(vname='CROWN_AREA_USTORY_BY_SZ', units = 'm2/ha',               &
          long='total crown area of understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_crown_area_understory_sz )

    call this%set_history_var(vname='LEAF_MD_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='LEAF_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_leaf_md_canopy_sz )
    
    call this%set_history_var(vname='ROOT_MD_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='ROOT_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_root_md_canopy_sz )

    call this%set_history_var(vname='BSTORE_MD_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='BSTORE_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bstore_md_canopy_sz )

    call this%set_history_var(vname='BDEAD_MD_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='BDEAD_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bdead_md_canopy_sz )

    call this%set_history_var(vname='BSW_MD_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='BSW_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bsw_md_canopy_sz )

    call this%set_history_var(vname='SEED_PROD_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='SEED_PROD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_seed_prod_canopy_sz )
    
   call this%set_history_var(vname='NPP_LEAF_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_LEAF for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=-999.9_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_leaf_canopy_sz )
    
   call this%set_history_var(vname='NPP_FROOT_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_FROOT for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_fnrt_canopy_sz )
    
   call this%set_history_var(vname='NPP_BSW_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_BSW for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_sapw_canopy_sz )
    
   call this%set_history_var(vname='NPP_BDEAD_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_BDEAD for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_dead_canopy_sz )
    
   call this%set_history_var(vname='NPP_BSEED_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_BSEED for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_seed_canopy_sz )
    
   call this%set_history_var(vname='NPP_STORE_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_STORE for canopy plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_stor_canopy_sz )
    
    call this%set_history_var(vname='LEAF_MR', units = 'kg C / m2 / yr',               &
          long='RDARK (leaf maintenance respiration)', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_leaf_mr )
    
    call this%set_history_var(vname='FROOT_MR', units = 'kg C / m2 / yr',               &
          long='fine root maintenance respiration)', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_froot_mr )
    
    call this%set_history_var(vname='LIVECROOT_MR', units = 'kg C / m2 / yr',               &
          long='live coarse root maintenance respiration)', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livecroot_mr )
    
    call this%set_history_var(vname='LIVESTEM_MR', units = 'kg C / m2 / yr',               &
          long='live stem maintenance respiration)', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livestem_mr )
    
    call this%set_history_var(vname='RDARK_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='RDARK for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_rdark_canopy_sz )
    
    call this%set_history_var(vname='LIVESTEM_MR_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='LIVESTEM_MR for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livestem_mr_canopy_sz )
    
    call this%set_history_var(vname='LIVECROOT_MR_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='LIVECROOT_MR for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livecroot_mr_canopy_sz )
    
    call this%set_history_var(vname='FROOT_MR_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='FROOT_MR for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_froot_mr_canopy_sz )
    
    call this%set_history_var(vname='RESP_G_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='RESP_G for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_g_canopy_sz )
    
    call this%set_history_var(vname='RESP_M_CANOPY_BY_SZ', units = 'kg C / ha / yr',               &
          long='RESP_M for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_m_canopy_sz )

    call this%set_history_var(vname='LEAF_MD_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='LEAF_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_leaf_md_understory_sz )
    
    call this%set_history_var(vname='ROOT_MD_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='ROOT_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_root_md_understory_sz )

    call this%set_history_var(vname='BSTORE_MD_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='BSTORE_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bstore_md_understory_sz )
    
    call this%set_history_var(vname='BDEAD_MD_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='BDEAD_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bdead_md_understory_sz )

    call this%set_history_var(vname='BSW_MD_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='BSW_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bsw_md_understory_sz )
    
    call this%set_history_var(vname='SEED_PROD_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='SEED_PROD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_seed_prod_understory_sz )
    
   call this%set_history_var(vname='NPP_LEAF_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_LEAF for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_leaf_understory_sz )
    
   call this%set_history_var(vname='NPP_FROOT_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_FROOT for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_fnrt_understory_sz )
    
   call this%set_history_var(vname='NPP_BSW_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_BSW for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_sapw_understory_sz )
    
   call this%set_history_var(vname='NPP_BDEAD_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_BDEAD for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_dead_understory_sz )
    
   call this%set_history_var(vname='NPP_BSEED_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_BSEED for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_seed_understory_sz )
    
   call this%set_history_var(vname='NPP_STORE_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
         long='NPP_STORE for understory plants by size class', use_default='inactive',   &
         avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_stor_understory_sz )
    
    call this%set_history_var(vname='RDARK_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='RDARK for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_rdark_understory_sz )
    
    call this%set_history_var(vname='STEM_MR_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='LIVESTEM_MR for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livestem_mr_understory_sz )
    
    call this%set_history_var(vname='CROOT_MR_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='LIVECROOT_MR for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livecroot_mr_understory_sz )
    
    call this%set_history_var(vname='FROOT_MR_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='FROOT_MR for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_froot_mr_understory_sz )
    
    call this%set_history_var(vname='RESP_G_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='RESP_G for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_g_understory_sz )
    
    call this%set_history_var(vname='RESP_M_USTORY_BY_SZ', units = 'kg C / ha / yr',               &
          long='RESP_M for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_m_understory_sz )


    ! CARBON BALANCE VARIABLES THAT DEPEND ON HLM BGC INPUTS

    call this%set_history_var(vname='NEP', units='gC/m^2/s', &
          long='net ecosystem production', use_default='active', &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_nep )

    call this%set_history_var(vname='Fire_CLOSS', units='gC/m^2/s', &
          long='ED/SPitfire Carbon loss to atmosphere', use_default='active', &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_fire_c_to_atm )
   
    call this%set_history_var(vname='FIRE_FLUX_BY_EL', units='g/m^2/s', &
          long='ED-spitfire loss to atmosphere of elements', use_default='active', &
          avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_burn_flux_elem )
   
    call this%set_history_var(vname='CBALANCE_ERROR_FATES', units='mgC/day',  &
         long='total carbon error, FATES', use_default='active', &
         avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cbal_err_fates )

    call this%set_history_var(vname='TOTALC_ERROR', units='mgC/day',  &
         long='total error, FATES mass-balance', use_default='active', &
         avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_err_fates )

    call this%set_history_var(vname='LITTER_FINES_AG_BY_EL', units='kg ha-1', &
          long='mass of above ground  litter in fines (leaves,nonviable seed)', use_default='active', &
          avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_fines_ag_elem )

    call this%set_history_var(vname='LITTER_FINES_BG_BY_EL', units='kg ha-1', &
          long='mass of below ground litter in fines (fineroots)', use_default='active', &
          avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_fines_bg_elem )

    call this%set_history_var(vname='LITTER_CWD_BG_BY_EL', units='kg ha-1', &
          long='mass of below ground litter in CWD (coarse roots)', use_default='active', &
          avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_bg_elem )

    call this%set_history_var(vname='LITTER_CWD_AG_BY_EL', units='kg ha-1', &
          long='mass of above ground litter in CWD (trunks/branches/twigs)', use_default='active', &
          avgflag='A', vtype=site_elem_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_ag_elem )

    call this%set_history_var(vname='LITTER_BY_ELCWD', units='kg ha-1', &
          long='total mass of litter in CWD', use_default='active', &
          avgflag='A', vtype=site_elcwd_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_elcwd )

    ! Mass states C/N/P SCPF dimensions
    ! CARBON
    call this%set_history_var(vname='TOTVEGC_BY_SCPF', units='kgC/ha', &
         long='total vegetation carbon mass in live plants by size-class x pft', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_totvegc_scpf )
    
    call this%set_history_var(vname='LEAFC_BY_SCPF', units='kgC/ha', &
         long='leaf carbon mass by size-class x pft', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_leafc_scpf )

    call this%set_history_var(vname='FNRTC_BY_SCPF', units='kgC/ha', &
         long='fine-root carbon mass by size-class x pft', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_fnrtc_scpf )
    
    call this%set_history_var(vname='SAPWC_BY_SCPF', units='kgC/ha', &
         long='sapwood carbon mass by size-class x pft', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_sapwc_scpf )
    
    call this%set_history_var(vname='STOREC_BY_SCPF', units='kgC/ha', &
         long='storage carbon mass by size-class x pft', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_storec_scpf )
    
    call this%set_history_var(vname='REPROC_BY_SCPF', units='kgC/ha', &
         long='reproductive carbon mass (on plant) by size-class x pft', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_reproc_scpf )
    
    call this%set_history_var(vname='CEFFLUX_BY_SCPF', units='kg/ha/day', &
         long='carbon efflux, root to soil, by size-class x pft', use_default='inactive', &
         avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
         upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cefflux_scpf )

    ! NITROGEN
    nitrogen_active_if2: if(any(element_list(:)==nitrogen_element)) then
       call this%set_history_var(vname='TOTVEGN_BY_SCPF', units='kgN/ha', &
            long='total (live) vegetation nitrogen mass by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_totvegn_scpf )

       call this%set_history_var(vname='LEAFN_BY_SCPF', units='kgN/ha', &
            long='leaf nitrogen mass by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_leafn_scpf )

       call this%set_history_var(vname='FNRTN_BY_SCPF', units='kgN/ha', &
            long='fine-root nitrogen mass by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_fnrtn_scpf )

       call this%set_history_var(vname='SAPWN_BY_SCPF', units='kgN/ha', &
            long='sapwood nitrogen mass by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_sapwn_scpf )

       call this%set_history_var(vname='STOREN_BY_SCPF', units='kgN/ha', &
            long='storage nitrogen mass by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_storen_scpf )

       call this%set_history_var(vname='REPRON_BY_SCPF', units='kgN/ha', &
            long='reproductive nitrogen mass (on plant) by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_repron_scpf )

       call this%set_history_var(vname='NUPTAKE_BY_SCPF', units='kgN d-1 ha-1', &
            long='nitrogen uptake, soil to root, by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nuptake_scpf )

       call this%set_history_var(vname='NEFFLUX_BY_SCPF', units='kgN d-1 ha-1', &
            long='nitrogen efflux, root to soil, by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nefflux_scpf )

       call this%set_history_var(vname='NNEEDGROW_BY_SCPF', units='kgN d-1 ha-1', &
            long='nitrogen needed to match growth, by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nneedgrow_scpf )

       call this%set_history_var(vname='NNEEDMAX_BY_SCPF', units='kgN d-1 ha-1', &
            long='nitrogen needed to reach max concentrations, by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nneedmax_scpf )
       
       
    end if nitrogen_active_if2

    ! PHOSPHORUS
    phosphorus_active_if2: if(any(element_list(:)==phosphorus_element))then
       call this%set_history_var(vname='TOTVEGP_BY_SCPF', units='kgP/ha', &
            long='total (live) vegetation phosphorus mass by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_totvegp_scpf )

       call this%set_history_var(vname='LEAFP_BY_SCPF', units='kgP/ha', &
            long='leaf phosphorus mass by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_leafp_scpf )

       call this%set_history_var(vname='FNRTP_BY_SCPF', units='kgP/ha', &
            long='fine-root phosphorus mass by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_fnrtp_scpf )

       call this%set_history_var(vname='SAPWP_BY_SCPF', units='kgP/ha', &
            long='sapwood phosphorus mass by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_sapwp_scpf )

       call this%set_history_var(vname='STOREP_BY_SCPF', units='kgP/ha', &
            long='storage phosphorus mass by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_storep_scpf )

       call this%set_history_var(vname='REPROP_BY_SCPF', units='kgP/ha', &
            long='reproductive phosphorus mass (on plant) by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_reprop_scpf )

       call this%set_history_var(vname='PUPTAKE_BY_SCPF', units='kg/ha/day', &
            long='phosphorus uptake, soil to root, by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_puptake_scpf )

       call this%set_history_var(vname='PEFFLUX_BY_SCPF', units='kg/ha/day', &
            long='phosphorus efflux, root to soil, by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_pefflux_scpf )

       call this%set_history_var(vname='PNEEDGROW_BY_SCPF', units='kg/ha/day', &
            long='phosphorus needed to match growth, by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_pneedgrow_scpf )

       call this%set_history_var(vname='PNEEDMAX_BY_SCPF', units='kg/ha/day', &
            long='phosphorus needed to reach max concentrations, by size-class x pft', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_pneedmax_scpf )
       
    end if phosphorus_active_if2

    ! organ-partitioned NPP / allocation fluxes
    call this%set_history_var(vname='NPP_LEAF', units='kgC/m2/yr',       &
          long='NPP flux into leaves', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_leaf )

    call this%set_history_var(vname='NPP_SEED', units='kgC/m2/yr',       &
          long='NPP flux into seeds', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_seed )

    call this%set_history_var(vname='NPP_STEM', units='kgC/m2/yr',       &
          long='NPP flux into stem', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_stem )

    call this%set_history_var(vname='NPP_FROOT', units='kgC/m2/yr',       &
          long='NPP flux into fine roots', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_froot )

    call this%set_history_var(vname='NPP_CROOT', units='kgC/m2/yr',       &
          long='NPP flux into coarse roots', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_croot )

    call this%set_history_var(vname='NPP_STOR', units='kgC/m2/yr',       &
          long='NPP flux into storage tissues', use_default='active',               &
          avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_stor )


    ! PLANT HYDRAULICS

    hydro_active_if: if(hlm_use_planthydro.eq.itrue) then
       
       call this%set_history_var(vname='ERRH2O_BY_SCPF', units='kg/indiv/s', &
             long='mean individual water balance error', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_errh2o_scpf )

       call this%set_history_var(vname='TRAN_BY_SCPF', units='kg/indiv/s', &
             long='mean individual transpiration rate', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_tran_scpf )

       call this%set_history_var(vname='SAPFLOW_BY_SCPF', units='kg/ha/s', &
             long='areal sap flow rate dimensioned by size x pft', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_sapflow_scpf )

       call this%set_history_var(vname='SAPFLOW', units='kg/ha/s', &
             long='areal sap flow rate', use_default='active', &
             avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_sapflow )

       
       call this%set_history_var(vname='ITERH1_BY_SCPF', units='count/indiv/step', &
             long='water balance error iteration diagnostic 1', &
             use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_iterh1_scpf )
       
       call this%set_history_var(vname='ITERH2_BY_SCPF', units='count/indiv/step', &
             long='water balance error iteration diagnostic 2', &
             use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_iterh2_scpf )
       
       call this%set_history_var(vname='ATH_BY_SCPF', units='m3 m-3', &
             long='absorbing root water content', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_ath_scpf )
       
       call this%set_history_var(vname='TTH_BY_SCPF', units='m3 m-3', &
             long='transporting root water content', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index =  ih_tth_scpf )
       
       call this%set_history_var(vname='STH_BY_SCPF', units='m3 m-3', &
             long='stem water contenet', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_sth_scpf )
       
       call this%set_history_var(vname='LTH_BY_SCPF', units='m3 m-3', &
             long='leaf water content', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_lth_scpf )

       call this%set_history_var(vname='AWP_BY_SCPF', units='MPa', &
             long='absorbing root water potential', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_awp_scpf )
       
       call this%set_history_var(vname='TWP_BY_SCPF', units='MPa', &
             long='transporting root water potential', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_twp_scpf )
       
       call this%set_history_var(vname='SWP_BY_SCPF', units='MPa', &
             long='stem water potential', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_swp_scpf )
       
       call this%set_history_var(vname='LWP_BY_SCPF', units='MPa', &
             long='leaf water potential', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_lwp_scpf )
 
       call this%set_history_var(vname='AFLC_BY_SCPF', units='fraction', &
             long='absorbing root fraction of condutivity', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_aflc_scpf )
       
       call this%set_history_var(vname='TFLC_BY_SCPF', units='fraction', &
             long='transporting root fraction of condutivity', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_tflc_scpf )
       
       call this%set_history_var(vname='SFLC_BY_SCPF', units='fraction', &
             long='stem water fraction of condutivity', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_sflc_scpf )
       
       call this%set_history_var(vname='LFLC_BY_SCPF', units='fraction', &
             long='leaf fraction of condutivity', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_lflc_scpf )
       
       call this%set_history_var(vname='BTRAN_BY_SCPF', units='unitless', &
             long='mean individual level btran', use_default='inactive', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_btran_scpf )
       
       call this%set_history_var(vname='ROOTWGT_SOILVWC', units='m3 m-3', &
            long='soil volumetric water content, weighted by root area', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootwgt_soilvwc )

       call this%set_history_var(vname='ROOTWGT_SOILVWCSAT', units='m3 m-3', &
            long='soil saturated volumetric water content, weighted by root area', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootwgt_soilvwcsat )
       
       call this%set_history_var(vname='ROOTWGT_SOILMATPOT', units='MPa', &
            long='soil matric potential, weighted by root area', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootwgt_soilmatpot )
       
       call this%set_history_var(vname='SOILMATPOT_BY_SL', units='MPa', &
            long='soil water matric potenial by soil layer', use_default='inactive', &
            avgflag='A', vtype=site_ground_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_soilmatpot_sl )
       
       call this%set_history_var(vname='SOILVWC_BY_SL', units='m3 m-3', &
            long='soil volumetric water content by soil layer', use_default='inactive', &
            avgflag='A', vtype=site_ground_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_soilvwc_sl )
       
       call this%set_history_var(vname='SOILVWCSAT_BY_SL', units='m3 m-3', &
            long='soil saturated volumetric water content by soil layer', use_default='inactive', &
            avgflag='A', vtype=site_ground_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_soilvwcsat_sl )
       
       call this%set_history_var(vname='ROOTUPTAKE', units='kg ha-1 s-1', &
            long='root water uptake rate', use_default='active', &
            avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake )
       
       call this%set_history_var(vname='ROOTUPTAKE_BY_SL', units='kg ha-1 s-1', &
            long='root water uptake rate by soil layer', use_default='inactive', &
            avgflag='A', vtype=site_ground_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake_sl )
       
       call this%set_history_var(vname='ROOTUPTAKE0_BY_SCPF', units='kg ha-1 m-1 s-1', &
            long='root water uptake from 0 to to 10 cm depth, by plant size x pft ', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake0_scpf )
       
       call this%set_history_var(vname='ROOTUPTAKE10_BY_SCPF', units='kg ha-1 m-1 s-1', &
            long='root water uptake from 10 to to 50 cm depth, by plant size x pft ', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake10_scpf )

       call this%set_history_var(vname='ROOTUPTAKE50_BY_SCPF', units='kg ha-1 m-1 s-1', &
            long='root water uptake from 50 to to 100 cm depth, by plant size x pft ', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake50_scpf )

       call this%set_history_var(vname='ROOTUPTAKE100_BY_SCPF', units='kg ha-1 m-1 s-1', &
            long='root water uptake below 100 cm depth, by plant size x pft ', use_default='inactive', &
            avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ELM', flushval=hlm_hio_ignore_val,    &
            upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake100_scpf )

       call this%set_history_var(vname='H2OVEG', units = 'kg/m2',               &
             long='water stored inside vegetation tissues (leaf, stem, roots)', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg )

       call this%set_history_var(vname='H2OVEG_DEAD', units = 'kg/m2',               &
             long='cumulative plant_stored_h2o in dead biomass due to mortality', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg_dead )

       call this%set_history_var(vname='H2OVEG_RECRUIT', units = 'kg/m2',               &
             long='amount of water in new recruits', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg_recruit )
    
       call this%set_history_var(vname='H2OVEG_GROWTURN_ERR', units = 'kg/m2',               &
             long='cumulative net borrowed (+) or lost (-) from plant_stored_h2o due to combined growth & turnover', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg_growturn_err )
    
       call this%set_history_var(vname='H2OVEG_PHENO_ERR', units = 'kg/m2',               &
             long='cumulative net borrowed (+) from plant_stored_h2o due to leaf emergence', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg_pheno_err )
     
       call this%set_history_var(vname='H2OVEG_HYDRO_ERR', units = 'kg/m2',               &
             long='cumulative net borrowed (+) from plant_stored_h2o due to plant hydrodynamics', use_default='inactive',   &
             avgflag='A', vtype=site_r8, hlms='CLM:ELM', flushval=0.0_r8,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_h2oveg_hydro_err )
    end if hydro_active_if 

    ! Must be last thing before return
    this%num_history_vars_ = ivar
    
  end subroutine define_history_vars


   ! ====================================================================================
   ! DEPRECATED, TRANSITIONAL OR FUTURE CODE SECTION
   ! ====================================================================================

   !subroutine set_fates_hio_str(tag,iotype_name, iostr_val)

!       ! Arguments
!       character(len=*), intent(in)           :: tag
!       character(len=*), optional, intent(in) :: iotype_name
!       integer, optional, intent(in)         :: iostr_val

!       ! local variables
!       logical              :: all_set
!       integer,  parameter  :: unset_int = -999
!       real(r8), parameter  :: unset_double = -999.9
!       integer              :: ityp, idim

!       select case (trim(tag))
!       case('flush_to_unset')
!          write(*, *) ''
!          write(*, *) 'Flushing FATES IO types prior to transfer from host'
!          do ityp=1,ubound(iovar_str, 1)
!             iovar_str(ityp)%dimsize = unset_int
!             iovar_str(ityp)%active  = .false.
!          end do

!       case('check_allset')
!          do ityp=1,ubound(iovar_str, 1)
!             write(*, *) 'Checking to see if ',iovar_str(ityp)%name, ' IO communicators were sent to FATES'
!             if(iovar_str(ityp)%active)then
!                if(iovar_str(ityp)%offset .eq. unset_int) then
!                   write(*, *) 'FATES offset information of IO type:', iovar_str(ityp)%name
!                   write(*, *) 'was never set'
!                   ! end_run('MESSAGE')
!                end if
!                do idim=1, iovar_str(ityp)%ndims
!                   if(iovar_str(ityp)%dimsize(idim) .eq. unset_int) then
!                      write(*, *) 'FATES dimension information of IO type:', iovar_str(ityp)%name
!                      write(*, *) 'was never set'
!                      ! end_run('MESSAGE')
!                   end if
!                end do
!             end if
!          end do
!          write(*, *) 'Checked. All history IO specifications properly sent to FATES.'
!       case default

!          ! Must have two arguments if this is not a check or flush
!          if(present(iostr_val) .and. present(iotype_name))then
!
!             ! Tag in this case is dimsize or offset
!             select case (trim(tag))
!
!             case('offset')
!                ityp=iotype_index(trim(iotype_name))
!                iovar_str(ityp)%offset = iostr_val
!                write(*, *) 'Transfering offset for IOTYPE',iotype_name, ' to FATES'

!             case('dimsize1')
!                ityp=iotype_index(trim(iotype_name))
!                iovar_str(ityp)%dimsize(1) = iostr_val
!                write(*, *) 'Transfering 1st dimension size for IOTYPE',iotype_name, ' to FATES'

!             case('dimsize2')
!                ityp=iotype_index(trim(iotype_name))
!                if(ubound(iovar_str(ityp)%dimsize, 1)==1)then
!                   write(fates_log(), *) 'Transfering second dimensional bound to unallocated space'
!                   write(fates_log(), *) 'type:', iotype_name
!                   ! end_run
!                end if
!                iovar_str(ityp)%dimsize(2) = iostr_val
!                write(*, *) 'Transfering 2nd dimension size for IOTYPE',iotype_name, ' to FATES'

!             case('dimsize3')
!                ityp=iotype_index(trim(iotype_name))
!                if(ubound(iovar_str(ityp)%dimsize, 1)<3)then
!                   write(fates_log(), *) 'Transfering third dimensional bound to unallocated space'
!                   write(fates_log(), *) 'type:', iotype_name
!                   ! end_run
!                end if
!                iovar_str(ityp)%dimsize(3) = iostr_val
!                write(*, *) 'Transfering 3rd dimension size for IOTYPE',iotype_name, ' to FATES'

!             case default
!                write(*, *) 'IO parameter not recognized:', trim(tag)
!                ! end_run
!             end select
!          else
!             write(*, *) 'no value was provided for the tag'
!          end if
!
!       end select
!       return
!     end subroutine set_fates_hio_str



end module FatesHistoryInterfaceMod
