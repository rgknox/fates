module PRTParametersMod


  use FatesConstantsMod,     only : r8 => fates_r8
  
  ! This module only holds the parameter definitions for PARTEH and allometry.
  ! This does not hold any of the code used for intiailizing and filling
  ! that data, for that is model dependent (ie FATES may have a different
  ! way than another TBM)
  ! This code does perform checks on parameters.

  type,public ::  prt_param_type

     ! The following three PFT classes 
     ! are mutually exclusive
     ! MLO: perhaps we should replace these three parameters with a single
     !      parameter (phenology(:)) that is assigned different indices?
     integer, allocatable :: stress_decid(:)        ! Is the plant stress deciduous?
                                                    ! 0 - No
                                                    ! 1 - Drought "hard" deciduous (i.e., PFT
                                                    !     sheds leaves all at once when stressed)
                                                    ! 2 - Drought semi-deciduous (i.e., PFT
                                                    !     sheds leaves gradually as drought
                                                    !     conditions deteriorate)
     integer, allocatable :: season_decid(:)        ! Is the plant seasonally deciduous (1=yes, 0=no)
     integer, allocatable :: evergreen(:)           ! Is the plant an evergreen (1=yes, 0=no)

     ! Drop fraction for tissues other than leaves (PFT-dependent)
     real(r8), allocatable :: phen_fnrt_drop_fraction(:) ! Abscission fraction of fine roots
     real(r8), allocatable :: phen_stem_drop_fraction(:) ! Abscission fraction of stems
     real(r8), allocatable :: phen_drought_threshold(:)   ! For obligate (hard) drought deciduous, this is the threshold
                                                          !    below which plants will abscise leaves, and 
                                                          !    above which plants will flush leaves. For semi-deciduous
                                                          !    plants, this is the threshold below which abscission will
                                                          !    be complete. This depends on the sign. If positive, these
                                                          !    are soil volumetric water content [m3/m3]. If
                                                          !    negative, the values are soil matric potential [mm].
                                                          !    Ignored for non-deciduous plants.
     real(r8), allocatable :: phen_moist_threshold(:)     ! For semi-deciduous, this is the threshold above which flushing 
                                                          !    will be complete.  This depends on the sign. If positive, these
                                                          !    are soil volumetric water content [m3/m3]. If
                                                          !    negative, the values are soil matric potential [mm].
                                                          !    Ignored for non-deciduous plants.
     real(r8), allocatable :: phen_doff_time(:)           ! Minimum number of days that plants must remain leafless before
                                                          !   flushing leaves again.


     ! Growth and Turnover Parameters
     real(r8), allocatable :: senleaf_long_fdrought(:)   ! Multiplication factor for leaf longevity of senescent 
                                                         ! leaves during drought( 1.0 indicates no change)
     real(r8), allocatable :: leaf_long(:,:)             ! Leaf turnover time (longevity) (pft x age-class)
                                                         ! If there is >1 class, it is the longevity from
                                                         ! one class to the next [yr]
                                                         !   For drought-deciduous PFTs, the sum of leaf
                                                         !   longevity across all leaf age classes is also
                                                         !   the maximum length of the growing (i.e., leaves on)
                                                         !   season.
     real(r8), allocatable :: root_long(:)               ! root turnover time (longevity) (pft)             [yr]
     real(r8), allocatable :: branch_long(:)             ! Turnover time for branchfall on live trees (pft) [yr]
     real(r8), allocatable :: turnover_nitr_retrans(:,:) ! nitrogen re-translocation fraction (pft x organ)
     real(r8), allocatable :: turnover_phos_retrans(:,:) ! phosphorus re-translocation fraction (pft x organ)
                                                         ! Parameters dimensioned by PFT and leaf age
     
     real(r8), allocatable :: grperc(:)                  ! Growth respiration per unit Carbon gained
                                                         ! One value for whole plant
     ! ONLY parteh_mode == 1  [kg/kg]
     !     real(r8), allocatable ::grperc_organ(:,:)     ! Unit growth respiration (pft x organ) [kg/kg]
     !                                                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !                                                        ! THIS IS NOT READ IN BY THE PARAMETER FILE
     !                                                        ! THIS IS JUST FILLED BY GRPERC.  WE KEEP THIS
     !                                                        ! PARAMETER FOR HYPOTHESIS TESTING (ADVANCED USE)
     !                                                        ! IT HAS THE PRT_ TAG BECAUSE THIS PARAMETER
     !                                                        ! IS USED INSIDE PARTEH, WHILE GRPERC IS APPLIED
     !                                                        ! IN THE LEAF BIOPHYSICS SCHEME
     !                                                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     real(r8), allocatable :: nitr_stoich_p1(:,:)        ! Parameter 1 for nitrogen stoichiometry (pft x organ) 
     real(r8), allocatable :: phos_stoich_p1(:,:)        ! Parameter 1 for phosphorus stoichiometry (pft x organ) 

     real(r8), allocatable :: nitr_store_ratio(:)        ! This is the ratio of the target nitrogen stored per
                                                         ! target nitrogen that is bound into the tissues
                                                         ! of leaves, fine-roots and sapwood
     
     
     real(r8), allocatable :: phos_store_ratio(:)        ! This is the ratio of the target phosphorus stored per
                                                         ! target phosphorus is bound into the tissues
                                                         ! of leaves, fine-roots and sapwood

     integer, allocatable :: organ_id(:)                 ! Mapping of the organ index in the parameter file, to the
                                                         ! global list of organs found in PRTGenericMod.F90
     real(r8), allocatable :: alloc_priority(:,:)        ! Allocation priority for each organ (pft x organ) [integer 0-6]
     real(r8), allocatable :: cushion(:)                 ! labile carbon storage target as multiple of leaf pool.
     real(r8), allocatable :: leaf_stor_priority(:)      ! leaf turnover vs labile carbon use prioritisation
                                                         ! (1 = lose  leaves, 0 = use store).
     real(r8), allocatable :: dbh_repro_threshold(:)     ! diameter at which mature plants shift allocation
     real(r8), allocatable :: seed_alloc_mature(:)       ! fraction of carbon balance allocated to 
                                                         ! clonal reproduction.
     real(r8), allocatable :: seed_alloc(:)              ! fraction of carbon balance allocated to seeds.
     real(r8), allocatable :: repro_alloc_a(:)           ! ahb added this; sigmoidal shape param relating dbh to seed allocation fraction
     real(r8), allocatable :: repro_alloc_b(:)           ! ahb added this; intercept param relating dbh to seed allocation fraction

     ! Derived parameters

     integer, allocatable :: organ_param_id(:)           ! This is the sparse reverse lookup index map. This is dimensioned
                                                         ! by all the possible organs in parteh, and each index
                                                         ! may point to the index in the parameter file, or will be -1
     



     

     ! PID controller parameters
     real(r8), allocatable :: pid_kp(:)                     ! proportion constant in the PID controller for fine-root biomass
     real(r8), allocatable :: pid_ki(:)                     ! integral constant in the PID controller for fine-root biomass
     real(r8), allocatable :: pid_kd(:)                     ! derivative constant in the PID controller for fine-root biomass
     
     real(r8), allocatable :: store_ovrflw_frac(:)          ! For a coupled nutrient enabled simulation with dynamic fine-root biomass,
                                                            ! there will be an excess of at least two of the three species C, N or P.
                                                            ! This specifies how much excess (overflow) is allowed to be retained in storage
                                                            ! beyond the target level before it is either burned (C) or exuded (N or P). The
                                                            ! maximum value is the target * (1+store_ovrflw_frac)
     

     real(r8), allocatable :: nfix_mresp_scfrac(:)            ! Surcharge (as a fraction) to add to maintentance respiration
                                                            ! that is used to pay for N-Fixation
     
  end type prt_param_type

  type(prt_param_type),public :: prt_params          ! Instantiation of the parameter object


  

  
  
end module PRTParametersMod

