program FatesTestLeafPhoto
  
  use FatesConstantsMod,           only : r8 => fates_r8
  use FATESPlantRespPhotosynthMod, only : LeafLayerPhotosynthesis
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  
  implicit none
  
  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader      ! param reader instance
  character(len=:), allocatable      :: param_file        ! input parameter file
  integer                            :: numpft            ! number of pfts (from parameter file)
  integer                            :: arglen            ! length of command line argument
  integer                            :: nargs             ! number of command line arguments
  real(r8)                           :: f_sun_lsl         ! fraction sun / shaded leaf area [0-1]
  real(r8)                           :: parsun_lsl        ! absorbed PAR in sunlit leaves []
  real(r8)                           :: parsha_lsl        ! absorved PAR in shaded leaves
  real(r8)                           :: laisun_lsl        ! LAI in sunlit leaves
  real(r8)                           :: laisha_lsl        ! LAI in shaded leaves
  real(r8)                           :: canopy_area_lsl   ! canopy area
  integer                            :: ft                ! plant functional type index
  real(r8)                           :: vcmax             ! maximum rate of carboxylation (umol co2/m**2/s)
  real(r8)                           :: jmax              ! maximum electron transport rate (umol electrons/m**2/s)
  real(r8)                           :: co2_rcurve_islope ! initial slope of CO2 response curve (C4 plants)
  real(r8)                           :: veg_tempk         ! vegetation temperature
  real(r8)                           :: veg_esat          ! saturation vapor pressure at veg_tempk (Pa)
  real(r8)                           :: can_press         ! air pressure NEAR the surface of the leaf (Pa)
  real(r8)                           :: can_co2_ppress    ! partial pressure of CO2 NEAR the leaf surface (Pa)
  real(r8)                           :: can_o2_ppress     ! partial pressure of O2 NEAR the leaf surface (Pa)
  real(r8)                           :: btran             ! transpiration wetness factor (0 to 1)
  real(r8)                           :: stom_int_btran    ! water-stressed minimum stomatal conductance (umol H2O/m**2/s)
  real(r8)                           :: cf                ! s m**2/umol -> s/m (ideal gas conversion) [umol/m3]
  real(r8)                           :: gb_mol            ! leaf boundary layer conductance (umol /m**2/s)
  real(r8)                           :: ceair             ! vapor pressure of airconstrained (Pa)
  real(r8)                           :: mm_kco2           ! Michaelis-Menten constant for CO2 (Pa)
  real(r8)                           :: mm_ko2            ! Michaelis-Menten constant for O2 (Pa)
  real(r8)                           :: co2_cpoint        ! CO2 compensation point (Pa)
  real(r8)                           :: lmr               ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
  real(r8)                           :: leaf_psi          ! Leaf water potential [MPa]
  real(r8)                           :: rb                ! Boundary Layer resistance of leaf [s/m]
  real(r8)                           :: psn_out           ! carbon assimilated in this leaf layer umolC/m2/s
  real(r8)                           :: rstoma_out        ! stomatal resistance (1/gs_lsl) (s/m)
  real(r8)                           :: anet_av_out       ! net leaf photosynthesis (umol CO2/m**2/s)
  real(r8)                           :: c13disc_z         ! carbon 13 in newly assimilated carbon
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'leaf_photo_out.nc' ! output file
  
  ! get parameter file from command-line argument
  nargs = command_argument_count()
  if (nargs /= 1) then
    write(*, '(ai2a)') "Incorrect number of arguments: ", nargs, ". Should be 1."
    stop
  else
    call get_command_argument(1, length=arglen)
    allocate(character(arglen) :: param_file)
    call get_command_argument(1, value=param_file)
  endif
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  
  print *, "hello photosynthesis"
  
end program FatesTestLeafPhoto