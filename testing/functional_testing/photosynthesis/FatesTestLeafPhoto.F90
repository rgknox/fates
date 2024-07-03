program FatesTestLeafPhoto
  
  use FatesConstantsMod,           only : r8 => fates_r8
  use FATESPlantRespPhotosynthMod, only : LeafLayerPhotosynthesis, LeafLayerBiophysicalRates
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use FatesParameterDerivedMod,    only : param_derived
  
  implicit none
  
  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader      ! param reader instance
  character(len=:), allocatable      :: param_file        ! input parameter file
  integer                            :: numpft            ! number of pfts (from parameter file)
  integer                            :: arglen            ! length of command line argument
  integer                            :: nargs             ! number of command line arguments
  real(r8)                           :: frac_sun          ! fraction sun / shaded leaf area [0-1]
  real(r8)                           :: par_sun           ! absorbed PAR in sunlit leaves []
  real(r8)                           :: par_sha           ! absorved PAR in shaded leaves
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
  real(r8)                           :: psn               ! carbon assimilated in this leaf layer umolC/m2/s
  real(r8)                           :: rstoma            ! stomatal resistance (1/gs_lsl) (s/m)
  real(r8)                           :: anet_av           ! net leaf photosynthesis (umol CO2/m**2/s)
  real(r8)                           :: c13disc           ! carbon 13 in newly assimilated carbon
  real(r8)                           :: kn 
  real(r8)                           :: cumulative_lai, nscaler
  real(r8)                           :: daylength_factor, kp
  
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
  
  ! initialize values
  frac_sun = 0.5_r8
  par_sun = 850.0_r8
  par_sha = 850.0_r8
  ft = 1
  kn = 0.11_r8
  cumulative_lai = 5.0
  nscaler = 1.0 !exp(-kn*cumulative_lai)
  daylength_factor = 1.0_r8
  veg_tempk = 25.0_r8 + 273.15_r8
  veg_esat = 
  btran = 1.0_r8

  ! need this to calculate vcmax, jmax, and kp
  call LeafLayerBiophysicalRates(50.0_r8, ft, EDPftvarcon_inst%vcmax25top(ft,1),         &
    param_derived%jmax25top(ft,1), param_derived%kp25top(ft,1), nscaler, veg_tempk,      &
    daylength_factor, veg_tempk, veg_tempk, btran, vcmax, jmax, kp)

  call LeafLayerPhotosynthesis(frac_sun, par_sun, par_sha, ft, vcmax, jmax,              &
    kp, veg_tempk, veg_esat, can_press, can_co2_ppress,                                  &
    can_o2_ppress, btran, stom_int_btran, cf, gb_mol, ceair, mm_kco2, mm_ko2, co2_cpoint, &  
    lmr, leaf_psi, rb, psn, rstoma, anet_av, c13disc)   
  
end program FatesTestLeafPhoto