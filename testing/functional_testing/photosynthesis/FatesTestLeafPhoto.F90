program FatesTestLeafPhoto
  
  use FatesConstantsMod,           only : r8 => fates_r8
  use FATESPlantRespPhotosynthMod, only : GetCanopyGasParameters
  use FATESPlantRespPhotosynthMod, only : LeafLayerBiophysicalRates
  use FATESPlantRespPhotosynthMod, only : LeafLayerPhotosynthesis
  use FatesTestPhotosynthesisMod,  only : CalcVaporPressure
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use FatesParameterDerivedMod,    only : param_derived
  use PRTParametersMod,            only : prt_params
  
  implicit none
  
  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader              ! param reader instance
  character(len=:), allocatable      :: param_file                ! input parameter file
  real(r8),         allocatable      :: veg_tempk(:)              ! vegetation temperature [K]
  real(r8),         allocatable      :: mm_kco2(:)                ! Michaelis-Menten constant for CO2 [Pa]
  real(r8),         allocatable      :: mm_ko2(:)                 ! Michaelis-Menten constant for O2 [Pa]
  real(r8),         allocatable      :: co2_compensation_pt(:)    ! CO2 compensation point [Pa]
  real(r8),         allocatable      :: cf(:)                     ! conversion factor between molar form and velocity form of conductance and resistance: [umol/m3]
  real(r8),         allocatable      :: leaf_bl_conductance(:)    ! leaf boundary layer conductance [umol H2O/m2/s]
  real(r8),         allocatable      :: constrained_air_vpress(:) ! vapor pressure of air, constrained [Pa]
  real(r8),         allocatable      :: vcmax(:,:)                ! maximum rate of carboxylation [umol CO2/m2/s]
  real(r8),         allocatable      :: jmax(:,:)                 ! maximum electron transport rate [umol electrons/m2/s]
  real(r8),         allocatable      :: co2_rcurve_islope(:,:)    ! initial slope of CO2 response curve (C4 plants)
  
  real(r8)                           :: can_air_vpress            ! canopy air vapor pressure [Pa]
  real(r8)                           :: veg_esat                  ! saturation vapor pressure at vegetation surface [Pa]
  integer                            :: num_temp                  ! array size
  integer                            :: numpft                    ! number of PFTs (from parameter file)
  integer                            :: i                         ! looping index
  integer                            :: ft                        ! plant functional type index
  integer                            :: nargs                     ! number of command-line arguments
  integer                            :: arglen                    ! command-line argument length
  
  ! real(r8)                           :: frac_sun             ! fraction sun / shaded leaf area [0-1]
  ! real(r8)                           :: par_sun              ! absorbed PAR in sunlit leaves []
  ! real(r8)                           :: par_sha              ! absorved PAR in shaded leaves

  
  ! real(r8)                           :: veg_tempk            ! vegetation temperature
  ! real(r8)                           :: veg_esat             ! saturation vapor pressure at veg_tempk (Pa)
  ! real(r8)                           :: can_press            ! air pressure NEAR the surface of the leaf (Pa)
  ! real(r8)                           :: can_co2_ppress       ! partial pressure of CO2 NEAR the leaf surface (Pa)
  ! real(r8)                           :: can_o2_ppress        ! partial pressure of O2 NEAR the leaf surface (Pa)
  ! real(r8)                           :: btran                ! transpiration wetness factor (0 to 1)
  ! real(r8)                           :: stom_int_btran       ! water-stressed minimum stomatal conductance (umol H2O/m**2/s)
  ! real(r8)                           :: cf                   ! s m**2/umol -> s/m (ideal gas conversion) [umol/m3]
  ! real(r8)                           :: gb_mol               ! leaf boundary layer conductance (umol /m**2/s)
  ! real(r8)                           :: ceair                ! vapor pressure of airconstrained (Pa)
  ! real(r8)                           :: lmr                  ! Leaf Maintenance Respiration  (umol CO2/m**2/s)
  ! real(r8)                           :: leaf_psi             ! Leaf water potential [MPa]
  ! real(r8)                           :: rb                   ! Boundary Layer resistance of leaf [s/m]
  ! real(r8)                           :: psn                  ! carbon assimilated in this leaf layer umolC/m2/s
  ! real(r8)                           :: rstoma               ! stomatal resistance (1/gs_lsl) (s/m)
  ! real(r8)                           :: anet_av              ! net leaf photosynthesis (umol CO2/m**2/s)
  ! real(r8)                           :: c13disc              ! carbon 13 in newly assimilated carbon
  ! real(r8)                           :: kn 
  ! real(r8)                           :: cumulative_lai, nscaler
  ! real(r8)                           :: daylength_factor, kp
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'leaf_photo_out.nc' ! output file
  real(r8),         parameter :: default_vpd = 1000.0_r8                   ! default vapor pressure deficit [Pa]
  real(r8),         parameter :: can_air_press = 101325.0_r8               ! default value for canopy air pressure [Pa]
  real(r8),         parameter :: can_o2_partial_press = 20900.0_r8         ! default value for O2 partial pressure [Pa]
  real(r8),         parameter :: default_veg_tempk = 25.0_r8 + 273.15_r8   ! default value for vegetation temperature [K]
  real(r8),         parameter :: default_can_air_tempk = default_veg_tempk ! default value for canopy air temperature [K]
  real(r8),         parameter :: leaf_bl_resistance = 125.0_r8             ! default value for leaf boundary layer resistance [umol H2O/m2/s]
  real(r8),         parameter :: min_temp = 279.15_r8                      ! minimum temperature to calculate [K]
  real(r8),         parameter :: max_temp = 313.15_r8                      ! maximum temperature to calculate [K]
  real(r8),         parameter :: temp_inc = 0.5_r8                         ! temperature increment to use [K]
  real(r8),         parameter :: default_nscaler = 1.0_r8                  ! default scaler for leaf nitrogen [0-1]
  real(r8),         parameter :: default_dayl_fact = 1.0_r8                ! default scaler for day length [0-1]
  real(r8),         parameter :: default_btran = 1.0_r8                    ! default scaler for BTRAN [0-1]
  
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
  
  ! calculate size of temperature array based on desired min/max/increment
  num_temp = int((max_temp - min_temp)/temp_inc + 1)
  numpft = size(prt_params%wood_density, dim=1)

  ! allocate arrays
  allocate(veg_tempk(num_temp))
  allocate(mm_kco2(num_temp))
  allocate(mm_ko2(num_temp))
  allocate(co2_compensation_pt(num_temp))
  allocate(cf(num_temp))
  allocate(leaf_bl_conductance(num_temp))
  allocate(constrained_air_vpress(num_temp))
  allocate(vcmax(num_temp, numpft))
  allocate(jmax(num_temp, numpft))
  allocate(co2_rcurve_islope(num_temp, numpft))
  
  ! initialize temperature array
  do i = 1, num_temp
    veg_tempk(i) = min_temp + temp_inc*(i-1)
  end do
  
  ! calculate photosynthesis as we scale temperature and hold everything else constant
  do i = 1, num_temp
    
    ! calculate saturation vapor pressure and vapor pressure
    call CalcVaporPressure(veg_tempk(i), default_vpd, veg_esat, can_air_vpress)
    
    ! calculate canopy gas parameters
    call GetCanopyGasParameters(can_air_press, can_o2_partial_press, veg_tempk(i),       &
      veg_tempk(i), can_air_vpress, veg_esat, leaf_bl_resistance, mm_kco2(i), mm_ko2(i), &
      co2_compensation_pt(i), cf(i), leaf_bl_conductance(i), constrained_air_vpress(i))
      
    do ft = 1, numpft 
      
      ! calculate leaf biophysical rates
      call LeafLayerBiophysicalRates(100.0_r8, ft, EDPftvarcon_inst%vcmax25top(ft,1),    &
      param_derived%jmax25top(ft,1), param_derived%kp25top(ft,1), default_nscaler,       &
      veg_tempk(i), default_dayl_fact, veg_tempk(i), veg_tempk(i), default_btran,        &
      vcmax(i,ft), jmax(i,ft), co2_rcurve_islope(i,ft))
      
      ! calculate leaf-level photosynthesis
      call LeafLayerPhotosynthesis()
    end do
      
  end do
  
end program FatesTestLeafPhoto

