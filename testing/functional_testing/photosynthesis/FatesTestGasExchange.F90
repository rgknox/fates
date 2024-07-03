program FatesTestGasExchange
  
  use FatesConstantsMod,           only : r8 => fates_r8
  use FATESPlantRespPhotosynthMod, only : GetCanopyGasParameters
  use FatesTestPhotosynthesisMod,  only : CalcVaporPressure
  
  implicit none
  
  ! LOCALS:
  real(r8)              :: can_air_vpress            ! canopy air vapor pressure [Pa]
  real(r8)              :: veg_esat                  ! saturation vapor pressure at vegetation surface [Pa]
  real(r8), allocatable :: vapor_press_deficit(:)    ! vapor pressure deficit [Pa]
  real(r8), allocatable :: veg_tempk(:)              ! vegetation temperature [K]
  real(r8), allocatable :: can_air_tempk(:)          ! canopy air temperature [K]
  real(r8), allocatable :: mm_kco2(:)                ! Michaelis-Menten constant for CO2 [Pa]
  real(r8), allocatable :: mm_ko2(:)                 ! Michaelis-Menten constant for O2 [Pa]
  real(r8), allocatable :: co2_compensation_pt(:)    ! CO2 compensation point [Pa]
  real(r8), allocatable :: cf(:)                     ! conversion factor between molar form and velocity form of conductance and resistance: [umol/m3]
  real(r8), allocatable :: leaf_bl_conductance(:)    ! leaf boundary layer conductance [umol H2O/m2/s]
  real(r8), allocatable :: constrained_air_vpress(:) ! vapor pressure of air, constrained [Pa]
  integer               :: num_temp                  ! array size
  integer               :: i                         ! looping index
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'gas_exchange_out.nc'             ! output file
  real(r8),         parameter :: default_vpd = 1000.0_r8                      ! default vapor pressure deficit [Pa]
  real(r8),         parameter :: can_air_press = 101325.0_r8                  ! default value for canopy air pressure [Pa]
  real(r8),         parameter :: can_o2_partial_press = 0.29_r8*can_air_press ! default value for O2 partial pressure [Pa]
  real(r8),         parameter :: default_veg_tempk = 25.0_r8 + 273.15_r8      ! default value for vegetation temperature [K]
  real(r8),         parameter :: default_can_air_tempk = default_veg_tempk    ! default value for canopy air temperature [K]
  real(r8),         parameter :: leaf_bl_resistance = 125.0_r8                ! default value for leaf boundary layer resistance [umol H2O/m2/s]
  real(r8),         parameter :: min_temp = 279.15_r8                         ! minimum temperature to calculate [K]
  real(r8),         parameter :: max_temp = 313.15_r8                         ! maximum temperature to calculate [K]
  real(r8),         parameter :: temp_inc = 0.5_r8                            ! temperature increment to use [K]
  
  numtemp = int((max_temp - min_temp)/temp_inc + 1)
  
  ! allocate arrays
  allocate(veg_tempk(numdbh, numpft))
  allocate(can_air_tempk(numdbh, numpft))
  allocate(mm_kco2(numdbh, numpft))
  allocate(mm_ko2(numdbh, numpft))
  allocate(co2_compensation_pt(numdbh, numpft))
  allocate(cf(numdbh, numpft))
  allocate(leaf_bl_conductance(numdbh, numpft))
  allocate(constrained_air_vpress(numdbh, numpft))
  
  ! initialize temperature array
  do i = 1, numtemp
    veg_tempk(i) = min_temp + temp_inc*(i-1)
  end do

  
  call GetCanopyGasParameters(can_air_press, can_o2_partial_press, default_veg_tempk, default_can_air_tempk, default_can_air_vpress, default_veg_esat,  &
    default_leaf_bl_res, mm_kco2, mm_ko2, co2_compensation_pt, cf,                       &
    leaf_bl_conductance, constrained_air_vpress)
      
end program FatesTestGasExchange