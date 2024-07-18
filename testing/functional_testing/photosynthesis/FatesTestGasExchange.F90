program FatesTestGasExchange
  
  use FatesConstantsMod,           only : r8 => fates_r8
  use FATESPlantRespPhotosynthMod, only : GetCanopyGasParameters
  use FatesTestPhotosynthesisMod,  only : CalcVaporPressure
  use FatesConstantsMod,           only : tfrz => t_water_freeze_k_1atm
  
  implicit none
  
  ! LOCALS:
  real(r8)              :: can_air_vpress            ! canopy air vapor pressure [Pa]
  real(r8)              :: veg_esat                  ! saturation vapor pressure at vegetation surface [Pa]
  real(r8), allocatable :: veg_tempk(:)              ! vegetation temperature [K]
  real(r8), allocatable :: mm_kco2(:)                ! Michaelis-Menten constant for CO2 [Pa]
  real(r8), allocatable :: mm_ko2(:)                 ! Michaelis-Menten constant for O2 [Pa]
  real(r8), allocatable :: co2_compensation_pt(:)    ! CO2 compensation point [Pa]
  real(r8), allocatable :: cf(:)                     ! conversion factor between molar form and velocity form of conductance and resistance: [umol/m3]
  real(r8), allocatable :: leaf_bl_conductance(:)    ! leaf boundary layer conductance [umol H2O/m2/s]
  real(r8), allocatable :: constrained_air_vpress(:) ! vapor pressure of air, constrained [Pa]
  real(r8)              :: veg_air_vpress            ! vegetation air vapor pressure [Pa]
  integer               :: num_temp                  ! array size
  integer               :: i                         ! looping index
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'gas_exchange_out.nc'          ! output file
  real(r8),         parameter :: default_RH = 80.0_r8                      ! default relative humidity [%]
  real(r8),         parameter :: can_air_press = 101325.0_r8               ! default value for canopy air pressure [Pa]
  real(r8),         parameter :: can_o2_partial_press = 20900.0_r8         ! default value for O2 partial pressure [Pa]
  real(r8),         parameter :: default_veg_tempk = 25.0_r8 + tfrz        ! default value for vegetation temperature [K]
  real(r8),         parameter :: default_can_air_tempk = 25.0_r8 + tfrz    ! default value for canopy air temperature [K]
  real(r8),         parameter :: leaf_bl_resistance = 42.3_r8              ! default value for leaf boundary layer resistance [umol H2O/m2/s]
  real(r8),         parameter :: min_temp = -5.0_r8                        ! minimum temperature to calculate [degC]
  real(r8),         parameter :: max_temp = 50.0_r8                        ! maximum temperature to calculate [degC]
  real(r8),         parameter :: temp_inc = 0.5_r8                         ! temperature increment to use [degC]
  
  interface

    subroutine WriteGasExchangeData(out_file, num_temp, veg_tempk, mm_kco2, mm_ko2,      &
      co2_compensation_pt, cf, leaf_bl_conductance, constrained_air_vpress)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar, RegisterVar
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file                 
      integer,          intent(in) :: num_temp                  
      real(r8),         intent(in) :: veg_tempk(:)             
      real(r8),         intent(in) :: mm_kco2(:)                
      real(r8),         intent(in) :: mm_ko2(:)                 
      real(r8),         intent(in) :: co2_compensation_pt(:)    
      real(r8),         intent(in) :: cf(:)                     
      real(r8),         intent(in) :: leaf_bl_conductance(:)   
      real(r8),         intent(in) :: constrained_air_vpress(:) 
      
    end subroutine WriteGasExchangeData

  end interface
  
  ! calculate size of temperature array based on desired min/max/increment
  num_temp = int((max_temp - min_temp)/temp_inc + 1)

  ! allocate arrays
  allocate(veg_tempk(num_temp))
  allocate(mm_kco2(num_temp))
  allocate(mm_ko2(num_temp))
  allocate(co2_compensation_pt(num_temp))
  allocate(cf(num_temp))
  allocate(leaf_bl_conductance(num_temp))
  allocate(constrained_air_vpress(num_temp))
  
  ! initialize temperature array
  do i = 1, num_temp
    veg_tempk(i) = (min_temp + temp_inc*(i-1)) + tfrz
  end do
  
  ! get canopy gas parameters as we scale temperature and hold everything else constant
  do i = 1, num_temp
    
    ! calculate saturation vapor pressure and vapor pressure
    call CalcVaporPressure(veg_tempk(i), default_RH, veg_esat, veg_air_vpress)
    can_air_vpress = (default_RH/100.0_r8)*veg_esat
    
    ! get canopy gas parameters
    call GetCanopyGasParameters(can_air_press, can_o2_partial_press, veg_tempk(i),       &
      default_can_air_tempk, can_air_vpress, veg_esat, leaf_bl_resistance, mm_kco2(i),   &
      mm_ko2(i), co2_compensation_pt(i), cf(i), leaf_bl_conductance(i),                  &
      constrained_air_vpress(i))

  end do
  
  call WriteGasExchangeData(out_file, num_temp, veg_tempk, mm_kco2, mm_ko2,              &
    co2_compensation_pt, cf, leaf_bl_conductance, constrained_air_vpress)
      
end program FatesTestGasExchange

! ----------------------------------------------------------------------------------------

subroutine WriteGasExchangeData(out_file, num_temp, veg_tempk, mm_kco2, mm_ko2,          &
  co2_compensation_pt, cf, leaf_bl_conductance, constrained_air_vpress)
  !
  ! DESCRIPTION:
  ! Writes out data from the gas exchange test
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar
  use FatesUnitTestIOMod, only : RegisterVar
  use FatesUnitTestIOMod, only : EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int

  implicit none

  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file                  ! output file name
  integer,          intent(in) :: num_temp                  ! size of temperature array
  real(r8),         intent(in) :: veg_tempk(:)              ! vegetation temperature [K]
  real(r8),         intent(in) :: mm_kco2(:)                ! Michaelis-Menten constant for CO2 [Pa]
  real(r8),         intent(in) :: mm_ko2(:)                 ! Michaelis-Menten constant for O2 [Pa]
  real(r8),         intent(in) :: co2_compensation_pt(:)    ! CO2 compensation point [Pa]
  real(r8),         intent(in) :: cf(:)                     ! conversion factor between molar form and velocity form of conductance and resistance: [umol/m3]
  real(r8),         intent(in) :: leaf_bl_conductance(:)    ! leaf boundary layer conductance [umol H2O/m2/s]
  real(r8),         intent(in) :: constrained_air_vpress(:) ! vapor pressure of air, constrained [Pa]

  ! LOCALS:
  integer          :: i              ! looping index
  integer          :: ncid           ! netcdf file id
  character(len=8) :: dim_names(1)   ! dimension names
  integer          :: dimIDs(1)      ! dimension IDs
  integer          :: tempID         ! variable IDs for dimensions
  integer          :: kco2ID, ko2ID
  integer          :: comp_ptID, cfID
  integer          :: blID, air_vpressID   

  ! dimension names
  dim_names = [character(len=12) :: 'temperature']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/num_temp/), 1, dimIDs)

  ! register temperature
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'K', 'vegetation temperature'], 2, tempID)

  ! register M-M CO2
  call RegisterVar(ncid, 'mm_co2', dimIDs(1:1), type_double,                             &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                         &
    [character(len=150) :: 'temperature', 'Pa', 'Michaelis-Menten constant for CO2'],    &
    3, kco2ID)

  ! register M-M O2
  call RegisterVar(ncid, 'mm_o2', dimIDs(1:1), type_double,                          &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                     &
    [character(len=150) :: 'temperature', 'Pa', 'Michaelis-Menten constant for O2'], &
    3, ko2ID)
    
  ! register CO2 compensation point
  call RegisterVar(ncid, 'co2_comp_pt', dimIDs(1:1), type_double,          &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],           &
    [character(len=150) :: 'temperature', 'Pa', 'CO2 compensation point'], &
    3, comp_ptID)
    
  ! register cf
    call RegisterVar(ncid, 'cf', dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],   &
    [character(len=150) :: 'temperature', 'umol m-3', 'conversion factor between molar form and velocity form of conductance and resistance'], &
    3, cfID)
    
  ! register boundary layer conductance
    call RegisterVar(ncid, 'bl_conductance', dimIDs(1:1), type_double,                            &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                                  &
    [character(len=150) :: 'temperature', 'umol H2O m-2 s-1', 'leaf boundary layer conductance'], &
    3, blID)
    
 ! register constrained vapor pressure of air
    call RegisterVar(ncid, 'air_vpress_constrained', dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                       &
    [character(len=150) :: 'temperature', 'Pa', 'vapor pressure of air, constrained'], &
    3, air_vpressID)

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, tempID, veg_tempk(:))
  call WriteVar(ncid, kco2ID, mm_kco2(:))
  call WriteVar(ncid, ko2ID, mm_ko2(:))
  call WriteVar(ncid, comp_ptID, co2_compensation_pt(:))
  call WriteVar(ncid, cfID, cf(:))
  call WriteVar(ncid, blID, leaf_bl_conductance(:))
  call WriteVar(ncid, air_vpressID, constrained_air_vpress(:))
  
  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteGasExchangeData