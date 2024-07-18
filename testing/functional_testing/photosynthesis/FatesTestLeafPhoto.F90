program FatesTestLeafPhoto
  
  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesTestPhotosynthesisMod,  only : CalcVaporPressure
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use PRTParametersMod,            only : prt_params
  use FatesConstantsMod,           only : tfrz => t_water_freeze_k_1atm
  
  implicit none
  
  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader           ! param reader instance
  character(len=:), allocatable      :: param_file             ! input parameter file
  real(r8),         allocatable      :: veg_tempk(:)           ! vegetation temperature [K]
  real(r8),         allocatable      :: PAR(:)                 ! PAR [umol/m2/s]
  real(r8),         allocatable      :: can_air_vpress(:)      ! canopy air vapor pressure [Pa]
  real(r8),         allocatable      :: leaf_photo_bytemp(:,:) ! net leaf photosynthesis with temperature [umol CO2/m2/s]
  real(r8),         allocatable      :: leaf_photo_byPAR(:,:)  ! net leaf photosynthesis with PAR [umol CO2/m2/s]
  real(r8),         allocatable      :: leaf_photo_byVPD(:,:)  ! net leaf photosynthesis with VPD [umol CO2/m2/s]
  real(r8),         allocatable      :: leaf_photo_byCO2(:,:)  ! net leaf photosynthesis with CO2 concentration [umol CO2/m2/s]
  real(r8)                           :: stomatal_resistance    ! stomatal resistance [s/m]
  real(r8)                           :: test_out               ! test output
  real(r8)                           :: veg_esat               ! saturation vapor pressure at vegetation surface [Pa]
  real(r8)                           :: veg_air_vpress         ! vapor pressure at vegetation surface [Pa]
  real(r8)                           :: default_can_air_vpress ! vapor pressure of canopy air [Pa]
  real(r8)                           :: default_veg_esat       ! default saturation vapor pressure at vegetation surface [Pa]
  integer                            :: num_temp, num_PAR      ! array sizes
  integer                            :: num_pft                ! number of PFTs (from parameter file)
  integer                            :: i                      ! looping index
  integer                            :: ft                     ! plant functional type index
  integer                            :: nargs                  ! number of command-line arguments
  integer                            :: arglen                 ! command-line argument length
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'leaf_photo_out.nc' ! output file
  
  ! DEFAULT VALUES ======================================================================
  real(r8), parameter :: default_can_air_press = 101325.0_r8  ! default value for canopy air pressure [Pa]
  real(r8), parameter :: default_can_co2_pp = 39.5_r8         ! default value for CO2 partial pressure [Pa]
  real(r8), parameter :: default_can_o2_pp = 20900.0_r8       ! default value for O2 partial pressure [Pa]
  real(r8), parameter :: default_can_tempk = 25.0_r8 + tfrz   ! default value for canopy air temperature [K]
  real(r8), parameter :: default_veg_tempk = 25.0_r8 + tfrz   ! default value for vegetation temperature [K]
  real(r8), parameter :: default_RH = 80.0_r8                 ! default value for relative humidity [%]  
  real(r8), parameter :: default_par = 2000.0_r8              ! default value for PAR [umol/m2/s]
  real(r8), parameter :: default_leaf_bl_resistance = 42.3_r8 ! default value for leaf boundary layer resistance [s/m]
  real(r8), parameter :: default_nscaler = 1.0_r8             ! default scaler for leaf nitrogen [0-1]
  real(r8), parameter :: default_dayl_fact = 1.0_r8           ! default scaler for day length [0-1]
  real(r8), parameter :: default_btran = 1.0_r8               ! default scaler for BTRAN [0-1]

  ! FOR CHANGING VALUES ==================================================================
  real(r8), parameter :: min_temp = -5.0_r8  ! minimum temperature to calculate [degC]
  real(r8), parameter :: max_temp = 50.0_r8  ! maximum temperature to calculate [degC]
  real(r8), parameter :: temp_inc = 0.5_r8   ! temperature increment to use [degC]
  real(r8), parameter :: min_PAR = 0.0_r8    ! minimum PAR to calculate [umol/m2/s]
  real(r8), parameter :: max_PAR = 1500.0_r8 ! maximum PAR to calculate [umol/m2/s]
  real(r8), parameter :: PAR_inc = 5.0_r8    ! PAR increment to use [umol/m2/s]
  real(r8), parameter :: min_VPD = 0.0_r8    ! minimum VPD to calculate [Pa]
  real(r8), parameter :: max_VPD = 4000.0_r8 ! maximum VPD to calculate [Pa]
  real(r8), parameter :: VPD_inc = 5.0_r8    ! VPD increment to use [Pa]
  real(r8), parameter :: min_co2 = 50.0_r8   ! minimum CO2 concentration to calculate [umol/mol]
  real(r8), parameter :: max_co2 = 600.0_r8  ! maximum CO2 concentration to calculate [umol/mol]
  real(r8), parameter :: co2_inc = 5.0_r8    ! CO2 concentration increment to use [umol/mol]
  
  interface

  subroutine WriteLeafPhotosynthesis(out_file, num_temp, numpft, veg_tempk, can_air_vpress, &
    leaf_photo_bytemp, leaf_photo_byPAR)  

    use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
    use FatesUnitTestIOMod, only : WriteVar, RegisterVar
    use FatesUnitTestIOMod, only : type_double, type_int
    use FatesConstantsMod,  only : r8 => fates_r8
    implicit none

    character(len=*), intent(in) :: out_file                
    integer,          intent(in) :: num_temp               
    integer,          intent(in) :: numpft                  
    real(r8),         intent(in) :: veg_tempk(:)
    real(r8),         intent(in) :: can_air_vpress(:)
    real(r8),         intent(in) :: leaf_photo_bytemp(:,:)
    real(r8),         intent(in) :: leaf_photo_byPAR(:,:)          
  end subroutine WriteLeafPhotosynthesis
  
  subroutine LeafLevelPhoto(can_air_press, can_co2_pp, can_o2_pp, veg_tempk, can_tempk,  &
    can_vpress, veg_esat, par, leaf_bl_resistance, nscaler, dayl_fact, btran, ft,        &
    leaf_photo, stomatal_resistance, test_out)
    
    use FatesConstantsMod,           only : r8 => fates_r8
    use FatesTestPhotosynthesisMod,  only : CalcVaporPressure
    use FATESPlantRespPhotosynthMod, only : GetCanopyGasParameters
    use EDPftvarcon,                 only : EDPftvarcon_inst
    use FatesParameterDerivedMod,    only : param_derived
    implicit none
    
    real(r8), intent(in)  :: can_air_press      
    real(r8), intent(in)  :: can_co2_pp          
    real(r8), intent(in)  :: can_o2_pp           
    real(r8), intent(in)  :: veg_tempk          
    real(r8), intent(in)  :: can_tempk          
    real(r8), intent(in)  :: can_vpress         
    real(r8), intent(in)  :: veg_esat           
    real(r8), intent(in)  :: par                 
    real(r8), intent(in)  :: leaf_bl_resistance
    real(r8), intent(in)  :: nscaler            
    real(r8), intent(in)  :: dayl_fact          
    real(r8), intent(in)  :: btran               
    integer,  intent(in)  :: ft                  
    real(r8), intent(out) :: leaf_photo         
    real(r8), intent(out) :: stomatal_resistance 
    real(r8), intent(out) :: test_out           
  end subroutine LeafLevelPhoto

end interface
  
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
  
  ! calculate sizes of arrays based on desired min/max/increments
  num_temp = int((max_temp - min_temp)/temp_inc + 1)
  num_PAR = int((max_PAR - min_PAR)/PAR_inc + 1)
  num_pft = size(prt_params%wood_density, dim=1)
  
  ! allocate arrays
  allocate(veg_tempk(num_temp))
  allocate(can_air_vpress(num_temp))
  allocate(leaf_photo_bytemp(num_temp, numpft))
  allocate(PAR(num_PAR))
  allocate(leaf_photo_byPAR(num_PAR, numpft))
  
  ! calculate saturation vapor pressure and vapor pressure for default values
  call CalcVaporPressure(default_veg_tempk, default_RH, default_veg_esat, veg_air_vpress)
  default_can_air_vpress = (default_RH/100.0_r8)*default_veg_esat
  
  !!! PAR --------------------------------------------------------------------------------

  ! calculate photosynthesis as we scale PAR and hold everything else constant
  do i = 1, num_PAR
    
    ! calculate PAR
    PAR(i) = min_PAR + PAR_inc*(i-1)
    
    do ft = 1, numpft 
      call LeafLevelPhoto(default_can_air_press, default_can_co2_pp, default_can_o2_pp,  &
      default_veg_tempk, default_can_tempk, default_can_air_vpress, default_veg_esat,    &
      par(i), default_leaf_bl_resistance, default_nscaler, default_dayl_fact,            &
      default_btran, ft, leaf_photo_byPAR(i,ft), stomatal_resistance, test_out)
    end do
  end do
  
  !!! Leaf temperature -------------------------------------------------------------------
  
  ! calculate photosynthesis as we scale temperature and hold everything else constant
  do i = 1, num_temp
          
    ! calculate temperature
    veg_tempk(i) = (min_temp + temp_inc*(i-1)) + tfrz
    
    ! calculate saturation vapor pressure and vapor pressure
    call CalcVaporPressure(veg_tempk(i), default_RH, veg_esat, veg_air_vpress)
    
    ! calculate canopy air vapor pressure so that RH and VPD are maintained
    can_air_vpress(i) = (default_RH/100.0_r8)*veg_esat
      
    do ft = 1, numpft 
      call LeafLevelPhoto(default_can_air_press, default_can_co2_pp, default_can_o2_pp,  &
      veg_tempk(i), default_can_tempk, can_air_vpress(i), veg_esat, default_par,         &
      default_leaf_bl_resistance, default_nscaler, default_dayl_fact, default_btran, ft, &
      leaf_photo_bytemp(i,ft), stomatal_resistance, test_out)
    end do
  end do
  
  !!! Write output and wrap up -----------------------------------------------------------
    
  call WriteLeafPhotosynthesis(out_file, num_temp, numpft, veg_tempk, can_air_vpress,    &
    leaf_photo) 

  ! deallocate arrays
  if (allocated(veg_tempk)) deallocate(veg_tempk)
  if (allocated(PAR)) deallocate(PAR)
  if (allocated(can_air_vpress)) deallocate(can_air_vpress)
  if (allocated(leaf_photo_bytemp)) deallocate(leaf_photo_bytemp)
  if (allocated(leaf_photo_byPAR)) deallocate(leaf_photo_byPAR)
  
end program FatesTestLeafPhoto

! ----------------------------------------------------------------------------------------

subroutine LeafLevelPhoto(can_air_press, can_co2_pp, can_o2_pp, veg_tempk, can_tempk,    &
    can_vpress, veg_esat, par, leaf_bl_resistance, nscaler, dayl_fact, btran, ft,        &
    leaf_photo, stomatal_resistance, test_out)
  !
  ! DESCRIPTION:
  ! Drives and individual photosynthesis calculation
  !
  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesTestPhotosynthesisMod,  only : CalcVaporPressure
  use FATESPlantRespPhotosynthMod, only : GetCanopyGasParameters
  use FATESPlantRespPhotosynthMod, only : LeafLayerBiophysicalRates
  use FATESPlantRespPhotosynthMod, only : LeafLayerPhotosynthesis
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use FatesParameterDerivedMod,    only : param_derived
  
  implicit none
  
  ! ARGUMENTS:
  real(r8), intent(in)  :: can_air_press       ! canopy air pressure [Pa]
  real(r8), intent(in)  :: can_co2_pp          ! canopy air CO2 partial pressure [Pa]
  real(r8), intent(in)  :: can_o2_pp           ! canopy air O2 partial pressure [Pa]
  real(r8), intent(in)  :: veg_tempk           ! leaf temperature [K]
  real(r8), intent(in)  :: can_tempk           ! canopy air temperature [K]
  real(r8), intent(in)  :: can_vpress          ! canopy air vapor pressure [Pa]
  real(r8), intent(in)  :: veg_esat            ! saturation vapor pressure at leaf surface [Pa]
  real(r8), intent(in)  :: par                 ! PAR [umol/m2/s]
  real(r8), intent(in)  :: leaf_bl_resistance  ! leaf boundary layer resistance [s/m]
  real(r8), intent(in)  :: nscaler             ! scaler for leaf nitrogen [0-1]
  real(r8), intent(in)  :: dayl_fact           ! scaler for day length factor [0-1]
  real(r8), intent(in)  :: btran               ! scaler for BTRAN [0-1]
  integer,  intent(in)  :: ft                  ! plant functional type index
  real(r8), intent(out) :: leaf_photo          ! net assimilation [umol C/m2/s]
  real(r8), intent(out) :: stomatal_resistance ! stomatal resistance [s/m]
  real(r8), intent(out) :: test_out            ! testing output
  
  ! LOCALS:
  real(r8) :: mm_kco2                ! Michaelis-Menten constant for CO2 [Pa]
  real(r8) :: mm_ko2                 ! Michaelis-Menten constant for O2 [Pa]
  real(r8) :: co2_compensation_pt    ! CO2 compensation point [Pa]
  real(r8) :: cf                     !  conversion factor between molar form and velocity form of conductance and resistance: [umol/m3]
  real(r8) :: leaf_bl_conductance    ! leaf boundary layer conductance [umol H2O/m2/s]
  real(r8) :: can_vpress_constrained ! canopy air vapor pressure, constrained [Pa]
  real(r8) :: vcmax                  ! maximum rate of carboxylation [umol CO2/m2/s]
  real(r8) :: jmax                   ! maximum electron transport rate [umol electrons/m2/s]
  real(r8) :: co2_rcurve_islope      ! initial slope of CO2 response curve (C4 plants)
  real(r8) :: stomatal_int_btran     ! water-stressed minimum stomatal conductance [umol H2O/m2/s]
  real(r8) :: leaf_maintenance_resp  ! leaf maintenance respiration [umol CO2/m2/s]
  real(r8) :: assimilation           ! CO2 assimilation [umol C/m2/s]
  real(r8) :: c13disc                ! C13 discrimination
  
  ! calculate canopy gas parameters
  call GetCanopyGasParameters(can_air_press, can_o2_pp, veg_tempk, can_tempk,            &
    can_vpress, veg_esat, leaf_bl_resistance, mm_kco2, mm_ko2, co2_compensation_pt,      &
    cf, leaf_bl_conductance, can_vpress_constrained)
    
  ! calculate leaf biophysical rates
  call LeafLayerBiophysicalRates(par, ft, EDPftvarcon_inst%vcmax25top(ft,1),             &
    param_derived%jmax25top(ft,1), param_derived%kp25top(ft,1), nscaler, veg_tempk,      &
    dayl_fact, can_tempk, can_tempk, btran, vcmax, jmax, co2_rcurve_islope)
    
  ! pulled out from model for now
  stomatal_int_btran = max(cf/2.E8_r8, EDPftvarcon_inst%stomatal_intercept(ft)*btran)
  
  if (ft == 12) then
    leaf_maintenance_resp = 0.025_r8*EDPftvarcon_inst%vcmax25top(ft,1)
  else 
    leaf_maintenance_resp = 0.015_r8*EDPftvarcon_inst%vcmax25top(ft,1)
  end if 
  
  ! calculate leaf-level photosynthesis
  call LeafLayerPhotosynthesis(1.0_r8, par, 0.0_r8, 10.0_r8, 0.0_r8, ft, vcmax, jmax,    &
    co2_rcurve_islope, veg_tempk, veg_esat, can_air_press, can_co2_pp, can_o2_pp, btran, &
    stomatal_int_btran, cf, leaf_bl_conductance, can_vpress_constrained, mm_kco2,        &
    mm_ko2, co2_compensation_pt, leaf_maintenance_resp, 999.0_r8, leaf_bl_resistance,    &
    assimilation, stomatal_resistance, leaf_photo, c13disc, test_out)
  
end subroutine LeafLevelPhoto

! ----------------------------------------------------------------------------------------

subroutine WriteLeafPhotosynthesis(out_file, num_temp, numpft, veg_tempk, can_air_vpress, &
  leaf_photo_bytemp, leaf_photo_byPAR)  
  !
  ! DESCRIPTION:
  ! Writes out data from the leaf photosynthesis test
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar
  use FatesUnitTestIOMod, only : RegisterVar
  use FatesUnitTestIOMod, only : EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int

  implicit none

  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file                 ! output file name
  integer,          intent(in) :: num_temp                 ! size of temperature array
  integer,          intent(in) :: numpft                   ! number of pfts
  real(r8),         intent(in) :: veg_tempk(:)             ! vegetation temperature [K]
  real(r8),         intent(in) :: can_air_vpress(:)        ! canopy air vapor pressure [Pa]
  real(r8),         intent(in) :: leaf_photo_bytemp(:,:)   ! net leaf photosynthesis [umol CO2/m2/s]
  real(r8),         intent(in) :: leaf_photo_byPAR(:,:)    ! net leaf photosynthesis [umol CO2/m2/s]

  ! LOCALS:
  integer, allocatable :: pft_indices(:) ! array of pft indices to write out
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=12)    :: dim_names(2)   ! dimension names
  integer              :: dimIDs(2)      ! dimension IDs
  integer              :: tempID, pftID  ! variable IDs for dimensions
  integer              :: stomaID
  integer              :: photoID
  integer              :: pressID
  integer              :: testID
  
  ! create pft indices
  allocate(pft_indices(numpft))
  do i = 1, numpft
    pft_indices(i) = i
  end do

  ! dimension names
  dim_names = [character(len=12) :: 'temperature', 'pft']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/num_temp, numpft/), 2, dimIDs)

  ! register temperature
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'K', 'vegetation temperature'], 2, tempID)
    
    ! register vapor pressure
  call RegisterVar(ncid, 'can_air_vpress', dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'Pa', 'vapor pressure'], 2, pressID)
    
  ! register pft
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_int,         &
    [character(len=20)  :: 'units', 'long_name'],                     &
    [character(len=150) :: '', 'plant functional type'], 2, pftID)
    
  ! register stomatal resistance
  call RegisterVar(ncid, 'stomatal_resistance', dimIDs(1:2), type_double,       &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                &
    [character(len=150) :: 'pft temperature', 's m-1]', 'stomatal resistance'], &
    3, stomaID)
    
  ! register leaf photosynthesis
  call RegisterVar(ncid, 'leaf_photo', dimIDs(1:2), type_double,                              &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                              &
    [character(len=150) :: 'pft temperature', 'umol CO2 m-2 s-1', 'net leaf photosynthesis'], &
    3, photoID)
    
  ! register test output
    call RegisterVar(ncid, 'test_output', dimIDs(1:2), type_double,                              &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                            &
    [character(len=150) :: 'pft temperature', '', 'for testing'], &
    3, testID)
    
  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, tempID, veg_tempk(:))
  call WriteVar(ncid, pftID, pft_indices(:))
  call WriteVar(ncid, pressID, can_air_vpress(:))
  call WriteVar(ncid, stomaID, stomatal_resistance(:,:))
  call WriteVar(ncid, photoID, leaf_photo(:,:))
  call WriteVar(ncid, testID, test_out(:,:))
  
  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteLeafPhotosynthesis

