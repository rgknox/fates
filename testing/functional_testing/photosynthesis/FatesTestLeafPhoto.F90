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
  real(r8),         allocatable      :: can_air_vpress(:)         ! canopy air vapor pressure [Pa]
  real(r8),         allocatable      :: veg_esat(:)               ! saturation vapor pressure at vegetation surface [Pa]
  real(r8),         allocatable      :: mm_kco2(:)                ! Michaelis-Menten constant for CO2 [Pa]
  real(r8),         allocatable      :: mm_ko2(:)                 ! Michaelis-Menten constant for O2 [Pa]
  real(r8),         allocatable      :: co2_compensation_pt(:)    ! CO2 compensation point [Pa]
  real(r8),         allocatable      :: cf(:)                     ! conversion factor between molar form and velocity form of conductance and resistance: [umol/m3]
  real(r8),         allocatable      :: leaf_bl_conductance(:)    ! leaf boundary layer conductance [umol H2O/m2/s]
  real(r8),         allocatable      :: constrained_air_vpress(:) ! vapor pressure of air, constrained [Pa]
  real(r8),         allocatable      :: vcmax(:,:)                ! maximum rate of carboxylation [umol CO2/m2/s]
  real(r8),         allocatable      :: jmax(:,:)                 ! maximum electron transport rate [umol electrons/m2/s]
  real(r8),         allocatable      :: co2_rcurve_islope(:,:)    ! initial slope of CO2 response curve (C4 plants)
  real(r8),         allocatable      :: assimilation(:,:)         ! carbon assimilation [umol C/m2/s]
  real(r8),         allocatable      :: stomatal_resistance(:,:)  ! stomatal resistance [s/m]
  real(r8),         allocatable      :: leaf_photo(:,:)           ! net leaf photosynthesis [umol CO2/m2/s]
  real(r8),         allocatable      :: c13disc(:,:)              ! carbon 13 in newly assimilated carbon
  real(r8),         allocatable      :: test_out(:,:)              ! test output
  real(r8)                           :: stomatal_intercept_btran  ! water-stressed minimum stomatal conductance [umol H2O/m2/s]
  real(r8)                           :: leaf_maintenance_resp     ! leaf maintenance respiration [umol CO2/m2/s]
  integer                            :: num_temp                  ! array size
  integer                            :: numpft                    ! number of PFTs (from parameter file)
  integer                            :: i                         ! looping index
  integer                            :: ft                        ! plant functional type index
  integer                            :: nargs                     ! number of command-line arguments
  integer                            :: arglen                    ! command-line argument length
  real(r8)                           :: vpd, can_air_esat, default_can_air_vpress
  real(r8)                           :: default_vpd
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'leaf_photo_out.nc' ! output file
  real(r8),         parameter :: default_RH = 15.0_r8                      ! default value for relative humidity [%]
  real(r8),         parameter :: can_air_press = 101325.0_r8               ! default value for canopy air pressure [Pa]
  real(r8),         parameter :: can_o2_partial_press = 20900.0_r8         ! default value for O2 partial pressure [Pa]
  real(r8),         parameter :: can_co2_partial_press = 39.5_r8           ! default value for CO2 partial pressure [Pa]
  real(r8),         parameter :: default_veg_tempk = 25.0_r8 + 273.15_r8   ! default value for vegetation temperature [K]
  real(r8),         parameter :: default_can_air_tempk = 25.0_r8 + 273.15_r8  ! default value for canopy air temperature [K]
  real(r8),         parameter :: leaf_bl_resistance = 42.3_r8             ! default value for leaf boundary layer resistance [s/m]
  real(r8),         parameter :: min_temp = 268.15                         ! minimum temperature to calculate [K]
  real(r8),         parameter :: max_temp = 323.15_r8                      ! maximum temperature to calculate [K]
  real(r8),         parameter :: temp_inc = 0.5_r8                         ! temperature increment to use [K]
  real(r8),         parameter :: default_nscaler = 1.0_r8                  ! default scaler for leaf nitrogen [0-1]
  real(r8),         parameter :: default_dayl_fact = 1.0_r8                ! default scaler for day length [0-1]
  real(r8),         parameter :: default_btran = 1.0_r8                    ! default scaler for BTRAN [0-1]
  real(r8),         parameter :: default_leaf_psi = 9999                   ! leaf water potential [MPa] - only effective if using plant hydro
  real(r8),         parameter :: default_par = 2000.0_r8                   ! default value for PAR [umol/m2/s]
  
  interface

  subroutine WriteLeafPhotosynthesis(out_file, num_temp, numpft, veg_tempk, veg_esat,    &
    can_air_vpress, assimilation, stomatal_resistance, leaf_photo, c13disc, test_out) 

    use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
    use FatesUnitTestIOMod, only : WriteVar, RegisterVar
    use FatesUnitTestIOMod, only : type_double, type_int
    use FatesConstantsMod,  only : r8 => fates_r8
    implicit none

    character(len=*), intent(in) :: out_file                
    integer,          intent(in) :: num_temp               
    integer,          intent(in) :: numpft                  
    real(r8),         intent(in) :: veg_tempk(:)
    real(r8),         intent(in) :: veg_esat(:)
    real(r8),         intent(in) :: can_air_vpress(:)            
    real(r8),         intent(in) :: assimilation(:,:)       
    real(r8),         intent(in) :: stomatal_resistance(:,:) 
    real(r8),         intent(in) :: leaf_photo(:,:)
    real(r8),         intent(in) :: c13disc(:,:)
    real(r8),         intent(in) :: test_out(:,:)  
  end subroutine WriteLeafPhotosynthesis

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
  
  ! calculate size of temperature array based on desired min/max/increment
  num_temp = int((max_temp - min_temp)/temp_inc + 1)
  numpft = size(prt_params%wood_density, dim=1)

  ! allocate arrays
  allocate(veg_tempk(num_temp))
  allocate(can_air_vpress(num_temp))
  allocate(veg_esat(num_temp))
  allocate(mm_kco2(num_temp))
  allocate(mm_ko2(num_temp))
  allocate(co2_compensation_pt(num_temp))
  allocate(cf(num_temp))
  allocate(leaf_bl_conductance(num_temp))
  allocate(constrained_air_vpress(num_temp))
  allocate(vcmax(num_temp, numpft))
  allocate(jmax(num_temp, numpft))
  allocate(co2_rcurve_islope(num_temp, numpft))
  allocate(assimilation(num_temp, numpft))
  allocate(stomatal_resistance(num_temp, numpft))
  allocate(leaf_photo(num_temp, numpft))
  allocate(c13disc(num_temp, numpft))
  allocate(test_out(num_temp, numpft))
  
  ! calculate saturation vapor pressure at canopy air (25)
  !call CalcVaporPressure(default_can_air_tempk, default_RH, default_vpd, can_air_esat,   &
  !  default_can_air_vpress)
  
  ! calculate vapor pressure with RH (80%)
  ! calculate vpd
  ! vpd_at_leaf_temp
  ! adjust canopy air so that VPD is the same

  ! calculate photosynthesis as we scale temperature and hold everything else constant
  do i = 1, num_temp
          
    ! initialize temperature array
    veg_tempk(i) = min_temp + temp_inc*(i-1)
    
    ! calculate saturation vapor pressure and vapor pressure
    call CalcVaporPressure(veg_tempk(i), default_rh, vpd, veg_esat(i), can_air_vpress(i))
    can_air_vpress(i) = (default_RH/100.0_r8)*veg_esat(i)
       
    ! calculate canopy gas parameters
    call GetCanopyGasParameters(can_air_press, can_o2_partial_press, veg_tempk(i),       &
      default_can_air_tempk, can_air_vpress(i), veg_esat(i), leaf_bl_resistance, mm_kco2(i), mm_ko2(i), &
      co2_compensation_pt(i), cf(i), leaf_bl_conductance(i), constrained_air_vpress(i))
      
    do ft = 1, numpft 
      
      ! calculate leaf biophysical rates
      call LeafLayerBiophysicalRates(100.0_r8, ft, EDPftvarcon_inst%vcmax25top(ft,1),    &
        param_derived%jmax25top(ft,1), param_derived%kp25top(ft,1), default_nscaler,     &
        veg_tempk(i), default_dayl_fact, veg_tempk(i), veg_tempk(i), default_btran,      &
        vcmax(i,ft), jmax(i,ft), co2_rcurve_islope(i,ft))
      
      ! pulled out from model for now
      stomatal_intercept_btran = max(cf(i)/2.E8_r8,                                      &
        EDPftvarcon_inst%stomatal_intercept(ft)*default_btran)
      
      if (ft == 12) then
        leaf_maintenance_resp = 0.025_r8*EDPftvarcon_inst%vcmax25top(ft,1)
      else 
        leaf_maintenance_resp = 0.015_r8*EDPftvarcon_inst%vcmax25top(ft,1)
      end if 
      
      ! calculate leaf-level photosynthesis
      call LeafLayerPhotosynthesis(1.0_r8, default_par, 0.0_r8, 10.0_r8, 0.0_r8, ft,     &
        vcmax(i,ft), jmax(i,ft), co2_rcurve_islope(i,ft), veg_tempk(i), veg_esat(i),        &
        can_air_press, can_co2_partial_press, can_o2_partial_press, default_btran,       &
        stomatal_intercept_btran, cf(i), leaf_bl_conductance(i),                         &
        constrained_air_vpress(i), mm_kco2(i), mm_ko2(i), co2_compensation_pt(i),        &
        leaf_maintenance_resp, default_leaf_psi, leaf_bl_resistance, assimilation(i,ft), &
        stomatal_resistance(i,ft), leaf_photo(i,ft), c13disc(i,ft), test_out(i,ft))
        
    end do
  end do
  
  call WriteLeafPhotosynthesis(out_file, num_temp, numpft, veg_tempk, veg_esat,          &
    constrained_air_vpress, assimilation, stomatal_resistance, leaf_photo, c13disc, test_out) 

end program FatesTestLeafPhoto

! ----------------------------------------------------------------------------------------

subroutine WriteLeafPhotosynthesis(out_file, num_temp, numpft, veg_tempk, veg_esat,      &
  can_air_vpress, assimilation, stomatal_resistance, leaf_photo, c13disc, test_out) 
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
  real(r8),         intent(in) :: veg_esat(:)              ! canopy air saturation vapor pressure [Pa]
  real(r8),         intent(in) :: can_air_vpress(:)        ! canopy air vapor pressure [Pa]
  real(r8),         intent(in) :: assimilation(:,:)        ! carbon assimilation [umol C/m2/s]
  real(r8),         intent(in) :: stomatal_resistance(:,:) ! stomatal resistance [s/m]
  real(r8),         intent(in) :: leaf_photo(:,:)          ! net leaf photosynthesis [umol CO2/m2/s]
  real(r8),         intent(in) :: c13disc(:,:)             ! carbon 13 in newly assimilated carbon
  real(r8),         intent(in) :: test_out(:,:)             ! test_output

  ! LOCALS:
  integer, allocatable :: pft_indices(:) ! array of pft indices to write out
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=12)    :: dim_names(2)   ! dimension names
  integer              :: dimIDs(2)      ! dimension IDs
  integer              :: tempID, pftID  ! variable IDs for dimensions
  integer              :: assimID, stomaID
  integer              :: photoID, c13ID
  integer              :: pressID, esatID
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
    
  ! register saturation vapor pressure
  call RegisterVar(ncid, 'veg_esat', dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'Pa', 'saturation vapor pressure'], 2, esatID)
    
    ! register vapor pressure
  call RegisterVar(ncid, 'can_air_vpress', dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'Pa', 'vapor pressure'], 2, pressID)
    
  ! register pft
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_int,         &
    [character(len=20)  :: 'units', 'long_name'],                     &
    [character(len=150) :: '', 'plant functional type'], 2, pftID)
    
  ! register assimilation
  call RegisterVar(ncid, 'assim', dimIDs(1:2), type_double,                                &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                           &
    [character(len=150) :: 'pft temperature', 'umol C m-2 s-1', 'carbon assimilation'],    &
    3, assimID)

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
    
  ! register c 13 descrimination
  call RegisterVar(ncid, 'c13_disc', dimIDs(1:2), type_double,                              &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                            &
    [character(len=150) :: 'pft temperature', '', 'carbon 13 in newly assimilated carbon'], &
    3, c13ID)
    
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
  call WriteVar(ncid, esatID, veg_esat(:))
  call WriteVar(ncid, pressID, can_air_vpress(:))
  call WriteVar(ncid, assimID, assimilation(:,:))
  call WriteVar(ncid, stomaID, stomatal_resistance(:,:))
  call WriteVar(ncid, photoID, leaf_photo(:,:))
  call WriteVar(ncid, c13ID, c13disc(:,:))
  call WriteVar(ncid, testID, test_out(:,:))
  
  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteLeafPhotosynthesis

