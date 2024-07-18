program FatesBiophysicalRates
  
  use FatesConstantsMod,           only : r8 => fates_r8
  use FATESPlantRespPhotosynthMod, only : LeafLayerBiophysicalRates
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use FatesParameterDerivedMod,    only : param_derived
  use PRTParametersMod,            only : prt_params
  
  implicit none
  
  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader ! param reader instance
  character(len=:), allocatable      :: param_file   ! input parameter file
  real(r8),         allocatable      :: veg_tempk(:) ! vegetation temperature [K]
  real(r8),         allocatable      :: vcmax(:,:)   ! maximum rate of carboxylation [umol CO2/m2/s]
  real(r8),         allocatable      :: jmax(:,:)    ! maximum electron transport rate [umol electrons/m2/s]
  real(r8),         allocatable      :: kp(:,:)      ! initial slope of CO2 response curve (C4 plants)
  integer                            :: num_temp     ! array size
  integer                            :: i            ! looping index
  integer                            :: ft           ! plant functional type index
  integer                            :: nargs        ! number of command-line arguments
  integer                            :: arglen       ! command-line argument length
  integer                            :: numpft       ! number of PFTs (from parameter file)
  
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'biophys_rates_out.nc' ! output file
  real(r8),         parameter :: par = 100.0_r8                    ! default PAR, just needs to be above 0.0 [W/m2]
  real(r8),         parameter :: default_nscaler = 1.0_r8          ! default scaler for leaf nitrogen [0-1]
  real(r8),         parameter :: default_dayl_fact = 1.0_r8        ! default scaler for day length [0-1]
  real(r8),         parameter :: default_btran = 1.0_r8            ! default scaler for BTRAN [0-1]
  real(r8),         parameter :: min_temp = 273.15_r8              ! minimum temperature to calculate [K]
  real(r8),         parameter :: max_temp = 323.15_r8              ! maximum temperature to calculate [K]
  real(r8),         parameter :: temp_inc = 0.5_r8                 ! temperature increment to use [K]
  
  interface

    subroutine WriteBiophysicalRates(out_file, num_temp, numpft, veg_tempk, vcmax, jmax, &
      kp)

      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar, RegisterVar
      use FatesUnitTestIOMod, only : type_double, type_int
      use FatesConstantsMod,  only : r8 => fates_r8
      implicit none

      character(len=*), intent(in) :: out_file                 
      integer,          intent(in) :: num_temp
      integer,          intent(in) :: numpft                 
      real(r8),         intent(in) :: veg_tempk(:)             
      real(r8),         intent(in) :: vcmax(:,:)                
      real(r8),         intent(in) :: jmax(:,:)                 
      real(r8),         intent(in) :: kp(:,:)
    end subroutine WriteBiophysicalRates

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
  
  ! pft array
  numpft = size(prt_params%wood_density, dim=1)

  ! allocate arrays
  allocate(veg_tempk(num_temp))
  allocate(vcmax(num_temp, numpft))
  allocate(jmax(num_temp, numpft))
  allocate(kp(num_temp, numpft))
  
  ! initialize temperature array
  do i = 1, num_temp
    veg_tempk(i) = min_temp + temp_inc*(i-1)
  end do
  
  ! get canopy gas parameters as we scale temperature and hold everything else constant
  do i = 1, num_temp
    do ft = 1, numpft
      call LeafLayerBiophysicalRates(par, ft, EDPftvarcon_inst%vcmax25top(ft,1),         &
        param_derived%jmax25top(ft,1), param_derived%kp25top(ft,1), default_nscaler,     &
        veg_tempk(i), default_dayl_fact, veg_tempk(i), veg_tempk(i), default_btran,      &
        vcmax(i,ft), jmax(i,ft), kp(i,ft))
    end do
  end do
  
  call WriteBiophysicalRates(out_file, num_temp, numpft, veg_tempk, vcmax, jmax, kp)
      
end program FatesBiophysicalRates

! ----------------------------------------------------------------------------------------

subroutine WriteBiophysicalRates(out_file, num_temp, numpft, veg_tempk, vcmax, jmax, kp) 
  !
  ! DESCRIPTION:
  ! Writes out data from the leaf biophysical rates test
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar
  use FatesUnitTestIOMod, only : RegisterVar
  use FatesUnitTestIOMod, only : EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int

  implicit none

  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file     ! output file name
  integer,          intent(in) :: num_temp     ! size of temperature array
  integer,          intent(in) :: numpft       ! number of pfts
  real(r8),         intent(in) :: veg_tempk(:) ! vegetation temperature [K]
  real(r8),         intent(in) :: vcmax(:,:)   ! maximum rate of carboxylation [umol co2/m2/s]
  real(r8),         intent(in) :: jmax(:,:)    ! maximum electron transport rate [umol electrons/m2/s]
  real(r8),         intent(in) :: kp(:,:)      ! initial slope of CO2 response curve (C4 plants)

  ! LOCALS:
  integer, allocatable :: pft_indices(:) ! array of pft indices to write out
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=12)    :: dim_names(2)   ! dimension names
  integer              :: dimIDs(2)      ! dimension IDs
  integer              :: tempID, pftID  ! variable IDs for dimensions
  integer              :: vcmaxID, jmaxID
  integer              :: kpID
  
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
    
  ! register pft
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_int,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: '', 'plant functional type'], 2, pftID)


  ! register vcmax
  call RegisterVar(ncid, 'vcmax', dimIDs(1:2), type_double,                                        &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                                   &
    [character(len=150) :: 'pft temperature', 'umol CO2 m-2 s-1', 'maximum rate of carboxylation'],    &
    3, vcmaxID)

  ! register jmax
  call RegisterVar(ncid, 'jmax', dimIDs(1:2), type_double,                        &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                  &
    [character(len=150) :: 'pft temperature', 'umol electrons m-2 s-1]', 'maximum electron transport rate'], &
    3, jmaxID)
    
  ! register kp
  call RegisterVar(ncid, 'kp', dimIDs(1:2), type_double,                              &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                      &
    [character(len=150) :: 'pft temperature', '', 'initial slope of CO2 response curve'], &
    3, kpID)
    
  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, tempID, veg_tempk(:))
  call WriteVar(ncid, pftID, pft_indices(:))
  call WriteVar(ncid, vcmaxID, vcmax(:,:))
  call WriteVar(ncid, jmaxID, jmax(:,:))
  call WriteVar(ncid, kpID, kp(:,:))
  
  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteBiophysicalRates