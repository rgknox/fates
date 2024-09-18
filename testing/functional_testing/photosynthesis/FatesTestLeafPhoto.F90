program FatesTestLeafPhoto
  
  use FatesConstantsMod,          only : r8 => fates_r8
  use FatesTestPhotosynthesisMod,  only : CalcVaporPressure
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use PRTParametersMod,            only : prt_params
  use FatesConstantsMod,           only : tfrz => t_water_freeze_k_1atm
  use FatesConstantsMod,           only : wm2_to_umolm2s
  
  implicit none
  
                                                                  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader              ! param reader instance
  character(len=:), allocatable      :: param_file                ! input parameter file
  real(r8),         allocatable      :: veg_tempk(:)              ! vegetation temperature [K]
  real(r8),         allocatable      :: PAR(:)                    ! PAR [umol/m2/s]
  real(r8),         allocatable      :: RH(:)                     ! relative humidity [%]
  real(r8),         allocatable      :: co2_pp(:)                 ! CO2 partial pressure [Pa]
  real(r8),         allocatable      :: nscaler(:)                ! scaler for leaf N [0-1]
  real(r8),         allocatable      :: btran(:)                  ! scaler for BTRAN [0-1]
  real(r8),         allocatable      :: can_air_vpress_bytemp(:)  ! canopy air vapor pressure (by temperature) [Pa]
  real(r8),         allocatable      :: can_air_vpress_byRH(:)    ! canopy air vapor pressure (by RH) [Pa]
  real(r8),         allocatable      :: veg_esat_bytemp(:)        ! saturation vapor pressure at vegetation surface - by temperature [Pa]
  real(r8),         allocatable      :: veg_esat_byRH(:)          ! saturation vapor pressure at vegetation surface - by RH [Pa]
  real(r8),         allocatable      :: leaf_photo_bytemp(:,:)    ! net leaf photosynthesis with temperature [umol CO2/m2/s]
  real(r8),         allocatable      :: leaf_photo_byPAR(:,:)     ! net leaf photosynthesis with PAR [umol CO2/m2/s]
  real(r8),         allocatable      :: leaf_photo_byRH(:,:)      ! net leaf photosynthesis with RH [umol CO2/m2/s]
  real(r8),         allocatable      :: leaf_photo_byCO2(:,:)     ! net leaf photosynthesis with CO2 concentration [umol CO2/m2/s]
  real(r8),         allocatable      :: leaf_photo_bynscaler(:,:) ! net leaf photosynthesis with nscaler [umol CO2/m2/s]
  real(r8),         allocatable      :: leaf_photo_bybtran(:,:)   ! net leaf photosynthesis with BTRAN [umol CO2/m2/s]
  real(r8)                           :: stomatal_conductance      ! stomatal conducance [umol h2o/m2/s]
  real(r8)                           :: net_photo                 ! net photosynthesis [umol h2o/m2/s]
  real(r8)                           :: c13disc                   ! c13 fraction of photosynthesis
  real(r8)                           :: veg_air_vpress            ! vapor pressure at vegetation surface [Pa]
  real(r8)                           :: default_can_air_vpress    ! default value for vapor pressure of canopy air [Pa]
  real(r8)                           :: default_veg_esat          ! default saturation vapor pressure at vegetation surface [Pa]
  real(r8)                           :: veg_temp                  ! vegetation temperature [K]
  integer                            :: num_temp, num_PAR         ! array sizes
  integer                            :: num_RH, num_co2           ! array sizes
  integer                            :: num_nscaler, num_btran    ! array sizes
  integer                            :: num_pft                   ! number of PFTs (from parameter file)
  integer                            :: i                         ! looping index
  integer                            :: ft                        ! plant functional type index
  integer                            :: nargs                     ! number of command-line arguments
  integer                            :: arglen                    ! command-line argument length
  
                                                                  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'leaf_photo_out.nc'   ! output file
  
                                                                  ! DEFAULT VALUES ======================================================================
  real(r8), parameter :: default_can_air_press = 101325.0_r8      ! default value for canopy air pressure [Pa]
  real(r8), parameter :: default_can_co2_pp = 38.5035_r8          ! default value for CO2 partial pressure [Pa]
  real(r8), parameter :: default_can_o2_pp = 20900.0_r8           ! default value for O2 partial pressure [Pa]
  real(r8), parameter :: default_can_tempk = 25.0_r8 + tfrz       ! default value for canopy air temperature [K]
  real(r8), parameter :: default_veg_tempk = 25.0_r8 + tfrz       ! default value for vegetation temperature [K]
  real(r8), parameter :: default_RH = 80.0_r8                     ! default value for relative humidity [%]  
  real(r8), parameter :: default_par = 2000.0_r8                  ! default value for PAR [umol/m2/s] ~ 430 w/m2
  real(r8), parameter :: default_leaf_bl_res_vel = 42.3_r8        ! default value for leaf boundary layer resistance [s/m] (velocity units)
  real(r8), parameter :: rgas = 8314.4598                         ! universal gas constant [J/K/kmol]
  real(r8), parameter :: umol_per_kmol = 1.0E9_r8                 ! number of umols in a kmol (many)
  real(r8), parameter :: default_leaf_bl_cond = (1._r8/default_leaf_bl_res_vel)*default_can_air_press/(rgas * default_can_tempk )*umol_per_kmol
  real(r8), parameter :: default_nscaler = 1.0_r8                 ! default scaler for leaf nitrogen [0-1]
  real(r8), parameter :: default_dayl_fact = 1.0_r8               ! default scaler for day length [0-1]
  real(r8), parameter :: default_btran = 1.0_r8                   ! default scaler for BTRAN [0-1]

                                                                  ! FOR CHANGING VALUES ==================================================================
  real(r8), parameter :: min_temp = -5.0_r8                       ! minimum temperature to calculate [degC]
  real(r8), parameter :: max_temp = 50.0_r8                       ! maximum temperature to calculate [degC]
  real(r8), parameter :: temp_inc = 0.5_r8                        ! temperature increment to use [degC]
  real(r8), parameter :: min_PAR = 0.0_r8                         ! minimum PAR to calculate [watts/m2/s]
  real(r8), parameter :: max_PAR = 1360.0_r8                      ! maximum PAR to calculate [watts/m2/s]
  real(r8), parameter :: PAR_inc = 5.0_r8                         ! PAR increment to use [umol/m2/s]
  real(r8), parameter :: min_RH = 0.0_r8                          ! minimum RH to calculate [%]
  real(r8), parameter :: max_RH = 100.0_r8                        ! maximum RH to calculate [%]
  real(r8), parameter :: RH_inc = 1.0_r8                          ! RH increment to use [%]
  real(r8), parameter :: min_co2 = 50.0_r8                        ! minimum CO2 concentration to calculate [umol/mol]
  real(r8), parameter :: max_co2 = 600.0_r8                       ! maximum CO2 concentration to calculate [umol/mol]
  real(r8), parameter :: co2_inc = 5.0_r8                         ! CO2 concentration increment to use [umol/mol]
  
  interface

  subroutine WriteLeafPhotosynthesis(out_file, num_temp, num_pft, num_par, num_RH,       &
    num_co2, num_nscaler, num_btran, veg_tempk, par, RH, co2, nscaler, btran,             &
    can_air_vpress_bytemp, veg_esat_bytemp, veg_esat_byRH, can_air_vpress_byRH,          &
    leaf_photo_bytemp, leaf_photo_byPAR, leaf_photo_byRH, leaf_photo_byCO2,              &
    leaf_photo_bynscaler, leaf_photo_bybtran)  
    use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
    use FatesUnitTestIOMod, only : WriteVar, RegisterVar
    use FatesUnitTestIOMod, only : type_double, type_int
    use FatesConstantsMod,  only : r8 => fates_r8

    character(len=*), intent(in) :: out_file                
    integer,          intent(in) :: num_temp               
    integer,          intent(in) :: num_pft    
    integer,          intent(in) :: num_par
    integer,          intent(in) :: num_RH
    integer,          intent(in) :: num_co2
    integer,          intent(in) :: num_nscaler     
    integer,          intent(in) :: num_btran           
    real(r8),         intent(in) :: veg_tempk(:)
    real(r8),         intent(in) :: par(:)
    real(r8),         intent(in) :: RH(:)
    real(r8),         intent(in) :: co2(:)
    real(r8),         intent(in) :: nscaler(:)
    real(r8),         intent(in) :: btran(:)
    real(r8),         intent(in) :: can_air_vpress_bytemp(:)
    real(r8),         intent(in) :: can_air_vpress_byRH(:)
    real(r8),         intent(in) :: veg_esat_bytemp(:)
    real(r8),         intent(in) :: veg_esat_byRH(:)
    real(r8),         intent(in) :: leaf_photo_bytemp(:,:)
    real(r8),         intent(in) :: leaf_photo_byPAR(:,:)   
    real(r8),         intent(in) :: leaf_photo_byRH(:,:)  
    real(r8),         intent(in) :: leaf_photo_byco2(:,:)
    real(r8),         intent(in) :: leaf_photo_bynscaler(:,:) 
    real(r8),         intent(in) :: leaf_photo_bybtran(:,:)     
  end subroutine WriteLeafPhotosynthesis
  
  subroutine LeafLevelPhoto(can_air_press, can_co2_pp, can_o2_pp, veg_tempk, can_tempk,  &
       can_vpress, veg_esat, par, leaf_bl_conductance, nscaler, dayl_fact, btran, ft,        &
       gross_photo, stomatal_cond, net_photo, c13disc)

    use FatesConstantsMod,          only : r8 => fates_r8
    
    real(r8), intent(in)  :: can_air_press      
    real(r8), intent(in)  :: can_co2_pp          
    real(r8), intent(in)  :: can_o2_pp           
    real(r8), intent(in)  :: veg_tempk          
    real(r8), intent(in)  :: can_tempk          
    real(r8), intent(in)  :: can_vpress         
    real(r8), intent(in)  :: veg_esat           
    real(r8), intent(in)  :: par                 
    real(r8), intent(in)  :: leaf_bl_conductance
    real(r8), intent(in)  :: nscaler            
    real(r8), intent(in)  :: dayl_fact          
    real(r8), intent(in)  :: btran               
    integer,  intent(in)  :: ft                  
    real(r8), intent(out) :: gross_photo
    real(r8), intent(out) :: stomatal_cond
    real(r8), intent(out) :: net_photo
    real(r8), intent(out) :: c13disc
  end subroutine LeafLevelPhoto

end interface
  
                                                                  ! get parameter file from command-line argument
  nargs = command_argument_count()
  if (nargs /= 1) then
    write(*, '(a,i2,a)') "Incorrect number of arguments: ", nargs, ". Should be 1."
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
  num_RH = int((max_RH - min_RH)/RH_inc + 1)
  num_co2 = int((max_co2 - min_co2)/co2_inc + 1)
  num_nscaler = int((1.0 - 0.0)/0.05 + 1)
  num_btran = int((1.0 - 0.0)/0.05 + 1)
  
                                                                  ! number of pfts from parameter file
  num_pft = size(prt_params%wood_density, dim=1)
  
                                                                  ! allocate arrays
  allocate(veg_tempk(num_temp))
  allocate(PAR(num_PAR))
  allocate(RH(num_RH))
  allocate(co2_pp(num_co2))
  allocate(nscaler(num_nscaler))
  allocate(btran(num_btran))
  allocate(can_air_vpress_bytemp(num_temp))
  allocate(can_air_vpress_byRH(num_RH))
  allocate(veg_esat_bytemp(num_temp))
  allocate(veg_esat_byRH(num_RH))
  allocate(leaf_photo_bytemp(num_temp, num_pft))
  allocate(leaf_photo_byPAR(num_PAR, num_pft))
  allocate(leaf_photo_byRH(num_RH, num_pft))
  allocate(leaf_photo_byco2(num_co2, num_pft))
  allocate(leaf_photo_bynscaler(num_nscaler, num_pft))
  allocate(leaf_photo_bybtran(num_btran, num_pft))
  
                                                                  ! calculate saturation vapor pressure and vapor pressure for default values
  call CalcVaporPressure(default_veg_tempk, default_RH, default_veg_esat, veg_air_vpress)
  default_can_air_vpress = (default_RH/100.0_r8)*default_veg_esat
  
                                                                  !!! PAR --------------------------------------------------------------------------------

                                                                  ! calculate photosynthesis as we scale PAR and hold everything else constant
  do i = 1, num_PAR
    
                                                                  ! calculate PAR and convert from watts/m2 to umols/m2/s
    PAR(i) = (min_PAR + PAR_inc*(i-1))*wm2_to_umolm2s
    
    do ft = 1, num_pft 
      call LeafLevelPhoto(default_can_air_press, default_can_co2_pp, default_can_o2_pp,  &
        default_veg_tempk, default_can_tempk, default_can_air_vpress, default_veg_esat,  &
        PAR(i), default_leaf_bl_cond, default_nscaler, default_dayl_fact,          &
        default_btran, ft, &
        leaf_photo_byPAR(i,ft), &
        stomatal_conductance, net_photo,c13disc)
    end do
  end do
  
  
                                                                  !!! CO2 concentration ------------------------------------------------------------------

                                                                  ! calculate photosynthesis as we scale CO2 ppm and hold everything else constant
  do i = 1, num_co2
    
                                                                  ! calculate PAR
    co2_pp(i) = ((min_co2 + co2_inc*(i-1))/1.0E6_r8)*default_can_air_press
    
    do ft = 1, num_pft 
      call LeafLevelPhoto(default_can_air_press, co2_pp(i), default_can_o2_pp,           &
        default_veg_tempk, default_can_tempk, default_can_air_vpress, default_veg_esat,  &
        default_par, default_leaf_bl_cond, default_nscaler, default_dayl_fact,     &
        default_btran, ft, &
        leaf_photo_byCO2(i,ft), stomatal_conductance, net_photo, c13disc)
    end do
  end do
  
  
                                                                  !!! RH --------------------------------------------------------------------------------

                                                                  ! calculate photosynthesis as we scale RH and hold everything else constant
  do i = 1, num_RH
    
                                                                  ! calculate PAR
    RH(i) = min_RH + RH_inc*(i-1)
    
                                                                  ! if (RH(i) >= 30.0_r8) then 
                                                                  !   veg_temp = default_veg_tempk
                                                                  ! else 
                                                                  !   veg_temp = 35.0_r8 + tfrz
                                                                  ! end if 
    
    veg_temp = default_veg_tempk
    
                                                                  ! calculate saturation vapor pressure and vapor pressure
    call CalcVaporPressure(veg_temp, RH(i), veg_esat_byRH(i), veg_air_vpress)
    
                                                                  ! calculate canopy air vapor pressure at this RH
    can_air_vpress_byRH(i) = (RH(i)/100.0_r8)*veg_esat_byRH(i)
    
    do ft = 1, num_pft 
      call LeafLevelPhoto(default_can_air_press, default_can_co2_pp, default_can_o2_pp,  &
        veg_temp, default_can_tempk, can_air_vpress_byRH(i), veg_esat_byRH(i),           &
        default_PAR, default_leaf_bl_cond, default_nscaler, default_dayl_fact,     &
        default_btran, ft, &
        leaf_photo_byRH(i,ft), stomatal_conductance, net_photo, c13disc)
    end do
  end do
  
                                                                  !!! Leaf temperature -------------------------------------------------------------------
  
                                                                  ! calculate photosynthesis as we scale temperature and hold everything else constant
  do i = 1, num_temp
          
                                                                  ! calculate temperature
    veg_tempk(i) = (min_temp + temp_inc*(i-1)) + tfrz
    
                                                                  ! calculate saturation vapor pressure and vapor pressure
    call CalcVaporPressure(veg_tempk(i), default_RH, veg_esat_bytemp(i), veg_air_vpress)
    
                                                                  ! calculate canopy air vapor pressure so that RH and VPD are maintained
    can_air_vpress_bytemp(i) = (default_RH/100.0_r8)*veg_esat_bytemp(i)
      
    do ft = 1, num_pft 
      call LeafLevelPhoto(default_can_air_press, default_can_co2_pp, default_can_o2_pp,  &
        veg_tempk(i), default_can_tempk, can_air_vpress_bytemp(i), veg_esat_bytemp(i),   &
        default_par, default_leaf_bl_cond, default_nscaler, default_dayl_fact,     &
        default_btran, ft, &
        leaf_photo_bytemp(i,ft), stomatal_conductance, net_photo, c13disc)
    end do
  end do
  
  
                                                                  !!! nscaler ----------------------------------------------------------------------------

                                                                  ! calculate photosynthesis as we scale the nscaler and hold everything else constant
  do i = 1, num_nscaler
    
                                                                  ! calculate nscaler
    nscaler(i) = 0.0 + 0.05*(i-1)
    
    do ft = 1, num_pft 
      call LeafLevelPhoto(default_can_air_press, default_can_co2_pp, default_can_o2_pp,  &
        default_veg_tempk, default_can_tempk, default_can_air_vpress, default_veg_esat,  &
        default_par, default_leaf_bl_cond, nscaler(i), default_dayl_fact,          &
        default_btran, ft, &
        leaf_photo_bynscaler(i,ft), stomatal_conductance, net_photo, c13disc)
    end do
  end do
  
  
                                                                  !!! btran ----------------------------------------------------------------------------

                                                                  ! calculate photosynthesis as we scale the BTRAN and hold everything else constant
  do i = 1, num_btran
    
                                                                  ! calculate btran
    btran(i) = 0.0 + 0.05*(i-1)
    
    do ft = 1, num_pft 
      call LeafLevelPhoto(default_can_air_press, default_can_co2_pp, default_can_o2_pp,  &
        default_veg_tempk, default_can_tempk, default_can_air_vpress, default_veg_esat,  &
        default_par, default_leaf_bl_cond, default_nscaler, default_dayl_fact,     &
        btran(i), ft, &
        leaf_photo_bybtran(i,ft), stomatal_conductance, net_photo, c13disc)
    end do
  end do
  
  
                                                                  !!! Write output -----------------------------------------------------------------------
    
  call WriteLeafPhotosynthesis(out_file, num_temp, num_pft, num_PAR, num_RH, num_co2,    &
    num_nscaler, num_btran, veg_tempk, par, RH, co2_pp, nscaler, btran,                  &
    can_air_vpress_bytemp, can_air_vpress_byRH, veg_esat_bytemp, veg_esat_byRH,          &
    leaf_photo_bytemp, leaf_photo_byPAR, leaf_photo_byRH, leaf_photo_byCO2,              &
    leaf_photo_bynscaler, leaf_photo_bybtran) 
    
    
                                                                  !!! Wrap up ----------------------------------------------------------------------------

                                                                  ! deallocate arrays
  if (allocated(veg_tempk)) deallocate(veg_tempk)
  if (allocated(PAR)) deallocate(PAR)
  if (allocated(RH)) deallocate(RH)
  if (allocated(co2_pp)) deallocate(co2_pp)
  if (allocated(can_air_vpress_bytemp)) deallocate(can_air_vpress_bytemp)
  if (allocated(leaf_photo_bytemp)) deallocate(leaf_photo_bytemp)
  if (allocated(leaf_photo_byPAR)) deallocate(leaf_photo_byPAR)
  if (allocated(leaf_photo_byRH)) deallocate(leaf_photo_byRH)
  if (allocated(leaf_photo_byco2)) deallocate(leaf_photo_byco2)
  
end program FatesTestLeafPhoto

                                                                  ! ----------------------------------------------------------------------------------------

subroutine LeafLevelPhoto(can_air_press, can_co2_pp, can_o2_pp, veg_tempk, can_tempk,    &
    can_vpress, veg_esat, par, leaf_bl_conductance, nscaler, dayl_fact, btran, ft,        &
    gross_photo, stomatal_cond, net_photo, c13disc) 
                                                                  !
                                                                  ! DESCRIPTION:
                                                                  ! Drives and individual photosynthesis calculation
                                                                  !
  use FatesConstantsMod,          only : r8 => fates_r8
  use FatesTestPhotosynthesisMod,  only : CalcVaporPressure
  use LeafBiophysicsMod,           only : LeafLayerBiophysicalRates
  use LeafBiophysicsMod,           only : LeafLayerPhotosynthesis
  use LeafBiophysicsMod,           only : GetCanopyGasParameters
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use FatesParameterDerivedMod,    only : param_derived
  
                                                                  ! ARGUMENTS:
  real(r8), intent(in)  :: can_air_press                          ! canopy air pressure [Pa]
  real(r8), intent(in)  :: can_co2_pp                             ! canopy air CO2 partial pressure [Pa]
  real(r8), intent(in)  :: can_o2_pp                              ! canopy air O2 partial pressure [Pa]
  real(r8), intent(in)  :: veg_tempk                              ! leaf temperature [K]
  real(r8), intent(in)  :: can_tempk                              ! canopy air temperature [K]
  real(r8), intent(in)  :: can_vpress                             ! canopy air vapor pressure [Pa]
  real(r8), intent(in)  :: veg_esat                               ! saturation vapor pressure at leaf surface [Pa]
  real(r8), intent(in)  :: par                                    ! PAR [umol/m2/s]
  real(r8), intent(in)  :: leaf_bl_conductance                    ! leaf boundary layer conductance [umol/m2/s]
  real(r8), intent(in)  :: nscaler                                ! scaler for leaf nitrogen [0-1]
  real(r8), intent(in)  :: dayl_fact                              ! scaler for day length factor [0-1]
  real(r8), intent(in)  :: btran                                  ! scaler for BTRAN [0-1]
  integer,  intent(in)  :: ft                                     ! plant functional type index

  real(r8), intent(out) :: gross_photo                            ! gross photosynthesis (umol CO2/m2/s)
  real(r8), intent(out) :: stomatal_cond                          ! leaf stomatal conductance (umol H2O/m2/s)
  real(r8), intent(out) :: net_photo                              ! net photosynthesis (umol CO2/m2/s)
  real(r8), intent(out) :: c13disc                                ! carbon 13 in newly assimilated carbon
  
  
                                                                  ! LOCALS:
  real(r8) :: mm_kco2                                             ! Michaelis-Menten constant for CO2 [Pa]
  real(r8) :: mm_ko2                                              ! Michaelis-Menten constant for O2 [Pa]
  real(r8) :: co2_compensation_pt                                 ! CO2 compensation point [Pa]

  real(r8) :: vcmax                                               ! maximum rate of carboxylation [umol CO2/m2/s]
  real(r8) :: jmax                                                ! maximum electron transport rate [umol electrons/m2/s]
  real(r8) :: co2_rcurve_islope                                   ! initial slope of CO2 response curve (C4 plants)
  real(r8) :: leaf_maintenance_resp                               ! leaf maintenance respiration [umol CO2/m2/s]


  real(r8), parameter :: leaf_psi_dummy = -999.9_r8               ! Dummy for leaf psi Mpa (not used)
  integer,  parameter :: c4_grass_pft = 12
  real(r8), parameter :: leaf_area_dummy = 1._r8                  ! This is only used in LeafLayerPhotosynthesis
                                                                  ! to test against zero, in the trivial
                                                                  ! case calculations are bypassed if there
                                                                  ! is no leaf
  
  
  ! calculate canopy gas parameters
  call GetCanopyGasParameters(can_air_press, can_o2_pp, veg_tempk, &
       mm_kco2,mm_ko2,co2_compensation_pt)

  
  ! calculate leaf biophysical rates
  call LeafLayerBiophysicalRates(ft, EDPftvarcon_inst%vcmax25top(ft,1),             &
    param_derived%jmax25top(ft,1),param_derived%kp25top(ft,1), nscaler, veg_tempk,      &
    dayl_fact, can_tempk, can_tempk, btran, vcmax, jmax, co2_rcurve_islope)
  
  if (ft == c4_grass_pft) then
    leaf_maintenance_resp = 0.025_r8*EDPftvarcon_inst%vcmax25top(ft,1)
  else 
    leaf_maintenance_resp = 0.015_r8*EDPftvarcon_inst%vcmax25top(ft,1)
  end if 
  
  ! calculate leaf-level photosynthesis
  call LeafLayerPhotosynthesis(par, leaf_area_dummy, ft, vcmax, jmax,    &
    co2_rcurve_islope, veg_tempk, veg_esat, can_air_press, can_co2_pp, can_o2_pp, btran, &
    leaf_bl_conductance, can_vpress, mm_kco2,        &
    mm_ko2, co2_compensation_pt, leaf_maintenance_resp, leaf_psi_dummy,    &
    gross_photo, stomatal_cond, net_photo, c13disc)

  
end subroutine LeafLevelPhoto

! ----------------------------------------------------------------------------------------

subroutine WriteLeafPhotosynthesis(out_file, num_temp, num_pft, num_par, num_RH,         &
  num_co2, num_nscaler, num_btran, veg_tempk, par, RH, co2, nscaler, btran,              &
  can_air_vpress_bytemp, can_air_vpress_byRH, veg_esat_bytemp, veg_esat_byRH,            &
  leaf_photo_bytemp, leaf_photo_byPAR, leaf_photo_byRH, leaf_photo_byCO2,                &
  leaf_photo_bynscaler, leaf_photo_bybtran)
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

  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file                 ! output file name
  integer,          intent(in) :: num_temp                 ! size of temperature array
  integer,          intent(in) :: num_pft                  ! number of pfts
  integer,          intent(in) :: num_par                  ! size of PAR array
  integer,          intent(in) :: num_RH                   ! size of RH array
  integer,          intent(in) :: num_co2                  ! size of co2 array
  integer,          intent(in) :: num_nscaler              ! size of nscaler array
  integer,          intent(in) :: num_btran                ! size of btran array
  real(r8),         intent(in) :: veg_tempk(:)             ! vegetation temperature [K]
  real(r8),         intent(in) :: par(:)                   ! PAR [umol/m2/s]
  real(r8),         intent(in) :: RH(:)                    ! relative humidity [%]
  real(r8),         intent(in) :: co2(:)                   ! CO2 partial pressure [Pa]
  real(r8),         intent(in) :: nscaler(:)               ! scaler for leaf N
  real(r8),         intent(in) :: btran(:)                 ! scaler for BTRAN
  real(r8),         intent(in) :: can_air_vpress_bytemp(:) ! canopy air vapor pressure by temperature [Pa]
  real(r8),         intent(in) :: can_air_vpress_byRH(:)   ! canopy air vapor pressure by RH [Pa]
  real(r8),         intent(in) :: veg_esat_bytemp(:)       ! vegetation saturation vapor pressure by temperature [Pa]
  real(r8),         intent(in) :: veg_esat_byRH(:)         ! vegetation saturation vapor pressure by RH [Pa]
  real(r8),         intent(in) :: leaf_photo_bytemp(:,:)   ! net leaf photosynthesis [umol CO2/m2/s]
  real(r8),         intent(in) :: leaf_photo_byPAR(:,:)    ! net leaf photosynthesis [umol CO2/m2/s]
  real(r8),         intent(in) :: leaf_photo_byRH(:,:)     ! net leaf photosynthesis [umol CO2/m2/s]
  real(r8),         intent(in) :: leaf_photo_byCO2(:,:)    ! net leaf photosynthesis [umol CO2/m2/s]
  real(r8),         intent(in) :: leaf_photo_bynscaler(:,:) ! net leaf photosynthesis [umol CO2/m2/s]
  real(r8),         intent(in) :: leaf_photo_bybtran(:,:)  ! net leaf photosynthesis [umol CO2/m2/s]
  

  ! LOCALS:
  integer, allocatable :: pft_indices(:) ! array of pft indices to write out
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=12)    :: dim_names(7)   ! dimension names
  integer              :: dimIDs(7)      ! dimension IDs
  integer              :: tempID, pftID  ! variable IDs for dimensions
  integer              :: parID, presstempID
  integer              :: pressRHID, esattempID
  integer              :: esatRHID, btranID
  integer              :: nscalerID
  integer              :: RHID, CO2ID
  integer              :: phototempID
  integer              :: photoPARID
  integer              :: photoRHID
  integer              :: photoCO2ID
  integer              :: photonscalerID
  integer              :: photobtranID
  
  ! create pft indices
  allocate(pft_indices(num_pft))
  do i = 1, num_pft
    pft_indices(i) = i
  end do

  ! dimension names
  dim_names = [character(len=12) :: 'temperature', 'par', 'RH', 'CO2', 'nscaler',        &
    'btran', 'pft']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/num_temp, num_par, num_RH, num_co2,             &
    num_nscaler, num_btran, num_pft/), 7, dimIDs)

  ! register temperature
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: 'K', 'vegetation temperature'], 2, tempID)
    
  ! register par
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_double,                       &
    [character(len=20)  :: 'units', 'long_name'],                                      &
    [character(len=150) :: 'umol m-2 s-1', 'photosynthestically active radiation'], 2, &
    parID)
    
  ! register RH
    call RegisterVar(ncid, dim_names(3), dimIDs(3:3), type_double, &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '%', 'relative humidity'], 2,           &
    RHID)
    
  ! register CO2
  call RegisterVar(ncid, dim_names(4), dimIDs(4:4), type_double, &
    [character(len=20)  :: 'units', 'long_name'],                &
    [character(len=150) :: 'Pa', 'CO2 partial pressure'], 2,     &
    CO2ID)
    
  ! register nscaler
  call RegisterVar(ncid, dim_names(5), dimIDs(5:5), type_double,      &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'scaler for leaf N'], 2, nscalerID)
    
  ! register btran
  call RegisterVar(ncid, dim_names(6), dimIDs(6:6), type_double,      &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'scaler for BTRAN'], 2, btranID)
    
  ! register pft
  call RegisterVar(ncid, dim_names(7), dimIDs(7:7), type_int,      &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'plant functional type'], 2, pftID)
    
  ! register vapor pressure (by temperature)
  call RegisterVar(ncid, 'can_air_vpress_bytemp', dimIDs(1:1), type_double,         &
    [character(len=20)  :: 'units', 'long_name'],                                   &
    [character(len=150) :: 'Pa', 'vapor pressure  - with changing temperature'], 2, &
    presstempID)
    
  ! register saturation vapor pressure (by temperature)
    call RegisterVar(ncid, 'veg_esat_bytemp', dimIDs(1:1), type_double,                        &
    [character(len=20)  :: 'units', 'long_name'],                                              &
    [character(len=150) :: 'Pa', 'saturation vapor pressure  - with changing temperature'], 2, &
    esattempID)
    
  ! register leaf photosynthesis (by temperature)
  call RegisterVar(ncid, 'leaf_photo_bytemp', (/dimIDs(1), dimIDs(7)/), type_double,                                      &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                                                          &
    [character(len=150) :: 'pft temperature', 'umol CO2 m-2 s-1', 'net leaf photosynthesis - with changing temperature'], &
    3, phototempID)
    
  ! register leaf photosynthesis (by PAR)
    call RegisterVar(ncid, 'leaf_photo_byPAR', (/dimIDs(2), dimIDs(7)/), type_double,                     &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                                          &
    [character(len=150) :: 'pft par', 'umol CO2 m-2 s-1', 'net leaf photosynthesis - with changing PAR'], &
    3, photoPARID)
    
  ! register leaf photosynthesis (by CO2)
  call RegisterVar(ncid, 'leaf_photo_byCO2', (/dimIDs(4), dimIDs(7)/), type_double,                       &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                                          &
    [character(len=150) :: 'pft CO2', 'umol CO2 m-2 s-1', 'net leaf photosynthesis - with changing CO2'], &
    3, photoCO2ID)
    
  ! register vapor pressure (by RH)
  call RegisterVar(ncid, 'can_air_vpress_byRH', dimIDs(3:3), type_double,    &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],             &
    [character(len=150) :: 'RH', 'Pa', 'vapor pressure - with changing RH'], &
    3, pressRHID)
    
  ! register saturation vapor pressure (by RH)
    call RegisterVar(ncid, 'veg_esat_byRH', dimIDs(3:3), type_double,                   &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                        &
    [character(len=150) :: 'RH', 'Pa', 'saturation vapor pressure - with changing RH'], &
    3, esatRHID)
  
  ! register leaf photosynthesis (by RH)
    call RegisterVar(ncid, 'leaf_photo_byRH', (/dimIDs(3), dimIDs(7)/), type_double,                    &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                                        &
    [character(len=150) :: 'pft RH', 'umol CO2 m-2 s-1', 'net leaf photosynthesis - with changing RH'], &
    3, photoRHID)
    
  ! register leaf photosynthesis (by nscaler)
    call RegisterVar(ncid, 'leaf_photo_bynscaler', (/dimIDs(5), dimIDs(7)/), type_double,                       &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                                          &
    [character(len=150) :: 'pft nscaler', 'umol CO2 m-2 s-1', 'net leaf photosynthesis - with changing nscaler'], &
    3, photonscalerID)
    
  ! register leaf photosynthesis (by btran)
    call RegisterVar(ncid, 'leaf_photo_bybtran', (/dimIDs(6), dimIDs(7)/), type_double,                       &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],                                          &
    [character(len=150) :: 'pft btran', 'umol CO2 m-2 s-1', 'net leaf photosynthesis - with changing btran'], &
    3, photobtranID)
  
  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, tempID, veg_tempk(:))
  call WriteVar(ncid, parID, par(:))
  call WriteVar(ncid, RHID, RH(:))
  call WriteVar(ncid, CO2ID, co2(:))
  call WriteVar(ncid, nscalerID, nscaler(:))
  call WriteVar(ncid, btranID, btran(:))
  call WriteVar(ncid, pftID, pft_indices(:))
  call WriteVar(ncid, presstempID, can_air_vpress_bytemp(:))
  call WriteVar(ncid, pressRHID, can_air_vpress_byRH(:))
  call WriteVar(ncid, esattempID, veg_esat_bytemp(:))
  call WriteVar(ncid, esatRHID, veg_esat_byRH(:))
  call WriteVar(ncid, phototempID, leaf_photo_bytemp(:,:))
  call WriteVar(ncid, photoPARID, leaf_photo_byPAR(:,:))
  call WriteVar(ncid, photoRHID, leaf_photo_byRH(:,:))
  call WriteVar(ncid, photoCO2ID, leaf_photo_byCO2(:,:))
  call WriteVar(ncid, photonscalerID, leaf_photo_bynscaler(:,:))
  call WriteVar(ncid, photobtranID, leaf_photo_bybtran(:,:))
  
  ! close the file
  call CloseNCFile(ncid)

end subroutine WriteLeafPhotosynthesis

