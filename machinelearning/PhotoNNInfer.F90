program PhotoNNInfer

  ! Import precision info from iso
  use, intrinsic :: iso_fortran_env, only : sp => real32
  use, intrinsic :: iso_fortran_env, only : dp => real64

  ! Import our library for interfacing with PyTorch
  use ftorch, only : torch_model, torch_tensor, torch_kCPU, torch_delete, &
       torch_tensor_from_array, torch_model_load, torch_model_forward

  ! Import our tools module for testing utils
  use ftorch_test_utils, only : assert_allclose

  use FatesConstantsMod, only  : tfrz => t_water_freeze_k_1atm
  use LeafBiophysicsMod, only  : QSat !(veg_tempk, can_press, qsat_alt, veg_esat)
  use LeafBiophysicsMod, only  : GetCanopyGasParameters
  use LeafBiophysicsMod, only  : LeafLayerBiophysicalRates
  use LeafBiophysicsMod, only  : LeafLayerPhotosynthesis
  use LeafBiophysicsMod, only  : LeafLayerMaintenanceRespiration_Ryan_1991
  use LeafBiophysSuppMod, only : AllocLeafParam
  use LeafBiophysSuppMod, only : SetLeafParam

  implicit none

  ! Set working precision for reals
  integer, parameter :: wp = sp
  integer :: num_args, ix
  character(len=128), dimension(:), allocatable :: args

  ! Set up Fortran data structures
  real(wp), dimension(13), target :: in_data
  real(wp), dimension(13), target :: in_data_norm
  real(wp), dimension(2),  target :: out_data
  real(wp), dimension(2),  target :: out_data_norm
  real(wp), dimension(2),  target :: expected

  ! Set up Torch data structures
  ! The net, a vector of input tensors (in this case we only have one), and the output tensor
  ! The size of 1 indicates a batch size of 1 (which is all we need for an inference step)
  ! Apparently the in_data and out_data will automatically allocate the correct size in
  ! the in_tensors

  type(torch_model) :: model
  type(torch_tensor), dimension(1) :: in_tensors
  type(torch_tensor), dimension(1) :: out_tensors

  ! Normalization

  ! c3psn_modelsd_i13_L16-Re-L32-Re-L32_c20250807-1112.pt
  
  real(wp), dimension(13) :: in_mean = [5.7742e+02, 1.6995e+01, 1.7551e+01, &
       2.9115e-01, 3.0815e+02, 3.9990e+02, &
       6.4835e+03, 2.7501e+06, 3.8906e+03, &
       1.8506e+02, 5.0139e+04, 7.9321e+00, &
       5.2570e-01]
  
  real(wp), dimension(13) :: in_std =[4.0494e+02, 2.9277e+01, 3.4885e+01, &
       3.7018e-01, 1.0610e+01, 1.4143e+02, &
       3.5785e+03, 1.5914e+06, 2.9987e+03, &
       1.6973e+02, 2.3459e+04, 3.8481e+00, &
       2.0702e-01]

  real(wp), dimension(2) :: out_mean = [1.8664e+00, 7.2744e+04]
  real(wp), dimension(2) :: out_std  = [4.6094e+00, 4.5607e+05]


  
  ! Environmental ranges for sampling a solution
  real(dp) :: leaf_tempc_min,leaf_tempc_max
  real(dp) :: rh_max,rh_min
  real(dp) :: co2_max,co2_min
  real(dp) :: par_abs_max,par_abs_min
  real(dp) :: gb_max,gb_min
  real(dp) :: btran_max,btran_min
  real(dp) :: vcmax25t_max,vcmax25t_min
  real(dp) :: air_press_max,air_press_min

  ! Local environmental variables
  real(dp) :: lnc_top      ! leaf N concentration top [gN/m2]
  real(dp) :: leaf_tempc   ! leaf temp C
  real(dp) :: rh           ! RH as a fraction
  real(dp) :: co2_ppm      ! atm CO2 concentration
  real(dp) :: par_abs      ! Absorbed PAR [umol/m2/s]
  real(dp) :: gb_umol      ! BL conductance [umol/m2/s]
  real(dp) :: btran        ! Transp wetness factor [0-1]
  real(dp) :: vcmax25t     ! VCmax @ 25C top of canopy [umol/m2/s]
  real(dp) :: jmax25t      ! JMax @ 25C top of canopy [umol/m2/s]
  real(dp) :: leaf_tempk  ! Leaf temperature in Kelvin
  real(dp) :: o2_pa        ! Atmospheric O2 partial pressure [Pa]
  real(dp) :: co2_pa       ! Atmospheric CO2 partial pressure [Pa]
  real(dp) :: co2_cpoint_pa ! CO2 compensation point [Pa]
  real(dp) :: mm_kco2       ! MM parameter for co2
  real(dp) :: mm_ko2        ! MM parameter for o2
  real(dp) :: qsat_sh       ! Saturation spec hum [g/kg?](unused)
  real(dp) :: qsat_pa       ! Saturation vapor pressure [Pa]
  real(dp) :: q_pa          ! Actual vapor pressure [Pa]
  real(dp) :: lmr           ! Leaf dark (maintenance) resp [umol/m2/s]
  real(dp) :: vcmax
  real(dp) :: jmax
  real(dp) :: kp_dummy
  real(dp) :: gs0
  real(dp) :: gs1
  real(dp) :: gs2

  ! Model output
  real(dp) :: agross
  real(dp) :: gs
  real(dp) :: anet
  real(dp) :: c13disc
  real(dp) :: ci
  integer  :: solve_iter
  
  integer, parameter :: ntests = 50
  integer :: i
  
  ! 1 standard atmosphere in [Pa]
  real(dp), parameter :: can_press_stdatm_pa = 101325.0
  ! Typical O2 concentration in atmosphere 209k ppm
  real(dp), parameter :: o2_ppp = 0.2095

  ! Simple conversion, number of micro-moles in a mole
  real(dp), parameter :: umol_per_mol = 1.e6
  real(dp), parameter :: mol_per_umol = 1.e-6

  ! Plant Traits and parameter constants (MATCH WITH TRAINING RUN!!!)
  real(dp), parameter :: fates_stoich_nitr = 0.033   ! N/C ratio of leaf [gN/gC]
  real(dp), parameter :: fates_leaf_slatop = 0.012   ! Specific Leaf Area at the top m2/gC
  real(dp), parameter :: kp25t_dummy = -9.e32_dp
  real(dp), parameter :: dayl_factor_full = 1._dp
  real(dp), parameter :: t_growth_dummy = 9.e32_dp
  real(dp), parameter :: t_home_dummy   = 9.e32_dp
  real(dp), parameter :: ci_tol = 0.1_dp
  
  real(dp), parameter :: nscaler_top = 1.0
  integer, parameter :: pft1 = 1   ! Token pft index for single type
  integer, parameter :: medlyn_model     = 2
  integer, parameter :: btran_on_gs_gs2  = 4
  integer, parameter :: btran_on_ag_vcmax_jmax = 2
  integer, parameter :: c3_path = 1
  integer, parameter :: daylen_on = 1
  integer, parameter :: fvcb1980 = 1
  integer, parameter :: jb2021 = 2
  integer, parameter :: net_assim_model = 1
  integer, parameter :: no_tempsense = 0

  
  ! Flag for testing
  logical :: test_pass

  interface
     function Normalize1DArray(x_in,x_mean,x_std) result(x_out)
       use, intrinsic :: iso_fortran_env, only : sp => real32
       implicit none
       real(sp), dimension(:), intent(in) :: x_in
       real(sp), dimension(:), intent(in) :: x_mean
       real(sp), dimension(:), intent(in) :: x_std
       real(sp), dimension(size(x_in)) :: x_out
     end function Normalize1DArray
     function DeNormalize1DArray(x_out,x_mean,x_std) result(x_in)
       use, intrinsic :: iso_fortran_env, only : sp => real32
       ! This is just the reverse process
       implicit none
       real(sp), dimension(:), intent(in) :: x_out
       real(sp), dimension(:), intent(in) :: x_mean
       real(sp), dimension(:), intent(in) :: x_std
      real(sp), dimension(size(x_out)) :: x_in
     end function DeNormalize1DArray
     function RandSamp(vmin,vmax) result(vout)
       use, intrinsic :: iso_fortran_env, only : dp => real64
       implicit none
       real(dp) :: vmin
       real(dp) :: vmax
       real(dp) :: vout
     end function RandSamp
  end interface

  ! Get TorchScript model file as a command line argument
  num_args = command_argument_count()
  allocate(args(num_args))
  do ix = 1, num_args
     call get_command_argument(ix,args(ix))
  end do

  ! Leaf temperature ranges [C]
  leaf_tempc_min = 20.0
  leaf_tempc_max = 50.0

  ! Relative Humidity Ranges
  rh_max = 1.00
  rh_min = 0.2

  ! CO2 concentration ranges [ppm]
  co2_max = 600.0
  co2_min = 200.0

  ! Absorbed PAR ranges [umol/m2/s]
  par_abs_min = 1.0
  par_abs_max = 800.0

  ! Boundary Conductance ranges [umol/m2/s]
  gb_min =  500000.0   ! Lower limit imposed by CLM/ELM 0.5 mol/m2/s
  gb_max = 5000000.0   ! 50% larger than  Roughly largest
  ! values seen at BCI (which are 2.5mol/m2/s)

  ! btran ranges
  btran_min = 0.1
  btran_max = 1.0

  ! vcmax25top ranges
  vcmax25t_min = 2.0
  vcmax25t_max = 125.0

  ! Leaf Nitrogen Concentration at canopy top gN/m2
  lnc_top  = fates_stoich_nitr/fates_leaf_slatop

  call AllocLeafParam(1)
  
  call SetLeafParam(real(fvcb1980,dp),0,'fates_electron_transport_model')
  call SetLeafParam(real(daylen_on,dp),0,'fates_daylength_factor_switch')
  call SetLeafParam(real(medlyn_model,dp),0,'fates_leaf_stomatal_model')
  call SetLeafParam(real(net_assim_model,dp),0,'fates_leaf_stomatal_assim_model')
  call SetLeafParam(real(no_tempsense,dp),0,'fates_leaf_photo_tempsens_model')
  call SetLeafParam(real(c3_path,dp),1,'fates_leaf_c3psn')
  call SetLeafParam(real(btran_on_gs_gs2,dp),1,'fates_leaf_stomatal_btran_model')
  call SetLeafParam(real(btran_on_ag_vcmax_jmax,dp),1,'fates_leaf_agross_btran_model')
  call SetLeafParam(0.15_dp,1,'fates_leaf_fnps')
  call SetLeafParam(4.1_dp,1,'fates_leaf_stomatal_slope_medlyn')
  call SetLeafParam(10000._dp,1,'fates_leaf_stomatal_intercept')
  call SetLeafParam(2.525e-6_dp,1,'fates_maintresp_leaf_ryan1991_baserate')
  call SetLeafParam(1.756_dp,1,'fates_maintresp_leaf_atkin2017_baserate')
  call SetLeafParam(65330._dp,1,'fates_leaf_vcmaxha')
  call SetLeafParam(43540._dp,1,'fates_leaf_jmaxha')
  call SetLeafParam(149250._dp,1,'fates_leaf_vcmaxhd')
  call SetLeafParam(149250._dp,1,'fates_leaf_jmaxhd')
  call SetLeafParam(485._dp,1,'fates_leaf_vcmaxse')
  call SetLeafParam(495._dp,1,'fates_leaf_jmaxse')

  call torch_model_load(model, args(1), torch_kCPU) !,requires_grad=.false.,is_training=.false.)

  print*,"MOD Ag, NN Ag, Mod gs, NN gs"
  
  do i = 1, ntests

     leaf_tempc = RandSamp(leaf_tempc_min,leaf_tempc_max)
     rh         = RandSamp(rh_min,rh_max)
     co2_ppm    = RandSamp(co2_min,co2_max)
     par_abs    = RandSamp(par_abs_min,par_abs_max)
     gb_umol    = RandSamp(gb_min,gb_max)
     btran      = RandSamp(btran_min,btran_max)
     vcmax25t   = RandSamp(vcmax25t_min,vcmax25t_max)
     
     ! Conversions
     jmax25t    = 1.67 * vcmax25t
     leaf_tempk = leaf_tempc + tfrz
     o2_pa      = o2_ppp * can_press_stdatm_pa
     co2_pa     = co2_ppm * 1.e-6 * can_press_stdatm_pa

     call QSat(leaf_tempk, can_press_stdatm_pa, qsat_sh, qsat_pa)

     q_pa = qsat_pa * rh
     
     call GetCanopyGasParameters(can_press_stdatm_pa, &
          o2_pa,     &
          leaf_tempk, &
          mm_kco2,   &
          mm_ko2,    &
          co2_cpoint_pa)

     call LeafLayerMaintenanceRespiration_Ryan_1991(lnc_top, &
          nscaler_top,       &
          pft1,            &
          leaf_tempk,     &
          lmr)

     !print*,"LMR: ",lnc_top,nscaler_top,pft1,leaf_tempk,lmr

     call LeafLayerBiophysicalRates(pft1,                 &
          vcmax25t, jmax25t, kp25t_dummy,                 &
          nscaler_top, leaf_tempk,                        &
          dayl_factor_full, t_growth_dummy, t_home_dummy, &
          btran, vcmax, jmax, kp_dummy, &
          gs0, gs1, gs2)

     in_data(1) = par_abs
     in_data(2) = vcmax
     in_data(3) = jmax
     in_data(4) = gs2
     in_data(5) = leaf_tempk
     in_data(6) = co2_ppm
     in_data(7) = qsat_pa
     in_data(8) = gb_umol
     in_data(9) = q_pa
     in_data(10) = mm_kco2
     in_data(11) = mm_ko2
     in_data(12) = co2_cpoint_pa
     in_data(13) = lmr

     in_data_norm = Normalize1DArray(in_data,in_mean,in_std)
     
     !print*,"----"
     !print*,"in: par_abs = ",par_abs,"vcmax = ",vcmax,"jmax = ",jmax
     !print*,"in: gs0 = ",gs0,"gs2 = ",gs2,"leaf_tempk = ",leaf_tempk
     !print*,"in: can_press = ",can_press_stdatm_pa,"co2_pa = ",co2_pa,"o2_pa = ",o2_pa
     !print*,"in: qsat_pa = ",qsat_pa,"gb_umol = ",gb_umol,"q_pa = ",q_pa
     !print*,"in: mm_kco2 = ",mm_kco2,"mm_ko2 = ",mm_ko2,"lmr = ",lmr
     !print*,"in: ci_tol = ",ci_tol
     
     call LeafLayerPhotosynthesis(         &
          par_abs, pft1,                   &  ! in
          vcmax, jmax, kp_dummy,           &  ! in
          gs0, gs1, gs2,                   &  ! in
          leaf_tempk, can_press_stdatm_pa, &  ! in
          co2_pa, o2_pa,                   &  ! in
          qsat_pa, gb_umol, q_pa,          &  ! in
          mm_kco2, mm_ko2,                 &  ! in
          co2_cpoint_pa,                   &  ! in
          lmr,                             &  ! in
          ci_tol,                          &  ! in
          agross,                          &  ! out
          gs,                              &  ! out
          anet,                            &  ! out
          c13disc,                         &  ! out
          ci,                              &  ! out
          solve_iter)                         ! out
     
     call torch_tensor_from_array(in_tensors(1), in_data_norm, torch_kCPU)
     call torch_tensor_from_array(out_tensors(1), out_data_norm, torch_kCPU)
     call torch_model_forward(model, in_tensors, out_tensors)

     out_data = DeNormalize1DArray(out_data_norm,out_mean,out_std)

     
     
     print*,agross,out_data(1),gs,out_data(2)
     !print*,"nn    out: ",out_data(1), out_data(2)
     
  end do

  !in_data = [4.60000000e+01, 3.28302733e-06, 1.37791026e-04, 1.00000000e-02, &
  !     2.23150000e+02, 3.00000000e+02, 3.93892034e+00, 5.00000000e+05, &
  !     2.16640618e+00, 8.63170186e-04, 2.03363077e+02, 2.57157011e-02, &
  !     1.52429652e-03]

  !write (*,*) in_data
  ! Initialise data
  ! Create Torch input/output tensors from the above arrays
  ! The (1) refers to a batch size of 1. I believe thta 

  write (*,*) out_data(:)
  
  ! Check output tensor matches expected value
  !   expected = [0.0_wp, 2.0_wp, 4.0_wp, 6.0_wp, 8.0_wp]
  !expected = [3.26423242e-06, 1.00000000e+04]
  expected = [3.27047356e-06, 1.00000000e+04]

  !!test_pass = assert_allclose(out_data, expected, test_name="SimpleNet", rtol=1e-5)

  ! Cleanup
  !!call torch_delete(model)
  !!call torch_delete(in_tensors)
  !!call torch_delete(out_tensors)

  !!if (.not. test_pass) then
  !!   stop 999
  !!end if

  write (*,*) "SimpleNet example ran successfully"

end program PhotoNNInfer

function Normalize1DArray(x_in,x_mean,x_std) result(x_out)
  use, intrinsic :: iso_fortran_env, only : sp => real32
  implicit none

  real(sp), dimension(:), intent(in) :: x_in
  real(sp), dimension(:), intent(in) :: x_mean
  real(sp), dimension(:), intent(in) :: x_std
  integer :: i
  real(sp), dimension(size(x_in)) :: x_out

  do i = 1,size(x_in,dim=1)
     x_out(i) = (x_in(i) - x_mean(i)) / x_std(i)
  end do

end function Normalize1DArray

function DeNormalize1DArray(x_out,x_mean,x_std) result(x_in)
  use, intrinsic :: iso_fortran_env, only : sp => real32
  ! This is just the reverse process
  implicit none
  real(sp), dimension(:), intent(in) :: x_out
  real(sp), dimension(:), intent(in) :: x_mean
  real(sp), dimension(:), intent(in) :: x_std
  integer :: i
  real(sp), dimension(size(x_out)) :: x_in

  do i = 1,size(x_out,dim=1)
     x_in(i) = x_out(i)*x_std(i)+x_mean(i)
  end do

end function DeNormalize1DArray

function RandSamp(vmin,vmax) result(vout)
  use, intrinsic :: iso_fortran_env, only : dp => real64
  implicit none
  real(dp) :: vmin
  real(dp) :: vmax
  real(dp) :: vout
  real(dp) :: x
  call random_number(x)
  vout = vmin + x*(vmax-vmin)
end function RandSamp

