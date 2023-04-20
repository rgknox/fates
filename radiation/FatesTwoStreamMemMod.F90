Module FatesTwoStreamMemMod

  ! ---------------------------------------------------------------------------
  ! This module is a space that holds data that defines how
  ! FATES in particular uses the two-stream radiation scheme.
  ! Alternatively, the TwoStreamMLPEMod is more agnostic.
  ! For instance, TwoStreamMLPEMod makes no assumptions about
  ! which or how many broad bands are used
  ! ---------------------------------------------------------------------------
  

  ! TODO: we use this cp_maxSWb only because we have a static array q(size=2) of
  ! land-ice abledo for vis and nir.  This should be a parameter, which would
  ! get us on track to start using multi-spectral or hyper-spectral (RGK 02-2017)
  
  integer, parameter, public :: num_swb = 2     ! Number of shortwave bands we use
                                                ! This needs to match what is used in the host model
                                                ! This is visible (1) and near-infrared (2)
  
  integer, parameter, public :: ivis = 1        ! This is the array index for short-wave
                                                ! radiation in the visible spectrum, as expected
                                                ! in boundary condition files and parameter
                                                ! files.  This will be compared with 
                                                ! the HLM's expectation in FatesInterfaceMod

  integer, parameter, public :: inir = 2        ! This is the array index for short-wave
                                                ! radiation in the near-infrared spectrum, as expected
                                                ! in boundary condition files and parameter
                                                ! files.  This will be compared with 
                                                ! the HLM's expectation in FatesInterfaceMod

  integer, parameter, public :: ipar = ivis     ! The photosynthetically active band
                                                ! can be approximated to be equal to the visible band

  integer, parameter  :: max_el_per_layer = 10
  
  real(r8), parameter :: init_max_vai_diff_per_elem = 0.2_r8
  
  ! albedo land ice by waveband (1=vis, 2=nir)
  real(r8), public  :: albice(num_swb)   =  (/ 0.80_r8, 0.55_r8 /)

  ! albedo land ice by waveband (1=vis, 2=nir)
  real(r8), public  :: rho_snow(num_swb) = (/ 0.80_r8, 0.55_r8 /)

  ! albedo land ice by waveband (1=vis, 2=nir)
  real(r8), public  :: tau_snow(num_swb) = (/ 0.01_r8, 0.01_r8 /)

  
end Module FatesTwoStreamMemMod
