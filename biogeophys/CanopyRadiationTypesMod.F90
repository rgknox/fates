Module CanopyRadiationTypesMod

  use FatesConstantsMod, only: r8 => fates_r8
  
  implicit none
  public

  integer,parameter :: num_swb = 2   ! Number of short-wave broadbands
                                     ! most likely 2, VIS & NIR

  integer, parameter :: num_rad_stream = 2  ! Direct and diffuse

  
  integer :: nclmax     = 3
  integer :: maxpft     = 20
  integer :: maxlevleaf = 50

  real(r8), public  :: albice(num_swb) = &       ! albedo land ice by waveband (1=vis, 2=nir)
       (/ 0.80_r8, 0.55_r8 /)
  real(r8), public  :: rho_snow(num_swb) = &     ! albedo land ice by waveband (1=vis, 2=nir)
       (/ 0.80_r8, 0.55_r8 /)
  real(r8), public  :: tau_snow(num_swb) = &     ! albedo land ice by waveband (1=vis, 2=nir)
       (/ 0.01_r8, 0.01_r8 /)


  
  
  ! Scratch arrays
  !real(r8), allocatable :: ftweight_scr(:,:,:)


  type rad_scratch
     real(r8) :: rho_layer(nclmax,maxpft,maxlevleaf,num_swb) ! Layer reflectance (0-1)
     real(r8) :: tau_layer(nclmax,maxpft,maxlevleaf,num_swb) ! Layer transmittance (0-1)
     real(r8) :: f_abs(nclmax,maxpft,maxlevleaf,num_swb)     ! Fraction absorbed
     real(r8) :: f_abs_leaf(nclmax,maxpft,maxlevleaf,num_swb)! Fraction absorbed by leaves

     ! These transmittance arrays are on the layer edges. The first edge (1) is above leaf layer
     ! 1, and the second edge is below
     real(r8) :: tr_dir_z(nclmax,maxpft,maxlevleaf+1)    ! Transmittance of direct beam radiation through a single layer
     real(r8) :: tr_dif_z(nclmax,maxpft,maxlevleaf+1)    ! Transmittance of diffuse radiation through a single layer
     
  end type rad_scratch

     
  type rad_params
     real(r8), allocatable :: rhol(:,:)     ! leaf reflectance: 1=vis, 2=nir    x pft
     real(r8), allocatable :: rhos(:,:)     ! stem reflectance: 1=vis, 2=nir    x pft
     real(r8), allocatable :: taul(:,:)     ! leaf transmittance: 1=vis, 2=nir  x pft
     real(r8), allocatable :: taus(:,:)     ! stem transmittance: 1=vis, 2=nir  x pft
     real(r8), allocatable :: xl(:)         ! leaf/stem orientation index, by pft
     real(r8), allocatable :: clumping_index(:) ! clumping index 0-1, when leaves stick together, by pft
  end type rad_params


  
  type scatter_type

     

     
  end type scatter_type

  

  
  type norman_type

     ! Information provided by the host model
     !integer               :: num_canl
     !integer               :: num_pft
     !integer, allocatable  :: num_vegl(:,:)  ! Number of vegetation scattering layers in each canopy and pft
     !real(r8)              :: fcansnow       ! Fraction of all vegetation covered by snow (uniformity assumption) 0-1


     ! Calculated by the radiation scheme
     real(r8), allocatable :: f_sun(:,:,:)   ! fraction of leaves in the sun in each (canopy-layer x  pft x leaf-layer)

     
  end type norman_type

  

contains

  

  
  
end Module CanopyRadiationTypesMod
