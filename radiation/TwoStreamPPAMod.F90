Module TwoStreamPPAMod
  
  ! This module holds the routines to calculate two-tream
  ! radiation scattering over a canopy represented with
  ! an imperfectly plastic approximation. In summary,
  ! there may various canopy layers, in each canopy layer,
  ! plant media for different functional types are grouped
  ! so that they inhabit their own exclusive footprint
  ! within the layer. Within these exclusive functional
  ! columns, there are further sub-layer discretizations,
  ! which are organized by top-down integrated vegetation
  ! area index.
  !
  ! Note that there is a separate allocation and call
  ! sequence for each broad band.  In other words, the
  ! two_stream_type is instantiated for each broad band.
  !
  ! 
  !
  ! Assumptions: band index 1 = visible (vis)
  !                         2 = near infrared (nir)
  !                         3 = thermal
  !

  integer, parameter :: r8 = selected_real_kind(12)
  
  
  type rad_params_type

     ! From the parameter file
     real(r8), allocatable :: rhol(:,:)         ! leaf material reflectance:   (band x pft)
     real(r8), allocatable :: rhos(:,:)         ! stem material reflectance:   (band x pft)
     real(r8), allocatable :: taul(:,:)         ! leaf material transmittance: (band x pft)
     real(r8), allocatable :: taus(:,:)         ! stem material transmittance: (band x pft)
     real(r8), allocatable :: xl(:)             ! leaf/stem orientation (pft)
     real(r8), allocatable :: clumping_index(:) ! clumping index 0-1, when
                                                ! leaves stick together (pft)

     ! Derived parameters
     real(r8), allocatable :: phi1(:)        ! intermediate term for kd and kb
     real(r8), allocatable :: phi2(:)        ! intermediate term for kd and kb
     real(r8), allocatable :: kd_leaf(:)     ! Mean optical depth per unit area leaves in diffuse
     real(r8), allocatable :: kd_stem(:)     ! Mean optical depth per unit area stems in diffuse
     real(r8), allocatable :: om_snow(:)     ! Snow material scattering albedo for snow (band)
     real(r8), allocatable :: om_leaf(:,:)   ! Leaf material scattering albedo (band x pft)
     real(r8), allocatable :: om_stem(:,:)   ! Stem material scattering albedo (band x pft)
     
  end type rad_params_type

  type(rad_params_type) :: rad_params

  ! This type describes the scattering elements, these are values
  ! that need to be defined by the ecosystem model, somewhat of
  ! an input to the solver. Since this is a Perfect Plasticity Approximation
  ! enabled system, we partition the scattering media into "columns" and "layers"
  ! Layers are canopy layers, think understory, mid-story and upper canopy. Columns
  ! are divisions of horizontal space, ie literal columns of space. The current
  ! implementation limits this space to media that has uniform scattering coefficients.
  ! So there could not be different PFTs in the same column, because they would undoubtedly
  ! have different joint scattering coefficients at different height levels in
  ! the column.  Therefore, every column is connected with a PFT.
  
  type scel_type
     integer  :: n_col
     integer, allocatable  :: col_pft(:)   ! pft index
     real(r8), allocatable :: col_area(:)  ! m2 col/m2 ground
     real(r8), allocatable :: col_lai(:)   ! m2 of leaf area / m2 col
     real(r8), allocatable :: col_sai(:)   ! m2 of stem area / m2 col
  end type scel_type


  ! This type contains the pre-processed scattering coefficients
  ! and routines.  This is the parent type that holds almost everything
  ! in the two-stream solver.
  type twostream_type
  
     type(scel_type), allocatable :: scel(:)     ! scattering elements for each layer
     integer                      :: n_scel      ! Number of scattering elements
     real(r8)                     :: grnd_albedo ! Albedo of the ground

     ! Element specific coefficients
     ! These contain weighted averages of element composition, stem, leaf and snow
     ! Ideally, we would like to spend time updating these

     real(r8), allocatable        :: eta_e(:)   ! mean element effective scattering area density
     real(r8), allocatable        :: kd_e(:)    ! mean element optical depth per unit scattering
                                                ! area of media
     real(r8), allocatable        :: om_e(:)    ! material scattering albedo (fraction)
     
     real(r8), allocatable        :: betad_e(:) ! Mean element backscatter fraction
                                                ! from diffuse radiation

     real(r8), allocatable        :: kb_e(:)    !
     real(r8), allocatable        :: betab_e(:) !

     real(r8)                     :: frac_snow     ! Current mean snow-fraction of the canopy
     real(r8)                     :: frac_snow_old ! Previous mean snow-fraction of the canopy
     
   contains

     procedure :: ZenithPrep     ! Update coefficients as zenith changes
     procedure :: CanopyPrep     ! Update coefficients as canopy changes
     procedure :: Solve          ! Perform the scattering solution
     procedure :: DecomposeFates ! Convert a FATES patch into scattering elements

     ! Note: Re-composing the result of scattering (such as absorbed PAR) back
     !       onto FATES cohorts happens outside this module, but will call
     !       this module's data structures.

     
  end type twostream_type


  ! Assumptions 1=visible (vis)
  !             2=near infrared (nir)
  !             3=thermal 
  

  ! Constants or reasonable approximations thereof
  real(r8), parameter :: snow_scatter_vis = 0.8_r8
  real(r8), parameter :: snow_scatter_nir = 0.4_r8
  real(r8), parameter :: snow_scatter_thm = 0.4_r8
  real(r8), parameter :: k_snow           = 1.0_r8  

  
contains


  subroutine ParamPrep(numpft,ib)

    real(r8) :: mu

    if(ib==1) then
       rad_params%om_snow(ib)    = snow_scatter_vis
       rad_params%betad_snow(ib) = 0.75_r8
    elseif(ib==2) then
       rad_params%om_snow(ib) = snow_scatter_nir
       rad_params%betad_snow(ib) = 0.75_r8
    else
       rad_params%om_snow(ib) = snow_scatter_thm
       rad_params%betad_snow(ib) = 0.75_r8
    end if
    
    do ft = 1,numpft

       ! The non-band specific parameters here will be re-derived for each
       ! band, which is inefficient, however this is an incredibly cheap
       ! routine to begin with, its only called during initialization, so
       ! just let it go, dont worry about it.
       
       ! There must be protections on xl to prevent div0 and other weirdness
       rad_params%phi1(ft) = 0.5_r8 - 0.633_r8*rad_params%xl(ft) - 0.330_r8*rad_params%xl(ft)*rad_params%xl(ft)
       rad_params%phi2(ft) = 0.877_r8 * (1._r8 - 2._r8*rad_params%phi1b(ft)) !0 = horiz leaves, 1 - vert leaves.

       mu = (1._r8/rad_params%phi2(ft))* &
            (1._r8-(rad_params%phi1(ft)/rad_params%phi2(ft))* &
            log((rad_params%phi2(ft)+rad_params%phi1(ft))/rad_params%phi1(ft)))
       
       rad_params%kd_leaf(ft) = 1/mu
       rad_params%kd_stem(ft) = 1._r8  ! Isotropic assumtion

       rad_params%om_leaf(ib,ft) = rad_params%rhol(ib,ft) + taul(ib,ft)
       rad_params%om_stem(ib,ft) = rad_params%rhos(ib,ft) + taus(ib,ft)
       
       rad_params%om_leaf(ib,ft) = rad_params%rhol(ib,ft) + taul(ib,ft)
       rad_params%om_stem(ib,ft) = rad_params%rhos(ib,ft) + taus(ib,ft)

    end do
    

  end subroutine ParamPrep

  
  subroutine CanopyPrep(this,ib,force)

    ! Pre-process things that change with canopy-geometry or snow cover
    ! We try to only run this when necessary. For instance we only
    ! run this when the canopy vegetation composition changes, or
    ! when the amount of snow-cover changes.  
    
    class(twostream_type) :: this
    integer               :: ib     ! The band of interest
    logical               :: force  ! We forcibly run this routine
                                    ! if the canopy composition changes
    ! But we check if the snow conditions
    ! change during the high frequency calls
    ! as well.
    
    integer :: ican  ! scattering element canopy layer index (top down)
    integer :: icol  ! scattering element column
    integer :: ilem  ! scattering element index
    real(r8) :: rho  ! element mean material reflectance
    real(r8) :: tau  ! element mean material transmittance

    if(.not.force) then
       if(abs(this%frac_snow-this%frac_snow_old)<nearzero) then
          this%frac_snow_old = this%frac_snow
          return
       end if
    end if
    
    this%frac_snow_old = this%frac_snow
    
    ilem = 0
    do ican = 1,ubound(this%scel,1)

       do icol = 1,this%scel(ican)%n_col

          ilem = ilem + 1
          ft = this%scel(ican)%col_pft(icol)
          
          this%eta_e(ilem) = this%scel(ican)%col_lai(icol) + this%scel(ican)%col_sai(icol)

          ! Mean element transmission coefficients w/o snow effects
          this%kd_e(ilem) =  (this%scel(ican)%col_lai(icol) * rad_params%kd_leaf(ft) + &
                            this%scel(ican)%col_sai(icol) * rad_params%kd_stem(ft))/this%eta_e(ilem)
          
          this%om_e(ilem) =  (this%scel(ican)%col_lai(icol)*rad_params%om_leaf(ib,ft) + &
                            this%scel(ican)%col_sai(icol)*rad_params%om_stem(ib,ft))/this%eta_e(ilem)

          ! Mean element transmission coefficients adding snow optical depth
          this%kd_e(ilem) = this%frac_snow*k_snow + this%kd_e(ilem)

          this%om_e(ilem) = this%frac_snow*rad_params%om_snow(ib) + (1._r8-this%frac_snow)*this%om_e(ilem)
          
          ! Diffuse backscatter, taken from G. Bonan's code

          rho = (this%scel(ican)%col_lai(icol) * rad_params%rho_leaf(ft) + &
                 this%scel(ican)%col_sai(icol) * rad_params%rho_stem(ft))/this%eta_e(ilem)
          tau = (this%scel(ican)%col_lai(icol) * rad_params%tau_leaf(ft) + &
                 this%scel(ican)%col_sai(icol) * rad_params%tau_stem(ft))/this%eta_e(ilem)
          
          this%betad_e(ilem)  = 0.5_r8 / this%om_e(ilem) * &
               ( this%om_e(ilem) + (rho-tau) * ((1._r8+rad_params%xl(ft))/2._r8)**2._r8 )

          this%betad_e(ilem) = this%frac_snow*rad_params%betad_snow(ib) + (1._r8-this%frac_snow)*this%betad_e(ilem)
          
       end do
    end do
    
    return
  end subroutine CanopyPrep

  ! =============================================================================
  
  subroutine ZenithPrep(this,cosz,ib)

    ! Pre-process things that change with the zenith angle
    ! i.e. the beam optical properties

    class(twostream_type) :: this
    integer :: ican  ! scattering element canopy layer index (top down)
    integer :: icol  ! scattering element column
    integer :: ilem  ! scattering element index

        
    ilem = 0
    do ican = 1,ubound(this%scel,1)

       do icol = 1,this%scel(ican)%n_col

          ilem = ilem + 1
          ft = this%scel(ican)%col_pft(icol)

          gdir = rad_params%phi1(ft) + rad_params%phi2(ft) * cosz
          
          !how much direct light penetrates a singleunit of lai?
          this%kb_e(ilem) = gdir / cosz

          ! RGK: The snow is adding area, so I don't see this as a weighted
          !      average, but purely adding optical depth
          this%kb_e(ilem) = this%frac_snow*k_snow + this%kb_e(ilem)

          ! betab - upscatter parameter for direct beam radiation, from G. Bonan

          avmu = (1._r8 - rad_params%phi1(ft)/rad_params%phi2(ft) * &
               log((rad_params%phi1(ft)+rad_params%phi2(ft))/rad_params%phi1(ft))) / rad_params%phi2(ft)

          tmp0 = gdir + phi2 * cosz
          tmp1 = phi1 * cosz
          tmp2 = 1._r8 - tmp1/tmp0 * log((tmp1+tmp0)/tmp1)
          asu = 0.5_r8 * this%om_e(ilem) * gdir / tmp0 * tmp2
          this%betab_e(ilem) = (1._r8 + avmu*this%kb_e(ilem)) / (this%om_e(ilem)*avmu*this%kb_e(ilem)) * asu
       
    end do
 end do
 
 return
end subroutine ZenithPrep


  function GetNSCel(this) result(n_scel)

    ! Simply return the total number
    ! of scattering elements from the
    ! multi-layer scattering element array
    
    class(twostream_type) :: this
    integer :: n_scel
    
    n_scel = 0
    do ic = 1,ubound(this%scel,1)
       n_scel = n_scel + this%scel(ic)%n_col
    end do
    
  end function GetNSCel
  
  
  subroutine Solve(this,ib,Rbeam0)

    class(twostream_type) :: this
    integer,intent(in)    :: ib       ! band index
    real(r8)              :: Rbeam0
    real(r8)              :: Rdiff0

    integer :: n_cl    ! Number of canopy layers
    integer :: n_scel  ! Number of scattering elements

    ! Two stream solution arrays
    ! Each of these are given generic names, because
    ! they are assemblages of many terms. But generally
    ! they fit the linear algebra formulation:
    !
    ! tau(:) = omega(:,:) * lambda(:)
    !
    ! Where, we invert to solve for the coefficients lambda

    
    real(r8),allocatable :: omega(:,:)
    real(r8),allocatable :: tau(:)
    real(r8),allocatable :: lambda(:)

    integer :: ic  ! Loop index for canopy layers
    integer :: ilem  ! Loop index for scattering elements

    allocate(omega(2*this%n_scel,2*this%n_scel))
    allocate(tau(2*this%n_scel,1))
    allocate(lambda(2*this%n_scel,1))

    !baf
    !om
    !Kb
    !etad
    !bd
    !bb
    !b1

    deallocate(omega)
    deallocate(tau)
    deallocate(lambda)
    
    return
  end subroutine Solve

  ! ============================================================
  

  subroutine DecomposeFates(this)

    ! In this routine, we decompose a FATES canopy
    ! into the data structures upon which we
    ! can perform the two-stream calculations
    
    class(twostream_type) :: this

    

    return
  end subroutine DecomposeFates


end Module TwoStreamPPAMod
