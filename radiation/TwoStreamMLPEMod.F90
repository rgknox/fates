Module TwoStreamMLPEMod
  
  ! This module holds the routines to calculate two-tream
  ! radiation scattering of vegetation in "M"ultiple "L"ayers
  ! with "P"arellel "E"lements "MLPE"
  !
  ! In summary,
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

  implicit none
  private
  
  integer, parameter :: r8 = selected_real_kind(12)
  real(r8),parameter :: nearzero = 1.e-20_r8
  logical, parameter :: debug=.true.
  logical, parameter :: use_derivation1 = .true.
  real(r8), parameter :: unset_r8 = 1.e-36_r8

  
  integer, parameter :: twostr_vis = 1         ! Named index of visible shortwave radiation
  integer, parameter :: twostr_nir = 2         ! Named index for near infrared shortwave radiation
  
  ! These are parameter constants, ie things that are specific to the plant material
  ! and radiation band.  Not all of these need to be used. 2-stream ultimately wants
  ! optical depth, material reflectance and backscatter fractions for diffuse and
  ! direct light. So there are various ways to get to these parameters, depending
  ! on the host model's available parameters.  The rho,tau,xl and clumping parameters
  ! are standard elm/clm parameters, and provided as a convenience.
  
  type, public :: rad_params_type

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
     real(r8), allocatable :: om_leaf(:,:)   ! Leaf material scattering albedo (band x pft)
     real(r8), allocatable :: om_stem(:,:)   ! Stem material scattering albedo (band x pft)

     real(r8), allocatable :: om_snow(:)     ! Snow material scattering albedo for snow (band)
     real(r8), allocatable :: betad_snow(:)  ! Diffuse backscatter fraction for snow (band)
     
  end type rad_params_type
  
  type(rad_params_type),public :: rad_params

  
  ! Information describing the scattering elements from a host model

  type scel_type
     integer  :: pft      ! pft index
     real(r8) :: area     ! m2 col/m2 ground
     real(r8) :: lai      ! m2 of leaf area / m2 col
     real(r8) :: sai      ! m2 of stem area / m2 col
  end type scel_type


  ! Information describing the scattering coefficients
  ! (this is allocated for each broad band)
  type scco_type
  
     ! Terms used in the final solution, also used for decomposing solution
     real(r8) :: Au       ! Compound intercept term
     real(r8) :: Ad       ! Compound intercept term
     real(r8) :: B1u       ! Compound term w/ lambdas
     real(r8) :: B2u       ! Compound term w/ lambdas
     real(r8) :: B1d       ! Compound term w/ lambdas
     real(r8) :: B2d       ! Compound term w/ lambdas
     real(r8) :: lambda1  ! Compount term w/ B1d and B1u
     real(r8) :: lambda2  ! Compound term w/ B2d and B2u
     real(r8) :: a        ! Complex term operating on veg area index
     real(r8) :: Kb       ! Optical depth of beam radiation
     real(r8) :: Kb_leaf  ! Optical depth of just leaves in beam radiation
     real(r8) :: Kd       ! Optical depth of diffuse radiation
     real(r8) :: om       ! material reflectance. ie portion that
                          ! is reflected from material assuming impact (fraction)
     real(r8) :: betad    ! backscatter fraction of diffuse radiation
     real(r8) :: betab    ! backscatter fraction of beam radiation
     
     real(r8) :: Rbeam0   ! Downwelling beam radiation at
                          ! top of the element [w/m2]

   contains 
       
     procedure :: GetRdUp
     procedure :: GetRdDn
     procedure :: GetRb
  end type scco_type


  type band_type

     type(scco_type), pointer :: scco(:,:)            ! array of scattering coefficients (layer, column)
                                                      ! can be sparse, will only solve indices up to

     integer                  :: ib                   ! band index, should be consistent with rad_params
     real(r8)                 :: albedo_grnd_diff     ! Ground albedo diffuse
     real(r8)                 :: albedo_grnd_beam     ! Ground albedo direct 

  end type band_type
  

  ! This type contains the pre-processed scattering coefficients
  ! and routines.  This is the parent type that holds almost everything
  ! in the two-stream solver.
  ! The scel structure describes the scattering elements, these are values
  ! that need to be defined by the ecosystem model, somewhat of
  ! an input to the solver. Since this is a Perfect Plasticity Approximation
  ! enabled system, we partition the scattering media into "columns" and "layers"
  ! Layers are canopy layers, think understory, mid-story and upper canopy. Columns
  ! are divisions of horizontal space, ie literal columns of space. The current
  ! implementation limits this space to media that has uniform scattering coefficients.
  ! So there could not be different PFTs in the same column, because they would undoubtedly
  ! have different joint scattering coefficients at different height levels in
  ! the column.  Therefore, every column is connected with a PFT.

  
  type, public :: twostream_type
  
     type(scel_type), pointer :: scel(:,:)            ! array of scattering elements (layer, column)
                                                      ! can be sparse, will only solve indices up to
                                                      ! n_lyr,n_col(n_lyr)

     type(band_type), pointer :: band(:)              ! Holds scattering coefficients for each band
                                                      ! vis,nir,etc (nothing that emits though, no thermal)

     integer                      :: n_lyr            ! number of (vertical) scattering element layers
     integer, allocatable         :: n_col(:)         ! number of (horizontal) scattering element columns per layer
     integer                      :: n_scel           ! total number of scattering elements
     logical                      :: force_prep       ! Some coefficients are only updated
                                                      ! when the canopy composition changes, ie
                                                      ! changes in leaf, stem or snow structure.
                                                      ! If so, this sets to true, signalling that diffuse
                                                      ! scattering coefficients should be updated.
                                                      ! Otherwise, we only updated zenith dependent
                                                      ! parameters on short sub-daily timesteps
     real(r8)                     :: frac_snow        ! Current mean snow-fraction of the canopy
     real(r8)                     :: frac_snow_old    ! Previous mean snow-fraction of the canopy

   contains
     
     procedure :: ZenithPrep     ! Update coefficients as zenith changes
     procedure :: CanopyPrep     ! Update coefficients as canopy changes
     procedure :: Solve          ! Perform the scattering solution
     procedure :: GetNSCel
     procedure :: AllocInitTwoStream
     procedure :: DeallocTwoStream

  end type twostream_type



  ! Assumptions 1=visible (vis)
  !             2=near infrared (nir)
  !             3=thermal 
  

  ! Constants or reasonable approximations thereof
  real(r8), parameter :: snow_scatter_vis = 0.85_r8  ! Tarboton 1995
  real(r8), parameter :: snow_scatter_nir = 0.75_r8 
  real(r8), parameter :: snow_scatter_thm = 0.65_r8  ! Tarboton 1995
  real(r8), parameter :: k_snow           = 1.0_r8  

  ! For air, use a nominal values to prevent div0s
  ! the key is that tai = 0
  real(r8), parameter :: k_air = 0.5_r8  
  real(r8), parameter :: om_air  = 0.5_r8
  real(r8), parameter :: beta_air = 0.5_r8

  integer, parameter :: air_ft = 0  ! 

  
  public :: ParamPrep
  public :: AllocateRadParams
  public :: GetAbsRad
  
contains

  subroutine AllocInitTwoStream(this,band_indices,ncan,ncol)

    class(twostream_type) :: this
    integer               :: band_indices(:)
    integer               :: ncan
    integer               :: ncol
    
    integer :: nbands
    integer :: ib

    nbands = ubound(band_indices,1)

    allocate(this%n_col(ncan)) 
    allocate(this%scel(ncan,ncol))
    allocate(this%band(nbands))
    
    this%frac_snow        = unset_r8
    this%frac_snow_old    = unset_r8

    do ib = 1,nbands
       
       allocate(this%band(ib)%scco(ncan,ncol))
       this%band(ib)%albedo_grnd_diff = unset_r8
       this%band(ib)%albedo_grnd_beam = unset_r8
       this%band(ib)%ib = band_indices(ib)

    end do
       
    return
  end subroutine AllocInitTwoStream

  ! ===============================================================================================
  
  subroutine DeallocTwoStream(this)

    class(twostream_type) :: this

    integer               :: nbands
    integer               :: ib
    
    nbands = ubound(this%band,1)
    
    deallocate(this%scel)
    deallocate(this%n_col)
    do ib = 1,nbands
       deallocate(this%band(ib)%scco)
    end do
    deallocate(this%band)
      
    return
  end subroutine DeallocTwoStream
  
  ! ===============================================================================================
  
  subroutine AllocateRadParams(n_pft,n_bands)

    integer,intent(in) :: n_pft
    integer,intent(in) :: n_bands
    
    allocate(rad_params%rhol(n_bands,n_pft))
    allocate(rad_params%rhos(n_bands,n_pft))
    allocate(rad_params%taul(n_bands,n_pft))
    allocate(rad_params%taus(n_bands,n_pft))
    allocate(rad_params%xl(n_pft))
    allocate(rad_params%clumping_index(n_pft))

    allocate(rad_params%phi1(n_pft))
    allocate(rad_params%phi2(n_pft))
    allocate(rad_params%kd_leaf(n_pft))
    allocate(rad_params%kd_stem(n_pft))
    allocate(rad_params%om_leaf(n_bands,n_pft))
    allocate(rad_params%om_stem(n_bands,n_pft))
    allocate(rad_params%om_snow(n_bands))
    allocate(rad_params%betad_snow(n_bands))
    
  end subroutine AllocateRadParams
  
  ! ================================================================================================

  function GetRdDn(this,vai) result(r_diff_dn)
    class(scco_type)    :: this
    real(r8),intent(in) :: vai
    real(r8)            :: r_diff_dn
    
    ! Rdn = Ad e−(Kbv) + Re + λ1 B1d e^(av) + λ2 B2d e^(−av)
    r_diff_dn = this%Ad*exp(-this%Kb*vai) + this%B1d*this%lambda1*exp(this%a*vai) + this%B2d*this%lambda2*exp(-this%a*vai)
  end function GetRdDn

  function GetRdUp(this,vai) result(r_diff_up)
    class(scco_type)    :: this
    real(r8),intent(in) :: vai
    real(r8)            :: r_diff_up
    
    ! Rup = Au e−(Kbv) + Re + λ1 B1u e^(av) + λ2 B2u e^(−av)
    r_diff_up = this%Au*exp(-this%Kb*vai) + this%B1u*this%lambda1*exp(this%a*vai) + this%B2u*this%lambda2*exp(-this%a*vai)
  end function GetRdUp
  
  function GetRb(this,vai) result(r_beam_dn)
    class(scco_type) :: this
    real(r8),intent(in)   :: vai
    real(r8)              :: r_beam_dn

    r_beam_dn = this%Rbeam0*exp(-this%Kb*vai)
  end function GetRb
  
  
  subroutine GetAbsRad(scel,scco,ib,vai_top,vai_bot,Rd_abs_leaf,Rb_abs_leaf,R_abs_stem)

    ! This routine is used to help decompose radiation scattering
    ! and return the amount of absorbed radiation.  The canopy layer and column
    ! index identify the element of interest. The other arguments are the upper and
    ! lower bounds within the element over which to evaluate absorbed radiation.
    ! The assumption is that the vegetation area index is zero at the top of the
    ! element, and increases going downwards.  As with all assumptions in this
    ! module, the scattering parameters are uniform within the element itself,
    ! which includes an assumption of the leaf/stem proportionality.
    ! ---------------------------------------------------------------------------
    ! Solution for radiative intensity of diffuse up and down at tai=v
    ! Rup = Au e−(Kbv) + Re + λ1 B1u e^(av) + λ2 B2u e^(−av)
    ! Rdn = Ad e−(Kbv) + Re + λ1 B1d e^(av) + λ2 B2d e^(−av)
    ! ---------------------------------------------------------------------------
    
    ! Arguments
    type(scel_type)      :: scel
    type(scco_type)      :: scco    
    integer,intent(in)    :: ib          ! broad band index
    real(r8),intent(in)   :: vai_top     ! veg area index (from the top of element) to start
    real(r8),intent(in)   :: vai_bot     ! veg area index (from the top of element) to finish
    real(r8), intent(out) :: Rb_abs_leaf ! Absorbed beam radiation from leaves [W/m2 ground]
    real(r8), intent(out) :: Rd_abs_leaf ! Absorbed diff radiation from leaves [W/m2 ground]
    real(r8), intent(out) :: R_abs_stem  ! Absorbed beam+diff radiation stems  [W/m2 ground]

    real(r8) :: Rd_abs                   ! Absorbed diffuse radiation [W/m2 ground]
    real(r8) :: Rb_abs                   ! Absorbed beam radiation 
    real(r8) :: vai_max                  ! total integrated (leaf+stem) area index of the current element
    real(r8) :: diff_wt_leaf             ! diffuse absorption weighting for leaves
    real(r8) :: diff_wt_stem             ! diffuse absorption weighting for stems
    real(r8) :: beam_wt_leaf             ! beam absorption weighting for leaves
    real(r8) :: beam_wt_stem             ! beam absorption weighting for stems
    
    associate( Rbeam0  => scco%Rbeam0, &
               a       => scco%a, &
               om      => scco%om, &
               lai     => scel%lai, &
               sai     => scel%sai, &
               ft      => scel%pft)

      ! The total vegetation area index of the element
      vai_max = lai + sai

      ! We have to disentangle the absorption between leaves and stems, we give them both
      ! a weighting fraction of total absorption of  area*K*(1-om)

      diff_wt_leaf = lai*(1._r8-rad_params%om_leaf(ft,ib))*rad_params%kd_leaf(ft)
      diff_wt_stem = sai*(1._r8-rad_params%om_stem(ft,ib))*rad_params%kd_stem(ft)

      beam_wt_leaf = lai*(1._r8-rad_params%om_leaf(ft,ib))*scco%Kb_leaf
      beam_wt_stem = sai*(1._r8-rad_params%om_stem(ft,ib))*1._r8

      if(debug) then
         if(vai_bot> vai_max)then
            print*,"During decomposition of the 2-stream radiation solution"
            print*,"A vegetation area index (VAI) was requested in GetAbsRad()"
            print*,"that is larger than the total integrated VAI of the "
            print*,"computation element of interest."
            print*,"vai_max: ",vai_max
            print*,"vai_bot: ",vai_bot
            stop
         end if
         if(vai_bot<vai_top) then
            print*,"During decomposition of the 2-stream radiation solution"
            print*,"the vegetation area index at the lower position was set"
            print*,"as greater than the value at the upper position."
            print*,"vai_max: ",vai_max
            print*,"vai_bot: ",vai_bot
            stop
         end if
      end if

      
      ! Amount of absorbed radiation is retrieved by doing an energy
      ! balance on this boundaries over the depth of interest
      ! Result is Watts / m2 of the element's area footprint NOT
      ! per m2 of tissue (at least not in this step)
      
      Rb_abs = scco%GetRb(vai_top)-scco%GetRb(vai_bot)*(1._r8-om)  
      Rd_abs = (scco%GetRdDn(vai_top) - scco%GetRdDn(vai_bot)) + &
               (scco%GetRdUp(vai_bot) - scco%GetRdUp(vai_top))
      
      ! Fraction of leaves exposed to direct sunlight
      Rb_abs_leaf = Rb_abs * beam_wt_leaf / (beam_wt_leaf+beam_wt_stem)
      Rd_abs_leaf = Rd_abs * diff_wt_leaf / (diff_wt_leaf+diff_wt_stem)
      
      R_abs_stem = Rb_abs*beam_wt_stem / (beam_wt_leaf+beam_wt_stem) + &
                   Rd_abs*diff_wt_stem / (diff_wt_leaf+diff_wt_stem)

      
    end associate
    return
  end subroutine GetAbsRad

  ! ================================================================================================
  
  subroutine ParamPrep(numpft,ib)

    integer,intent(in) :: numpft   ! Number of pfts
    integer,intent(in) :: ib       ! Band index
    
    real(r8) :: mu
    integer  :: ft

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
       rad_params%phi2(ft) = 0.877_r8 * (1._r8 - 2._r8*rad_params%phi1(ft)) !0 = horiz leaves, 1 - vert leaves.

       mu = (1._r8/rad_params%phi2(ft))* &
            (1._r8-(rad_params%phi1(ft)/rad_params%phi2(ft))* &
            log((rad_params%phi2(ft)+rad_params%phi1(ft))/rad_params%phi1(ft)))
       
       rad_params%kd_leaf(ft) = 1/mu
       rad_params%kd_stem(ft) = 1._r8  ! Isotropic assumtion

       rad_params%om_leaf(ib,ft) = rad_params%rhol(ib,ft) + rad_params%taul(ib,ft)
       rad_params%om_stem(ib,ft) = rad_params%rhos(ib,ft) + rad_params%taus(ib,ft)
       
    end do
    
    return
  end subroutine ParamPrep

  ! ================================================================================================
  
  subroutine CanopyPrep(this,ib)

    ! Pre-process things that change with canopy-geometry or snow cover
    ! We try to only run this when necessary. For instance we only
    ! run this when the canopy vegetation composition changes, or
    ! when the amount of snow-cover changes.  
    
    class(twostream_type) :: this
    integer               :: ib     ! The band of interest

    ! But we check if the snow conditions
    ! change during the high frequency calls
    ! as well.
    
    integer :: ican  ! scattering element canopy layer index (top down)
    integer :: icol  ! scattering element column
    real(r8) :: rho  ! element mean material reflectance
    real(r8) :: tau  ! element mean material transmittance
    real(r8) :: vai  ! vegetation area index lai+sai
    !type(scel_type),pointer :: scelp   ! Pointer to the scel data structure
    !type(scco_type),pointer :: sccop   ! Pointer to the scco data structure
    
    if(.not.this%force_prep) then
       if(abs(this%frac_snow-this%frac_snow_old)<nearzero) then
          this%frac_snow_old = this%frac_snow
          return
       end if
    end if
    
    this%frac_snow_old = this%frac_snow
    
    do ican = 1,this%n_lyr
       do icol = 1,this%n_col(ican)

          associate(lai => this%scel(ican,icol)%lai, &
                    sai => this%scel(ican,icol)%sai, &
                    ft  => this%scel(ican,icol)%pft, &
                    sclep => this%scel(ican,icol),   &
                    sccop => this%band(ib)%scco(ican,icol))
            
            if(ft==0)then
               ! Simple provisions for a ghost element (air)
               sccop%kd = k_air
               sccop%om = om_air
               sccop%betad = beta_air
            else

               vai = lai + sai

               ! Mean element transmission coefficients w/o snow effects
               sccop%kd =  (lai * rad_params%kd_leaf(ft) + &
                    sai * rad_params%kd_stem(ft))/vai

               sccop%om =  (lai*rad_params%om_leaf(ib,ft) + &
                    sai*rad_params%om_stem(ib,ft))/vai

               ! Mean element transmission coefficients adding snow optical depth
               !!sccop%Kd = this%frac_snow*k_snow + sccop%Kd

               !!sccop%om = this%frac_snow*rad_params%om_snow(ib) + (1._r8-this%frac_snow)*sccop%om

               ! Diffuse backscatter, taken from G. Bonan's code

               rho = (lai * rad_params%rhol(ib,ft) + &
                    sai * rad_params%rhos(ib,ft))/vai
               tau = (lai * rad_params%taul(ib,ft) + &
                    sai * rad_params%taus(ib,ft))/vai

               sccop%betad  = 0.5_r8 / sccop%om * &
                    ( sccop%om + (rho-tau) * ((1._r8+rad_params%xl(ft))/2._r8)**2._r8 )

               !!sccop%betad  = this%frac_snow*rad_params%betad_snow(ib) + (1._r8-this%frac_snow)*sccop%betad 
            end if
          end associate
       end do
    end do
    
    return
  end subroutine CanopyPrep

  ! ================================================================================================
  
  subroutine ZenithPrep(this,ib,cosz)

    ! Pre-process things that change with the zenith angle
    ! i.e. the beam optical properties

    class(twostream_type) :: this
    integer               :: ib      ! band index, matches indexing of rad_params
    real(r8),intent(in)   :: cosz    ! Cosine of the zenith angle

    integer :: ican  ! scattering element canopy layer index (top down)
    integer :: icol  ! scattering element column
    real(r8) :: asu  ! single scattering albedo
    real(r8) :: avmu ! Average inverse diffuse optical depth per unit leaf area
    real(r8) :: gdir
    real(r8) :: tmp0,tmp1,tmp2
    type(scco_type), pointer :: sccop
    
    do ican = 1,this%n_lyr
       do icol = 1,this%n_col(ican)

          sccop =>  this%band(ib)%scco(ican,icol)
          
          associate(ft => this%scel(ican,icol)%pft, &
               scelp => this%scel(ican,icol))

            if(ft==0)then

               ! Simple provisions for a ghost element (air)
               sccop%Kb_leaf = k_air
               sccop%Kb = k_air
               sccop%betab = beta_air

            else
               gdir = rad_params%phi1(ft) + rad_params%phi2(ft) * cosz

               !how much direct light penetrates a singleunit of lai?
               sccop%Kb_leaf = gdir / cosz

               sccop%Kb = (scelp%lai*sccop%Kb_leaf + scelp%sai*1.0)/(scelp%lai+scelp%sai)

               ! RGK: The snow is adding area, so I don't see this as a weighted
               !      average, but purely adding optical depth
               !!sccop%kb = this%frac_snow*k_snow + sccop%kb

               ! betab - upscatter parameter for direct beam radiation, from G. Bonan

               avmu = (1._r8 - rad_params%phi1(ft)/rad_params%phi2(ft) * &
                    log((rad_params%phi1(ft)+rad_params%phi2(ft))/rad_params%phi1(ft))) / rad_params%phi2(ft)

               tmp0 = gdir +  rad_params%phi2(ft) * cosz
               tmp1 =  rad_params%phi1(ft) * cosz
               tmp2 = 1._r8 - tmp1/tmp0 * log((tmp1+tmp0)/tmp1)
               asu = 0.5_r8 * sccop%om * gdir / tmp0 * tmp2
               sccop%betab = (1._r8 + avmu*sccop%kb) / (sccop%om*avmu*sccop%kb) * asu

            end if

            !!sccop%betab = this%frac_snow*beta_snow + sccop%betab

            !this%albedo_grnd_beam = 1.e-36  ! Must fill this in
            
          end associate
       end do
    end do

    return
  end subroutine ZenithPrep
  
  ! ================================================================================================

  subroutine GetNSCel(this)

    ! Simply return the total number
    ! of scattering elements from the
    ! multi-layer scattering element array
    
    class(twostream_type) :: this
    integer :: ican
    
    this%n_scel = 0
    do ican = 1,this%n_lyr
       this%n_scel = this%n_scel + this%n_col(ican)
    end do
    return
  end subroutine GetNSCel

  ! ===============================================================
  
  subroutine Solve(this,ib,Rbeam_atm,Rdiff_atm)

    class(twostream_type) :: this
    integer               :: ib         ! Band of interest, matches indexing of rad_params
    real(r8)              :: Rbeam_atm  ! Intensity of beam radiation at top of canopy [W/m2]
    real(r8)              :: Rdiff_atm  ! Intensity of diffuse radiation at top of canopy [W/m2]

    ! Two stream solution arrays
    ! Each of these are given generic names, because
    ! they are assemblages of many terms. But generally
    ! they fit the linear algebra formulation:
    !
    ! TAU(:) = OMEGA(:,:) * LAMBDA(:)
    !
    ! Where, we invert to solve for the coefficients LAMBDA
    
    real(r8),allocatable :: OMEGA(:,:)
    real(r8),allocatable :: TAU(:)
    real(r8),allocatable :: LAMBDA(:)

    integer :: ican  ! Loop index for canopy layers
    integer :: ibot  ! layer index for top side of layer divide
    integer :: itop  ! layer index for bottom side of layer divide
    integer :: icol  ! Loop index for canopy columns
    integer :: jcol  ! Another loop index for canopy columns
    integer :: ilem  ! Index for scattering elements
    integer :: k1,k2 ! Indices for the lambda terms in the OMEGA and LAMBDA array
    integer :: qp    ! Equation position index
    
    integer :: ilem_off ! Offset, or total number of elements above layer of interest
    real(r8) :: b1,b2,a2,nu_sqrd,nu ! intermediate terms, see documentation
    real(r8) :: Rbeam_top           ! Mean beam radiation at top of layer      [W/m2]
    real(r8) :: Rbeam_bot           ! Mean beam radiation at bottom of layer   [W/m2]
    real(r8) :: Rbeam_exp_ground    ! Mean beam radiation incident on exposed ground   [W/m2]
    real(r8) :: frac_exp_ground     ! Fraction of ground that is not occupied by columns [m2/m2]
    real(r8) :: vai                 ! Vegetation area index [m2 vegetation / m2 ground]

    type(scel_type),pointer :: scelp   ! Pointer to the scel data structure
    type(scco_type),pointer :: sccop   ! Pointer to the scco data structure
    
    ! Parameters for solving via LAPACK DGELS()
    character(1),parameter :: trans = 'N'           ! Input matrix is not transposed
    integer, parameter :: workmax = 100             ! Maximum iterations to minimize work
    real(r8) :: work(workmax)                       ! Work array
    integer  :: lwork                               ! Dimension of work array
    integer :: info                                 ! Procedure diagnostic ouput

    ! Testing switch
    ! If true, then allow elements
    ! of different layers, but same row, to have priority
    ! flux into the other element, instead of a mix
    logical, parameter :: continuity_on = .true.    
    
    ! ------------------------------------------------------------------------------------
    ! Example system of equations for 2 parallel columns in each of two canopy
    ! layers.  Each line is one of the balanc equations. And the x's are
    ! the unknown coefficients used in those equations.  2 coefficients
    ! map to each element, and read left to right.
    ! EL1 is the element in top layer left column.
    ! EL2 is the element in the top layer, right column
    ! EL3 is the element in the bottom layer, left column
    ! EL4 is the element in the bottom layer, right column
    !
    !                                              EL1 EL2 EL3 EL4
    ! EQ: Idn balance with upper BC can1, col 1:   x x 
    ! EQ: Idn balance with upper BC can1, col 2:       x x
    ! EQ: Idn balance between upper & lower        x x x x x x
    ! EQ: Idn balance between upper & lower        x x x x     x x
    ! EQ: Iup balance between lower & upper        x x x x x x x x
    ! EQ: Iup balance between lower & upper        x x x x x x x x
    ! EQ: Iup/Idn balance with ground, 1st col:            x x
    ! EQ: Iup/Idn Balance with ground, 2nd lower col:          x x  
    !
    !     Note: The Iup balance between layers requires ALL
    !     terms, because light comes out of both
    !     upper canopy elements and reflects off soil
    !     AND, light upwells from both lower elements.
    !
    ! --------------------------------------------------------------------------
    
    ! --------------------------------------------------------------------------
    ! Beam Scattering
    ! First do the direct beam stuff.  It is a trivial solution
    ! and is required as a boundary condition to the diffuse solver
    ! All parallel layers recieve downwelling form the
    ! atmosphere.
    ! Rbeam0 is the upper boundary condition provided by data or another
    ! model.
    ! Rbeam() is the incident beam radiation at the top of each layer
    ! upper canopy.
    ! --------------------------------------------------------------------------

    Rbeam_top = Rbeam_atm

    do ican = 1,this%n_lyr

       Rbeam_bot = 0._r8
       do icol = 1,this%n_col(ican)
          scelp => this%scel(ican,icol)
          sccop => this%band(ib)%scco(ican,icol)
          sccop%Rbeam0 = Rbeam_top
          Rbeam_bot = Rbeam_bot + &
               Rbeam_top*scelp%area*exp(-sccop%Kb*(scelp%lai+scelp%sai))
       end do

       ! Save the beam intensity that would be imposed on any exposed earth
       if(ican==this%n_lyr) then
          frac_exp_ground=1._r8
          do icol = 1,this%n_col(ican)
             frac_exp_ground = frac_exp_ground - this%scel(ican,icol)%area
          end do
          
          Rbeam_exp_ground = Rbeam_top
          
       end if
       Rbeam_top = Rbeam_bot
    end do

    ! Calculate element-level intermediate terms to the solve
    ! These are dependent on leaf level scattering and beam scattering
    ! These values will be used to populate the matrix solve
    ! =====================================================================

    do ican = 1,this%n_lyr
       do icol = 1,this%n_col(ican)

          scelp => this%scel(ican,icol)
          sccop => this%band(ib)%scco(ican,icol)
          
          a2 = sccop%Kd*sccop%Kd*(sccop%om-1._r8)*(sccop%om-1._r8-2._r8*sccop%om*sccop%betad)
          
          if(a2<0._r8) then
             print*,'a^2 is less than zero'
             stop
          end if
          
          sccop%a  = sqrt(a2)
       
          b1 = (sccop%Kd*(1._r8-sccop%om)*(1._r8-2._r8*sccop%betab)+sccop%Kb) * &
               sccop%om*sccop%Kb*sccop%Rbeam0
          b2 = (sccop%Kd*(sccop%om-1._r8-2._r8*sccop%om*sccop%betad) - &
               (1._r8-2._r8*sccop%betab)*sccop%Kb) * &
               sccop%om*sccop%Kb*sccop%Rbeam0

          if(use_derivation1) then
             
             nu_sqrd = (1._r8-sccop%om+2._r8*sccop%om*sccop%betad)/(1._r8-sccop%om)
             
             if(nu_sqrd<0._r8)then
                print*,'nu_sqrd is less than zero'
                stop
             end if
             
             ! B_1 up term from documentation:
             sccop%B1u  = 0.5_r8*(1._r8+sqrt(nu_sqrd))
             
             ! B_2 up term from documentation
             sccop%B2u = 0.5_r8*(1._r8-sqrt(nu_sqrd))

             ! B_1 down term from documentation:
             sccop%B1d  = -0.5_r8*(1._r8-sqrt(nu_sqrd))
             
             ! B_2 down term from documentation
             sccop%B2d = -0.5_r8*(1._r8+sqrt(nu_sqrd))
             
             ! A_2 term from documentation
             sccop%Ad    = -0.5_r8*(b2-b1)/(sccop%a*sccop%a-sccop%Kb*sccop%Kb)   ! aka half b2 minus b1
             
             ! A_1 term from documentation
             sccop%Au    = -0.5_r8*(b2+b1)/(sccop%a*sccop%a-sccop%Kb*sccop%Kb)   ! aka half b1 plus b2

          else

             nu_sqrd = (sccop%om-1._r8)/(sccop%om - 1._r8-2._r8*sccop%om*sccop%betad)

             nu = (sccop%Kd*(sccop%om-1._r8))/sccop%a
             
             b1 = -b1
             
             ! B 1 up term from documentation
             !sccop%B1u  = 0.5_r8*(1._r8-nu)
             sccop%B1u  = 0.5_r8*(1._r8-sqrt(nu_sqrd))
             
             ! B_2 term from documentation
             !sccop%B2u = 0.5_r8*(1._r8+nu)
             sccop%B2u  = 0.5_r8*(1._r8+sqrt(nu_sqrd))
             
             ! B 1 up term from documentation
             !sccop%B1d  = 0.5_r8*(1._r8+nu)
             sccop%B1d  = 0.5_r8*(1._r8+sqrt(nu_sqrd))
             
             ! B_2 term from documentation
             !sccop%B2d = 0.5_r8*(1._r8-nu)
             sccop%B2d  = 0.5_r8*(1._r8-sqrt(nu_sqrd))
             
             ! A_2 term from documentation
             sccop%Ad    = -0.5_r8*(b2+b1)/(sccop%a*sccop%a-sccop%Kb*sccop%Kb)   ! aka half b2 minus b1
             
             ! A_1 term from documentation
             sccop%Au    = -0.5_r8*(b2-b1)/(sccop%a*sccop%a-sccop%Kb*sccop%Kb)   ! aka half b1 plus b2
             
          end if
          
          
       end do
    end do
    
    ! =====================================================================
    ! Set up the linear systems solver
    !
    ! [TAU] = [OMEGA]*[LAMBDA]
    ! OMEGA(n_equations,n_coefficients)
    ! TAU(n_equations)
    ! LAMBDA (n_coefficients) (the solution)
    !
    ! Indexing Variables
    ! ilem : element position
    ! k1 : coefficient 1 position
    ! k2 : coefficient 2 position
    ! qp : equation position, this continues to increment
    ! =====================================================================

    ! TO-DO: MAKE THIS SCRATCH SPACE AT THE SITE LEVEL?
    allocate(OMEGA(2*this%n_scel,2*this%n_scel))
    allocate(TAU(2*this%n_scel))
    allocate(LAMBDA(2*this%n_scel))

    OMEGA  = 0._r8
    !TAU    = 0._r8 ! Dont need to zero thsese
    !LAMBDA = 0._r8
    
    ! --------------------------------------------------------------------
    ! I. Flux equations with the atmospheric boundary
    ! These balance with all elements in the upper
    ! canopy, only.  The upper canopy is layer 1.
    ! --------------------------------------------------------------------
    
    qp   = 0    ! e"Q"uation "P"osition

    do icol = 1,this%n_col(1)

       scelp => this%scel(1,icol)
       sccop => this%band(ib)%scco(1,icol)
       ilem = icol
       qp   = qp   + 1
       
       k1 = 2*(ilem-1)+1
       k2 = k1+1
       
       TAU(qp)      =  Rdiff_atm - sccop%Ad
       OMEGA(qp,k1) =  sccop%B1d
       OMEGA(qp,k2) =  sccop%B2d

    end do
    
    
    if_understory: if(this%n_lyr>1) then


       ! -------------------------------------------------------------------
       ! II. Flux equations between canopy layers, DOWNWELLING
       ! We only perform flux balancing between layers
       ! if we have any understory, this is true if ican>1
       ! -------------------------------------------------------------------
       ! Refer to Equation X in technical document
       ! ------------------------------------------------------------
       
       ! This is the index offset for the layer above the
       ! current layer of interest. We start by evaluating
       ! Layer 2, so the offset refers to layer 1, and a
       ! value of 0
       ilem_off = 0
       
       do_dn_ican: do ican = 2,this%n_lyr

          itop = ican-1  ! Top layer of the balance
          ibot = ican    ! Bottom layer of the balance
          
          ! Downwelling, includes all members from top for
          ! each independant member below
          
          do jcol = 1,this%n_col(ibot)

             qp = qp + 1
             ilem = ilem_off + this%n_col(itop) + jcol
             k1 = 2*(ilem-1)+1
             k2 = k1 + 1
              
             ! Include the self terms for the current element
             ! This term is at v=0

             TAU(qp) = this%band(ib)%scco(ibot,jcol)%Ad
             OMEGA(qp,k1) = OMEGA(qp,k1) - this%band(ib)%scco(ibot,jcol)%B1d
             OMEGA(qp,k2) = OMEGA(qp,k2) - this%band(ib)%scco(ibot,jcol)%B2d

             ! We need to include the terms from
             ! all elements above the current element of interest
             ! (this can be moved out of jcol loop for efficiency)
             do icol = 1,this%n_col(itop)

                ilem = ilem_off + icol
                k1 = 2*(ilem-1)+1
                k2 = k1 + 1
                
                scelp => this%scel(itop,icol)
                sccop => this%band(ib)%scco(itop,icol)
                
                vai = scelp%lai + scelp%sai

                TAU(qp) = TAU(qp) - scelp%area * sccop%Ad *exp(-sccop%Kb*vai)
                OMEGA(qp,k1) = OMEGA(qp,k1) + scelp%area * sccop%B1d*exp(sccop%a*vai)
                OMEGA(qp,k2) = OMEGA(qp,k2) + scelp%area * sccop%B2d*exp(-sccop%a*vai)

             end do
          
          end do

          ilem_off = ilem_off + this%n_col(itop)
          
       end do do_dn_ican


       ! -------------------------------------------------------------------
       ! III. Flux equations between canopy layers, UPWELLING
       ! -------------------------------------------------------------------
       ! Refer to equation X in the technical documentation.
       ! Note the upwelling balance is performed on the upper layer,
       ! one equation for each element in the upper layer.
       ! Note that since we use "ghost elements" or air elements
       ! we don't have to factor in reflections from exposed ground.
       ! These effects will be mediated through the ghost elements
       ! -------------------------------------------------------------------
       
       ilem_off = 0
       
       do_up_ican: do ican = 2,this%n_lyr

          itop = ican-1
          ibot = ican

          do icol = 1,this%n_col(itop)
             
             qp = qp + 1
             
             ! Self terms (ie the upwelling evaluated at the bottom edge of each top element)
             ilem = ilem_off + icol
             k1   = 2*(ilem-1)+1
             k2   = k1 + 1
             scelp => this%scel(itop,icol)
             sccop => this%band(ib)%scco(itop,icol)
             
             vai = scelp%lai + scelp%sai
             TAU(qp) = sccop%Au*exp(-sccop%Kb*vai)
             OMEGA(qp,k1) = OMEGA(qp,k1) - sccop%B1u*exp(sccop%a*vai)
             OMEGA(qp,k2) = OMEGA(qp,k2) - sccop%B2u*exp(-sccop%a*vai)

             ! Terms for mean diffuse exiting lower elements (move out of this loop for efficiency)
             do jcol = 1,this%n_col(ibot)
                ilem = ilem_off + this%n_col(itop) + jcol
                k1 = 2*(ilem-1)+1
                k2 = k1 + 1
                scelp => this%scel(ibot,jcol)
                sccop => this%band(ib)%scco(ibot,jcol)
                
                TAU(qp) = TAU(qp) - scelp%area*sccop%Au
                OMEGA(qp,k1) = OMEGA(qp,k1) + scelp%area*sccop%B1u
                OMEGA(qp,k2) = OMEGA(qp,k2) + scelp%area*sccop%B2u
             end do
                
          end do

          ilem_off = ilem_off + this%n_col(itop)
       end do do_up_ican

       
    end if if_understory


    ! Flux balance equations between the understory elements, and
    ! the ground below them
    ilem_off = 0
    do ican=1,this%n_lyr-1
       ilem_off = ilem_off + this%n_col(ican)
    end do

    do jcol = 1,this%n_col(this%n_lyr)
       
       ilem = ilem_off + jcol
       qp = qp + 1
       k1 = 2*(ilem-1)+1
       k2 = k1 + 1

       scelp => this%scel(this%n_lyr,jcol)
       sccop => this%band(ib)%scco(this%n_lyr,jcol)
       
       vai = scelp%lai + scelp%sai

       TAU(qp) = sccop%Au*exp(-sccop%Kb*vai)  &
            - this%band(ib)%albedo_grnd_diff*sccop%Ad*exp(-sccop%Kb*vai) &
            - this%band(ib)%albedo_grnd_beam*sccop%Rbeam0*exp(-sccop%Kb*vai)
       
       OMEGA(qp,k1) = OMEGA(qp,k1) - sccop%B1u*exp(sccop%a*vai)
       OMEGA(qp,k2) = OMEGA(qp,k2) - sccop%B2u*exp(-sccop%a*vai)

       OMEGA(qp,k1) = OMEGA(qp,k1) + this%band(ib)%albedo_grnd_diff*sccop%B1d*exp(sccop%a*vai)
       OMEGA(qp,k2) = OMEGA(qp,k2) + this%band(ib)%albedo_grnd_diff*sccop%B2d*exp(-sccop%a*vai)
       
    end do
    
    LAMBDA = TAU
    ! Solution borrowed from G Lemieux's usage during FATES canopy trimming:
    ! Compute the optimum size of the work array
    
    lwork = -1 ! Ask dgels to compute optimal number of entries for work
    call dgels(trans, this%n_scel*2, this%n_scel*2, 1, OMEGA, this%n_scel*2, LAMBDA, this%n_scel*2, work, lwork, info)
    lwork = int(work(1)) ! Pick the optimum.  TBD, can work(1) come back with greater than work size?
    
    ! Compute the minimum of 2-norm of of the least squares fit to solve for X
    ! Note that dgels returns the solution by overwriting the LAMBDA array.
    ! The result has the form: X = [b; m]
    call dgels(trans, this%n_scel*2, this%n_scel*2, 1, OMEGA, this%n_scel*2, LAMBDA, this%n_scel*2, work, lwork, info)

    ! Save the solution terms
    ilem_off = 0
    do ican = 1,this%n_lyr
       do icol = 1,this%n_col(ican)
          ilem = ilem_off + icol
          k1 = 2*(ilem-1)+1
          k2 = k1 + 1
          scelp => this%scel(ican,icol)
          sccop => this%band(ib)%scco(ican,icol)
          sccop%lambda1 = LAMBDA(k1)
          sccop%lambda2 = LAMBDA(k2)
       end do
       ilem_off = ilem_off + this%n_col(ican)
    end do
    
    deallocate(OMEGA)
    deallocate(TAU)
    deallocate(LAMBDA)
    
    return
  end subroutine Solve


end Module TwoStreamMLPEMod
