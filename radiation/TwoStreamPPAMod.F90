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

  implicit none
  private
  
  integer, parameter :: r8 = selected_real_kind(12)
  real(r8),parameter :: nearzero = 1.e-20_r8
  
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
  type, public :: twostream_type
  
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
  real(r8), parameter :: snow_scatter_vis = 0.85_r8  ! Tarboton 1995
  real(r8), parameter :: snow_scatter_nir = 0.75_r8 
  real(r8), parameter :: snow_scatter_thm = 0.65_r8  ! Tarboton 1995
  real(r8), parameter :: k_snow           = 1.0_r8  
  
  public :: ParamPrep
  
contains


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
       
       rad_params%om_leaf(ib,ft) = rad_params%rhol(ib,ft) + rad_params%taul(ib,ft)
       rad_params%om_stem(ib,ft) = rad_params%rhos(ib,ft) + rad_params%taus(ib,ft)

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
    integer :: ft    ! Plant functional type index
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
          !!this%kd_e(ilem) = this%frac_snow*k_snow + this%kd_e(ilem)

          !!this%om_e(ilem) = this%frac_snow*rad_params%om_snow(ib) + (1._r8-this%frac_snow)*this%om_e(ilem)
          
          ! Diffuse backscatter, taken from G. Bonan's code

          rho = (this%scel(ican)%col_lai(icol) * rad_params%rhol(ib,ft) + &
                 this%scel(ican)%col_sai(icol) * rad_params%rhos(ib,ft))/this%eta_e(ilem)
          tau = (this%scel(ican)%col_lai(icol) * rad_params%taul(ib,ft) + &
                 this%scel(ican)%col_sai(icol) * rad_params%taus(ib,ft))/this%eta_e(ilem)
          
          this%betad_e(ilem)  = 0.5_r8 / this%om_e(ilem) * &
               ( this%om_e(ilem) + (rho-tau) * ((1._r8+rad_params%xl(ft))/2._r8)**2._r8 )

          !!this%betad_e(ilem) = this%frac_snow*rad_params%betad_snow(ib) + (1._r8-this%frac_snow)*this%betad_e(ilem)
          
       end do
    end do
    
    return
  end subroutine CanopyPrep

  ! =============================================================================
  
  subroutine ZenithPrep(this,cosz,ib)

    ! Pre-process things that change with the zenith angle
    ! i.e. the beam optical properties

    class(twostream_type) :: this
    real(r8),intent(in)   :: cosz    ! Cosine of the zenith angle
    integer,intent(in)    :: ib
    
    integer :: ican  ! scattering element canopy layer index (top down)
    integer :: icol  ! scattering element column
    integer :: ilem  ! scattering element index
    integer :: ft    ! (plant) functional type index
    real(r8) :: asu  ! single scattering albedo
    real(r8) :: avmu ! Average inverse diffuse optical depth per unit leaf area
    real(r8) :: gdir
    real(r8) :: tmp0,tmp1,tmp2
    
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
          !!this%kb_e(ilem) = this%frac_snow*k_snow + this%kb_e(ilem)

          ! betab - upscatter parameter for direct beam radiation, from G. Bonan

          avmu = (1._r8 - rad_params%phi1(ft)/rad_params%phi2(ft) * &
               log((rad_params%phi1(ft)+rad_params%phi2(ft))/rad_params%phi1(ft))) / rad_params%phi2(ft)

          tmp0 = gdir +  rad_params%phi2(ft) * cosz
          tmp1 =  rad_params%phi1(ft) * cosz
          tmp2 = 1._r8 - tmp1/tmp0 * log((tmp1+tmp0)/tmp1)
          asu = 0.5_r8 * this%om_e(ilem) * gdir / tmp0 * tmp2
          this%betab_e(ilem) = (1._r8 + avmu*this%kb_e(ilem)) / (this%om_e(ilem)*avmu*this%kb_e(ilem)) * asu
          
          !!this%betab_e(ilem) = this%frac_snow*beta_snow + this%betab_e(ilem)
          
          
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
    integer :: ican
    
    n_scel = 0
    do ican = 1,ubound(this%scel,1)
       n_scel = n_scel + this%scel(ican)%n_col
    end do
    
  end function GetNSCel

  ! ===============================================================
  
  subroutine Solve(this,ib,Rbeam0,Rdiff0,omega_g)

    class(twostream_type) :: this
    integer,intent(in)    :: ib       ! band index
    real(r8)              :: Rbeam0   ! Intensity of beam radiation at top of canopy [W/m2]
    real(r8)              :: Rdiff0   ! Intensity of diffuse radiation at top of canopy [W/m2]
    real(r8)              :: omega_g  ! Albedo of the ground (ie fraction of radiation reflected
                                      ! in this band

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

    integer :: n_can ! number of canopy layers
    integer :: ican  ! Loop index for canopy layers
    integer :: icol  ! Loop index for canopy columns
    integer :: jcol  ! Another loop index for canopy columns
    integer :: ilem  ! Index for scattering elements

    real(r8) :: Rbeam_coltop  ! Mean beam radiation at top of columns in layer      [W/m2]
    real(r8) :: Rbeam_colbot  ! Mean beam radiation at bottom of columns in layer   [W/m2]
    real(r8) :: open_area     ! Fraction of ground that is not occupied by columns [m2/m2]

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
    ! RbeamV is the area weighted mean radiation that passes through the
    ! upper canopy.
    ! --------------------------------------------------------------------------
    
    Rbeam_coltop = Rbeam0
    ilem = 0
    do ican = 1,ubound(this%scel,1)
       Rbeam_colbot = 0._r8
       do icol = 1,this%scel(ican)%n_col
          ilem=ilem+1
          Rbeam(ilem) = Rbeam_coltop
          Rbeam_colbot = Rbeam_colbot + &
               Rbeam_coltop*this%scel(ican)%col_area(icol)*exp(-this%kb_e(ilem)*this%eta_e(ilem))
       end do
       open_area = 1._r8-sum(this%scel(ican)%col_area(1:this%scel(ican)%n_col),1)
       if(open_area>nearzero)then
          RbeamV = RbeamV + RbeamTop*open_area
       end if
       Rbeamtop = Rbeam_colbot
    end do

    ! Beam radiation at the surface
    Rbeam_sfc = Rbeam_colbot

    ! Calculate element-level intermediate terms to the solve
    ! These are dependent on leaf level scattering and beam scattering
    ! These values will be used to populate the inversion matrixes
    ! =====================================================================
    
    do ilem = 1, this%n_scel

       baf(ilem) = exp(-Kb(ilem)*this%etad_e(ilem))

       a2 = Kd_e(ilem)*Kd_e(ilem)*(om_e(ilem)-1._r8)*(om_e(ilem)-1._r8-2._r8*om_e(ilem)*bd_e(ilem))
       if(a2<0._r8) then
          print*,'a^2 is less than zero'
          stop
       end if

       a(ilem)  = sqrt(a2)
       b1 = (Kd_e(ilem)*(1._r8-om_e(ilem))* &
            (1._r8-2._r8*bb(ilem))+Kb_e(ilem))*om(k,ib)*Kb(ilem)*Rbeam(ilem)
       b2 = (Kd_e(ilem)*(om_e(ilem)-1._r8-2._r8*om_e(ilem)*betad_e(ilem))- &
            (1._r8-2._r8*betab_e(ilem))*Kb_e(ilem))*om_e(ilem)*Kb_e(ilem)*Rbeam(ilem)
       eta = a(ilem)*a(ilem) - Kb_e(ilem)*Kb_e(ilem)

       dafp(ilem) = exp(a(ilem)*etad_e(ilem))
       dafm(ilem) = exp(-a(ilem)*etad_e(ilem))

       nu_sqrd = (1._r8-om_e(ilem)+2._r8*om_e(ilem)*bd_e(ilem))/(1._r8-om_e(ilem))

       if(nu_sqrd<0._r8)then
          print('nu_sqrd is less than zero')
          stop
       end if

       
       ! B_1 term from documentation:
       B1(ilem)  = 0.5_r8*(1._r8+sqrt(nu_sqrd)) ! aka half (1 plus nu)

       ! B_2 term from documentation
       B2(ilem) = 0.5_r8*(1._r8-sqrt(nu_sqrd)) ! aka half (1 minus nu)

       ! A_2 term from documentation
       A2(ilem)    = -0.5_r8*(b2-b1)/eta       ! aka half b2 minus b1

       ! A_1 term from documentation
       A1(ilem)    = -0.5_r8*(b2+b1)/eta       ! aka half b1 plus b2

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

    
    allocate(OMEGA(2*this%n_scel,2*this%n_scel))
    allocate(TAU(2*this%n_scel))
    allocate(LAMBDA(2*this%n_scel))


    n_can = ubound(this%scel,1)
    
    ! --------------------------------------------------------------------
    ! I. Flux equations with the atmospheric boundary
    ! These balance with all elements in the upper
    ! canopy, only.  The upper canopy is layer 1.
    ! --------------------------------------------------------------------
    
    ilem = 0
    qp   = 0
    
    do icol = 1,this%scel(1)%n_col
       ilem = ilem + 1
       qp   = qp   + 1

       k1 = 2*(ilem-1)+1
       k2 = k1+1
       
       TAU(qp)      =  Rdiff0 - A2(ilem)
       OMEGA(qp,k1) =  -B2(ilem)
       OMEGA(qp,k2) =  -B1(ilem)
       
    end do

    ! -------------------------------------------------------------------
    ! II. Flux equations between canopy layers.
    ! We only perform flux balancing between layers
    ! if we have any understory, this is true if ican>1
    ! -------------------------------------------------------------------
    
    if_understory: if(n_can>1) then
       
       ! Perform downwelling flux balances
       ! Refer to Equation X in technical document
       ! ------------------------------------------------------------
       
       ! This is the index offset for the layer above the
       ! current layer of interest. We start by evaluating
       ! Layer 2, so the offset refers to layer 1, and a
       ! value of 0
       ilem_off = 0
       
       do_dn_ican: do ican = 2,n_can

          itop = ican-1  ! Top layer of the balance
          ibot = ican    ! Bottom layer of the balance
          
          ! Downwelling, includes all members from top for
          ! each independant member below
          
          do jcol = 1,this%scel(ibot)%n_col
             qp = qp + 1
             jlem = ilem_off + this%scel(itop)%n_col + jcol
             
             ! Include the self terms for the current element
             ! This term is at v=0

             TAU(qp) = TAU(qp) + A2(jlem)
             k1 = 2*(jlem-1)+1
             k2 = k1 + 1
             OMEGA(qp,k1) = OMEGA(qp,k1) + B2(jlem)
             OMEGA(qp,k2) = OMEGA(qp,k2) + B1(jlem)
             
             ! We need to include the terms from
             ! all elements above the current element of interest
             do icol = 1,this%scel(itop)%n_col
                ilem = ilem_off + icol
                TAU(qp) = TAU(qp) - &
                     this%scel(itop)%col_area(icol) * &
                     A2(ilem)*exp(-Kb(ilem)*etad(ilem))

                k1 = 2*(ilem-1)+1
                k2 = k1 + 1
                
                OMEGA(qp,k1) = OMEGA(qp,k1) - this%scel(itop)%col_area(icol) * B2(ilem)*exp(a*etad(ilem))
                OMEGA(qp,k2) = OMEGA(qp,k2) - this%scel(itop)%col_area(icol) * B1(ilem)*exp(-a*etad(ilem))

             end do
          
          end do

          ilem_off = ilem_off + this%scel(ican-1)%n_col
          
       end do do_dn_ican

       ! Perform flux balancing for upwelling radiation
       ! between layers. Refer to equation X in the
       ! technical documentation.

       
       ilem_off = 0
       do_up_ican: do ican = 2,ubound(this%scel,1)

          

          

          
          ilem_off = ilem_off + this%scel(ican-1)%n_co
       end do do_up_ican

       
    end if if_understory


    

    
    deallocate(OMEGA)
    deallocate(TAU)
    deallocate(LAMBDA)
    
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
