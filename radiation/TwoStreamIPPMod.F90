Module TwoStreamIPAMod

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
  ! Assumptions: band index 1 = visible (vis)
  !                         2 = near infrared (nir)
  !                         3 = thermal
  !
  
  type rad_params_type

     ! From the parameter file
     real(r8), allocatable :: rhol(:,:)         ! leaf reflectance:   (band x pft)
     real(r8), allocatable :: rhos(:,:)         ! stem reflectance:   (band x pft)
     real(r8), allocatable :: taul(:,:)         ! leaf transmittance: (band x pft)
     real(r8), allocatable :: taus(:,:)         ! stem transmittance: (band x pft)
     real(r8), allocatable :: xl(:)             ! leaf/stem orientation (pft)
     real(r8), allocatable :: clumping_index(:) ! clumping index 0-1, when
                                                ! leaves stick together (pft)

     ! Derived
     real(r8), allocatable :: phi1(:)        ! intermediate term for kd and kb
     real(r8), allocatable :: phi2(:)        ! intermediate term for kd and kb
     real(r8), allocatable :: kd_leaf(:)     ! Mean optical depth per unit area leaves in diffuse
     real(r8), allocatable :: kd_stem(:)     ! Mean optical depth per unit area stems in diffuse
     real(r8), allocatable :: betad_leaf(:)  ! Diffuse backscatter fraction for leaves
     real(r8), allocatable :: betad_stem(:)  ! Diffuse backscatter fraction for stems
     real(r8), allocatable :: om_snow(:)     ! Snow scattering albedo for snow (band)
     real(r8), allocatable :: om_leaf(:,:)   ! Leaf scattering albedo (band x pft)
     real(r8), allocatable :: om_stem(:,:)   ! Stem scattering albedo (band x pft)

     
     !real(r8), allocatable :: betab(:,:)
     

     
  end type rad_params_type

  type(rad_params_type) :: rad_params

  
  type scel_type
     integer  :: n_col
     integer, allocatable  :: col_pft(:)   ! pft index
     real(r8), allocatable :: col_area(:)  ! m2 col/m2 ground
     real(r8), allocatable :: col_lai(:)   ! m2 of leaf area / m2 col
     real(r8), allocatable :: col_sai(:)   ! m2 of stem area / m2 col
  end type scel_type


  type twostream_type
  
     type(scel_type), allocatable :: scel(:)  ! scattering elements for each layer
     integer                      :: n_scel
     real(r8)                     :: grnd_albedo 

     real(r8), allocatable        :: eta_e(:)   ! effective scattering area density
     real(r8), allocatable        :: kd_e(:)    ! mean optical depth per unit scattering
                                                ! area of media in the element
                                                ! weighted combination of kd_leaf and kd_stem
     real(r8), allocatable        :: bd(:)    
     real(r8), allocatable        :: om_e(:)    ! single scattering albedo
     
     
     
   contains

     !procedure : PrepFast
     !procedure : PrepSlow
     procedure : Solve
     
  end type twostream_type


  ! Assumptions 1=visible (vis)
  !             2=near infrared (nir)
  !             3=thermal 
  

  ! Constants or reasonable approximations thereof
  real(r8), parameter :: snow_scatter_vis = 0.8_r8
  real(r8), parameter :: snow_scatter_nir = 0.4_r8
  real(r8), parameter :: snow_scatter_thm = 0.4_r8
  real(r8), parameter :: kd_snow          = 1.0_r8  

  
contains


  subroutine ParamPrep(numpft)

    real(r8) :: mu

    do ft = 1,numpft

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

       !rad_params%betad_leaf(ib,ft) = rad_params%rhol(ib,ft)/rad_params%om_leaf(ib,ft)
       !rad_params%betad_stem(ib,ft) = rad_params%rhos(ib,ft)/rad_params%om_stem(ib,ft)

      


       
    end do
    

  end subroutine ParamPrep

  
  subroutine CanopyPrep(this,ib)

    ! Pre-process things that change with canopy-geometry
    
    class(twostream_type) :: this

    integer :: ie  ! scattering element index
    integer :: ic  ! canopy layer index (top down)
    
    ie = 0
    do ican = 1,ubound(this%scel,1)

       do icol = 1,this%scel(ican)%n_col

          ie = ie + 1
          ft = this%scel(ican)%col_pft(icol)
          
          this%eta_e(ie) = this%scel(ican)%col_lai(icol) + this%scel(ican)%col_sai(icol)

          ! Mean element transmission coefficients w/o snow effects
          this%kd_e(ie) =  (this%scel(ican)%col_lai(icol) * rad_params%kd_leaf(ft) + &
                            this%scel(ican)%col_sai(icol) * rad_params%kd_stem(ft))/this%eta_e(ie)
          
          this%om_e(ie) =  (this%scel(ican)%col_lai(icol)*rad_params%om_leaf(ib,ft) + &
                            this%scel(ican)%col_sai(icol)*rad_params%om_stem(ib,ft))/this%eta_e(ie)

          !betad_e(ie) = (this%scel(ican)%col_lai(icol)*rad_params%betad_leaf(ib,ft) + &
          !               this%scel(ican)%col_sai(icol)*rad_params%betad_stem(ib,ft))/this%eta_e(ie)

          ! Mean element transmission coefficients with snow effects
          this%kd_e(ie) = fsnow*kd_snow + (1._r8-fsnow)*this%kd_e(ie)

          this%om_e(ie) = fsnow*rad_params%om_snow(ib) + (1._r8-fsnow)*this%om_e(ie)
          
          ! Diffuse backscatter, taken from G. Bonan's code

          rho = (this%scel(ican)%col_lai(icol) * rad_params%rho_leaf(ft) + &
                 this%scel(ican)%col_sai(icol) * rad_params%rho_stem(ft))/this%eta_e(ie)
          tau = (this%scel(ican)%col_lai(icol) * rad_params%tau_leaf(ft) + &
                 this%scel(ican)%col_sai(icol) * rad_params%tau_stem(ft))/this%eta_e(ie)
          
          this%betad_e(ie)  = 0.5_r8 / this%om_e(ie) * &
               ( this%om_e(ie) + (rho-tau) * ((1._r8+rad_params%xl(ft))/2._r8)**2._r8 )

          this%betad_e(ie) = fsnow*rad_params%betad_snow(ib) + (1._r8-fsnow)*this%betad_e(ie)
          
          
          
       end do
    end do
    
    return
  end subroutine CanopyPrep

  ! =============================================================================
  
  subroutine ZenithPrep(cosz,ib)

    ! Pre-process things that change with the zenith angle
    ! i.e. the beam optical properties
        
    ie = 0
    do ican = 1,ubound(this%scel,1)

       do icol = 1,this%scel(ican)%n_col

          ie = ie + 1
          ft = this%scel(ican)%col_pft(icol)

          gdir = rad_params%phi1(ft) + rad_params%phi2(ft) * cosz
          
          !how much direct light penetrates a singleunit of lai?
          this%kb_e(ie) = gdir / cosz

          ! betab - upscatter parameter for direct beam radiation, from G. Bonan

          avmu = (1._r8 - rad_params%phi1(ft)/rad_params%phi2(ft) * &
               log((rad_params%phi1(ft)+rad_params%phi2(ft))/rad_params%phi1(ft))) / rad_params%phi2(ft)

          tmp0 = gdir + phi2 * cosz
          tmp1 = phi1 * cosz
          tmp2 = 1._r8 - tmp1/tmp0 * log((tmp1+tmp0)/tmp1)
          asu = 0.5_r8 * this%om_e(ie) * gdir / tmp0 * tmp2
          this%betab_e(ie) = (1._r8 + avmu*this%kb_e(ie)) / (this%om_e(ie)*avmu*this%kb_e(ie)) * asu
       
    end do !FT


  
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
    integer :: ie  ! Loop index for scattering elements

    associate(n_scel => this%n_scel)
    
    allocate(omega(2*this%n_scel,2*this%n_scel))
    allocate(tau(2*this%n_scel,1))
    allocate(lambda(2*this%n_scel,1))

    baf
    om
    Kb
    etad
    bd
    bb
    b1
    
    

    do ie = 1,n_scel
       baf(k) = exp(-Kb(k)*etad(k))    ! Beam attenuation fraction
       a2 = Kd(k)*Kd(k) *(om(k,ib)-1)*(om(k,ib)-1-2*om(k,ib)*bd(k,ib));
       if(a2<0)
       display('a^2 is less than zero');
       pause;
       return;
    end do
        a(k)  = sqrt( a2  );
        b1 = (Kd(k)*(1-om(k,ib))*(1-2*bb(k))+Kb(k))*om(k,ib)*Kb(k)*Rbeam(k);
        b2 = (Kd(k)*(om(k,ib)-1-2*om(k,ib)*bd(k,ib))-(1-2*bb(k))*Kb(k))*om(k,ib)*Kb(k)*Rbeam(k);
        eta = a(k).^2 - Kb(k).^2;
        dafp(k) = exp(a(k)*etad(k));
        dafm(k) = exp(-a(k)*etad(k));
%        b1v(k) = b1;
%        b2v(k) = b2;
%        etav(k) = eta;
        nu_sqrd = (1-om(k,ib)+2.*om(k,ib)*bd(k,ib))./(1-om(k,ib));
        if(numel(nu_sqrd)~=1)
            display('NU^2 is not scalar?');
            display(size(om(k,ib)));
            display(size(bd(k,ib)));
            pause;
            return;
        end
        if(nu_sqrd(1)<0)
            display('nu_sqrd is less than zero');
            pause;
            return;
        end
        nu_hplus(k)  = 0.5*(1+sqrt(nu_sqrd)); % aka half (1 plus nu)
        nu_hminus(k) = 0.5*(1-sqrt(nu_sqrd)); % aka half (1 minus nu)
        hb2mb1(k)    = 0.5*(b2-b1)/eta;       % aka half b2 minus b1
        hb2pb1(k)    = 0.5*(b2+b1)/eta;       % aka half b1 plus b2
       

    end do

    


    deallocate(omega)
    deallocate(tau)
    deallocate(lambda)
    
    
    real(r8) :: Rbeam(this%n_element)

    
    ! Beam Scattering
    ! =====================================================================
    ! First do the direct beam stuff.  It is a trivial solution
    ! and is required as a boundary condition to the diffuse solver
    ! All parallel layers plus 1 (the mixed) recieve downwelling from the
    ! atmosphere.  Rbeam0 is the upper boundary condition provided by data
    ! or another model.
    ! Rbeam() is the incident beam radiation at the top of each layer
    ! RbeamV is the area weighted mean radiation that passes through the
    ! upper canopy.
    ! =====================================================================


    Rbeam(:) = spval
    RbeamV    = 0
    sumfrac   = 0
    ak        = 0

    do ic = 1,numcl

       do ip = 

       

       
    do k = 1,this%nparlayer
       ak = k
       Rbeam(ak)   = Rbeam0
       RbeamV      = RbeamV + this%area_frac(ak)*Rbeam(ak)*exp(-Kb(ak)*etad(ak))
       sumfrac     = sumfrac + afrac(ak)
    end do
    
    

    ! =====================================================================
    ! First lets see if we can auto-generate the shape
    ! The rows are the different equations, the columns are the different
    ! Coefficients, for two-stream equations, there are two
    ! =====================================================================
    
    ! Indexing Variables
    ! ep : element position
    ! k1 : coefficient 1 position
    ! k2 : coefficient 2 position
    ! qp : equation position, this continues to increment
    
    
    ! Upper boundary condition for the parallel layers
    ! Known downwelling diffuse radiation
    ! Upper boundary condition
    ! Downwelling diffuse equals that provided by atm model or data
    ! =====================================================================
    
    qp=1
    do k=1,nparlayer
        ep = k
        k1 = 2*(ep-1)+1
        k2 = k1+1
        TAU(qp)      =  -hb2mb1(ep) + R_emiss(ep,ib) - Rdiff0(ib)
        OMEGA(qp,k1) =  nu_hminus(ep)
        OMEGA(qp,k2) =  nu_hplus(ep)
        qp=qp+1
     end do

    


    return
  end subroutine TwoStreamDrive

  ! ============================================================
  

  subroutine CanopyToScatterEl()

    ! In this routine, we decompose a FATES canopy
    ! into the data structures upon which we
    ! can perform the two-stream calculations
    
    

    

    return
  end subroutine CanopyToScatterEl



  
  

  

  


  
end Module TwoStreamIPAMod
