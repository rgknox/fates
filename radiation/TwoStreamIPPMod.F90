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

  type rad_params_type
     real(r8), allocatable :: rhol(:,:)     ! leaf reflectance: 1=vis, 2=nir    x pft
     real(r8), allocatable :: rhos(:,:)     ! stem reflectance: 1=vis, 2=nir    x pft
     real(r8), allocatable :: taul(:,:)     ! leaf transmittance: 1=vis, 2=nir  x pft
     real(r8), allocatable :: taus(:,:)     ! stem transmittance: 1=vis, 2=nir  x pft
     real(r8), allocatable :: xl(:)         ! leaf/stem orientation index, by pft
     real(r8), allocatable :: clumping_index(:) ! clumping index 0-1, when leaves stick together, by pft
  end type rad_params_type

  type(rad_params_type), allocatable :: rad_params(:)

  
  type scel_type
     integer  :: n_col
     integer, allocatable  :: col_pft(:)   ! pft index
     real(r8), allocatable :: col_area(:)  ! m2 col/m2 ground
     real(r8), allocatable :: col_lai(:)   ! m2 of leaf area / m2 col
     real(r8), allocatable :: col_sai(:)   ! m2 of stem area / m2 col
  end type scel_type


  type twostream_type
  
     type(scel_type), allocatable :: scel(:)  ! scattering elements for each layer
     real(r8)                     :: grnd_albedo 

     real(r8), allocatable        :: etad(:)  ! effective scattering area density
     real(r8), allocatable        :: Kd(:)    ! diffuse scattering area density
     real(r8), allocatable        :: bd(:)    
     real(r8), allocatable        :: om(:)    ! single scattering albedo
     
     
     
   contains

     !procedure : PrepFast
     !procedure : PrepSlow
     procedure : Solve
     
  end type twostream_type

  
contains

  subroutine ParmPrep()


    
    
    etad = zeros(nlayer,1);  % Effective Scattering Area Density
    Kd   = zeros(nlayer,1);
    bd   = zeros(nlayer,3);
    om   = zeros(nlayer,3);

    
    return
  end subroutine ParmPrep

  
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
    
    n_cl = ubound(scel,1)

    n_scel = 0
    do ic = 1,ncl
       n_scel = n_scel + scel(ic)%n_col
    end do


    allocate(omega(2*n_scel,2*n_scel))
    allocate(tau(2*n_scel,1))
    allocate(lambda(2*n_scel,1))

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
