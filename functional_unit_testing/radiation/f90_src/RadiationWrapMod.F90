module RadiationWrapMod

  use TwoStreamMLPEMod
  use iso_c_binding, only : c_char
  use iso_c_binding, only : c_int
  use iso_c_binding, only : r8 => c_double
  
  implicit none
  public
  save
  
  integer(kind=c_int), parameter :: param_string_length = 32
  
  type(twostream_type) :: twostream
  
  
contains
  
  subroutine InitAllocate(n_layer,n_column)

    integer(kind=c_int), intent(in) :: n_layer
    integer(kind=c_int), intent(in) :: n_column
    
    integer(kind=c_int) :: ican


    call twostream%AllocInitTwoStream((/1,2/),n_layer,n_column)

    
    twostream%n_lyr = n_layer

    do ican = 1,n_layer
       twostream%n_col(ican) = n_column
    end do

    twostream%force_prep = .true.

    call twostream%GetNScel()

    twostream%frac_snow = 0._r8
    twostream%frac_snow_old = 1._r8
    
    print*,"Allocated twostream instance"
    print*," with ",twostream%n_scel," elements"
    
    return
  end subroutine InitAllocate


  subroutine Dealloc()

    call twostream%DeallocTwoStream()
    
  end subroutine Dealloc

  
  subroutine SetRadParam(val,pft,ib,pname)

    real(r8), intent(in)            :: val
    character(kind=c_char,len=*), intent(in)    :: pname
    integer(kind=c_int), intent(in) :: pft
    integer(kind=c_int), intent(in) :: ib

    select case(trim(pname))
    case('rhol')
       rad_params%rhol(ib,pft) = val
    case('rhos')
       rad_params%rhos(ib,pft) = val
    case('taul')
       rad_params%taul(ib,pft) = val
    case('taus')
       rad_params%taus(ib,pft) = val
    case('xl')
       rad_params%xl(pft) = val     
    case('clumping_index')
       rad_params%clumping_index(pft) = val
    case default
       print*,"An unknown parameter name was sent to the parameter"
       print*,"initialization function."
       print*,"name:--",trim(pname),"--"
       stop
    end  select
       
  end subroutine SetRadParam

  ! =============================================================================
  
  subroutine SetGroundSnow(ib,val,pname)

    real(r8), intent(in)            :: val
    integer,  intent(in)            :: ib
    character(kind=c_char,len=*), intent(in)    :: pname

    select case(trim(pname))
    case('albedo_grnd_diff')
       twostream%band(ib)%albedo_grnd_diff = val
    case('albedo_grnd_beam')
       twostream%band(ib)%albedo_grnd_beam = val
       case default
       print*,"An unknown parameter name was sent to ground/snow"
       print*,"initialization function."
       print*,"name:--",trim(pname),"--"
       stop
    end  select
  end subroutine SetGroundSnow

  ! =============================================================================
  
  subroutine SetupCanopy(ican,icol,pft,area,lai,sai)

    integer(kind=c_int), intent(in) :: ican  ! Canopy layer index
    integer(kind=c_int), intent(in) :: icol  ! Column (pft) position index
    integer(kind=c_int), intent(in) :: pft   ! PFT index
    real(r8), intent(in)            :: area  ! columns fraction of the ground
    real(r8), intent(in)            :: lai   ! LAI
    real(r8), intent(in)            :: sai
    
    
    twostream%scel(ican,icol)%pft  = pft
    twostream%scel(ican,icol)%area = area
    twostream%scel(ican,icol)%lai  = lai
    twostream%scel(ican,icol)%sai  = sai

    return
  end subroutine SetupCanopy
  
  subroutine WrapCanopyPrep(ib)

    integer(kind=c_int),intent(in) :: ib

    call twostream%CanopyPrep(ib)
    
  end subroutine WrapCanopyPrep
  
  subroutine WrapZenithPrep(ib,cosz)

    integer,intent(in)       :: ib
    real(kind=r8),intent(in) :: cosz
    
    call twostream%ZenithPrep(ib,cosz)

    return
  end subroutine WrapZenithPrep

  subroutine WrapSolve(ib,Rbeam_atm,Rdiff_atm)

    integer(c_int),intent(in) :: ib
    real(kind=r8),intent(in)  :: Rbeam_atm
    real(kind=r8),intent(in)  :: Rdiff_atm

    call twostream%Solve(ib,Rbeam_atm,Rdiff_atm)

    return
  end subroutine WrapSolve

  subroutine WrapGetIntensity(ican,icol,ib,vai,r_diff_dn,r_diff_up,r_beam)

    integer(c_int) :: ican, icol
    integer(c_int) :: ib
    real(r8)    :: vai
    real(r8)    :: r_diff_dn
    real(r8)    :: r_diff_up
    real(r8)    :: r_beam
    
    r_diff_dn = twostream%band(ib)%scco(ican,icol)%GetRdDn(vai)
    r_diff_up = twostream%band(ib)%scco(ican,icol)%GetRdUp(vai)
    r_beam    = twostream%band(ib)%scco(ican,icol)%GetRb(vai)

    return
  end subroutine WrapGetIntensity

  subroutine WrapGetAbsRad(ican,icol,ib,vai_top,vai_bot,Rd_abs_leaf,Rb_abs_leaf,R_abs_stem)

    integer(c_int) :: ican, icol
    integer(c_int) :: ib
    real(r8)    :: vai_top,vai_bot
    real(r8)    :: Rd_abs_leaf,Rb_abs_leaf,R_abs_stem

    associate(sccop => twostream%band(ib)%scco(ican,icol), &
              scelp => twostream%scel(ican,icol) )
      
      call GetAbsRad(scelp,sccop,ib,vai_top,vai_bot,Rd_abs_leaf,Rb_abs_leaf,R_abs_stem)
      
    end associate
    return
  end subroutine WrapGetAbsRad
  
  subroutine WrapGetParams(ican,icol,ib,Kb,Kd,om,betad,betab)

    integer(c_int) :: ican, icol
    integer(c_int) :: ib
    real(r8)    :: Kb,Kd,om,betad,betab

    Kb    = twostream%band(ib)%scco(ican,icol)%Kb
    Kd    = twostream%band(ib)%scco(ican,icol)%Kd
    om    = twostream%band(ib)%scco(ican,icol)%om
    betad = twostream%band(ib)%scco(ican,icol)%betad
    betab = twostream%band(ib)%scco(ican,icol)%betab
    
    return
  end subroutine WrapGetParams

  subroutine WrapForceParams(ican,icol,ib,val,pname)

    ! This will overwrite the 2-stream parameters
    ! that are derived from the fates params

    integer(c_int) :: ican, icol
    integer(c_int) :: ib
    real(r8), intent(in)            :: val
    character(kind=c_char,len=*), intent(in)    :: pname
    
    select case(trim(pname))
    case('Kb')
       twostream%band(ib)%scco(ican,icol)%Kb = val
    case('Kd')
       twostream%band(ib)%scco(ican,icol)%Kd = val
    case('om')
       twostream%band(ib)%scco(ican,icol)%om = val
    case('betab')
       twostream%band(ib)%scco(ican,icol)%betab = val
    case('betad')
       twostream%band(ib)%scco(ican,icol)%betad = val
    case default
       print*,"An unknown parameter name was sent to the parameter"
       print*,"initialization function."
       print*,"name:--",trim(pname),"--"
       stop
    end  select
  
  end subroutine WrapForceParams
  
end module RadiationWrapMod
