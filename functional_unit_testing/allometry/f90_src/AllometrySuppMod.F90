module AllometrySuppMod

  use iso_c_binding, only : c_char
  use iso_c_binding, only : c_int
  use iso_c_binding, only : r8 => c_double

  use PRTParametersMod, only : prt_params
  implicit none
  public
  save
  
  integer(kind=c_int), parameter :: param_string_length = 32
  logical, parameter :: debug = .true.

  ! We do not use the param_derived_type from the main
  ! FATES codebase, ie main/FatesParameterDerivedMod.F90
  ! This is because it has non-trivial dependencies

contains

  ! =============================================================
  
  subroutine AllocAllomParam(numpft)

    integer, intent(in) :: numpft
    
    allocate(prt_params%allom_d2h1(numpft))
    allocate(prt_params%allom_d2h2(numpft))
    allocate(prt_params%allom_d2h3(numpft))
    allocate(prt_params%allom_hmode(numpft))
    allocate(prt_params%allom_dbh_maxheight(numpft))
    allocate(prt_params%allom_agb1(numpft))
    allocate(prt_params%allom_agb2(numpft))
    allocate(prt_params%allom_agb3(numpft))
    allocate(prt_params%allom_agb4(numpft))
    allocate(prt_params%wood_density(numpft))
    allocate(prt_params%c2b(numpft))
    allocate(prt_params%allom_agb_frac(numpft))
    allocate(prt_params%allom_amode(numpft))
    allocate(prt_params%slatop(numpft))
    allocate(prt_params%allom_lmode(numpft))
    allocate(prt_params%allom_d2bl1(numpft))
    allocate(prt_params%allom_d2bl2(numpft))
    allocate(prt_params%allom_d2bl3(numpft))
    allocate(prt_params%allom_blca_expnt_diff(numpft))
    allocate(prt_params%allom_d2ca_coefficient_min(numpft))
    allocate(prt_params%allom_d2ca_coefficient_max(numpft))
    allocate(prt_params%leafn_vert_scaler_coeff1(numpft))
    allocate(prt_params%leafn_vert_scaler_coeff2(numpft))
    allocate(prt_params%slamax(numpft))
    allocate(prt_params%allom_sai_scaler(numpft))
    allocate(prt_params%allom_smode(numpft))
    allocate(prt_params%allom_cmode(numpft))
    allocate(prt_params%allom_fmode(numpft))
    allocate(prt_params%allom_stmode(numpft))
    allocate(prt_params%cushion(numpft)) 
    allocate(prt_params%allom_la_per_sa_int(numpft))
    allocate(prt_params%allom_la_per_sa_slp(numpft))
    allocate(prt_params%allom_h2cd1(numpft))
    allocate(prt_params%allom_h2cd2(numpft))
    allocate(prt_params%allom_dmode(numpft))
    allocate(prt_params%fnrt_prof_mode(numpft))
    allocate(prt_params%fnrt_prof_a(numpft))
    allocate(prt_params%fnrt_prof_b(numpft))
    allocate(prt_params%woody(numpft))

    return
  end subroutine AllocAllomParam

  subroutine DeallocAllomParam()

    deallocate(prt_params%allom_d2h1)
    deallocate(prt_params%allom_d2h2)
    deallocate(prt_params%allom_d2h3)
    deallocate(prt_params%allom_hmode)
    deallocate(prt_params%allom_dbh_maxheight)
    deallocate(prt_params%allom_agb1)
    deallocate(prt_params%allom_agb2)
    deallocate(prt_params%allom_agb3)
    deallocate(prt_params%allom_agb4)
    deallocate(prt_params%wood_density)
    deallocate(prt_params%c2b)
    deallocate(prt_params%allom_agb_frac)
    deallocate(prt_params%allom_amode)
    deallocate(prt_params%slatop)
    deallocate(prt_params%allom_lmode)
    deallocate(prt_params%allom_d2bl1)
    deallocate(prt_params%allom_d2bl2)
    deallocate(prt_params%allom_d2bl3)
    deallocate(prt_params%allom_blca_expnt_diff)
    deallocate(prt_params%allom_d2ca_coefficient_min)
    deallocate(prt_params%allom_d2ca_coefficient_max)
    deallocate(prt_params%leafn_vert_scaler_coeff1)
    deallocate(prt_params%leafn_vert_scaler_coeff2)
    deallocate(prt_params%slamax)
    deallocate(prt_params%allom_sai_scaler)
    deallocate(prt_params%allom_smode)
    deallocate(prt_params%allom_cmode)
    deallocate(prt_params%allom_fmode)
    deallocate(prt_params%allom_stmode)
    deallocate(prt_params%cushion) 
    deallocate(prt_params%allom_la_per_sa_int)
    deallocate(prt_params%allom_la_per_sa_slp)
    deallocate(prt_params%allom_h2cd1)
    deallocate(prt_params%allom_h2cd2)
    deallocate(prt_params%allom_dmode)
    deallocate(prt_params%fnrt_prof_mode)
    deallocate(prt_params%fnrt_prof_a)
    deallocate(prt_params%fnrt_prof_b)
    deallocate(prt_params%woody)
    
  end subroutine DeallocAllomParam

  ! =====================================================================================

  function IsAllomParamAllocated() result(is_allocated)
    
    logical :: is_allocated 

    ! use c3psn to test allocated/deallocated status
    if(allocated(prt_params%woody))then
       is_allocated = .true.
    else
       is_allocated = .false.
    end if
       
  end function IsAllomParamAllocated
  
  ! =====================================================================================
  
  subroutine SetAllomParam(val,pft,pname)

    real(r8), intent(in)            :: val
    character(kind=c_char,len=*), intent(in)    :: pname
    integer(kind=c_int), intent(in) :: pft

    select case(trim(pname))
    case('fates_allom_d2h1')
       prt_params%allom_d2h1(pft) = val
    case('fates_allom_d2h2')
       prt_params%allom_d2h2(pft) = val
    case('fates_allom_d2h3')
       prt_params%allom_d2h3(pft) = val
    case('fates_allom_hmode')
       prt_params%allom_hmode(pft) = nint(val)
    case('fates_allom_dbh_maxheight')
       prt_params%allom_dbh_maxheight(pft) = val
    case('fates_allom_agb1')
       prt_params%allom_agb1(pft) = val
    case('fates_allom_agb2')
       prt_params%allom_agb2(pft) = val
    case('fates_allom_agb3')
       prt_params%allom_agb3(pft) = val
    case('fates_allom_agb4')
       prt_params%allom_agb4(pft) = val
    case('fates_wood_density')
       prt_params%wood_density(pft) = val
    case('fates_c2b')
       prt_params%c2b(pft) = val
    case('fates_allom_agb_frac')
       prt_params%allom_agb_frac(pft) = val
    case('fates_allom_amode')
       prt_params%allom_amode(pft) = nint(val)
    case('fates_leaf_slatop')
       prt_params%slatop(pft) = val
    case('fates_leaf_slamax')
       prt_params%slamax(pft) = val
    case('fates_allom_lmode')
       prt_params%allom_lmode(pft) = nint(val)
    case('fates_allom_d2bl1')
       prt_params%allom_d2bl1(pft) = val
    case('fates_allom_d2bl2')
       prt_params%allom_d2bl2(pft) = val
    case('fates_allom_d2bl3')
       prt_params%allom_d2bl3(pft) = val
    case('fates_allom_blca_expnt_diff')
       prt_params%allom_blca_expnt_diff(pft) = val
    case('fates_allom_d2ca_coefficient_min')
       prt_params%allom_d2ca_coefficient_min(pft) = val
    case('fates_allom_d2ca_coefficient_max')
       prt_params%allom_d2ca_coefficient_max(pft) = val
    case('fates_leafn_vert_scaler_coeff1')
       prt_params%leafn_vert_scaler_coeff1(pft) = val
    case('fates_leafn_vert_scaler_coeff2')
       prt_params%leafn_vert_scaler_coeff2(pft) = val
    case('fates_allom_sai_scaler')
       prt_params%allom_sai_scaler(pft) = val
    case('fates_allom_smode')
       prt_params%allom_smode(pft) = nint(val)
    case('fates_allom_cmode')
       prt_params%allom_cmode(pft) = nint(val)
    case('fates_allom_fmode')
       prt_params%allom_fmode(pft) = nint(val)
    case('fates_allom_stmode')
       prt_params%allom_stmode(pft) = nint(val)
    case('fates_alloc_storage_cushion')
       prt_params%cushion(pft) = val 
    case('fates_allom_la_per_sa_int')
       prt_params%allom_la_per_sa_int(pft) = val
    case('fates_allom_la_per_sa_slp')
       prt_params%allom_la_per_sa_slp(pft) = val
    case('fates_allom_h2cd1')
       prt_params%allom_h2cd1(pft) = val
    case('fates_allom_h2cd2')
       prt_params%allom_h2cd2(pft) = val
    case('fates_allom_dmode')
       prt_params%allom_dmode(pft) = nint(val)
    case('fates_allom_fnrt_prof_mode')
       prt_params%fnrt_prof_mode(pft) = nint(val)
    case('fates_allom_fnrt_prof_a')
       prt_params%fnrt_prof_a(pft) = val
    case('fates_allom_fnrt_prof_b')
       prt_params%fnrt_prof_b(pft) = val
    case('fates_woody')
       prt_params%woody(pft) = val
    case default
       print*,"An unknown parameter name was sent to the parameter"
       print*,"initialization function."
       print*,"name:--",trim(pname),"--"
       stop
    end  select
       
  end subroutine SetAllomParam
  
  subroutine DumpAllomParams()

    print*,'woody:',prt_params%woody(:)
    
  end subroutine DumpAllomParams
  
end module AllometrySuppMod
