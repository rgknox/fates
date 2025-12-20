module EDParamsDerivedSuppMod

  use iso_c_binding, only : c_char
  use iso_c_binding, only : c_int
  use iso_c_binding, only : r8 => c_double

  use FatesParameterDerivedMod, only : param_derived
  
  implicit none

  private
  public :: SetDerivedParam
  
contains
  
  ! =====================================================================================
  
  subroutine SetDerivedParam(val,pft,pname)

    real(r8), intent(in)            :: val
    character(kind=c_char,len=*), intent(in)    :: pname
    integer(kind=c_int), intent(in) :: pft

    select case(trim(pname))
    case('jmax25top')
       param_derived%jmax25top(pft,1) = val
    case('tpu25top')
       param_derived%tpu25top(pft,1) = val
    case('kp25top')
       param_derived%kp25top(pft,1) = val
    case('branch_fraction')
       param_derived%branch_frac(pft) = val   !sum(SF_val_CWD_frac(1:3))
    case default
       print*,"An unknown parameter name was sent to the parameter"
       print*,"initialization function."
       print*,"name:--",trim(pname),"--"
       stop
    end select

  end subroutine SetDerivedParam
  
end module EDParamsDerivedSuppMod

