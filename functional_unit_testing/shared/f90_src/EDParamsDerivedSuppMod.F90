module FatesParameterDerivedMod

  use iso_c_binding, only : c_char
  use iso_c_binding, only : c_int
  use iso_c_binding, only : r8 => c_double

  implicit none
  
  type, public :: param_derived_type
     real(r8), allocatable :: branch_frac(:)  ! fraction of aboveground woody biomass in branches (as
                                              ! oppose to stems) - for use in damage allometries
  end type param_derived_type

  type(param_derived_type),public :: param_derived          ! Instantiation of the parameter object

  private
  public :: AllocDerivedParam
  public :: IsDerivedParamAllocated
  public :: SetDerivedParam
  
contains
  
  subroutine AllocDerivedParam(numpft)

    integer, intent(in) :: numpft
    allocate(param_derived%branch_frac(numpft))
    
  end subroutine AllocDerivedParam

  function IsDerivedParamAllocated() result(is_allocated)
    
    logical :: is_allocated 

    if(allocated(param_derived%branch_frac))then
       is_allocated = .true.
    else
       is_allocated = .false.
    end if
    
  end function IsDerivedParamAllocated
  
  ! =====================================================================================
  
  subroutine SetDerivedParam(val,pft,pname)

    real(r8), intent(in)            :: val
    character(kind=c_char,len=*), intent(in)    :: pname
    integer(kind=c_int), intent(in) :: pft

    select case(trim(pname))
    case('branch_fraction')
       param_derived%branch_frac(pft) = val   !sum(SF_val_CWD_frac(1:3))
    case default
       print*,"An unknown parameter name was sent to the parameter"
       print*,"initialization function."
       print*,"name:--",trim(pname),"--"
       stop
    end select

  end subroutine SetDerivedParam
  
end module FatesParameterDerivedMod



module EDParamsMod

  use iso_c_binding, only : c_char
  use iso_c_binding, only : c_int
  use iso_c_binding, only : r8 => c_double
  use FatesConstantsMod, only : min_vai_bin_sum
  
  implicit none

  integer                      :: nlevleaf 
  real(r8), allocatable,public :: dinc_vai(:)
  real(r8), allocatable,public :: dlower_vai(:)
  real(r8), allocatable,public :: ED_val_history_damage_bin_edges(:)

  private
  public :: AllocateDincVai
  public :: SetDincVai
  
contains

  subroutine AllocateDincVai(nlevleaf_in)

    integer, intent(in) :: nlevleaf_in

    nlevleaf = nlevleaf_in
    allocate(dinc_vai(nlevleaf))
    allocate(dlower_vai(nlevleaf))

  end subroutine AllocateDincVai

  ! =====================================================================================
  
  subroutine SetDincVai(vai_top_bin_width,vai_width_increase_factor)

    real(r8),intent(in) :: vai_top_bin_width
    real(r8),intent(in) :: vai_width_increase_factor
    integer :: i
    
    do i = 1,nlevleaf
       dinc_vai(i) = vai_top_bin_width * vai_width_increase_factor ** (i-1)
    end do
    
    if (sum(dinc_vai) < min_vai_bin_sum ) then
       print*,'You messed up your dinc_vai'
       stop
    end if
    
    ! lower edges of VAI bins
    dlower_vai(1) = 0._r8
    do i = 2,nlevleaf
       dlower_vai(i) =  dlower_vai(i-1) + dinc_vai(i-1)
    end do

  end subroutine SetDincVai
    
end module EDParamsMod
