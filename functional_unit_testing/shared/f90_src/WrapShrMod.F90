module shr_log_mod
   use iso_c_binding, only : c_char
   use iso_c_binding, only : c_int
    
   public :: shr_log_errMsg
   
 contains
   function shr_log_errMsg(source, line) result(ans)
     character(kind=c_char,len=*), intent(in) :: source
     integer(c_int), intent(in) :: line
     character(kind=c_char,len=4) :: cline  ! character version of int
     character(kind=c_char,len=128) :: ans

     write(cline,'(I4)') line
     ans = "source: " // trim(source) // " line: "// trim(cline)
     
   end function shr_log_errMsg
   
end module shr_log_mod

module shr_sys_mod

  public :: shr_sys_abort

contains

  subroutine shr_sys_abort(msg)
    character(len=*),optional :: msg
    if(present(msg))then
       write(*,*)msg
    end if
    call exit(0)
  end subroutine shr_sys_abort
  
end module shr_sys_mod

module shr_infnan_mod
  use iso_c_binding, only : c_double
  use, intrinsic :: ieee_arithmetic
  private
  real(c_double), public :: shr_infnan_nan
  public :: init_nan
  
contains
  subroutine init_nan
    shr_infnan_nan = ieee_value(shr_infnan_nan,ieee_quiet_nan)
  end subroutine init_nan
  
end module shr_infnan_mod

!use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)
