module shr_log_mod
  public :: shr_log_errMsg
contains

  function shr_log_errMsg(file, line)
    
    ! !INPUT/OUTPUT PARAMETERS:
    
    character(len=512)   :: shr_log_errMsg
    character(len=*), intent(in) :: file
    integer         , intent(in) :: line
    character(len=6) :: line_str

    write (line_str,"(I5)") line
    shr_log_errMsg = 'ERROR in '//trim(file)//' at line '//trim(line_str)
    
  end function shr_log_errMsg
end module shr_log_mod


module shr_sys_mod
  public :: shr_sys_abort
contains
  
  subroutine shr_sys_abort(string,rc)
    
    !----- arguments -----
    character(len=*)    , intent(in), optional :: string  ! error message string
    integer, intent(in), optional :: rc      ! error code
    character(len=512) :: local_string
    
    if (present(string)) then
       local_string = trim(string)
    else
       local_string = "Unknown error submitted to shr_abort_abort."
    end if
    call abort()
  end subroutine shr_sys_abort
  
end module shr_sys_mod
