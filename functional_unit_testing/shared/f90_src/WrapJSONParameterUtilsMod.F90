module WrapJSONParameterUtilsMod

  ! This is a helper module that assists in communication
  ! between python modules and FATES routines.

  use JSONParameterUtilsMod
  use FatesParametersInterface, only : pstruct
  use iso_c_binding, only : c_char
  use iso_c_binding, only : c_int
  use iso_c_binding, only : r8 => c_double
  use shr_sys_mod   , only: shr_sys_abort
  
  implicit none

  private
  public :: WrapJSONSetParameter
  public :: WrapJSONRead
  
contains
  
  ! =====================================================================================

  subroutine WrapJSONRead(filename)

    character(len=*),intent(in) :: filename

    call JSONRead(filename,pstruct)

  end subroutine WrapJSONRead

  ! =====================================================================================

  subroutine WrapJSONSetParameter(param_name,index,rval,ival,sval)

    ! This is mostly used as an override, the typical assumption
    ! is that the parameters are read through a file, and that
    ! file is easily manipulated via python tools

    character(len=*),intent(in) :: param_name
    integer,intent(in)          :: index
    real(r8), optional          :: rval
    integer,  optional          :: ival
    character(len=*),optional   :: sval

    type(param_type), pointer :: param
    logical :: found_param
    logical :: type_error
    integer :: i,j,k
    integer :: nrows
    
    type_error  = .false.
    found_param = .false.
    loop_params: do k = 1,size(pstruct%parameters)
       param=>pstruct%parameters(k)
       if(trim(param_name)==trim(param%name))then
          select case(param%dtype)
          case(r_scalar_type)
             if(.not.present(rval) .or. index.ne.0)then
                type_error = .true.
                exit loop_params
             end if
             param%r_data_scalar = rval
          case(r_1d_type)
             if(.not.present(rval).or.index>size(param%r_data_1d))then
                type_error = .true.
                exit loop_params
             end if
             param%r_data_1d(index) = rval
          case(r_2d_type)
             if(.not.present(rval).or.index>size(param%r_data_2d))then
                type_error = .true.
                exit loop_params
             end if
             nrows = size(param%r_data_2d,dim=1)
             i = mod(index-1,nrows) + 1
             j = (index-1)/nrows + 1
             param%r_data_2d(i,j) = rval
          case(i_scalar_type)
             if(.not.present(ival) .or. index.ne.0)then
                type_error = .true.
                exit loop_params
             end if
             param%i_data_scalar = ival
          case(i_1d_type)
             if(.not.present(ival).or.index>size(param%i_data_1d))then
                type_error = .true.
                exit loop_params
             end if
             param%i_data_1d(index) = ival
          case(i_2d_type)
             if(.not.present(ival).or.index>size(param%i_data_2d))then
                type_error = .true.
                exit loop_params
             end if
             nrows = size(param%i_data_2d,dim=1)
             i = mod(index-1,nrows) + 1
             j = (index-1)/nrows + 1
             param%i_data_2d(i,j) = ival
          case(c_solo_type)
             if(.not.present(sval) .or. index.ne.0)then
                type_error = .true.
                exit loop_params
             end if
             param%c_data = trim(sval)
          case(c_1d_type)
             if(.not.present(sval).or.index>size(param%c_data_1d))then
                type_error = .true.
                exit loop_params
             end if
             param%c_data_1d(index) = trim(sval)
          case(c_2d_type)
             if(.not.present(sval).or.index>size(param%c_data_2d))then
                type_error = .true.
                exit loop_params
             end if
             nrows = size(param%c_data_2d,dim=1)
             i = mod(index-1,nrows) + 1
             j = (index-1)/nrows + 1
             param%c_data_2d(i,j) = trim(sval)
          case default
             write(*,*) 'WrapJSON: unknown datatype for ',trim(param_name)
             write(*,*) 'dtype: ',param%dtype
             call shr_sys_abort()
          end select
          found_param = .true.
          exit loop_params
       end if
    end do loop_params

    if(.not.found_param)then
       write(*,*) 'WrapJSON: could not find the desired parameter name:',trim(param_name)
       write(*,*) 'in the dataset'
       call shr_sys_abort()
    end if
    
    if(type_error)then
       write(*,*) 'WrapJSON: type error, wrong argument present'
       write(*,*) ' or, incorrect indexing'
       write(*,*) 'param_name:',trim(param_name)
       write(*,*) 'index:',index
       write(*,*) 'present rval: ',present(rval)
       write(*,*) 'present ival: ',present(ival)
       write(*,*) 'present sval: ',present(sval)
       call shr_sys_abort()
    end if
    
    return
  end subroutine WrapJSONSetParameter
  
end module WrapJSONParameterUtilsMod

