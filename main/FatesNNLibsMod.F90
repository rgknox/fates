module FatesNNLibsMod

  use FatesConstantsMod            , only : r8 => fates_r8
  use, intrinsic :: iso_fortran_env, only : sp => real32
  
  implicit none
  
contains
  
  function Normalize1DArray(x_in,x_mean,x_std) result(x_out)

    ! FTorch typically uses single precision tensors
    ! so we process the input/output arrays as singles too
    
    real(sp), dimension(:), intent(in) :: x_in
    real(sp), dimension(:), intent(in) :: x_mean
    real(sp), dimension(:), intent(in) :: x_std
    integer :: i
    real(sp), dimension(size(x_in)) :: x_out
    
    do i = 1,size(x_in,dim=1)
       x_out(i) = (x_in(i) - x_mean(i)) / x_std(i)
    end do
    
  end function Normalize1DArray
  
  function DeNormalize1DArray(x_out,x_mean,x_std) result(x_in)

    ! This is just the reverse process

    real(sp), dimension(:), intent(in) :: x_out
    real(sp), dimension(:), intent(in) :: x_mean
    real(sp), dimension(:), intent(in) :: x_std
    integer :: i
    real(sp), dimension(size(x_out)) :: x_in
    
    do i = 1,size(x_out,dim=1)
       x_in(i) = x_out(i)*x_std(i)+x_mean(i)
    end do
    
  end function DeNormalize1DArray
  
  function RandSamp(vmin,vmax) result(vout)
    
    real(r8) :: vmin
    real(r8) :: vmax
    real(r8) :: vout
    real(r8) :: x
    call random_number(x)
    vout = vmin + x*(vmax-vmin)
  end function RandSamp
  
end module FatesNNLibsMod
