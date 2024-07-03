module FatesTestPhotosynthesisMod
  ! Helper methods for testing the FATES photosynthesis routines
  
  use FatesConstantsMod, only : r8 => fates_r8
  
  implicit none
  
  contains 
  
  subroutine CalcVaporPressure(temperature, vapor_pressure_deficit, sat_vapor_pressure,  &
      vapor_pressure)
    !
    ! DESCRIPTION:
    ! Calculates saturation vapor pressure and actual vapor pressure based on input
    ! temperature and vapor pressure deficit
    !
    
    ! ARGUMENTS:
    real(r8), intent(in)  :: temperature            ! temperature [K]
    real(r8), intent(in)  :: vapor_pressure_deficit ! vapor pressure deficit [Pa]
    real(r8), intent(out) :: sat_vapor_pressure     ! saturation vapor pressure [Pa]
    real(r8), intent(out) :: vapor_pressure         ! vapor pressure [Pa]
    
    sat_vapor_pressure = sat_vapor_press(temperature - 273.15_r8) ! function uses temperature in Celsius
    
    vapor_pressure = max(0.0_r8, sat_vapor_pressure - vapor_pressure_deficit)
    
  end subroutine CalcVaporPressure
  
  ! --------------------------------------------------------------------------------------
  
  real(r8) function sat_vapor_press(temperature)
    !
    ! DESCRIPTION:
    ! Calculates saturation vapor pressure and actual vapor pressure based on input
    ! temperature and vapor pressure deficit
    ! Equation and constants taken from Bonan 2019, Climate Change and Terrestrial 
    ! Ecosystem Modeling, pg 47-48
    !
    
    ! ARGUMENTS:
    real(r8), intent(in) :: temperature ! temperature [degrees C]
    
    ! LOCALS:
    real(r8), pointer :: a(0:8)
    
    ! CONSTANTS:
    real(r8), parameter :: a_water(0:8) =                                                &
      (/6.11213476_r8, 0.444007856_r8, 0.143064234E-1_r8, 0.264461437E-3_r8,             &
      0.305903558E-5_r8, 0.196237241E-7_r8, 0.892344772E-10_r8, -0.373208410E-12_r8,     &
      0.209339997E-15_r8/)
    
    real(r8), parameter :: a_ice(0:8) =                                                  &
      (/6.11123516_r8, 0.503109514_r8, 0.188369801E-1_r8, 0.420547422E-3_r8,             &
      0.614396778E-5_r8, 0.602780717E-7_r8, 0.387940929E-9_r8, 0.149436277E-11_r8,       &
      0.262655803E-14_r8/)
      
    if (temperature >= 0.0_r8 .and. temperature <= 100.0_r8) then 
      a => a_water
    else if (temperature >= -75.0_r8 .and. temperature < 0.0_r8) then
      a => a_ice
    else 
      print *, "Temperature ", temperature, " is outside the range that we can calculate."
      stop
    end if
    
    sat_vapor_press = 100.0*(a(0) + temperature*(a(1) + temperature*(a(2) +                &
      temperature*(a(3) + temperature*(a(4) + temperature*(a(5) + temperature *(a(6) +   &
      temperature*(a(7) + temperature*a(8)))))))))
  
  end function sat_vapor_press
  
  
end module FatesTestPhotosynthesisMod