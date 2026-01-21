module FatesKineticsMod

  use FatesConstantsMod, only : r8 => fates_r8
  use shr_log_mod,       only : errMsg => shr_log_errMsg
  use FatesGlobals,      only : endrun => fates_endrun
  use FatesGlobals,      only : fates_log
  
  implicit none
  private
  
  public :: AqueousUptakeMM
  public :: SolveConcNewtRaphMM
  
contains

  subroutine AqueousUptakeMM (ncomp, s0, v_max, k_m, surf_exp, dt, theta, s_final, mass_flux)

    ! -------------------------------------------------------------------------
    ! This module represents competitive uptake of a nutrient substrate
    ! for an arbitrary number of competitors. Uptake
    ! is based on Michaelis-Menten kinetics. We perform an approximation
    ! of this system by creating a single set of effective uptake
    ! parameters vmax and k_m.  We integrate the system to find the total
    ! change in nutrient concentration, and then partition the flux
    ! to the competitors based on their uptake proportions at the
    ! beginning of the time-step.
    ! --------------------------------------------------------------------------
    
    integer,  intent(in) :: ncomp                ! Number of competitors
    real(r8), intent(in) :: s0                   ! Initial nutrient concentration (g/m3 h2o)
    real(r8), intent(in) :: v_max(ncomp)         ! Max uptake rate at saturation (g/m3 h2o/t)
    real(r8), intent(in) :: k_m(ncomp)           ! Half-saturation constants (g/m3 h2o)
    real(r8), intent(in) :: surf_exp(ncomp)      ! Fraction of the competitor surface that is 
                                                 ! exposed to the substrate (0-1)
    real(r8), intent(in) :: dt                   ! Time step (s)
    real(r8), intent(in) :: theta                ! water content (m3 h20/m3 soil)
    real(r8), intent(out) :: s_final             ! Final concentration (g/m3 h2o)
    real(r8), intent(out) :: mass_flux(ncomp)    ! Grams moved per competitor (g/m3 soil/s)
    
    real(r8) :: v_max_eff           ! Effective vmax of all competitors combined (sum)
                                    ! (g/m3 h2o/s)
    real(r8) :: k_m_eff             ! Effective K_m of all competitors (harmonic mean)
    real(r8) :: mass_flux_tot       ! grams of nutrient moved (g/m3 soil/s)
    real(r8) :: initial_flux_total  ! total substrate flux at t0
    real(r8) :: initial_flux(ncomp) ! competitor substrate flux at t0
    integer  :: i                   ! competitor loop index
    
    ! 1. Calculate Effective Parameters
    v_max_eff = sum(surf_exp * v_max) 
    k_m_eff = v_max_eff / sum(v_max * surf_exp / k_m)
    
    ! 2. Solve for S_final using Newton's method
    call SolveConcNewtRaphMM(s0, v_max_eff, k_m_eff, dt, s_final)
    
    ! 3. Calculate total mass flux rate into all competitors
    mass_flux_tot = (s0 - s_final) * theta / dt
    
    ! 4. Partition based on initial proportions
    initial_flux = v_max * surf_exp * s0 / (k_m + s0)
    initial_flux_total = sum(initial_flux)

    ! 5. Apply fraction from initial flux to the time integrated flux
    mass_flux = mass_flux_tot * initial_flux/initial_flux_total

    return
  end subroutine AqueousUptakeMM

  ! ================================================================

  subroutine SolveConcNewtRaphMM(s0, vmax, km, dt, s_next)

    ! ---------------------------------------------------
    ! This routine finds the concentration of a substrate
    ! at the end of a time-step dt.
    ! We do this by integrating the differential equation
    ! that defines the change in substrate concentration
    ! following a classic Michaeles-Menten formulation.
    ! For change in concentration S, with max uptake
    ! and half-saturation constants: Vmax and Km:
    !
    ! dS/dt = - Vmax*S/(Km+S)
    !
    ! (Km/S + 1) dS = -Vmax * dt
    !
    !  Integrate both sides, from S0 to Sf (solving for S "final"):
    !  int_S0^Sf (Km/S + 1) ds = int_0^t_f - Vmax dt
    !
    !  Km * ln(Sf) + Sf - (Km ln(S0) + S0 ) = -Vmax * deltat
    !  Km * ln(Sf/S0) + Sf - S0 = -Vmax * deltat
    !
    ! Lets find the root:
    !  f(Sf) = 0 = (S0 - Sf) + Km * ln(S0/Sf) - Vmax*deltat
    !
    ! Derivative:
    !  f'(Sf) = -1 - (Km/Sf)
    ! ----------------------------------------------------------------------
    
    implicit none
    
    ! Input variables
    real(r8), intent(in)  :: s0      ! Initial concentration [mass/volume]
    real(r8), intent(in)  :: vmax    ! Max uptake rate [mass/volume/time]
    real(r8), intent(in)  :: km      ! Half-saturation constant [mass/volume]
    real(r8), intent(in)  :: dt      ! Time step [time]
    
    ! Output variable
    real(r8), intent(out) :: s_next
    
    ! Internal parameters
    integer, parameter :: max_iter = 20
    real(r8), parameter :: tol = 1.0d-10
    real(r8)            :: f, df, s_curr
    integer            :: i

    ! 1. Check if concentration is already zero or negligible
    if (s0 <= 1.0d-15) then
        s_next = 0.0d0
        return
    endif

    ! 2. Check for total depletion (preventing log of negative)
    ! If total potential uptake > current substrate + some buffer,
    ! it's likely S will hit zero.
    if ( (vmax * dt) > (s0 + km) ) then
        ! Optional: More complex check, but usually safety-first for soil:
        s_curr = s0 * 0.1d0 ! Start with a much smaller guess
    else
        s_curr = s0 - (vmax * s0 / (km + s0)) * dt ! Forward Euler guess
    endif

    ! 3. Newton-Raphson Iteration
    do i = 1, max_iter
        ! f(s) = (s0 - s) + km * ln(s0/s) - vmax*dt
        f = (s0 - s_curr) + km * log(s0 / s_curr) - (vmax * dt)
        
        ! f'(s) = -1 - (km/s)
        df = -1.0d0 - (km / s_curr)
        
        ! Newton step: s_new = s_old - f/f'
        s_next = s_curr - (f / df)
        
        ! Safety: Substrate cannot be negative or larger than s0
        if (s_next <= 0.0d0) s_next = 1.0d-18 
        if (s_next > s0) s_next = s0

        ! Convergence check
        if (abs(s_next - s_curr) < tol) return
        
        s_curr = s_next
    end do

    ! If it reaches here, it didn't converge (rare for this function)
    ! s_next is already set to the last best guess.
  end subroutine SolveConcNewtRaphMM
  
end module FatesKineticsMod
