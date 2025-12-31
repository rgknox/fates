module FatesKineticsMod
  
  implicit none

  
contains

  subroutine SolveMultipleCompetitors(s0, v_max, k_m, dt, vol, ncomp, s_final, total_mass, mass_indiv)

    ! -------------------------------------------------------------------------
    ! This module represents competitive uptake of a nutrient substrate
    ! for an arbitrary number of competitors. It is assumed that uptake
    ! is based on Michaelis-Menten kinetics. We perform an approximation
    ! of this system by creating a single set of effective uptake
    ! parameters vmax and k_m.  We integrate the system to find the total
    ! change in nutrient concentration, and then partition the flux
    ! to the competitors based on their uptake proportions at the
    ! beginning of the time-step.
    ! --------------------------------------------------------------------------
    
    integer, intent(in) :: ncomp               ! Number of competitors
    real(r8), intent(in) :: s0                 ! Initial nutrient concentration (g/m3 h2o)
    real(r8), intent(in) :: v_max(ncomp)       ! Max uptake rates (g/m3 h2o/t)
    real(r8), intent(in) :: k_m(ncomp)         ! Half-saturation constants (g/m3 h2o)
    real(r8), intent(in) :: dt                 ! Time step (s)
    real(r8), intent(in) :: vol                ! Water volume (m3)
    
    real(r8), intent(out) :: s_final           ! Final concentration (g/m3 h2o)
    real(r8), intent(out) :: total_mass        ! Total grams moved   (g)
    real(r8), intent(out) :: mass_indiv(ncomp) ! Grams moved per competitor
    
    real(r8) :: v_max_eff  ! Effective vmax of all competitors combined (sum)
    real(r8) :: k_m_eff    ! Effective K_m of all competitors (harmonic mean)
    real(r8) :: initial_flux_total   ! total substrate flux at t0
    real(r8) :: initial_flux(ncomp)  ! competitor substrate flux at t0
    integer  :: i                    ! competitor loop index
    
    ! 1. Calculate Effective Parameters
    v_max_eff = sum(v_max)
    k_m_eff = v_max_eff / sum(v_max / k_m)
    
    ! 2. Solve for S_final using Newton's method (Analytic solution to MM)
    ! The implicit equation: S0 - St + Km*ln(S0/St) = Vmax*dt
    s_final = SolveConc(s0, v_max_eff, k_m_eff, dt)
    
    ! 3. Calculate Mass Fluxes
    total_mass = (s0 - s_final) * vol
    
    ! 4. Partition based on initial proportions
    initial_flux = v_max * s0 / (k_m + s0)
    initial_flux_total = sum(initial_flux)
    
    do i = 1, ncomp
       mass_indiv(i) = (initial_flux(i) / initial_flux_total) * total_mass
    end do
    
  end subroutine SolveMultipleCompetitors

  ! ================================================================
  
  function SolveConc(s0, vmax, km, dt) result(st)

    ! Internal Newton-Raphson solver to find the root of the MM integration
    real(8), intent(in) :: s0, vmax, km, dt
    real(8) :: st, f, df, targ
    integer :: iter

    targ = vmax * dt
    st = s0  ! Initial guess
    
    do iter = 1, 20
       ! f(s) = (s0 - s) + km * ln(s0/s) - vmax*dt
       ! We want f(s) = 0
       f = (s0 - st) + km * log(s0/st) - targ
       ! df/ds = -1 - km/s
       df = -1.0d0 - (km / st)
       
       st = st - f/df
       
       if (abs(f) < 1.0d-10) exit
    end do
  end function SolveConc

end module FatesKineticsMod
