module FatesCNDDMod

  use FatesConstantsMod, only: r8 => fates_r8
  use FatesConstantsMod, only: nearzero
  use FatesGlobals,      only: fates_log
  use FatesCohortMod,    only: fates_cohort_type
  use FatesPatchMod,     only: fates_patch_type
  use PRTGenericMod,     only: leaf_organ
  use PRTGenericMod,     only: carbon12_element
  use PRTGenericMod,     only: element_pos
  use EDPftvarcon,       only: pft_param => EDPftvarcon_inst
  use FatesInterfaceTypesMod, only: numpft
  use FatesLitterMod,    only: litter_type
  use EDParamsMod,       only: dev_arbitrary

  implicit none
  private
  
  ! Turnover Scaling Function, max and min
  real(r8), parameter :: cndd_ts_min = 0.2_r8
  real(r8), parameter :: cndd_ts_max = 5.0_r8

  public :: CNDDPatch
  
  contains
  
  subroutine CNDDPatch(patch)

    type(fates_patch_type),intent(inout)  :: patch
    type(fates_cohort_type),pointer :: cohort
    type(litter_type),pointer       :: littc
    integer  :: p
    real(r8) :: pft_share(1:numpft)
    
    integer, parameter :: height_filter = 1
    integer, parameter :: layer_filter = 2
    integer, parameter :: no_filter = 3
    integer, parameter :: cndd_filter = no_filter
    real(r8),parameter :: height_crit = 10._r8

    ! The fraction of seed biomass that can be extrapolated
    ! to unrealized saplings below the minimum recruit size,
    ! but only considering the organs that serve as a 
    ! predation source for herbivores and parasitic fungi/bacteria.
    
    real(r8),parameter :: sapling_pred_frac = 0.2_r8
    
    ! First, lets loop through the PFTs, and the cohorts
    ! that comprise those PFTs and quantify the presence of the 
    ! PFT. This may based on leaf mass, leaf area, seed mass, or
    ! some combination of both. We will call this score the
    ! pft_share. We will go on to normalize this.

    littc => patch%litter(element_pos(carbon12_element))
    
    pft_share(:) = 0._r8
    do p = 1, numpft
       cohort => patch%tallest
       do while (associated(cohort))
          if(cohort%pft==p)then
             select case(cndd_filter)
             case(no_filter)
                pft_share(p) = pft_share(p)+cohort%prt%GetState(leaf_organ,carbon12_element)*cohort%n
                !pft_share(p) = pft_share(p)+cohort%tree_lai*cohort%c_area
             case(layer_filter)
                if(cohort%canopy_layer>1)then
                   pft_share(p) = pft_share(p)+cohort%prt%GetState(leaf_organ,carbon12_element)*cohort%n
                   !pft_share(p) = pft_share(p)+cohort%tree_lai*cohort%c_area
                end if
             case(height_filter)
                if(cohort%height<height_crit)then
                   pft_share(p) = pft_share(p)+cohort%prt%GetState(leaf_organ,carbon12_element)*cohort%n
                   !pft_share(p) = pft_share(p)+cohort%tree_lai*cohort%c_area
                end if
             end select
          end if
          cohort => cohort%shorter
       enddo

       ! Include un-realized saplings in the presence of the PFT
       ! The factor involved should consider both the fraction of relevent biomass
       ! 
       pft_share(p) = pft_share(p) + littc%seed(p)*patch%area*sapling_pred_frac
       
    end do
    
    ! Second, normalize the pft_share vector, and
    ! use those shares to create an objective function
    ! that will affect turnover
    
    if(sum(pft_share)>nearzero)then
       pft_share = pft_share / sum(pft_share)
       do p = 1, numpft
          if(pft_param%cndd_frac(p) > nearzero) then
             patch%cndd_ts(p) = CNDDRFunc1(pft_share(p),pft_param%cndd_frac(p))
             !patch%cndd_ts(p) = CNDDRFunc2(pft_share(p),pft_param%cndd_frac(p))
          else
             patch%cndd_ts(p) = 1.0_r8
          end if
       end do
       
    else
       patch%cndd_ts(1:numpft) = 1.0_r8
    end if

  end subroutine CNDDPatch

  ! ============================================================================
  
  pure function CNDDMM(cndd_ts) result(cndd_ts_protected)

    real(r8),intent(in) :: cndd_ts
    real(r8) :: cndd_ts_protected
    
    ! Protect the response function in a min/max
    cndd_ts_protected = min(max(cndd_ts,cndd_ts_min),cndd_ts_max)
    
  end function CNDDMM

  ! ============================================================================
  
  pure function CNDDRFunc1(cndd_share, cndd_target) result(cndd_ts)

    real(r8),intent(in)  :: cndd_share
    real(r8),intent(in)  :: cndd_target
    real(r8) :: cndd_ts
    real(r8), parameter :: alpha = 1._r8
    
    !cndd_ts = CNDDMM(exp(alpha*(cndd_share-cndd_target)))
    cndd_ts = CNDDMM(exp(dev_arbitrary*(cndd_share-cndd_target)))
    
  end function CNDDRFunc1
  
  ! ============================================================================
  
  pure function CNDDRFunc2(cndd_share, cndd_target) result(cndd_ts)

    real(r8),intent(in)  :: cndd_share
    real(r8),intent(in)  :: cndd_target
    real(r8) :: cndd_ts
    
    cndd_ts = CNDDMM(cndd_share/cndd_target)
    
  end function CNDDRFunc2

  
  
end module FatesCNDDMod
