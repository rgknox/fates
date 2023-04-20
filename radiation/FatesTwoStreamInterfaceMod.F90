Module FatesTwoStreamInterfaceMod

  ! This module holds routines that are specific to connecting FATES with
  ! the two-stream radiation module. These routines are used to
  ! describe the scattering elements from cohort and patch data, and are
  ! used to decompose the scattering elements to return values
  ! at the cohort, or patch-pft scale.
  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesConstantsMod     , only : ifalse
  use FatesConstantsMod     , only : itrue
  use FatesConstantsMod     , only : nearzero
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use FatesGlobals          , only : fates_log
  use FatesGlobals          , only : endrun => fates_endrun
  use shr_infnan_mod        , only : nan => shr_infnan_nan, assignment(=)
  use FatesInterfaceTypesMod, only : numpft,hlm_numSWb
  use FatesRadiationMemMod  , only : ivis, inir
  use TwoStreamMLPEMod      , only : air_ft, AllocateRadParams, rad_params
  use EDTypesMod            , only : ed_patch_type, ed_cohort_type
  use EDTypesMod            , only : nclmax
  use TwoStreamMLPEMod      , only : twostream_type
  use EDPftvarcon           , only : EDPftvarcon_inst

  implicit none
  
  logical, parameter :: debug  = .false. ! local debug flag
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  
  contains

    
    subroutine FatesConstructElements(patch)


      type(ed_patch_type) :: patch

      type(ed_cohort_type), pointer :: cohort
      integer :: n_col(nclmax) ! Number of parallel column elements per layer
      integer :: ican,ft,icol
      type(twostream_type), pointer :: twostr
      real(r8), parameter :: canopy_open_frac = 0.0_r8

      integer :: maxcol
      real(r8) :: canopy_frac
      
      ! These parameters are not used yet
      !real(r8) :: max_vai_diff_per_elem ! The maximum vai difference in any element
      !                                  ! between the least and most vai of constituting
      !                                  ! cohorts.  THe objective is to reduce this.
      !integer, parameter  :: max_el_per_layer = 10
      !real(r8), parameter :: init_max_vai_diff_per_elem = 0.2_r8
      !type(ed_cohort_type), pointer :: elem_co_ptrs(ncl*max_el_per_layer,100)
      

      associate(twostr => patch%twostr)
      

      ! Identify how many elements we need, and possibly consolidate
      ! cohorts into elements where they are very similar (LAI and PFT)
      ! -------------------------------------------------------------------------------------------
      
      !max_vai_diff_per_elem = init_max_vai_diff_per_elem
      !iterate_count_do: do while(iterate_element_count)then
         
         ! Identify how many elements we need
         n_col(1:nclmax) = 0
         cohort => patch%tallest
         do while (associated(cohort))
            
            ft = cohort%pft
            ican = cohort%canopy_layer

            n_col(ican) = n_col(ican) + 1
            
            cohort => cohort%shorter
         enddo

      !   iterate_element_count = .false.
      !end do iterate_count_do

      ! Determine if we need an air element in each layer
      ! -------------------------------------------------------------------------------------------

      do ican = 1,patch%ncl_p-1
         if(canopy_open_frac>nearzero)then
            n_col(ican) = n_col(ican) + 1
         end if
      end do
      n_col(patch%ncl_p) = n_col(patch%ncl_p) + 1

      maxcol = -9
      do ican = 1,patch%ncl_p
         if (n_col(ican)>maxcol) maxcol=n_col(ican)
      end do

         
      ! Handle memory
      ! If the two-stream object is not large enough
      ! or if it is way larger than what is needed
      ! re-allocate the object
      ! -------------------------------------------------------------------------------------------
      
      if(.not.associated(twostr%scel)) then

         call twostr%AllocInitTwoStream((/ivis,inir/),patch%ncl_p,maxcol+2)
         
      else
         if(ubound(twostr%scel,2) <  maxcol .or. &
            ubound(twostr%scel,2) > (maxcol+4) .or. &
            ubound(twostr%scel,1) < patch%ncl_p ) then

            call twostr%DeallocTwoStream()
            
            ! Add a little more space than necessary so
            ! we don't have to keep allocating/deallocating
            call twostr%AllocInitTwoStream((/ivis,inir/),patch%ncl_p,maxcol+2)
            
         end if
      end if


      ! Fill the elements with their basic data and
      ! reference the cohort to the elements
      ! -------------------------------------------------------------------------------------------
      
      n_col(1:nclmax) = 0
      cohort => patch%tallest
      do while (associated(cohort))
         
         ft = cohort%pft
         ican = cohort%canopy_layer
         
         ! Every cohort gets its own element right now
         n_col(ican) = n_col(ican)+1

         if(ican.ne.patch%ncl_p) then
            canopy_frac = 1._r8 - canopy_open_frac
         else
            canopy_frac = 1._r8
         end if
         
         twostr%scel(ican,n_col(ican))%pft = ft
         twostr%scel(ican,n_col(ican))%area = canopy_frac*cohort%c_area/patch%total_canopy_area
         twostr%scel(ican,n_col(ican))%lai  = cohort%treelai
         twostr%scel(ican,n_col(ican))%sai  = cohort%treesai

         ! Cohort needs to know which column its in
         cohort%twostr_col = n_col(ican)
         
         cohort => cohort%shorter
      enddo


      ! Add the air (open) elements
      do ican = 1,patch%ncl_p

         canopy_frac = 0._r8
         do icol=1,n_col(ican)
            canopy_frac = canopy_frac + twostr%scel(ican,icol)%area
         end do
         
         if( (1._r8-canopy_frac)>nearzero ) then
            
            n_col(ican) = n_col(ican) + 1
            
            twostr%scel(ican,n_col(ican))%pft  = air_ft
            twostr%scel(ican,n_col(ican))%area = 1._r8-canopy_frac
            twostr%scel(ican,n_col(ican))%lai  = 0._r8
            twostr%scel(ican,n_col(ican))%sai  = 0._r8
         end if
            
      end do

      twostr%n_col(1:patch%ncl_p) = n_col(1:patch%ncl_p)

      if(debug) then
         do ican = 1,patch%ncl_p

            canopy_frac = 0._r8
            do icol=1,n_col(ican)
               canopy_frac = canopy_frac + twostr%scel(ican,icol)%area
            end do
            if( abs(1-canopy_frac) > nearzero ) then
               write(fates_log(),*) 'Element areas do not add up to 1: ',ican
               write(fates_log(),*) 'canopy_frac = ',canopy_frac
               do icol=1,n_col(ican)
                  write(fates_log(),*) 'area: ',twostr%scel(ican,icol)%area
               end do
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
         end do
      end if

      
      ! Set up some non-element parameters
      ! -------------------------------------------------------------------------------------------

      twostr%n_lyr = patch%ncl_p   ! Number of layers

      call twostr%GetNSCel()       ! Total number of elements
      
      twostr%force_prep = .true.   ! This signals that two-stream scattering coefficients
                                   ! that are dependent on geometry need to be updated

      
    end associate
    return
  end subroutine FatesConstructElements
  
  ! =============================================================================================
  
  subroutine FatesGetCohortAbsRad(patch,cohort,ib,vaitop,vaibot,rd_abs_leaf,rb_abs_leaf,r_abs_stem)

    ! This subroutine retrieves the absorbed radiation on
    ! leaves and stems, as well as the leaf sunlit fraction
    ! over a specified interval of VAI (vegetation area index)
    ! VAI is exposed leaf + stem area index

    type(ed_patch_type)  :: patch
    type(ed_cohort_type) :: cohort
    integer,intent(in)   :: ib
    real(r8),intent(in)  :: vaitop
    real(r8),intent(in)  :: vaibot
    real(r8),intent(out) :: rb_abs_leaf
    real(r8),intent(out) :: rd_abs_leaf
    real(r8),intent(out) :: r_abs_stem

    real(r8) :: vai_top_el
    real(r8) :: vai_bot_el
    real(r8) :: rd_abs_leaf_el
    real(r8) :: rb_abs_leaf_el
    real(r8) :: r_abs_stem_el

    associate(scelp => patch%twostr%scel(cohort%canopy_layer,cohort%twostr_col), &
              sccop => patch%twostr%band(ib)%scco(cohort%canopy_layer,cohort%twostr_col) )

      ! Convert the vai coordinate from the cohort to the element
      vai_top_el = vaitop * (scelp%lai+scelp%sai)/(cohort%treelai+cohort%treesai)
      vai_bot_el = vaibot * (scelp%lai+scelp%sai)/(cohort%treelai+cohort%treesai)

      ! Return the absorbed radiation for the element over that band
      call GetAbsRad(scelp,sccop,ib,vai_top_el,vai_bot_el,rd_abs_leaf_el,rb_abs_leaf_el,r_abs_stem_el)

      ! Scale the absorbed rad back to the cohort        
      rd_abs_leaf = rd_abs_leaf_el * cohort%treelai/scelp%lai
      rb_abs_leaf = rb_abs_leaf_el * cohort%treelai/scelp%lai
      r_abs_stem  = r_abs_stem_el * cohort%treesai/scelp%sai

    end associate
  end subroutine FatesGetCohortAbsRad

  ! =============================================================================================

  subroutine TransferRadParams()

    integer :: ft,ib  ! loop indices

    call AllocateRadParams(numpft,hlm_numSWb)

    do ft = 1,numpft
       do ib = 1,hlm_numSWb

          rad_params%rhol(ib,ft) = EDPftvarcon_inst%rhol(ib,ft)
          rad_params%rhos(ib,ft) = EDPftvarcon_inst%rhos(ib,ft)
          rad_params%taul(ib,ft) = EDPftvarcon_inst%taul(ib,ft)
          rad_params%taus(ib,ft) = EDPftvarcon_inst%taus(ib,ft)

       end do
       rad_params%xl(ft) = EDPftvarcon_inst%xl(ft)
       rad_params%clumping_index(ft) = EDPftvarcon_inst%clumping_index(ft)
    end do

  end subroutine TransferRadParams


end Module FatesTwoStreamInterfaceMod
