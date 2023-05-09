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
  use FatesRadiationMemMod  , only : rho_snow,tau_snow
  use TwoStreamMLPEMod      , only : air_ft, AllocateRadParams, rad_params
  use EDTypesMod            , only : ed_patch_type, ed_cohort_type, ed_site_type
  use EDTypesMod            , only : nclmax
  use TwoStreamMLPEMod      , only : twostream_type
  use TwoStreamMLPEMod      , only : ParamPrep
  use TwoStreamMLPEMod      , only : AllocateRadParams
  use EDPftvarcon           , only : EDPftvarcon_inst
  use FatesRadiationMemMod  , only : rad_solver,twostr_solver
  use FatesAllometryMod     , only : VegAreaLayer
  
  
  implicit none

  logical, parameter :: debug  = .false. ! local debug flag
  character(len=*), parameter, private :: sourcefile = &
       __FILE__


  public :: FatesConstructRadElements
  public :: FatesGetCohortAbsRad
  public :: FatesPatchFsun
  
contains


  subroutine FatesConstructRadElements(site,fcansno_pa,coszen_pa)

    type(ed_site_type)  :: site
    type(ed_patch_type),pointer :: patch
    real(r8)                    :: fcansno_pa(:)
    real(r8)                    :: coszen_pa(:)
    
    type(ed_cohort_type), pointer :: cohort
    integer :: n_col(nclmax) ! Number of parallel column elements per layer
    integer :: ican,ft,icol
    type(twostream_type), pointer :: twostr


    ! DO NOT MAKE CANOPY_OPEN_FRAC >0 UNTIL LAI COMPRESSION
    ! HAS BEEN THOUGHT THROUGH. WE CANT JUST DECREASE THE
    ! AREA WITHOUT CONSERVING TOTAL LEAF AND STEM AREA
    real(r8), parameter :: canopy_open_frac = 0.00_r8

    integer :: maxcol
    real(r8) :: canopy_frac
    integer  :: ifp
    ! Area indices for the cohort [m2 media / m2 crown footprint]
    real(r8) :: elai_cohort,tlai_cohort,esai_cohort,tsai_cohort
    real(r8) :: vai_top,vai_bot  ! veg area index at top and bottom of cohort (dummy vars)

    ! These parameters are not used yet
    !real(r8) :: max_vai_diff_per_elem ! The maximum vai difference in any element
    !                                  ! between the least and most vai of constituting
    !                                  ! cohorts.  THe objective is to reduce this.
    !integer, parameter  :: max_el_per_layer = 10
    !real(r8), parameter :: init_max_vai_diff_per_elem = 0.2_r8
    !type(ed_cohort_type), pointer :: elem_co_ptrs(ncl*max_el_per_layer,100)

    
    if(rad_solver.ne.twostr_solver)return

    ifp=0
    patch => site%oldest_patch
    do while (associated(patch))
       ifp=ifp+1
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

         ! If there is only one layer, then we don't
         ! need to add an air element to the only
         ! layer. This is because all non-veg
         ! area will be attributed to a ground patch
         ! But if there is more than one layer, then
         ! an air element is needed for all the non
         ! occupied space, even if the canopy_open_frac
         ! is zero.
            
         if(patch%ncl_p>1)then
            canopy_frac = 0._r8
            do icol=1,n_col(patch%ncl_p)
               canopy_frac = canopy_frac + twostr%scelg(patch%ncl_p,icol)%area
            end do
            if( (1._r8-canopy_frac)>nearzero ) then
               n_col(patch%ncl_p) = n_col(patch%ncl_p) + 1
            end if
         end if


         ! Handle memory
         ! If the two-stream object is not large enough
         ! or if it is way larger than what is needed
         ! re-allocate the object
         ! -------------------------------------------------------------------------------------------

         maxcol = 0
         do ican = 1,patch%ncl_p
            if (n_col(ican)>maxcol) maxcol=n_col(ican)
         end do

         if(.not.associated(twostr%scelg)) then

            call twostr%AllocInitTwoStream((/ivis,inir/),patch%ncl_p,maxcol+2)

         else

            print*,maxcol,ubound(twostr%scelg,2)
            
            if(ubound(twostr%scelg,2) <  maxcol .or. &
               ubound(twostr%scelg,2) > (maxcol+4) .or. &
               ubound(twostr%scelg,1) < patch%ncl_p ) then

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

            patch%canopy_mask(ican,ft) = 1
            
            ! Every cohort gets its own element right now
            n_col(ican) = n_col(ican)+1

            if(ican.ne.patch%ncl_p) then
               canopy_frac = 1._r8 - canopy_open_frac
            else
               canopy_frac = 1._r8
            end if

            ! If we pass layer index 0 to this routine
            ! it will return the total plant LAIs and SAIs
            call VegAreaLayer(cohort%treelai, &
                              cohort%treesai, &
                              cohort%hite,    &
                              0,                     &
                              cohort%nv,      &
                              cohort%pft,     &
                              site%snow_depth,       &
                              vai_top, vai_bot,      & 
                              elai_cohort,esai_cohort)
            
            twostr%scelg(ican,n_col(ican))%pft = ft
            twostr%scelg(ican,n_col(ican))%area = canopy_frac*cohort%c_area/patch%total_canopy_area
            twostr%scelg(ican,n_col(ican))%lai  = elai_cohort
            twostr%scelg(ican,n_col(ican))%sai  = esai_cohort

            ! Cohort needs to know which column its in
            cohort%twostr_col = n_col(ican)

            cohort => cohort%shorter
         enddo

         ! Add the air (open) elements

         if(patch%ncl_p>1)then
            canopy_frac = 0._r8
            ican = patch%ncl_p
            do icol=1,n_col(ican)
               canopy_frac = canopy_frac + twostr%scelg(ican,icol)%area
            end do
            if( (1._r8-canopy_frac)>nearzero ) then
               n_col(ican) = n_col(ican) + 1
               twostr%scelg(ican,n_col(ican))%pft  = air_ft
               twostr%scelg(ican,n_col(ican))%area = 1._r8-canopy_frac
               twostr%scelg(ican,n_col(ican))%lai  = 0._r8
               twostr%scelg(ican,n_col(ican))%sai  = 0._r8
            end if
         end if

         twostr%n_col(1:patch%ncl_p) = n_col(1:patch%ncl_p)

         if(debug) then
            do ican = 1,patch%ncl_p
               canopy_frac = 0._r8
               do icol=1,n_col(ican)
                  canopy_frac = canopy_frac + twostr%scelg(ican,icol)%area
               end do
               if( abs(1-canopy_frac) > nearzero ) then
                  write(fates_log(),*) 'Element areas do not add up to 1: ',ican
                  write(fates_log(),*) 'canopy_frac = ',canopy_frac
                  do icol=1,n_col(ican)
                     write(fates_log(),*) 'area: ',twostr%scelg(ican,icol)%area
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
         call twostr%CanopyPrep(fcansno_pa(ifp))
         call twostr%ZenithPrep(coszen_pa(ifp))
         
     end associate

      

       
       patch => patch%younger
    end do

    return
  end subroutine FatesConstructRadElements

  ! =============================================================================================

  subroutine FatesPatchFSun(patch,fsun,laisun,laisha)

    type(ed_patch_type) :: patch
    real(r8)            :: fsun    ! Patch average sunlit fraction
    real(r8)            :: laisun  ! Patch average LAI of leaves in sun
    real(r8)            :: laisha  ! Patch average LAI of leaves in shade

    integer :: ican, icol  ! Canopy vertical and horizontal element index

    ! Dummy variables
    real(r8)            :: Rb_abs,Rd_abs,Rd_abs_leaf,Rb_abs_leaf,R_abs_stem,R_abs_snow
    real(r8)            :: leaf_sun_frac  ! Element specific sunlit fraction of leaf

    laisun = 0._r8
    laisha = 0._r8

    associate(twostr => patch%twostr, &
              scelg => patch%twostr%scelg)
    
      do ican = 1,twostr%n_lyr
         do icol = 1,twostr%n_col(ican)
            
            call two_str%GetAbsRad(ican,icol,ivis,0._r8,scelg%lai+scelg%sai, &
                 Rb_abs,Rd_abs,Rd_abs_leaf,Rb_abs_leaf,R_abs_stem,R_abs_snow,leaf_sun_frac)

            laisun = laisun + scelg%area*scelg%lai*leaf_sun_frac
            laisha = laisha + scelg%area*scelg%lai*(1._r8-leaf_sun_frac)
            
         end do
      end do

      if(laisun+laisha>nearzero)then
         fsun = laisun / (laisun+laisha)
      else
         fsun = 0.5_r8  ! Nominal value, should not affect results if no leaves or light!
      end if
         
    end associate
    return
  end subroutine FatesPatchFSun

   ! =============================================================================================
  subroutine FatesGetCohortAbsRad(patch,cohort,ib,vaitop,vaibot,cohort_elai,cohort_esai,rd_abs_leaf,rb_abs_leaf,leaf_sun_frac )

    ! This subroutine retrieves the absorbed radiation on
    ! leaves and stems, as well as the leaf sunlit fraction
    ! over a specified interval of VAI (vegetation area index)
    ! VAI is exposed leaf + stem area index

    type(ed_patch_type)  :: patch
    type(ed_cohort_type) :: cohort
    integer,intent(in)   :: ib
    real(r8),intent(in)  :: vaitop
    real(r8),intent(in)  :: vaibot
    real(r8),intent(in)  :: cohort_elai
    real(r8),intent(in)  :: cohort_esai
    real(r8),intent(out) :: rb_abs_leaf
    real(r8),intent(out) :: rd_abs_leaf
    real(r8),intent(out) :: leaf_sun_frac

    real(r8) :: rb_abs,rd_abs
    real(r8) :: rd_abs_el,rb_abs_el
    real(r8) :: vai_top_el
    real(r8) :: vai_bot_el
    real(r8) :: rd_abs_leaf_el
    real(r8) :: rb_abs_leaf_el
    real(r8) :: r_abs_stem_el
    real(r8) :: r_abs_snow_el
    real(r8) :: diff_wt_leaf,diff_wt_elem
    real(r8) :: beam_wt_leaf,beam_wt_elem
    real(r8) :: evai_cvai  ! element VAI / cohort VAI
    
    associate(scelg => patch%twostr%scelg(cohort%canopy_layer,cohort%twostr_col), &
         scelb => patch%twostr%band(ib)%scelb(cohort%canopy_layer,cohort%twostr_col) )

      evai_cvai = (scelg%lai+scelg%sai)/(cohort_elai+cohort_esai)
      
      if(abs(evai_cvai-1._r8)>1.e-8_r8)then
         print*,"EVAI_CVAI: ",evai_cvai
         stop
      end if
      
      ! Convert the vai coordinate from the cohort to the element
      vai_top_el = vaitop * evai_cvai
      vai_bot_el = vaibot * evai_cvai

      ! Return the absorbed radiation for the element over that band
      call patch%twostr%GetAbsRad(cohort%canopy_layer,cohort%twostr_col,ib,vai_top_el,vai_bot_el, & 
           Rb_abs_el,Rd_abs_el,rd_abs_leaf_el,rb_abs_leaf_el,r_abs_stem_el,r_abs_snow_el,leaf_sun_frac)

      rd_abs = rd_abs_el / evai_cvai
      rb_abs = rb_abs_el / evai_cvai
      
      diff_wt_leaf = (1._r8-patch%twostr%frac_snow)*cohort_elai*(1._r8-rad_params%om_leaf(ib,cohort%pft))*rad_params%Kd_leaf(cohort%pft)
      diff_wt_elem = (cohort_elai+cohort_esai)*(1._r8-scelb%om)*scelg%Kd

      beam_wt_leaf = (1._r8-patch%twostr%frac_snow)*cohort_elai*(1._r8-rad_params%om_leaf(ib,cohort%pft))*scelg%Kb_leaf
      beam_wt_elem = (cohort_elai+cohort_esai)*(1._r8-scelb%om)*scelg%Kb

      print*,"---"
      print*,diff_wt_leaf,diff_wt_elem
      print*,cohort_elai,cohort_esai
      print*,rad_params%om_leaf(ib,cohort%pft),rad_params%Kd_leaf(cohort%pft)
      print*,scelb%om,scelg%Kd

      
      rd_abs_leaf = rd_abs * min(1.0_r8,diff_wt_leaf / diff_wt_elem)
      rb_abs_leaf = rb_abs * min(1.0_r8,beam_wt_leaf / beam_wt_elem)

      
    end associate
  end subroutine FatesGetCohortAbsRad

  ! =============================================================================================

  subroutine TransferRadParams()

    integer :: ft,ib  ! loop indices

    call AllocateRadParams(numpft,hlm_numSWb)

    do ft = 1,numpft
       do ib = 1,hlm_numSWb

          rad_params%rhol(ib,ft) = EDPftvarcon_inst%rhol(ft,ib)
          rad_params%rhos(ib,ft) = EDPftvarcon_inst%rhos(ft,ib)
          rad_params%taul(ib,ft) = EDPftvarcon_inst%taul(ft,ib)
          rad_params%taus(ib,ft) = EDPftvarcon_inst%taus(ft,ib)

       end do
       rad_params%xl(ft) = EDPftvarcon_inst%xl(ft)
       rad_params%clumping_index(ft) = EDPftvarcon_inst%clumping_index(ft)
    end do

    call ParamPrep()

    return
  end subroutine TransferRadParams


end Module FatesTwoStreamInterfaceMod
