Module NormanRadiationMod

  use CanopyRadiationTypesMod

  use FatesConstantsMod, only: fates_unset_r8
  use FatesConstantsMod, only: r8 => fates_r8
  use FatesConstantsMod, only: pi_const
  
  subroutine PatchNormanRadiation (nrad, &
       numcl,         &
       numpft,        &
       fcansnow,      &   !
       szen_angle,    &   ! 
       elai_profile,  &   ! (ic,ft,iv)
       esai_profile,  &   ! (ic,ft,iv)
       carea_profile, &   ! (ic,ft,iv)
       albd_parb_out, &   ! (ft,ib)
       albi_parb_out, &   ! (ft,ib)
       fabd_parb_out, &   ! (ft,ib)
       fabi_parb_out, &   ! (ft,ib)
       ftdd_parb_out, &   ! (ft,ib)
       ftid_parb_out, &   ! (ft,ib)
       ftii_parb_out)     ! (ft,ib)
    
    ! -----------------------------------------------------------------------------------
    !
    ! This routine performs the Norman Radiation scattering for each patch.
    !
    ! -----------------------------------------------------------------------------------
    
    ! -----------------------------------------------------------------------------------
    ! !ARGUMENTS:
    ! -----------------------------------------------------------------------------------

    !type(ed_patch_type), intent(inout), target :: currentPatch

    integer,  intent(in)    :: numcl                ! number of canopy layers
    integer,  intent(in)    :: numpft               ! number of pfts
    integer,  intent(in)    :: nrad(maxpft,maxpft)  ! Number of vegetation layers with
                                                    ! any scattering media canopy layer x pft
    type(rad_scratch)       :: sc                     ! Scratch arrays
    real(r8), intent(in)    :: fcansnow               ! Fraction of canopy with snow on it (0-1)
    real(r8), intent(in)    :: szen_angle             ! solar zenith angle (radians)

    real(r8), intent(in)    :: elai_profile(numcl,numpft,maxlevleaf)
    real(r8), intent(in)    :: esai_profile(numcl,numpft,maxlevleaf)
    real(r8), intent(in)    :: carea_profile(numcl,numpft,maxlevleaf)
    
    real(r8), intent(inout) :: albd_parb_out(num_swb)
    real(r8), intent(inout) :: albi_parb_out(num_swb)
    real(r8), intent(inout) :: fabd_parb_out(num_swb)
    real(r8), intent(inout) :: fabi_parb_out(num_swb)
    real(r8), intent(inout) :: ftdd_parb_out(num_swb)
    real(r8), intent(inout) :: ftid_parb_out(num_swb)
    real(r8), intent(inout) :: ftii_parb_out(num_swb)

    
    ! Locals
    ! -----------------------------------------------------------------------------------

    integer  :: ic         ! canopy layer index
    integer  :: iv         ! vegetation layer index
    integer  :: ft         ! pft index (can't do "if", for obvious reasons)
    integer  :: ib         ! band index
    integer  :: iter                                          ! Iteration index
    integer  :: irep                                          ! Flag to exit iteration loop


    ! Scratch arrays
    real(r8) :: rho_layer(nclmax,maxpft,maxlevleaf,num_swb)  ! Weighted verage reflectance of layer
    real(r8) :: tau_layer(nclmax,maxpft,maxlevleaf,num_swb)  ! Weighted average transmittance of layer
    real(r8) :: tr_dir(nclmax,maxpft,maxlevleaf+1)         ! Transmittance of direct beam radiation through a single layer
    real(r8) :: tr_dif(nclmax,maxpft,maxlevleaf+1)         ! Transmittance of diffuse radiation through a single layer
    real(r8) :: f_abs(nclmax,maxpft,maxlevleaf,num_swb)      ! Fraction absorbed
    real(r8) :: f_abs_leaf(nclmax,maxpft,maxlevleaf,num_swb) ! Fraction absorbed by leaves
    real(r8) :: lai_change(nclmax,maxpft,nlevleaf-1)         ! Crown area change between layers

    
    ! For scattering elements
    real(r8) :: frac_lai                                ! Fraction of lai in each layer
    real(r8) :: frac_sai                                ! Fraction of sai in each layer
    real(r8) :: laisum                                  ! cumulative lai+sai for canopy layer (at middle of layer)

    
    ! for direct beam extinction coefficient
    real(r8) :: sb
    real(r8) :: cosz             ! 0.001 <= coszen <= 1.000
    real(r8) :: chil
    real(r8) :: gdir
    real(r8) :: k_dir(maxpft)                              ! Direct beam extinction coefficient
    real(r8) :: phi1b(maxpft)                                 ! Radiation transmitted to the soil surface.
    real(r8) :: phi2b(maxpft)
    real(r8) :: angle

    
    real(r8) :: error                                         ! Error check
    real(r8) :: down_rad, up_rad                              ! Iterative solution do Dif_dn and Dif_up

    

   

    real(r8) :: weighted_dir_tr(nclmax)
    real(r8) :: weighted_fsun(nclmax)
    real(r8) :: weighted_dif_ratio(nclmax,maxSWb)
    real(r8) :: weighted_dif_down(nclmax)
    real(r8) :: weighted_dif_up(nclmax)
    real(r8) :: refl_dif(nclmax,maxpft,nlevleaf,maxSWb)  ! Term for diffuse radiation reflected by laye
    real(r8) :: tran_dif(nclmax,maxpft,nlevleaf,maxSWb)  ! Term for diffuse radiation transmitted by layer
    real(r8) :: dif_ratio(nclmax,maxpft,nlevleaf,maxSWb) ! Ratio of upward to forward diffuse fluxes
    real(r8) :: Dif_dn(nclmax,maxpft,nlevleaf)           ! Forward diffuse flux onto canopy layer J (W/m**2 ground area)
    real(r8) :: Dif_up(nclmax,maxpft,nlevleaf)           ! Upward diffuse flux above canopy layer J (W/m**2 ground area)
    

    


    
    real(r8) :: Abs_dir_z(maxpft,nlevleaf)
    real(r8) :: Abs_dif_z(maxpft,nlevleaf)
    real(r8) :: abs_rad(maxSWb)                               !radiation absorbed by soil
    real(r8) :: tr_soili                                      ! Radiation transmitted to the soil surface.
    real(r8) :: tr_soild                                      ! Radiation transmitted to the soil surface.
    

    real(r8),parameter :: tolerance = 0.000000001_r8

    real(r8),parameter :: zenith_tolerance = 1.e-6_r8

    integer, parameter :: max_diag_nlevleaf = 4
    integer, parameter :: diag_nlevleaf = min(nlevleaf,max_diag_nlevleaf)  ! for diagnostics, write a small number of leaf layers

    real(r8) :: denom
    real(r8) :: lai_reduction(nclmax)

    integer  :: fp,iv,s      ! array indices
    


    real(r8), parameter :: forc_dir(n_rad_stream_types) = (/ 1.0_r8, 0.0_r8 /)   ! These are binary switches used
    real(r8), parameter :: forc_dif(n_rad_stream_types) = (/ 0.0_r8, 1.0_r8 /)   ! to turn off and on radiation streams


    ! Initialize local arrays
    
    weighted_dir_tr(:)   = 0._r8
    weighted_dif_down(:) = 0._r8
    weighted_dif_up(:)   = 0._r8
    
    Dif_up(:,:,:)        = 0._r8
    Dif_dn(:,:,:)        = 0._r8
    refl_dif(:,:,:,:)    = 0.0_r8
    tran_dif(:,:,:,:)    = 0.0_r8
    dif_ratio(:,:,:,:)   = 0.0_r8
    
    
    ! Initialize the ouput arrays
    ! ---------------------------------------------------------------------------------
    albd_parb_out(1:num_swb) = 0.0_r8
    albi_parb_out(1:num_swb) = 0.0_r8
    fabd_parb_out(1:num_swb) = 0.0_r8
    fabi_parb_out(1:num_swb) = 0.0_r8
    ftdd_parb_out(1:num_swb) = 1.0_r8
    ftid_parb_out(1:num_swb) = 1.0_r8
    ftii_parb_out(1:num_swb) = 1.0_r8
    

    
    if(debug)then
       if (sum(carea_profile(1,:,1))<0.999_r8)then
          write(fates_log(),*) 'canopy not full',carea_profile(1,:,1)
          stop
       endif
       if (sum(carea_profile(1,:,1))>1.0001_r8)then
          write(fates_log(),*) 'canopy too full',carea_profile(1,:,1)
          stop
       endif
    end if

    
    ! ==========================================================================
    ! PART I  -  GET LEAF LAYER REFLECTANCE, TRANSMISSION AND ABSORPTION
    ! ==========================================================================

    if(debug)then
       rho_layer(:,:,:,:)  = fates_unset_r8
       tau_layer(:,:,:,:)  = fates_unset_r8
       f_abs(:,:,:,:)      = fates_unset_r8
       f_abs_leaf(:,:,:,:) = fates_unset_r8
       
       tr_dif(:,:,:)     = fates_unset_r8
    end if
    
    do ic = 1,numcl
       do ft = 1,numpft
          do  iv = 1,nrad(ic,ft)

             ! veg layer index below current
             ivb = iv + 1
             
             ! area_change shows how much the crown area changes between layers.
             ! If there is more crown area in one layer than the layer below
             ! some of the radiation should move through that next layer unobstructed
             ! This change is on boundaries BETWEEN layers. Index 1 is between veg layer 1 and 2
             
             carea_change(ic,ft,iv) = 0._r8   ! default assumption, no change
             if (iv<nrad(ic,ft)) then
                if (( carea_profile(ic,ft,iv+1) >  0.0_r8 ) .and. ( carea_profile(ic,ft,iv+1)  <  carea_profile(ic,ft,iv) ))then
                   carea_change(ic,ft,iv) = carea_profile(ic,ft,iv)-carea_profile(ic,ft,iv+1)
                endif
             end if

                
             if(elai_profile(ic,ft,iv)+ esai_profile(ic,ft,iv).gt.0.0_r8) then
                frac_lai = elai_profile(ic,ft,iv)/&
                     (elai_profile(ic,ft,iv)+ esai_profile(ic,ft,iv))
             else
                frac_lai = 1.0_r8
             endif

             frac_sai = 1.0_r8 - frac_lai
                  
             ! layer level reflectance qualities
             do ib = 1,num_swb 

                rho_layer(ic,ft,iv,ib)=frac_lai*rad_params%rhol(ft,ib)+frac_sai*rad_params%rhos(ft,ib)
                tau_layer(ic,ft,iv,ib)=frac_lai*rad_params%taul(ft,ib)+frac_sai*rad_params%taus(ft,ib)
                
                ! adjust reflectance and transmittance for canopy snow
                rho_layer(ic,ft,iv,ib)=rho_layer(ic,ft,iv,ib)*(1.0_r8-fcansnow) &
                     + rad_params%rho_snow(ib) * fcansnow
                tau_layer(ic,ft,iv,ib)=tau_layer(ic,ft,iv,ib)*(1.0_r8-fcansnow) &
                     + rad_params%tau_snow(ib) * fcansnow
                
                ! fraction of incoming light absorbed by leaves or stems.
                f_abs(ic,ft,iv,ib) = 1.0_r8 - tau_layer(ic,ft,iv,ib) - rho_layer(ic,ft,iv,ib)
                
                ! the fraction of the vegetation absorbed light which is absorbed by leaves
                f_abs_leaf(ic,ft,iv,ib) = (1.0_r8- fcansnow) * frac_lai* &
                     (1.0_r8 - rad_params%rhol(ft,ib) - rad_params%taul(ft,ib))/f_abs(ic,ft,iv,ib)
                
             end do !ib

          end do !iv
       end do !ft
    end do !ic
    
    ! ==============================================================================
    ! PART II.  Direct beam extinction coefficient, k_dir. PFT specific. (HFREQ)
    !===============================================================================
    
    ! The zenith angle is measured as the deviation from straight above
    ! cos(0)    = 1, which implies direct sun
    ! cos(pi/2) or cos(-pi/2) = 0, which implies no direct sun
    
    ! Make sure the zenith angle is not perfectly at +/ 90 degrees, will cause math issues
    cosz = max(0.001_r8, cosz_in )
    
    do ft = 1,numpft
       sb = max(0.001_r8,(90._r8 - (acos(cosz)*180._r8/pi_const)) * (pi_const / 180._r8))
       chil = rad_params%xl(ft) !min(max(xl(ft), -0.4_r8), 0.6_r8 )
       if ( abs(chil) <= 0.01_r8) then
          chil  = 0.01_r8
       end if
       phi1b(ft) = 0.5_r8 - 0.633_r8*chil - 0.330_r8*chil*chil
       phi2b(ft) = 0.877_r8 * (1._r8 - 2._r8*phi1b(ft)) !0 = horiz leaves, 1 - vert leaves.
       gdir = phi1b(ft) + phi2b(ft) * sin(sb)
       !how much direct light penetrates a singleunit of lai?
       k_dir(ft) = clumping_index(ft) * gdir / sin(sb)
    end do !FT

    !------------------------------------------------------------------------------!
    ! Diffuse transmittance, tr_dif, do each layer with thickness elai_z.
    ! Estimated do nine sky angles in increments of 10 degrees
    ! PFT specific...  
    !------------------------------------------------------------------------------!
    
    do ic = 1,numcl !start at the top canopy layer (1 is the top layer.)
       do ft = 1,numpft
          tr_dif(ic,ft,:) = 0._r8
          do iv = 1,nrad(ic,ft)
             do j = 1,9
                angle = (5._r8 + real(j - 1,r8) * 10._r8) * pi_const / 180._r8
                gdir = phi1b(ft) + phi2b(ft) * sin(angle)
                tr_dif(ic,ft,iv) = tr_dif(ic,ft,iv) + exp(-rad_params%clumping_index(ft) * &
                     gdir / sin(angle) * &
                     (elai_profile(ic,ft,iv)+esai_profile(ic,ft,iv))) * &
                     sin(angle)*cos(angle)
             end do
             tr_dif(ic,ft,iv) = tr_dif(ic,ft,iv) * 2._r8 * (10._r8 * pi_const / 180._r8)
          end do
       end do
    end do
    
    ! ==============================================================================
    ! PART III.  
    ! ==============================================================================

    ! Get direct beam occlusion and sun/shade fractions
    call NormanDirTransFSun(nrad, &   ! IN
         numcl,                   &   ! IN
         numpft,                  &   ! IN
         elai_profile,            &   ! IN
         esai_profile,            &   ! IN
         carea_profile,           &   ! IN
         tr_dir,                  &   ! OUT
         f_sun,                   &   ! OUT
         weighted_dir_tr,         &   ! OUT
         weighted_fsun)               ! OUT

       
    ! do this for diffuse and direct beam
    do radtype = 1, num_rad_stream

         do ic = numcl,1, -1 !start at the bottom and work up.
            do ft = 1,numpft

                  !==============================================================================!
                  ! Iterative solution do scattering
                  !==============================================================================!

                  do ib = 1,num_swb !vis, nir
                     !------------------------------------------------------------------------------!
                     ! Leaf scattering coefficient and terms do diffuse radiation reflected
                     ! and transmitted by a layer
                     !------------------------------------------------------------------------------!

                     do iv = 1,nrad(ic,ft)
                        !How much diffuse light is intercepted and then reflected?
                        refl_dif(ic,ft,iv,ib) = (1._r8 - tr_dif(ic,ft,iv)) * rho_layer(ic,ft,iv,ib)
                        !How much diffuse light in this layer is transmitted?
                        tran_dif(ic,ft,iv,ib) = (1._r8 - tr_dif(ic,ft,iv)) * &
                             tau_layer(ic,ft,iv,ib) + tr_dif(ic,ft,iv)
                     end do

                     !------------------------------------------------------------------------------!
                     ! Ratio of upward to forward diffuse fluxes, dif_ratio
                     !------------------------------------------------------------------------------!
                     ! Soil diffuse reflectance (ratio of down to up radiation).
                     iv = nrad(ic,ft) + 1
                     if (ic  == numcl)then !nearest the soil
                        dif_ratio(ic,ft,iv,ib) = currentPatch%gnd_alb_dif(ib)  !bc_in(s)%albgr_dif_rb(ib)
                     else
                        dif_ratio(ic,ft,iv,ib) = weighted_dif_ratio(ic+1,ib)
                     end if
                     ! Canopy layers, working upwardfrom soil with dif_ratio(iv+1) known
                     ! FIX(RF,032414) ray tracing eqution - need to find derivation of this...
                     ! for each unit going down, there are x units going up.
                     do iv = nrad(ic,ft),1, -1
                        dif_ratio(ic,ft,iv,ib) = dif_ratio(ic,ft,iv+1,ib) * &
                             tran_dif(ic,ft,iv,ib)*tran_dif(ic,ft,iv,ib) / &
                             (1._r8 - dif_ratio(ic,ft,iv+1,ib) * refl_dif(ic,ft,iv,ib)) &
                             + refl_dif(ic,ft,iv,ib)
                        dif_ratio(ic,ft,iv,ib) = dif_ratio(ic,ft,iv,ib) * &
                             carea_profile(ic,ft,iv)/carea_profile(ic,ft,1)
                        dif_ratio(ic,ft,iv,ib) = dif_ratio(ic,ft,iv,ib) + dif_ratio(ic,ft,iv+1,ib) * &
                             (carea_profile(ic,ft,1)-carea_profile(ic,ft,iv))/carea_profile(ic,ft,1)
                     end do
                     weighted_dif_ratio(ic,ib) = weighted_dif_ratio(ic,ib) + &
                          dif_ratio(ic,ft,1,ib) * carea_profile(ic,ft,1)
                     !instance where the first layer carea_profile is used a proxy for the whole column. FTWA
                  end do!num_swb

                  
            end do!ft
         end do!ic

         ! Zero out the radiation error for the current patch before conducting the conservation check
         currentPatch%radiation_error = 0.0_r8

         do ib = 1,num_swb
            Dif_dn(:,:,:) = 0.00_r8
            Dif_up(:,:,:) = 0.00_r8
            do ic = 1, numcl !work down from the top of the canopy.
               weighted_dif_down(ic) = 0._r8
               do ft = 1, numpft


                     !------------------------------------------------------------------------------!
                     ! First estimates do downward and upward diffuse flux
                     !
                     ! Dif_dn =  forward diffuse flux onto layer J
                     ! Dif_up =  Upward diffuse flux above layer J
                     !
                     ! Solved here without direct beam radiation and using dif_ratio = Dif_up / Dif_dn
                     !------------------------------------------------------------------------------!
                     ! downward diffuse flux onto the top surface of the canopy

                     if (ic == 1)then
                        Dif_dn(ic,ft,1) = forc_dif(radtype)
                     else
                        Dif_dn(ic,ft,1) = weighted_dif_down(ic-1)
                     end if
                     ! forward diffuse flux within the canopy and at soil, working forward through canopy
                     do iv = 1,nrad(ic,ft)
                        denom = refl_dif(ic,ft,iv,ib) *  dif_ratio(ic,ft,iv,ib)
                        denom = 1._r8 - denom
                        Dif_dn(ic,ft,iv+1) = Dif_dn(ic,ft,iv) * tran_dif(ic,ft,iv,ib) / &
                             denom *carea_profile(ic,ft,iv)/carea_profile(ic,ft,1)
                        if (iv > 1)then
                           if (carea_change(ic,ft,iv-1) > 0.0_r8)then
                              !here we are thinking about whether the layer above had an laichange,
                              !but calculating the flux onto the layer below.
                              Dif_dn(ic,ft,iv+1) = Dif_dn(ic,ft,iv+1)+ Dif_dn(ic,ft,iv)* &
                                   carea_change(ic,ft,iv-1)/carea_profile(ic,ft,1)
                              Dif_dn(ic,ft,iv+1) = Dif_dn(ic,ft,iv+1)+ Dif_dn(ic,ft,iv-1)* &
                                   (carea_profile(ic,ft,1)-carea_profile(ic,ft,iv-1)/carea_profile(ic,ft,1))
                           else
                              Dif_dn(ic,ft,iv+1) = Dif_dn(ic,ft,iv+1) + Dif_dn(ic,ft,iv) * &
                                   (carea_profile(ic,ft,1)-carea_profile(ic,ft,iv))/carea_profile(ic,ft,1)
                           endif
                        else
                           Dif_dn(ic,ft,iv+1)    = Dif_dn(ic,ft,iv+1) + Dif_dn(ic,ft,iv) * &
                                (carea_profile(ic,ft,1)-carea_profile(ic,ft,iv))/carea_profile(ic,ft,1)
                        endif
                     end do

                     weighted_dif_down(ic) = weighted_dif_down(ic) + Dif_dn(ic,ft,nrad(ic,ft)+1) * &
                          carea_profile(ic,ft,1)

                     !instance where the first layer carea_profile is used a proxy for the whole column. FTWA

               end do !ft
               if (ic == currentPatch%NCic_p.and.currentPatch%NCic_p > 1)then !is the the (incomplete) understorey?
                  !Add on the radiation going through the canopy gaps.
                  weighted_dif_down(ic) = weighted_dif_down(ic) + weighted_dif_down(ic-1)*(1.0-sum(carea_profile(ic,:,1)))
                  !instance where the first layer carea_profile is used a proxy for the whole column. FTWA
               endif
            end do !ic

            do ic = numcl,1 ,-1 !work up from the bottom.
               weighted_dif_up(ic) = 0._r8
               do ft = 1, numpft

                     !Bounce diffuse radiation off soil surface.
                     iv = nrad(ic,ft) + 1
                     if (ic==numcl)then !is this the bottom layer ?
                        Dif_up(ic,ft,iv) = currentPatch%gnd_alb_dif(ib) * Dif_dn(ic,ft,iv)
                     else
                        Dif_up(ic,ft,iv) = weighted_dif_up(ic+1)
                     end if
                     ! Upward diffuse flux within the canopy and above the canopy, working upward through canopy

                     do iv = nrad(ic,ft), 1, -1
                        if (carea_change(ic,ft,iv) > 0.0_r8)then
                           Dif_up(ic,ft,iv) = dif_ratio(ic,ft,iv,ib) * Dif_dn(ic,ft,iv) * &
                                carea_profile(ic,ft,iv) / carea_profile(ic,ft,1)
                           Dif_up(ic,ft,iv) = Dif_up(ic,ft,iv) + Dif_up(ic,ft,iv+1) * &
                                tran_dif(ic,ft,iv,ib) * carea_change(ic,ft,iv)/carea_profile(ic,ft,1)
                           Dif_up(ic,ft,iv) = Dif_up(ic,ft,iv) + Dif_up(ic,ft,iv+1) * &
                                (carea_profile(ic,ft,1)-carea_profile(ic,ft,iv))/carea_profile(ic,ft,1)
                           !nb is this the right constuction?
                           ! the radiation that hits the empty space is not reflected.
                        else
                           Dif_up(ic,ft,iv) = dif_ratio(ic,ft,iv,ib) * Dif_dn(ic,ft,iv) * carea_profile(ic,ft,iv)
                           Dif_up(ic,ft,iv) = Dif_up(ic,ft,iv) + Dif_up(ic,ft,iv+1) * (1.0_r8-carea_profile(ic,ft,iv))
                        endif
                     end do

                     weighted_dif_up(ic) = weighted_dif_up(ic) + Dif_up(ic,ft,1) * carea_profile(ic,ft,1)
                     !instance where the first layer carea_profile is used a proxy for the whole column. FTWA

               end do !ft
               if (ic == numcl.and.numcl > 1)then !is this the (incomplete) understorey?
                  !Add on the radiation coming up through the canopy gaps.
                  !diffuse to diffuse
                  weighted_dif_up(ic) = weighted_dif_up(ic) +(1.0_r8-sum(carea_profile(ic,1:numpft,1))) * &
                       weighted_dif_down(ic-1) * currentPatch%gnd_alb_dif(ib)
                  !direct to diffuse
                  weighted_dif_up(ic) = weighted_dif_up(ic) + forc_dir(radtype) * &
                       weighted_dir_tr(ic-1) * (1.0_r8-sum(carea_profile(ic,1:numpft,1))) * currentPatch%gnd_alb_dir(ib)
               endif
            end do !ic

            !------------------------------------------------------------------------------!
            ! 3. Iterative calculation of forward and upward diffuse fluxes, iNCL_puding
            !    scattered direct beam
            !------------------------------------------------------------------------------!

            ! Flag to exit iteration loop: 0 = exit and 1 = iterate
            irep = 1
            ! Iteration loop
            iter = 0
            do while(irep ==1 .and. iter<50)

               iter = iter + 1
               irep = 0
               do ic = 1,numcl !working from the top down
                  weighted_dif_down(ic) = 0._r8
                  do ft =1,numpft

                        ! forward diffuse flux within the canopy and at soil, working forward through canopy
                        ! with Dif_up -from previous iteration-. Dif_dn(1) is the forward diffuse flux onto the canopy.
                        ! Note: down = forward flux onto next layer
                        if (ic == 1)then !is this the top layer?
                           Dif_dn(ic,ft,1) = forc_dif(radtype)
                        else
                           Dif_dn(ic,ft,1) = weighted_dif_down(ic-1)
                        end if
                        down_rad = 0._r8

                        do iv = 1, nrad(ic,ft)
                           ! down rad'n is the sum of the down and upwards reflected diffuse fluxes...
                           down_rad = Dif_dn(ic,ft,iv) * tran_dif(ic,ft,iv,ib) + &
                                Dif_up(ic,ft,iv+1) * refl_dif(ic,ft,iv,ib)

                           !... plus the direct beam intercepted and intransmitted by this layer.
                           down_rad = down_rad + forc_dir(radtype) * tr_dir(ic,ft,iv) * (1.00_r8 - &
                                exp(-k_dir(ft) *  (elai_profile(ic,ft,iv)+ &
                                esai_profile(ic,ft,iv))  )) * tau_layer(ic,ft,iv,ib)


                           !... plus the direct beam intercepted and intransmitted by this layer.
                           ! modified to spread it out over the whole of incomplete layers.

                           down_rad = down_rad *(carea_profile(ic,ft,iv)/carea_profile(ic,ft,1))

                           if (iv > 1)then
                              if (carea_change(ic,ft,iv-1) > 0.0_r8)then
                                 down_rad = down_rad + Dif_dn(ic,ft,iv)   * carea_change(ic,ft,iv-1)/carea_profile(ic,ft,1)
                                 down_rad = down_rad + Dif_dn(ic,ft,iv-1) * (carea_profile(ic,ft,1)-carea_profile(ic,ft,iv-1))/ &
                                      carea_profile(ic,ft,1)
                              else
                                 down_rad = down_rad + Dif_dn(ic,ft,iv)   * (carea_profile(ic,ft,1)-carea_profile(ic,ft,iv))/ &
                                      carea_profile(ic,ft,1)
                              endif
                           else
                              down_rad = down_rad + Dif_dn(ic,ft,iv)   * (carea_profile(ic,ft,1)-carea_profile(ic,ft,iv))/ &
                                   carea_profile(ic,ft,1)
                           endif

                           !this is just Dif down, plus refl up, plus dir intercepted and turned into dif... ,
                           if (abs(down_rad - Dif_dn(ic,ft,iv+1)) > tolerance)then
                              irep = 1
                           end if
                           Dif_dn(ic,ft,iv+1) = down_rad

                        end do !iv

                        weighted_dif_down(ic) = weighted_dif_down(ic) + Dif_dn(ic,ft,nrad(ic,ft)+1) * &
                             carea_profile(ic,ft,1)

                     
                  end do!ft
                  if (ic == numcl.and.numcl > 1)then !is this the (incomplete) understorey?
                     weighted_dif_down(ic) = weighted_dif_down(ic) + weighted_dif_down(ic-1) * &
                          (1.0_r8-sum(carea_profile(ic,1:numpft,1)))
                  end if
               end do ! do ic loop

               do ic = 1, numcl ! working from the top down.
                  weighted_dif_up(ic) = 0._r8
                  do ft =1,numpft

                        ! Upward diffuse flux at soil or from lower canopy (forward diffuse and unscattered direct beam)
                        iv = nrad(ic,ft) + 1
                        if (ic==currentPatch%NCic_p)then  !In the bottom canopy layer, reflect off the soil
                           Dif_up(ic,ft,iv) = Dif_dn(ic,ft,iv) * currentPatch%gnd_alb_dif(ib) + &
                                forc_dir(radtype) * tr_dir(ic,ft,iv) * currentPatch%gnd_alb_dir(ib)
                        else      !In the other canopy layers, reflect off the underlying vegetation.
                           Dif_up(ic,ft,iv) =  weighted_dif_up(ic+1)
                        end if

                        ! Upward diffuse flux within and above the canopy, working upward through canopy
                        ! with Dif_dn from previous interation.  Note: up = upward flux above current layer
                        do iv = nrad(ic,ft),1,-1
                           !this is radiation up, by layer transmittance, by

                           !reflection of the lower layer,
                           up_rad = Dif_dn(ic,ft,iv) * refl_dif(ic,ft,iv,ib)
                           up_rad = up_rad + forc_dir(radtype) * tr_dir(ic,ft,iv) * (1.00_r8 - exp(-k_dir(ft) * &
                                (elai_profile(ic,ft,iv)+esai_profile(ic,ft,iv))))* &
                                rho_layer(ic,ft,iv,ib)
                           up_rad = up_rad + Dif_up(ic,ft,iv+1) * tran_dif(ic,ft,iv,ib)
                           up_rad = up_rad * carea_profile(ic,ft,iv)/carea_profile(ic,ft,1)
                           up_rad = up_rad + Dif_up(ic,ft,iv+1) *(carea_profile(ic,ft,1)-carea_profile(ic,ft,iv))/carea_profile(ic,ft,1)
                           ! THE LOWER LAYER FLUX IS HOMOGENIZED, SO WE DON"T CONSIDER THE CAREA_CHANGE HERE...

                           if (abs(up_rad - Dif_up(ic,ft,iv)) > tolerance) then !are we close to the tolerance level?
                              irep = 1
                           end if
                           Dif_up(ic,ft,iv) = up_rad

                        end do  !iv
                        weighted_dif_up(ic) = weighted_dif_up(ic) + Dif_up(ic,ft,1) * carea_profile(ic,ft,1)

                  end do!ft

                  if (ic == numcl.and.numcl > 1)then  !is this the (incomplete) understorey?
                     !Add on the radiation coming up through the canopy gaps.
                     weighted_dif_up(ic) = weighted_dif_up(ic) +(1.0_r8-sum(carea_profile(ic,1:numpft,1))) * &
                          weighted_dif_down(ic-1) * currentPatch%gnd_alb_dif(ib)
                     weighted_dif_up(ic) = weighted_dif_up(ic) + forc_dir(radtype) * &
                          weighted_dir_tr(ic-1) * (1.0_r8-sum(carea_profile(ic,1:numpft,1)))*currentPatch%gnd_alb_dir(ib)
                  end if
               end do!ic
            end do ! do while over iter

            abs_rad(ib) = 0._r8
            tr_soili = 0._r8
            tr_soild = 0._r8

            do ic = 1, numcl !working from the top down.
               abs_dir_z(:,:) = 0._r8
               abs_dif_z(:,:) = 0._r8
               do ft =1,numpft

                     !==============================================================================!
                     ! Compute absorbed flux densities
                     !==============================================================================!

                     ! Absorbed direct beam and diffuse do leaf layers
                     do iv = 1, nrad(ic,ft)
                        Abs_dir_z(ft,iv) = carea_profile(ic,ft,iv)* forc_dir(radtype) * tr_dir(ic,ft,iv) * &
                             (1.00_r8 - exp(-k_dir(ft) *  (elai_profile(ic,ft,iv)+ &
                             esai_profile(ic,ft,iv)) )) * f_abs(ic,ft,iv,ib)
                        Abs_dif_z(ft,iv) = carea_profile(ic,ft,iv)* ((Dif_dn(ic,ft,iv) + &
                             Dif_up(ic,ft,iv+1)) * (1.00_r8 - tr_dif(ic,ft,iv)) * f_abs(ic,ft,iv,ib))
                     end do

                     ! Absorbed direct beam and diffuse do soil
                     if (ic == numcl)then
                        iv = nrad(ic,ft) + 1
                        Abs_dif_z(ft,iv) = carea_profile(ic,ft,1)*Dif_dn(ic,ft,iv) * (1.0_r8 - currentPatch%gnd_alb_dif(ib) )
                        Abs_dir_z(ft,iv) = carea_profile(ic,ft,1)*forc_dir(radtype) * &
                             tr_dir(ic,ft,iv) * (1.0_r8 - currentPatch%gnd_alb_dir(ib)  )
                        tr_soild = tr_soild + carea_profile(ic,ft,1)*forc_dir(radtype) * tr_dir(ic,ft,iv)
                        tr_soili = tr_soili + carea_profile(ic,ft,1)*Dif_dn(ic,ft,iv)
                     end if

                     ! Absorbed radiation, shaded and sunlit portions of leaf layers
                     !here we get one unit of diffuse radiation... how much of
                     !it is absorbed?
                     if (ib == ivis) then ! only set the absorbed PAR for the visible light band.
                        do iv = 1, nrad(ic,ft)
                           if (radtype==idirect) then
                              if ( debug ) then
                                 write(fates_log(),*) 'EDsurfAlb 730 ',Abs_dif_z(ft,iv),currentPatch%f_sun(ic,ft,iv)
                                 write(fates_log(),*) 'EDsurfAlb 731 ', currentPatch%fabd_sha_z(ic,ft,iv), &
                                      currentPatch%fabd_sun_z(ic,ft,iv)
                              endif
                              currentPatch%fabd_sha_z(ic,ft,iv) = Abs_dif_z(ft,iv) * &
                                   (1._r8 - currentPatch%f_sun(ic,ft,iv))*f_abs_leaf(ic,ft,iv,ib)
                              currentPatch%fabd_sun_z(ic,ft,iv) =( Abs_dif_z(ft,iv) * &
                                   currentPatch%f_sun(ic,ft,iv) + &
                                   Abs_dir_z(ft,iv))*f_abs_leaf(ic,ft,iv,ib)
                           else
                              currentPatch%fabi_sha_z(ic,ft,iv) = Abs_dif_z(ft,iv) * &
                                   (1._r8 - currentPatch%f_sun(ic,ft,iv))*f_abs_leaf(ic,ft,iv,ib)
                              currentPatch%fabi_sun_z(ic,ft,iv) = Abs_dif_z(ft,iv) * &
                                   currentPatch%f_sun(ic,ft,iv)*f_abs_leaf(ic,ft,iv,ib)
                           endif
                           if ( debug ) then
                              write(fates_log(),*) 'EDsurfAlb 740 ', currentPatch%fabd_sha_z(ic,ft,iv), &
                                   currentPatch%fabd_sun_z(ic,ft,iv)
                           endif
                        end do
                     endif ! ib


                     !==============================================================================!
                     ! Sum fluxes
                     !==============================================================================!
                     ! Solar radiation absorbed by ground
                     iv = nrad(ic,ft) + 1
                     if (ic==numcl)then
                        abs_rad(ib) = abs_rad(ib) +  (Abs_dir_z(ft,iv) + Abs_dif_z(ft,iv))
                     end if
                     ! Solar radiation absorbed by vegetation and sunlit/shaded leaves
                     do iv = 1,nrad(ic,ft)
                        if (radtype == idirect)then
                           currentPatch%fabd(ib) = currentPatch%fabd(ib) + &
                                Abs_dir_z(ft,iv)+Abs_dif_z(ft,iv)
                           ! bc_out(s)%fabd_parb_out(ib) = currentPatch%fabd(ib)
                        else
                           currentPatch%fabi(ib) = currentPatch%fabi(ib) + Abs_dif_z(ft,iv)
                           ! bc_out(s)%fabi_parb_out(ib) = currentPatch%fabi(ib)
                        endif
                     end do

                     ! Albefor
                     if (ic==1)then !top canopy layer.
                        if (radtype == idirect)then
                           albd_parb_out(ib) = albd_parb_out(ib) + &
                                Dif_up(ic,ft,1) * carea_profile(ic,ft,1)
                        else
                           albi_parb_out(ib) = albi_parb_out(ib) + &
                                Dif_up(ic,ft,1) * carea_profile(ic,ft,1)
                        end if
                     end if

                     ! pass normalized PAR profiles for use in diagnostic averaging for history fields
                     if (ib == ivis) then ! only diagnose PAR profiles for the visible band
                        do iv = 1, nrad(ic,ft)
                           currentPatch%nrmlzd_parprof_pft_dir_z(radtype,ic,ft,iv) = &
                                forc_dir(radtype) * tr_dir(ic,ft,iv)
                           currentPatch%nrmlzd_parprof_pft_dif_z(radtype,ic,ft,iv) = &
                                Dif_dn(ic,ft,iv) + Dif_up(ic,ft,iv)
                           !
                           currentPatch%nrmlzd_parprof_dir_z(radtype,ic,iv) = &
                                currentPatch%nrmlzd_parprof_dir_z(radtype,ic,iv) + &
                                (forc_dir(radtype) * tr_dir(ic,ft,iv)) * &
                                (carea_profile(ic,ft,iv) / sum(carea_profile(ic,1:numpft,iv)))
                           currentPatch%nrmlzd_parprof_dif_z(radtype,ic,iv) = &
                                currentPatch%nrmlzd_parprof_dif_z(radtype,ic,iv) + &
                                (Dif_dn(ic,ft,iv) + Dif_up(ic,ft,iv)) * &
                                (carea_profile(ic,ft,iv) / sum(carea_profile(ic,1:numpft,iv)))
                        end do
                     end if ! ib = visible

               end do !ft
               if (radtype == idirect)then
                  fabd_parb_out(ib) = currentPatch%fabd(ib)
               else
                  fabi_parb_out(ib) = currentPatch%fabi(ib)
               endif


               !radiation absorbed from fluxes through unfilled part of lower canopy.
               if (numcl > 1.and.ic == numcl)then
                  abs_rad(ib) = abs_rad(ib) + weighted_dif_down(ic-1) * &
                       (1.0_r8-sum(carea_profile(ic,1:numpft,1)))*(1.0_r8-currentPatch%gnd_alb_dif(ib) )
                  abs_rad(ib) = abs_rad(ib) + forc_dir(radtype) * weighted_dir_tr(ic-1) * &
                       (1.0_r8-sum(carea_profile(ic,1:numpft,1)))*(1.0_r8-currentPatch%gnd_alb_dir(ib) )
                  tr_soili = tr_soili + weighted_dif_down(ic-1) * (1.0_r8-sum(carea_profile(ic,1:numpft,1)))
                  tr_soild = tr_soild + forc_dir(radtype) * weighted_dir_tr(ic-1) * (1.0_r8-sum(carea_profile(ic,1:numpft,1)))
               endif

               if (radtype == idirect)then
                  currentPatch%tr_soil_dir(ib) = tr_soild
                  currentPatch%tr_soil_dir_dif(ib) = tr_soili
                  currentPatch%sabs_dir(ib)     = abs_rad(ib)
                  ftdd_parb_out(ib)  = tr_soild
                  ftid_parb_out(ib) =  tr_soili
               else
                  currentPatch%tr_soil_dif(ib) = tr_soili
                  currentPatch%sabs_dif(ib)     = abs_rad(ib)
                  ftii_parb_out(ib) =  tr_soili
               end if

            end do!l


            !==============================================================================!
            ! Conservation check
            !==============================================================================!
            ! Total radiation balance: absorbed = incoming - outgoing

            if (radtype == idirect)then
               error = abs(currentPatch%sabs_dir(ib) - (currentPatch%tr_soil_dir(ib) * &
                    (1.0_r8-currentPatch%gnd_alb_dir(ib) ) + &
                    currentPatch%tr_soil_dir_dif(ib) * (1.0_r8-currentPatch%gnd_alb_dif(ib)     )))

               if(debug)then
                  if ( abs(error) > 0.0001)then
                     write(fates_log(),*)'dir ground absorption error',error,currentPatch%sabs_dir(ib), &
                          currentPatch%tr_soil_dir(ib)* &
                          (1.0_r8-currentPatch%gnd_alb_dir(ib)  ),numcl,ib,sum(carea_profile(1,1:numpft,1))
                     write(fates_log(),*) 'albedos',currentPatch%sabs_dir(ib) ,currentPatch%tr_soil_dir(ib), &
                          (1.0_r8-currentPatch%gnd_alb_dir(ib)  )
                     do ft =1,numpft
                        iv = nrad(1,ft) + 1
                        write(fates_log(),*) 'abs soil fluxes', Abs_dir_z(ft,iv),Abs_dif_z(ft,iv)
                     end do
                  end if
               end if

            else
               if (debug) then
                  if ( abs(currentPatch%sabs_dif(ib)-(currentPatch%tr_soil_dif(ib) * &
                       (1.0_r8-currentPatch%gnd_alb_dif(ib)  ))) > 0.0001_r8)then
                     write(fates_log(),*)'dif ground absorption error',currentPatch%sabs_dif(ib) , &
                          (currentPatch%tr_soil_dif(ib)* &
                          (1.0_r8-currentPatch%gnd_alb_dif(ib)  )),numcl,ib,sum(carea_profile(1,1:numpft,1))
                  endif
               end if
            endif

            if (radtype == idirect)then
               error = (forc_dir(radtype) + forc_dif(radtype)) - &
                    (fabd_parb_out(ib)  + albd_parb_out(ib) + currentPatch%sabs_dir(ib))
            else
               error = (forc_dir(radtype) + forc_dif(radtype)) - &
                    (fabi_parb_out(ib)  + albi_parb_out(ib) + currentPatch%sabs_dif(ib))
            endif

            ! ignore the current patch radiation error if the veg-covered fraction of the patch is really small
            if ( (currentPatch%total_canopy_area / currentPatch%area) .gt. tolerance ) then
               ! normalize rad error by the veg-covered fraction of the patch because that is
               ! the only part that this code applies to
               currentPatch%radiation_error = currentPatch%radiation_error + error &
                    * currentPatch%total_canopy_area / currentPatch%area
            endif

            lai_reduction(:) = 0.0_r8
            do ic = 1, numcl
               do ft =1,numpft
                     do iv = 1, nrad(ic,ft)
                        if (carea_change(ic,ft,iv) > 0.0_r8)then
                           lai_reduction(ic) = max(lai_reduction(ic),carea_change(ic,ft,iv))
                        endif
                     enddo
               enddo
            enddo

            if (radtype == idirect)then
               !here we are adding a within-ED radiation scheme tolerance, and then adding the diffrence onto the albedo
               !it is important that the lower boundary for this is ~1000 times smaller than the tolerance in surface albedo.
               if (abs(error)  >  1.e-9_r8 .and. abs(error) < 0.15_r8)then
                  albd_parb_out(ib) = albd_parb_out(ib) + error
                  !this terms adds the error back on to the albedo. While this is partly inexcusable, it is
                  ! in the medium term a solution that
                  ! prevents the model from crashing with small and occasional energy balances issues.
                  ! These are extremely difficult to debug, many have been solved already, leading
                  ! to the complexity of this code, but where the system generates occasional errors, we
                  ! will deal with them for now.
               end if

               if (abs(error)  >  0.15_r8)then
                  if(debug)then
                     write(fates_log(),*) 'Large Dir Radn consvn error',error ,ib
                     write(fates_log(),*) 'diags', albd_parb_out(ib), ftdd_parb_out(ib), &
                          ftid_parb_out(ib), fabd_parb_out(ib)
                     write(fates_log(),*) 'elai',currentpatch%elai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                     write(fates_log(),*) 'esai',currentpatch%esai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                     write(fates_log(),*) 'carea_profile',carea_profile(1,1:numpft,1:diag_nlevleaf)
                     write(fates_log(),*) 'cp',currentPatch%area, currentPatch%patchno
                     write(fates_log(),*) 'ground albedo diffuse (ib)', currentPatch%gnd_alb_dir(ib)
                  end if
                  albd_parb_out(ib) = albd_parb_out(ib) + error
               end if
            else

               if (abs(error)  >  1.e-9_r8 .and. abs(error) < 0.15_r8)then
                  albi_parb_out(ib) = albi_parb_out(ib) + error
               end if

               if (abs(error)  >  0.15_r8)then
                  if(debug)then
                     write(fates_log(),*)  'lg Dif Radn consvn error',error ,ib
                     write(fates_log(),*) 'diags', albi_parb_out(ib), ftii_parb_out(ib), &
                          fabi_parb_out(ib)
                     !write(fates_log(),*) 'carea_change',carea_change(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                     !write(fates_log(),*) 'elai',currentpatch%elai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                     !write(fates_log(),*) 'esai',currentpatch%esai_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                     !write(fates_log(),*) 'carea_profile',carea_profile(currentpatch%ncl_p,1:numpft,1:diag_nlevleaf)
                     write(fates_log(),*) 'cp',currentPatch%area, currentPatch%patchno
                     write(fates_log(),*) 'ground albedo diffuse (ib)', currentPatch%gnd_alb_dir(ib)
                     !write(fates_log(),*) 'rhol',rhol(1:numpft,:)
                     !write(fates_log(),*) 'ftw',sum(carea_profile(1,1:numpft,1)),carea_profile(1,1:numpft,1)
                     !write(fates_log(),*) 'CAP',currentPatch%canopy_area_profile(1,1:numpft,1)
                  end if
                  albi_parb_out(ib) = albi_parb_out(ib) + error
               end if

               if (radtype == idirect)then
                  error = (forc_dir(radtype) + forc_dif(radtype)) - &
                       (fabd_parb_out(ib)  + albd_parb_out(ib) + currentPatch%sabs_dir(ib))
               else
                  error = (forc_dir(radtype) + forc_dif(radtype)) - &
                       (fabi_parb_out(ib)  + albi_parb_out(ib) + currentPatch%sabs_dif(ib))
               endif

               if(debug) then
                  if (abs(error)  >  0.00000001_r8)then
                     write(fates_log(),*)  'there is still error after correction',error ,ib
                  end if
               end if

            end if
         end do !num_swb

      enddo ! rad-type


    end associate
    return
  end subroutine PatchNormanRadiation

  ! ================================================================================

  
  subroutine NormanDirTransFSun(nrad,              &   ! IN
                                numcl,             &   ! IN
                                numpft,            &   ! IN
                                elai_profile,      &   ! IN
                                esai_profile,      &   ! IN
                                carea_profile,     &   ! IN
                                carea_change,      &   ! IN
                                tr_dir,            &   ! OUT
                                f_sun,             &   ! OUT
                                weighted_dir_tr,   &   ! OUT
                                weighted_fsun)         ! OUT

    ! -------------------------------------------------------------------------------
    ! This routine returns, at different veg layers:
    ! 1) the direct beam free path transmission (ie not including forward scattering)
    ! 2) The sun/shade fractions
    ! -------------------------------------------------------------------------------

    integer,  intent(in)    :: numcl                ! number of canopy layers
    integer,  intent(in)    :: numpft               ! number of pfts
    integer,  intent(in)    :: nrad(maxpft,maxpft)  ! Number of vegetation layers with
                                                    ! any scattering media canopy layer x pft

    real(r8), intent(in)    :: elai_profile(numcl,numpft,maxlevleaf)
    real(r8), intent(in)    :: esai_profile(numcl,numpft,maxlevleaf)
    real(r8), intent(in)    :: carea_profile(numcl,numpft,maxlevleaf)
    real(r8), intent(in)    :: carea_change(numcl,numpft,maxlevleaf-1)

    real(r8), intent(out) :: tr_dir(nclmax,maxpft,maxlevleaf+1)
    real(r8), intent(out) :: f_sun(nclmax,maxpft,maxlevleaf)

    real(r8), intent(out) :: weighted_dir_tr(nclmax)
    real(r8), intent(out) :: weighted_fsun(nclmax)

    real(r8) :: vaisum       ! Vertically integrated vegetation area index

    ! Methods of calculating direct interception
    
    integer, parameter :: dir_int_base = 1
    integer, parameter :: dir_int_simp = 2
    integer, parameter :: dir_int = dir_int_base

    
    integer :: ic,ft,iv
    
    ! Based on the leaf+stem profiles, calculate transmittance of direct radation
    ! and fsun

    if(debug) then
       tr_dir(:,:,:) = fates_unset_r8
       f_sun(:,:,:)  = fates_unset_r8
    end if

    
    select case(dir_int)
    case(dir_int_simp) then

       do ic = 1,numcl !start at the top canopy layer (1 is the top layer.)
          
          ! etai_above:  mean total (stem+leaf) area index in the
          ! canopy layer above current (per m2 of total canopy area)
          if(ic==1)then
             vai_above = 0._r8
          else
             do ft = 1,numpft
                do iv = 1,nrad(ic,ft)
                   vai_above=vai_above + &
                        (elai_profile(ic,ft,iv)+esai_profile(ic,ft,iv))*carea_profile(ic,ft,iv)
                end do
             end if
          end if

          
          do ft = 1,numpft
             
             !------------------------------------------------------------------------------!
             ! Direct beam transmittance, tr_dir, uses cumulative VAI above layer iv
             ! (ie iv+1) to give unscattered direct beam onto layer iv. do each PFT section.
             ! This is just an  decay curve based on k_dir. (leaf & sun angle)
             !------------------------------------------------------------------------------!
             if (ic==1)then
                tr_dir(ic,ft,1) = 1._r8
             else
                tr_dir(ic,ft,1) = weighted_dir_tr(ic-1)
             endif
             
             vaisum = vai_above

             sum_tr_dir: do iv = 1,nrad(ic,ft)

                ! VAI of the current layer, in m2 per m2 of the PFT's footprint
                mean_pft_vai = (elai_profile(ic,ft,iv)+esai_profile(ic,ft,iv)) * &
                                carea_profile(ic,ft,iv)/carea_profile(ic,ft,1)
                
                ! Integrated VAI to mid-point of the bin
                vaisum_cent = vaisum + 0.5*mean_pft_vai
                
                ! Integrated LAI to the bottom of the bin
                vaisum = vaisum + mean_pft_vai

                ! Fraction of the beam radiation attenuated 
                frac_atten = exp(-k_dir(ft) * vaisum) - exp(-k_dir(ft) * (vaisum+mean_pft_vai))

                tr_dir(ic,ft,iv+1) = tr_dir(ic,ft,iv)-frac_atten
                
                f_sun(ic,ft,iv) = 

                
                int_frac_intercept = 
                
                tr_dir(ic,ft,iv+1) = exp(-k_dir(ft) * vaisum)

                if (ic == 1)then !top canopy layer
                   f_sun(ic,ft,iv) = exp(-k_dir(ft) * laisum)* &
                        (carea_profile(ic,ft,iv)/carea_profile(ic,ft,1))
                else
                   f_sun(ic,ft,iv) = weighted_fsun(ic-1)* exp(-k_dir(ft) * laisum)* &
                        (carea_profile(ic,ft,iv)/carea_profile(ic,ft,1))
                endif
                
             end do sum_tr_dir
          end do
       end do
    case(dir_int_base) then

       do ic = 1,numcl !start at the top canopy layer (1 is the top layer.)
          
          weighted_dir_tr(ic)              = 0._r8
          weighted_fsun(ic)                = 0._r8
          
          ! Each canopy layer (canopy, understorey) has multiple 'parallel' pft's
          
          do ft = 1,numpft
             
             !------------------------------------------------------------------------------!
             ! Direct beam transmittance, tr_dir, uses cumulative LAI above layer iv
             ! (ie iv+1) to give unscattered direct beam onto layer iv. do each PFT section.
             ! This is just an  decay curve based on k_dir. (leaf & sun angle)
             !------------------------------------------------------------------------------!
             if (ic==1)then
                tr_dir(ic,ft,1) = 1._r8
             else
                tr_dir(ic,ft,1) = weighted_dir_tr(ic-1)
             endif
                   
             laisum = 0._r8
             !total direct beam getting to the bottom of the top canopy.
             sum_tr_dir: do iv = 1,nrad(ic,ft)

                ! This is light coming straight through the canopy.
                ! leaf+stem area per m2 of the PFT's ground area
                laisum = laisum + elai_profile(ic,ft,iv)+esai_profile(ic,ft,iv)
                if (ic==1)then
                   tr_dir(ic,ft,iv+1) = exp(-k_dir(ft) * laisum) * &
                        (carea_profile(ic,ft,iv)/carea_profile(ic,ft,1))
                else
                   tr_dir(ic,ft,iv+1) = weighted_dir_tr(ic-1)*exp(-k_dir(ft) * laisum)* &
                        (carea_profile(ic,ft,iv)/carea_profile(ic,ft,1))
                endif

             
             
             if (iv == 1)then
                !this is the top layer.
                tr_dir(ic,ft,iv+1) = tr_dir(ic,ft,iv+1) + tr_dir(ic,ft,iv) * &
                     ((carea_profile(ic,ft,1)-carea_profile(ic,ft,iv))/carea_profile(ic,ft,1))
             else
                ! the carea_change(iv) affects the light incident on layer iv+2 not iv+1
                ! light coming from the layer above (iv-1) goes through iv and onto iv+1.
                if (carea_change(ic,ft,iv-1) > 0.0_r8)then
                   tr_dir(ic,ft,iv+1) = tr_dir(ic,ft,iv+1) + tr_dir(ic,ft,iv)* &
                        carea_change(ic,ft,iv-1) / carea_profile(ic,ft,1)
                   tr_dir(ic,ft,iv+1) = tr_dir(ic,ft,iv+1) + tr_dir(ic,ft,iv-1)* &
                        (carea_profile(ic,ft,1)-carea_profile(ic,ft,iv-1))/carea_profile(ic,ft,1)
                else
                   !account fot the light that comes striaght down from unfilled layers above.
                   tr_dir(ic,ft,iv+1) = tr_dir(ic,ft,iv+1) + tr_dir(ic,ft,iv) * &
                        ((carea_profile(ic,ft,1)-carea_profile(ic,ft,iv))/carea_profile(ic,ft,1))
                endif
             endif
             
          end do sum_tr_dir
             
          !add up all the weighted contributions from the different PFT columns.
          weighted_dir_tr(ic) = weighted_dir_tr(ic) + tr_dir(ic,ft,nrad(ic,ft)+1)*carea_profile(ic,ft,1)
          
          !------------------------------------------------------------------------------!
          ! Sunlit and shaded fraction of leaf layer
          !------------------------------------------------------------------------------!
          
          !laisum = 0._r8
          do_fsun: do iv = 1,nrad(ic,ft)

             ! Cumulative leaf area. Original code uses cumulative lai do layer.
             ! Now use cumulative lai at center of layer.
             ! Same as tr_dir calcualtions, but in the middle of the layer? FIX(RF,032414)-WHY?
             
             if (iv  ==  1) then
                laisum = 0.5_r8 * (elai_profile(ic,ft,iv)+esai_profile(ic,ft,iv))
             else
                laisum = laisum + elai_profile(ic,ft,iv)+esai_profile(ic,ft,iv)
             end if
             
             if (ic == 1)then !top canopy layer
                f_sun(ic,ft,iv) = exp(-k_dir(ft) * laisum)* &
                     (carea_profile(ic,ft,iv)/carea_profile(ic,ft,1))
             else
                f_sun(ic,ft,iv) = weighted_fsun(ic-1)* exp(-k_dir(ft) * laisum)* &
                     (carea_profile(ic,ft,iv)/carea_profile(ic,ft,1))
             endif
             
             if ( iv > 1 ) then  ! becasue we are looking at this layer (not the next)
                ! we only ever add fluxes if iv>1
                if (carea_change(ic,ft,iv-1) > 0.0_r8)then
                   f_sun(ic,ft,iv) = f_sun(ic,ft,iv) + &
                        f_sun(ic,ft,iv) * &
                        carea_change(ic,ft,iv-1)/carea_profile(ic,ft,1)
                   f_sun(ic,ft,iv) = f_sun(ic,ft,iv) + &
                        f_sun(ic,ft,iv-1) * &
                        (carea_profile(ic,ft,1)-carea_profile(ic,ft,iv-1))/carea_profile(ic,ft,1)
                else
                   f_sun(ic,ft,iv) = f_sun(ic,ft,iv) + &
                        f_sun(ic,ft,iv-1) * &
                        (carea_profile(ic,ft,1)-carea_profile(ic,ft,iv))/carea_profile(ic,ft,1)
                endif
             endif
             
          end do do_fsun
          
          weighted_fsun(ic) = weighted_fsun(ic) + f_sun(ic,ft,nrad(ic,ft))* &
               carea_profile(ic,ft,1)
          
          ! instance where the first layer carea_profile is used a proxy for the whole column. FTWA
          ! this is possibly a source of slight error. If we use the carea_profile at the top of the PFT column,
          ! then we willl underestimate fsun, but if we use carea_profile at the bottom of the column, we will
          ! underestimate it. Really, we should be tracking the release of direct light from the column as it tapers
          ! towards the ground. Is that necessary to get energy closure? It would be quite hard...
          
          
          
       end do!pft loop
       
    end do !ic
    
    return
  end subroutine NormanDirTransFSun
  
end Module NormanRadiationMod
