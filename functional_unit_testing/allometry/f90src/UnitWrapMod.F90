
! =======================================================================================
!
! This file is an alternative to key files in the fates
! filesystem. Noteably, we replace fates_r8 and fates_in
! with types that work with "ctypes".  This is
! a key step in working with python
!
! We also wrap FatesGlobals to reduce the dependancy
! cascade that it pulls in from shr_log_mod.
!
! =======================================================================================

module shr_log_mod

   use iso_c_binding, only : c_char
   use iso_c_binding, only : c_int

   contains

   function shr_log_errMsg(source, line) result(ans)
      character(kind=c_char,len=*), intent(in) :: source
      integer(c_int), intent(in) :: line
      character(kind=c_char,len=128) :: ans

      ans = "source: " // trim(source) // " line: "
   end function shr_log_errMsg

end module shr_log_mod


module FatesGlobals

   contains

   integer function fates_log()
      fates_log = 11
   end function fates_log

   subroutine fates_endrun(msg)

      implicit none
      character(len=*), intent(in) :: msg    ! string to be printed

      stop

   end subroutine fates_endrun

end module FatesGlobals


module PRTParametersMod
  use iso_c_binding, only : r8 => c_double
  use iso_c_binding, only : i4 => c_int

  

  type,public ::  prt_param_type

  ! VARIABLE-DEFINITIONS-HERE

	 real(r8), pointer :: allom_zroot_k(:)
	 real(r8), pointer :: allom_zroot_min_z(:)
	 real(r8), pointer :: allom_zroot_min_dbh(:)
	 real(r8), pointer :: allom_zroot_max_z(:)
	 real(r8), pointer :: allom_zroot_max_dbh(:)
	 real(r8), pointer :: allom_agb4(:)
	 real(r8), pointer :: allom_agb3(:)
	 real(r8), pointer :: allom_agb2(:)
	 real(r8), pointer :: allom_agb1(:)
	 real(r8), pointer :: allom_d2ca_coefficient_min(:)
	 real(r8), pointer :: allom_d2ca_coefficient_max(:)
	 real(r8), pointer :: allom_blca_expnt_diff(:)
	 real(r8), pointer :: allom_d2bl3(:)
	 real(r8), pointer :: allom_d2bl2(:)
	 real(r8), pointer :: allom_d2bl1(:)
	 real(r8), pointer :: allom_d2h3(:)
	 real(r8), pointer :: allom_d2h2(:)
	 real(r8), pointer :: allom_d2h1(:)
	 real(r8), pointer :: allom_agb_frac(:)
	 real(r8), pointer :: allom_l2fr(:)
	 real(r8), pointer :: allom_la_per_sa_slp(:)
	 real(r8), pointer :: allom_la_per_sa_int(:)
	 real(r8), pointer :: allom_stmode(:)
	 real(r8), pointer :: allom_smode(:)
	 real(r8), pointer :: allom_cmode(:)
	 real(r8), pointer :: allom_amode(:)
	 real(r8), pointer :: allom_fmode(:)
	 real(r8), pointer :: allom_lmode(:)
	 real(r8), pointer :: allom_hmode(:)
	 real(r8), pointer :: allom_dbh_maxheight(:)
	 real(r8), pointer :: allom_sai_scaler(:)
	 real(r8), pointer :: slatop(:)
	 real(r8), pointer :: slamax(:)
	 real(r8), pointer :: crown(:)
	 real(r8), pointer :: woody(:)
	 real(r8), pointer :: wood_density(:)
	 real(r8), pointer :: c2b(:)
	 real(r8), pointer :: fnrt_prof_b(:)
	 real(r8), pointer :: fnrt_prof_a(:)
	 real(r8), pointer :: fnrt_prof_mode(:)
	 real(r8), pointer :: seed_alloc(:)
	 real(r8), pointer :: seed_alloc_mature(:)
	 real(r8), pointer :: dbh_repro_threshold(:)
	 real(r8), pointer :: leaf_stor_priority(:)
	 real(r8), pointer :: cushion(:)
	 integer, pointer :: organ_id(:)
	 real(r8), pointer :: phos_store_ratio(:)
	 real(r8), pointer :: nitr_store_ratio(:)
	 real(r8), pointer :: phos_stoich_p2(:,:)
	 real(r8), pointer :: phos_stoich_p1(:,:)
	 real(r8), pointer :: nitr_stoich_p2(:,:)
	 real(r8), pointer :: nitr_stoich_p1(:,:)
	 real(r8), pointer :: grperc(:)
	 real(r8), pointer :: turnover_phos_retrans(:,:)
	 real(r8), pointer :: turnover_nitr_retrans(:,:)
	 real(r8), pointer :: turnover_carb_retrans(:,:)
	 real(r8), pointer :: turnover_retrans_mode(:)
	 real(r8), pointer :: branch_long(:)
	 real(r8), pointer :: root_long(:)
	 real(r8), pointer :: leaf_long(:)
	 real(r8), pointer :: senleaf_long_fdrought(:)
	 integer, pointer :: evergreen(:)
	 integer, pointer :: season_decid(:)
	 integer, pointer :: stress_decid(:)
  end type prt_param_type

  type(prt_param_type),public :: prt_params

end module PRTParametersMod


module PRTParamsGeneric

   use iso_c_binding, only : r8 => c_double
   use iso_c_binding, only : i4 => c_int
   use iso_c_binding, only : c_char
   use PRTParametersMod,only: prt_params
   integer,parameter :: shr_kind_cs = 80                     ! short char
  type ptr_var1
     real(r8), dimension(:), pointer :: var_rp
     integer(i4), dimension(:), pointer :: var_ip
     character(len=shr_kind_cs) :: var_name
     integer :: vtype
  end type ptr_var1

  type ptr_var2
     real(r8), dimension(:,:), pointer :: var_rp
     integer(i4), dimension(:,:), pointer :: var_ip
     character(len=shr_kind_cs) :: var_name
     integer :: vtype
  end type ptr_var2


  type prt_params_ptr_type
     type(ptr_var1), allocatable :: var1d(:)
     type(ptr_var2), allocatable :: var2d(:)
  end type prt_params_ptr_type


  type(prt_params_ptr_type),  public :: prt_params_ptr


  integer :: numparm1d ! Number of different PFT parameters
  integer :: numparm2d
  integer :: numpft
  logical, parameter ::  debug = .false.

contains


   subroutine PRTParamsPySet(rval,ival,indx1,indx2,name)

      implicit none
      ! Arguments
      character(kind=c_char,len=*), intent(in) :: name
      real(r8),intent(in) :: rval
      integer(i4),intent(in) :: ival
	   integer(i4),intent(in) :: indx1
	   integer(i4),intent(in) :: indx2
      ! Locals
      logical :: npfound
      integer :: ip
      integer :: namelen

      namelen = len(trim(name))

      if(debug) print*,"F90: ARGS: ",trim(name),"INDX1: ",indx1,"INDX2: ",indx2," RVAL: ",rval," IVAL: ",ival
      ip=0
      npfound = .false.
      do ip=1,numparm1d
         
         if (trim(name) == trim(prt_params_ptr%var1d(ip)%var_name ) ) then
            if(debug) print*,"F90: Found ",trim(name)," in lookup table as index:",ip

            npfound = .true.
            if(prt_params_ptr%var1d(ip)%vtype == 1) then ! real
               prt_params_ptr%var1d(ip)%var_rp(indx1) = rval
            elseif(prt_params_ptr%var1d(ip)%vtype == 2) then ! integer
               prt_params_ptr%var1d(ip)%var_ip(indx1) = ival
            else
               print*,"F90: STRANGE TYPE"
               stop
            end if
         end if
      end do

      do ip=1,numparm2d
         if (trim(name) == trim(prt_params_ptr%var2d(ip)%var_name ) ) then
            if(debug) print*,"F90: Found ",trim(name)," in lookup table"
            npfound = .true.
            if(prt_params_ptr%var2d(ip)%vtype == 1) then ! real
               prt_params_ptr%var2d(ip)%var_rp(indx1,indx2) = rval
            elseif(prt_params_ptr%var2d(ip)%vtype == 2) then ! integer
               prt_params_ptr%var2d(ip)%var_ip(indx1,indx2) = ival
            else
               print*,"F90: STRANGE TYPE"
               stop
            end if
         end if
      end do



      if(.not.npfound)then
         print*,"F90: The parameter you loaded DNE: ",name(:)
         stop
      end if


	  ! Perform a check to see if the target array is being filled
      if (trim(name) == 'fates_allom_d2h1') then
         if (prt_params%allom_d2h1(indx1) == rval) then
            if(debug) print*,"F90: POINTER CHECK PASSES:",rval," = ",prt_params%allom_d2h1(indx1)
         else
            print*,"F90: POINTER CHECK FAILS:",rval," != ",prt_params%allom_d2h1(indx1)
            stop
         end if
      end if

      if (trim(name) == 'fates_wood_density' ) then
         if (prt_params%wood_density(indx1) == rval) then
            if(debug) print*,"F90: POINTER CHECK PASSES:",rval," = ",prt_params%wood_density(indx1)
         else
            print*,"F90: POINTER CHECK FAILS:",rval," != ",prt_params%wood_density(indx1)
            stop
         end if
      end if

      return
  end subroutine PRTParamsPySet


  subroutine PRTParamsAlloc(fates_NCWD, & 
   fates_history_age_bins, & 
   fates_history_height_bins, & 
   fates_history_size_bins, & 
   fates_history_coage_bins, & 
   fates_hydr_organs, & 
   fates_leafage_class, & 
   fates_litterclass, & 
   fates_pft, & 
   fates_prt_organs, & 
   fates_string_length, & 
   fates_hlm_pftno)
    !
    use FatesConstantsMod, only : fates_unset_r8
	 use FatesConstantsMod, only : fates_unset_int

	implicit none

   integer,intent(in) :: fates_NCWD
   integer,intent(in) :: fates_history_age_bins
   integer,intent(in) :: fates_history_height_bins
   integer,intent(in) :: fates_history_size_bins
   integer,intent(in) :: fates_history_coage_bins
   integer,intent(in) :: fates_hydr_organs
   integer,intent(in) :: fates_leafage_class
   integer,intent(in) :: fates_litterclass
   integer,intent(in) :: fates_pft
   integer,intent(in) :: fates_prt_organs
   integer,intent(in) :: fates_string_length
   integer,intent(in) :: fates_hlm_pftno


    ! LOCALS:
    integer                    :: iv1   ! The parameter incrementer
    integer                    :: iv2
    !------------------------------------------------------------------------

    allocate( prt_params_ptr%var1d(100)) ! Make this plenty large
	 allocate( prt_params_ptr%var2d(100))
	 iv1=0
	 iv2=0

	! POINTER-SPECIFICATION-HERE (DO NOT REMOVE THIS LINE, OR MOVE IT)

	 allocate(prt_params%stress_decid(fates_pft))
	 prt_params%stress_decid(:) = fates_unset_int
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_phen_stress_decid"
	 prt_params_ptr%var1d(iv1)%var_ip   => prt_params%stress_decid
	 prt_params_ptr%var1d(iv1)%vtype    = 2

	 allocate(prt_params%season_decid(fates_pft))
	 prt_params%season_decid(:) = fates_unset_int
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_phen_season_decid"
	 prt_params_ptr%var1d(iv1)%var_ip   => prt_params%season_decid
	 prt_params_ptr%var1d(iv1)%vtype    = 2

	 allocate(prt_params%evergreen(fates_pft))
	 prt_params%evergreen(:) = fates_unset_int
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_phen_evergreen"
	 prt_params_ptr%var1d(iv1)%var_ip   => prt_params%evergreen
	 prt_params_ptr%var1d(iv1)%vtype    = 2

	 allocate(prt_params%senleaf_long_fdrought(fates_pft))
	 prt_params%senleaf_long_fdrought(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_senleaf_long_fdrought"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%senleaf_long_fdrought
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%leaf_long(fates_pft))
	 prt_params%leaf_long(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_senleaf_long_fdrought"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%leaf_long
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%root_long(fates_pft))
	 prt_params%root_long(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_root_long"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%root_long
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%branch_long(fates_pft))
	 prt_params%branch_long(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_branch_turnover"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%branch_long
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%turnover_retrans_mode(fates_pft))
	 prt_params%turnover_retrans_mode(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_turnover_retrans_mode"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%turnover_retrans_mode
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%turnover_carb_retrans(fates_pft,fates_prt_organs))
	 prt_params%turnover_carb_retrans(:,:) = fates_unset_r8
	 iv2 = iv2 + 1
	 prt_params_ptr%var2d(iv2)%var_name = "fates_turnover_carb_retrans"
	 prt_params_ptr%var2d(iv2)%var_rp   => prt_params%turnover_carb_retrans
	 prt_params_ptr%var2d(iv2)%vtype    = 1

	 allocate(prt_params%turnover_nitr_retrans(fates_pft,fates_prt_organs))
	 prt_params%turnover_nitr_retrans(:,:) = fates_unset_r8
	 iv2 = iv2 + 1
	 prt_params_ptr%var2d(iv2)%var_name = "fates_turnover_nitr_retrans"
	 prt_params_ptr%var2d(iv2)%var_rp   => prt_params%turnover_nitr_retrans
	 prt_params_ptr%var2d(iv2)%vtype    = 1

	 allocate(prt_params%turnover_phos_retrans(fates_pft,fates_prt_organs))
	 prt_params%turnover_phos_retrans(:,:) = fates_unset_r8
	 iv2 = iv2 + 1
	 prt_params_ptr%var2d(iv2)%var_name = "fates_turnover_phos_retrans"
	 prt_params_ptr%var2d(iv2)%var_rp   => prt_params%turnover_phos_retrans
	 prt_params_ptr%var2d(iv2)%vtype    = 1

	 allocate(prt_params%grperc(fates_pft))
	 prt_params%grperc(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_grperc"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%grperc
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%nitr_stoich_p1(fates_pft,fates_prt_organs))
	 prt_params%nitr_stoich_p1(:,:) = fates_unset_r8
	 iv2 = iv2 + 1
	 prt_params_ptr%var2d(iv2)%var_name = "fates_prt_nitr_stoich_p1"
	 prt_params_ptr%var2d(iv2)%var_rp   => prt_params%nitr_stoich_p1
	 prt_params_ptr%var2d(iv2)%vtype    = 1

	 allocate(prt_params%nitr_stoich_p2(fates_pft,fates_prt_organs))
	 prt_params%nitr_stoich_p2(:,:) = fates_unset_r8
	 iv2 = iv2 + 1
	 prt_params_ptr%var2d(iv2)%var_name = "fates_prt_nitr_stoich_p2"
	 prt_params_ptr%var2d(iv2)%var_rp   => prt_params%nitr_stoich_p2
	 prt_params_ptr%var2d(iv2)%vtype    = 1

	 allocate(prt_params%phos_stoich_p1(fates_pft,fates_prt_organs))
	 prt_params%phos_stoich_p1(:,:) = fates_unset_r8
	 iv2 = iv2 + 1
	 prt_params_ptr%var2d(iv2)%var_name = "fates_prt_phos_stoich_p1"
	 prt_params_ptr%var2d(iv2)%var_rp   => prt_params%phos_stoich_p1
	 prt_params_ptr%var2d(iv2)%vtype    = 1

	 allocate(prt_params%phos_stoich_p2(fates_pft,fates_prt_organs))
	 prt_params%phos_stoich_p2(:,:) = fates_unset_r8
	 iv2 = iv2 + 1
	 prt_params_ptr%var2d(iv2)%var_name = "fates_prt_phos_stoich_p2"
	 prt_params_ptr%var2d(iv2)%var_rp   => prt_params%phos_stoich_p2
	 prt_params_ptr%var2d(iv2)%vtype    = 1

	 allocate(prt_params%nitr_store_ratio(fates_pft))
	 prt_params%nitr_store_ratio(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_nitr_store_ratio"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%nitr_store_ratio
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%phos_store_ratio(fates_pft))
	 prt_params%phos_store_ratio(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_phos_store_ratio"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%phos_store_ratio
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%organ_id(fates_prt_organs))
	 prt_params%organ_id(:) = fates_unset_int
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_prt_organ_id"
	 prt_params_ptr%var1d(iv1)%var_ip   => prt_params%organ_id
	 prt_params_ptr%var1d(iv1)%vtype    = 2

	 allocate(prt_params%cushion(fates_pft))
	 prt_params%cushion(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_alloc_storage_cushion"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%cushion
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%leaf_stor_priority(fates_pft))
	 prt_params%leaf_stor_priority(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_leaf_stor_priority"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%leaf_stor_priority
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%dbh_repro_threshold(fates_pft))
	 prt_params%dbh_repro_threshold(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_seed_dbh_repro_threshold"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%dbh_repro_threshold
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%seed_alloc_mature(fates_pft))
	 prt_params%seed_alloc_mature(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_seed_alloc_mature"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%seed_alloc_mature
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%seed_alloc(fates_pft))
	 prt_params%seed_alloc(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_seed_alloc_mature"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%seed_alloc
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%fnrt_prof_mode(fates_pft))
	 prt_params%fnrt_prof_mode(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_fnrt_prof_mode"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%fnrt_prof_mode
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%fnrt_prof_a(fates_pft))
	 prt_params%fnrt_prof_a(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_fnrt_prof_a"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%fnrt_prof_a
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%fnrt_prof_b(fates_pft))
	 prt_params%fnrt_prof_b(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_fnrt_prof_b"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%fnrt_prof_b
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%c2b(fates_pft))
	 prt_params%c2b(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_c2b"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%c2b
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%wood_density(fates_pft))
	 prt_params%wood_density(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_wood_density"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%wood_density
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%woody(fates_pft))
	 prt_params%woody(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_woody"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%woody
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%crown(fates_pft))
	 prt_params%crown(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_fire_crown_depth_frac"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%crown
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%slamax(fates_pft))
	 prt_params%slamax(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_leaf_slamax"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%slamax
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%slatop(fates_pft))
	 prt_params%slatop(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_leaf_slatop"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%slatop
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_sai_scaler(fates_pft))
	 prt_params%allom_sai_scaler(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_sai_scaler"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_sai_scaler
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_dbh_maxheight(fates_pft))
	 prt_params%allom_dbh_maxheight(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_dbh_maxheight"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_dbh_maxheight
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_hmode(fates_pft))
	 prt_params%allom_hmode(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_hmode"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_hmode
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_lmode(fates_pft))
	 prt_params%allom_lmode(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_lmode"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_lmode
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_fmode(fates_pft))
	 prt_params%allom_fmode(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_fmode"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_fmode
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_amode(fates_pft))
	 prt_params%allom_amode(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_amode"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_amode
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_cmode(fates_pft))
	 prt_params%allom_cmode(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_cmode"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_cmode
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_smode(fates_pft))
	 prt_params%allom_smode(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_smode"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_smode
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_stmode(fates_pft))
	 prt_params%allom_stmode(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_stmode"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_stmode
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_la_per_sa_int(fates_pft))
	 prt_params%allom_la_per_sa_int(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_la_per_sa_int"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_la_per_sa_int
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_la_per_sa_slp(fates_pft))
	 prt_params%allom_la_per_sa_slp(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_la_per_sa_slp"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_la_per_sa_slp
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_l2fr(fates_pft))
	 prt_params%allom_l2fr(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_l2fr"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_l2fr
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_agb_frac(fates_pft))
	 prt_params%allom_agb_frac(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_agb_frac"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_agb_frac
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_d2h1(fates_pft))
	 prt_params%allom_d2h1(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_d2h1"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_d2h1
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_d2h2(fates_pft))
	 prt_params%allom_d2h2(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_d2h2"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_d2h2
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_d2h3(fates_pft))
	 prt_params%allom_d2h3(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_d2h3"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_d2h3
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_d2bl1(fates_pft))
	 prt_params%allom_d2bl1(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_d2bl1"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_d2bl1
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_d2bl2(fates_pft))
	 prt_params%allom_d2bl2(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_d2bl2"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_d2bl2
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_d2bl3(fates_pft))
	 prt_params%allom_d2bl3(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_d2bl3"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_d2bl3
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_blca_expnt_diff(fates_pft))
	 prt_params%allom_blca_expnt_diff(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_blca_expnt_diff"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_blca_expnt_diff
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_d2ca_coefficient_max(fates_pft))
	 prt_params%allom_d2ca_coefficient_max(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_d2ca_coefficient_max"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_d2ca_coefficient_max
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_d2ca_coefficient_min(fates_pft))
	 prt_params%allom_d2ca_coefficient_min(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_d2ca_coefficient_min"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_d2ca_coefficient_min
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_agb1(fates_pft))
	 prt_params%allom_agb1(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_agb1"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_agb1
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_agb2(fates_pft))
	 prt_params%allom_agb2(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_agb2"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_agb2
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_agb3(fates_pft))
	 prt_params%allom_agb3(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_agb3"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_agb3
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_agb4(fates_pft))
	 prt_params%allom_agb4(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_agb4"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_agb4
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_zroot_max_dbh(fates_pft))
	 prt_params%allom_zroot_max_dbh(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_zroot_max_dbh"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_zroot_max_dbh
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_zroot_max_z(fates_pft))
	 prt_params%allom_zroot_max_z(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_zroot_max_z"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_zroot_max_z
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_zroot_min_dbh(fates_pft))
	 prt_params%allom_zroot_min_dbh(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_zroot_min_dbh"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_zroot_min_dbh
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_zroot_min_z(fates_pft))
	 prt_params%allom_zroot_min_z(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_zroot_min_z"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_zroot_min_z
	 prt_params_ptr%var1d(iv1)%vtype    = 1

	 allocate(prt_params%allom_zroot_k(fates_pft))
	 prt_params%allom_zroot_k(:) = fates_unset_r8
	 iv1 = iv1 + 1
	 prt_params_ptr%var1d(iv1)%var_name = "fates_allom_zroot_k"
	 prt_params_ptr%var1d(iv1)%var_rp   => prt_params%allom_zroot_k
	 prt_params_ptr%var1d(iv1)%vtype    = 1


   numparm1d = iv1
	numparm2d = iv2

   if(debug) print*,"F90: ALLOCATED ",numparm1d," 1D PARAMETERS"
   if(debug) print*,"F90: ALLOCATED ",numparm2d," 2D PARAMETERS"
   if(debug) print*,"FOR ",fates_pft," PFTs"

	numpft = fates_pft

    return
 end subroutine PRTParamsAlloc

end module PRTParamsGeneric