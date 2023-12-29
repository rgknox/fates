

program DriveKmean

  real(r8),parameter :: vais = (/0.2,0.9,1.1,3.0,3.5,4.5,4.6/)
  integer,parameter  :: max_groups = 3

  integer :: group(:)

  allocate(group(ubound(vais,1)))
  
  call ClusterCohortRadKmean(vais,max_groups,group)

  print*,"VAIS: ",vais(:)
  print*,"GROUP: ",group(:)

  
  deallocate(group)
  
end program DriveKmean




subroutine ClusterCohortRadKmean(points,max_groups,groups)

  real(r8),intent(in)  :: points(:)
  integer,intent(in)   :: max_groups  ! maximum allowable number of groups
  integer,intent(out)  :: groups(:)

  integer :: n_points  ! number of cohorts providing VAI
  integer :: n_groups  ! this is "k" the number of groups

  real(r8) :: centroids(200)
  real(r8) :: variances(200)
  integer :: init_groups(200)

  real(r8) :: converge_var
  real(r8), parameter :: converge_tol = 1.e-2  ! relative change in variance
  
  
  n_points = ubound(points,dim=1)

  ! Determine the ideal number of groups

  ! Initialize the groups 

  ! if the ideal number of groups is less than max_groups
  if(n_points <= max_groups) then

     ! Trivial solution, each cohort is it's own element
     do i = 1,n_points
        groups(i) = i
     end do

  else

     n_groups = max_groups

     ! Quasi-uniform initialization matched to closest point
     do i = 1,n_groups
        icent = ceiling(real(n_points,r8)*real(i,r8)/real(n_groups,r8))
        centroids(i) = points(icent)
     end do
     
     varsum = 1.e10_r8
     varsum_del = 1.e10_r8
     do while(varsum_del>varsum_tol)

        ! Update the group list
        groupcount(1:n_groups) = 0
        do i = 1,n_vai
           do j = 1,n_groups
              dist_cent(j) = abs(centroid(j)-points(i))
           end do
           icent = minloc(dist_cent(1:n_groups), 1)
           group(i) = icent
           groupcount(i) = groupcount(i) + 1
        end do

        if(vasum_del<varsum_tol) exit
        
        ! Update the centroids
        centroid(:) = 0._r8
        do i = 1,n_vai
           centroid(group(i)) = centroid(group(i)) + points(i)/real(groupcount(i,r8))
        end do



     end do

  end if

  return
end subroutine ClusterCohortRadKmean
