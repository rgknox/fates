module ClusteringMod

  use shr_log_mod      , only: errMsg => shr_log_errMsg
  use shr_sys_mod      , only: shr_sys_abort
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals     , only : fates_log
  
contains
  
  subroutine ClusterKMeanNaive1D(points,max_clusters,max_iter,conv_tol,clusterids,n_clusters,out_code)

    ! -----------------------------------------------------------------------------------------
    ! This is a simple naive k-means clustering algorithm. The algorithm was written based on
    ! the concepts explained in wikipedia and are not based on any pre-existing code.
    !
    ! This algorithm will evaluate a 1D vector of points, and based on a total number of clusters,
    ! will assign those points to a cluster, with the objective of minimizing the distance between
    ! all points and the closest cluster centroid. It does this by iteratively updating
    ! the cluster centroid positions based on the mean of the cluster constituent points.
    ! This method initializes the cluster centroids based, by matching the centroid with
    ! one of the existing points in fractional rank order. ie the 5th out of 10 clusters, will
    ! be roughly lined up wih the point half-way down the list. This algorithm is used
    ! for fates, and therefore the points will most likely be ordered.
    ! Ryan Knox 1/24
    ! -----------------------------------------------------------------------------------------

   
    
    implicit none

    ! Input Arguments
    real(r8),intent(in)  :: points(:)                ! the values of the points we are clustering
    integer,intent(in)   :: max_clusters             ! maximum allowable number of clusters
    integer,intent(in)   :: max_iter                 ! maximum number of iterations to allow
    real(r8),intent(in)  :: conv_tol                 ! convergence tolerance

    ! Output Arguments
    integer,intent(out)  :: n_clusters               ! The number of clusters generated
    ! this could be lower than max in the future
    integer,intent(out)  :: clusterids(:)             ! the cluster index associated with each point
    integer,intent(out)  :: out_code                 ! output code:  -1 = trival solution (few points)
    !                0 = solution found
    !                1 = (poot) max iterations exceeded

    integer :: n_points                              ! number of points to be clustered
    integer,parameter :: abs_max_clusters = 500      ! don't allow numbers of clusters larger than this
    real(r8) :: centroids(abs_max_clusters)          ! the cluster centroids (means)
    real(r8) :: dist_cent(abs_max_clusters)          ! distance between point and centroid
    integer  :: n_clust_pts(abs_max_clusters)        ! number of points in each cluster 
    integer  :: i                                    ! Point index
    integer  :: icent                                ! Point index for choosing a starting centroid
    integer  :: j                                    ! cluster index

    ! Convergence criteria, two types:
    ! 1) Either use a variance tolerance, or
    ! 2) Check if the cluster associations have
    ! changed over some number of trys
    integer, parameter :: max_conv_mem = 10
    integer :: clusterid_sum
    integer :: clusterid_sums(max_conv_mem) ! Unique cluster mapping sums from current and previous iterations
    integer  :: n_conv_mem             ! Number of convergence indices to remember when evaluating
    real(r8) :: varsum                 ! Sum of variances between points and centers
    real(r8) :: varsum_old             ! Sum of variances between points and centers from previous iteration
    integer  :: conv_iter              ! convergence iterator

    character(len=*), parameter, private :: sourcefile = &
         __FILE__
    
    n_points = ubound(points,dim=1)

    ! Determine the ideal number of clusters (COMING SOON?)
    ! We are just using the max clusters for the number
    ! of clusters. Alternatively, an analysis could determine
    ! if the max clusters is unnecessary and propose a smaller
    ! number of clusters. This will be available in a next version
    ! potentially.

    n_clusters = max_clusters

    ! If the ideal number of clusters is less than max_clusters
    if_trivial: if(n_points <= n_clusters) then

       ! Trivial solution, not many points so
       ! each point is it's own cluster
       do i = 1,n_points
          clusterids(i) = i
       end do
       n_clusters = n_points
       out_code = 0

    else

       if(n_clusters>abs_max_clusters)then
          write(fates_log(),*) "The K-means algorithm is capped at ",abs_max_clusters,"clusters"
          write(fates_log(),*) "You specified: ",n_clusters
          write(fates_log(),*) "Stopping Run:"
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       ! Quasi-uniform initialization matched to closest point
       do j = 1,n_clusters
          i = ceiling(real(n_points,r8)*real(j,r8)/real(n_clusters,r8))
          centroids(j) = points(i)
       end do

       ! If a positive variance tolerance is provided, then
       ! we are using variance tolerance as a convergence metric.
       ! if a negative value is provided, we are assuming convergence
       ! is determined if assignments do not change

       if(abs(conv_tol)<1.e-10)then
          write(fates_log(),*)"The K-means clustering algorithm needs a more"
          write(fates_log(),*)"realistic convergence criteria than: ",conv_tol
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       if(conv_tol<0._r8)then
          ! Number of memory slots to compare
          n_conv_mem = nint(abs(conv_tol))
          ! Initialize the previous indices as a negative
          clusterid_sums(1:n_conv_mem) = -100
       end if

       conv_iter = 0
       adjust_centroids: do

          ! There are two different convergence criteria
          ! which we use to exit this loop, no need to
          ! specify a while. We just keep going until we dont

          ! I.
          ! -------------------------------------------------------------------
          ! Update the current cluster list and the sum of the
          ! variances which may be used as goodness of fit
          ! of goodness of fit

          n_clust_pts(1:n_clusters) = 0
          varsum = 0._r8
          clusterid_sum = 0
          do i = 1,n_points

             do j = 1,n_clusters
                dist_cent(j) = (centroids(j)-points(i))**2._r8
             end do

             j = minloc(dist_cent(1:n_clusters), 1)
             clusterids(i)  = j  ! This is the cluster assignment
             n_clust_pts(j) = n_clust_pts(j) + 1

             ! Create a unique integer that is associated
             ! with the layout of the groupings
             varsum = varsum + dist_cent(j)
             clusterid_sum = clusterid_sum + j*(10000+i)
          end do

          ! II.
          ! If too many iterations, leave, but report with
          ! a substandard output code

          if(conv_iter==max_iter)then
             out_code = -2
             exit
          end if

          ! III
          ! Test the fitness of the cluster by either of the
          ! two methods, using normalized variance or unique ids

          if(conv_tol>0._r8) then

             if( abs(varsum_old-varsum)/(varsum+varsum_old)<conv_tol) then
                out_code = conv_iter
                exit
             end if
             varsum_old = varsum

          else

             ! Test if assignments are different
             if( all(clusterid_sum==clusterid_sums(1:n_conv_mem))   ) then
                out_code = conv_iter
                exit
             else
                do j = n_conv_mem,2,-1
                   clusterid_sums(j) = clusterid_sums(j-1)
                end do
                clusterid_sums(1) = clusterid_sum
             end if
          end if

          ! IV.
          ! If the point-to-cluster associations have not stabilized, 
          ! update the centroids and then try again

          centroids(1:n_clusters) = 0._r8
          do i = 1,n_points
             j = clusterids(i)
             centroids(j) = centroids(j) + points(i)/real(n_clust_pts(j),r8)
          end do
          conv_iter = conv_iter + 1
       end do adjust_centroids

    end if if_trivial

    return
  end subroutine ClusterKMeanNaive1D
end module ClusteringMod



program DriveKmean

  use ClusteringMod,only: ClusterKMeanNaive1d

  integer, parameter :: r8 = selected_real_kind(12)

  integer, parameter :: max_pts = 1000
  real(r8) :: vais(max_pts)

  !    (/0.2,0.9,1.1,3.0,3.5,4.5,4.6/)
  integer, parameter :: max_clusters = 3
  integer, parameter :: max_iter = 200
  real(r8),parameter :: conv_tol = -3._r8   ! five iterations must match

  integer             :: clusterids(max_pts)
  integer             :: n_clusters
  integer             :: out_code
  integer             :: n_vai


  n_vai = 20
  call random_number(vais(1:n_vai))

  call ClusterKMeanNaive1D(vais(1:n_vai),max_clusters,max_iter,conv_tol, &
       clusterids(1:n_vai),n_clusters,out_code)

  print,"OUT CODE: ",out_code

  do j = 1,n_clusters
     print*,"Cluster: ",j
     do i = 1,n_vai
        if(clusterids(i)==j)then
           print*,"point:",vais(i)
        end if
     end do
  end do

end program DriveKmean

