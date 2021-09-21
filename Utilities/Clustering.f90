module Clustering

    use Constants, only: wp

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: kMeans

contains

pure subroutine kMeans(k, x, xAvg, clusters, info)
    !* Arguments *!
    integer,           intent(in)  :: k
    real(wp),          intent(in)  :: x(:,:)
    real(wp),          intent(out) :: xAvg(:,:)
    integer,           intent(out) :: clusters(:)
    integer, optional, intent(out) :: info

    !* Local variables *!
    integer, parameter   :: max_iters = 50
    integer, allocatable :: prevClusters(:)
    integer              :: nDims
    integer              :: nPts
    integer              :: ind
    integer              :: i


    nDims = size(x, dim=1)
    nPts  = size(x, dim=2)
    allocate(prevClusters(nPts))
 
    ! Initial centroids
    do i = 1, k
        ind = floor( 1.0_wp + real(i-1, kind=wp)/(k-1) * (nPts-1) )
        xAvg(:,i) = x(:,ind)
    end do
    call kMeans_compute_clusters(k, nDims, nPts, x, xAvg, clusters)
 
    ! Loop until convergence
    do i = 1, max_iters
        prevClusters = clusters
        call kMeans_compute_centroids(k, nDims, nPts, x, xAvg, clusters)
        call kMeans_compute_clusters(k, nDims, nPts, x, xAvg, clusters)
        if (all(prevClusters == clusters)) exit
    end do
 
    ! Check convergence
    if (i > max_iters) clusters = -1
    if (present(info)) info = merge(i-1, i, i > max_iters)

end subroutine kMeans

pure subroutine kMeans_compute_clusters(k, nDims, nPts, x, xAvg, clusters)
    !* Arguments *!
    integer,  intent(in)  :: k
    integer,  intent(in)  :: nDims
    integer,  intent(in)  :: nPts
    real(wp), intent(in)  :: x(nDims,nPts)
    real(wp), intent(in)  :: xAvg(nDims,k)
    integer,  intent(out) :: clusters(nPts)

    !* Local variables *!
    integer  :: i
    integer  :: j
    real(wp) :: dist
    real(wp) :: minDist


    do concurrent (i = 1: nPts)
        minDist = norm2( x(:,i) - xAvg(:,1) )
        clusters(i) = 1
        do j = 2, k
            dist = norm2( x(:,i) - xAvg(:,j) )
            if (dist < minDist) then
                minDist = dist
                clusters(i) = j
            end if
        end do
    end do

end subroutine kMeans_compute_clusters

pure subroutine kMeans_compute_centroids(k, nDims, nPts, x, xAvg, clusters)
    !* Arguments *!
    integer,  intent(in)  :: k
    integer,  intent(in)  :: nDims
    integer,  intent(in)  :: nPts
    real(wp), intent(in)  :: x(nDims,nPts)
    real(wp), intent(out) :: xAvg(nDims,k)
    integer,  intent(in)  :: clusters(nPts)

    !* Local variables *!
    integer  :: i
    real(wp) :: ptsInCluster(k)


    xAvg = 0.0_wp
    ptsInCluster = 0.0_wp
    do i = 1, nPts
        xAvg(:,clusters(i)) = xAvg(:,clusters(i)) + x(:,i)
        ptsInCluster(clusters(i)) = ptsInCluster(clusters(i)) + 1.0_wp
    end do
    do i = 1, k
        xAvg(:,i) = xAvg(:,i) / ptsInCluster(clusters(i))
    end do

end subroutine kMeans_compute_centroids

end module Clustering
