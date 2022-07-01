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
    integer,           intent(in)    :: k
    real(wp),          intent(in)    :: x(:,:)
    real(wp),          intent(inout) :: xAvg(:,:)
    integer,           intent(out)   :: clusters(:)
    integer, optional, intent(out)   :: info

    !* Local variables *!
    integer, parameter   :: max_iters = 50
    integer, allocatable :: prevClusters(:)
    integer              :: nDims
    integer              :: nPts
    integer              :: i


    nPts  = size(x, dim=1)
    nDims = size(x, dim=2)

    ! Initial centroids
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
    real(wp), intent(in)  :: x(nPts,nDims)
    real(wp), intent(in)  :: xAvg(k,nDims)
    integer,  intent(out) :: clusters(nPts)

    !* Local variables *!
    integer  :: i
    integer  :: j
    real(wp) :: dist
    real(wp) :: minDist


    do i = 1, nPts
        minDist = norm2( x(i,:) - xAvg(1,:) )
        clusters(i) = 1
        do j = 2, k
            dist = norm2( x(i,:) - xAvg(j,:) )
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
    real(wp), intent(in)  :: x(nPts,nDims)
    real(wp), intent(out) :: xAvg(k,nDims)
    integer,  intent(in)  :: clusters(nPts)

    !* Local variables *!
    integer :: i
    integer :: ptsInCluster(k)


    xAvg = 0.0_wp
    ptsInCluster = 0
    do i = 1, nPts
        xAvg(clusters(i),:) = xAvg(clusters(i),:) + x(i,:)
        ptsInCluster(clusters(i)) = ptsInCluster(clusters(i)) + 1
    end do
    do i = 1, k
        if (ptsInCluster(i) == 0) cycle
        xAvg(i,:) = xAvg(i,:) / real(ptsInCluster(i), kind=wp)
    end do

end subroutine kMeans_compute_centroids

end module Clustering