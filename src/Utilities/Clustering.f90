module Clustering

    use Constants, only: wp, PI2, FLT_EPS
    use Utilities, only: almost_zero, matinv, matdet, sort_indices

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: kMeans
    public :: GMM

    type GaussianList_t
        integer               :: d              ! Space dimension
        integer               :: n              ! # of components
        real(wp), allocatable :: tau(:)         ! [n]: weights
        real(wp), allocatable :: mu(:,:)        ! [d, n]: centroids
        real(wp), allocatable :: cov(:,:,:)     ! [d, d, n]: convariance matrices
        real(wp), allocatable :: covinv(:,:,:)  ! [d, d, n]: inverse covariance matrices
    end type GaussianList_t

contains

pure subroutine kMeans(k, x, xavg, clusters, info)
    !* Arguments *!
    integer,           intent(in)    :: k
    real(wp),          intent(in)    :: x(:,:)
    real(wp),          intent(inout) :: xavg(:,:)
    integer,           intent(out)   :: clusters(:)
    integer, optional, intent(out)   :: info

    !* Local variables *!
    integer, parameter   :: max_iters = 50
    integer, allocatable :: prevClusters(:)
    integer              :: ndims
    integer              :: npts
    integer              :: i


    ndims = size(x, dim=1)
    npts  = size(x, dim=2)

    ! Initial centroids
    call kMeans_compute_clusters(k, ndims, npts, x, xavg, clusters)

    ! Loop until convergence
    do i = 1, max_iters
        prevClusters = clusters
        call kMeans_compute_centroids(k, ndims, npts, x, xavg, clusters)
        call kMeans_compute_clusters(k, ndims, npts, x, xavg, clusters)
        if (all(prevClusters == clusters)) exit
    end do

    ! Check convergence
    if (present(info)) info = merge(-1, i, i > max_iters)

    ! Sort clusters
    call kMeans_sort_clusters(k, ndims, npts, xavg, clusters)

end subroutine kMeans

pure subroutine kMeans_compute_clusters(k, ndims, npts, x, xavg, clusters)
    !* Arguments *!
    integer,  intent(in)  :: k
    integer,  intent(in)  :: ndims
    integer,  intent(in)  :: npts
    real(wp), intent(in)  :: x(ndims,npts)
    real(wp), intent(in)  :: xavg(ndims,k)
    integer,  intent(out) :: clusters(npts)

    !* Local variables *!
    integer  :: i
    integer  :: j
    real(wp) :: dist
    real(wp) :: minDist


    do i = 1, npts
        minDist = huge(1.0_wp)
        do j = 1, k
            dist = norm2( x(:,i) - xavg(:,j) )
            if (dist < minDist) then
                minDist = dist
                clusters(i) = j
            end if
        end do
    end do

end subroutine kMeans_compute_clusters

pure subroutine kMeans_compute_centroids(k, ndims, npts, x, xavg, clusters)
    !* Arguments *!
    integer,  intent(in)  :: k
    integer,  intent(in)  :: ndims
    integer,  intent(in)  :: npts
    real(wp), intent(in)  :: x(ndims,npts)
    real(wp), intent(out) :: xavg(ndims,k)
    integer,  intent(in)  :: clusters(npts)

    !* Local variables *!
    integer :: i
    integer :: ptsInCluster(k)


    xavg = 0.0_wp
    ptsInCluster = 0
    do i = 1, npts
        xavg(:,clusters(i)) = xavg(:,clusters(i)) + x(:,i)
        ptsInCluster(clusters(i)) = ptsInCluster(clusters(i)) + 1
    end do
    do i = 1, k
        if (ptsInCluster(i) == 0) cycle
        xavg(:,i) = xavg(:,i) / real(ptsInCluster(i), kind=wp)
    end do

end subroutine kMeans_compute_centroids

pure subroutine kMeans_sort_clusters(k, ndims, npts, xavg, clusters)
    !* Arguments *!
    integer,  intent(in)    :: k
    integer,  intent(in)    :: ndims
    integer,  intent(in)    :: npts
    real(wp), intent(inout) :: xavg(ndims,k)
    integer,  intent(inout) :: clusters(npts)

    !* Local variables *!
    integer  :: i, j
    real(wp) :: normvec(k)
    real(wp) :: xavg_tmp(ndims,k)
    integer  :: indices(k)
    integer  :: clustermap(k)


    ! Sort the clusters by their distance to the origin
    normvec = norm2(xavg, dim=1)
    call sort_indices(normvec, indices)
    do i = 1, k
        clustermap(indices(i)) = i
    end do

    ! Assign nodes to clusters
    do i = 1, npts
        clusters(i) = clustermap(clusters(i))
    end do

    xavg_tmp = xavg
    do i = 1, k
        j = clustermap(i)
        xavg(:,j) = xavg_tmp(:,i)
    end do

end subroutine kMeans_sort_clusters

pure subroutine GMM(init_nclusters, nclusters, x, xavg, clusters, info)
    !* Arguments *!
    integer,           intent(in)    :: init_nclusters
    integer,           intent(out)   :: nclusters
    real(wp),          intent(in)    :: x(:,:)
    real(wp),          intent(inout) :: xavg(:,:)
    integer,           intent(out)   :: clusters(:)
    integer, optional, intent(out)   :: info

    !* Local variables *!
    real(wp),        parameter :: tol = 1e-2_wp             ! Magic!!
    real(wp),        parameter :: small = 1e4_wp * FLT_EPS  ! Magic!!
    integer,         parameter :: max_iters = 50

    type(GaussianList_t)  :: g
    integer               :: info_
    integer               :: ndims
    integer               :: npts
    integer               :: iter
    real(wp)              :: ll
    real(wp)              :: llprev
    real(wp), allocatable :: R(:,:)


    nclusters = init_nclusters

    ! Avoid computations for the trivial case
    if (nclusters == 1) then
        clusters = 1
        xavg(:,1) = (minval(x, dim=2) + maxval(x, dim=2)) * 0.5_wp
        return
    end if

    ndims = size(x, dim=1)
    npts  = size(x, dim=2)

    ! Initialize
    call GMM_init(g, ndims, nclusters, xavg)

    ! Update gaussians
    allocate(R(npts, nclusters))
    ll = huge(1.0_wp)
    do iter = 1, max_iters
        llprev = ll
        call GMM_Estep(ndims, npts, nclusters, g, x, R, small, ll)
        call GMM_Mstep(ndims, npts, nclusters, g, x, R, small)
        nclusters = g%n
        if (abs((ll - llprev) / ll) <= tol .or. nclusters < 1) exit
    end do

    ! Check convergence
    if (iter > max_iters) then
        call kMeans(nclusters, x, xavg, clusters, info_)
        if (info_ == -1) then
            info_ = -2
        end if
    else
        info_ = 0
    end if
    if (present(info)) info = info_

    ! Reorder clusters
    call GMM_sort_clusters(ndims, npts, nclusters, g, x, xavg, clusters)

end subroutine GMM

pure subroutine GMM_init(g, d, n, xavg)
    !* Arguments *!
    type(GaussianList_t), intent(out) :: g
    integer,              intent(in)  :: d
    integer,              intent(in)  :: n
    real(wp),             intent(in)  :: xavg(d,n)

    !* Local variables *!
    integer  :: i, j


    g%d = d
    g%n = n
    allocate(g%tau(n))
    allocate(g%mu(d,n))
    allocate(g%cov(d,d,n)); g%cov = 0.0_wp
    allocate(g%covinv, source=g%cov)

    g%tau = 1.0_wp / n
    g%mu = xavg
    do i = 1, n
        do j = 1, d
            g%cov(j,j,i) = 1.0_wp
            g%covinv(j,j,i) = 1.0_wp
        end do
    end do

end subroutine GMM_init

pure function GMM_pdf(d, x, mu, S, Sinv)
    !* Arguments *!
    integer,  intent(in) :: d
    real(wp), intent(in) :: x(d)
    real(wp), intent(in) :: mu(d)
    real(wp), intent(in) :: S(d,d)
    real(wp), intent(in) :: Sinv(d,d)
    real(wp)             :: GMM_pdf

    !* Local variables *!
    real(wp) :: xmu(d)
    real(wp) :: det

    xmu = x - mu
    det = matdet(d, S)
    GMM_pdf = exp(-0.5_wp * dot_product(xmu, matmul(Sinv, xmu))) / &
              sqrt(PI2 ** d * det)

end function GMM_pdf

pure subroutine GMM_Estep(ndims, npts, nclusters, g, x, R, small, logL)
    !* Arguments *!
    integer,              intent(in)    :: ndims
    integer,              intent(in)    :: npts
    integer,              intent(in)    :: nclusters
    type(GaussianList_t), intent(in)    :: g
    real(wp),             intent(in)    :: x(:,:)
    real(wp),             intent(inout) :: R(:,:)
    real(wp),             intent(in)    :: small
    real(wp),             intent(out)   :: logL

    !* Local variables *!
    integer  :: i, j
    real(wp) :: den

    logL = 0.0_wp
    do i = 1, npts
        den = 0.0_wp
        do j = 1, nclusters
            R(i,j) = g%tau(j) * GMM_pdf(ndims, x(:,i), g%mu(:,j), g%cov(:,:,j), g%covinv(:,:,j))
            den = den + R(i,j)
        end do

        ! The point might be far from all clusters
        if (almost_zero(den)) then
            R(i,1:nclusters) = 1.0_wp / nclusters
            den = small
        else
            R(i,1:nclusters) = R(i,1:nclusters) / den
        end if
        logL = logL + log(den)
    end do

end subroutine GMM_Estep

pure subroutine GMM_Mstep(ndims, npts, nclusters, g, x, R, tol)
    !* Arguments *!
    integer,              intent(in)    :: ndims
    integer,              intent(in)    :: npts
    integer,              intent(inout) :: nclusters
    type(GaussianList_t), intent(inout) :: g
    real(wp),             intent(in)    :: x(:,:)
    real(wp),             intent(in)    :: R(:,:)
    real(wp),             intent(in)    :: tol

    !* Local variables *!
    integer  :: i, j, k, l, m
    integer  :: info_
    logical  :: deleteit
    real(wp) :: ns
    real(wp) :: xmu(ndims)

    k = 1
    do j = 1, nclusters
        ! "Number of points" in the cluster
        ns = sum(R(:,j))

        ! Delete empty clusters...
        if (ns < tol) then
            deleteit = .true.

        ! ...and update the rest
        else
            g%tau(k) = ns / npts
            g%mu(:,k) = matmul(x, R(:,j)) / ns
            g%cov(:,:,k) = 0.0_wp
            do i = 1, npts
                xmu = x(:,i) - g%mu(:,k)
                do m = 1, ndims
                    do l = 1, ndims
                        g%cov(l,m,k) = g%cov(l,m,k) + R(i,j) * xmu(l) * xmu(m)
                    end do
                end do
            end do
            g%cov(:,:,k) = g%cov(:,:,k) / ns
            call matinv(ndims, g%cov(:,:,k), g%covinv(:,:,k), tol=tol, info=info_)

            deleteit = .false.

            ! If a cluster collapses, limit its size
            if (info_ < 0) then
                do l = 1, ndims
                    g%cov(:,l,k) = 0.0_wp
                    g%cov(l,l,k) = tol**(1.0_wp/ndims) * 1e1_wp
                end do
                call matinv(ndims, g%cov(:,:,k), g%covinv(:,:,k), tol=tol)
            end if

            ! Also look for overlapping clusters
            do l = 1, k-1
                if (all(almost_zero(g%mu(:,k) - g%mu(:,l), tol))) then
                    deleteit = .true.
                    exit
                end if
            end do

        end if

        ! Delete the cluster and step forward onto the next iteration
        if (deleteit) then
            do l = k+1, g%n
                g%tau(l-1) = g%tau(l)
                g%mu(:,l-1) = g%mu(:,l)
                g%cov(:,:,l-1) = g%cov(:,:,l)
                g%covinv(:,:,l-1) = g%covinv(:,:,l)
            end do
            g%n = g%n - 1
        else
            k = k + 1
        end if
    end do

end subroutine GMM_Mstep

pure subroutine GMM_sort_clusters(ndims, npts, nclusters, g, x, xavg, clusters)
    !* Arguments *!
    integer,              intent(in)  :: ndims
    integer,              intent(in)  :: npts
    integer,              intent(in)  :: nclusters
    type(GaussianList_t), intent(in)  :: g
    real(wp),             intent(in)  :: x(:,:)
    real(wp),             intent(out) :: xavg(:,:)
    integer,              intent(out) :: clusters(:)

    !* Local variables *!
    integer  :: i, j
    real(wp) :: val
    real(wp) :: tmp
    real(wp) :: normvec(nclusters)
    integer  :: indices(nclusters)
    integer  :: clustermap(nclusters)


    ! Sort the clusters by their distance to the origin
    normvec = norm2(g%mu(:,1:nclusters), dim=1)
    call sort_indices(normvec, indices)
    do i = 1, nclusters
        clustermap(indices(i)) = i
    end do

    ! Assign nodes to clusters
    do i = 1, npts
        clusters(i) = 0
        val = 0.0_wp
        tmp = 0.0_wp
        do j = 1, nclusters
            tmp = GMM_pdf(ndims, x(:,i), g%mu(:,j), g%cov(:,:,j), g%covinv(:,:,j))
            if (tmp > val) then
                val = tmp
                clusters(i) = clustermap(j)
            end if
        end do
    end do

    do i = 1, nclusters
        j = clustermap(i)
        xavg(:,j) = g%mu(:,i)
    end do

end subroutine GMM_sort_clusters

end module Clustering
