!*******************************************************************************
!  MODULE: WENO
!
!> @author
!> Andres Mateo
!
!> @brief
!> Implements a generic, one-dimensional WENO scheme, based on the works:
!>   - Shu, C.-W. (1997). Essentially Non-Oscillatory and Weighted Essentially
!>     Non-Oscillatory Schemes for Hyperbolic Conservation Laws.
!>   - Fisher, T. C., & Carpenter, M. H. (2013). High-order entropy stable
!>     finite difference schemes for nonlinear conservation laws: Finite
!>     domains. Journal of Computational Physics, 252, 518–557.
!*******************************************************************************

module WENO

    use Constants, only: wp, eLeft, eRight, fNone
    use Physics, only: NEQS
    use ExceptionsAndMessages

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: WENO_interpolation
    public :: compute_WENO_coeffs
    public :: crj_coeffs
    public :: crj_coeffs_derivatives
    public :: target_coeffs

    ! Explicitly define public types
    public :: WENO_stencil_t


    type WENO_stencil_t
        integer  :: k        = 0   !< size of the stencils
        integer  :: n        = 0   !< number of nodes in the element
        integer  :: leftElem = -1  !< left-most element. of the stencil
        integer  :: numElems = 0   !< number of elements in the stencil
        integer  :: nL       = 0   !< number of nodes to the left
        integer  :: nLL      = 0   !< number of nodes of the left element
        integer  :: nRR      = 0   !< number of nodes of the right element
        integer  :: nSt      = 0   !< number of nodes in the stencil
        real(wp) :: leftInt        !< coeff. for the left bound node
        real(wp) :: rightInt       !< coeff. for the right bound node
        real(wp), allocatable :: cjsir(:,:,:,:)       !< interpolation coefficients
        real(wp), allocatable :: cjsirpl(:,:,:,:,:,:) !< derivative coefficients
        real(wp), allocatable :: dsir(:,:,:)          !< target coefficients
        real(wp), allocatable :: dx(:)                !< width of the cells
    contains
        procedure :: construct => WENO_stencil_constructor
    end type WENO_stencil_t

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Constructor for the WENO_stencil_t class. 'dx' must be set explicitly when
!> all the nodes of the stencil are known.
!
!> @param[in]  n          number of nodes in the element that "owns" the stencil
!> @param[in]  k          size of the WENO stencil
!> @param[in]  hasBounds  .true. if the nodes may be at the element bounds
!> @param[in]  x          coords of the main nodes
!> @param[in]  xc         coords of the cell interfaces
!···············································································
subroutine WENO_stencil_constructor(this, n, k, hasBounds, x, xc)
    !* Arguments *!
    integer,  intent(in) :: n
    integer,  intent(in) :: k
    logical,  intent(in) :: hasBounds
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: xc(:)
    ! Derived types
    class(WENO_stencil_t), intent(inout) :: this

    !* Local variables *!
    real(wp) :: xInt


    this%k = k
    this%n = n

    ! WENO interpolation coefficients
    if (k > 0) then
        allocate(this%cjsir(k, 2, 0:n, 0:k), source=0.0_wp)
        allocate(this%cjsirpl(k, 2, 0:n, 0:k, n, k-1), source=0.0_wp)
        allocate(this%dsir(2, 0:n, 0:k), source=0.0_wp)
    end if

    ! Linear interpolation factors for the boundary nodes
    if (.not. hasBounds) then
        this%leftInt  = 0.0_wp
        this%rightInt = 1.0_wp
    else
        xInt = (xc(2) + xc(1)) / 2.0_wp
        this%leftInt = (xInt-x(1))/(x(2)-x(1))
        xInt = (xc(n+1) + xc(n)) / 2.0_wp
        this%rightInt = (xInt-x(n-1))/(x(n)-x(n-1))
    end if

end subroutine WENO_stencil_constructor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Compute the WENO fluxes at the interfaces of the finite volume cells with
!> widths 'dx'. The interpolation coefficients and target weights must be
!> known at the time this subroutine is called.
!
!> @param[in]   vj       values at the cell centres
!> @param[in]   il       index of the left-most node to be interpolated
!> @param[in]   ir       index of the right-most node to be interpolated
!> @param[in]   side     side of the face to be interpolated
!> @param[in]   stencil  stencil object with information about the interpolation
!> @param[out]  vbar     values at the interfaces
!···············································································
subroutine WENO_interpolation(vj, il, ir, side, stencil, vbar)
    !* Arguments *!
    real(wp), intent(in)  :: vj(:, :)
    integer,  intent(in)  :: il
    integer,  intent(in)  :: ir
    integer,  intent(in)  :: side
    real(wp), intent(out) :: vbar(0:, :)
    ! Derived types
    type(WENO_stencil_t), intent(in) :: stencil

    !* Local variables *!
    real(wp) :: EPS

    integer :: n
    integer :: k
    integer :: iLeft
    integer :: iRight
    integer :: idx
    integer :: r_min
    integer :: r_max
    integer :: r_lb
    integer :: r_ub
    integer :: li
    integer :: i
    integer :: r
    integer :: i_eq

    integer  :: dr
    integer  :: iLoc
    real(wp) :: beta

    real(wp), allocatable :: vst(:,:)
    real(wp), allocatable :: alpha(:,:)
    real(wp), allocatable :: sum_alpha(:)
    real(wp), allocatable :: omega(:)


    ! Initialization
    EPS = 1e-30_wp
    n   = stencil%nSt
    k   = stencil%k

    ! Set the limits depending on the interpolation type and physical boundaries
    ! If FD WENO is used, the faces are associated to their left cell
    iLeft  = iL
    iRight = iR
    iLoc   = 0
    idx    = iLeft
    r_min  = 0
    r_max  = k

    ! Update some parameters to avoid reading "outside" the mesh
    if (side == eLeft) then
        if (iLeft == 0) then
            iLeft = iLeft+1
            iLoc  = iLoc+1
        end if
        r_max = k-1

    else if (side == eRight) then
        if (iRight == n) then
            iRight = iRight-1
        end if
        idx   = idx+1
        r_min = 1

    ! Finite Differences cases
    else if (iLeft == 0) then
        iLeft = iLeft+1
        iLoc  = iLoc+1

    else if (iRight == n) then
        iRight = iRight-1

    end if

    allocate(vst(0:k, NEQS))
    allocate(alpha(0:k, NEQS))
    allocate(sum_alpha(NEQS))
    allocate(omega(NEQS))

    ! Loop over all the nodes of the mesh given by 'vj'
    do i = iLeft, iRight

        ! Range of 'r' usable for this node
        r_lb = max(r_min, k-i)
        r_ub = min(r_max, n-i)
        dr   = r_ub - r_lb + 1

        ! Compute the sensors for each stencil
        do i_eq = 1, NEQS

            li = i - (k-1) + r_lb

            do r = r_lb, r_ub

                li = i - (k-1) + r
                vst(r,i_eq) = dot_product(stencil%cjsir(:,side,iLoc,r), &
                                          vj(li:li+k-1,i_eq))

                beta = compute_beta(vj(li:li+k-1,i_eq), stencil%dx(idx), &
                                    stencil%cjsirpl(:,side,iLoc,r,:,:))

                alpha(r,i_eq) = stencil%dsir(side,iLoc,r) / (EPS + beta)**2

            end do
        end do

        ! Denominator of 'omega'
        sum_alpha = sum(alpha(r_lb:r_ub,:), dim=1)

        ! Final WENO weights and values at the interfaces
        vbar(iLoc,:) = 0.0_wp
        do r = r_lb, r_ub

            omega = alpha(r,:) / sum_alpha
            vbar(iLoc,:) = vbar(iLoc,:) + omega * vst(r,:)

        end do

        iLoc = iLoc + 1
        idx  = idx  + 1

    end do

end subroutine WENO_interpolation

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Compute all the coefficients required to implement a WENO scheme.
!
!> @param[in]   x        nodes where the derivative will be interpolated
!> @param[in]   xc       interfaces between finite volumes
!> @param[in]   i        index of the cell inside the given stencil (xc)
!> @param[in]   k        size of the stencils
!> @param[in]   r_lower  lower bound of the r-index
!> @param[in]   r_upper  upper bound of the r-index
!> @param[out]  cjr      coefficients of the interpolating polynomial
!> @param[out]  cjrpl    coefficients of the derivative polynomial
!> @param[out]  d        target weights of the WENO scheme
!> @param[out]  c2j      coeff. of the wider interpolating polynomial (optional)
!···············································································
subroutine compute_WENO_coeffs(x, xc, i, k, r_lower, r_upper, cjr, cjrpl, d, c2j)
    !* Arguments *!
    real(wp),           intent(in)  :: x(:)
    real(wp),           intent(in)  :: xc(:)
    integer,            intent(in)  :: i
    integer,            intent(in)  :: k
    integer,            intent(in)  :: r_lower
    integer,            intent(in)  :: r_upper
    real(wp),           intent(out) :: cjr(:, 0:)
    real(wp),           intent(out) :: cjrpl(:, 0:, :, :)
    real(wp),           intent(out) :: d(0:)
    real(wp), optional, intent(out) :: c2j(:)

    !* Local variables *!
    integer :: r
    integer :: l

    integer :: dr
    integer :: k2
    integer :: li

    real(wp) :: c2kj(2*k)


    ! First, set outputs to 0
    cjr   = 0.0_wp
    cjrpl = 0.0_wp
    d     = 0.0_wp

    ! Range of available stencils
    dr = r_upper - r_lower + 1
    k2 = k + dr - 1

    ! For each stencil
    do r = r_lower, r_upper

        ! Get interpolation coefficients
        call crj_coeffs(xc=xc, k=k, i=i, r=r, left_index=li, &
                        cj=cjr(:,r))

        ! Get derivative coefficients
        do l = 1, k-1
            call crj_coeffs_derivatives(x=x, xc=xc,             &
                                        order=l, k=k, i=i, r=r, &
                                        left_index=li,          &
                                        cjp=cjrpl(:,r,:,l))
        end do

    end do

    ! Get greater interpolation and target coefficients
    call crj_coeffs(xc=xc, k=k2, i=i, r=r_upper, left_index=li, cj=c2kj)
    call target_coeffs(cjr(:,r_lower:r_upper), c2kj, d(r_lower:r_upper))

    ! Update the output if required
    if (present(c2j)) then
        c2j = 0.0_wp
        c2j(r_lower:r_lower+k2-1) = c2kj(r_lower:r_lower+k2-1)
    end if

end subroutine compute_WENO_coeffs

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Compute the derivative of order 'order' of the Lagrange polynomial for the
!> solution interpolator \f$ p(x) \f$. 'xc' are the interfaces of the cells,
!> whereas 'x' is the node where the interpolation is actually performed.
!
!> @param[in]   x           nodes where the coefficients are computed
!> @param[in]   xc          interfaces between finite volumes
!> @param[in]   order       order of the derivative
!> @param[in]   k           size of the stencils
!> @param[in]   i           index of the 'main' node of the stencils
!> @param[in]   left_index  index of the left-most node of the stencil
!> @param[out]  cjp         coefficients of the interpolating polynomial
!···············································································
subroutine crj_coeffs_derivatives(x, xc, order, k, i, r, left_index, cjp)
    !* Arguments *!
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: xc(:)
    integer,  intent(in)  :: order
    integer,  intent(in)  :: k
    integer,  intent(in)  :: i
    integer,  intent(in)  :: r
    integer,  intent(out) :: left_index
    real(wp), intent(out) :: cjp(:,:)  ! j, p

    !* Local variables *!
    real(wp) :: stencil_xc(k+1)
    real(wp) :: cm(size(x),k)
    integer  :: j
    integer  :: m

    ! Index of the left-most 'face'
    left_index = i - (k-1) + r

    ! Check bounds first (not so obvious with this reconstruction)
    if ((left_index < 1) .or. (left_index+k-1 > ubound(xc, dim=1))) then
        call printError("WENO", "r must lie in the interval [0, k]")
    end if

    ! Points of the stencil 'r'
    stencil_xc = xc(left_index:left_index+k)

    do m = 1, k
        cm(:,m) = expansion(0, order+1, m, [m], x, stencil_xc)
    end do

    cjp = 0.0_wp
    do j = 1, k

        do m = j, k
            cjp(j,:) = cjp(j,:) + cm(:,m)
        end do

        cjp(j,:) = cjp(j,:) * (xc(left_index+j)-xc(left_index+j-1))

    end do

end subroutine crj_coeffs_derivatives

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Utility to compute the inner sums and products of the derivatives of the
!> Lagrange polynomials. To be called only from 'crj_coeffs_derivatives'.
!
!> @param[in]  l      current level of recursion
!> @param[in]  depth  max. depth of recursion
!> @param[in]  j      index of the first sum
!> @param[in]  inds   all the indices to be avoided
!> @param[in]  x      interpolated nodes
!> @param[in]  xc     interfaces between cells
!
!> @return            sum of levels lower than 'l'
!···············································································
recursive function expansion(l, depth, j, inds, x, xc) result(next)
    !* Arguments *!
    integer  :: l
    integer  :: depth
    integer  :: j
    integer  :: inds(:)
    real(wp) :: x(:)
    real(wp) :: xc(0:)

    !* Return value *!
    real(wp) :: next(size(x))

    !* Local variables *!
    integer  :: i
    integer  :: k
    real(wp) :: next_level(size(x))

    ! Number of cells in the stencil
    k = size(xc) - 1

    ! The last step is always a product
    if (depth == l) then

        next = 1.0_wp
        do i = 0, k

            if (any(inds == i)) cycle
            next = next * (x - xc(i)) / (xc(j) - xc(i))

        end do

    ! Sum the values from lower steps
    else

        next = 0.0_wp
        do i = 0, k

            if (any(inds == i)) cycle
            next_level = expansion(l+1, depth, j, [inds, i], x, xc)
            next = next + next_level / (xc(j) - xc(i))

        end do

    end if

end function expansion

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Interpolate the value of the interface \f$ i+1/2 \f$ from the stencil
!> \f$ r \f$ of size \f$ k \f$.
!
!> @param[in]   xc          interfaces between finite volumes
!> @param[in]   k           size of the stencils
!> @param[in]   i           index of the 'main' node of the stencils
!> @param[in]   r           left offset of the stencil
!> @param[in]   left_index  index of the left-most node of the stencil
!> @param[out]  cj          coefficients of the interpolating polynomial
!···············································································
subroutine crj_coeffs(xc, k, i, r, left_index, cj)
    !* Arguments *!
    real(wp), intent(in)  :: xc(:)
    integer,  intent(in)  :: k
    integer,  intent(in)  :: i
    integer,  intent(in)  :: r
    integer,  intent(out) :: left_index
    real(wp), intent(out) :: cj(:)

    !* Local variables *!
    real(wp) :: cj2(size(cj), 1)


    ! Just a wrapper around the routine for the derivatives, with \f$ l=0 \f$
    call crj_coeffs_derivatives(x=[xc(i+1)], xc=xc, order=0, k=k, i=i, r=r, &
                                left_index=left_index, cjp=cj2)

    cj = cj2(:,1)

end subroutine crj_coeffs

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Determine the target weights for the stencils of 'ck' so that their addition
!> yields the same polynomial as that of 'c2k'.
!
!> @param[in]   ck   coefficients of the stencils
!> @param[in]   c2k  coefficients of the greater stencil
!> @param[out]  d    target weights
!···············································································
subroutine target_coeffs(ck, c2k, d)
    !* Arguments *!
    real(wp), intent(in)  :: ck(:,0:)  ! j, r
    real(wp), intent(in)  :: c2k(:)    ! j
    real(wp), intent(out) :: d(0:)     ! r

    !* Local variables *!
    real(wp), allocatable :: b(:)
    integer               :: i
    integer               :: k
    integer               :: st
    integer               :: nst

    ! Some initializations
    nst = size(ck, dim=2)
    k   = size(ck, dim=1)
    st  = k-1

    ! Copy of the larger stencil
    allocate(b(0:size(c2k)-1))
    b = c2k

    do i = 0, nst-1
        d(i) = b(i) / ck(1,i)
        b(i:i+st) = b(i:i+st) - d(i) * ck(:,i)
    end do

end subroutine target_coeffs

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Determine the target weights for the stencils of 'ck' so that their addition
!> yields the same polynomial as that of 'c2k'.
!
!> @param[in]  vj    values of the solution at the nodes of the stencil
!> @param[in]  dx    size of the 'main' cell
!> @param[in]  cjpl  coefficients of the 'k-1' derivatives (modified)
!
!> @return         value of the sensor
!···············································································
function compute_beta(vj, dx, cjpl) result(beta)
    !* Arguments *!
    real(wp), intent(in) :: vj(:)
    real(wp), intent(in) :: dx
    real(wp), intent(in) :: cjpl(:,:,:)  ! j, p, l

    !* Return value *!
    real(wp) :: beta

    !* Local variables *!
    integer  :: k
    integer  :: l


    k    = size(vj)
    beta = 0.0_wp

    do l = 1, k-1
        beta = beta + sum( matmul(vj, cjpl(:,:,l))**2 ) * dx**(2*l-1)
    end do

end function compute_beta

end module WENO
