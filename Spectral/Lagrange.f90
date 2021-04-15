!*******************************************************************************
!  MODULE: Legendre
!
!> @author
!> Andres Mateo
!
!> @brief
!> Routines related to the Lagrange polynomials.
!*******************************************************************************

module Lagrange

    use Constants
    use Utilities

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public methods
    public :: barycentric_weights
    public :: lagrange_interpolation
    public :: polynomial_interpolation_matrix
    public :: lagrange_interpolating_polynomials
    public :: lagrange_interpolant_derivative
    public :: polynomial_derivative_matrix
    public :: mth_order_polynomial_derivative_matrix

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 30
!
!> @param[in]  n  number of nodes
!> @param[in]  x  coordinates of the nodes
!> @param[out] w  barycentric weights for the Lagrange interpolator
!···············································································
subroutine barycentric_weights(n, x, w)
    !* Arguments *!
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: x(:)
    real(wp), intent(out) :: w(:)

    !* Local variables *!
    integer :: i
    integer :: j

    ! Initialise weights
    w = 1.0_wp

    ! And iterate over all the nodes
    do i = 2, n
        do j = 1, i-1
            w(j) = w(j) * (x(j) - x(i))
            w(i) = w(i) * (x(i) - x(j))
        end do
    end do
    w = 1.0_wp / w
end subroutine barycentric_weights

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 31
!
!> @param[in]  n  number of nodes
!> @param[in]  x  coordinate of the evaluation point
!> @param[in]  xn nodes of the Lagrange interpolator
!> @param[in]  fn values of the function at the nodes
!> @param[in]  wn barycentric weights for the Lagrange interpolator
!> @param[out] p  interpolated value at 'x'
!···············································································
subroutine lagrange_interpolation(n, x, xn, fn, wn, p)
    !* Arguments *!
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: x
    real(wp), intent(in)  :: xn(:)
    real(wp), intent(in)  :: fn(:)
    real(wp), intent(in)  :: wn(:)
    real(wp), intent(out) :: p

    !* Local variables *!
    integer  :: i
    real(wp) :: num
    real(wp) :: den
    real(wp) :: t

    ! Initialise numerator and denominator
    num = 0.0_wp
    den = 0.0_wp

    ! Sweep all the nodes
    do i = 1, n

        ! Stop if 'x' is one of the nodes
        if (almost_zero(x-xn(i))) then
            p = fn(i)
            return
        end if

        t   = wn(i) / (x - xn(i))
        num = num + t*fn(i)
        den = den + t
    end do

    ! Interpolated value
    p = num / den
end subroutine lagrange_interpolation

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 32
!
!> @param[in]  n  number of nodes (number of columns)
!> @param[in]  xn nodes of the Lagrange interpolator
!> @param[in]  wn barycentric weights for the Lagrange interpolator
!> @param[in]  m  number of evaluation points (number of rows)
!> @param[in]  x  evaluation points ($\xi$)
!> @param[out] T  interpolation matrix
!···············································································
subroutine polynomial_interpolation_matrix(n, xn, wn, m, x, T)
    !* Arguments *!
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: xn(:)
    real(wp), intent(in)  :: wn(:)
    integer,  intent(in)  :: m
    real(wp), intent(in)  :: x(:)
    real(wp), intent(out) :: T(:,:)

    !* Local variables *!
    logical  :: rowHasMatch
    integer  :: i
    integer  :: j
    real(wp) :: r
    real(wp) :: s

    ! For each node 'x'
    do i = 1, m
        rowHasMatch = .false.

        ! Check if the interpolated 'x' equals one of the nodes 'xn'
        do j = 1, n
            T(i,j) = 0.0_wp
            if (almost_zero(x(i)-xn(j))) then
                rowHasMatch = .true.
                T(i,j) = 1.0_wp
            end if
        end do

        ! If there is no coincidence, calculate 'T'
        if (.not. rowHasMatch) then
            s = 0.0_wp
            do j = 1, n
                r      = wn(j) / (x(i) - xn(j))
                T(i,j) = r
                s      = s + r
            end do
            T(i,:) = T(i,:) / s
        end if
    end do
end subroutine polynomial_interpolation_matrix

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 34
!
!> @param[in]  n  number of nodes
!> @param[in]  x  point to evaluate the Lagrange polynomials
!> @param[in]  xn nodes for the Lagrange polynomials
!> @param[in]  wn barycentric weights for the Lagrange interpolator
!> @param[out] l  values of the 'n' Lagrange polynomials at 'x'
!···············································································
subroutine lagrange_interpolating_polynomials(n, x, xn, wn, l)
    !* Arguments *!
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: x
    real(wp), intent(in)  :: xn(:)
    real(wp), intent(in)  :: wn(:)
    real(wp), intent(out) :: l(:)

    !* Local variables *!
    logical  :: xMatchesNode
    integer  :: i
    real(wp) :: r
    real(wp) :: s

    ! Check if the point 'x' equals any node 'xn'
    xMatchesNode = .false.
    do i = 1, n
        l(i) = 0.0_wp
        if (almost_zero(x-xn(i))) then
            l(i) = 1.0_wp
            xMatchesNode = .true.
        end if
    end do

    ! Return in that case
    if (xMatchesNode) return

    ! Determine the polynomials otherwise
    s = 0.0_wp
    do i = 1, n
        r = wn(i) / (x - xn(i))
        l(i) = r
        s = s + r
    end do
    l = l / s
end subroutine lagrange_interpolating_polynomials

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 36
!
!> @param[in]  n  number of nodes
!> @param[in]  x  point to evaluate the Lagrange interpolant derivative
!> @param[in]  xn nodes for the Lagrange polynomials
!> @param[in]  fn values of the function at the nodes
!> @param[in]  wn barycentric weights for the Lagrange interpolator
!> @param[out] pp value of Lagrange polynomial derivative at 'x'
!···············································································
subroutine lagrange_interpolant_derivative(n, x, xn, fn, wn, pp)
    !* Arguments *!
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: x
    real(wp), intent(in)  :: xn(:)
    real(wp), intent(in)  :: fn(:)
    real(wp), intent(in)  :: wn(:)
    real(wp), intent(out) :: pp

    !* Local variables *!
    logical  :: atNode
    integer  :: i
    integer  :: ind
    real(wp) :: num
    real(wp) :: den
    real(wp) :: p
    real(wp) :: r

    ! Initialisation
    atNode = .false.
    num    = 0.0_wp

    ! Check if the point 'x' equals a node 'xn'
    do i = 1, n
        if (almost_zero(x-xn(i))) then
            atNode = .true.
            p   = fn(i)
            den = -wn(i)
            ind = i
        end if
    end do

    ! Faster calculation when $x \equiv x_n$
    if (atNode) then
        do i = 1, n
            if (i /= ind) then
                num = num + wn(i) * (p - fn(i)) / (x - xn(i))
            end if
        end do

    ! Loop over all the nodes otherwise
    else
        den = 0.0_wp
        call lagrange_interpolation(n, x, xn, fn, wn, p)
        do i = 1, n
            r = wn(i) / (x - xn(i))
            num = num + r * (p - fn(i)) / (x - xn(i))
            den = den + r
        end do
    end if

    pp = num /den
end subroutine lagrange_interpolant_derivative

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 37
!
!> @param[in]  n  number of nodes
!> @param[in]  xn nodes for the Lagrange polynomials
!> @param[out] D  First derivative approximation matrix
!···············································································
subroutine polynomial_derivative_matrix(n, xn, D)
    !* Arguments *!
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: xn(:)
    real(wp), intent(out) :: D(:,:)

    !* Local variables *!
    integer  :: i
    integer  :: j
    real(wp) :: wn(n)

    call barycentric_weights(n, xn, wn)
    do i = 1, n
        D(i,i) = 0.0_wp
        do j = 1, n
            if (j /= i) then
                D(i,j) = wn(j) / ( wn(i) * (xn(i) - xn(j)) )
                D(i,i) = D(i,i) - D(i,j)
            end if
        end do
    end do
end subroutine polynomial_derivative_matrix

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 38 (not working yet)
!
!> @param[in]  n  number of nodes
!> @param[in]  m  order of the derivative
!> @param[in]  xn nodes for the Lagrange polynomials
!> @param[out] D  First derivative approximation matrix
!···············································································
subroutine mth_order_polynomial_derivative_matrix(n, m, xn, Dm)
    !* Arguments *!
    integer,  intent(in)  :: n
    integer,  intent(in)  :: m
    real(wp), intent(in)  :: xn(:)
    real(wp), intent(out) :: Dm(:,:)

    !* Local variables *!
    integer  :: i
    integer  :: j
    integer  :: k
    real(wp) :: wn(n)

    ! Initialisation
    call barycentric_weights(n, xn, wn)
    call polynomial_derivative_matrix(n, xn, Dm)

    ! Nothing else to do if $m = 1$
    if (m == 1) return

    ! $m > 1$
    do k = 2, m
        do i = 1, n
            Dm(i,i) = 0.0_wp
            do j = 1, n
                if (j /= i) then
                    Dm(i,j) = k / (xn(i) - xn(j)) * &
                             (wn(j)/wn(i) * Dm(i,i) - Dm(i,j))
                    Dm(i,i) = Dm(i,i) - Dm(i,j)
                end if
            end do
        end do
    end do
end subroutine mth_order_polynomial_derivative_matrix

end module Lagrange
