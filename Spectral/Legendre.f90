!*******************************************************************************
!  MODULE: Legendre
!
!> @author
!> Andres Mateo
!
!> @brief
!> Routines related to the Legendre polynomials.
!*******************************************************************************

module Legendre

    use Constants
    use Utilities
    use ExceptionsAndMessages

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public methods
    public :: legendre_polynomial
    public :: legendre_polynomial_derivative
    public :: legendre_gauss_nodes_weights
    public :: legendre_gauss_lobatto_nodes_weights

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 20
!
!> @param[in]  k  degree of the Legendre polynomial
!> @param[in]  x  node to calculate 'l'
!> @param[out] l  value of the Legendre polynomial
!···············································································
subroutine legendre_polynomial(k, x, l)
    !* Arguments *!
    integer,  intent(in)  :: k
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: l

    !* Local variables *!
    real(wp) :: l1
    real(wp) :: l2
    integer  :: i

    ! Recursion only if k>1
    select case (k)

    case (0)
        l = 1.0_wp

    case (1)
        l = x

    case default
        l2 = 1.0_wp
        l1 = x
        do i = 2, k
            l = (2.0_wp*i-1.0_wp)/i*x*l1 + (1.0_wp-i)/i*l2
            l2 = l1
            l1 = l
        end do

    end select
end subroutine legendre_polynomial

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 22
!
!> @param[in]  k  degree of the Legendre polynomial
!> @param[in]  x  node to calculate 'l'
!> @param[out] l  value of the Legendre polynomial
!> @param[out] lp value of the Legendre polynomial derivative
!···············································································
subroutine legendre_polynomial_derivative(k, x, l, lp)
    !* Arguments *!
    integer,  intent(in)  :: k
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: l
    real(wp), intent(out) :: lp

    !* Local variables *!
    real(wp) :: l1
    real(wp) :: l2
    real(wp) :: lp1
    real(wp) :: lp2
    integer  :: i

    ! Recursion only if k>1
    select case (k)

    case (0)
        l  = 1.0_wp
        lp = 0.0_wp

    case(1)
        l  = x
        lp = 1.0_wp

    case default
        l2  = 1.0_wp
        l1  = x
        lp2 = 0.0_wp
        lp1 = 1.0_wp
        do i = 2, k
            l  = (2.0_wp*i-1.0_wp)/i*x*l1 + (1.0_wp-i)/i*l2
            lp = lp2 + (2.0_wp*i-1)*l1
            l2  = l1
            l1  = l
            lp2 = lp1
            lp1 = lp
        end do

    end select
end subroutine legendre_polynomial_derivative

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 23
!
!> @param[in]  n  number of zeros
!> @param[out] z  zeros of the n-th Legendre polynomial
!> @param[out] w  corresponding weights for the Gauss quadrature
!···············································································
subroutine legendre_gauss_nodes_weights(n, z, w)
    !* Arguments *!
    integer,  intent(in)  :: n
    real(wp), intent(out) :: z(:)
    real(wp), intent(out) :: w(:)

    !* Local variables *!
    integer,  parameter :: maxIters = 50
    real(wp), parameter :: relTol   = 4.0_wp*epsilon(1.0_wp)
    integer  :: i
    integer  :: j
    real(wp) :: l
    real(wp) :: lp
    real(wp) :: delta

    ! Iterative method only if n>2
    select case (n)

    case (1)
        z(1) = 0.0_wp
        w(1) = 2.0_wp

    case (2)
        z(1) = -sqrt(1.0_wp/3.0_wp)
        z(2) = -z(1)
        w(1) = 1.0_wp
        w(2) = w(1)

    case (3:)

        do i = 1, n/2
            ! Initial guess
            z(i) = -cos((2.0_wp*i-1.0_wp)/(2.0_wp*n)*PI)

            ! Newton-Raphson algorithm
            do j = 1, maxIters
                call legendre_polynomial_derivative(n, z(i), l, lp)
                delta = -l/lp
                z(i) = z(i) + delta
                if (abs(delta)<relTol*abs(z(i))) then
                    exit
                end if
            end do

            call legendre_polynomial_derivative(n, z(i), l, lp)
            z(n+1-i) = -z(i)
            w(i)     = 2.0_wp / ((1.0_wp-z(i)**2)*lp**2)
            w(n+1-i) = w(i)
        end do

        ! In case 'n' is odd, add z=0 as a zero
        if (mod(n, 2) == 1) then
            call legendre_polynomial_derivative(n, 0.0_wp, l, lp)
            z(n/2+1) = 0.0_wp
            w(n/2+1) = 2.0_wp/lp**2
        end if

    case default

        call printError("Legendre", "There must be at least 1 Gauss node.")

    end select
end subroutine legendre_gauss_nodes_weights

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 24
!
!> @param[in]  k  order of the Legendre polynomial
!> @param[in]  x  node of evaluation
!> @param[out] q  difference $L_{k+1}-L_{k-1}$
!> @param[out] qp difference $L'_{k+1}-L'_{k-1}$
!> @param[out] l  value of the Legendre polynomial
!···············································································
subroutine q_L_evaluation(k, x, q, qp, l)
    !* Arguments *!
    integer,  intent(in)  :: k
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: q
    real(wp), intent(out) :: qp
    real(wp), intent(out) :: l

    !* Local variables *!
    real(wp) :: l1
    real(wp) :: l2
    real(wp) :: lp
    real(wp) :: lp1
    real(wp) :: lp2
    integer  :: i

    ! Iterate until $L_{k+1}$
    l2  = 1.0_wp
    l1  = x
    lp2 = 0.0_wp
    lp1 = 1.0_wp
    do i = 2, k
        l  = (2.0_wp*i-1.0_wp)/i*x*l1 + (1.0_wp-i)/i*l2
        lp = lp2 + (2.0_wp*i-1.0_wp)*l1
        l2  = l1
        l1  = l
        lp2 = lp1
        lp1 = lp
    end do
    l  = (2.0_wp*k+1.0_wp)/(k+1.0_wp)*x*l1 - k/(k+1.0_wp)*l2
    lp = lp2 + (2.0_wp*k+1.0_wp)*l1
    q  = l - l2
    qp = lp - lp2
    l  = l1
end subroutine q_L_evaluation

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 25
!
!> @param[in]  n  number of zeros
!> @param[out] z  zeros of the n-th legendre polynomial
!> @param[out] w  corresponding weights for the Gauss quadrature
!···············································································
subroutine legendre_gauss_lobatto_nodes_weights(n, z, w)
    !* Arguments *!
    integer,  intent(in)  :: n
    real(wp), intent(out) :: z(:)
    real(wp), intent(out) :: w(:)

    !* Local variables *!
    integer,  parameter :: maxIters = 20
    real(wp), parameter :: relTol   = 2.0_wp*epsilon(1.0_wp)
    integer  :: i
    integer  :: j
    real(wp) :: q
    real(wp) :: qp
    real(wp) :: l
    real(wp) :: delta

    ! Iterative method only if n>2
    if (n == 2) then
        z(1) = -1.0_wp
        z(2) = 1.0_wp
        w(1) = 1.0_wp
        w(2) = w(1)

    else if (n > 2) then
        z(1) = -1.0_wp
        z(n) = 1.0_wp
        w(1) = 2.0_wp/(n*(n-1))
        w(n) = w(1)
        do i = 2, n/2
            ! Initial guess
            z(i) = -cos(((i-0.75_wp)*PI)/(n-1) &
                 - 3.0_wp/(8.0_wp*(n-1)*PI)/(i-0.75_wp))

            ! Newton-Raphson algorithm
            do j = 1, maxIters
                call q_L_evaluation(n-1, z(i), q, qp, l)
                delta = -q/qp
                z(i) = z(i) + delta
                if (abs(delta)<relTol*abs(z(i))) then
                    exit
                end if
            end do

            call legendre_polynomial(n-1, z(i), l)
            z(n+1-i) = -z(i)
            w(i)     = 2.0_wp / (n*(n-1)*l**2)
            w(n+1-i) = w(i)
        end do

    else

        call printError("Legendre", &
                        "There must be at least 2 Gauss-Lobatto nodes.")

    end if

    ! In case 'n' is odd, add z=0 as a zero
    if (mod(n, 2) == 1) then
        call legendre_polynomial(n-1, 0.0_wp, l)
        z(n/2+1) = 0.0_wp
        w(n/2+1) = 2.0_wp / (n*(n-1)*l**2)
    end if
end subroutine legendre_gauss_lobatto_nodes_weights

end module Legendre
