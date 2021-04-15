!*******************************************************************************
!  MODULE: NodalDG_GaussLobattoClass
!
!> @author
!> Andres Mateo
!
!> @brief
!> Defines the NodalDG_GLL class.
!*******************************************************************************

module NodalDG_GaussLobattoClass

    use Constants
    use StdExpansionClass
    use Utilities
    use ExceptionsAndMessages
    use Legendre
    use Lagrange
    use WENO

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public types
    public :: NodalDG_GLL

!···············································································
!> @class NodalDG_GLL
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 58.
!···············································································
    type, extends(Stdexp_t) :: NodalDG_GLL
        real(wp), allocatable :: Q(:,:)       !< SBP matrix operator
        real(wp), allocatable :: Ds(:,:)      !< \f$D^{\sharp}\f$ matrix
        real(wp), allocatable :: inpC(:,:)    !< main to complementary grid
    contains
        procedure, private :: nodalDG_GLL_assignment
        procedure :: construct     => nodalDG_GLL_constructor
        procedure :: extrapolate   => nodalDG_GLL_extrapolator
        generic   :: assignment(=) => nodalDG_GLL_assignment
    end type NodalDG_GLL

contains


!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 59.
!
!> @param[in]  n  number of nodes
!···············································································
subroutine nodalDG_GLL_constructor(this, n)
    !* Arguments *!
    integer, intent(in) :: n
    ! Derived types
    class(NodalDG_GLL), intent(inout) :: this

    !* Local variables *!
    integer               :: i
    integer               :: j
    real(wp), allocatable :: xn(:)
    real(wp), allocatable :: avgTmp(:,:)

    ! Set the number of nodes
    this%hasBounds = .true.
    this%n         = n
    this%nt        = eGaussLobatto

    ! Allocate attributes
    ! 'inp' is not allocated here, it is set when needed by adaptation rutines
    ! 'fineRes' is allocated when saving the mesh with a 'Printer' object
    allocate(this%wb(n), source=0.0_wp)
    allocate(this%D(n, n), source=0.0_wp)
    allocate(this%Q(n, n), source=0.0_wp)
    allocate(this%Dh(n, n), source=0.0_wp)
    allocate(this%Ds(n ,n), source=0.0_wp)
    allocate(this%bdMode(2, n), source=0.0_wp)
    allocate(this%x(n), source=0.0_wp)
    allocate(this%xc(0:n), source=0.0_wp)
    allocate(this%avgInp(n, n), source=0.0_wp)
    allocate(this%w(n), source=0.0_wp)
    allocate(this%iw(n), source=0.0_wp)
    allocate(this%modeN2(n), source=0.0_wp)
    allocate(this%Fwd(n, n), source=0.0_wp)
    allocate(this%Bwd(n, n), source=0.0_wp)

    ! Nodes and weights of the quadrature
    call legendre_gauss_lobatto_nodes_weights(n, this%x, this%w)
    call barycentric_weights(n, this%x, this%wb)
    this%iw = 1.0_wp / this%w

    ! Calculate the interpolating polynomials at the boundaries
    this%bdMode = 0.0_wp
    this%bdMode(eLeft,1)  = 1.0_wp
    this%bdMode(eRight,n) = 1.0_wp

    ! Calculate the derivative matrices
    call polynomial_derivative_matrix(n, this%x, this%D)
    do i = 1, n
        do j = 1, n
            this%Q(j,i)  = this%D(j,i) * this%w(j)
            this%Dh(i,j) = -this%Q(j,i) * this%iw(i)
        end do
    end do

    ! And the split-form derivative matrix
    this%Ds      = 2.0_wp * this%D
    this%Ds(1,1) = this%Ds(1,1) + this%iw(1)
    this%Ds(n,n) = this%Ds(n,n) - this%iw(n)

    ! Variables for the telescopic fluxes. Refer to:
    ! Fisher, T. C., & Carpenter, M. H. (2013). High-order entropy stable
    ! finite difference schemes for nonlinear conservation laws: Finite domains.
    ! Journal of Computational Physics, 252, 518–557.

    ! Complementary grid...
    this%xc(0) = -1.0_wp
    this%xc(n) = 1.0_wp

    do i = 1, (n-1)/2
        this%xc(i)   = this%xc(i-1)   + this%w(i)
        this%xc(n-i) = this%xc(n+1-i) - this%w(n+1-i)
    end do

    if (mod(n-1, 2) /= 0) then
        this%xc(n/2) = this%xc(n/2-1) + this%w(n/2)
    end if

    ! ...and the corresponding interpolation coefficients
    allocate(avgTmp(n,n))
    do i = 1, n
        xn = this%xc(i-1) + (this%x / 2.0_wp + 0.5_wp) * this%w(i)
        call polynomial_interpolation_matrix(n, this%x, this%wb, &
                                             n, xn, avgTmp)
        this%avgInp(i,:) = matmul(this%w, avgTmp) / 2.0_wp
    end do

    ! Values of the Lagrange-Legendre transforms and their norms
    ! Loop over the polynomial orders
    do j = 1, n

        ! Loop over the quadrature nodes
        do i = 1, n

            ! Backwards transformation matrix...
            call legendre_polynomial(j-1, this%x(i), this%Bwd(i,j))

            ! ...and discrete norm of the $P_j$ Legendre polynomial
            this%modeN2(j) = this%modeN2(j) + this%Bwd(i,j)**2 * this%w(i)

            ! Numerator of the forward transformation matrix
            this%Fwd(j,i) = this%Bwd(i,j) * this%w(i)

        end do

        ! Division by the squared norm of the base element
        this%Fwd(j,:) = this%Fwd(j,:) / this%modeN2(j)

    end do

end subroutine nodalDG_GLL_constructor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Approximate the value at the boundaries from the inner nodes.
!
!> @param[in]   u   values inside the element
!> @param[out]  uL  values at the left boundary
!> @param[out]  uR  values at the right boundary
!···············································································
subroutine nodalDG_GLL_extrapolator(this, u, uL, uR)
    !* Arguments *!
    real(wp), intent(in)  :: u(:,:)
    real(wp), intent(out) :: uL(:)
    real(wp), intent(out) :: uR(:)
    ! Derived types
    class(NodalDG_GLL), intent(in) :: this

    uL = u(1,:)
    uR = u(this%n,:)

end subroutine nodalDG_GLL_extrapolator

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Overloaded assignment between 'NodalDG' objects.
!
!> @param[out]  lhs  NodalDG object to which the values are copied
!> @param[in]   rhs  NodalDG object from which the values are copied
!···············································································
subroutine nodalDG_GLL_assignment(lhs, rhs)
    !* Arguments *!
    class(NodalDG_GLL), intent(out) :: lhs
    class(NodalDG_GLL), intent(in)  :: rhs

    ! Base class
    lhs%n         = rhs%n
    lhs%nt        = rhs%nt
    lhs%hasBounds = rhs%hasBounds
    lhs%x         = rhs%x
    lhs%xc        = rhs%xc
    lhs%avgInp    = rhs%avgInp
    lhs%D         = rhs%D
    lhs%bdMode    = rhs%bdMode
    lhs%wb        = rhs%wb
    lhs%Dh        = rhs%Dh
    lhs%w         = rhs%w
    lhs%iw        = rhs%iw
    lhs%modeN2    = rhs%modeN2
    lhs%Fwd       = rhs%Fwd
    lhs%Bwd       = rhs%Bwd

    if (allocated(rhs%fineInp)) lhs%fineInp = rhs%fineInp
    if (allocated(rhs%inp))     lhs%inp = rhs%inp

    ! This class
    lhs%Q    = rhs%Q
    lhs%Ds   = rhs%Ds
    lhs%inpC = rhs%inpC

end subroutine nodalDG_GLL_assignment

end module NodalDG_GaussLobattoClass
