!*******************************************************************************
!  MODULE: StdExpansionClass
!
!> @author
!> Andres Mateo
!
!> @brief
!> Defines an abstract base class for the discretization of a 1D standard
!> element, and another one for the DGSEM.
!*******************************************************************************

module StdExpansionClass

    use Constants, only: wp, fNone, eLeft, eRight
    use ExceptionsAndMessages, only: printError
    use Legendre
    use Lagrange

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public types
    public :: StdExp_t

!···············································································
!> @class StdExp_t
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Base class for the discretization of a 1D element.
!···············································································
    type :: StdExp_t
        integer  :: n  = 0                    !< number of modes
        integer  :: nt = fNone                !< type of discretization
        logical  :: hasBounds                 !< nodes 1 and n are at the boundary
        real(wp), allocatable :: x(:)         !< nodes
        real(wp), allocatable :: xc(:)        !< nodes of the complementary grid
        real(wp), allocatable :: lh(:,:)      !< interpolation polynomials at the complementary grid
        real(wp), allocatable :: avgInp(:,:)  !< coefficients for avg. interpolation
        real(wp), allocatable :: D(:,:)       !< derivative matrix
        real(wp), allocatable :: bdMode(:,:)  !< interpolation polynomial at the bounds
        real(wp), allocatable :: wb(:)        !< barycentric weights
        real(wp), allocatable :: Dh(:,:)      !< \f$\hat{D}\f$ matrix
        real(wp), allocatable :: w(:)         !< quadrature weights
        real(wp), allocatable :: iw(:)        !< quadrature weights (inverse)
        real(wp), allocatable :: modeN2(:)    !< norm^2 of the interpolation polynomial
        real(wp), allocatable :: Fwd(:,:)     !< nodal to modal transform
        real(wp), allocatable :: Bwd(:,:)     !< modal to nodal transform
        real(wp), allocatable :: inp(:,:)     !< interpolation matrix
        real(wp), allocatable :: fineInp(:,:) !< interpolation matrix to fine mesh
        real(wp), allocatable :: Hsvv(:,:)    !< SVV filtering matrix H
    contains
        procedure :: StdExp_allocator
        procedure :: StdExp_assign
        procedure :: init    => StdExp_constructor
        procedure :: project => StdExp_projector
    end type StdExp_t

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Base class constructor, do not call directly (Algorithm 59).
!
!> @param[in]  n     number of nodes
!> @param[in]  Psvv  exponent of the SVV filtering kernel (optional)
!···············································································
subroutine StdExp_constructor(this, n, Psvv)
    !* Arguments *!
    integer,           intent(in) :: n
    integer, optional, intent(in) :: Psvv
    ! Derived types
    class(StdExp_t), intent(inout) :: this

    !* Local variables *!
    integer               :: i
    integer               :: j
    real(wp), allocatable :: xn(:)
    real(wp), allocatable :: avgTmp(:,:)
    real(wp), allocatable :: coeff(:)


    ! Nodes and weights of the quadrature
    call barycentric_weights(n, this%x, this%wb)
    this%iw = 1.0_wp / this%w

    ! Calculate the interpolating polynomials at the boundaries
    call lagrange_interpolating_polynomials(n, -1.0_wp, this%x, this%wb, &
                                            this%bdMode(eLeft,:))
    call lagrange_interpolating_polynomials(n, 1.0_wp, this%x, this%wb,  &
                                            this%bdMode(eRight,:))

    ! Calculate the derivative matrices
    call polynomial_derivative_matrix(n, this%x, this%D)
    do i = 1, n
        do j = 1, n
            this%Dh(i,j) = -this%D(j,i) * this%w(j) * this%iw(i)
        end do
    end do

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

    ! Interpolation coefficients for subelement averaging
    allocate(avgTmp(n,n))
    do i = 1, n
        xn = this%xc(i-1) + (this%x / 2.0_wp + 0.5_wp) * this%w(i)
        call polynomial_interpolation_matrix(n, this%x, this%wb, n, xn, avgTmp)
        this%avgInp(i,:) = matmul(this%w, avgTmp) / 2.0_wp
    end do

    ! Interpolation coefficients at the nodes of the complementary grid
    do i = 0, n
        call lagrange_interpolating_polynomials(n, this%xc(i), this%x, &
                                                this%wb, this%lh(:,i))
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

    ! SVV filtering coefficients and final H matrix: Bwd.diag(F).Fwd
    if (present(Psvv)) then

        coeff = [((real(i, kind=wp)/(n-1))**Psvv, i = 0, n-1)]

        do i = 1, n
            this%Hsvv(i,:) = coeff(i) * this%Fwd(i,:)
        end do
        this%Hsvv = matmul(this%Bwd, this%Hsvv)

    end if

end subroutine StdExp_constructor

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
subroutine StdExp_projector(this, u, uL, uR)
    !* Arguments *!
    real(wp), intent(in)  :: u(:,:)
    real(wp), intent(out) :: uL(:)
    real(wp), intent(out) :: uR(:)
    ! Derived types
    class(StdExp_t), intent(in) :: this


    ! Avoid warnings (at least with intel)
    uL = 0.0_wp
    uR = 0.0_wp

    call printError("StdExpansionClass.f90", &
                    "Base class projector is not implemented.")

end subroutine StdExp_projector

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Allocate space for allocatables of the base class. Do not call directly.
!
!> @param[in]  n          number of quadrature points.
!> @param[in]  allocHsvv  .true. if Hsvv must be allocated
!···············································································
subroutine StdExp_allocator(this, n, allocHsvv)
    !* Arguments *!
    integer, intent(in) :: n
    logical, intent(in) :: allocHsvv
    ! Derived types
    class(StdExp_t), intent(inout) :: this


    this%n = n

    ! 'inp' is not allocated here, it is set when needed by adaptation rutines
    ! 'fineRes' is allocated when constructing a 'Printer' object
    allocate(this%wb(n), source=0.0_wp)
    allocate(this%D(n, n), source=0.0_wp)
    allocate(this%Dh(n ,n), source=0.0_wp)
    allocate(this%bdMode(2, n), source=0.0_wp)
    allocate(this%x(n), source=0.0_wp)
    allocate(this%xc(0:n), source=0.0_wp)
    allocate(this%lh(n, 0:n), source=0.0_wp)
    allocate(this%avgInp(n, n), source=0.0_wp)
    allocate(this%w(n), source=0.0_wp)
    allocate(this%iw(n), source=0.0_wp)
    allocate(this%modeN2(n), source=0.0_wp)
    allocate(this%Fwd(n, n), source=0.0_wp)
    allocate(this%Bwd(n, n), source=0.0_wp)

    if (allocHsvv) then
        allocate(this%Hsvv(n, n), source=0.0_wp)
    end if

end subroutine StdExp_allocator

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Overloaded assignment between 'StdExp_t' objects. Do not call directly.
!
!> @param[out]  lhs  StdExp_t object to which the values are copied
!> @param[in]   rhs  StdExp_t object from which the values are copied
!···············································································
subroutine StdExp_assign(lhs, rhs)
    !* Arguments *!
    class(StdExp_t), intent(out) :: lhs
    class(StdExp_t), intent(in)  :: rhs

    lhs%n         = rhs%n
    lhs%nt        = rhs%nt
    lhs%hasBounds = rhs%hasBounds
    lhs%x         = rhs%x
    lhs%xc        = rhs%xc
    lhs%lh        = rhs%lh
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
    lhs%Hsvv      = rhs%Hsvv

    if (allocated(rhs%fineInp)) lhs%fineInp = rhs%fineInp
    if (allocated(rhs%inp))     lhs%inp     = rhs%inp

end subroutine StdExp_assign

end module StdExpansionClass
