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

    use Constants, only: wp, eGaussLobatto
    use StdExpansionClass
    use Legendre

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
    contains
        procedure, private :: nodalDG_GLL_assignment
        procedure :: init          => nodalDG_GLL_constructor
        procedure :: project       => nodalDG_GLL_projector
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
!> @param[in]  Psvv  exponent of the SVV filtering kernel (optional)
!···············································································
subroutine nodalDG_GLL_constructor(this, n, Psvv)
    !* Arguments *!
    integer,           intent(in) :: n
    integer, optional, intent(in) :: Psvv
    ! Derived types
    class(NodalDG_GLL), intent(inout) :: this

    !* Local variables *!
    integer               :: i
    integer               :: j

    ! Set implementation-specific attributes
    this%hasBounds = .true.
    this%nt        = eGaussLobatto

    ! Allocate dynamic storage
    if (present(Psvv)) then
        call this%StdExp_t%StdExp_allocator(n, allocHsvv=.true.)
    else
        call this%StdExp_t%StdExp_allocator(n, allocHsvv=.false.)
    end if
    allocate(this%Q(n, n), source=0.0_wp)
    allocate(this%Ds(n ,n), source=0.0_wp)

    ! Nodes and weights of the quadrature
    call legendre_gauss_lobatto_nodes_weights(n, this%x, this%w)

    ! Initialize the base attributes
    if (present(Psvv)) then
        call this%StdExp_t%init(n, Psvv)
    else
        call this%StdExp_t%init(n)
    end if

    ! SBP matrix so that: MD = Q
    do i = 1, n
        do j = 1, n
            this%Q(j,i)  = this%D(j,i) * this%w(j)
        end do
    end do

    ! Split-form derivative matrix
    this%Ds      = 2.0_wp * this%D
    this%Ds(1,1) = this%Ds(1,1) + this%iw(1)
    this%Ds(n,n) = this%Ds(n,n) - this%iw(n)

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
subroutine nodalDG_GLL_projector(this, u, uL, uR)
    !* Arguments *!
    real(wp), intent(in)  :: u(:,:)
    real(wp), intent(out) :: uL(:)
    real(wp), intent(out) :: uR(:)
    ! Derived types
    class(NodalDG_GLL), intent(in) :: this

    uL = u(1,:)
    uR = u(this%n,:)

end subroutine nodalDG_GLL_projector

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
    call lhs%StdExp_t%StdExp_assign(rhs)

    ! This class
    lhs%Q    = rhs%Q
    lhs%Ds   = rhs%Ds

end subroutine nodalDG_GLL_assignment

end module NodalDG_GaussLobattoClass
