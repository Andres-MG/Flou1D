!*******************************************************************************
!  MODULE: NodalDG_GaussClass
!
!> @author
!> Andres Mateo
!
!> @brief
!> Defines the NodalDG_GL class.
!*******************************************************************************

module NodalDG_GaussClass

    use Constants, only: wp, eGauss, eLeft, eRight
    use StdExpansionClass
    use Legendre

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public types
    public :: NodalDG_GL

!···············································································
!> @class NodalDG_GL
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 58.
!···············································································
    type, extends(Stdexp_t) :: NodalDG_GL
    contains
        procedure, private :: nodalDG_GL_assignment
        procedure :: init          => nodalDG_GL_constructor
        procedure :: project       => nodalDG_GL_projector
        generic   :: assignment(=) => nodalDG_GL_assignment
    end type NodalDG_GL

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
subroutine nodalDG_GL_constructor(this, n, Psvv)
    !* Arguments *!
    integer,           intent(in) :: n
    integer, optional, intent(in) :: Psvv
    ! Derived types
    class(NodalDG_GL), intent(inout) :: this


    ! Set implementation-specific attributes
    this%hasBounds = .false.
    this%nt        = eGauss

    ! Allocate dynamic storage
    if (present(Psvv)) then
        call this%StdExp_t%StdExp_allocator(n, allocHsvv=.true.)
    else
        call this%StdExp_t%StdExp_allocator(n, allocHsvv=.false.)
    end if

    ! Nodes and weights of the quadrature
    call legendre_gauss_nodes_weights(n, this%x, this%w)

    ! And initialize the base attributes
    if (present(Psvv)) then
        call this%StdExp_t%init(n, Psvv)
    else
        call this%StdExp_t%init(n)
    end if

end subroutine nodalDG_GL_constructor

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
subroutine nodalDG_GL_projector(this, u, uL, uR)
    !* Arguments *!
    real(wp), intent(in)  :: u(:,:)
    real(wp), intent(out) :: uL(:)
    real(wp), intent(out) :: uR(:)
    ! Derived types
    class(NodalDG_GL), intent(in) :: this

    uL = matmul(this%bdMode(eLeft,:), u)
    uR = matmul(this%bdMode(eRight,:), u)

end subroutine nodalDG_GL_projector

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
subroutine nodalDG_GL_assignment(lhs, rhs)
    !* Arguments *!
    class(NodalDG_GL), intent(out) :: lhs
    class(NodalDG_GL), intent(in)  :: rhs


    ! Base class assignment
    call lhs%StdExp_t%StdExp_assign(rhs)

end subroutine nodalDG_GL_assignment

end module NodalDG_GaussClass
