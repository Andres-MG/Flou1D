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

    use Constants

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
    type, abstract :: StdExp_t
        integer  :: n  = 0                    !< number of modes
        integer  :: nt = fNone                !< type of discretization
        logical  :: hasBounds                 !< nodes 1 and n are at the boundary
        real(wp), allocatable :: x(:)         !< nodes
        real(wp), allocatable :: xc(:)        !< nodes of the complementary grid
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
    contains
        procedure(exp_constructor),  deferred :: construct
        procedure(exp_extrapolator), deferred :: extrapolate
    end type StdExp_t

    abstract interface
        ! Interface for the constructor subroutine
        subroutine exp_constructor(this, n)
            import StdExp_t
            integer,         intent(in)    :: n    !< number of modes
            class(StdExp_t), intent(inout) :: this
        end subroutine exp_constructor

        ! Interface for the extrapolator subroutine
        subroutine exp_extrapolator(this, u, uL, uR)
            use Constants, only: wp
            import StdExp_t
            real(wp), intent(in)  :: u(:,:)
            real(wp), intent(out) :: uL(:)
            real(wp), intent(out) :: uR(:)
            class(StdExp_t), intent(in) :: this
        end subroutine exp_extrapolator
    end interface

end module StdExpansionClass
