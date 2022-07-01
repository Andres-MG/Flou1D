module Physics_WE

    use Constants

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public variables
    integer, parameter, public :: NEQS = 2

    real(wp), public, protected :: pC

    ! Explicitly define public functions
    public :: initPhysics
    public :: flux

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initialises the especific values of the compressible Euler equations
!···············································································
subroutine initPhysics(c)
    !* Arguments *!
    real(wp), intent(in) :: c

    pC = c

end subroutine initPhysics

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 121
!···············································································
function flux(Phi)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:)

    !* Return values *!
    real(wp) :: flux(NEQS)

    flux(1) = pC**2 * Phi(2)
    flux(2) = Phi(1)

end function flux

end module Physics_WE
