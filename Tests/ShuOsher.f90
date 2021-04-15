!*******************************************************************************
!  MODULE: ShuOsher
!
!> @author
!> Andres Mateo
!
!> @brief
!> Boundary and initial conditions for the 1D Shu-Osher test case.
!*******************************************************************************

module ShuOsher

    use Constants
    use Physics

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: ShuOsherIC

    ! Generic interfaces
    interface ShuOsherIC
        module procedure ShuOsherICscalar
        module procedure ShuOsherICvector
    end interface ShuOsherIC

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initial conditions for the Shu-Osher test case (1D).
!
!> @param[in]   x      coordinates of the node
!> @param[in]   gamma  gamma constant
!> @param[out]  Phi    values at the node
!···············································································
subroutine ShuOsherICscalar(x, gamma, Phi)
    !* Arguments *!
    real(wp), intent(in)  :: x
    real(wp), intent(in)  :: gamma
    real(wp), intent(out) :: Phi(:)

    !* Local variables *!
    real(wp) :: invGammaOne

    ! Set up
    invGammaOne = 1.0_wp / (gamma - 1.0_wp)

    if (x < 0.1_wp) then
        Phi(1) = 3.857143_wp
        Phi(2) = 3.857143_wp * 2.629369_wp
        Phi(3) = 10.33333_wp*invGammaOne + 3.857143_wp*2.629369_wp**2 * 0.5_wp
    else
        Phi(1) = 1.0_wp + 0.2_wp * sin(50.0_wp * x)
        Phi(2) = 0.0_wp
        Phi(3) = invGammaOne
    end if

end subroutine ShuOsherICscalar

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initial conditions for the Shu-Osher test case (1D).
!
!> @param[in]   x      coordinates of the nodes
!> @param[in]   gamma  gamma constant
!> @param[out]  Phi    values at the nodes
!···············································································
subroutine ShuOsherICvector(x, gamma, Phi)
    !* Arguments *!
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: gamma
    real(wp), intent(out) :: Phi(:,:)

    !* Local variables *!
    real(wp) :: invGammaOne

    ! Set up
    invGammaOne = 1.0_wp / (gamma - 1.0_wp)

    where (x < 0.1_wp)
        Phi(:,1) = 3.857143_wp
        Phi(:,2) = 3.857143_wp * 2.629369_wp
        Phi(:,3) = 10.33333_wp*invGammaOne + 3.857143_wp*2.629369_wp**2 * 0.5_wp
    elsewhere
        Phi(:,1) = 1.0_wp + 0.2_wp * sin(50.0_wp * x)
        Phi(:,2) = 0.0_wp
        Phi(:,3) = invGammaOne
    end where

end subroutine ShuOsherICvector

end module ShuOsher
