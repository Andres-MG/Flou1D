!*******************************************************************************
!  MODULE: SodsTube
!
!> @author
!> Andres Mateo
!
!> @brief
!> Boundary and initial conditions for the 1D Sod's tube test case.
!*******************************************************************************

module SodsTube

    use Constants
    use Physics

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: SodsTubeIC

    ! Generic interfaces
    interface SodsTubeIC
        module procedure SodsTubeICscalar
        module procedure SodsTubeICvector
    end interface SodsTubeIC

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initial conditions for the Sod's tube test case (1D).
!
!> @param[in]   x      coordinates of the node
!> @param[out]  Phi    values at the node
!···············································································
subroutine SodsTubeICscalar(x, Phi)
    !* Arguments *!
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: Phi(:)

    if (x < 0.5_wp) then
        Phi(IRHO)  = 1.0_wp
        Phi(IRHOU) = 0.0_wp
        Phi(IRHOE) = 2.5_wp
   else
        Phi(IRHO)  = 0.125_wp
        Phi(IRHOU) = 0.0_wp
        Phi(IRHOE) = 0.25_wp
    end if

end subroutine SodsTubeICscalar

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initial conditions for the Sod's tube test case (1D).
!
!> @param[in]   x      coordinates of the nodes
!> @param[out]  Phi    values at the nodes
!···············································································
subroutine SodsTubeICvector(x, Phi)
    !* Arguments *!
    real(wp), intent(in)  :: x(:)
    real(wp), intent(out) :: Phi(:,:)

    where (x < 0.5_wp)
        Phi(:,IRHO)  = 1.0_wp
        Phi(:,IRHOU) = 0.0_wp
        Phi(:,IRHOE) = 2.5_wp
    elsewhere
        Phi(:,IRHO)  = 0.125_wp
        Phi(:,IRHOU) = 0.0_wp
        Phi(:,IRHOE) = 0.25_wp
    end where

end subroutine SodsTubeICvector

end module SodsTube
