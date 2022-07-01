module RiemannSolvers_WE

    use Constants
    use Physics_WE

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: setRiemannSolver
    public :: RiemannSolver
    public :: destructRiemannSolver

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initialises the Riemann solver selected by the user
!···············································································
subroutine setRiemannSolver()

    ! Empty for the 1D wave equation

end subroutine setRiemannSolver

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 88: modified for a 1D case
!
!> @param[in]   PhiL  values on the left of the face
!> @param[in]   PhiR  values on the right of the face
!> @param[in]   nHat  unit vector pointing outwards orthogonally
!> @param[out]  F     Flux at the interface
!···············································································
subroutine RiemannSolver(PhiL, PhiR, nHat, F)
    !* Arguments *!
    real(wp), intent(in)  :: PhiL(:)
    real(wp), intent(in)  :: PhiR(:)
    real(wp), intent(in)  :: nHat
    real(wp), intent(out) :: F(:)

    !* Local variables *!
    real(wp) :: wr  ! Wave moving to the right
    real(wp) :: wl  ! Wave moving to the left

    ! Intermediate values
    wr = PhiL(1) + pC * PhiL(2)
    wl = PhiR(1) - pC * PhiR(2)

    ! Riemann fluxes
    F(1) = pC * ( wr - wl ) * 0.5_wp
    F(2) =      ( wr + wl ) * 0.5_wp
    F    = F * nHat

end subroutine RiemannSolver

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Releases the used memory
!···············································································
subroutine destructRiemannSolver()

    ! Empty subroutine

end subroutine destructRiemannSolver

end module RiemannSolvers_WE
