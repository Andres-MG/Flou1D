!*******************************************************************************
!  MODULE: NormalShock
!
!> @author
!> Andres Mateo
!
!> @brief
!> Boundary and initial conditions for a 1D test case with a normal shock wave.
!*******************************************************************************

module NormalShock

    use Constants
    use Physics

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: afterShockValues

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Calculates the value of the conservative variables after a normal shock wave
!> given the conservative values of the incoming flow.
!
!> @param[in]   gamma   gamma constant
!> @param[in]   PhiIn   inflow conservative values
!> @param[out]  PhiOut  outflow conservative values
!···············································································
subroutine afterShockValues(gamma, PhiIn, PhiOut)
    !* Arguments *!
    real(wp), intent(in)  :: gamma
    real(wp), intent(in)  :: PhiIn(:)
    real(wp), intent(out) :: PhiOut(:)

    !* Local variables *!
    real(wp) :: p
    real(wp) :: a
    real(wp) :: M1
    real(wp) :: M2
    real(wp) :: pR
    real(wp) :: rR
    real(wp) :: TR
    real(wp) :: p0R

    ! Initialisation
    p  = (gamma-1.0_wp) * (PhiIn(IRHOE) - 0.5_wp*PhiIn(IRHOU)**2/PhiIn(IRHO))
    a  = sqrt( gamma * p / PhiIn(IRHO) )
    M1 = PhiIn(IRHOU)/PhiIn(IRHO) / a

    ! Get variables after the shock
    call NormalShockValues(M1, gamma, M2, pR, rR, TR, p0R)

    ! Switch back to conservative variables
    ! Rho
    PhiOut(IRHO) = rR * PhiIn(IRHO)

    ! Rho * u
    p = pR * p
    a = sqrt( gamma * p / PhiOut(IRHO) )
    PhiOut(IRHOU) = PhiOut(IRHO) * M2 * a

    ! Rho * e
    PhiOut(IRHOE) = p / (gamma-1.0_wp) + 0.5_wp*PhiOut(IRHOU)**2/PhiOut(IRHO)

end subroutine afterShockValues

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Computes different ratios between thermodynamical variables before and after
!> a normal shock wave.
!
!> @param[in]   M1    incoming Mach number of the flow
!> @param[in]   gam   gamma constant
!> @param[out]  M2    outgoing Mach number
!> @param[out]  pR    ratio of pressures (\$p_2/p_1\$)
!> @param[out]  rhoR  ratio of densities (\$\rho_2/\rho_1\$)
!> @param[out]  TR    ratio of temperatures (\$T_2/T_1\$)
!> @param[out]  p0R   ratio of stagnation pressures (\$p_{0,2}/p_{0,1}\$)
!···············································································
subroutine NormalShockValues(M1, gam, M2, pR, rhoR, TR, p0R)
    !* Arguments *!
    real(wp), intent(in)  :: M1
    real(wp), intent(in)  :: gam
    real(wp), intent(out) :: M2
    real(wp), intent(out) :: pR
    real(wp), intent(out) :: rhoR
    real(wp), intent(out) :: TR
    real(wp), intent(out) :: p0R

    !* Local variables *!
    real(wp) :: gm1
    real(wp) :: gp1
    real(wp) :: gm12
    real(wp) :: gp12

    ! Initialisation
    gm1  = gam - 1.0_wp
    gp1  = gam + 1.0_wp
    gm12 = 0.5_wp * gm1
    gp12 = 0.5_wp * gp1

    ! Normal shock equations
    M2   = sqrt( (M1**2 * gm1 + 2.0_wp) / (2.0_wp * gam * M1**2 - gm1) )
    pR   = gam / gp12 * M1**2 - gm1 / gp1
    rhoR = gp1 * M1**2 / (gm1 * M1**2 + 2.0_wp)
    TR   = (1.0_wp+gm12*M1**2)*(gam/gm12*M1**2-1.0_wp) / (M1**2*(gam/gm12+gm12))
    p0R  = (gp12*M1**2 / (1.0_wp + gm12*M1**2)) ** (gam/gm1) * &
           (1.0_wp / (gam/gp12*M1**2 - gm1/gp1)) ** (1.0_wp/gm1)

end subroutine NormalShockValues

end module NormalShock
