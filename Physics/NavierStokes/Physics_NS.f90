!*******************************************************************************
!  MODULE: Physics_NS
!
!> @author
!> Andres Mateo
!
!> @brief
!> Defines some functions related to the NS equations.
!*******************************************************************************

module Physics_NS

    use Constants
    use PhysicsStorage_NS
    use ArtificialViscosity_NS

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: getPressure
    public :: getPhysicalEntropy
    public :: EulerFlux
    public :: ViscousFlux
    public :: computeMaxEigenvalue
    public :: getEntropyVars

    interface getPressure
        module procedure getPressure_scalar
        module procedure getPressure_vector
    end interface getPressure

    interface getPhysicalEntropy
        module procedure getPhysicalEntropy_scalar
        module procedure getPhysicalEntropy_vector
    end interface getPhysicalEntropy

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Calculates the pressure from the physical variables \f$\Phi\f$.
!
!> @param[in]  Phi  variables of the compressible 1D flow
!
!> @return          value of the pressure
!···············································································
function getPressure_scalar(Phi) result(p)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:)

    !* Return values *!
    real(wp) :: p

    !* Local variables *!
    real(wp) :: u

    u = Phi(IRHOU) / Phi(IRHO)
    p = Phys%GammaMinusOne * ( Phi(IRHOE) - 0.5_wp * Phi(IRHOU) * u )

end function getPressure_scalar

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Calculates the pressure from the physical variables \f$\Phi\f$.
!
!> @param[in]  Phi  variables of the compressible 1D flow
!
!> @return          value of the pressure
!···············································································
function getPressure_vector(Phi) result(p)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:,:)

    !* Return values *!
    real(wp) :: p(size(Phi, dim=1))

    !* Local variables *!
    real(wp) :: u(size(Phi, dim=1))

    u = Phi(:,IRHOU) / Phi(:,IRHO)
    p = Phys%GammaMinusOne * ( Phi(:,IRHOE) - 0.5_wp * Phi(:,IRHOU) * u )

end function getPressure_vector

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Calculates the physical entropy from the physical variables \f$\Phi\f$.
!
!> @param[in]  Phi  variables of the compressible 1D flow
!
!> @return          value of the entropy
!···············································································
function getPhysicalEntropy_scalar(Phi) result(s)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:)

    !* Return values *!
    real(wp) :: s

    !* Local variables *!
    real(wp) :: p

    p = getPressure(Phi)
    s = log(p) - Phys%Gamma * log(Phi(IRHO))

end function getPhysicalEntropy_scalar

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Calculates the physical entropy from the physical variables \f$\Phi\f$.
!
!> @param[in]  Phi  variables of the compressible 1D flow
!
!> @return          value of the entropy
!···············································································
function getPhysicalEntropy_vector(Phi) result(s)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:,:)

    !* Return values *!
    real(wp) :: s(size(Phi, dim=1))

    !* Local variables *!
    real(wp) :: p(size(Phi, dim=1))

    p = getPressure(Phi)
    s = log(p) - Phys%Gamma * log(Phi(:,IRHO))

end function getPhysicalEntropy_vector

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Calculates some useful variables for a compressible flow
!
!> @param[in]   Phi   variables of the compressible 1D flow
!> @param[out]  vars  entropy variables
!···············································································
subroutine getEntropyVars(Phi, vars)
    !* Arguments *!
    real(wp), intent(in)  :: Phi(:)
    real(wp), intent(out) :: vars(:)

    !* Local variables *!
    real(wp) :: p
    real(wp) :: invRT
    real(wp) :: s
    real(wp) :: u

    ! Intermediate thermodynamic variables
    u     = Phi(IRHOU) / Phi(IRHO)
    p     = getPressure(Phi)
    invRT = Phi(IRHO) / p
    s     = getPhysicalEntropy(Phi)

    ! Final values
    vars(IRHO)  = (Phys%Gamma - s) / Phys%GammaMinusOne - 0.5_wp * invRT * u**2
    vars(IRHOU) = u * invRT
    vars(IRHOE) = -invRT

end subroutine getEntropyVars

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Flux of the compressible Euler equation:
!>  \f[ Phi_1 \quad \rho   \f]
!>  \f[ Phi_2 \quad \rho u \f]
!>  \f[ Phi_3 \quad \rho e \f]
!
!> @param[in]  Phi  values of the conservative variables
!
!> @return          flux of the Euler equation
!···············································································
function EulerFlux(Phi)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:)

    !* Return values *!
    real(wp) :: EulerFlux(NEQS)

    !* Local variables *!
    real(wp) :: p

    ! Calculate pressure
    p = getPressure(Phi)

    EulerFlux(IRHO)  = Phi(IRHOU)
    EulerFlux(IRHOU) = p + Phi(IRHOU)**2 / Phi(IRHO)
    EulerFlux(IRHOE) = (Phi(IRHOE) + p) * Phi(IRHOU) / Phi(IRHO)

end function EulerFlux

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Viscous flux correspoinding to the Navier-Stokes equations.
!
!> @param[in]  Phi    values of the conservative variables
!> @param[in]  Grad   gradients of the conservative/entropy variables
!> @param[in]  alpha  artificial viscosity constant (optional)
!
!> @return            viscous flux
!···············································································
function ViscousFlux(Phi, Grad, alpha)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:)
    real(wp), intent(in) :: Grad(:)
    real(wp), optional, intent(in) :: alpha

    !* Return values *!
    real(wp) :: ViscousFlux(NEQS)

    !* Local variables *!
    real(wp) :: rho
    real(wp) :: u
    real(wp) :: e
    real(wp) :: p
    real(wp) :: ux
    real(wp) :: eix
    real(wp) :: gam
    real(wp) :: Bm(NEQS,NEQS)


    if (Phys%Mu /= 0) then

        if (Phys%WithEntropyVars) then

            ! Intermediate variables
            rho  = Phi(IRHO)
            u    = Phi(IRHOU) / rho
            p    = Phys%GammaMinusOne * ( Phi(IRHOE) - rho * u**2 * 0.5_wp )
            gam  = Phys%Gamma/Phys%GammaMinusOne / Phys%Pr * p / rho

            ! Viscous matrix
            Bm(1,:) = [ 0.0_wp, 0.0_wp,          0.0_wp                ]
            Bm(2,:) = [ 0.0_wp, 4.0_wp/3.0_wp,   4.0_wp/3.0_wp*u          ]
            Bm(3,:) = [ 0.0_wp, 4.0_wp/3.0_wp*u, 4.0_wp/3.0_wp*u**2 + gam ]

            ViscousFlux = Phys%Mu * p/rho * matmul( Bm, Grad )

        else

            ! Intermediate values
            u = Phi(IRHOU) / Phi(IRHO)
            e = Phi(IRHOE) / Phi(IRHO)

            ! Calculate derivatives
            ux  = ( Grad(IRHOU) - Grad(IRHO)*u ) / Phi(IRHO)
            eix = ( Grad(IRHOE) - Grad(IRHO)*e ) / Phi(IRHO) - u * ux

            ! Viscous fluxes
            ViscousFlux(IRHO)  = 0.0_wp
            ViscousFlux(IRHOU) = 4.0_wp/3.0_wp * ux
            ViscousFlux(IRHOE) = 4.0_wp/3.0_wp * u * ux &
                                      + Phys%Gamma/Phys%Pr * eix

            ViscousFlux = ViscousFlux * Phys%Mu

        end if

    else

        ViscousFlux = 0.0_wp

    end if

    ! Add artificial viscosity
    if (present(alpha)) then
    if (alpha /= 0.0_wp) then

        ViscousFlux = ViscousFlux + ArtViscousFlux(Phi, Grad, alpha &
                                  * [1.0_wp, Phys%Alpha(2), Phys%Alpha(3)])

    end if
    end if

end function ViscousFlux

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Calculates the eigenvalue of the linearized NS equations (1D) with the
!> highest absolute value.
!
!> @param[in]  Phi  values of the conservative variables
!
!> @return          absolute value of the maximum eigenvalue (1D)
!···············································································
function computeMaxEigenvalue(Phi) result(lambda)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:)

    !* Return values *!
    real(wp) :: lambda

    !* Local variables *!
    real(wp) :: u
    real(wp) :: p
    real(wp) :: a

    u = Phi(IRHOU) / Phi(IRHO)
    p = getPressure(Phi)
    a = sqrt(Phys%Gamma * p / Phi(IRHO))

    ! Find the greater lambda number
    lambda = a + abs(u)

end function computeMaxEigenvalue

end module Physics_NS
