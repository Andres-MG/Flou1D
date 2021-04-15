!*******************************************************************************
!  MODULE: ArtificialViscosity_CE
!
!> @author
!> Andres Mateo
!
!> @brief
!> Definition of the implemented artificial viscosity formulations.
!*******************************************************************************

module ArtificialViscosity_NS

    use Constants
    use Utilities
    use ExceptionsAndMessages
    use PhysicsStorage_NS

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: initArtificialViscosity
    public :: ArtViscousFlux
    public :: destructArtificialViscosity

    abstract interface
        function ArtViscFlux_Int(Phi, Grad, alpha)
            import wp
            import NEQS
            real(wp), intent(in) :: Phi(:)
            real(wp), intent(in) :: Grad(:)
            real(wp), intent(in) :: alpha(:)
            real(wp) :: ArtViscFlux_Int(NEQS)
        end function ArtViscFlux_Int
    end interface

    procedure(ArtViscFlux_Int), pointer :: ArtViscousFlux => null()

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initialises the function pointer to the selected viscous flux implementation.
!
!> @param[in]  artViscType  type of artificial viscosity
!···············································································
subroutine initArtificialViscosity(artViscType)
    !* Arguments *!
    integer,  intent(in) :: artViscType

    ! Select artificial viscosity flux fomrulation
    select case (artViscType)

    case (eLaplacianVisc)
        ArtViscousFlux => LaplacianViscosity
        ! Phys%WithEntropyVars = .false.
        Phys%IsViscous = .true.

    case (eGuermondPhysical)
        ArtViscousFlux => GuermondPhysical
        ! Phys%WithEntropyVars = .false.
        Phys%IsViscous = .true.

    case (eGuermondEntropy)
        ArtViscousFlux => GuermondEntropy
        Phys%WithEntropyVars = .true.
        Phys%IsViscous = .true.

    case (fNone)
        ArtViscousFlux => noArtificialViscosity
        ! Phys%WithEntropyVars = .false.
        ! Phys%IsViscous = .false.

    case default
        call printError("Physics_CE.f90", &
                        "The selected artificial viscosity is not available.")

    end select

end subroutine initArtificialViscosity

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> 'Empty' artificial viscosity flux. Use only for inviscid cases.
!>  \f[ Phi_1 \quad \rho   \f]
!>  \f[ Phi_2 \quad \rho u \f]
!>  \f[ Phi_3 \quad \rho e \f]
!
!> @param[in]  Phi     values of the conservative variables
!> @param[in]  Grad    gradients of the conservative variables
!> @param[in]  alphas  artificial viscosity constants
!
!> @return            zero viscous flux
!···············································································
function noArtificialViscosity(Phi, Grad, alphas)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:)
    real(wp), intent(in) :: Grad(:)
    real(wp), intent(in) :: alphas(:)

    !* Return values *!
    real(wp) :: noArtificialViscosity(NEQS)

    noArtificialViscosity = 0.0_wp

end function noArtificialViscosity

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Viscous flux defined as the gradient of the conservative variables times a
!> numerical constant 'pAlpha'.
!>  \f[ Phi_1 \quad \rho   \f]
!>  \f[ Phi_2 \quad \rho u \f]
!>  \f[ Phi_3 \quad \rho e \f]
!
!> @param[in]  Phi     values of the conservative variables (not used)
!> @param[in]  Grad    gradients of the conservative variables
!> @param[in]  alphas  artificial viscosity constants
!
!> @return             viscous flux
!···············································································
function LaplacianViscosity(Phi, Grad, alphas)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:)
    real(wp), intent(in) :: Grad(:)
    real(wp), intent(in) :: alphas(:)

    !* Return values *!
    real(wp) :: LaplacianViscosity(NEQS)

    ! Viscous fluxes
    LaplacianViscosity = alphas(1) * Grad

end function LaplacianViscosity

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Density regularisation based on the physical variables \f$ (\rho,u,e) \f$,
!> where a \f$ \lambda \f$ term is added to the reference to get a Navier-
!> Stokes-like viscous flux. Refer to Guermond, J. L., & Popov, B. (2014).
!> Viscous regularization of the Euler equations and entropy principles. SIAM
!> Journal on Applied Mathematics, 74(2), 284–305.
!
!>  \f[ Phi_1  \quad \rho   \f]
!>  \f[ Phi_2  \quad \rho u \f]
!>  \f[ Phi_3  \quad \rho e \f]
!>  \f[ Grad_1 \quad \nabla \rho   \f]
!>  \f[ Grad_2 \quad \nabla (\rho u) \f]
!>  \f[ Grad_3 \quad \nabla (\rho e) \f]
!
!> @param[in]  Phi     values of the conservative variables
!> @param[in]  Grad    gradients of the conservative variables
!> @param[in]  alphas  artificial viscosity constants
!
!> @return             viscous flux
!···············································································
function GuermondPhysical(Phi, Grad, alphas)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:)
    real(wp), intent(in) :: Grad(:)
    real(wp), intent(in) :: alphas(:)

    !* Return values *!
    real(wp) :: GuermondPhysical(NEQS)

    !* Local variables *!
    real(wp) :: u
    real(wp) :: e
    real(wp) :: ux
    real(wp) :: eix
    real(wp) :: alpha
    real(wp) :: beta
    real(wp) :: lambda

    ! Intermediate values
    u   = Phi(IRHOU) / Phi(IRHO)
    e   = Phi(IRHOE) / Phi(IRHO)
    ux  = ( Grad(IRHOU) - Grad(IRHO)*u ) / Phi(IRHO)
    eix = ( Grad(IRHOE) - Grad(IRHO)*e ) / Phi(IRHO) - u * ux

    ! Density regularisation
    GuermondPhysical(IRHO)  = Grad(IRHO)
    GuermondPhysical(IRHOU) = u * Grad(IRHO)
    GuermondPhysical(IRHOE) = Grad(IRHOE) - Phi(IRHO) * u * ux

    ! Viscous term
    alpha  = alphas(1)
    beta   = alphas(2)
    lambda = alphas(3)
    GuermondPhysical = alpha * GuermondPhysical          &
                     + beta * ux * [ 0.0_wp, 1.0_wp, u ] &
                     + lambda * [0.0_wp, 0.0_wp, eix ]

end function GuermondPhysical

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Density regularisation based on the entropy variables \f$ (w_1,w_2,w_3) \f$,
!> where a \f$ \lambda \f$ term is added to the reference to get a Navier-
!> Stokes-like viscous flux. Refer to Guermond, J. L., & Popov, B. (2014).
!> Viscous regularization of the Euler equations and entropy principles. SIAM
!> Journal on Applied Mathematics, 74(2), 284–305.
!
!>  \f[ Phi_1  \quad \rho   \f]
!>  \f[ Phi_2  \quad \rho u \f]
!>  \f[ Phi_3  \quad \rho e \f]
!>  \f[ Grad_1 \quad \nabla w_1 \f]
!>  \f[ Grad_2 \quad \nabla w_2 \f]
!>  \f[ Grad_3 \quad \nabla w_3 \f]
!
!> @param[in]  Phi     values of the conservative variables
!> @param[in]  Grad    gradients of the entropy variables
!> @param[in]  alphas  artificial viscosity constants
!
!> @retur              viscous flux
!···············································································
function GuermondEntropy(Phi, Grad, alphas)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:)
    real(wp), intent(in) :: Grad(:)
    real(wp), intent(in) :: alphas(:)

    !* Return values *!
    real(wp) :: GuermondEntropy(NEQS)

    !* Local variables *!
    real(wp) :: alpha
    real(wp) :: beta
    real(wp) :: lambda
    real(wp) :: rho
    real(wp) :: u
    real(wp) :: e
    real(wp) :: p
    real(wp) :: lam2
    real(wp) :: tmp(NEQS,NEQS)
    real(wp) :: Bk(NEQS,NEQS)
    real(wp) :: Bm(NEQS,NEQS)
    real(wp) :: Bl(NEQS,NEQS)

    ! Intermediate variables
    rho   = Phi(IRHO)
    u     = Phi(IRHOU) / rho
    e     = Phi(IRHOE) / rho
    p     = Phys%GammaMinusOne * ( Phi(IRHOE) - rho * u**2 * 0.5_wp )
    lam2  = ( p / rho )**2 / Phys%GammaMinusOne

    ! Density regularisation matrix
    Bk(1,:) = [ 1.0_wp, u,    e           ]
    Bk(2,:) = [ u,      u**2, u*e         ]
    Bk(3,:) = [ e,      u*e,  e**2 + lam2 ]

    ! Viscous matrix (momentum)
    Bm(1,:) = [ 0.0_wp, 0.0_wp, 0.0_wp  ]
    Bm(2,:) = [ 0.0_wp, 1.0_wp, u       ]
    Bm(3,:) = [ 0.0_wp, u,      u**2    ]

    ! Viscous matrix (energy)
    Bl(1,:) = [ 0.0_wp, 0.0_wp, 0.0_wp ]
    Bl(2,:) = [ 0.0_wp, 0.0_wp, 0.0_wp ]
    Bl(3,:) = [ 0.0_wp, 0.0_wp, lam2   ]

    ! Guermond flux
    alpha  = alphas(1)
    beta   = alphas(2)
    lambda = alphas(3)

    ! Avoid compilation error with Intel Fortran
    tmp = alpha * rho * Bk + beta * p / rho * Bm + lambda * Bl
    GuermondEntropy = matmul(tmp, Grad)

end function GuermondEntropy

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Destructs the artificial viscous flux function pointer.
!···············································································
subroutine destructArtificialViscosity()

    if (associated(artViscousFlux)) then
        nullify(artViscousFlux)
    end if

end subroutine destructArtificialViscosity

end module ArtificialViscosity_NS
