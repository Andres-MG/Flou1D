!*******************************************************************************
!  MODULE: RiemannSolvers_NS
!
!> @author
!> Andres Mateo
!
!> @brief
!> Implementation of different Riemann Solvers for the NS equations.
!*******************************************************************************

module RiemannSolvers_NS

    use Constants
    use Utilities
    use Physics_NS
    use PhysicsStorage_NS
    use ExceptionsAndMessages

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: initRiemannSolver
    public :: TwoPointFlux
    public :: DissipativeFlux
    public :: RiemannSolver
    public :: destructRiemannSolver

    abstract interface
        function Flux_Int(PhiL, PhiR)
            import wp
            import NEQS
            real(wp), intent(in) :: PhiL(:)
            real(wp), intent(in) :: PhiR(:)
            real(wp) :: Flux_Int(NEQS)
        end function Flux_Int
    end interface

    ! Decomposition of the numerical flux
    procedure(Flux_Int), pointer :: TwoPointFlux    => null()
    procedure(Flux_Int), pointer :: DissipativeFlux => null()

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initialises the Riemann solver selected by the user. 'avgType' and
!> 'dissipationType' set the different algorithms for the two terms of the
!> Lax-Friedrichs numerical flux:
!>  \f[ f^*(u_L,u_R,\hat{n}) = \hat{n} \cdot (avg(u_L,u_R) - diss(u_L,u_R)) /f]
!
!> @param[in]  avgType          averaging algorithm
!> @param[in]  dissipationType  dissipation method
!···············································································
subroutine initRiemannSolver(avgType, dissipationType)
    !* Arguments *!
    integer, intent(in) :: avgType
    integer, intent(in) :: dissipationType

    ! Set the averaging algorithm
    select case (avgType)

    case (eStdAvg)
        TwoPointFlux => standardAvgFlux

    case (eDucros)
        TwoPointFlux => DucrosFlux

    case (eKennedyGruber)
        TwoPointFlux => KennedyGruberFlux

    case (ePirozzoli)
        TwoPointFlux => PirozzoliFlux

    case (eIsmailRoe)
        TwoPointFlux => IsmailRoeFlux

    case (eChandrashekar)
        TwoPointFlux => ChandrashekarFlux

    case (fNone)
        TwoPointFlux => nullFlux

    case default
        call printError("RiemannSolvers_NS.f90", &
                        "The selected averaging algorithm is not available.")

    end select

    ! Set the dissipation algorithm
    select case (dissipationType)

    case (eLocalLxF)
        DissipativeFlux => localAvgDissipation

    case (eMatrixDiss)
        DissipativeFlux => matrixDissipation

    case (eIRdiss)
        DissipativeFlux => IsmailRoeDissipation

    case (eCHdiss)
        DissipativeFlux => ChandrashekarDissipation

    case (fNone)
        DissipativeFlux => nullFlux

    case default
        call printError("RiemannSolvers_NS.f90", &
                        "The selected dissipation algorithm is not available.")

    end select

end subroutine initRiemannSolver

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Numerical flux
!>  \f[ f^*(u_L,u_R,\hat{n}) = \hat{n} \cdot (avg(u_L,u_R) - diss(u_L,u_R)) /f]
!
!> @param[in]   PhiL  values at the left of the boundary
!> @param[in]   PhiR  values at the right of the boundary
!> @param[in]   nHat  direction of the normal (pointing outwards)
!
!> @return            numerical flux at the boundary
!···············································································
function RiemannSolver(PhiL, PhiR, nHat) result(F)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(:)
    real(wp), intent(in) :: PhiR(:)
    real(wp), intent(in) :: nHat

    !* Return values *!
    real(wp) :: F(NEQS)

    !* Local variables *!
    real(wp) :: TPflux(NEQS)
    real(wp) :: DissFlux(NEQS)

    TPflux = TwoPointFlux(PhiL, PhiR)
    DissFlux = DissipativeFlux(PhiL, PhiR)
    F = nHat * ( TPflux - DissFlux )

end function RiemannSolver

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Returns an empty vector.
!
!> @param[in]  PhiL  values at the left of the boundary
!> @param[in]  PhiR  values at the right of the boundary
!
!> @return           numerical flux at the boundary
!···············································································
function nullFlux(PhiL, PhiR) result(F)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(:)
    real(wp), intent(in) :: PhiR(:)

    !* Return values *!
    real(wp) :: F(NEQS)

    F = 0.0_wp

end function nullFlux

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Symmetric part of the numerical flux as standard average
!
!> @param[in]  PhiL  values at the left of the boundary
!> @param[in]  PhiR  values at the right of the boundary
!
!> @result           numerical flux at the boundary
!···············································································
function standardAvgFlux(PhiL, PhiR) result(F)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(:)
    real(wp), intent(in) :: PhiR(:)

    !* Return values *!
    real(wp) :: F(NEQS)

    F = 0.5_wp * ( EulerFlux(PhiL) + EulerFlux(PhiR) )

end function standardAvgFlux

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Volume flux of F. Ducros et al., High-order fluxes for conservative
!  skew-symmetric-like schemes in structured meshes: application to
!  compressible flows, J. Comput. Phys. 161 (2000) 114–139
!
!> @param[in]  PhiL  values at the left of the boundary
!> @param[in]  PhiR  values at the right of the boundary
!
!> @return           numerical flux at the boundary
!···············································································
function DucrosFlux(PhiL, PhiR) result(F)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(:)
    real(wp), intent(in) :: PhiR(:)

    !* Return values *!
    real(wp) :: F(NEQS)

    !* Local variables *!
    real(wp) :: pL
    real(wp) :: pR
    real(wp) :: rho
    real(wp) :: rhoU
    real(wp) :: u
    real(wp) :: p
    real(wp) :: rhoE

    ! Get pressure values
    pL = Phys%GammaMinusOne * ( PhiL(IRHOE) - 0.5_wp*PhiL(IRHOU)**2/PhiL(IRHO) )
    pR = Phys%GammaMinusOne * ( PhiR(IRHOE) - 0.5_wp*PhiR(IRHOU)**2/PhiR(IRHO) )

    ! Averaged values
    rho  = 0.5_wp * ( PhiL(IRHO)             + PhiR(IRHO)             )
    rhoU = 0.5_wp * ( PhiL(IRHOU)            + PhiR(IRHOU)            )
    u    = 0.5_wp * ( PhiL(IRHOU)/PhiL(IRHO) + PhiR(IRHOU)/PhiR(IRHO) )
    p    = 0.5_wp * ( pL                     + pR                     )
    rhoE = 0.5_wp * ( PhiL(IRHOE)            + PhiR(IRHOE)            )

    ! Ducros flux
    F(IRHO)  = rho * u
    F(IRHOU) = rhoU * u + p
    F(IRHOE) = (rhoE + p) * u

end function DucrosFlux

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Volume flux of C.A. Kennedy, A. Gruber, Reduced aliasing formulations of
!  the convective terms within the Navier–Stokes equations for a compressible
!  fluid, J. Comput. Phys. 227 (2008) 1676–1700
!
!> @param[in]  PhiL  values at the left of the boundary
!> @param[in]  PhiR  values at the right of the boundary
!
!> @return           numerical flux at the boundary
!···············································································
function KennedyGruberFlux(PhiL, PhiR) result(F)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(:)
    real(wp), intent(in) :: PhiR(:)

    !* Return values *!
    real(wp) :: F(NEQS)

    !* Local variables *!
    real(wp) :: pL
    real(wp) :: pR
    real(wp) :: rho
    real(wp) :: u
    real(wp) :: p
    real(wp) :: e
    real(wp) :: rhoU

    ! Get pressure values
    pL = Phys%GammaMinusOne * ( PhiL(IRHOE) - 0.5_wp*PhiL(IRHOU)**2/PhiL(IRHO) )
    pR = Phys%GammaMinusOne * ( PhiR(IRHOE) - 0.5_wp*PhiR(IRHOU)**2/PhiR(IRHO) )

    ! Averaged values
    rho  = 0.5_wp * ( PhiL(IRHO)             + PhiR(IRHO)             )
    u    = 0.5_wp * ( PhiL(IRHOU)/PhiL(IRHO) + PhiR(IRHOU)/PhiR(IRHO) )
    p    = 0.5_wp * ( pL                     + pR                     )
    e    = 0.5_wp * ( PhiL(IRHOE)/PhiL(IRHO) + PhiR(IRHOE)/PhiR(IRHO) )
    rhoU = rho * u

    ! Kennedy & Gruber flux
    F(IRHO)  = rhoU
    F(IRHOU) = rhoU * u + p
    F(IRHOE) = rhoU * e + p * u

end function KennedyGruberFlux

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Volume flux of S. Pirozzoli, Numerical methods for high-speed flows, Annu.
!  Rev. Fluid Mech. 43 (2011) 163–1
!
!> @param[in]  PhiL  values at the left of the boundary
!> @param[in]  PhiR  values at the right of the boundary
!
!> @return           numerical flux at the boundary
!···············································································
function PirozzoliFlux(PhiL, PhiR) result(F)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(:)
    real(wp), intent(in) :: PhiR(:)

    !* Return values *!
    real(wp) :: F(NEQS)

    !* Local variables *!
    real(wp) :: pL
    real(wp) :: pR
    real(wp) :: hL
    real(wp) :: hR
    real(wp) :: rho
    real(wp) :: u
    real(wp) :: p
    real(wp) :: h
    real(wp) :: rhoU

    ! Get pressure values
    pL = Phys%GammaMinusOne * ( PhiL(IRHOE) - 0.5_wp*PhiL(IRHOU)**2/PhiL(IRHO) )
    pR = Phys%GammaMinusOne * ( PhiR(IRHOE) - 0.5_wp*PhiR(IRHOU)**2/PhiR(IRHO) )

    ! Get enthalpy values
    hL = ( PhiL(IRHOE) + pL ) / PhiL(IRHO)
    hR = ( PhiR(IRHOE) + pR ) / PhiR(IRHO)

    ! Averaged values
    rho  = 0.5_wp * ( PhiL(IRHO)             + PhiR(IRHO)             )
    u    = 0.5_wp * ( PhiL(IRHOU)/PhiL(IRHO) + PhiR(IRHOU)/PhiR(IRHO) )
    p    = 0.5_wp * ( pL                     + pR                     )
    h    = 0.5_wp * ( hL                     + hR                     )
    rhoU = rho * u

    ! Pirozzoli flux
    F(IRHO)  = rhoU
    F(IRHOU) = rhoU * u + p
    F(IRHOE) = rhoU * h

end function PirozzoliFlux

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Volume flux of F. Ismail, P.L. Roe, Affordable, entropy-consistent Euler
!  flux functions II: entropy production at shocks, J. Comput. Phys. 228 (15)
!  (2009) 5410–5436. [20]
!
!> @param[in]  PhiL  values at the left of the boundary
!> @param[in]  PhiR  values at the right of the boundary
!
!> @return           numerical flux at the boundary
!···············································································
function IsmailRoeFlux(PhiL, PhiR) result(F)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(:)
    real(wp), intent(in) :: PhiR(:)

    !* Return values *!
    real(wp) :: F(NEQS)

    !* Local variables *!
    real(wp) :: pL
    real(wp) :: pR
    real(wp) :: z(3)
    real(wp) :: zL(3)
    real(wp) :: zR(3)
    real(wp) :: zLN(3)
    real(wp) :: rho
    real(wp) :: u
    real(wp) :: p1
    real(wp) :: p2
    real(wp) :: h

    ! Get pressure values
    pL = Phys%GammaMinusOne * ( PhiL(IRHOE) - 0.5_wp*PhiL(IRHOU)**2/PhiL(IRHO) )
    pR = Phys%GammaMinusOne * ( PhiR(IRHOE) - 0.5_wp*PhiR(IRHOU)**2/PhiR(IRHO) )

    ! Get Z parameter vector
    zL  = [ sqrt(PhiL(IRHO)/pL), PhiL(IRHOU)/sqrt(PhiL(IRHO)*pL), &
            sqrt(PhiL(IRHO)*pL) ];
    zR  = [ sqrt(PhiR(IRHO)/pR), PhiR(IRHOU)/sqrt(PhiR(IRHO)*pR), &
            sqrt(PhiR(IRHO)*pR) ];
    zLN = logarithmicMean(zL, zR)

    ! Averaged values
    z   = 0.5_wp * ( zL + zR )
    rho = z(1) * zLN(3)
    u   = z(2) / z(1)
    p1  = z(3) / z(1)
    p2  = Phys%GammaPlusOne / (2.0_wp*Phys%Gamma) * zLN(3)/zLN(1) + &
          Phys%GammaMinusOne / (2.0_wp*Phys%Gamma) * z(3)/z(1)
    h   = Phys%Gamma * p2 / (rho * Phys%GammaMinusOne) + 0.5_wp * u**2

    ! Ismail and Roe flux
    F(IRHO)  = rho * u
    F(IRHOU) = rho * u**2 + p1
    F(IRHOE) = rho * u * h

end function IsmailRoeFlux

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Volume flux of P. Chandrashekar, Kinetic energy preserving and entropy
!  stable finite volume schemes for compressible Euler and Navier–Stokes
!  equations, Commun. Comput. Phys. 14 (5) (2013) 1252–1286.
!
!> @param[in]  PhiL  values at the left of the boundary
!> @param[in]  PhiR  values at the right of the boundary
!
!> @return           numerical flux at the boundary
!···············································································
function ChandrashekarFlux(PhiL, PhiR) result(F)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(:)
    real(wp), intent(in) :: PhiR(:)

    !* Return values *!
    real(wp) :: F(NEQS)

    !* Local variables *!
    real(wp) :: uL
    real(wp) :: uR
    real(wp) :: RTL
    real(wp) :: RTR
    real(wp) :: betaL
    real(wp) :: betaR
    real(wp) :: rhoLN
    real(wp) :: betaLN
    real(wp) :: rho
    real(wp) :: u
    real(wp) :: u2
    real(wp) :: beta
    real(wp) :: p
    real(wp) :: h

    ! Some initial values
    uL    = PhiL(IRHOU) / PhiL(IRHO)
    uR    = PhiR(IRHOU) / PhiR(IRHO)
    RTL   = Phys%GammaMinusOne * ( PhiL(IRHOE)/PhiL(IRHO) - 0.5_wp * uL**2 )
    RTR   = Phys%GammaMinusOne * ( PhiR(IRHOE)/PhiR(IRHO) - 0.5_wp * uR**2 )
    betaL = 1.0_wp / (2.0_wp * RTL)
    betaR = 1.0_wp / (2.0_wp * RTR)

    ! Logarithmic values
    rhoLN  = logarithmicMean(PhiL(IRHO), PhiR(IRHO))
    betaLN = logarithmicMean(betaL,      betaR     )

    ! Standard average
    rho  = 0.5_wp * ( PhiL(IRHO) + PhiR(IRHO) )
    u    = 0.5_wp * ( uL         + uR         )
    u2   = 0.5_wp * ( uL * uL    + uR * uR    )
    beta = 0.5_wp * ( betaL      + betaR      )
    p    = rho / (2.0_wp * beta)
    h    = 0.5_wp / (betaLN*Phys%GammaMinusOne) - 0.5_wp * u2 + p/rhoLN + u**2

    ! KEPEC values
    F(IRHO)  = rhoLN * u
    F(IRHOU) = rhoLN * u**2 + p
    F(IRHOE) = rhoLN * u * h

end function ChandrashekarFlux

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Dissipation for the local Lax-Friedrichs flux
!
!> @param[in]  PhiL  values at the left of the boundary
!> @param[in]  PhiR  values at the right of the boundary
!
!> @return           numerical flux at the boundary
!···············································································
function localAvgDissipation(PhiL, PhiR) result(F)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(:)
    real(wp), intent(in) :: PhiR(:)

    !* Return values *!
    real(wp) :: F(NEQS)

    !* Local variables *!
    real(wp) :: lambdaL
    real(wp) :: lambdaR
    real(wp) :: lambda

    ! Get values of lambda on both sides
    lambdaL = computeMaxEigenvalue(PhiL)
    lambdaR = computeMaxEigenvalue(PhiR)

    ! Find the greater lambda number
    lambda = max(lambdaL, lambdaR)

    ! Stabilisation term
    F = 0.5_wp * lambda * ( PhiR - PhiL )

end function localAvgDissipation

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Matrix dissipation algorithm
!
!> @param[in]  PhiL  values at the left of the boundary
!> @param[in]  PhiR  values at the right of the boundary
!
!> @return           numerical flux at the boundary
!···············································································
function matrixDissipation(PhiL, PhiR) result(F)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(:)
    real(wp), intent(in) :: PhiR(:)

    !* Return values *!
    real(wp) :: F(NEQS)

    !* Local variables *!
    integer  :: i
    integer  :: j
    integer  :: k
    real(wp) :: wL(NEQS)
    real(wp) :: wR(NEQS)
    real(wp) :: uL
    real(wp) :: uR
    real(wp) :: pL
    real(wp) :: pR
    real(wp) :: betaL
    real(wp) :: betaR
    real(wp) :: rho
    real(wp) :: beta
    real(wp) :: p
    real(wp) :: u
    real(wp) :: u2
    real(wp) :: rhoLN
    real(wp) :: betaLN
    real(wp) :: a
    real(wp) :: h
    real(wp) :: R(NEQS, NEQS)
    real(wp) :: lambda(NEQS)
    real(wp) :: S(NEQS)

    call getEntropyVars(PhiL, wL)
    call getEntropyVars(PhiR, wR)

    ! Left / right values
    uL    = PhiL(IRHOU) / PhiL(IRHO)
    uR    = PhiR(IRHOU) / PhiR(IRHO)
    pL    = Phys%GammaMinusOne * ( PhiL(IRHOE) - 0.5_wp * uL * PhiL(IRHOU) )
    pR    = Phys%GammaMinusOne * ( PhiR(IRHOE) - 0.5_wp * uR * PhiR(IRHOU) )
    betaL = PhiL(IRHO) / ( 2.0_wp * pL )
    betaR = PhiR(IRHO) / ( 2.0_wp * pR )

    ! Mean values
    rho  = 0.5_wp * ( PhiL(IRHO) + PhiR(IRHO) )
    beta = 0.5_wp * ( betaL      + betaR      )
    u    = 0.5_wp * ( uL         + uR         )
    u2   = 0.5_wp * ( uL**2      + uR**2      )
    p    = rho / ( 2.0_wp * beta )

    ! Logarithmic values
    rhoLN  = logarithmicMean(PhiL(IRHO), PhiR(IRHO))
    betaLN = logarithmicMean(betaL,      betaR     )

    ! Some thermodynamic variables
    a = sqrt( Phys%Gamma * p / rhoLN )
    h = 0.5_wp / betaLN  * Phys%Gamma/Phys%GammaMinusOne + u**2 - 0.5_wp * u2

    ! More intermediate values
    R(:,1) = [ 1.0_wp, u - a, h - u * a        ]
    R(:,2) = [ 1.0_wp, u    , u**2 - 0.5_wp*u2 ]
    R(:,3) = [ 1.0_wp, u + a, h + u * a        ]

    lambda = abs([ u - a, u, u + a ])

    S(1)   = rhoLN / ( 2.0_wp * Phys%Gamma )
    S(2)   = 2.0_wp * S(1) * Phys%GammaMinusOne
    S(3)   = S(1)

    ! Iteratively construct the dissipative term
    F = 0.0_wp
    do i = 1, NEQS;     do j = 1, NEQS;     do k = 1, NEQS
        F(i) = F(i) + R(i,k) * lambda(k) * S(k) * R(j,k) * ( wR(j) - wL(j) )
    end do;             end do;             end do
    F = 0.5_wp * F

end function matrixDissipation

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Dissipation for the Imail & Roe flux
!  M. Carpenter, T. Fisher, E. Nielsen, S. Frankel, Entropy stable spectral
!  collocation schemes for the Navier–Stokes equations: discontinuous
!  interfaces, SIAM J. Sci. Comput. 36 (5) (2014) B835–B867
!
!> @param[in]  PhiL  values at the left of the boundary
!> @param[in]  PhiR  values at the right of the boundary
!
!> @return           numerical flux at the boundary
!···············································································
function IsmailRoeDissipation(PhiL, PhiR) result(F)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(:)
    real(wp), intent(in) :: PhiR(:)

    !* Return values *!
    real(wp) :: F(NEQS)

    !* Local variables *!
    real(wp) :: pL
    real(wp) :: pR
    real(wp) :: z(3)
    real(wp) :: zL(3)
    real(wp) :: zR(3)
    real(wp) :: zLN(3)
    real(wp) :: rho
    real(wp) :: u
    real(wp) :: p1
    real(wp) :: p2
    real(wp) :: h
    real(wp) :: e
    real(wp) :: wL(3)
    real(wp) :: wR(3)
    real(wp) :: aL
    real(wp) :: aR
    real(wp) :: a
    real(wp) :: lambda
    real(wp) :: R(3,3)

    ! Get pressure values
    pL = Phys%GammaMinusOne * ( PhiL(IRHOE) - 0.5_wp*PhiL(IRHOU)**2/PhiL(IRHO) )
    pR = Phys%GammaMinusOne * ( PhiR(IRHOE) - 0.5_wp*PhiR(IRHOU)**2/PhiR(IRHO) )

    ! Get Z parameter vector
    zL  = [ sqrt(PhiL(IRHO)/pL), PhiL(IRHOU)/sqrt(PhiL(IRHO)*pL), &
            sqrt(PhiL(IRHO)*pL) ];
    zR  = [ sqrt(PhiR(IRHO)/pR), PhiR(IRHOU)/sqrt(PhiR(IRHO)*pR), &
            sqrt(PhiR(IRHO)*pR) ];
    zLN = logarithmicMean(zL, zR)

    ! Averaged values
    z   = 0.5_wp * ( zL + zR )
    rho = z(1) * zLN(3)
    u   = z(2) / z(1)
    p1  = z(3) / z(1)
    p2  = Phys%GammaPlusOne / (2.0_wp*Phys%Gamma) * zLN(3)/zLN(1) + &
          Phys%GammaMinusOne / (2.0_wp*Phys%Gamma) * z(3)/z(1)
    h   = Phys%Gamma * p2 / (rho * Phys%GammaMinusOne) + 0.5_wp * u**2
    e   = h - p1 / rho
    a   = sqrt( Phys%Gamma * p1 / rho )

    ! Find the greater lambda number
    aL = sqrt(Phys%Gamma * pL / PhiL(IRHO))
    aR = sqrt(Phys%Gamma * pR / PhiR(IRHO))
    lambda = max((aL + abs(PhiL(IRHOU) / PhiL(IRHO))), &
                 (aR + abs(PhiR(IRHOU) / PhiR(IRHO))))

    ! Entropy variables
    call getEntropyVars(PhiL, wL)
    call getEntropyVars(PhiR, wR)

    ! Stabilisation term
    R(1,:) = [ rho,    rho * u,         rho * e                               ]
    R(2,:) = [ R(1,2), rho * u**2 + p1, rho * h * u                           ]
    R(3,:) = [ R(1,3), R(2,3),          rho*h**2 - a**2*p1/Phys%GammaMinusOne ]
    F = 0.5_wp * lambda * matmul(R, (wR - wL))

end function IsmailRoeDissipation

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Dissipation for the Chandrashekar flux
!  Gassner, G. J., Winters, A. R., & Kopriva, D. A. (2016). Split form nodal
!  discontinuous Galerkin schemes with summation-by-parts property for the
!  compressible Euler equations. Journal of Computational Physics, 327, 39–66.
!
!> @param[in]  PhiL  values at the left of the boundary
!> @param[in]  PhiR  values at the right of the boundary
!
!> @return           numerical flux at the boundary
!···············································································
function ChandrashekarDissipation(PhiL, PhiR) result(F)
    !* Arguments *!
    real(wp), intent(in)  :: PhiL(:)
    real(wp), intent(in)  :: PhiR(:)

    !* Return values *!
    real(wp) :: F(NEQS)

    !* Local variables *!
    real(wp) :: uL
    real(wp) :: uR
    real(wp) :: aL
    real(wp) :: aR
    real(wp) :: pL
    real(wp) :: pR
    real(wp) :: lambda
    real(wp) :: betaL
    real(wp) :: betaR
    real(wp) :: betaLN
    real(wp) :: rho
    real(wp) :: u

    ! Get pressures and calculate the two sound velocities
    uL = PhiL(IRHOU) / PhiL(IRHO)
    uR = PhiR(IRHOU) / PhiR(IRHO)
    pL = Phys%GammaMinusOne * (PhiL(IRHOE) - PhiL(IRHOU) * uL * 0.5_wp)
    pR = Phys%GammaMinusOne * (PhiR(IRHOE) - PhiR(IRHOU) * uR * 0.5_wp)
    aL = sqrt(Phys%Gamma * pL / PhiL(IRHO))
    aR = sqrt(Phys%Gamma * pR / PhiR(IRHO))

    ! Find the greater lambda number
    lambda = max((aL + abs(PhiL(IRHOU) / PhiL(IRHO))), &
                 (aR + abs(PhiR(IRHOU) / PhiR(IRHO))))

    ! Energy stabilisation term
    betaL = PhiL(IRHO) / ( 2.0_wp * pL)
    betaR = PhiR(IRHO) / ( 2.0_wp * pR)
    ! Averaged values
    betaLN = logarithmicMean(betaL, betaR)
    rho = 0.5_wp * ( PhiL(IRHO) + PhiR(IRHO) )
    u   = 0.5_wp * ( uL         + uR         )

    ! Stabilisation term
    F(IRHO)  = 0.5_wp * lambda * ( PhiR(IRHO)  - PhiL(IRHO)  )
    F(IRHOU) = 0.5_wp * lambda * ( PhiR(IRHOU) - PhiL(IRHOU) )
    F(IRHOE) = 0.5_wp * lambda * ( (0.5_wp/(Phys%GammaMinusOne*betaLN)  &
             + uR*uL) * (PhiR(IRHO) - PhiL(IRHO)) + rho*u * (uR - uL)   &
             + 0.5_wp*rho/Phys%GammaMinusOne * (1.0_wp/betaR - 1.0_wp/betaL) )

end function ChandrashekarDissipation

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Releases the used memory.
!···············································································
subroutine destructRiemannSolver()

    if (associated(TwoPointFlux))    nullify(TwoPointFlux)
    if (associated(DissipativeFlux)) nullify(DissipativeFlux)

end subroutine destructRiemannSolver

end module RiemannSolvers_NS
