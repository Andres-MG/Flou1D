!*******************************************************************************
!  MODULE: ManufacturedSolutions_CE
!
!> @author
!> Andres Mateo
!
!> @brief
!> Implementation of an oscillatory, 1D manufactured solution for the
!> Navier-Stokes equations with the artificial viscosity proposed by
!> Guermond & Popov.
!*******************************************************************************

module ManufacturedSolutions_NS

    use Constants

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: MMSsource
    public :: MMSexact

    ! Generic interfaces
    interface MMSsource
        module procedure MMSeulerSource
        module procedure MMSviscousSource
        module procedure MMSguermondSource
    end interface MMSsource

    interface MMSexact
        module procedure MMSexactScalar
        module procedure MMSexactVector
    end interface MMSexact

    ! Some local variables
    integer,  parameter :: IRHO  = 1
    integer,  parameter :: IRHOU = 2
    integer,  parameter :: IRHOE = 3
    real(wp), parameter :: Cr0   = 2.0_wp
    real(wp), parameter :: Cr1   = 0.1_wp
    real(wp), parameter :: kr    = PI
    real(wp), parameter :: wr    = 2.0_wp * PI
    real(wp), parameter :: Cu0   = -1.0_wp
    real(wp), parameter :: Cu1   = 1.0_wp
    real(wp), parameter :: Ce1   = 1.0_wp

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Source term of an inviscid flow with the following steady-state solution:
!
!>  \f[ \rho   = C_{r0} + C_{r1} \sin \left( k_r x - \omega_r t \right) \f]
!>  \f[ u      = C_{u0} + C_{u1} \rho \f]
!>  \f[ \rho e = C_{e1} \rho^2 \f]
!
!> @param[in]   gamma  specific heat ratio
!> @param[in]   x      coordinates of the nodes
!> @param[in]   t      time instant
!> @param[out]  S      values at the nodes
!···············································································
subroutine MMSeulerSource(gamma, x, t, S)
    !* Arguments *!
    real(wp), intent(in)  :: gamma
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: t
    real(wp), intent(out) :: S(:,:)

    !* Local variables *!
    real(wp), allocatable :: rho(:)
    real(wp), allocatable :: u(:)
    real(wp), allocatable :: e(:)
    real(wp), allocatable :: ei(:)
    real(wp), allocatable :: p(:)

    real(wp), allocatable :: rho_t(:)
    real(wp), allocatable :: u_t(:)
    real(wp), allocatable :: e_t(:)

    real(wp), allocatable :: rho_x(:)
    real(wp), allocatable :: u_x(:)
    real(wp), allocatable :: e_x(:)
    real(wp), allocatable :: ei_x(:)
    real(wp), allocatable :: p_x(:)

    real(wp) :: Phi_dot(size(x), 3)
    real(wp) :: Fe_x(size(x), 3)

    ! Variables and their derivatives
    rho = Cr0 + Cr1 * sin(kr*x - wr*t)
    u   = Cu0 + Cu1 * rho
    e   = Ce1 * rho
    ei  = e - 0.5_wp * u**2
    p   = (gamma - 1.0_wp) * rho * ei

    rho_t = -Cr1 * wr * cos(kr*x - wr*t)
    u_t   = Cu1 * rho_t
    e_t   = Ce1 * rho_t

    rho_x = Cr1 * kr * cos(kr*x - wr*t)
    u_x   = Cu1 * rho_x
    e_x   = Ce1 * rho_x
    ei_x  = e_x - u * u_x
    p_x   = (gamma - 1.0_wp) * (rho_x*ei + rho*ei_x)

    ! Time derivatives of the conservative variables
    Phi_dot(:, IRHO)  = rho_t
    Phi_dot(:, IRHOU) = rho_t*u + rho*u_t
    Phi_dot(:, IRHOE) = rho_t*e + rho*e_t

    ! Divergence of the inviscid flux
    Fe_x(:, IRHO)  = rho_x*u + rho*u_x
    Fe_x(:, IRHOU) = rho_x*u**2 + 2.0_wp*rho*u*u_x + p_x
    Fe_x(:, IRHOE) = (rho_x*e + rho*e_x + p_x) * u + (rho*e + p) * u_x

    ! Source term
    S = Phi_dot + Fe_x

end subroutine MMSeulerSource

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Source term of a viscous flow with the following steady-state solution:
!
!>  \f[ \rho   = C_{r0} + C_{r1} \sin \left( k_r x - \omega_r t \right) \f]
!>  \f[ u      = C_{u0} + C_{u1} \rho \f]
!>  \f[ \rho e = C_{e1} \rho^2 \f]
!
!> @param[in]   Pr     Prandtl number
!> @param[in]   mu     NS viscosity
!> @param[in]   gamma  specific heat ratio
!> @param[in]   x      coordinates of the nodes
!> @param[in]   t      time instant
!> @param[out]  S      values at the nodes
!···············································································
subroutine MMSviscousSource(Pr, mu, gamma, x, t, S)
    !* Arguments *!
    real(wp), intent(in)  :: Pr
    real(wp), intent(in)  :: mu
    real(wp), intent(in)  :: gamma
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: t
    real(wp), intent(out) :: S(:,:)

    !* Local variables *!
    real(wp), allocatable :: rho(:)
    real(wp), allocatable :: u(:)
    real(wp), allocatable :: e(:)
    real(wp), allocatable :: ei(:)
    real(wp), allocatable :: p(:)

    real(wp), allocatable :: rho_t(:)
    real(wp), allocatable :: u_t(:)
    real(wp), allocatable :: e_t(:)

    real(wp), allocatable :: rho_x(:)
    real(wp), allocatable :: u_x(:)
    real(wp), allocatable :: e_x(:)
    real(wp), allocatable :: ei_x(:)
    real(wp), allocatable :: p_x(:)

    real(wp), allocatable :: rho_xx(:)
    real(wp), allocatable :: u_xx(:)
    real(wp), allocatable :: e_xx(:)
    real(wp), allocatable :: ei_xx(:)

    real(wp) :: Phi_dot(size(x), 3)
    real(wp) :: Fe_x(size(x), 3)
    real(wp) :: Fv_x(size(x), 3)

    ! Variables and their derivatives
    rho = Cr0 + Cr1 * sin(kr*x - wr*t)
    u   = Cu0 + Cu1 * rho
    e   = Ce1 * rho
    ei  = e - 0.5_wp * u**2
    p   = (gamma - 1.0_wp) * rho * ei

    rho_t = -Cr1 * wr * cos(kr*x - wr*t)
    u_t   = Cu1 * rho_t
    e_t   = Ce1 * rho_t

    rho_x = Cr1 * kr * cos(kr*x - wr*t)
    u_x   = Cu1 * rho_x
    e_x   = Ce1 * rho_x
    ei_x  = e_x - u * u_x
    p_x   = (gamma - 1.0_wp) * (rho_x*ei + rho*ei_x)

    rho_xx = -Cr1 * kr**2 * sin(kr*x - wr*t)
    u_xx   = Cu1 * rho_xx
    e_xx   = Ce1 * rho_xx
    ei_xx  = e_xx - u_x**2 - u*u_xx

    ! Time derivatives of the conservative variables
    Phi_dot(:, IRHO)  = rho_t
    Phi_dot(:, IRHOU) = rho_t*u + rho*u_t
    Phi_dot(:, IRHOE) = rho_t*e + rho*e_t

    ! Divergence of the inviscid flux
    Fe_x(:, IRHO)  = rho_x*u + rho*u_x
    Fe_x(:, IRHOU) = rho_x*u**2 + 2.0_wp*rho*u*u_x + p_x
    Fe_x(:, IRHOE) = (rho_x*e + rho*e_x + p_x) * u + (rho*e + p) * u_x

    ! Divergence of the NS viscous flux
    Fv_x(:, IRHO)  = 0.0_wp
    Fv_x(:, IRHOU) = 4.0_wp/3.0_wp * mu*u_xx
    Fv_x(:, IRHOE) = 4.0_wp/3.0_wp * mu*(u_x**2 + u*u_xx) + mu*gamma/Pr * ei_xx

    ! Source term
    S = Phi_dot + Fe_x - Fv_x

end subroutine MMSviscousSource

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Source term of an NS flow with artificial viscosity according to
!> Guermond, J. L., & Popov, B. (2014). Viscous regularization of the Euler
!> equations and entropy principles. SIAM Journal on Applied Mathematics,
!> 74(2), 284–305. https://doi.org/10.1137/120903312
!
!>  \f[ f_{GP,1}^a = \alpha \rho_x \f],
!>  \f[ f_{GP,2}^a = \alpha u \rho_x + \beta u_x \f],
!>  \f[ f_{GP,3}^a = \alpha (\rho e_i)_x + \frac{u^2}{2}\rho_x + \beta u u_x \\
!>      + \lambda (ei)_x \f]
!
!> and the following steady-state solution:
!
!>  \f[ \rho   = C_{r0} + C_{r1} \sin \left( k_r x - \omega_r t \right) \f]
!>  \f[ u      = C_{u0} + C_{u1} \rho \f]
!>  \f[ \rho e = C_{e1} \rho^2 \f]
!
!> @param[in]   Pr      Prandtl number
!> @param[in]   mu      NS viscosity
!> @param[in]   alpha   artificial viscosity coeff.
!> @param[in]   beta    artificial viscosity coeff.
!> @param[in]   lambda  artificial viscosity coeff.
!> @param[in]   gamma   specific heat ratio
!> @param[in]   x       coordinates of the nodes
!> @param[in]   t       time instant
!> @param[out]  S       values at the nodes
!···············································································
subroutine MMSguermondSource(Pr, mu, alpha, beta, lambda, gamma, x, t, S)
    !* Arguments *!
    real(wp), intent(in)  :: Pr
    real(wp), intent(in)  :: mu
    real(wp), intent(in)  :: alpha
    real(wp), intent(in)  :: beta
    real(wp), intent(in)  :: lambda
    real(wp), intent(in)  :: gamma
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: t
    real(wp), intent(out) :: S(:,:)

    !* Local variables *!
    real(wp), allocatable :: rho(:)
    real(wp), allocatable :: u(:)
    real(wp), allocatable :: e(:)
    real(wp), allocatable :: ei(:)
    real(wp), allocatable :: p(:)

    real(wp), allocatable :: rho_t(:)
    real(wp), allocatable :: u_t(:)
    real(wp), allocatable :: e_t(:)

    real(wp), allocatable :: rho_x(:)
    real(wp), allocatable :: u_x(:)
    real(wp), allocatable :: e_x(:)
    real(wp), allocatable :: ei_x(:)
    real(wp), allocatable :: p_x(:)

    real(wp), allocatable :: rho_xx(:)
    real(wp), allocatable :: u_xx(:)
    real(wp), allocatable :: e_xx(:)
    real(wp), allocatable :: ei_xx(:)
    real(wp), allocatable :: rho_ei_xx(:)

    real(wp) :: Phi_dot(size(x), 3)
    real(wp) :: Fe_x(size(x), 3)
    real(wp) :: Fv_x(size(x), 3)
    real(wp) :: Fgp_x(size(x), 3)

    ! Variables and their derivatives
    rho = Cr0 + Cr1 * sin(kr*x - wr*t)
    u   = Cu0 + Cu1 * rho
    e   = Ce1 * rho
    ei  = e - 0.5_wp * u**2
    p   = (gamma - 1.0_wp) * rho * ei

    rho_t = -Cr1 * wr * cos(kr*x - wr*t)
    u_t   = Cu1 * rho_t
    e_t   = Ce1 * rho_t

    rho_x = Cr1 * kr * cos(kr*x - wr*t)
    u_x   = Cu1 * rho_x
    e_x   = Ce1 * rho_x
    ei_x  = e_x - u * u_x
    p_x   = (gamma - 1.0_wp) * (rho_x*ei + rho*ei_x)

    rho_xx = -Cr1 * kr**2 * sin(kr*x - wr*t)
    u_xx   = Cu1 * rho_xx
    e_xx   = Ce1 * rho_xx
    ei_xx  = e_xx - u_x**2 - u*u_xx

    rho_ei_xx = rho_xx*ei + 2.0_wp*rho_x*ei_x + rho*ei_xx

    ! Time derivatives of the conservative variables
    Phi_dot(:, IRHO)  = rho_t
    Phi_dot(:, IRHOU) = rho_t*u + rho*u_t
    Phi_dot(:, IRHOE) = rho_t*e + rho*e_t

    ! Divergence of the inviscid flux
    Fe_x(:, IRHO)  = rho_x*u + rho*u_x
    Fe_x(:, IRHOU) = rho_x*u**2 + 2.0_wp*rho*u*u_x + p_x
    Fe_x(:, IRHOE) = (rho_x*e + rho*e_x + p_x) * u + (rho*e + p) * u_x

    ! Divergence of the NS viscous flux
    Fv_x(:, IRHO)  = 0.0_wp
    Fv_x(:, IRHOU) = 4.0_wp/3.0_wp * mu*u_xx
    Fv_x(:, IRHOE) = 4.0_wp/3.0_wp * mu*(u_x**2 + u*u_xx) + mu*gamma/Pr * ei_xx

    ! Divergence of the GP artificial flux
    Fgp_x(:, IRHO)  = alpha*rho_xx
    Fgp_x(:, IRHOU) = alpha*(u_x*rho_x + u*rho_xx) + beta*u_xx
    Fgp_x(:, IRHOE) = alpha*(rho_ei_xx + 0.5_wp*rho_xx*u**2 + u*u_x*rho_x) &
                    + beta*(u_x**2 + u*u_xx) + lambda * ei_xx

    ! Source term
    S = Phi_dot + Fe_x - Fv_x - Fgp_x

end subroutine MMSguermondSource

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Exact solution of the form:
!
!>  \f[ \rho   = C_{r0} + C_{r1} \sin \left( k_r x - \omega_r t \right) \f]
!>  \f[ u      = C_{u0} + C_{u1} \rho \f]
!>  \f[ \rho e = C_{e1} \rho^2 \f]
!
!> @param[in]   x      coordinates of the nodes
!> @param[in]   t      time instant
!> @param[out]  Phi    values at the nodes
!> @param[out]  Grad   gradients at the nodes (optional)
!···············································································
subroutine MMSexactScalar(x, t, Phi, Grad)
    !* Arguments *!
    real(wp),           intent(in)  :: x
    real(wp),           intent(in)  :: t
    real(wp),           intent(out) :: Phi(:)
    real(wp), optional, intent(out) :: Grad(:)

    Phi(IRHO)  = Cr0 + Cr1 * sin(kr*x - wr*t)
    Phi(IRHOU) = Phi(IRHO) * ( Cu0 + Cu1 * Phi(IRHO) )
    Phi(IRHOE) = Ce1 * Phi(IRHO) ** 2

    if (present(Grad)) then
        Grad(IRHO)  = Cr1 * kr * cos(kr*x - wr*t)
        Grad(IRHOU) = Cu0 * Grad(IRHO) + 2.0_wp * Cu1 * Phi(IRHO) * Grad(IRHO)
        Grad(IRHOE) = 2.0_wp * Ce1 * Phi(IRHO) * Grad(IRHO)
    end if

end subroutine MMSexactScalar

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Exact solution of the form:
!
!>  \f[ \rho   = C_{r0} + C_{r1} \sin \left( k_r x - \omega_r t \right) \f]
!>  \f[ u      = C_{u0} + C_{u1} \rho \f]
!>  \f[ \rho e = C_{e1} \rho^2 \f]
!
!> @param[in]   x      coordinates of the nodes
!> @param[in]   t      time instant
!> @param[out]  Phi    values at the nodes
!> @param[out]  Grad   gradients at the nodes (optional)
!···············································································
subroutine MMSexactVector(x, t, Phi, Grad)
    !* Arguments *!
    real(wp),           intent(in)  :: x(:)
    real(wp),           intent(in)  :: t
    real(wp),           intent(out) :: Phi(:,:)
    real(wp), optional, intent(out) :: Grad(:,:)

    Phi(:,IRHO)  = Cr0 + Cr1 * sin(kr*x - wr*t)
    Phi(:,IRHOU) = Phi(:,IRHO) * ( Cu0 + Cu1 * Phi(:,IRHO) )
    Phi(:,IRHOE) = Ce1 * Phi(:,IRHO) ** 2

    if (present(Grad)) then
        Grad(:,IRHO)  = Cr1 * kr * cos(kr*x - wr*t)
        Grad(:,IRHOU) = Cu0*Grad(:,IRHO) + 2.0_wp*Cu1*Phi(:,IRHO)*Grad(:,IRHO)
        Grad(:,IRHOE) = 2.0_wp * Ce1 * Phi(:,IRHO) * Grad(:,IRHO)
    end if

end subroutine MMSexactVector

end module ManufacturedSolutions_NS
