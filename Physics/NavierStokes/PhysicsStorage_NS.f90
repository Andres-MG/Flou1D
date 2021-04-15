!*******************************************************************************
!  MODULE: PhysicsStorage_NS
!
!> @author
!> Andres Mateo
!
!> @brief
!> Module defining a class to store all the data related to the physics of the
!> compressible Navier-Stokes equations.
!*******************************************************************************

module PhysicsStorage_NS

    use Constants

    implicit none

    ! Explicitly define public types
    public :: PhysicsNS_t

    ! Explicitly define public variables
    public :: Phys

    ! Explicitly define public parameters
    integer, parameter, public :: NEQS  = 3
    integer, parameter, public :: IRHO  = 1
    integer, parameter, public :: IRHOU = 2
    integer, parameter, public :: IRHOE = 3
    character(len=5), parameter, public :: VarNames(3) = &
                                           ["rho  ", "rho*u", "rho*e"]

!···············································································
!> @class PhysicsNS_t
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Class representing the physical properties of the problem.
!···············································································
    type :: PhysicsNS_t
        real(wp) :: Gamma
        real(wp) :: GammaPlusOne
        real(wp) :: GammaMinusOne
        real(wp) :: Mu
        real(wp) :: Pr
        real(wp) :: Alpha(3)
        real(wp) :: SSFVc
        logical  :: WithEntropyVars
        logical  :: IsViscous
    contains
        procedure :: construct => NS_constructor
    end type PhysicsNS_t

    type(PhysicsNS_t) :: Phys

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initialises some characteristic constants of the compressible NS equations.
!
!> @param[in]  gam           specific heat ratio
!> @param[in]  mu            dynamic viscosity [Pa·s]
!> @param[in]  Pr            Prandtl number
!> @param[in]  alpha         max. artificial viscosity
!> @param[in]  alpha2beta    ratio between the 2-1 art. visc. coefficients
!> @param[in]  alpha2lambda  ratio between the 3-1 art. visc. coefficients
!> @param[in]  SSFVblending  SSFV blending constant
!···············································································
subroutine NS_constructor(this, gam, mu, Pr, alpha, alpha2beta, alpha2lambda, &
                          SSFVblending)
    !* Arguments *!
    real(wp), intent(in) :: gam
    real(wp), intent(in) :: mu
    real(wp), intent(in) :: Pr
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: alpha2beta
    real(wp), intent(in) :: alpha2lambda
    real(wp), intent(in) :: SSFVblending
    ! Derived types
    class(PhysicsNS_t), intent(inout) :: this

    ! Euler equation constants
    this%Gamma         = gam
    this%GammaPlusOne  = gam + 1.0_wp
    this%GammaMinusOne = gam - 1.0_wp

    ! Navier-Stokes equations constants
    if (mu == 0) then
        this%IsViscous = .false.
    else
        this%IsViscous = .true.
    end if
    this%Mu = mu
    this%Pr = Pr

    ! Default to conservative variables
    this%WithEntropyVars = .false.

    ! Artificial viscosity constants
    this%Alpha(1) = alpha
    this%Alpha(2) = alpha2beta
    this%Alpha(3) = alpha2lambda

    this%SSFVc = SSFVblending

end subroutine NS_constructor

end module PhysicsStorage_NS
