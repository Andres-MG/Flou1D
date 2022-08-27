!*******************************************************************************
!  MODULE: Setup_params
!
!> @author
!> Andres Mateo
!
!> @brief
!> Definition of the input parameters for the simulation .
!*******************************************************************************

module Setup_params

    use Constants

    implicit none

    ! Parameters of the simulation
    integer,  parameter :: P        = 9                 ! Initial expansion order (P)
    integer,  parameter :: K        = 20                ! Number of elements
    real(wp), parameter :: MAXRES   = 1e-8_wp           ! Max residual allowed
    real(wp), parameter :: TSPAN(2) = [0.0_wp, 0.18_wp] ! Time span
    real(wp), parameter :: TSTEP    = 1e-4_wp           ! Time step
    logical,  parameter :: PERIODIC = .false.           ! Periodic BCs
    integer,  parameter :: SAVEINT  = 100               ! Steps between prints
    integer,  parameter :: SAVERES  = 0                 ! Points per element
    integer,  parameter :: CFLSTEP  = 400               ! Steps between CFL

    !-> Solution file name
    !       - Tecplot format: .dat
    !       - Gnuplot format: .gp
    character(len=CHAR_LEN) :: FILENAME = "./results/test.gp"

    !-> Mathematical discretization
    !       Advection term
    !           - Weak:          eWeakAdvection
    !           - Strong:        eStrongAdvection
    !           - Split form:    eSplitAdvection
    !           - WENO:          eWENOadvection
    !           - SSWENO:        eSSWENOadvection
    !           - FV:            eFVadvection
    !           - SSFV:          eSSFVadvection
    integer, parameter :: ADVECTION = eSplitAdvection

    !-> Advection of elements marked by the sensor
    !           - Same as above: fNone
    !           - Weak:          eWeakAdvection
    !           - ...
    integer, parameter :: SENSEDADVECTION = fNone

    !       Gradients
    !           - Weak:   eWeakGradient
    !           - Strong: eStrongGradient
    integer, parameter :: GRADIENTS = eWeakGradient

    !       Viscous terms
    !           - Bassi-Rebay 1 (BR1): eBR1
    integer, parameter :: VISCOSITY = eBR1

    !       Quadrature nodes type
    !           - Gauss:         eGauss
    !           - Gauss-Lobatto: eGaussLobatto
    integer, parameter :: NODETYPE = eGaussLobatto

    !       Split form & Riemann solver
    !           Two-point flux type
    !               - Standard average: eStdAvg
    !               - Ducros:           eDucros
    !               - Kenedy & Gruber:  eKennedyGruber
    !               - Pirozzoli:        ePirozzoli
    !               - Ismail & Roe:     eIsmailRoe
    !               - Chandrashekar:    eChandrashekar
    integer, parameter :: TWOPOINTTYPE = eChandrashekar

    !           Stabilisation term for the Riemann Solver
    !               - Local LxF:          eLocalLxF
    !               - Matrix dissipation: eMatrixDiss
    !               - Ismail & Roe:       eIRdiss
    !               _ Chandrashekar:      eCHdiss
    integer, parameter :: DISSTYPE = eMatrixDiss

    !-> Time integration method
    !       - Runge-Kutta 3: eRK3
    !       - Runge-Kutta 5: eRK5
    integer,  parameter :: INT_METHOD = eRK3

    !-> Physical parameters
    !       Euler equations
    real(wp), parameter :: GAMMA = 1.4_wp

    !       Viscous terms
    !           - Navier-Stokes viscosity
    real(wp), parameter :: MU = 0.0_wp
    real(wp), parameter :: PR = 0.71_wp

    !-> Artificial viscosity formulation
    !       - No added viscosity:            fNone
    !       - Laplacian flux:                eLaplacianVisc
    !       - GP flux w/ physical variables: eGuermondPhysical
    !       - GP flux w/ entropy variables:  eGuermondEntropy
    integer,  parameter :: ARTVISC      = eGuermondEntropy
    real(wp), parameter :: ALPHAMAX     = 0.0003_wp * (P+1) * K
    real(wp), parameter :: ALPHA2BETA   = 1.0_wp
    real(wp), parameter :: ALPHA2LAMBDA = 0.0_wp

    !-> Entropy-stable FV-DGSEM blending constant
    !       - Higher is "more" DGSEM
    !       - Lower is "more" FV
    real(wp), parameter :: SSFVBLEND = 1e+2_wp

    !-> Spectral Vanishing Viscosity (SVV)
    !       Type of SVV
    !           - No SVV:                       fNone
    !           - GP flux w/ entropy variables: eGuermondEntropySVV
    !
    !       - SVVPOW: Exponent of the filtering law (i/P)**SVVPOW
    integer,  parameter :: SVVTYPE         = fNone
    integer,  parameter :: SVVPOW          = 1
    real(wp), parameter :: ALPHASVV        = 0.0001_wp
    real(wp), parameter :: ALPHA2BETASVV   = 1.0_wp
    real(wp), parameter :: ALPHA2LAMBDASVV = 0.0_wp

    !-> Sensor definition
    !       Type of sensor
    !           - Aliasing error:   eAliasingSensor
    !           - Truncation error: eTruncationError
    !           - Modal sensor:     eModalSensor
    !           - Density sensor:   eDensitySensor
    !           - Jump sensor:      eJumpSensor
    !           - kMeans sensor:    eKmeansSensor
    !           - GMM sensor:       eGMMSensor
    integer, parameter :: SENSORTYPE    = eGMMSensor
    integer, parameter :: SECSENSORTYPE = eJumpSensor

    !       Shape of the activation function
    !           - RAMPTOP/BOTTOM for elements of orden P>0
    !           - SECRAMPTOP/BOTTOM for elements of order P=0
    !           - Bottom values serve as a threshold to trigger SVV and SSFV
    !
    !  Scaled value |    .  .  ._____
    !               |    .  . /.
    !               |    .  ./ .
    !               |    .  .  .
    !               |    . /.  .
    !               |____./ .  .
    !               |____.__.__._____ Sensor
    !
    real(wp), parameter :: RAMPTOP         = -8.0_wp !-10.0_wp
    real(wp), parameter :: RAMPBOTTOM      = -15.0_wp !-20.0_wp
    real(wp), parameter :: SECRAMPTOP      = -1.0_wp
    real(wp), parameter :: SECRAMPBOTTOM   = -2.0_wp
    integer,  parameter :: SENSORWINDOW(2) = [0, huge(1)]
    integer,  parameter :: SENSORSTEP      = 1

    !       Number of clusters for clustering sensors
    integer, parameter :: NUMCLUSTERS = 4

    !       Sensed variable
    !           - rho*p: 0
    !           - rho:   1
    !           - rho*u: 2
    !           - rho*e: 3
    integer, parameter :: SENSEDVAR = 0

    !       Unsteady correction of the truncation error
    logical, parameter :: TECORRECTION = .true.

    !       Truncation error estimator
    !           - No TE estimation: fNone
    !           - Non-isolated TE:  eLocalTrunc
    !           - Isolated TE:      eIsolatedLocalTrunc
    integer, parameter :: TRUNCERRORTYPE = eIsolatedLocalTrunc

    !       Compute the truncation error map
    logical, parameter :: TEMAP      = .false.
    integer, parameter :: LOWERLIMIT = 1

    ! -> WENO stencil size
    integer,  parameter :: KWENO = 3

    !-> h-adaptation
    !       Type of adaptation
    !           - No adaptation:      fNone
    !           - Dynamic adaptation: eDynamicAdaptation
    integer, parameter :: ADAPTTYPE = fNone

    !       Threshold value of the sensor
    real(wp), parameter :: ADAPTHIGHTHRES   = RAMPTOP
    real(wp), parameter :: ADAPTLOWTHRES    = -1.5_wp
    real(wp), parameter :: SECADAPTLOWTHRES = -2.0_wp

    !       Limit values of the approximation order (P)
    integer, parameter :: ADAPTMINORDER = 2
    integer, parameter :: ADAPTMAXORDER = 6

    !       Limit sizes of the elements
    real(wp), parameter :: MINSIZE = 0.005_wp
    real(wp), parameter :: MAXSIZE = 0.05_wp

    !       First and last dynamic adaptation timesteps
    integer, parameter :: FIRSTDYNSTEP = 0
    integer, parameter :: LASTDYNSTEP  = huge(1)

    !       Timesteps between mesh adaptations
    integer, parameter :: DYNAMICSTEP = 30

    !-> Stop at NaN (for debugging purposes):
    logical, parameter :: STOPNAN = .true.

end module Setup_params

!*******************************************************************************
!  MODULE: Setup_routines
!
!> @author
!> Andres Mateo
!
!> @brief
!> Definition of some user defined functions, such as boundary conditions, the
!> initial state and the exact solution.
!*******************************************************************************

module Setup_routines

    use Constants
    use Physics
    use Tests
    use Utilities

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: meshDefinition
    public :: initialCondition
    public :: externalState
    public :: sourceTerm
    public :: exactSolution

    ! Local variables
    real(wp), parameter :: xL    = 0.0_wp
    real(wp), parameter :: xR    = 1.0_wp

    real(wp), parameter :: Pr    = 0.71_wp
    real(wp), parameter :: gamma = 1.4_wp
    real(wp), parameter :: mu    = 0.0_wp

    real(wp), parameter :: alpha  = 0.0_wp
    real(wp), parameter :: beta   = 0.0_wp
    real(wp), parameter :: lambda = 0.0_wp

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initialisation of the mesh nodes
!
!> @param[in]   k   number of elements in the mesh
!> @param[out]  xk  nodes of the mesh
!···············································································
subroutine meshDefinition(k, xk)
    !* Arguments *!
    integer,  intent(in)  :: k
    real(wp), intent(out) :: xk(:)

    !* Local variables *!
    integer  :: i
    real(wp) :: x

    ! Construct mesh
    !call random_number(xk)
    !xk = 2.0_wp * (xk+0.2_wp) / sum(xk+0.2_wp)
    !x = -1.0_wp
    do i = 1, k+1
        x = 2.0_wp*(i-1)/k-1.0_wp
        !x = x + xk(i)
        call linear_mapping(xL, xR, x, xk(i))
    end do

end subroutine meshDefinition

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Definition of the initial condition
!
!> @param[in]   x    coordinates to evaluate the IC
!> @param[out]  Phi  prescribed value at the beginning
!···············································································
subroutine initialCondition(x, Phi)
    !* Arguments *!
    real(wp), intent(in)  :: x(:)
    real(wp), intent(out) :: Phi(:,:)

    !call MMSexact(x, 0.0_wp, Phi)
    call SodsTubeIC(x, Phi)
    ! call ShuOsherIC(x, gamma, Phi)

end subroutine initialCondition

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Definition of boundary conditions
!
!> @param[in]   x          coordinates to evaluate the BC
!> @param[in]   t          time instant
!> @param[in]   PhiInt     interpolated values at the boundary
!> @param[in]   gradInt    interpolated values of the grad. at the boundary
!> @param[out]  PhiExt     prescribed value at the boundary
!> @param[out]  gradExt    prescribed value of the gradient at the boundary
!···············································································
subroutine externalState(x, t, PhiInt, gradInt, PhiExt, gradExt)
    !* Arguments *!
    real(wp), intent(in)    :: x
    real(wp), intent(in)    :: t
    real(wp), intent(inout) :: PhiInt(:)
    real(wp), intent(inout) :: gradInt(:)
    real(wp), intent(out)   :: PhiExt(:)
    real(wp), intent(out)   :: gradExt(:)

    !call MMSexact(x, t, PhiExt)
    call SodsTubeIC(x, PhiExt)
    ! call ShuOsherIC(x, gamma, PhiExt)

    gradExt = gradInt

end subroutine externalState

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Definition of the source term.
!
!> @param[in]   x  coordinates to evaluate the source term
!> @param[in]   t  time instant
!> @param[out]  S  value of the source term at x
!···············································································
subroutine sourceTerm(x, t, S)
    !* Arguments *!
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: t
    real(wp), intent(out) :: S(:,:)

    !call MMSsource(Pr, mu, alpha, beta, lambda, gamma, x, t, S)
    S = 0.0_wp

end subroutine sourceTerm

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Exact solution of the problem. Returns .false. if it is unknown
!
!> @param[in]   x    coordinates to evaluate the IC
!> @param[in]   t    time instant
!> @param[out]  Phi  prescribed value at the beginning
!
!> @return .false. if the results are not valid, i.e. the solution is unknown
!···············································································
function exactSolution(x, t, Phi)
    !* Arguments *!
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: t
    real(wp), intent(out) :: Phi(:,:)

    !* Return value *!
    logical :: exactSolution

    !call MMSexact(x, t, Phi)
    exactSolution = .false.
    Phi = 0.0_wp

end function exactSolution

end module Setup_routines
