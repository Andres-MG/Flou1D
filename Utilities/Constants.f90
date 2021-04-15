module Constants

    use iso_fortran_env, only: stdin  => input_unit,  &
                               stdout => output_unit, &
                               stderr => error_unit,  &
                               r64    => real64,      &
                               r32    => real32

    implicit none

    !! Global precision
#ifdef SINGLE_PRECISION
    integer,          parameter :: wp = r32
    character(len=*), parameter :: REAL_FMT = "e16.6"
#else
    integer,          parameter :: wp = r64
    character(len=*), parameter :: REAL_FMT = "e26.16"
#endif

    !! Formats for outputting time
    character(len=*), parameter :: TIME0_FMT = "f7.5"
    character(len=*), parameter :: TIME1_FMT = "f0.5"

    !! Numerical constants
    real(wp), parameter :: PI       = acos(-1.0_wp)
    real(wp), parameter :: FLT_EPS  = epsilon(1.0_wp)
    integer,  parameter :: CHAR_LEN = 100

    !! Flags and ennumerations

    ! Enum: types of nodes
    integer, parameter :: eGauss        = 10
    integer, parameter :: eGaussLobatto = 11

    ! Enum: Face extrapolation type
    integer, parameter :: eWENO5 = 20

    ! Enum: symmetric numerical flux
    integer, parameter :: eStdAvg        = 40
    integer, parameter :: eDucros        = 41
    integer, parameter :: eKennedyGruber = 42
    integer, parameter :: ePirozzoli     = 43
    integer, parameter :: eIsmailRoe     = 44
    integer, parameter :: eChandrashekar = 45

    ! Enum: skew numerical flux (numerical stabilisation)
    integer, parameter :: eLocalLxF   = 50
    integer, parameter :: eMatrixDiss = 51
    integer, parameter :: eIRdiss     = 52
    integer, parameter :: eCHdiss     = 53

    ! Enum: advection formulation
    integer, parameter :: eWeakAdvection   = 60
    integer, parameter :: eStrongAdvection = 61
    integer, parameter :: eSplitAdvection  = 62
    integer, parameter :: eWENOadvection   = 63
    integer, parameter :: eSSWENOadvection = 64
    integer, parameter :: eFVadvection     = 65
    integer, parameter :: eSSFVadvection   = 66

    ! Enum: gradient formulation
    integer, parameter :: eWeakGradient   = 70
    integer, parameter :: eStrongGradient = 71

    ! Enum: viscous term formulation
    integer, parameter :: eBR1 = 80

    ! Enum: artificial viscosity flux types
    integer, parameter :: eLaplacianVisc    = 90
    integer, parameter :: eGuermondPhysical = 91
    integer, parameter :: eGuermondEntropy  = 92

    ! Enum: time integration schemes
    integer, parameter :: eInvalid = 100
    integer, parameter :: eRK3     = 101
    integer, parameter :: eRK5     = 102

    ! Enum: sensor types & h/p adaption
    integer, parameter :: eDynamicAdaptation  = 110
    integer, parameter :: eLocalTrunc         = 111
    integer, parameter :: eIsolatedLocalTrunc = 112
    integer, parameter :: eTruncationError    = 113
    integer, parameter :: eModalSensor        = 114
    integer, parameter :: eJumpSensor         = 115

    ! Enum: type of polynomial expansion
    integer, parameter :: eNodal = 120
    integer, parameter :: eModal = 121

    ! Enum: export file extensions
    integer, parameter :: eGnuplot = 130
    integer, parameter :: eTecPlot = 131

    ! Enum: boundaries
    integer, parameter :: eLeft  = 1
    integer, parameter :: eRight = 2

    ! Flag: empty element
    integer, parameter :: fNone = -1

end module Constants
