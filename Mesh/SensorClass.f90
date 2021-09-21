!*******************************************************************************
!  MODULE: SensorsClass
!
!> @author
!> Andres Mateo
!
!> @brief
!> Implementation of different shock-capturing sensors.
!*******************************************************************************

module SensorClass

    use Constants
    use Utilities
    use PDEclass
    use ElementClass
    use TruncationErrorClass
    use ExceptionsAndMessages
    use Physics

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public types
    public :: Sensor_t

    ! Explicitly define public variables
    public :: Sensor

    ! Private variables
    real(wp), parameter :: MIN_SENSOR = -100.0_wp

!···············································································
!> @class Sensor_t
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Wraps the different implementations of the shock-capturing sensors.
!···············································································
    type Sensor_t
        integer  :: SensorType
        integer  :: SecSensorType
        integer  :: SensedVariable
        real(wp) :: mainMin
        real(wp) :: mainMax
        real(wp) :: mainS0
        real(wp) :: mainDelta
        real(wp) :: secMin
        real(wp) :: secMax
        real(wp) :: secS0
        real(wp) :: secDelta
        procedure(Sensor_Int), private, &
            pointer, nopass :: mainSensor => null()
        procedure(Sensor_Int), private, &
            pointer, nopass :: secSensor  => null()
    contains
        procedure :: construct  => Sensor_constructor
        procedure :: sense      => Sensor_compute_sensor
        procedure :: updateMesh => Sensor_update_sensor_values
        final     :: Sensor_destructor
    end type Sensor_t

    abstract interface
        function Sensor_Int(elem, time)
            import wp, Elem_t
            type(Elem_t), intent(inout) :: elem
            real(wp),     intent(in)    :: time
            real(wp) :: Sensor_Int
        end function Sensor_Int
    end interface

    type(Sensor_t) :: Sensor

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initialise the sensor derived type.
!
!> @param[in]  sensType     sensor type to be used
!> @param[in]  secSensType  secondary sensor type to be used
!> @param[in]  sensVar      variable used to compute the sensor
!> @param[in]  mainBottom   lower limit of the main sensor ramp
!> @param[in]  mainTop      upper limit of the main sensor ramp
!> @param[in]  secBottom    lower limit of the secondary sensor ramp
!> @param[in]  secTop       upper limit of the secondary sensor ramp
!···············································································
subroutine Sensor_constructor(this, sensType, secSensType, sensVar, &
                              mainBottom, mainTop, secBottom, secTop)
    !* Arguments *!
    integer,  intent(in) :: sensType
    integer,  intent(in) :: secSensType
    integer,  intent(in) :: sensVar
    real(wp), intent(in) :: mainBottom
    real(wp), intent(in) :: mainTop
    real(wp), intent(in) :: secBottom
    real(wp), intent(in) :: secTop
    ! Derived types
    class(Sensor_t), intent(inout) :: this

    ! Set the global flags
    this%SensorType     = sensType
    this%SecSensorType  = secSensType
    this%SensedVariable = sensVar

    this%mainMin   = mainBottom
    this%mainMax   = mainTop
    this%mainS0    = (mainTop + mainBottom) / 2.0_wp
    this%mainDelta = mainTop - mainBottom
    this%secMin    = secBottom
    this%secMax    = secTop
    this%secS0     = (secTop + secBottom) / 2.0_wp
    this%secDelta  = secTop - secBottom

    ! Point to the selected sensor function
    select case (sensType)

    case (eTruncationError)
        this%mainSensor => TEsensor

    case (eModalSensor)
        this%mainSensor => LegendreSensor

    case (eDensitySensor)
        this%mainSensor => densitySensor

    case (eJumpSensor)
        this%mainSensor => jumpSensor

    case default
        call printError("SensorClass.f90", &
                        "The selected sensor is not available.")

    end select

    ! Check the selected secondary sensor
    select case (secSensType)

    case (eTruncationError)
        this%secSensor => TEsensor

    case (eModalSensor)
        this%secSensor => LegendreSensor

    case (eDensitySensor)
        this%secSensor => densitySensor

    case (eJumpSensor)
        this%secSensor => jumpSensor

    case default
        call printError("SensorClass.f90", &
                        "The selected secondary sensor is not available.")

    end select

end subroutine Sensor_constructor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Wrapper over the different sensor implementations. It applies the main sensor
!> to the elements with \f$ P>0 \f$ and the secondary one, to those with
!> \f$ P=0 \f$.
!
!> @param[in]   elem       element where the sensor is computed
!> @param[in]   time       time instant
!> @param[out]  sensedVal  value of the sensor
!> @param[out]  scaledVal  value of the sensor after scaling (optional)
!···············································································
subroutine Sensor_compute_sensor(this, elem, time, sensedVal, scaledVal)
    !* Arguments *!
    real(wp),           intent(in)  :: time
    real(wp),           intent(out) :: sensedVal
    real(wp), optional, intent(out) :: scaledVal
    ! Derived types
    class(Sensor_t), intent(in)    :: this
    type(Elem_t),    intent(inout) :: elem

    ! Call the main sensor on elements with P > 1
    if (elem%std%n > 1) then
        sensedVal = this%mainSensor(elem, time)

        if (present(scaledVal)) then
            scaledVal = sinRamp(sensedVal-this%mainS0, this%mainDelta)
        end if

    ! Apply the secondary sensor to the constant elements
    else
        sensedVal = this%secSensor(elem, time)

        if (present(scaledVal)) then
            scaledVal = sinRamp(sensedVal-this%secS0, this%secDelta)
        end if

    end if

end subroutine Sensor_compute_sensor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Updates the value of the sensor in all the elements of the mesh.
!
!> @param[in]  time    time instant
!> @param[in]  scaled  .true. if the scaled value must be computed
!···············································································
subroutine Sensor_update_sensor_values(this, time, scaled)
    !* Arguments *!
    real(wp), intent(in) :: time
    logical,  intent(in) :: scaled
    ! Derived types
    class(Sensor_t), intent(in) :: this

    !* Local variables *!
    integer :: i
    type(Elem_t), pointer :: elem

    ! Loop over all the elements
    call PDE%mesh%elems%reset_last(0)
    do i = 1, PDE%mesh%elems%size()

        elem => PDE%mesh%elems%next()

        if (scaled) then
            call this%sense(elem, time, elem%sens, elem%aVis)
        else
            call this%sense(elem, time, elem%sens)
        end if

        ! Mark elements detected by the sensor
        if (elem%std%n > 1) then
            if (elem%sens >= this%mainMin) then
                elem%sensed = .true.
            else
                elem%sensed = .false.
            end if

        else
            if (elem%sens >= this%secMin) then
                elem%sensed = .true.
            else
                elem%sensed = .false.
            end if

        end if

    end do

    nullify(elem)

end subroutine Sensor_update_sensor_values

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Sensor based on the truncation error estimator from the solution residual.
!
!> @param[inout]  elem  element where the sensor is computed
!> @param[in]     time  time instant
!···············································································
function TEsensor(elem, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    type(Elem_t), intent(inout) :: elem

    !* Return values *!
    real(wp) :: TEsensor
    integer  :: orderDec

        ! Defaults to 0 if not the main sensor
        if (elem%std%n == 1) then

            TEsensor = 0.0_wp

        ! Maximize the difference between the coarse and the fine mesh
        else

            select case (elem%std%n)
            case (2:)
                orderDec = 1
            case default
                orderDec = 0   ! Should be useless
            end select

            ! Temporarily store the error in the output variable
            TEsensor = TruncError%estimate(elem, orderDec, time, &
                                           Sensor%SensedVariable)

        end if

    ! Apply the decimal logarithm to the error
    if (TEsensor == 0.0_wp) then
        TEsensor = MIN_SENSOR  ! Avoid -inf values
    else
        TEsensor = log10(TEsensor)
    end if

end function TEsensor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Sensor based on the norm of the higher mode of the expansion, according to
!> Huerta, A., Casoni, E., & Peraire, J. (2012). A simple shock-capturing
!> technique for high-order discontinuous Galerkin methods. International
!> Journal for Numerical Methods in Fluids, 69(10), 1614–1632.
!
!> @param[in]  elem  element where the sensor is computed
!> @param[in]  time  time instant
!···············································································
function LegendreSensor(elem, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    type(Elem_t), intent(inout) :: elem

    !* Return values *!
    real(wp) :: LegendreSensor

    !* Local variables *!
    integer  :: i
    real(wp) :: sensedVariable
    real(wp) :: num
    real(wp) :: den

    ! Shortcut calculations if possible
    if (elem%std%n <= 1) then
        LegendreSensor = MIN_SENSOR
    end if

    ! Switch to modal space...
    call elem%toLegendre()

    ! ...and compute the error
    den = 0.0_wp
    do i = 1, elem%std%n
        sensedVariable = getSensedVar(elem%Phi(i,:))
        den = den + elem%std%modeN2(i) * sensedVariable**2
    end do
    sensedVariable = getSensedVar(elem%Phi(elem%std%n,:))
    num = elem%std%modeN2(elem%std%n) * sensedVariable**2

    ! Final expression of the sensor
    if (num == 0.0_wp) then
        LegendreSensor = MIN_SENSOR
    else
        LegendreSensor = log10(num/den)
    end if

    ! Switch back to nodal space
    call elem%fromLegendre()

end function LegendreSensor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Sensor based on the norm of the density gradient.
!
!> @param[in]  elem  element where the sensor is computed
!> @param[in]  time  time instant
!···············································································
function densitySensor(elem, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    type(Elem_t), intent(inout) :: elem

    !* Return value *!
    real(wp) :: densitySensor

    !* Local variables *!
    real(wp) :: gradRho2
    integer  :: i


    densitySensor = 0.0_wp

    do i = 1, elem%std%n

        if (Phys%WithEntropyVars) then
            gradRho2 = sum(elem%Phi(i,:), elem%Grad(i,:))**2
        else
            gradRho2 = elem%Grad(i,IRHO)**2
        end if

        densitySensor = densitySensor + gradRho2 * elem%std%w(i)

    end do

    densitySensor = sqrt(densitySensor * elem%jac)

end function densitySensor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Sensor based on the jump between elements: Jameson, A., & Mavriplis, D.
!> (1986). Finite Volume Dolution of the Two-Dimensional Euler Equations on a
!> Regular Triangular Mesh. AIAA Journal, 24(4), 611–618.
!
!> @param[in]  elem  element where the sensor is computed
!> @param[in]  time  time instant
!···············································································
function jumpSensor(elem, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    type(Elem_t), intent(inout) :: elem

    !* Return values *!
    real(wp) :: jumpSensor

    !* Local variables *!
    real(wp) :: sVarLL
    real(wp) :: sVarLR
    real(wp) :: sVarRL
    real(wp) :: sVarRR
    real(wp) :: num
    real(wp) :: den

    !---- Begin associate ----!
    associate(faceLeft  => PDE%mesh%faces%at(elem%faceLeft), &
              faceRight => PDE%mesh%faces%at(elem%faceRight))

        sVarLL = getSensedVar(faceLeft%PhiL)
        sVarLR = getSensedVar(faceLeft%PhiR)
        sVarRL = getSensedVar(faceRight%PhiL)
        sVarRR = getSensedVar(faceRight%PhiR)

        num = abs(sVarLL + sVarRR - sVarLR - sVarRL)
        den = abs(sVarLL + sVarRR + sVarLR + sVarRL)

        ! Avoid -inf result from log10(0)
        if (num == 0.0_wp) then
            jumpSensor = MIN_SENSOR
        else
            jumpSensor = log10(num / den)
        end if

    end associate
    !----- End associate -----!

end function jumpSensor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Compute the value of the sensed variable according to 'pSensedVar'.
!
!> @param[in]  Phi  function values
!
!> @return          variable to be used as an input for the sensor
!···············································································
function getSensedVar(Phi) result(var)
    !* Arguments *!
    real(wp), intent(in) :: Phi(:)

    !* Return values *!
    real(wp) :: var

    if (Sensor%SensedVariable > 0) then
        var = Phi(Sensor%SensedVariable)

    else
        if (Phi(IRHO) /= 0.0_wp) then
            var = getPressure(Phi) * Phi(IRHO)
        else
            var = 0.0_wp
        end if

    end if

end function getSensedVar

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Finalises the pointers defined in 'Sensors.f90'.
!···············································································
subroutine Sensor_destructor(this)
    !* Arguments *!
    type(Sensor_t), intent(inout) :: this

    if (associated(this%mainSensor)) then
        nullify(this%mainSensor)
    end if

    if (associated(this%secSensor)) then
        nullify(this%secSensor)
    end if

end subroutine Sensor_destructor

end module SensorClass
