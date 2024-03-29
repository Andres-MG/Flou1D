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
    use MeshClass
    use ElementClass
    use TruncationErrorClass
    use ExceptionsAndMessages
    use Physics
    use Clustering, only: kMeans, GMM

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
        integer  :: nClusters
        real(wp) :: mainMin
        real(wp) :: mainMax
        real(wp) :: mainS0
        real(wp) :: mainDelta
        real(wp) :: secMin
        real(wp) :: secMax
        real(wp) :: secS0
        real(wp) :: secDelta
        procedure(SensorGlobal_Int), private, pointer :: globalSensor => null()
        procedure(Sensor_Int),       private, pointer :: mainSensor   => null()
        procedure(Sensor_Int),       private, pointer :: secSensor    => null()
    contains
        procedure :: construct  => Sensor_constructor
        procedure :: sense      => Sensor_compute_sensor
        procedure :: updateMesh => Sensor_update_sensor_values
        final     :: Sensor_destructor
    end type Sensor_t

    abstract interface
        subroutine SensorGlobal_Int(self, mesh, time)
            import wp, Sensor_t, Mesh_t
            class(Sensor_t), intent(in)    :: self
            type(Mesh_t),    intent(inout) :: mesh
            real(wp),        intent(in)    :: time
        end subroutine SensorGlobal_Int
        function Sensor_Int(self, elem, time)
            import wp, Sensor_t, Elem_t
            class(Sensor_t), intent(in)    :: self
            type(Elem_t),    intent(inout) :: elem
            real(wp),        intent(in)    :: time
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
!> @param[in]  numClusters  number of clusters to used by clustering algorithms
!> @param[in]  mainBottom   lower limit of the main sensor ramp
!> @param[in]  mainTop      upper limit of the main sensor ramp
!> @param[in]  secBottom    lower limit of the secondary sensor ramp
!> @param[in]  secTop       upper limit of the secondary sensor ramp
!···············································································
subroutine Sensor_constructor(this, sensType, secSensType, sensVar, numClusters, &
                              mainBottom, mainTop, secBottom, secTop)
    !* Arguments *!
    integer,  intent(in) :: sensType
    integer,  intent(in) :: secSensType
    integer,  intent(in) :: sensVar
    integer,  intent(in) :: numClusters
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
    this%nClusters      = numClusters

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

    case (eAliasingSensor)
        this%mainSensor => aliasingSensor

    case (eTruncationError)
        this%mainSensor => TEsensor

    case (eModalSensor)
        this%mainSensor => LegendreSensor

    case (eDensitySensor)
        this%mainSensor => densitySensor

    case (eJumpSensor)
        this%mainSensor => jumpSensor

    case (eKmeansSensor)
        this%globalSensor => kMeansSensor

    case (eGMMSensor)
        this%globalSensor => GMMSensor

    case default
        call printError("SensorClass.f90", &
                        "The selected sensor is not available.")

    end select

    ! Check the selected secondary sensor
    select case (secSensType)

    case (eAliasingSensor)
        this%secSensor => aliasingSensor

    case (eTruncationError)
        this%secSensor => TEsensor

    case (eModalSensor)
        this%secSensor => LegendreSensor

    case (eDensitySensor)
        this%secSensor => densitySensor

    case (eJumpSensor)
        this%secSensor => jumpSensor

    case (eKmeansSensor)
        call printError("SensorClass.f90", &
                        "The k-means sensor cannot be selected as a secondary sensor.")

    case (eGMMSensor)
        call printError("SensorClass.f90", &
                        "The GMM sensor cannot be selected as a secondary sensor.")

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
!> Map the "raw" sensor value to the [0, 1] region.
!
!> @param[in]  elem       element where the sensor was computed
!> @param[in]  sensedVal  value of the sensor
!···············································································
function Sensor_rescale_element(this, elem, sensedVal) result(scaledVal)
    !* Arguments *!
    real(wp), intent(in) :: sensedVal
    ! Derived types
    type(Sensor_t), intent(in) :: this
    type(Elem_t),   intent(in) :: elem

    !* Return values *!
    real(wp) :: scaledVal

    if (elem%std%n > 1) then
        scaledVal = sinRamp(sensedVal-this%mainS0, this%mainDelta)
    else
        scaledVal = sinRamp(sensedVal-this%secS0, this%secDelta)
    end if

end function Sensor_rescale_element

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
!···············································································
function Sensor_compute_sensor(this, elem, time) result(sensedVal)
    !* Arguments *!
    real(wp), intent(in)  :: time
    ! Derived types
    class(Sensor_t), intent(in)    :: this
    type(Elem_t),    intent(inout) :: elem

    !* Return values *!
    real(wp) :: sensedVal

    if (elem%std%n > 1 .and. associated(this%mainSensor)) then
        sensedVal = this%mainSensor(elem, time)
    else if (elem%std%n == 1 .and. associated(this%secSensor)) then
        sensedVal = this%secSensor(elem, time)
    end if

end function Sensor_compute_sensor

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
    integer  :: i
    real(wp) :: rescaled
    type(Elem_t), pointer :: elem

    ! Global sensor
    if (associated(this%globalSensor)) then
        call this%globalSensor(PDE%mesh, time)
    end if

    ! Loop over all the elements
    call PDE%mesh%elems%reset_last(0)
    do i = 1, PDE%mesh%elems%size()

        elem => PDE%mesh%elems%next()

        if (associated(this%globalSensor)) then
            rescaled = elem%sens
        else
            elem%sens = this%sense(elem, time)
            rescaled = Sensor_rescale_element(this, elem, elem%sens)
        end if
        if (scaled) elem%aVis = rescaled

        ! Mark elements detected by the sensor
        elem%sensed    = merge(.true., .false., rescaled >  0.0_wp)
        elem%saturated = merge(.true., .false., rescaled >= 1.0_wp)

    end do

    nullify(elem)

end subroutine Sensor_update_sensor_values

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Sensor based on the aliasing error.
!
!> @param[inout]  elem  element where the sensor is computed
!> @param[in]     time  time instant
!···············································································
function aliasingSensor(this, elem, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(Sensor_t), intent(in) :: this
    type(Elem_t), intent(inout) :: elem

    !* Return values *!
    real(wp) :: aliasingSensor

    !* Local variables *!
    integer               :: i, j
    integer               :: n
    real(wp), allocatable :: Fs(:,:,:)
    real(wp), allocatable :: aliasing(:,:)

    n = elem%std%n

    allocate(Fs(n, n, NEQS))
    do i = 1, n
        Fs(i,i,:) = EulerFlux(elem%Phi(i,:))
        do j = i+1, n
            Fs(i,j,:) = TwoPointFlux(elem%Phi(i,:), elem%Phi(j,:))
            Fs(j,i,:) = Fs(i,j,:)
        end do
    end do

    allocate(aliasing(n, NEQS), source=0.0_wp)
    do j = 1, n
        do i = 1, n
            aliasing(i,:) = aliasing(i,:) + elem%std%D(i,j) * (2.0_wp*Fs(i,j,:) - Fs(j,j,:))
            aliasing(i,:) = abs(aliasing(i,:))
        end do
    end do

    ! Aliasing error estimation
    aliasingSensor = maxval(aliasing)
    if (aliasingSensor == 0.0_wp) then
        aliasingSensor = MIN_SENSOR
    else
        aliasingSensor = log10(aliasingSensor)
    end if

end function aliasingSensor

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
function TEsensor(this, elem, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(Sensor_t), intent(in) :: this
    type(Elem_t), intent(inout) :: elem

    !* Return values *!
    real(wp) :: TEsensor

    !* Local variables *!
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
function LegendreSensor(this, elem, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(Sensor_t), intent(in) :: this
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
function densitySensor(this, elem, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(Sensor_t), intent(in) :: this
    type(Elem_t), intent(inout) :: elem

    !* Return value *!
    real(wp) :: densitySensor

    !* Local variables *!
    real(wp) :: gradRho2
    integer  :: i


    densitySensor = 0.0_wp

    do i = 1, elem%std%n

        if (Phys%WithEntropyVars) then
            gradRho2 = sum(elem%Phi(i,:) * elem%Grad(i,:))**2
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
function jumpSensor(this, elem, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(Sensor_t), intent(in) :: this
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
!> Sensor using the kMeans clustering algorithm
!> (based on the degree thesis of Carlos Piqueras).
!
!> @param[inout]  mesh  mesh with all the elements to be sensed
!> @param[in]     time  time instant
!···············································································
subroutine kMeansSensor(this, mesh, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(Sensor_t), intent(in)    :: this
    type(Mesh_t),    intent(inout) :: mesh

    !* Local variables *!
    integer :: ie, i
    type(Elem_t), pointer :: elem

    integer  :: nnodes
    integer  :: cnt
    integer  :: n
    integer  :: info
    real(wp) :: u2
    real(wp) :: higherCluster
    real(wp) :: centroids(2,this%nClusters)
    real(wp), allocatable :: derivs(:,:)
    integer,  allocatable :: clusters(:)

    ! Compute derivatives
    nnodes = 0
    call mesh%elems%reset_last(0)
    do ie = 1, mesh%elems%size()
        elem => mesh%elems%next()
        nnodes = nnodes + elem%std%n
    end do
    allocate(derivs(2,nnodes))  ! (rhox, px)
    allocate(clusters(nnodes))

    cnt = 0
    call mesh%elems%reset_last(0)
    do ie = 1, mesh%elems%size()
        elem => mesh%elems%next()
        n = elem%std%n
        do i = 1, n
            cnt = cnt + 1
            u2 = (elem%Phi(i,IRHOU) / elem%Phi(i,IRHO))**2
            derivs(1,cnt) = elem%PhiD(i,IRHO)
            derivs(2,cnt) = Phys%GammaMinusOne * (elem%PhiD(i,IRHOE) - &
                            0.5_wp * u2 * elem%PhiD(i,IRHO) -          &
                            elem%Phi(i,IRHOU) * elem%PhiD(i,IRHOU))
        end do
    end do

    ! Rescale the derivatives for clustering
    call getClusterVariables(this%nClusters, derivs, centroids)

    ! Clustering
    call kMeans(this%nClusters, derivs, centroids, clusters, info)
    if (info < 0) then
        call printWarning("SensorClass.f90", "kMeans did not converge")
    end if

    ! Assign sensor value according to clusters
    cnt = 0
    call mesh%elems%reset_last(0)
    do ie = 1, mesh%elems%size()
        elem => mesh%elems%next()
        n = elem%std%n
        higherCluster = maxval(clusters(cnt+1:cnt+n))
        elem%sens = real(higherCluster - 1, kind=wp) / (this%nClusters - 1)
        cnt = cnt + n
    end do

    nullify(elem)

end subroutine kMeansSensor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Sensor using the Gaussian Mixture Model
!> (based on the degree thesis of Carlos Piqueras).
!
!> @param[inout]  mesh  mesh with all the elements to be sensed
!> @param[in]     time  time instant
!···············································································
subroutine GMMSensor(this, mesh, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(Sensor_t), intent(in)    :: this
    type(Mesh_t),    intent(inout) :: mesh

    !* Local variables *!
    integer :: ie, i
    type(Elem_t), pointer :: elem

    integer  :: nclusters
    integer  :: nnodes
    integer  :: cnt
    integer  :: n
    integer  :: info
    real(wp) :: u2
    real(wp) :: higherCluster
    real(wp) :: centroids(2,this%nClusters)
    real(wp), allocatable :: derivs(:,:)
    integer,  allocatable :: clusters(:)


    ! Compute derivatives
    nnodes = 0
    call mesh%elems%reset_last(0)
    do ie = 1, mesh%elems%size()
        elem => mesh%elems%next()
        nnodes = nnodes + elem%std%n
    end do
    allocate(derivs(2,nnodes))  ! (rhox, px)
    allocate(clusters(nnodes))

    cnt = 0
    call mesh%elems%reset_last(0)
    do ie = 1, mesh%elems%size()
        elem => mesh%elems%next()
        n = elem%std%n
        do i = 1, n
            cnt = cnt + 1
            u2 = (elem%Phi(i,IRHOU) / elem%Phi(i,IRHO))**2
            derivs(1,cnt) = elem%PhiD(i,IRHO)
            derivs(2,cnt) = Phys%GammaMinusOne * (elem%PhiD(i,IRHOE) - &
                            0.5_wp * u2 * elem%PhiD(i,IRHO) -          &
                            elem%Phi(i,IRHOU) * elem%PhiD(i,IRHOU))
        end do
    end do

    ! Rescale the derivatives for clustering
    call getClusterVariables(this%nClusters, derivs, centroids)

    ! Clustering
    call GMM(this%nClusters, nclusters, derivs, centroids, clusters, info)
    if (info == -1) then
        call printWarning("SensorClass.f90", "GMM used kMeans")
    elseif (info == -2) then
        call printWarning("SensorClass.f90", "GMM did not converge")
    end if

    ! Assign sensor value according to clusters
    cnt = 0
    call mesh%elems%reset_last(0)
    do ie = 1, mesh%elems%size()
        elem => mesh%elems%next()
        n = elem%std%n
        if (nclusters <= 1) then
            elem%sens = 0.0_wp
        else
            higherCluster = maxval(clusters(cnt+1:cnt+n))
            elem%sens = real(higherCluster - 1, kind=wp) / (nclusters - 1)
        end if
        cnt = cnt + n
    end do

    nullify(elem)

end subroutine GMMSensor

subroutine getClusterVariables(nclusters, x, xavg)
    !* Arguments *!
    integer,  intent(in)    :: nclusters
    real(wp), intent(inout) :: x(:,:)
    real(wp), intent(out)   :: xavg(:,:)

    !* Local variables *!
    integer  :: i
    real(wp) :: minimum(2)
    real(wp) :: maximum(2)
    real(wp) :: diff(2)


    x = abs(x)
    minimum = minval(x, dim=2)
    maximum = maxval(x, dim=2)

    if (almost_zero(maximum(1) - minimum(1))) then
        if (maximum(1) > 0.0_wp) then
            x(1,:) = 1.0_wp
            minimum(1) = 1.0_wp
            maximum(1) = 1.0_wp
        else
            x(1,:) = 0.0_wp
            minimum(1) = 0.0_wp
            maximum(1) = 0.0_wp
        end if
    else
        x(1,:) = (x(1,:) - minimum(1)) / (maximum(1) - minimum(1))
        minimum(1) = 0.0_wp
        maximum(1) = 1.0_wp
    end if

    if (almost_zero(maximum(2) - minimum(2))) then
        if (maximum(2) > 0.0_wp) then
            x(2,:) = 20.0_wp
            minimum(2) = 20.0_wp
            maximum(2) = 20.0_wp
        else
            x(2,:) = 0.0_wp
            minimum(2) = 0.0_wp
            maximum(2) = 0.0_wp
        end if
    else
        x(2,:) = (x(2,:) - minimum(2)) / (maximum(2) - minimum(2)) * 20.0_wp
        minimum(2) = 0.0_wp
        maximum(2) = 20.0_wp
    end if

    ! Make clusters
    xavg(:,1) = minimum
    xavg(:,nClusters) = maximum
    diff = (maximum - minimum) / (nClusters - 1)
    do i = 2, nClusters-1
        xavg(:,i) = minimum + diff * (i - 1)
    end do

end subroutine getClusterVariables

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
        var = getPressure(Phi) * Phi(IRHO)

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
