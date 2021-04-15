!*******************************************************************************
!  MODULE: TimeIntegratorClass
!
!> @author
!> Andres Mateo
!
!> @brief
!> Definition of the 'TimeIntegrator_t' derived type.
!*******************************************************************************

module TimeIntegratorClass

    use Constants
    use Utilities
    use ExceptionsAndMessages
    use ElementClass
    use PDEclass
    use ExplicitMethods
    use TruncationErrorClass
    use Physics
    use SensorClass
    use PrinterClass
    !use hpAdaptation

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public types
    public :: TimeIntegrator_t

    ! Explicitly define public variables
    public :: TimeIntegrator

!···············································································
!> @class TimeIntegrator_t
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Handles the time integration of the equations.
!···············································································
    type TimeIntegrator_t
        integer  :: method       = eInvalid  !< time-stepping algorithm
        integer  :: printInt     = 0         !< time steps between prints
        integer  :: numSteps     = 0         !< no. of steps so far
        integer  :: dynStep      = 0         !< time steps between dyn. adapts.
        integer  :: firstDynStep = 0         !< first step for dynamic adapt.
        integer  :: lastDynStep  = 0         !< last step for dynamic adapt.
        integer  :: aViscSteps   = 0         !< no. of steps between vis. calc.
        integer  :: aViscLowStep = 0         !< first step of art. viscosity
        integer  :: aViscTopStep = 0         !< last step of art. viscosity
        integer  :: cflSteps     = 0         !< no. of steps between cfl calc.
        logical  :: completed    = .false.   !< .true. when the end is reached
        logical  :: stopNaN      = .true.    !< .true. to avoid NaN values
        real(wp) :: maxRes       = 0.0_wp    !< max. value of the residual
        real(wp) :: tInit        = 0.0_wp    !< initial time of the simulation
        real(wp) :: tFinal       = 0.0_wp    !< final time of the simulation
        real(wp) :: dt           = 0.0_wp    !< fixed time step
        real(wp) :: t            = 0.0_wp    !< 'present' time of the sim.
        procedure(Integrator_Int), &
            pointer, nopass, private :: &
                integrator => null()         !< Pointer to the int. routine
    contains
        procedure :: construct      => TimeIntegrator_constructor
        procedure :: stepInTime     => TimeIntegrator_time_step
        procedure :: integrateUntil => TimeIntegrator_integrate_until
        procedure :: integrate      => TimeIntegrator_integrate
        final     :: TimeIntegrator_destructor
    end type TimeIntegrator_t

    abstract interface
        subroutine Integrator_Int(t, dt)
            import wp
            real(wp), intent(in) :: t
            real(wp), intent(in) :: dt
        end subroutine Integrator_Int
    end interface

    type(TimeIntegrator_t) :: TimeIntegrator

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Constructor of the 'TimeIntegrator_t' derived type.
!
!> @param[in]  maxResidual       max. value of the residual (steady cases)
!> @param[in]  t0                initial time of the simulation
!> @param[in]  tf                final time of the simulation
!> @param[in]  dt                time step during the simulation
!> @param[in]  intMethod         time integration method chosen by the user
!> @param[in]  printStep         number of time steps between prints
!> @param[in]  dynamicStep       no. of time steps between dynamic adaptations
!> @param[in]  firstDynamicStep  first step where dynamic adapt. is performed
!> @param[in]  lastDynamicStep   last step where dynamic adapt. is performed
!> @param[in]  artViscSteps      no. of steps between art. visc. calculations
!> @param[in]  viscWindow        window of time steps where a. visc. is active
!> @param[in]  cflStep           no. of steps between CFL number calculations
!> @param[in]  stopAtNaN         .true. if the integration should avoid NaNs
!···············································································
subroutine timeIntegrator_constructor(this, maxResidual, t0, tf, dt,     &
                                      intMethod, printStep, dynamicStep, &
                                      firstDynamicStep, lastDynamicStep, &
                                      artViscSteps, viscWindow, cflStep, &
                                      stopAtNaN)
    !* Arguments *!
    real(wp),     intent(in) :: maxResidual
    real(wp),     intent(in) :: t0
    real(wp),     intent(in) :: tf
    real(wp),     intent(in) :: dt
    integer,      intent(in) :: intMethod
    integer,      intent(in) :: printStep
    integer,      intent(in) :: dynamicStep
    integer,      intent(in) :: firstDynamicStep
    integer,      intent(in) :: lastDynamicStep
    integer,      intent(in) :: artViscSteps
    integer,      intent(in) :: viscWindow(2)
    integer,      intent(in) :: cflStep
    logical,      intent(in) :: stopAtNaN
    ! Derived types
    class(TimeIntegrator_t), intent(inout) :: this

    ! Set some time constants
    this%completed    = .false.
    this%stopNaN      = stopAtNaN
    this%maxRes       = maxResidual
    this%tInit        = t0
    this%tFinal       = tf
    this%dt           = dt
    this%printInt     = printStep
    this%dynStep      = dynamicStep
    this%firstDynStep = firstDynamicStep
    this%lastDynStep  = lastDynamicStep
    this%aViscSteps   = artViscSteps
    this%aViscLowStep = viscWindow(1)
    this%aViscTopStep = viscWindow(2)
    this%cflSteps     = cflStep

    ! Select the time integration scheme
    select case (intMethod)

    case (eRK3)
        this%integrator => DGstepByRK3

    case (eRK5)
        this%integrator => DGstepByRK5

    case default
        call printError("TimeIntegratorClass.f90", &
                        "A valid time integration scheme must be provided.")

    end select

end subroutine timeIntegrator_constructor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Integrates the equations one step in time.
!
!> @param[out]  realDt  $\Delta t$ used in this time step
!> @param[in]   dt      size of the time step (optional)
!···············································································
subroutine timeIntegrator_time_step(this, realDt, dt)
    !* Arguments *!
    real(wp),           intent(out) :: realDt
    real(wp), optional, intent(in)  :: dt
    ! Derived types
    class(TimeIntegrator_t), intent(inout) :: this

    !* Local variables *!
    character(len=CHAR_LEN) :: errorMsg
    integer                 :: i
    type(Elem_t), pointer   :: elem

    ! Check optional arguments
    if (present(dt)) then
        realDt = dt
    else
        ! Do not go over 'tFinal'
        if (this%tFinal > 0.0_wp) then
            realDt = merge(this%dt, this%tFinal-this%t, &
                           this%tFinal-this%t > this%dt)
        ! When tFinal <= 0, we suppose a steady simulation
        else
            realDt = this%dt
        end if
    end if

    ! Step forward
    call this%integrator(this%t, realDt)

    ! Update 'TimeIntegrator_t' attributes
    this%t        = this%t + realDt
    this%numSteps = this%numSteps + 1

    ! Look for unconverged values
    call PDE%mesh%elems%reset_last(0)
    do i = 1, PDE%mesh%elems%size()

        elem => PDE%mesh%elems%next()

        ! Error handling (save unconverged solution and throw the error)
        if (this%stopNaN .and. any(isnan(elem%Phi))) then
            call Printer%saveMesh(this%t)
            write(errorMsg, '(a,f0.5)') "NaN found at time t = ", this%t
            call printError("TimeIntegration", errorMsg)
        end if

    end do

    nullify(elem)

    ! If there were no errors, update everything (1 more "useless" computation)
    call PDE%globalTimeDerivative(this%t)

    ! Check if the final time was reached (unsteady case)
    if (this%tFinal > 0.0_wp) then
        if (this%t >= this%tFinal) then
            this%completed = .true.
        else
            this%completed = .false.
        end if

        ! ...or the residual is below the threshold (steady case)
    else
        call PDE%computeMaxResidual()
        if (PDE%maxResidual < this%maxRes) then
            this%completed = .true.
        else
            this%completed = .false.
        end if
    end if

end subroutine timeIntegrator_time_step

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Integrates the equations until the instant specified by 'tf' is reached.
!
!> @param[in]  tf  time instant until which the integration is performed
!···············································································
subroutine timeIntegrator_integrate_until(this, tf)
    !* Arguments *!
    real(wp), intent(in) :: tf
    ! Derived types
    class(TimeIntegrator_t), intent(inout) :: this

    !* Local variables *!
    real(wp) :: dt

    ! Integration loop
    do while (this%t < tf)

        ! Determine the maximum posible 'dt'
        dt = min(this%dt, tf-this%t)

        ! Make one time step
        call this%stepInTime(dt)

    end do

end subroutine timeIntegrator_integrate_until

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Integrates the equations until final time specified at creation.
!
!> @param[in]  TEmap  .true. if the TE map has to be computed
!···············································································
subroutine timeIntegrator_integrate(this, TEmap)
    !* Arguments *!
    logical, intent(in) :: TEmap
    ! Derived types
    class(TimeIntegrator_t), intent(inout) :: this

    !* Local variables *!
    character(len=CHAR_LEN) :: infoMsg
    real(wp)                :: dt
    logical                 :: adapt
    logical                 :: sensorUpdated
    logical                 :: viscUpdated
    integer                 :: i
    type(Elem_t), pointer   :: elem

    ! Initialisation
    adapt         = .false.
    sensorUpdated = .false.

    ! Compute sensors before first save
    call PDE%globalTimeDerivative(this%t)
    call Sensor%updateMesh(this%t, scaled=.true.)

    ! Assign a fixed value for the artificial viscosity coeff. if requested
    if (Phys%Alpha(1) /= 0.0_wp) then

        if (this%aViscSteps <= 0) then

            call PDE%mesh%elems%reset_last(0)
            do i = 1, PDE%mesh%elems%size()

                elem => PDE%mesh%elems%next()
                elem%aVis = Phys%Alpha(1)

            end do

            nullify(elem)

            call printWarning("TimeIntegratorClass.f90", &
                "Artificial viscosity coeffs. have been set as constants.")

        else if (this%numSteps >= this%aViscLowStep) then

            call PDE%mesh%updateArtViscosity()

        end if

    end if

    !! Print initial mesh
    call Printer%saveMesh(this%t)

    !! Main loop of integration
    do while (.not. this%completed)

        ! Take step
        call this%stepInTime(dt)

        ! Update the sensor and the artificial viscous coefficients
        if (this%numSteps >= this%aViscLowStep .and. &
            this%numSteps <= this%aViscTopStep) then

            if (Phys%Alpha(1) /= 0.0_wp .and. this%aViscSteps > 0) then

                if (mod(this%numSteps, this%aViscSteps) == 0) then

                    call Sensor%updateMesh(this%t, scaled=.true.)
                    call PDE%mesh%updateArtViscosity()
                    sensorUpdated = .true.

                end if

            end if

        end if

        !if (pAdaptationType == eDynamicAdaptation) then

            !! Adapt if it is the time
            !if (this%numSteps == this%firstDynStep) then

                !adapt = .true.

            !else if (this%numSteps == this%lastDynStep) then

                !adapt = .true.

            !else if (this%numSteps > this%firstDynStep .and. &
                     !this%numSteps < this%lastDynStep) then

                !if (mod(this%numSteps-this%firstDynStep, this%dynStep)==0) then
                    !adapt = .true.
                !end if

            !end if

        !end if

        ! Computation of the CFL condition
        if (this%cflSteps > 0) then
            if (mod(this%numSteps, this%cflSteps) == 0) then

                call PDE%computeCFLnumber(dt)
                call PDE%computeMaxResidual()

                write(infoMsg, '(a,e15.8,2a,f0.5,a,i0,a,f0.5,a)')             &
                    "Max. residual: ", PDE%maxResidual, ". ",                 &
                    "Max. CFL: ", PDE%maxCFL, " in element ", PDE%maxCFLelem, &
                    " at time ", this%t, " s."
                call printInfo("TimeIntegratorClass.f90", infoMsg)

            end if
        end if

        ! Print if it is time for it (always before adaptation)
        if (this%printInt > 0) then
            if (mod(this%numSteps, this%printInt) == 0) then

                if (.not. sensorUpdated) then
                    call Sensor%updateMesh(this%t, scaled=.false.)
                    sensorUpdated = .true.
                end if

                call Printer%saveMesh(this%t)

            end if
        end if

        ! Perform h-adaptation
        if (adapt) then

            ! Recall that the sensor might have been calculated before
            if (.not. sensorUpdated) then
                call Sensor%updateMesh(this%t, scaled=.false.)
                sensorUpdated = .true.
            end if

            ! Adapt
            !call hAdaptation(mesh)
            adapt = .false.

            ! Notify
            write(infoMsg, '(a,a,f0.5)') "Dynamic adaptation performed at ", &
                                         "time ", this%t
            call printInfo("TimeIntegratorClass.f90", infoMsg)

        end if

        ! Reset the flags
        sensorUpdated = .false.
        viscUpdated   = .false.

    end do

    !! Print final solution after reaching the end
    if (this%printInt <= 0) then

        call Sensor%updateMesh(this%t, scaled=.false.)
        call Printer%saveMesh(this%t)

    else if (mod(this%numSteps, this%printInt) /= 0) then

        call Sensor%updateMesh(this%t, scaled=.false.)
        call Printer%saveMesh(this%t)

    end if

    !! Truncation error map
    if (TEmap) then

        ! Compute the truncation error map
        call TruncError%mapError(this%t)

        ! Print the mesh with the error map
        call Printer%saveTEmap(this%t)

        call printInfo("TimeIntegratorClass.f90", &
                       "Truncation error map computed.")

    end if

end subroutine timeIntegrator_integrate

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Finalises the timeIntegrator object.
!···············································································
subroutine TimeIntegrator_destructor(this)
    !* Arguments *!
    type(TimeIntegrator_t), intent(inout) :: this

    if (associated(this%integrator)) then
        nullify(this%integrator)
    end if

end subroutine TimeIntegrator_destructor

end module TimeIntegratorClass
