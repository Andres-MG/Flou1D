program main

    use Constants, only: wp
    use PrinterClass
    use PDEclass
    use Physics
    use TimeIntegratorClass
    use TruncationErrorClass
    use SensorClass
    !use hpAdaptation
    use Setup_params

    implicit none

    !* Some local variables *!
    real(wp) :: startTime
    real(wp) :: endTime

    ! Printer
    call Printer%construct(FILENAME, SAVERES)

    ! Physics-related parameters
    call Phys%construct(GAMMA, MU, PR, ALPHAMAX, ALPHA2BETA, ALPHA2LAMBDA, &
                        SVVPOW, ALPHASVV, ALPHA2BETASVV, ALPHA2LAMBDASVV,  &
                        SSFVBLEND)
    call initArtificialViscosity(ARTVISC, SVVTYPE)
    call initRiemannSolver(TWOPOINTTYPE, DISSTYPE)

    ! Case initialization
    call PDE%construct(K, P+1, NODETYPE, PERIODIC, ADVECTION, &
                       GRADIENTS, VISCOSITY, SENSEDADVECTION, KWENO)

    ! Time integrator
    call TimeIntegrator%construct(MAXRES, TSPAN(1), TSPAN(2), TSTEP, &
                                  INT_METHOD, SAVEINT, DYNAMICSTEP,  &
                                  FIRSTDYNSTEP, LASTDYNSTEP,         &
                                  SENSORSTEP, SENSORWINDOW, CFLSTEP,  &
                                  STOPNAN)

    ! Truncation error estimator (for maps and sensors)
    call TruncError%construct(TRUNCERRORTYPE, TECORRECTION, LOWERLIMIT)

    ! Shock-capturing sensor
    call Sensor%construct(SENSORTYPE, SECSENSORTYPE, SENSEDVAR, NUMCLUSTERS, &
                          RAMPBOTTOM, RAMPTOP, SECRAMPBOTTOM, SECRAMPTOP)

    ! h/p adaptator
    !call initAdaptationHP(ADAPTTYPE, ADAPTLOWTHRES, ADAPTHIGHTHRES,       &
                          !SECADAPTLOWTHRES, ADAPTMINORDER, ADAPTMAXORDER, &
                          !MINSIZE, MAXSIZE)

    ! Time integration loop
    call cpu_time(startTime)
    call TimeIntegrator%integrate(TEMAP)
    call cpu_time(endTime)

    ! Print summary of the computation
    call Printer%endSummary(TimeIntegrator%t, endTime-startTime)

    ! Finalisation
    call destructArtificialViscosity
    call destructRiemannSolver

end program main
