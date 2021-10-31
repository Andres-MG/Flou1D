!*******************************************************************************
!  MODULE: PDEclass
!
!> @author
!> Andres Mateo
!
!> @brief
!> Definition of the 'PDE_t' derived type.
!*******************************************************************************

module PDEclass

    use Constants
    use Utilities
    use ExceptionsAndMessages
    use MeshClass
    use ElementClass
    use AdvectionOperators
    use GradientOperators
    use ViscousOperators
    use Physics
    use Setup_routines

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public types
    public :: PDE_t

    ! Explicitly define public variables
    public :: PDE

!···············································································
!> @class PDE_t
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Class containing the information of the problem discretization.
!···············································································
    type PDE_t
        real(wp), allocatable :: L2error(:)   !< max. \f$L_2\f$ error
        real(wp), allocatable :: LinfError(:) !< max. \f$L_\infty\f$ error
        real(wp) :: maxResidual    = 0.0_wp   !< max. residual
        real(wp) :: maxCFL         = 0.0_wp   !< max CFL value
        integer  :: maxCFLelem     = 0        !< element where the CFL is higher
        integer  :: kWENO          = 0        !< WENO stencil size
        logical  :: calculateGrads = .false.  !< gradients flag
        logical  :: viscous        = .false.  !< .false. if inviscid
        logical  :: useEntropyVars = .false.  !< variables for the visc. flux
        type(Mesh_t) :: mesh                  !< description of the mesh
        procedure(DGSEMadvection_Int), &
            pointer, nopass :: advect => null()      !< advection formulation
        procedure(DGSEMadvection_Int), &
            pointer, nopass :: sensAdvect => null()  !< advection with sensor
        procedure(DGSEMgradient_Int), &
            pointer, nopass :: gradient => null()    !< gradient formulation
        procedure(Viscous_Int), &
            pointer, nopass :: viscousTerm => null() !< viscous formulation
    contains
        procedure, private :: PDE_assignment
        procedure :: construct            => PDE_constructor
        procedure :: globalTimeDerivative => PDE_global_time_derivative
        procedure :: elemTimeDerivative   => PDE_element_time_derivative
        procedure :: computeMaxResidual   => PDE_compute_max_residual
        procedure :: computeCFLnumber     => PDE_compute_CFL_number
        procedure :: computeL2error       => PDE_L2_error
        procedure :: computeLinfError     => PDE_Linf_error
        procedure :: printEndSummary      => PDE_print_final_summary
        generic   :: assignment(=)        => PDE_assignment
        final     :: PDE_destruct
    end type PDE_t

    type(PDE_t) :: PDE

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initialises the mesh by calling the 'meshDefinition' in Setup
!
!> @param[in]  k          number of elements in the mesh
!> @param[in]  n          number of quadrature nodes in each element
!> @param[in]  nodeType   type of nodes (Gauss or Gauss-Lobatto)
!> @param[in]  periodic   flag for periodic boundary conditions
!> @param[in]  advection  formulation of the advection operator
!> @param[in]  gradient   formulation of the gradient operator
!> @param[in]  viscosity  formulation of the viscous term
!> @param[in]  sensorAdv  advection operator for the sensed elements
!> @param[in]  kWENO      size of the WENO stencil
!···············································································
subroutine PDE_constructor(this, k, n, nodeType, periodic, advection, &
                           gradient, viscosity, sensorAdv, kWENO)
    !* Arguments *!
    integer, intent(in) :: k
    integer, intent(in) :: n
    integer, intent(in) :: nodeType
    logical, intent(in) :: periodic
    integer, intent(in) :: advection
    integer, intent(in) :: gradient
    integer, intent(in) :: viscosity
    integer, intent(in) :: sensorAdv
    integer, intent(in) :: kWENO
    ! Derived types
    class(PDE_t), intent(inout) :: this

    !* Local variables *!
    integer :: i
    real(wp)              :: xk(k+1)
    real(wp), allocatable :: x(:)
    type(Elem_t), pointer :: elem

    ! Get nodes from Setup file
    call meshDefinition(k, xk)

    ! Construct the mesh
    call this%mesh%construct(n, xk, nodeType, kWENO, periodic)

    ! Initialise the elements
    call this%mesh%elems%reset_last(0)
    do i = 1, this%mesh%elems%size()

        ! Set the initial values
        elem => this%mesh%elems%next()
        x = elem%getCoords()
        call initialCondition(x, elem%Phi)

    end do

    ! Set some attributes
    this%maxResidual = 0.0_wp
    this%maxCFL      = 0.0_wp
    this%kWENO       = kWENO
    this%maxCFLelem  = 0
    allocate(this%L2error(NEQS))
    allocate(this%LinfError(NEQS))

    ! Set the advection type
    select case (advection)

    case (eWeakAdvection)
        this%advect => AdvectionDGSEMweak

    case (eStrongAdvection)
        this%advect => AdvectionDGSEMstrong

    case (eSplitAdvection)
        if (.not. elem%std%hasBounds) then  ! Probe one element
            call printError("PDEclass.f90", &
                "The DGSEM with split-forms requires Gauss-Lobatto nodes.")
        end if
        this%advect => AdvectionDGSEMsplit

    case (eWENOadvection)
        this%advect => AdvectionWENO

    case (eSSWENOadvection)
        if (.not. elem%std%hasBounds) then  ! Probe one element
            call printError("PDEclass.f90", &
                "The DGSEM with SSWENO requires Gauss-Lobatto nodes.")
        end if
        this%advect => AdvectionSSWENO

    case (eFVadvection)
        this%advect => AdvectionFV

    case (eSSFVadvection)
        if (.not. elem%std%hasBounds) then  ! Probe one element
            call printError("PDEclass.f90", &
                "The DGSEM with entropy-stable FV requires Gauss-Lobatto nodes.")
        end if
        this%advect => AdvectionSSFV

    case default
        call printError("PDEclass.f90", &
            "The requested advection operator is not implemented.")

    end select

    ! Set the advection type in sensed elements
    select case (sensorAdv)

    case (eWeakAdvection)
        this%sensAdvect => AdvectionDGSEMweak

    case (eStrongAdvection)
        this%sensAdvect => AdvectionDGSEMstrong

    case (eSplitAdvection)
        if (.not. elem%std%hasBounds) then  ! Probe one element
            call printError("PDEclass.f90", &
                "The DGSEM with split-forms requires Gauss-Lobatto nodes.")
        end if
        this%sensAdvect => AdvectionDGSEMsplit

    case (eWENOadvection)
        this%sensAdvect => AdvectionWENO

    case (eSSWENOadvection)
        if (.not. elem%std%hasBounds) then  ! Probe one element
            call printError("PDEclass.f90", &
                "The DGSEM with SSWENO requires Gauss-Lobatto nodes.")
        end if
        this%sensAdvect => AdvectionSSWENO

    case (eFVadvection)
        this%sensAdvect => AdvectionFV

    case (eSSFVadvection)
        if (.not. elem%std%hasBounds) then  ! Probe one element
            call printError("PDEclass.f90", &
                "The DGSEM with entropy-stable FV requires Gauss-Lobatto nodes.")
        end if
        this%sensAdvect => AdvectionSSFV

    case default
        this%sensAdvect => this%advect

    end select

    ! Set the gradient type
    select case (gradient)

    case (eWeakGradient)
        this%gradient => GradientDGSEMweak

    case (eStrongGradient)
        this%gradient => GradientDGSEMstrong

    case default
        call printError("PDEclass.f90", &
            "Available gradient formulation are: weak and strong")

    end select

    ! Set the viscous term
    select case (viscosity)

    case (eBR1)
        this%viscousTerm => ViscousBR1

    case default
        call printError("PDEclass.f90", &
            "Only the BR1 viscous scheme is implemented so far.")

    end select

    nullify(elem)

end subroutine PDE_constructor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 122
!
!> @param[in]  time  time instant when the derivative is computed
!···············································································
subroutine PDE_global_time_derivative(this, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(PDE_t), intent(inout) :: this

    !* Local variables *!
    integer               :: i
    real(wp), allocatable :: x(:)
    type(Elem_t), pointer :: elem
    type(Face_t), pointer :: face


    ! Reset the value of the derivative
    call this%mesh%elems%reset_last(0)
    do i = 1, this%mesh%elems%size()

        elem => this%mesh%elems%next()
        elem%PhiD = 0.0_wp

    end do

    ! Approximate values at the boundaries of the elements
    call this%mesh%projectToFaces()

    ! Get external boundary values (wrong gradients here!!)
    call this%mesh%updateBDs(time)

    ! Calculate face fluxes
    call this%mesh%faces%reset_last(0)
    do i = 1, this%mesh%faces%size()

        ! Euler flux
        face    => this%mesh%faces%next()
        face%EF  = RiemannSolver(face%PhiL, face%PhiR, 1.0_wp)
        face%EFL = EulerFlux(face%PhiL)
        face%EFR = EulerFlux(face%PhiR)

    end do

    ! Compute numerical viscous fluxes
    if (Phys%IsViscous) then
        call this%viscousTerm(this%mesh, time, this%gradient)
    end if

    do i = 1, this%mesh%elems%size()

        !---- Begin associate ----!
        elem => this%mesh%elems%next()
        associate(faceLeft  => this%mesh%faces%at(elem%faceLeft), &
                  faceRight => this%mesh%faces%at(elem%faceRight))

        ! Calculate the advective term
        if (elem%saturated) then
            call this%sensAdvect(this%mesh, elem%ID)
        else
            call this%advect(this%mesh, elem%ID)
        end if

        ! Scale the advective and viscous terms
        elem%PhiD = elem%invJ * elem%PhiD

        ! And add the source terms
        x = elem%getCoords()
        call SourceTerm(x, time, elem%Src)
        elem%PhiD = elem%PhiD + elem%Src

        end associate
        !----- End associate -----!

    end do

    nullify(elem)
    nullify(face)

end subroutine PDE_global_time_derivative

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Return the time derivative of the solution inside the element 'ind'.
!
!> @param[in]  time      time instant when the derivative is computed
!> @param[in]  ind       index of the element
!> @param[in]  isolated  .true. if the element is 'isolated' from the mesh
!···············································································
subroutine PDE_element_time_derivative(this, time, ind, isolated)
    !* Arguments *!
    real(wp), intent(in) :: time
    integer,  intent(in) :: ind
    logical,  intent(in) :: isolated
    ! Derived types
    class(PDE_t) :: this

    !* Local variables *!
    real(wp), allocatable :: x(:)
    type(Elem_t), pointer :: elem
    type(Face_t), pointer :: faceLeft
    type(Face_t), pointer :: faceRight

    !---- Begin associate ----!
    elem      => this%mesh%elems%at(ind)
    faceLeft  => this%mesh%faces%at(elem%faceLeft)
    faceRight => this%mesh%faces%at(elem%faceRight)

    ! Reset the value of the derivative
    elem%PhiD = 0.0_wp

    ! Approximate values at the boundaries of the element
    call elem%projectToFaces(elem%Phi, faceLeft%PhiR, faceRight%PhiL)

    ! Update BDs if the element is at a physical boundary
    if (ind == this%mesh%leftBound .or. ind == this%mesh%rightBound) then
        call this%mesh%updateBDs(time)
    end if

    ! Calculate face fluxes
    faceLeft%EFR  = EulerFlux(faceLeft%PhiR)
    faceRight%EFL = EulerFlux(faceRight%PhiL)
    if (isolated) then
        faceLeft%EF  = faceLeft%EFR
        faceRight%EF = faceRight%EFL
    else
        faceLeft%EF  = RiemannSolver(faceLeft%PhiL,  faceLeft%PhiR,  1.0_wp)
        faceRight%EF = RiemannSolver(faceRight%PhiL, faceRight%PhiR, 1.0_wp)
    end if

    ! Compute numerical viscous fluxes
    if (Phys%IsViscous) then
        call this%viscousTerm(this%mesh, time, this%gradient, ind, isolated)
    end if

    ! Calculate the advective term
    if (elem%saturated) then
        call this%sensAdvect(this%mesh, elem%ID)
    else
        call this%advect(this%mesh, elem%ID)
    end if

    ! Scale the advective and viscous terms
    elem%PhiD = elem%invJ * elem%PhiD

    ! And add the source term
    x = elem%getCoords()
    call SourceTerm(x, time, elem%Src)
    elem%PhiD = elem%PhiD + elem%Src

    nullify(elem)
    nullify(faceLeft)
    nullify(faceRight)

end subroutine PDE_element_time_derivative

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Norm \$L_{\infty}\$ of the steady-state residual. Calculated as
!> \$ R^N(u^N) = \left|\left| \dot{\Phi}_i \right|\right|_{\infty} \$
!> NOTE: PhiD might need to be updated before the call.
!···············································································
subroutine PDE_compute_max_residual(this)
    !* Arguments *!
    class(PDE_t), target, intent(inout) :: this

    !* Local variables *!
    integer               :: i
    type(Elem_t), pointer :: elem => null()

    !--- Begin associate ---!
    associate(mesh => this%mesh)

    ! Compute the local residual and choose the highest one
    this%maxResidual = 0.0_wp
    call mesh%elems%reset_last(0)
    do i = 1, mesh%elems%size()

        elem => mesh%elems%next()
        this%maxResidual = max(this%maxResidual, &
                               maxval(abs(elem%PhiD)))

    end do

    nullify(elem)

    end associate
    !----- End associate -----!

end subroutine PDE_compute_max_residual

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Computes the maximum CFL number of the mesh.
!
!> @param[in]  dt  time step of the integration
!···············································································
subroutine PDE_compute_CFL_number(this, dt)
    !* Arguments *!
    real(wp), intent(in) :: dt
    ! Derived types
    class(PDE_t), intent(inout) :: this

    !* Local variables *!
    integer               :: i
    real(wp)              :: tmpCFL
    type(Elem_t), pointer :: elem => null()

    !--- Begin associate ---!
    associate(mesh => this%mesh)

    ! Initialise
    this%maxCFL = -1.0_wp

    ! Loop over the elements
    call mesh%elems%reset_last(0)
    do i = 1, mesh%elems%size()

        elem => mesh%elems%next()
        tmpCFL = elem%computeCFL(dt)

        ! Select the maximum CFL number
        if (tmpCFL > this%maxCFL) then
            this%maxCFL     = tmpCFL
            this%maxCFLelem = i
        end if

    end do

    nullify(elem)

    end associate
    !----- End associate -----!

end subroutine PDE_compute_CFL_number

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Update the \f$L^2\f$ error between the computed solution and the exact
!> solution (if no exact solution is provided, the error will be 0).
!
!> @param[in]  time  time instant
!···············································································
subroutine PDE_L2_error(this, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(PDE_t), intent(inout) :: this

    !* Local variables *!
    integer               :: i
    integer               :: nDOFs
    logical               :: available_exact
    real(wp)              :: tmpErr(NEQS)
    real(wp), allocatable :: x(:)
    type(Elem_t), pointer :: elem

    !--- Begin associate ---!
    associate(mesh => this%mesh)

    ! Initialisation
    nDOFs        = 0
    this%L2error = 0.0_wp

    ! Add the error of the elements to the total error
    call mesh%elems%reset_last(0)
    do i = 1, mesh%elems%size()

        elem => mesh%elems%next()

        ! Compute the exact solution
        x = elem%getCoords()
        available_exact = exactSolution(x, time, elem%exact)

        ! Add to the error if the exact solution exists
        if (available_exact) then

            tmpErr = sum( (elem%Phi - elem%exact) ** 2, dim=1)
            this%L2error = this%L2error + tmpErr
            nDOFs = nDOFs + elem%std%n

        end if

    end do

    nullify(elem)

    if (nDOFs /= 0) then
        this%L2error = sqrt(this%L2error/nDOFs)
    end if

    end associate
    !----- End associate -----!

end subroutine PDE_L2_error

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Update the \f$L^{\infty}\f$ error between the computed solution and the
!> exact solution (if no exact solution is provided, the error will be 0).
!
!> @param[in]  time  time instant when the error is calculated
!···············································································
subroutine PDE_Linf_error(this, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(PDE_t), intent(inout) :: this

    !* Local variables *!
    integer               :: i
    logical               :: available_exact
    real(wp)              :: tmpErr(NEQS)
    real(wp), allocatable :: x(:)
    type(Elem_t), pointer :: elem

    !--- Begin associate ---!
    associate(mesh => this%mesh)

    ! Initialisation
    this%LinfError = 0.0_wp

    ! Find the maximum error
    call mesh%elems%reset_last(0)
    do i = 1, mesh%elems%size()

        elem => mesh%elems%next()

        ! Compute the exact solution
        x = elem%getCoords()
        available_exact = exactSolution(x, time, elem%exact)

        ! Add to the error if the exact solution exists
        if (available_exact) then

            tmpErr = maxval(abs(elem%Phi - elem%exact), dim=1)
            this%LinfError = max(this%LinfError, tmpErr)

        end if

    end do

    !---- End associate ----!
    end associate

    nullify(elem)

end subroutine PDE_Linf_error

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Print a summary of the case results once the computation is finished.
!
!> @param[in]  time  time instant when the error is calculated
!···············································································
subroutine PDE_print_final_summary(this, time, elapsed_time)
    !* Arguments *!
    real(wp), intent(in) :: time
    real(wp), intent(in) :: elapsed_time
    ! Derived types
    class(PDE_t), intent(inout) :: this

    !* Local variables *!
    integer                 :: i
    character(len=CHAR_LEN) :: fmat_res
    character(len=CHAR_LEN) :: fmat_err

    ! Get values
    call this%computeMaxResidual()
    call this%computeL2error(time)
    call this%computeLinfError(time)

    ! Set output format
#ifdef SINGLE_PRECISION
        fmat_res = "(a,e12.6)"
        fmat_err = "(a,a,a,e12.6)"
#else
        fmat_res = "(a,e22.16)"
        fmat_err = "(a,a,a,e22.16)"
#endif

    ! Print results
    write(stdout, "(a)") ""
    write(stdout, "(a)") "****************************************"
    write(stdout, "(a,f0.3,a)") " Exec. time: ", elapsed_time, " s"
    write(stdout, "(a)") ""

    write(stdout, fmat_res) " Max. residual: ", this%maxResidual
    write(stdout, "(a)") ""

    write(stdout, "(a)") " L2 error:"
    do i = 1, NEQS
        write(stdout, fmat_err) "    + ", VarNames(i), ": ", this%L2error(i)
    end do
    write(stdout, "(a)") ""

    write(stdout, "(a)") " Linf error:"
    do i = 1, NEQS
        write(stdout, fmat_err) "    + ", VarNames(i), ": ", this%LinfError(i)
    end do

    write(stdout, "(a)") "****************************************"

end subroutine PDE_print_final_summary

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Assignment operator overloading. Copies the values in 'rhs' into 'lhs'.
!
!> @param[out]  lhs  PDE_t object to which the values are copied
!> @param[in]   rhs  PDE_t object from which the values are copied
!···············································································
subroutine PDE_assignment(lhs, rhs)
    !* Arguments *!
    class(PDE_t), target, intent(out) :: lhs
    class(PDE_t),         intent(in)  :: rhs

    lhs%L2error        = rhs%L2error
    lhs%LinfError      = rhs%LinfError
    lhs%maxResidual    = rhs%maxResidual
    lhs%maxCFL         = rhs%maxCFL
    lhs%kWENO          = rhs%kWENO
    lhs%maxCFLelem     = rhs%maxCFLelem
    lhs%calculateGrads = rhs%calculateGrads
    lhs%viscous        = rhs%viscous
    lhs%useEntropyVars = rhs%useEntropyVars
    lhs%mesh           = rhs%mesh
    lhs%advect        => rhs%advect
    lhs%sensAdvect    => rhs%sensAdvect
    lhs%gradient      => rhs%gradient
    lhs%viscousTerm   => rhs%viscousTerm

end subroutine PDE_assignment

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Releases the memory used by the instance of 'PDE_t'.
!···············································································
subroutine PDE_destruct(this)
    !* Arguments *!
    type(PDE_t), intent(inout) :: this

    if (associated(this%advect)) then
        nullify(this%advect)
    end if
    if (associated(this%sensAdvect)) then
        nullify(this%sensAdvect)
    end if
    if (associated(this%gradient)) then
        nullify(this%gradient)
    end if
    if (associated(this%viscousTerm)) then
        nullify(this%viscousTerm)
    end if

end subroutine PDE_destruct

end module PDEclass
