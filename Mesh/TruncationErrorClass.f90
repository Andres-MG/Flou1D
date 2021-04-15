!*******************************************************************************
!  MODULE: TruncationErrorClass
!
!> @author
!> Andres Mateo
!
!> @brief
!> Implements a truncation error estimator for unsteady and steady cases with
!> converged or unconverged solutions.
!*******************************************************************************

module TruncationErrorClass

    use Constants
    use ExceptionsAndMessages
    use Utilities
    use PDEclass
    use ElementClass
    use StdExpansions
    use StdExpansionClass
    use Physics

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public types
    public :: TruncError_t

    ! Explicitly define public variables
    public :: TruncError

!···············································································
!> @class TruncError_t
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Wraps the different implementations of the truncation error estimators.
!···············································································
    type TruncError_t
        logical :: Active
        logical :: Isolated
        logical :: Corrected
        integer :: MapLimit
    contains
        procedure :: construct => TE_constructor
        procedure :: mapError  => TE_map_error
        procedure :: estimate  => TE_local_trunc_error
    end type TruncError_t

    type(TruncError_t) :: TruncError

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Initialises some global variables related to the truncation error estimation
!> module.
!
!> @param[in]  TEtype      truncation error estimation method
!> @param[in]  correction  .true. if the unsteady correction should be applied
!> @param[in]  lowerLimit  minimum order of the estimated TE map
!···············································································
subroutine TE_constructor(this, TEtype, correction, lowerLimit)
    !* Arguments *!
    integer, intent(in) :: TEtype
    logical, intent(in) :: correction
    integer, intent(in) :: lowerLimit
    ! Derived types
    class(TruncError_t), intent(inout) :: this

    ! Set global variables
    this%Corrected = correction
    this%MapLimit  = merge(lowerLimit, 0, lowerLimit >= 0)

    ! Select the estimator chosen by the user
    select case (TEtype)

    case (fNone)
        this%Isolated = .false.
        this%Active   = .false.

    case (eLocalTrunc)
        this%Isolated = .false.
        this%Active   = .true.

    case (eIsolatedLocalTrunc)
        this%Isolated = .true.
        this%Active   = .true.

    case default
        call printError("TruncationError.f90", &
            "The selected truncation error estimator is not implemented.")

    end select

end subroutine TE_constructor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Creates a map of the dependency of the truncation error as a function of the
!> approximation order in each element of the 'mesh'. The mesh fluxes must be
!> updated before calling this function, and the error is computed by using the
!> variable indicated in 'variable' or the first variable if none is provided.
!
!> @param[in]  time  time instant
!···············································································
subroutine TE_map_error(this, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(TruncError_t) :: this

    !* Local variables *!
    integer  :: i
    integer  :: j
    integer  :: k
    real(wp) :: error
    type(Elem_t), pointer :: elem

    ! Do nothing if the truncation error estimator is set to 'none'
    if (.not. this%Active) then
        return
    end if

    ! Loop over the elements...
    call PDE%mesh%elems%reset_last(0)
    do i = 1, PDE%mesh%elems%size()

        elem => PDE%mesh%elems%next()

        ! ...and over the lower order meshes to get the error map
        do j = 1, elem%std%n-1 - this%MapLimit

            ! Proceed only if the auxiliary P>=0
            if (elem%std%n-j <= 0) then
                exit
            end if

            ! Compute the error
            do k = 1, NEQS

                error = this%estimate(elem, j, time, k)

                ! Sort the errors in the direction of increasing order
                elem%Terr(elem%std%n-j, k) = error

            end do

        end do

    end do

    nullify(elem)

end subroutine TE_map_error

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Estimates the truncation error for a converged solution. The values
!> \f$ \tau_P^N \f$ are returned in 'error', where \$ P-N = \text{orderDecrease}
!> \f$.
!
!> @param[in]  elem            element to use
!> @param[in]  orderDecrease   amount to decrease the order of the elements
!> @param[in]  time            time instant of the simulation
!> @param[in]  sensedVariable  variable used to compute the error
!
!> @return                      estimated truncation error in element 'eInd'
!···············································································
function TE_local_trunc_error(this, elem, orderDecrease, time, &
                              sensedVariable) result(error)
    !* Arguments *!
    integer,  intent(in) :: orderDecrease
    real(wp), intent(in) :: time
    integer,  intent(in) :: sensedVariable
    ! Derived types
    class(TruncError_t), intent(in)    :: this
    type(Elem_t),        intent(inout) :: elem

    !* Return values *!
    real(wp) :: error

    !* Local variables *!
    integer :: i
    integer :: nt
    logical :: isolatedElem

    real(wp), allocatable :: sensedVariableTime(:)
    real(wp), allocatable :: origResidual(:)
    real(wp), allocatable :: residual(:)

    class(Stdexp_t), allocatable :: auxStd
    type(Elem_t)                 :: auxElem
    type(Face_t)                 :: tmpFaceLeft
    type(Face_t)                 :: tmpFaceRight


    ! Backup the affected element and faces
    auxElem      = elem
    tmpFaceLeft  = PDE%mesh%faces%val_at(auxElem%faceLeft)
    tmpFaceRight = PDE%mesh%faces%val_at(auxElem%faceRight)

    ! Check if the requested computaton is 'doable'
    if (auxElem%std%n - orderDecrease <= 0) then
        error = huge(1.0)
        return
    end if

    ! Get a reference to a new std. element with the desired order
    nt = merge(eGauss, auxElem%std%nt, auxElem%std%n-orderDecrease == 1)

    ! Interpolate to the new approximation order
    select case (nt)
    case (eGauss)
        allocate(NodalDG_GL :: auxStd)
    case (eGaussLobatto)
        allocate(NodalDG_GLL :: auxStd)
    end select

    call auxStd%construct(auxElem%std%n-orderDecrease)
    call elem%interpTo(PDE%mesh%faces%at(elem%faceLeft),  &
                       PDE%mesh%faces%at(elem%faceRight), &
                       auxStd, project=.true.)

    ! Determine the type of interpolation for the faces
    if (this%Isolated .and. auxStd%n /= 1) then
        isolatedElem = .true.

    ! Non isolated error also when P=0
    else
        isolatedElem = .false.

    end if

    ! And compute the time derivative in the temporary element
    call PDE%elemTimeDerivative(time, elem%ID, isolated=isolatedElem)

    ! Finally, compute the error
    allocate(sensedVariableTime(elem%std%n))
    do i = 1, size(sensedVariableTime)
        sensedVariableTime(i) = getSensedVarDerivative(elem%Phi(i,:),  &
                                                       elem%PhiD(i,:), &
                                                       sensedVariable)
    end do

    residual = sensedVariableTime

    if (this%Corrected) then

        call reallocate(auxElem%std%n, sensedVariableTime)
        do i = 1, size(sensedVariableTime)
            sensedVariableTime(i) = getSensedVarDerivative(auxElem%Phi(i,:),  &
                                                           auxElem%PhiD(i,:), &
                                                           sensedVariable)
        end do

        ! Project to the LOWER order mesh
        allocate(origResidual, mold=residual)
        call auxElem%interpAt(elem%std%x, sensedVariableTime, origResidual)

        residual = residual - origResidual

    end if

    ! And return the maximum estimation of the error
    residual = elem%jac * elem%std%w * residual
    error = maxval(abs(residual))

    ! Reset everything before leaving
    call PDE%mesh%faces%update(tmpFaceLeft,  elem%faceLeft)
    call PDE%mesh%faces%update(tmpFaceRight, elem%faceRight)
    elem = auxElem

end function TE_local_trunc_error

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Computes the value of the time derivative of the sensed variable.
!
!> @param[in]  Phi   function values
!> @param[in]  PhiD  time derivative of the function values
!> @param[in]  var   index of the variable to be sensed (<0 to use \f$p\rho\f$)
!
!> @return           time derivative of the sensed variable
!···············································································
function getSensedVarDerivative(Phi, PhiD, variable) result(varD)
    !* Arguments *!
    real(wp), intent(in) :: Phi(NEQS)
    real(wp), intent(in) :: PhiD(NEQS)
    integer,  intent(in) :: variable

    !* Return values *!
    real(wp) :: varD

    !* Local variables *!
    real(wp) :: u
    real(wp) :: p
    real(wp) :: ut
    real(wp) :: pt

    if (variable > 0) then
        varD = PhiD(variable)

    else
        u = Phi(IRHOU) / Phi(IRHO)
        p = getPressure(Phi)
        ut = PhiD(IRHOU)/Phi(IRHO) - u/Phi(IRHO) * PhiD(IRHO)
        pt = Phys%GammaMinusOne * (PhiD(IRHOE) - 0.5_wp &
           * (PhiD(IRHOU)*u + Phi(IRHOU)*ut))

        varD = pt * Phi(IRHO) + p * PhiD(IRHO)

    end if

end function getSensedVarDerivative

end module TruncationErrorClass
