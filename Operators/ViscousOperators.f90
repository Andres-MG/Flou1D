module ViscousOperators

    use Constants, only: wp
    use Utilities
    use Physics
    use MeshClass
    use ElementClass
    use StdExpansionClass
    use GradientOperators

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: ViscousBR1
    public :: Viscous_Int

    abstract interface
        subroutine Viscous_Int(mesh, time, gradient, elemInd, isolated)
            import wp, Mesh_t, DGSEMgradient_Int
            type(Mesh_t),      intent(inout) :: mesh
            real(wp),          intent(in)    :: time
            integer, optional, intent(in)    :: elemInd
            logical, optional, intent(in)    :: isolated
            procedure(DGSEMgradient_Int)     :: gradient
        end subroutine Viscous_Int
    end interface

contains

subroutine ViscousBR1(mesh, time, gradient, elemInd, isolated)
    !* Arguments *!
    real(wp),          intent(in) :: time
    integer, optional, intent(in) :: elemInd
    logical, optional, intent(in) :: isolated
    ! Derived types
    type(Mesh_t), intent(inout) :: mesh
    ! Procedures
    procedure(DGSEMgradient_Int) :: gradient

    !* Local variables *!
    logical  :: isolatedElem
    integer  :: iElemMin
    integer  :: iElemMax
    integer  :: iLeft
    integer  :: iRight
    integer  :: iStep
    integer  :: i
    integer  :: k
    integer  :: n
    real(wp) :: visc
    real(wp) :: viscL(NEQS)
    real(wp) :: viscR(NEQS)
    real(wp) :: PhiL(NEQS)
    real(wp) :: PhiR(NEQS)
    real(wp) :: flux(NEQS)
    real(wp), allocatable :: viscFlux(:,:)
    type(Elem_t), pointer :: elem
    type(Face_t), pointer :: face
    type(Face_t), pointer :: faceLeft
    type(Face_t), pointer :: faceRight


    if (present(isolated)) then
        isolatedElem = isolated
    else
        isolatedElem = .false.
    end if

    ! Global or one element
    if (present(elemInd)) then
        elem => mesh%elems%at(elemInd)
        iElemMin = elemInd
        iElemMax = elemInd
        iLeft    = elem%faceLeft
        iRight   = elem%faceRight
        iStep    = iRight - iLeft
    else
        iElemMin = 1
        iElemMax = mesh%elems%size()
        iLeft    = 1
        iRight   = mesh%faces%size()
        iStep    = 1
    end if

    ! Update the gradients at the interfaces
    do i = iLeft, iRight, iStep

        face => mesh%faces%at(i)

        ! Interface flux
        if (Phys%WithEntropyVars) then
            call getEntropyVars(face%PhiL, PhiL)
            call getEntropyVars(face%PhiR, PhiR)
        else
            PhiL = face%PhiL
            PhiR = face%PhiR
        end if

        ! Riemann solver
        if (isolatedElem) then

            face%GradL = PhiL
            face%GradR = PhiR

        else

            ! BR1 scheme
            flux = 0.5_wp * ( PhiL + PhiR )
            face%GradL = flux
            face%GradR = flux

        end if

    end do

    ! Calculate gradients
    call mesh%elems%reset_last(iElemMin-1)
    do i = iElemMin, iElemMax

        elem      => mesh%elems%next()
        faceLeft  => mesh%faces%at(elem%faceLeft)
        faceRight => mesh%faces%at(elem%faceRight)

        ! Update the volume values
        if (Phys%WithEntropyVars) then

            ! Entropy vars if required
            do k = 1, elem%std%n
                call getEntropyVars(elem%Phi(k,:), elem%Grad(k,:))
            end do

        else

            elem%Grad = elem%Phi

        end if

        ! And calculate the gradient
        elem%Grad = gradient(mesh, elem%Grad, faceLeft%GradR, faceRight%GradL, &
                             faceLeft%PhiR, faceRight%PhiL, elem%ID)
        elem%Grad = elem%invJ * elem%Grad

        call elem%projectToFaces(elem%Grad, faceLeft%GradR, faceRight%GradL)

    end do

    ! Reload the BDs in case they depend on the inner gradients
    if (.not. present(elemInd)) then
        call mesh%updateBDs(time)
    else if (elemInd == mesh%leftBound .or. elemInd == mesh%rightBound) then
        call mesh%updateBDs(time)
    end if

    ! Compute SVV interior and face fluxes
    if (Phys%IsSVV) then

        call mesh%elems%reset_last(iElemMin-1)
        do i = iElemMin, iElemMax

            elem      => mesh%elems%next()
            faceLeft  => mesh%faces%at(elem%faceLeft)
            faceRight => mesh%faces%at(elem%faceRight)

            ! When the sensor is in (0,1)
            if (elem%sensed .and. .not. elem%saturated) then
                elem%Fsvv = SVVflux(elem%Phi, elem%Grad, elem%std%Hsvv)
                call elem%projectToFaces(elem%Fsvv, faceLeft%SVVR, faceRight%SVVL)

            ! Let the more dissipative schemes work otherwise
            else
                elem%Fsvv      = 0.0_wp
                faceLeft%SVVR  = 0.0_wp
                faceRight%SVVL = 0.0_wp

            end if

        end do

    end if

    ! Face fluxes loop
    do i = iLeft, iRight, iStep

        face => mesh%faces%at(i)

        ! NS viscosity
        viscL = ViscousFlux(face%PhiL, face%GradL)
        viscR = ViscousFlux(face%PhiR, face%GradR)

        ! Artificial viscosity
        ! At the physical BDs, use the interior viscosities
        if (face%elemLeft >= 0) then
            elem => mesh%elems%at(face%elemLeft)
            visc = elem%aVis
        else
            elem => mesh%elems%at(face%elemRight)
            visc = elem%aVis
        end if

        if (visc > 0.0_wp) then
            viscL = viscL + ArtViscousFlux(face%PhiL, face%GradL, visc)
        end if

        if (face%elemRight >= 0) then
            elem => mesh%elems%at(face%elemRight)
            visc = elem%aVis
        else
            elem => mesh%elems%at(face%elemLeft)
            visc = elem%aVis
        end if

        if (visc > 0.0_wp) then
            viscR = viscR + ArtViscousFlux(face%PhiR, face%GradR, visc)
        end if

        ! SVV
        if (Phys%IsSVV) then
            viscL = viscL + face%SVVL
            viscR = viscR + face%SVVR
        end if

        ! Riemann solver
        if (isolatedElem) then
            face%VFL = viscL
            face%VFR = viscR
        ! BR1 scheme
        else
            flux = 0.5_wp * (viscL + viscR)
            face%VFL = flux
            face%VFR = flux
        end if

    end do

    ! Compute the integrals, interior and boundary terms
    call mesh%elems%reset_last(iElemMin-1)
    do i = iElemMin, iElemMax

        !---- Begin associate ----!
        elem      => mesh%elems%next()
        faceLeft  => mesh%faces%at(elem%faceLeft)
        faceRight => mesh%faces%at(elem%faceRight)
        associate(std        => elem%std,     &
                  FstarLeft  => faceLeft%VFR, &
                  FstarRight => faceRight%VFL)

        n = std%n

        ! There are three contributions to the viscous fluxes:
        !   - NS
        !   - Artificial viscosity activated by the sensor
        !   - Artificial viscosity deactivated by the sensor (SVV)
        call reallocate(n, NEQS, viscFlux)
        do k = 1, n
            viscFlux(k,:) = ViscousFlux(elem%Phi(k,:), elem%Grad(k,:))
            viscFlux(k,:) = viscFlux(k,:) + ArtViscousFlux(elem%Phi(k,:),  &
                                                           elem%Grad(k,:), &
                                                           elem%aVis)
        end do
        if (Phys%IsSVV .and. .not. elem%sensed) viscFlux = viscFlux + elem%Fsvv

        ! Interior values
        viscFlux = matmul(std%Dh, viscFlux)

        ! Add boundary flux
        if (std%hasBounds) then

            viscFlux(1,:) = viscFlux(1,:) - FstarLeft  * std%iw(1)
            viscFlux(n,:) = viscFlux(n,:) + FstarRight * std%iw(n)

        else

            do k = 1, NEQS
                viscFlux(:,k) = viscFlux(:,k) &
                              + ( FstarRight(k) * std%bdMode(eRight,:) &
                              - FstarLeft(k) * std%bdMode(eLeft,:) ) * std%iw
            end do

        end if


        ! Finally, add the viscous term to the time derivative
        elem%PhiD = elem%PhiD + viscFlux

        end associate
        !----- End associate -----!

    end do

    nullify(elem)
    nullify(face)
    nullify(faceLeft)
    nullify(faceRight)

end subroutine ViscousBR1

end module ViscousOperators
