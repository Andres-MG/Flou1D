!*******************************************************************************
!  MODULE: MeshClass
!
!> @author
!> Andres Mateo
!
!> @brief
!> Definition of the mesh object that stores all the elements of the domain.
!*******************************************************************************

module MeshClass

    use Constants
    use Utilities, only: linear_mapping
    use StdExpansionClass
    use ElementClass
    use StorageClass
    use Setup_routines
    use Physics, only: Phys, NEQS, EulerFlux
    use WENO, only: WENO_stencil_t, compute_WENO_coeffs

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public types
    public :: Mesh_t

!···············································································
!> @class Mesh_t
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 120
!···············································································
    type Mesh_t
        integer  :: nodeType       = 0           !< initial type of nodes
        integer  :: leftBound      = 0           !< leftmost elements index
        integer  :: rightBound     = 0           !< rightmost elements index
        real(wp) :: width          = 0.0_wp      !< maximum width of the mesh
        logical  :: periodicBCs    = .false.     !< periodic boundary conditions
        class(ElemSt_t), allocatable :: elems    !< array of DGSEM elements
        class(FaceSt_t), allocatable :: faces    !< array of faces
    contains
        procedure, private :: Mesh_assignment
        procedure :: construct          => Mesh_constructor
        procedure :: projectToFaces     => Mesh_project_faces
        procedure :: updateArtViscosity => Mesh_update_artificial_viscosity
        procedure :: updateBDs          => Mesh_get_boundary_values
        procedure :: updateWENOcoeffs   => Mesh_update_WENO_coefficients
        procedure :: computeAvgVals     => Mesh_avg_stencil_values
        generic   :: assignment(=)      => Mesh_assignment
    end type Mesh_t

    ! Generic interface for structured meshes in different dimensions
    interface initConnectivities
        module procedure init_connectivities_1D
    end interface initConnectivities

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 122
!
!> @param[in]  n         number of quadrature nodes
!> @param[in]  xk        coordinates of the mesh nodes
!> @param[in]  nt        type of nodes (Gauss or Gauss-Lobatto)
!> @param[in]  WENOsize  size of the WENO stencil
!> @param[in]  periodic  flag for periodic boundary conditions
!···············································································
subroutine Mesh_constructor(this, n, xk, nt, WENOsize, periodic)
    !* Arguments *!
    integer,  intent(in) :: n
    real(wp), intent(in) :: xk(:)
    integer,  intent(in) :: nt
    integer,  intent(in) :: WENOsize
    logical,  intent(in) :: periodic
    ! Derived types
    class(Mesh_t), intent(inout) :: this

    !* Local variables *!
    integer  :: k
    integer  :: i
    type(Elem_t), allocatable :: aux_elem
    type(Face_t), allocatable :: aux_face

    ! Get the size of the mesh
    k = size(xk) - 1

    ! Construct elements (same order in all of them) and faces
    ! Only static arrays so far
    allocate(ElemVector_t :: this%elems)
    allocate(FaceVector_t :: this%faces)

    ! The first item goes with the constructor
    call this%elems%construct(k)
    call this%faces%construct(k+1)

    do i = 1, k

        allocate(aux_elem)
        allocate(aux_face)

        call aux_elem%construct(xk(i), xk(i+1), n, WENOsize, nt, i)
        call aux_face%construct(i)

        call this%elems%push_back(aux_elem)
        call this%faces%push_back(aux_face)

        deallocate(aux_elem)
        deallocate(aux_face)

    end do

    ! Last face
    allocate(aux_face)
    call aux_face%construct(k+1)
    call this%faces%push_back(aux_face)

    ! Set the rest of the attributes
    this%nodeType    = nt
    this%leftBound   = 1
    this%rightBound  = k
    this%width       = xk(k+1) - xk(1)
    this%periodicBCs = periodic

    ! Initialise connectivities
    call initConnectivities(this)

    ! Compute all the coefficients related to the (SS)WENO scheme
    if (WENOsize > 0) then
        do i = 1, k
            call this%updateWENOcoeffs(WENOsize, i)
        end do
    end if

end subroutine Mesh_constructor

subroutine Mesh_avg_stencil_values(this, elemID, Phi)
    !* Arguments *!
    integer,               intent(in)  :: elemID
    real(wp), allocatable, intent(out) :: Phi(:,:)
    ! Derived types
    class(Mesh_t), target, intent(inout) :: this

    !* Local variables *!
    integer               :: i
    integer               :: n
    integer               :: nodeCnt
    type(Elem_t), pointer :: elem_ptr


    !---- Begin associate ----!
    associate(element => this%elems%at(elemID))
    associate(stencil => element%WENOst)

    ! Allocate the output
    allocate(Phi(stencil%nSt, NEQS))

    elem_ptr => this%elems%at(stencil%leftElem)
    n = elem_ptr%std%n

    ! First element first
    call elem_ptr%computeAvgVals((n+1)-stencil%nLL, n, Phi(:stencil%nLL,:))
    nodeCnt = stencil%nLL

    ! Intermediate elements
    elem_ptr => this%elems%at(elem_ptr%elemRight)
    do i = 2, stencil%numElems-1

        n = elem_ptr%std%n
        call elem_ptr%computeAvgVals(1, n, Phi(nodeCnt+1:nodeCnt+n,:))

        nodeCnt = nodeCnt + n
        elem_ptr => this%elems%at(elem_ptr%elemRight)

    end do

    ! Last element
    call elem_ptr%computeAvgVals(1, stencil%nRR, Phi(nodeCnt+1:,:))

    end associate
    end associate
    !----- End associate -----!

    nullify(elem_ptr)

end subroutine Mesh_avg_stencil_values

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Computes and sets the value of the artificial viscosity coefficient for
!> the current values of 'Phi', 'Grad' and 'sensor'. The value of the scaling
!> factor must be stored in the 'aVis' attribute of the elements.
!
!> @param[in]  this  mesh object to modify
!···············································································
subroutine Mesh_update_artificial_viscosity(this)
    !* Arguments *!
    ! Derived types
    class(Mesh_t), intent(inout) :: this

    !* Local variables *!
    integer :: i
    type(Elem_t), pointer :: elem

    ! Update all the elements, scaling with 'h/p'
    call this%elems%reset_last(0)
    do i = 1, this%elems%size()
        elem => this%elems%next()
        elem%aVis = elem%hp * elem%aVis*Phys%Alpha(1)
    end do

    nullify(elem)

end subroutine Mesh_update_artificial_viscosity

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Calculates the values of 'Phi' at both sides of the interfaces between
!> elements if 'gradient' is .false. or not present. Otherwise, projects the
!> gradients of 'Phi' or the entropy variables.
!
!> @param[inout]  this      mesh object to update
!> @param[in]     gradient  .true. if the gradients are to be projected
!···············································································
subroutine Mesh_project_faces(this, gradient)
    !* Arguments *!
    logical, optional, intent(in)    :: gradient
    class(Mesh_t),     intent(inout) :: this

    !* Local variables *!
    integer :: i
    type(Elem_t), pointer :: elem
    type(Face_t), pointer :: faceL
    type(Face_t), pointer :: faceR

    ! Project the gradients if requested
    if (present(gradient)) then

        if (gradient) then

            call this%elems%reset_last(0)
            do i = 1, this%elems%size()

                elem  => this%elems%next()
                faceL => this%faces%at(elem%faceLeft)
                faceR => this%faces%at(elem%faceRight)

                !---- Begin associate ----!
                associate(grad  => elem%Grad,   &
                          gradL => faceL%GradR, &
                          gradR => faceR%GradL)

                    call elem%projectToFaces(grad, gradL, gradR)

                end associate
                !----- End associate -----!

            end do

            ! Do NOT extrapolate Phi
            return

        end if

    end if

    call this%elems%reset_last(0)
    do i = 1, this%elems%size()

        !---- Begin associate ----!
        elem  => this%elems%next()
        faceL => this%faces%at(elem%faceLeft)
        faceR => this%faces%at(elem%faceRight)
        associate(Phi  => elem%Phi, PhiL => faceL%PhiR, PhiR => faceR%PhiL)

            call elem%projectToFaces(Phi, PhiL, PhiR)

        end associate
        !----- End associate -----!

    end do

    nullify(elem)
    nullify(faceL)
    nullify(faceR)

end subroutine Mesh_project_faces

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Returns the values at the boundaries of the domain.
!
!> @param[in]  time  time instant of the simulation
!···············································································
subroutine Mesh_get_boundary_values(this, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(Mesh_t), intent(inout) :: this

    !* Local variables *!
    type(Elem_t), pointer :: elemLeft
    type(Elem_t), pointer :: elemRight
    type(Face_t), pointer :: faceLeft
    type(Face_t), pointer :: faceRight


    elemLeft  => this%elems%at(this%leftBound)
    faceLeft  => this%faces%at(elemLeft%faceLeft)
    elemRight => this%elems%at(this%rightBound)
    faceRight => this%faces%at(elemRight%faceRight)

    ! Periodic BCs do not require the definition of external values
    if (this%periodicBCs) then

        faceLeft%PhiL   = faceRight%PhiL
        faceLeft%GradL  = faceRight%GradL
        faceRight%PhiR  = faceLeft%PhiR
        faceRight%GradR = faceLeft%GradR

    else

        call externalState(elemLeft%affineMap(-1.0_wp), time, &
                           faceLeft%PhiR, faceLeft%GradR,     &
                           faceLeft%PhiL, faceLeft%GradL)

        call externalState(elemRight%affineMap(1.0_wp), time, &
                           faceRight%PhiL, faceRight%GradL,   &
                           faceRight%PhiR, faceRight%GradR)

    end if

    nullify(elemLeft)
    nullify(elemRight)
    nullify(faceLeft)
    nullifY(faceRight)

end subroutine Mesh_get_boundary_values

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Computes the values of the WENO coefficients with global stencils.
!
!> @param[in]  kWENO   size of the WENO stencil
!> @param[in]  elemID  ID of the element to be updated
!···············································································
subroutine Mesh_update_WENO_coefficients(this, kWENO, elemID)
    !* Arguments *!
    integer, intent(in) :: kWENO
    integer, intent(in) :: elemID
    ! Derived types
    class(Mesh_t), target, intent(inout) :: this

    !* Local variables *!
    integer :: i
    integer :: p
    integer :: iLoc
    integer :: n
    integer :: nL
    integer :: nLL
    integer :: nR
    integer :: nRR
    integer :: nPrev
    integer :: ntmp
    integer :: elemCount
    integer :: elemCountLeft
    integer :: rLow
    integer :: rHigh
    integer :: leftBD
    integer :: rightBD

    real(wp), allocatable :: xc(:)
    real(wp), allocatable :: dx(:)
    real(wp), allocatable :: xtmp(:)
    real(wp), allocatable :: xn(:)

    type(Elem_t),         pointer :: elem
    type(Elem_t),         pointer :: elem_ptr
    class(StdExp_t),      pointer :: std
    type(WENO_stencil_t), pointer :: stencil


    elem    => this%elems%at(elemID)
    std     => elem%std
    stencil => elem%WENOst

    n   = std%n
    nL  = 0
    nR  = 0
    nLL = n    ! If there are no elems. to the left  -> nLL = n
    nRR = n    ! If there are no elems. to the right -> nRR = n
    leftBD  = 0
    rightBD = huge(0) / 10

    ! Go right first
    elem_ptr => elem
    elemCount = 1
    do

        ! Check if the element is at the boundary before moving forward
        if (this%periodicBCs .and. elem_ptr%ID == this%rightBound) then
            rightBD = elemCount
        end if

        ! Check if there are more elements to the right
        if (elem_ptr%elemRight < 1) exit

        elem_ptr => this%elems%at(elem_ptr%elemRight)
        elemCount = elemCount + 1

        ! Count the number of nodes added in this pass
        if (nR + elem_ptr%std%n < kWENO) then
            nRR = elem_ptr%std%n
            nR  = nR + nRR
        else
            nRR = kWENO - nR
            nR  = kWENO
            exit
        end if

    end do

    ! +1 to account for the right-most FV boundary
    nRR = nRR + 1
    nR  = nR  + 1

    ! And then left
    elem_ptr => elem
    elemCountLeft = 0
    do

        ! Check if the element is at the boundary before moving forward
        if (this%periodicBCs .and. elem_ptr%ID == this%leftBound) then
            leftBD = elemCount
        end if

        ! Check if there are more elements to the left
        if (elem_ptr%elemLeft < 1) exit

        elem_ptr => this%elems%at(elem_ptr%elemLeft)
        elemCount = elemCount + 1
        elemCountLeft = elemCountLeft + 1

        ! Count the number of nodes added in this pass
        if (nL + elem_ptr%std%n < kWENO) then
            nLL = elem_ptr%std%n
            nL  = nL + nLL
        else
            nLL = kWENO - nL
            nL  = kWENO
            exit
        end if

    end do

    ! If there are no elems. to the right -> nR = 1
    if (nR == 0) nR = 1

    ! Set the correct values of the 'left/rightBD' variables
    if (this%periodicBCs) then
        rightBD = elemCountLeft + rightBD
        leftBD  = elemCount - leftBD + 1
    end if

    ! Allocate the vector of complementary nodes
    n = nL + n + nR
    allocate(xc(n))

    ! Update the element structure (for main nodes, not complementary ones)
    elem%WENOst%leftElem = elem_ptr%ID
    elem%WENOst%numElems = elemCount
    elem%WENOst%nL       = nL
    elem%WENOst%nLL      = nLL
    elem%WENOst%nRR      = nRR-1
    elem%WENOst%nSt      = n-1

    ! Remember, the loop starts from the left
    xtmp = elem_ptr%getCoords(complementary=.true.)
    if (1 < leftBD) xtmp = xtmp - this%width
    ntmp = size(xtmp)
    xc(:nLL) = xtmp(ntmp-nLL:ntmp-1)
    nPrev = nLL

    do i = 2, elemCount-1

        elem_ptr => this%elems%at(elem_ptr%elemRight)
        xtmp = elem_ptr%getCoords(complementary=.true.)
        if (i < leftBD) then
            xtmp = xtmp - this%width
        else if (i > rightBD) then
            xtmp = xtmp + this%width
        end if
        ntmp = elem_ptr%std%n
        xc(nPrev+1:nPrev+ntmp) = xtmp(:ntmp)
        nPrev = nPrev + ntmp

    end do

    ! The last element goes out of the loop
    elem_ptr => this%elems%at(elem_ptr%elemRight)
    xtmp = elem_ptr%getCoords(complementary=.true.)
    if (elemCount > rightBD) xtmp = xtmp + this%width
    xc(nPrev+1:) = xtmp(:nRR)

    ! Once the vector of nodes is filled, store the cell sizes...
    ! Be careful with the physical BDs of the mesh
    dx = xc(2:) - xc(:n-1)
    if (nL < kWENO) then
        if (nR < kWENO+1) then
            allocate(stencil%dx(0:n))
            stencil%dx(n) = dx(n-1)
        else
            allocate(stencil%dx(0:n-1))
        end if

        stencil%dx(0) = dx(1)

    else if (nR < kWENO+1) then
        allocate(stencil%dx(1:n))
        stencil%dx(n) = dx(n-1)

    else
        allocate(stencil%dx(1:n-1))

    end if
    stencil%dx(1:n-1) = dx

    !... and call the WENO routine
    allocate(xn(std%n))
    do iLoc = 0, std%n

        i = iLoc + nL

        ! Limits for r (left side)
        rLow  = max(0, kWENO-i)
        rHigh = min(kWENO-1, n-1-i)

        ! Skip left boundary of the mesh
        if (rLow <= rHigh) then

            ! Nodes for cell integration
            call linear_mapping(xc(i), xc(i+1), std%x, xn)

            call compute_WENO_coeffs(x=xn, xc=xc,                               &
                                     i=i, k=kWENO,                              &
                                     r_lower=rLow, r_upper=rHigh,               &
                                     cjr=stencil%cjsir(:,eLeft,iLoc,:),         &
                                     cjrpl=stencil%cjsirpl(:,eLeft,iLoc,:,:,:), &
                                     d=stencil%dsir(eLeft,iLoc,:))

            ! Scaling of the derivative coefficient for the sensor integral
            do p = 1, std%n
                stencil%cjsirpl(:,eLeft,iLoc,:,p,:) = stencil%cjsirpl(:,eLeft,iLoc,:,p,:) &
                                                    * sqrt(std%w(p) * stencil%dx(i) / 2.0_wp)
            end do

         end if

        ! Limits for r (right side)
        rLow  = max(1, kWENO-i)
        rHigh = min(kWENO, n-1-i)

        ! Skip right boundary
        if (rLow <= rHigh) then

            ! Nodes for cell integration
            call linear_mapping(xc(i+1), xc(i+2), std%x, xn)

            call compute_WENO_coeffs(x=xn, xc=xc,                                &
                                     i=i, k=kWENO,                               &
                                     r_lower=rLow, r_upper=rHigh,                &
                                     cjr=stencil%cjsir(:,eRight,iLoc,:),         &
                                     cjrpl=stencil%cjsirpl(:,eRight,iLoc,:,:,:), &
                                     d=stencil%dsir(eRight,iLoc,:))

            ! Scaling of the derivative coefficient for the sensor integral
            do p = 1, std%n
                stencil%cjsirpl(:,eRight,iLoc,:,p,:) = stencil%cjsirpl(:,eRight,iLoc,:,p,:) &
                                                     * sqrt(std%w(p) * stencil%dx(i+1) / 2.0_wp)
            end do

        end if

    end do

    nullify(elem)
    nullify(elem_ptr)
    nullify(std)
    nullify(stencil)

end subroutine Mesh_update_WENO_coefficients

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Set the connectivities for an element and its surrounding faces.
!
!> @param[inout]  mesh  mesh where the connectivities are set
!···············································································
subroutine init_connectivities_1D(mesh)
    !* Arguments *!
    type(Mesh_t), intent(inout) :: mesh

    !* Local variables *!
    integer :: i
    type(Elem_t), pointer :: elem
    type(Face_t), pointer :: face

    ! Loop over all the elements
    call mesh%faces%reset_last(0)
    call mesh%elems%reset_last(0)
    do i = 1, mesh%elems%size()

        face => mesh%faces%next()
        face%elemLeft  = i-1
        face%elemRight = i

        elem => mesh%elems%next()
        elem%elemLeft  = i-1
        elem%elemRight = i+1
        elem%faceLeft  = i
        elem%faceRight = i+1

    end do

    ! Left boundary
    face => mesh%faces%front()
    elem => mesh%elems%front()
    if (mesh%periodicBCs) then
        face%elemLeft = mesh%rightBound
        elem%elemLeft = mesh%rightBound
    else
        face%elemLeft = -1
        elem%elemLeft = -1
    end if

    ! Right boundary
    face => mesh%faces%back()
    elem => mesh%elems%back()
    face%elemLeft = mesh%rightBound
    if (mesh%periodicBCs) then
        face%elemRight = mesh%leftBound
        elem%elemRight = mesh%leftBound
    else
        face%elemRight = -1
        elem%elemRight = -1
    end if

    nullify(elem)
    nullify(face)

end subroutine init_connectivities_1D

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Assignment operator overloading. Copies the values in 'rhs' into 'lhs'.
!
!> @param[out]  lhs  Mesh_t object to which the values are copied
!> @param[in]   rhs  Mesh_t object from which the values are copied
!···············································································
subroutine Mesh_assignment(lhs, rhs)
    !* Arguments *!
    class(Mesh_t), target, intent(out) :: lhs
    class(Mesh_t),         intent(in)  :: rhs

    ! Assign numerical attributes
    lhs%nodeType       = rhs%nodeType
    lhs%leftBound      = rhs%leftBound
    lhs%rightBound     = rhs%rightBound
    lhs%width          = rhs%width
    lhs%periodicBCs    = lhs%periodicBCs
    lhs%elems          = rhs%elems
    lhs%faces          = rhs%faces

end subroutine Mesh_assignment

end module MeshClass
