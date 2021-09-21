!*******************************************************************************
!  MODULE: ElementClass
!
!> @author
!> Andres Mateo
!
!> @brief
!> DGSEM elements class for 1D cases.
!*******************************************************************************

module ElementClass

    use Constants
    use ExceptionsAndMessages
    use Utilities
    use Lagrange
    use Legendre
    use Physics
    use StdExpansionClass
    use StdExpansions
    use WENO, only: WENO_stencil_t

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public types
    public :: Elem_t
    public :: Face_t

!···············································································
!> @class Elem_t
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 120
!···············································································
    type Elem_t
        real(wp) :: dx   = 2.0_wp  !< size of the element
        real(wp) :: xL   = -1.0_wp !< position of the left boundary
        real(wp) :: xR   = 1.0_wp  !< position of the right boundary
        real(wp) :: sens = 0.0_wp  !< shock-capturing sensor
        real(wp) :: aVis = 0.0_wp  !< artificial viscosity
        real(wp) :: hp   = 0.0_wp  !< value of dx / (dG%n-1)
        real(wp) :: jac  = 0.0_wp  !< jacobian (constant)
        real(wp) :: invJ = 0.0_wp  !< jacobian (inverse)
        integer  :: pSpace = 0     !< polynomial space type
        integer  :: ID             !< position in the global list
        integer  :: faceLeft       !< index of the left face in faces
        integer  :: faceRight      !< index of the right face in faces
        integer  :: elemLeft       !< index of the left element in elems
        integer  :: elemRight      !< index of the right element in elems
        logical  :: sensed         !< .true. if the sensor detected this element
        real(wp), allocatable :: Terr(:,:)   !< truncation error
        real(wp), allocatable :: PhiD(:,:)   !< time derivative at the nodes
        real(wp), allocatable :: Src(:,:)    !< source term at the nodes
        real(wp), allocatable :: Phi(:,:)    !< function values at the nodes
        real(wp), allocatable :: Fsvv(:,:)   !< SVV fluxes at the nodes
        real(wp), allocatable :: exact(:,:)  !< exact values at the nodes
        real(wp), allocatable :: Grad(:,:)   !< gradient at the nodes
        real(wp), allocatable :: G(:,:)      !< vectors for low storage RK
        type(WENO_stencil_t)  :: WENOst      !< stencil info for WENO methods
        class(StdExp_t), allocatable :: std  !< standard element discretization
    contains
        procedure, private :: s_map => Element_scalar_affine_map
        procedure, private :: v_map => Element_vector_affine_map
        procedure, private :: Element_assignment
        procedure, private :: interpVector => Element_interpolate_values_vector
        procedure, private :: interpMatrix => Element_interpolate_values_matrix
        generic   :: affineMap      => s_map, v_map
        generic   :: interpAt       => interpVector, interpMatrix
        procedure :: construct      => Element_constructor
        procedure :: projectToFaces => Element_face_projector
        procedure :: getCoords      => Element_get_coordinates
        procedure :: interpTo       => Element_interpolate_to_std
        procedure :: toLegendre     => Element_to_Legendre
        procedure :: fromLegendre   => Element_from_Legendre
        procedure :: computeAvgVals => Element_avg_values
        procedure :: computeCFL     => Element_CFL_number
        generic   :: assignment(=)  => Element_assignment
    end type Elem_t

!···············································································
!> @class Face_t
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 116 (modified: only for faces)
!···············································································
    type Face_t
        integer  :: ID            !< position in the list
        integer  :: elemLeft      !< element to the left
        integer  :: elemRight     !< element to the right
        real(wp) :: PhiL(NEQS)    !< function values at the left side
        real(wp) :: PhiR(NEQS)    !< function values at the right side
        real(wp) :: GradL(NEQS)   !< function gradients at the left side
        real(wp) :: GradR(NEQS)   !< function gradients at the right side
        real(wp) :: EFL(NEQS)     !< Euler fluxes at the left side
        real(wp) :: EFR(NEQS)     !< Euler fluxes at the right side
        real(wp) :: VFL(NEQS)     !< viscous fluxes at the left side
        real(wp) :: VFR(NEQS)     !< viscous fluxes at the right side
        real(wp) :: SVVL(NEQS)    !< SVV fluxes at the left side
        real(wp) :: SVVR(NEQS)    !< SVV fluxes at the right side
        real(wp) :: EF(NEQS)      !< Riemann flux
    contains
        procedure, private :: Face_assignment
        procedure :: construct     => Face_constructor
        generic   :: assignment(=) => Face_assignment
    end type Face_t

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Constructor for the Elem_t class.
!
!> @param[in]  xL        position of the right node
!> @param[in]  xR        position of the left node
!> @param[in]  nNodes    number of discretization nodes
!> @param[in]  WENOsize  size of the WENO stencil
!> @param[in]  nType     type of discretization
!> @param[in]  pos       position in the list of elements
!···············································································
subroutine Element_constructor(this, xL, xR, nNodes, WENOsize, nType, pos)
    !* Arguments *!
    real(wp), intent(in) :: xL
    real(wp), intent(in) :: xR
    integer,  intent(in) :: nNodes
    integer,  intent(in) :: WENOsize
    integer,  intent(in) :: nType
    integer,  intent(in) :: pos
    ! Derived types
    class(Elem_t), intent(inout) :: this


    ! First, construct the local standard space
    select case (nType)
    case (eGauss)
        allocate(NodalDG_GL :: this%std)

    case (eGaussLobatto)
        allocate(NodalDG_GLL :: this%std)

    case default
        call printError("ElementClass.f90", &
                        "Only Gauss and Gauss-Lobatto points are allowed.")
    end select

    ! Pass Psvv if the case includes SVV
    if (Phys%IsSVV) then
        call this%std%init(nNodes, Phys%Psvv)
    else
        call this%std%init(nNodes)
    end if

    call this%WENOst%construct(nNodes, WENOsize, this%std%hasBounds, &
                               this%std%x, this%std%xc)

    ! Set some values
    this%xL        = xL
    this%xR        = xR
    this%dx        = xR - xL
    this%sens      = 0.0_wp
    this%aVis      = 0.0_wp
    this%hp        = this%dx / this%std%n
    this%jac       = 0.5_wp * this%dx
    this%invJ      = 1.0_wp / this%jac
    this%pSpace    = eNodal
    this%ID        = pos
    this%faceLeft  = -1
    this%faceRight = -1
    this%elemLeft  = -1
    this%elemRight = -1
    this%sensed    = .false.

    ! Allocate the rest of the attributes
    allocate(this%Terr(this%std%n-1, NEQS), source=0.0_wp)
    allocate(this%PhiD(this%std%n, NEQS), source=0.0_wp)
    allocate(this%Src(this%std%n, NEQS), source=0.0_wp)
    allocate(this%Phi(this%std%n, NEQS), source=0.0_wp)
    allocate(this%exact(this%std%n, NEQS), source=0.0_wp)
    allocate(this%Grad(this%std%n, NEQS), source=0.0_wp)
    allocate(this%G(this%std%n, NEQS), source=0.0_wp)

    if (Phys%IsSVV) allocate(this%Fsvv(this%std%n, NEQS), source=0.0_wp)

    ! Just for printing
    this%Terr = 0.0_wp
    this%G    = 0.0_wp

end subroutine Element_constructor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Algorithm 121
!
!> @param[in]   u   values to project into the faces
!> @param[out]  uL  values at the left face
!> @param[out]  uR  values at the right face
!···············································································
subroutine Element_face_projector(this, u, uL ,uR)
    !* Arguments *!
    real(wp), intent(in)  :: u(:,:)
    real(wp), intent(out) :: uL(:)
    real(wp), intent(out) :: uR(:)
    ! Derived types
    class(Elem_t), intent(in) :: this

    call this%std%project(u, uL, uR)

end subroutine Element_face_projector

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Affine mapping between the reference domain and the physical domain
!
!> @param[in]  xi  coordinate in reference space \f$ [-1,1] \f$
!> @return         coordinate in physical space \f$ [x_l,x_r] \f$
!···············································································
function Element_scalar_affine_map(this, xi)
    !* Arguments *!
    real(wp), intent(in) :: xi
    ! Derived types
    class(Elem_t), intent(in) :: this

    !* Return value *!
    real(wp) :: Element_scalar_affine_map

    Element_scalar_affine_map = this%xL + (xi + 1.0_wp) * this%jac

end function Element_scalar_affine_map

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Affine mapping between the reference domain and the physical domain
!
!> @param[in]  xi  coordinates in reference space \f$ [-1,1] \f$
!> @return         coordinates in physical space \f$ [x_l,x_r] \f$
!···············································································
function Element_vector_affine_map(this, xi)
    !* Arguments *!
    real(wp), intent(in) :: xi(:)
    ! Derived types
    class(Elem_t), intent(in) :: this

    !* Return value *!
    real(wp) :: Element_vector_affine_map(size(xi))

    Element_vector_affine_map = this%xL + (xi + 1.0_wp) * this%jac

end function Element_vector_affine_map

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Returns the position and the values at the nodes (or the complementary grid
!> specified.
!
!> @param[in]   complementary  .true. if 'xc' is to be returned (optional)
!
!> @return                     nodes in physical space
!···············································································
function Element_get_coordinates(this, complementary) result(x)
    !* Arguments *!
    logical, optional, intent(in)  :: complementary
    ! Derived types
    class(Elem_t),       intent(in) :: this

    !* Return values *!
    real(wp), allocatable :: x(:)

    ! Map the Gauss nodes to physical space
    if (present(complementary)) then
        if (complementary) then
            x = this%affineMap(this%std%xc)
        else
            x = this%affineMap(this%std%x)
        end if
    else
        x = this%affineMap(this%std%x)
    end if

end function Element_get_coordinates

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Interpolates the values of the quadrature nodes to 'xp'.
!
!> @param[in]   xp       coordinates in local space where the values are interp.
!> @param[in]   PhiIn    values to interpolate at the nodes of the element
!> @param[out]  PhiOut   interpolated values at the specified nodes
!···············································································
subroutine Element_interpolate_values_matrix(this, xp, PhiIn, PhiOut)
    !* Arguments *!
    real(wp), intent(in)  :: xp(:)
    real(wp), intent(in)  :: PhiIn(:,:)
    real(wp), intent(out) :: PhiOut(:,:)
    ! Derived types
    class(Elem_t), intent(in) :: this

    !* Local variables *!
    integer               :: n
    integer               :: m
    real(wp), allocatable :: T(:,:)

    ! Sizes and allocations
    n = this%std%n
    m = size(xp)
    allocate(T(m,n))

    ! Interpolation
    call polynomial_interpolation_matrix(n, this%std%x, this%std%wb, m, xp, T)

    PhiOut = matmul(T, PhiIn)

end subroutine Element_interpolate_values_matrix

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Interpolates the values of the quadrature nodes to 'xp'.
!
!> @param[in]   xp       coordinates in local space where the values are interp.
!> @param[in]   PhiIn    values to interpolate at the nodes of the element
!> @param[out]  PhiOut   interpolated values at the specified nodes
!···············································································
subroutine Element_interpolate_values_vector(this, xp, PhiIn, PhiOut)
    !* Arguments *!
    real(wp), intent(in)  :: xp(:)
    real(wp), intent(in)  :: PhiIn(:)
    real(wp), intent(out) :: PhiOut(:)
    ! Derived types
    class(Elem_t), intent(in) :: this

    !* Local variables *!
    integer               :: n
    integer               :: m
    real(wp), allocatable :: T(:,:)

    ! Sizes and allocations
    n = this%std%n
    m = size(xp)
    allocate(T(m,n))

    ! Interpolation
    call polynomial_interpolation_matrix(n, this%std%x, this%std%wb, m, xp, T)

    PhiOut = matmul(T, PhiIn)

end subroutine Element_interpolate_values_vector

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Interpolates the values in the quadrature nodes from the quadrature nodes
!> of another element 'elem'. Only Phi is interpolated, so all the remaining
!> allocatable attributes should be considered undefined.
!
!> @param[inout]  faceLeft   face to the left of the element
!> @param[inout]  faceRight  face to the right of the element
!> @param[in]     stdNew     standard element to where the values are interp.
!> @param[in]     project    .true. if the values at the faces must be projected
!···············································································
subroutine Element_interpolate_to_std(this, faceLeft, faceRight, stdNew, project)
    !* Arguments *!
    logical, intent(in) :: project
    ! Derived types
    class(Elem_t),   intent(inout) :: this
    type(Face_t),    intent(inout) :: faceLeft
    type(Face_t),    intent(inout) :: faceRight
    class(StdExp_t), intent(in)    :: stdNew

    !* Local variables *!
    logical               :: matrixChanged
    real(wp), allocatable :: Phi(:,:)

    ! Initialise
    matrixChanged = .false.

    ! Make sure the interpolation matrix is correct
    if (.not. allocated(this%std%inp)) then
        allocate(this%std%inp(stdNew%n, this%std%n))
        matrixChanged = .true.
    else if (size(this%std%inp, 1) /= stdNew%n .or. &
             size(this%std%inp, 2) /= this%std%n) then
        deallocate(this%std%inp)
        allocate(this%std%inp(stdNew%n, this%std%n))
        matrixChanged = .true.
    end if

    ! Compute interpolation matrix
    if (matrixChanged .or. this%std%nt /= stdNew%nt) then
        call polynomial_interpolation_matrix(this%std%n, this%std%x, this%std%wb, &
                                             stdNew%n, stdNew%x, this%std%inp)
    end if

    ! Copy the old values before reallocation
    Phi = this%Phi

    ! Reallocate some arrays
    call reallocate(stdNew%n-1, NEQS, this%Terr)
    call reallocate(stdNew%n, NEQS, this%PhiD)
    call reallocate(stdNew%n, NEQS, this%Src)
    call reallocate(stdNew%n, NEQS, this%Phi)
    call reallocate(stdNew%n, NEQS, this%Grad)
    call reallocate(stdNew%n, NEQS, this%G)

    ! And interpolate values
    this%sens = 0.0_wp
    this%hp   = this%dx/stdNew%n
    this%Terr = 0.0_wp                   !< New order, new errors
    this%Phi  = matmul(this%std%inp, Phi)

    ! Set new standard, reference element
    this%std = stdNew

    ! Finally, extrapolate if requested
    if (project) then
        call this%projectToFaces(this%Phi,  faceLeft%PhiR, faceRight%PhiL)
        call this%projectToFaces(this%Grad, faceLeft%GradR, faceRight%GradL)
    end if

end subroutine Element_interpolate_to_std

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Project the solution from the N-Lagrange basis into the N-Legendre
!> polynomial basis (from nodal to modal form).
!···············································································
subroutine Element_to_Legendre(this)
    !* Arguments *!
    class(Elem_t),       intent(inout) :: this

    ! Only if the current polynomial space is different
    if (this%pSpace == eNodal) then
        this%Phi = matmul(this%std%Fwd, this%Phi)
    end if

    ! Update the polynomial expansion type
    this%pSpace = eModal

end subroutine Element_to_Legendre

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Projects the modal coefficients of the N-Legendre basis back into the
!> N-Lagrange polynomial basis (from modal to nodal form).
!···············································································
subroutine Element_from_Legendre(this)
    !* Arguments *!
    class(Elem_t), intent(inout) :: this

    ! Only if the current polynomial space is different
    if (this%pSpace == eModal) then
        this%Phi = matmul(this%std%Bwd, this%Phi)
    end if

    ! Update the polynomial expansion type
    this%pSpace = eNodal

end subroutine Element_from_Legendre

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Compute the average values of the subcell FV.
!
!> @param[in]   iLeft    first cell to be computed
!> @param[in]   iRight   last cell to be computed
!> @param[out]  Phi      average values of \f$\Phi\f$ at the FV
!···············································································
subroutine Element_avg_values(this, iLeft, iRight, Phi)
    !* Arguments *!
    integer,  intent(in)  :: iLeft
    integer,  intent(in)  :: iRight
    real(wp), intent(out) :: Phi(:,:)
    ! Derived types
    class(Elem_t), intent(in) :: this

    !* Local variables *!
    integer :: i


    do i = iLeft, iRight
        Phi(i,:) = matmul(this%std%avgInp(i,:), this%Phi)
    end do

end subroutine Element_avg_values

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Computes the CFL number of the element.
!
!> @param[in]  dt  time step of the integration method
!> @return         CFL number of the element
!···············································································
function Element_CFL_number(this, dt) result(cfl)
    !* Arguments *!
    real(wp), intent(in) :: dt
    ! Derived types
    class(Elem_t), intent(in) :: this

    !* Return values *!
    real(wp) :: cfl

    !* Local variables *!
    integer :: i

    ! Value at the first node
    cfl = computeMaxEigenvalue(this%Phi(1,:))

    ! Loop over the rest of the nodes
    do i = 2, size(this%Phi, dim=1)
        cfl = max(cfl, computeMaxEigenvalue(this%Phi(i,:)))
    end do

    cfl = cfl * dt / this%dx

end function Element_CFL_number

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Overloaded assignment between 'Elem_t' objects.
!
!> @param[out]  lhs  Elem_t object to which the values are copied
!> @param[in]   rhs  Elem_t object from which the values are copied
!···············································································
subroutine Element_assignment(lhs, rhs)
    !* Arguments *!
    class(Elem_t), intent(out) :: lhs
    class(Elem_t), intent(in)  :: rhs

    lhs%dx        = rhs%dx
    lhs%xL        = rhs%xL
    lhs%xR        = rhs%xR
    lhs%sens      = rhs%sens
    lhs%aVis      = rhs%aVis
    lhs%hp        = rhs%hp
    lhs%jac       = rhs%jac
    lhs%invJ      = rhs%invJ
    lhs%pSpace    = rhs%pSpace
    lhs%ID        = rhs%ID
    lhs%faceLeft  = rhs%faceLeft
    lhs%faceRight = rhs%faceRight
    lhs%elemLeft  = rhs%elemLeft
    lhs%elemRight = rhs%elemRight
    lhs%sensed    = rhs%sensed
    lhs%Terr      = rhs%Terr
    lhs%PhiD      = rhs%PhiD
    lhs%Src       = rhs%Src
    lhs%Phi       = rhs%Phi
    lhs%exact     = rhs%exact
    lhs%Grad      = rhs%Grad
    lhs%G         = rhs%G
    lhs%std       = rhs%std
    lhs%WENOst    = rhs%WENOst

    if (Phys%IsSVV) lhs%Fsvv = rhs%Fsvv

end subroutine Element_assignment

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Constructor of the Face_t class.
!
!> @param[in]  pos  position in the faces list
!···············································································
subroutine Face_constructor(this, pos)
    !* Arguments *!
    integer, intent(in) :: pos
    ! Derived types
    class(Face_t), intent(inout) :: this

    this%ID        = pos
    this%elemLeft  = -1
    this%elemRight = -1

    ! Initialize to 0 just in case
    this%SVVL = 0.0_wp
    this%SVVR = 0.0_wp

end subroutine Face_constructor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Overloaded assignment between 'Face_t' objects.
!
!> @param[out]  lhs  faceClass object to which the values are copied
!> @param[in]   rhs  faceClass  object from which the values are copied
!···············································································
subroutine Face_assignment(lhs, rhs)
    !* Arguments *!
    class(Face_t), intent(out) :: lhs
    class(Face_t), intent(in)  :: rhs

    lhs%ID        = rhs%ID
    lhs%elemLeft  = rhs%elemLeft
    lhs%elemRight = rhs%elemRight
    lhs%PhiL      = rhs%PhiL
    lhs%PhiR      = rhs%PhiR
    lhs%GradL     = rhs%GradL
    lhs%GradR     = rhs%GradR
    lhs%EFL       = rhs%EFL
    lhs%EFR       = rhs%EFR
    lhs%VFL       = rhs%VFL
    lhs%VFR       = rhs%VFR
    lhs%SVVL      = rhs%SVVL
    lhs%SVVR      = rhs%SVVR
    lhs%EF        = rhs%EF

end subroutine Face_assignment

end module ElementClass
