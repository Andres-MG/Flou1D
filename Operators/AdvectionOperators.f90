module AdvectionOperators

    use Constants, only: wp, eLeft, eRight, eGaussLobatto
    use Physics
    use StdExpansionClass
    use StdExpansions
    use ElementClass
    use MeshClass
    use WENO, only: WENO_interpolation
    use Clustering, only: kMeans

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: AdvectionDGSEMweak
    public :: AdvectionDGSEMstrong
    public :: AdvectionDGSEMsplit
    public :: AdvectionWENO
    public :: AdvectionSSWENO
    public :: AdvectionSSFV
    public :: AdvectionFV
    public :: DGSEMadvection_Int

    abstract interface
        subroutine DGSEMadvection_Int(mesh, elemID)
            import wp, Mesh_t
            type(Mesh_t), intent(inout) :: mesh
            integer,      intent(in)    :: elemID
        end subroutine DGSEMadvection_Int
    end interface

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Advection operator (weak) for the DGSEM.
!
!> @param[inout]  mesh    mesh with all the elements of the domain
!> @param[in]     elemID  ID of the element where the operator acts
!···············································································
subroutine AdvectionDGSEMweak(mesh, elemID)
    !* Arguments *!
    integer, intent(in) :: elemID
    ! Derived types
    type(Mesh_t), intent(inout) :: mesh

    !* Local variabes *!
    integer :: i
    integer :: n
    real(wp), allocatable :: DivF(:,:)
    real(wp), allocatable :: F(:,:)
    type(Elem_t), pointer :: element


    !---- Begin associate ----!
    element => mesh%elems%at(elemID)
    associate(std       => element%std,                     &
              faceLeft  => mesh%faces%at(element%faceLeft), &
              faceRight => mesh%faces%at(element%faceRight))

    n = std%n

    ! Calculate fluxes at the nodes
    allocate(F, mold=element%Phi)
    do i = 1, n
        F(i,:) = EulerFlux(element%Phi(i,:))
    end do

    ! Interior values
    DivF = matmul(std%Dh, F)

    ! Add boundary flux
    if (std%hasBounds) then

        DivF(1,:) = DivF(1,:) - faceLeft%EF  * std%iw(1)
        DivF(n,:) = DivF(n,:) + faceRight%EF * std%iw(n)

    else

        do i = 1, NEQS
            DivF(:,i) = DivF(:,i) + ( faceRight%EF(i) * std%bdMode(eRight,:) &
                      - faceLeft%EF(i) * std%bdMode(eLeft,:) ) * std%iw
        end do

    end if

    ! Add the advection contribution to the time derivative
    element%PhiD = element%PhiD - DivF

    end associate
    !----- End associate -----!

    nullify(element)

end subroutine AdvectionDGSEMweak

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Advection operator (strong) for the DGSEM.
!
!> @param[inout]  mesh    mesh with all the elements of the domain
!> @param[in]     elemID  ID of the element where the operator acts
!···············································································
subroutine AdvectionDGSEMstrong(mesh, elemID)
    !* Arguments *!
    integer, intent(in) :: elemID
    ! Derived types
    type(Mesh_t), intent(inout) :: mesh

    !* Local variabes *!
    integer :: i
    integer :: n
    real(wp), allocatable :: DivF(:,:)
    real(wp), allocatable :: F(:,:)
    type(Elem_t), pointer :: element


    !---- Begin associate ----!
    element => mesh%elems%at(elemID)
    associate(std       => element%std,                     &
              faceLeft  => mesh%faces%at(element%faceLeft), &
              faceRight => mesh%faces%at(element%faceRight))

    n = std%n

    ! Calculate fluxes at the nodes
    allocate(F, mold=element%Phi)
    do i = 1, n
        F(i,:) = EulerFlux(element%Phi(i,:))
    end do

    ! Interior values
    DivF = matmul(std%D, F)

    ! Add boundary flux
    if (std%hasBounds) then

        DivF(1,:) = DivF(1,:) - (faceLeft%EF-F(1,:))  * std%iw(1)
        DivF(n,:) = DivF(n,:) + (faceRight%EF-F(n,:)) * std%iw(n)

    else

        do i = 1, NEQS
            DivF(:,i) = DivF(:,i) + ( (faceRight%EF(i)-faceRight%EFL(i))        &
                      * std%bdMode(eRight,:) - (faceLeft%EF(i)-faceLeft%EFR(i)) &
                      * std%bdMode(eLeft,:) ) * std%iw
        end do

    end if


    ! Add the advection contribution to the time derivative
    element%PhiD = element%PhiD - DivF

    end associate
    !----- End associate -----!

    nullify(element)

end subroutine AdvectionDGSEMstrong

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Advection operator (split form) for the DGSEM. Only with GLL nodes!!
!
!> @param[inout]  mesh    mesh with all the elements of the domain
!> @param[in]     elemID  ID of the element where the operator acts
!···············································································
subroutine AdvectionDGSEMsplit(mesh, elemID)
    !* Arguments *!
    integer, intent(in) :: elemID
    ! Derived types
    type(Mesh_t), intent(inout) :: mesh

    !* Local variables *!
    integer :: i
    integer :: j
    integer :: n
    real(wp), allocatable :: DivF(:,:)
    real(wp), allocatable :: F(:,:)
    real(wp), allocatable :: Fs(:,:,:)
    type(Elem_t), pointer :: element

    !---- Begin associate ----!
    element => mesh%elems%at(elemID)
    associate(faceLeft  => mesh%faces%at(element%faceLeft), &
              faceRight => mesh%faces%at(element%faceRight))

    ! Requires Gauss-Lobatto nodes to work
    select type (std => element%std)
    type is (NodalDG_GLL)

    n = std%n

    ! Fill the principal diagonal of the \f$F^{\sharp}\f$ matrix
    allocate(F(n, NEQS))
    allocate(Fs(n, n, NEQS))
    do i = 1, n
        F(i,:)    = EulerFlux(element%Phi(i,:))
        Fs(i,i,:) = F(i,:)
    end do

    ! And then the rest of the elements
    do i = 1, n
        do j = i+1, n
            Fs(i,j,:) = TwoPointFlux(element%Phi(i,:), element%Phi(j,:))
            Fs(j,i,:) = Fs(i,j,:)
        end do
    end do

    ! Skew-symmetric derivative
    allocate(DivF, mold=element%PhiD)
    DivF = 0.0_wp

    do i = 1, n
        do j = 1, n
            DivF(i,:) = DivF(i,:) + std%Ds(i,j) * Fs(i,j,:)
        end do
    end do

    ! Add boundary flux
    DivF(1,:) = DivF(1,:) - faceLeft%EF  * std%iw(1)
    DivF(n,:) = DivF(n,:) + faceRight%EF * std%iw(n)

    ! Add the advection contribution to the time derivative
    element%PhiD = element%PhiD - DivF

    end select

    end associate
    !----- End associate -----!

    nullify(element)

end subroutine AdvectionDGSEMsplit

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Advection operator (WENO) for the DGSEM.
!
!> @param[inout]  mesh    mesh with all the elements of the domain
!> @param[in]     elemID  ID of the element where the operator acts
!···············································································
subroutine AdvectionWENO(mesh, elemID)
    !* Arguments *!
    integer, intent(in) :: elemID
    ! Derived types
    type(Mesh_t), intent(inout) :: mesh

    !* Local variables *!
    integer  :: n
    integer  :: i
    real(wp), allocatable :: Phi(:,:)
    real(wp), allocatable :: PhiLeft(:,:)
    real(wp), allocatable :: PhiRight(:,:)
    real(wp), allocatable :: Fbar(:,:)
    real(wp), allocatable :: DivF(:,:)
    type(Elem_t), pointer :: element
    type(Face_t), pointer :: faceL
    type(Face_t), pointer :: faceR


    !---- Begin associate ----!
    element => mesh%elems%at(elemID)
    associate(std       => element%std,                     &
              stencil   => element%WENOst,                  &
              faceLeft  => mesh%faces%at(element%faceLeft), &
              faceRight => mesh%faces%at(element%faceRight))

    n = std%n
    allocate(DivF, mold=element%PhiD)
    allocate(Fbar(0:n, NEQS))
    allocate(PhiLeft(0:n, NEQS))
    allocate(PhiRight(0:n, NEQS))

    ! "Averages" in the FV cells
    allocate(Phi(stencil%nSt, NEQS))
    call mesh%computeAvgVals(elemID, Phi)

    !---- Begin associate ----!
    faceL => mesh%faces%at(element%faceLeft)
    faceR => mesh%faces%at(element%faceRight)
    associate(PhiExtLeft  => faceL%PhiL, &
              PhiExtRight => faceR%PhiR)

    ! Get left-side Phi- values...
    if (.not. mesh%periodicBCs .and. elemID == mesh%leftBound) then
        PhiLeft(0,:) = PhiExtLeft
    end if
    call WENO_interpolation(Phi, stencil%nL, stencil%nL+n, eLeft, stencil, &
                            PhiLeft)

    ! ...and right-side Phi+ values
    if (.not. mesh%periodicBCs .and. elemID == mesh%rightBound) then
        PhiRight(n,:) = PhiExtRight
    end if
    call WENO_interpolation(Phi, stencil%nL, stencil%nL+n, eRight, stencil, &
                            PhiRight)

    end associate
    !----- End associate -----!

    nullify(faceL)
    nullify(faceR)

    ! Compute fluxes from Phi- and Phi+ values
    do i = 0, n
        Fbar(i,:) = RiemannSolver(PhiLeft(i,:), PhiRight(i,:), 1.0_wp)
    end do

    ! Finite volume derivative
    do i = 1, n
        DivF(i,:) = (Fbar(i,:) - Fbar(i-1,:)) / stencil%dx(stencil%nL+i)
    end do

    ! Add the advection contribution to the time derivative
    ! Scale with the Jacobian to avoid the inverse Jacobian of the DGSEM
    element%PhiD = element%PhiD - DivF * element%jac

    end associate
    !----- End associate -----!

    nullify(element)

end subroutine AdvectionWENO

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Advection operator (SSWENO) for the DGSEM.
!
!> @param[inout]  mesh    mesh with all the elements of the domain
!> @param[in]     elemID  ID of the element where the operator acts
!···············································································
subroutine AdvectionSSWENO(mesh, elemID)
    !* Arguments *!
    integer, intent(in) :: elemID
    ! Derived types
    type(Mesh_t), intent(inout) :: mesh

    !* Local variables *!
    integer  :: n
    integer  :: i
    integer  :: j
    integer  :: k
    real(wp), allocatable :: FSbar(:,:)
    real(wp), allocatable :: FWbar(:,:)
    real(wp), allocatable :: Fbar(:,:)
    real(wp), allocatable :: Ftmp(:)
    real(wp), allocatable :: PhiAvg(:,:)
    real(wp), allocatable :: PhiLeft(:,:)
    real(wp), allocatable :: PhiRight(:,:)
    real(wp), allocatable :: wEntropy(:,:)
    real(wp), allocatable :: c
    real(wp), allocatable :: b
    real(wp), allocatable :: d
    real(wp), allocatable :: DivF(:,:)
    type(Elem_t), pointer :: element


    !---- Begin associate ----!
    element   => mesh%elems%at(elemID)
    associate(Phi       => element%Phi,                     &
              stencil   => element%WENOst,                  &
              faceLeft  => mesh%faces%at(element%faceLeft), &
              faceRight => mesh%faces%at(element%faceRight))

    ! ONLY for Gauss-Lobatto nodes
    select type (std => element%std)
    type is (NodalDG_GLL)

    n = std%n
    allocate(FSbar(n-1, NEQS), source=0.0_wp)

    ! Entropy-conservative flux on the complementary grid
    do i = 1, n-1
        do j = i+1, n
            do k = 1, i
                Ftmp = TwoPointFlux(Phi(k,:), Phi(j,:))
                FSbar(i,:) = FSbar(i,:) + std%Q(k,j) * Ftmp
            end do
        end do
    end do
    FSbar = 2.0_wp * FSbar

    ! WENO flux
    allocate(PhiLeft(n-1, NEQS))
    allocate(PhiRight(n-1, NEQS))

    ! "Averages" in the FV cells
    allocate(PhiAvg(stencil%nSt, NEQS))
    call mesh%computeAvgVals(elemID, PhiAvg)

    ! Get left-side Phi- values...
    call WENO_interpolation(PhiAvg, stencil%nL+1, stencil%nL+n-1, eLeft, &
                            stencil, PhiLeft)

    ! ...and right-side Phi+ values
    call WENO_interpolation(PhiAvg, stencil%nL+1, stencil%nL+n-1, eRight, &
                            stencil, PhiRight)

    ! Compute fluxes from Phi- and Phi+ values
    allocate(FWbar(n-1, NEQS))
    do i = 1, n-1
        FWbar(i,:) = RiemannSolver(PhiLeft(i,:), PhiRight(i,:), 1.0_wp)
    end do

    ! Linear combination of both fluxes
    c = 1e-12_wp
    allocate(wEntropy(2, NEQS))
    allocate(Fbar(0:n, NEQS))
    do i = 1, n-1

        call getEntropyVars(Phi(i,:),   wEntropy(eLeft,:))
        call getEntropyVars(Phi(i+1,:), wEntropy(eRight,:))

        b = dot_product(wEntropy(eRight,:) - wEntropy(eLeft,:), &
                        FSbar(i,:) - FWbar(i,:))
        d = sqrt(b**2 + c**2)
        d = (d-b) / d

        Fbar(i,:) = FWbar(i,:) + d * (FSbar(i,:) - FWbar(i,:))

    end do

    Fbar(0,:) = faceLeft%EF
    Fbar(n,:) = faceRight%EF

    ! Derivative using the telescopic approach
    allocate(DivF, mold=element%PhiD)
    do i = 1, n
        DivF(i,:) = std%iw(i) * (Fbar(i,:) - Fbar(i-1,:))
    end do

    ! Add boundary flux
    DivF(1,:) = DivF(1,:) - (faceLeft%EF-faceLeft%EFR)   * std%iw(1)
    DivF(n,:) = DivF(n,:) + (faceRight%EF-faceRight%EFL) * std%iw(n)

    ! Add the advection contribution to the time derivative
    element%PhiD = element%PhiD - DivF

    end select

    end associate
    !----- End associate -----!

    nullify(element)

end subroutine AdvectionSSWENO

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Advection operator (Finite Volumes) for the DGSEM.
!
!> @param[inout]  mesh    mesh with all the elements of the domain
!> @param[in]     elemID  ID of the element where the operator acts
!···············································································
subroutine AdvectionFV(mesh, elemID)
    !* Arguments *!
    integer, intent(in) :: elemID
    ! Derived types
    type(Mesh_t), intent(inout) :: mesh

    !* Local variables *!
    integer  :: i
    integer  :: n
    real(wp), allocatable :: Fbar(:,:)
    real(wp), allocatable :: PhiAvg(:,:)
    real(wp), allocatable :: DivF(:,:)
    type(Elem_t), pointer :: element


    !---- Begin associate ----!
    element => mesh%elems%at(elemID)
    associate(std       => element%std,                     &
              Phi       => element%Phi,                     &
              faceLeft  => mesh%faces%at(element%faceLeft), &
              faceRight => mesh%faces%at(element%faceRight))

    n = std%n

    ! FV flux on the complementary grid
    allocate(Fbar(0:n, NEQS))
    allocate(PhiAvg(n, NEQS))

    call element%computeAvgVals(1, n, PhiAvg)

    do i = 1, n-1
        Fbar(i,:) = RiemannSolver(PhiAvg(i,:), PhiAvg(i+1,:), 1.0_wp)
    end do

    Fbar(0,:) = faceLeft%EFR
    Fbar(n,:) = faceRight%EFL

    ! FV derivative (dx = w in the complementary grid)
    allocate(DivF(n, NEQS))
    do i = 1, n
        DivF(i,:) = std%iw(i) * (Fbar(i,:) - Fbar(i-1,:))
    end do

    ! Add boundary flux
    DivF(1,:) = DivF(1,:) - (faceLeft%EF-faceLeft%EFR)   * std%iw(1)
    DivF(n,:) = DivF(n,:) + (faceRight%EF-faceRight%EFL) * std%iw(n)

    ! Add the advection contribution to the time derivative
    element%PhiD = element%PhiD - DivF

    end associate
    !----- End associate -----!

    nullify(element)

end subroutine AdvectionFV

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Advection operator (Entropy-stable FV) for the DGSEM.
!
!> @param[inout]  mesh    mesh with all the elements of the domain
!> @param[in]     elemID  ID of the element where the operator acts
!···············································································
subroutine AdvectionSSFV(mesh, elemID)
    !* Arguments *!
    integer, intent(in) :: elemID
    ! Derived types
    type(Mesh_t), intent(inout) :: mesh

    !* Local variables *!
    integer  :: i
    integer  :: j
    integer  :: k
    integer  :: n
    real(wp) :: b
    real(wp) :: d
    real(wp) :: Ftmp(NEQS)
    real(wp) :: PhiLeft(NEQS)
    real(wp) :: PhiRight(NEQS)
    real(wp) :: jump(NEQS)
    integer,  allocatable :: side(:)
    real(wp), allocatable :: PhiAvg(:,:)
    real(wp), allocatable :: PhiBar(:,:)
    real(wp), allocatable :: PhiV(:,:)
    real(wp), allocatable :: FSbar(:,:)
    real(wp), allocatable :: FVbar(:,:)
    real(wp), allocatable :: Fbar(:,:)
    real(wp), allocatable :: wEntropy(:,:)
    real(wp), allocatable :: DivF(:,:)
    type(Elem_t), pointer :: element


    !---- Begin associate ----!
    element => mesh%elems%at(elemID)
    associate(Phi       => element%Phi,                     &
              faceLeft  => mesh%faces%at(element%faceLeft), &
              faceRight => mesh%faces%at(element%faceRight))

    ! ONLY for Gauss-Lobatto nodes
    select type (std => element%std)
    type is (NodalDG_GLL)

    n = std%n
    allocate(FSbar(n-1, NEQS), source=0.0_wp)

    ! Entropy-conservative flux on the complementary grid
    do i = 1, n-1
        do j = i+1, n
            do k = 1, i
                Ftmp = TwoPointFlux(Phi(k,:), Phi(j,:))
                FSbar(i,:) = FSbar(i,:) + std%Q(k,j) * Ftmp
            end do
        end do
    end do
    FSbar = 2.0_wp * FSbar

    ! Interpolate the complementary grid values
    ! Clusters
    ! allocate(PhiV, source=Phi)
    ! allocate(PhiAvg(2, NEQS))
    ! allocate(side(n))
    ! call kMeans(2, PhiV, PhiAvg, side)

    ! Remove the offset to get less oscillatory polynomials
    ! jump = PhiAvg(2,:) - PhiAvg(1,:)
    ! do i = 1, n
    !     if (side(i) == 2) then
    !         PhiV(i,:) = PhiV(i,:) - jump
    !     end if
    ! end do

    ! Interpolate into the complementary grid
    ! allocate(PhiBar(0:n, NEQS))
    ! PhiBar = matmul(std%lh, PhiV)

    ! FV flux on the complementary grid
    allocate(FVbar(n-1, NEQS))
    do i = 1, n-1

        ! if (side(i) == side(i+1)) then
            PhiLeft  = Phi(i,:)
            PhiRight = Phi(i+1,:)
        ! else
        !     PhiLeft  = PhiBar(i,:)
        !     PhiRight = PhiBar(i,:)
        !     if (side(i) == 2)   PhiLeft  = PhiLeft  + jump
        !     if (side(i+1) == 2) PhiRight = PhiRight + jump
        ! end if

        FVbar(i,:) = SSFVdissipativeFlux(PhiLeft, PhiRight)

    end do

    ! Linear combination of both fluxes
    allocate(wEntropy(2, NEQS))
    allocate(Fbar(0:n, NEQS))
    do i = 1, n-1

        call getEntropyVars(Phi(i,:),   wEntropy(eLeft,:))
        call getEntropyVars(Phi(i+1,:), wEntropy(eRight,:))

        b = dot_product(wEntropy(eRight,:) - wEntropy(eLeft,:), &
                        FSbar(i,:) - FVbar(i,:))
        d = sqrt(b**2 + Phys%SSFVc**2)
        d = (d-b) / d

        Fbar(i,:) = FVbar(i,:) + d * (FSbar(i,:) - FVbar(i,:))

    end do

    Fbar(0,:) = faceLeft%EFR
    Fbar(n,:) = faceRight%EFL

    ! Derivative using the telescopic approach
    allocate(DivF(n, NEQS))
    do i = 1, n
        DivF(i,:) = std%iw(i) * (Fbar(i,:) - Fbar(i-1,:))
    end do

    ! Add boundary flux
    DivF(1,:) = DivF(1,:) - (faceLeft%EF-faceLeft%EFR)   * std%iw(1)
    DivF(n,:) = DivF(n,:) + (faceRight%EF-faceRight%EFL) * std%iw(n)

    ! Add the advection contribution to the time derivative
    element%PhiD = element%PhiD - DivF

    end select

    end associate
    !----- End associate -----!

    nullify(element)

end subroutine AdvectionSSFV

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Dissipative flux (LxF) to be added to the SC flux in the SSFV formulation.
!
!> @param[in]  PhiL  Values on the left side of the complementary face
!> @param[in]  PhiR  Values on the right side of the complementary face
!···············································································
function SSFVdissipativeFlux(PhiL, PhiR) result(Fv)
    !* Arguments *!
    real(wp), intent(in) :: PhiL(NEQS)
    real(wp), intent(in) :: PhiR(NEQS)

    !* Return values *!
    real(wp) :: Fv(NEQS)

    !* Local variables *!
    real(wp) :: Fl(NEQS)
    real(wp) :: Fr(NEQS)
    real(wp) :: lL, lR, l


    Fl = EulerFlux(PhiL)
    Fr = EulerFlux(PhiR)
    Fv = (Fl + Fr) / 2.0_wp

    lL = computeMaxEigenvalue(PhiL)
    lR = computeMaxEigenvalue(PhiR)
    l  = max(lL, lR)

    Fv = Fv - l/2.0_wp * (PhiR-PhiL)

end function SSFVdissipativeFlux

end module AdvectionOperators
