!*******************************************************************************
!  MODULE: hpAdaptation
!
!> @author
!> Andres Mateo
!
!> @brief
!> h- and p-adaptation routines.
!*******************************************************************************

module hpAdaptation

    !use Constants
    !use MeshClass
    !use NodalDGclass
    !use DGSEM1Dclass
    !use Physics, only: NEQS
    !use ExceptionsAndMessages

    !implicit none

    !! Defaults to private
    !private

    !! Explicitly define public functions
    !public :: initAdaptationHP
    !public :: hAdaptation

    !! Local variables
    !integer,  public, protected :: pAdaptationType
    !integer,  public, protected :: pHadaptMinNodes
    !integer,  public, protected :: pHadaptMaxNodes
    !real(wp), public, protected :: pHadaptLowThres
    !real(wp), public, protected :: pHadaptHighThres
    !real(wp), public, protected :: pHadaptSecLowThres
    !real(wp), public, protected :: pHadaptMinSize
    !real(wp), public, protected :: pHadaptMaxSize

!contains

!!···············································································
!!> @author
!!> Andres Mateo
!!
!!  DESCRIPTION
!!> @brief
!!> Initialises some local variables to the user defined values.
!!
!!> @param[in]  adaptType
!!> @param[in]  thresholdLow     min. sensor value before joining elements
!!> @param[in]  thresholdHigh    max. sensor value before splitting an element
!!> @param[in]  secThresholdLow  min. sec. sensor value before joining elements
!!> @param[in]  minOrder         min. order of the approximating polynomials
!!> @param[in]  maxOrder         max. order of the approximating polynomials
!!> @param[in]  minSize          min. width of an element
!!> @param[in]  maxSize          max. width of an element
!!···············································································
!subroutine initAdaptationHP(adaptType, thresholdLow, thresholdHigh, &
                            !secThresholdLow, minOrder, maxOrder,    &
                            !minsize, maxsize)
    !!* Arguments *!
    !integer,  intent(in) :: adaptType
    !real(wp), intent(in) :: thresholdLow
    !real(wp), intent(in) :: thresholdHigh
    !real(wp), intent(in) :: secThresholdLow
    !integer,  intent(in) :: minOrder
    !integer,  intent(in) :: maxOrder
    !real(wp), intent(in) :: minSize
    !real(wp), intent(in) :: maxSize

    !pAdaptationType    = adaptType
    !pHadaptLowThres    = thresholdLow
    !pHadaptHighThres   = thresholdHigh
    !pHadaptSecLowThres = secThresholdLow
    !pHadaptMinNodes    = minOrder+1
    !pHadaptMaxNodes    = maxOrder+1
    !pHadaptMinSize     = minSize
    !pHadaptMaxSize     = maxSize

!end subroutine initAdaptationHP

!!···············································································
!!> @author
!!> Andres Mateo
!!
!!  DESCRIPTION
!!> @brief
!!> Splits the elements where the sensor value is higher than 'HadaptThres'. The
!!> last (n-1) new elements are appended to the list of elements, while the old
!!> element (before splitting) is substituted by the first new element.
!!
!!> @param[inout]  mesh  mesh object to h-adapt
!!···············································································
!subroutine hAdaptation(mesh)
    !!* Arguments *!
    !type(Mesh1D), target, intent(inout) :: mesh

    !!* Local variables *!
    !integer  :: leftElem
    !integer  :: elemInd
    !real(wp) :: newSize

    !! Initialise
    !leftElem = -1

    !! Loop over the elements of the mesh
    !elemInd = mesh%leftBound
    !do while (elemInd >= 0)

        !!---- Begin associate ----!
        !associate(dG   => mesh%dGst(mesh%elems(elemInd)%dG), &
                  !elem => mesh%elems(elemInd))

            !! Only adapt if the sensor is over the thresholds
            !! Join elements if the sensor is too low
            !if (dG%n <= pHadaptMaxNodes) then

                !if ((dG%n >  1 .and. elem%sens < pHadaptLowThres) .or. &
                    !(dG%n <= 1 .and. elem%sens < pHadaptSecLowThres)) then

                    !if (leftElem < 0) then

                        !leftElem = elemInd

                    !else

                        !newSize = mesh%elems(leftElem)%dx + elem%dx

                        !if ( newSize <= pHadaptMaxSize) then
                            !call joinElements(mesh, [leftElem, elemInd])
                        !end if

                        !leftElem = -1

                    !end if

                !! Split elements if the sensor is too high
                !else if (elem%sens > pHadaptHighThres .and. &
                         !dG%n >= pHadaptMinNodes) then

                    !newSize = 0.5_wp * elem%dx

                    !if (newSize >= pHadaptMinSize) then
                        !call checkMeshSize(mesh)    ! Breaks the associate :/
                        !call splitElement(mesh, elemInd)
                    !end if

                    !! Reset the counter if we split the element
                    !leftElem = -1

                !else

                    !! Reset the counter if we do nothing as well
                    !leftElem = -1

                !end if

            !else

                !! Reset the counter in any other case (?)
                !leftElem = -1

            !end if

        !end associate
        !!----- End associate -----!

        !! Next one
        !elemInd = mesh%elems(elemInd)%elemRight

    !end do

!end subroutine hAdaptation

!!···············································································
!!> @author
!!> Andres Mateo
!!
!!  DESCRIPTION
!!> @brief
!!> Subdivides the element in position 'elemInd' into two subelements. The new
!!> element is appended to the list of elements, while the old element (before
!!> splitting) is substituted by the first new element.
!!
!!> @param[inout]  mesh     mesh object to h-adapt
!!> @param[in]     elemInd  index of the element to be subdivided
!!···············································································
!subroutine splitElement(mesh, elemInd)
    !!* Arguments *!
    !integer, intent(in) :: elemInd
    !! Derived types
    !type(Mesh1D), intent(inout) :: mesh

    !!* Local variables *!
    !integer               :: oldN
    !real(wp)              :: oldXl
    !real(wp)              :: oldXr
    !real(wp)              :: newVis
    !real(wp)              :: newDx
    !integer               :: newN
    !integer               :: nodeType
    !integer               :: dGind
    !integer               :: newElemInd
    !integer               :: newFaceInd
    !real(wp), allocatable :: interpNodes(:)
    !real(wp), allocatable :: tmpPhi(:,:)
    !! Derived types
    !type(DGSEM1D) :: tmpElem

    !! Find a place for the new element
    !if (mesh%nElemHoles == 0) then

        !newElemInd    = mesh%lastElem + 1
        !mesh%lastElem = newElemInd

    !else

        !newElemInd = mesh%elemHoles(mesh%nElemHoles)
        !mesh%elemHoles(mesh%nElemHoles) = 0
        !mesh%nElemHoles = mesh%nElemHoles - 1

    !end if

    !! And the new face
    !if (mesh%nFaceHoles == 0) then

        !newFaceInd    = mesh%lastFace + 1
        !mesh%lastFace = newFaceInd

    !else

        !newFaceInd = mesh%faceHoles(mesh%nFaceHoles)
        !mesh%faceHoles(mesh%nFaceHoles) = 0
        !mesh%nFaceHoles = mesh%nFaceHoles - 1

    !end if

    !!---- Begin associate ----!
    !associate(oldElem => mesh%elems(elemInd),             &
              !newElem => mesh%elems(newElemInd),          &
              !oldDG   => mesh%dGst(mesh%elems(elemInd)%dG))

        !! Keep some attributes
        !oldN  = oldDG%n
        !oldXl = oldElem%xl
        !oldXr = oldElem%xr

        !! Initialisation
        !newDx = 0.5_wp * oldElem%dx
        !newN  = ceiling(0.5_wp * oldN)
        !newN  = max(pHadaptMinNodes, newN)

        !! Add standard dG elmt. to the mesh
        !if (newN > 1) then
            !nodeType = oldDG%nt
        !else
            !newN     = 1       ! Never less than one
            !nodeType = eGauss
        !end if
        !call mesh%updateDGlist(newN, nodeType, dGind)

        !! Append the new element
        !! Create a new element
        !call newElem%construct(oldXl, oldXr-newDx, dGind, mesh%dGst, &
                               !mesh%derivativeType, newElemInd)

        !! New value for the artificial viscosity
        !newVis = oldElem%aVis * newElem%hp / oldElem%hp

        !! Interpolate the old element to the new element
        !interpNodes = (mesh%dGst(newElem%dG)%xg - 1.0_wp) * 0.5_wp
        !call oldElem%interpAt(mesh%dGst, interpNodes, oldElem%Phi, newElem%Phi)

        !! Fill the new element
        !newElem%aVis = newVis

        !! And add the face between the two elements
        !call mesh%faces(newFaceInd)%construct(newFaceInd)

        !! Update the old element to be the second new element
        !! Interpolate the values of 'Phi'
        !allocate(tmpPhi, mold=newElem%Phi)
        !interpNodes = (mesh%dGst(newElem%dG)%xg + 1.0_wp) * 0.5_wp
        !call oldElem%interpAt(mesh%dGst, interpNodes, oldElem%Phi, tmpPhi)

        !! Backup the connectivities
        !tmpElem = oldElem

        !! And update the element
        !call oldElem%reset()
        !call oldElem%construct(oldXl+newDx, oldXr, newElem%dG, mesh%dGst, &
                               !mesh%derivativeType, oldElem%ID)
        !oldElem%Phi  = tmpPhi
        !oldElem%aVis = newVis

        !! Reset its previous connectivities
        !oldElem%faceLeft  = tmpElem%faceLeft
        !oldElem%faceRight = tmpElem%faceRight
        !oldElem%elemLeft  = tmpElem%elemLeft
        !oldElem%elemRight = tmpElem%elemRight

    !end associate
    !!----- End associate -----!

    !! Update mesh
    !mesh%nElems = mesh%nElems + 1
    !mesh%nFaces = mesh%nFaces + 1
    !call updateConnectivites(mesh, elemInd, newElemInd, newFaceInd)

!end subroutine splitElement

!!···············································································
!!> @author
!!> Andres Mateo
!!
!!  DESCRIPTION
!!> @brief
!!> Merge the elements indicated by 'elemInds'. The new element gets the ID of
!!> the second one, while the first one is eliminated from the list.
!!
!!> @param[inout]  mesh      mesh object to h-adapt
!!> @param[in]     elemInds  indices of the elements to be merged
!!···············································································
!subroutine joinElements(mesh, elemInds)
    !!* Arguments *!
    !integer, intent(in) :: elemInds(2)
    !! Derived types
    !type(Mesh1D), target, intent(inout) :: mesh

    !!* Local variables *!
    !integer               :: nNew
    !integer               :: nLeft
    !integer               :: nRight
    !integer               :: dGind
    !integer               :: tmpRightElem
    !integer               :: tmpRightFace
    !real(wp)              :: limit
    !real(wp)              :: prevVis(2)
    !integer,  allocatable :: tmpHoles(:)
    !real(wp), allocatable :: xLeft(:)
    !real(wp), allocatable :: xRight(:)
    !real(wp), allocatable :: PhiLeft(:,:)
    !real(wp), allocatable :: PhiRight(:,:)

    !!---- Begin associate ----!
    !associate(elemLeft  => mesh%elems(elemInds(1)), &
              !elemRight => mesh%elems(elemInds(2)))

        !! Scaled artificial viscosities
        !prevVis = [ elemLeft%aVis / elemLeft%hp, elemRight%aVis / elemRight%hp ]

        !! Limit the order of the polynomial
        !nNew = mesh%dGst(elemLeft%dG)%n + mesh%dGst(elemRight%dG)%n
        !if (nNew > 2) then
            !nNew = nNew - 1
        !end if
        !nNew = min(pHadaptMaxNodes, nNew)

        !! Add the new dG element to the list
        !call mesh%updateDGlist(nNew, mesh%nodeType, dGind)

        !! Perform the interpolation
        !limit  = elemLeft%dx / (elemRight%xR - elemLeft%xL) * 2.0_wp - 1.0_wp
        !xLeft  = pack(mesh%dGst(dGind)%xg, mesh%dGst(dGind)%xg <= limit)
        !xRight = pack(mesh%dGst(dGind)%xg, mesh%dGst(dGind)%xg >= limit)
        !xLeft  = 2.0_wp * (xLeft+1.0_wp)  / (1.0_wp+limit) - 1.0_wp
        !xRight = 2.0_wp * (xRight-1.0_wp) / (1.0_wp-limit) + 1.0_wp
        !nLeft  = size(xLeft)
        !nRight = size(xRight)

        !allocate(PhiLeft(nLeft, NEQS))
        !allocate(PhiRight(nRight, NEQS))
        !call elemLeft%interpAt(mesh%dGst, xLeft, elemLeft%Phi, PhiLeft)
        !call elemRight%interpAt(mesh%dGst, xRight, elemRight%Phi, PhiRight)

        !! Backup some data
        !tmpRightElem = elemRight%elemRight
        !tmpRightFace = elemRight%faceRight

        !! Reconstruct the right element -> new element
        !call elemRight%reset()
        !call elemRight%construct(elemLeft%xL, elemRight%xR, dGind, &
                                 !mesh%dGst, mesh%derivativeType,   &
                                 !elemInds(2))

        !! Reset connectivities
        !elemRight%elemLeft  = elemLeft%elemLeft
        !elemRight%faceLeft  = elemLeft%faceLeft
        !elemRight%elemRight = tmpRightElem
        !elemRight%faceRight = tmpRightFace

        !! New artificial viscosity
        !elemRight%aVis = sum(prevVis) * 0.5_wp * elemRight%hp

        !! Average the value that falls between the two elements
        !if (nLeft + nRight > nNew) then

            !elemRight%Phi(:nLeft-1,:) = PhiLeft(:nLeft-1,:)
            !elemRight%Phi(nLeft+1:,:) = PhiRight(2:,:)
            !elemRight%Phi(nLeft,:)    = (PhiLeft(nLeft,:)+PhiRight(1,:))/2.0_wp

        !else

            !elemRight%Phi(:nLeft,:)   = PhiLeft
            !elemRight%Phi(nLeft+1:,:) = PhiRight

        !end if

        !! Reconstruct the left face
        !mesh%faces(elemLeft%faceLeft)%elemRight = elemInds(2)  ! :/

        !! Set some properties in the next element to the left
        !if (elemInds(1) == mesh%leftBound) then
            !mesh%leftBound = elemInds(2)
        !else 
            !mesh%elems(elemLeft%elemLeft)%elemRight = elemInds(2)  ! :/
        !end if

        !! Update the 'holes' vectors of the mesh
        !if (size(mesh%elemHoles) == mesh%nElemHoles) then

            !allocate(tmpHoles(2*mesh%nElemHoles))  ! Magic!!
            !tmpHoles(:size(mesh%elemHoles)) = mesh%elemHoles
            !call move_alloc(tmpHoles, mesh%elemHoles)

        !end if

        !if (size(mesh%faceHoles) == mesh%nfaceHoles) then

            !allocate(tmpHoles(2*mesh%nfaceHoles))  ! Magic!!
            !tmpHoles(:size(mesh%faceHoles)) = mesh%faceHoles
            !call move_alloc(tmpHoles, mesh%faceHoles)

        !end if

        !mesh%elemHoles(mesh%nElemHoles+1) = elemInds(1)
        !mesh%faceHoles(mesh%nFaceHoles+1) = elemLeft%faceRight
        !mesh%nElemHoles = mesh%nElemHoles + 1
        !mesh%nFaceHoles = mesh%nFaceHoles + 1

        !! And update the mesh
        !mesh%nElems = mesh%nElems - 1
        !mesh%nFaces = mesh%nFaces - 1
        !call mesh%faces(elemLeft%faceRight)%reset()
        !call elemLeft%reset()

    !end associate
    !!----- End associate -----!

!end subroutine joinElements

!!···············································································
!!> @author
!!> Andres Mateo
!!
!!  DESCRIPTION
!!> @brief
!!> Sets the connectivities of the subelements and their contiguous elements and
!!> faces.
!!
!!> @param[inout]  mesh        mesh object to h-adapt
!!> @param[in]     elemInd     index of the element to be subdivided (right)
!!> @param[in]     newElemInd  index of the new element (left)
!!> @param[in]     newFaceInd  index of the new face (in between)
!!···············································································
!subroutine updateConnectivites(mesh, elemInd, newElemInd, newFaceInd)
    !!* Arguments *!
    !integer, intent(in) :: elemInd
    !integer, intent(in) :: newElemInd
    !integer, intent(in) :: newFaceInd
    !! Derived types
    !type(Mesh1D), target, intent(inout) :: mesh

    !!* Local variables *!
    !type(faceClass), pointer :: face => null()
    !type(DGSEM1D),   pointer :: elem => null()

    !!---- Begin associate ----!
    !associate(elemOrigin => mesh%elems(elemInd),    &
              !elemNew    => mesh%elems(newElemInd), &
              !faceNew    => mesh%faces(newFaceInd))

        !! Update the right face and its contiguous elements
        !face => mesh%faces(elemOrigin%faceLeft)
        !face%elemRight = newElemInd
        !if (face%elemLeft >= 0) then
            !elem => mesh%elems(face%elemLeft)
            !elem%elemRight = newElemInd
        !else
            !! The leftmost element is now different
            !mesh%leftBound = newElemInd
        !end if

        !! New face
        !faceNew%elemLeft  = newElemInd
        !faceNew%elemRight = elemInd

        !! New element
        !elemNew%elemLeft  = elemOrigin%elemLeft
        !elemNew%elemRight = elemInd
        !elemNew%faceLeft  = elemOrigin%faceLeft
        !elemNew%faceRight = newFaceInd

        !! Update the connectivities of the second element of the subset
        !elemOrigin%elemLeft = newElemInd
        !elemOrigin%faceLeft = newFaceInd

    !end associate
    !!----- End associate -----!

    !! Release memory
    !nullify(face)
    !if (associated(elem)) then
        !nullify(elem)
    !end if

!end subroutine updateConnectivites

!!···············································································
!!> @author
!!> Andres Mateo
!!
!!  DESCRIPTION
!!> @brief
!!> Checks that the size of the elements and faces arrays are large enough to
!!> hold all the data of the mesh. If the arrays are full, reallocates them with
!!> a predefined size margin.
!!
!!> @param[inout]  mesh  mesh object containing the arrays to be allocated
!!··············································································
!subroutine checkMeshSize(mesh)
    !!* Arguments *!
    !type(Mesh1D), intent(inout) :: mesh

    !!* Local variables *!
    !real(wp), parameter :: allocationMargin = 2.0  !< Magic!!
    !integer             :: newSize
    !! Derived types
    !type(DGSEM1D),   allocatable :: tmpElems(:)
    !type(faceClass), allocatable :: tmpFaces(:)

    !! Reallocate the elements storage if it is full
    !if (mesh%nElemHoles == 0 .and. size(mesh%elems) == mesh%nElems) then

        !! Size margin
        !newSize = ceiling(mesh%nElems * allocationMargin)

        !allocate(tmpElems(newSize))
        !tmpElems(:size(mesh%elems)) = mesh%elems
        !call move_alloc(tmpElems, mesh%elems)

    !end if

    !! Same thing for the faces
    !if (mesh%nFaceHoles == 0 .and. size(mesh%faces) == mesh%nFaces) then

        !! Size margin
        !newSize = ceiling(mesh%nFaces * allocationMargin)

        !allocate(tmpFaces(newSize))
        !tmpFaces(:size(mesh%faces)) = mesh%faces
        !call move_alloc(tmpFaces, mesh%faces)

    !end if

!end subroutine checkMeshSize
 
end module hpAdaptation
