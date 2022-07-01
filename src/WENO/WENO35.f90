!*******************************************************************************
!  MODULE: WENO35
!
!> @author
!> Andres Mateo
!
!> @brief
!> Implementation of the WENO35 algorithm for equispaced meshes only.
!*******************************************************************************

module WENO35

    use Constants
    use ElementClass
    use ExceptionsAndMessages

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: WENO35reconstruction

    ! Local variables
    real(wp), parameter :: gam(3) = [ 1.0_wp/10.0_wp, 3.0_wp/5.0_wp, &
                                      3.0_wp/10.0_wp ]
    real(wp), parameter :: eps    = 1e-6_wp

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> WENO approximation of the function values at the interafaces, according to
!> Shu, Chi-Wang (2009). High order weighted essentially nonoscillatory schemes
!> for convection dominated problems. SIAM Review, 51(1), 82–126.
!
!> @param[in]     nCons     number of constraints (equations)
!> @param[in]     elems     list of elements in the mesh
!> @param[in]     nFaces    number of faces
!> @param[inout]  faces     list of faces in the mesh
!> @param[in]     periodic  .true. if the BCs are periodic
!···············································································
subroutine WENO35reconstruction(nCons, elems, nFaces, faces, periodic)
    !* Arguments *!
    integer, intent(in) :: nCons
    integer, intent(in) :: nFaces
    logical, intent(in) :: periodic
    ! Derived types
    type(Elem_t), intent(in)    :: elems(:)
    type(Face_t), intent(inout) :: faces(:)

    !* Local variables *!
    integer  :: i
    integer  :: j
    integer  :: k
    integer  :: nLeft
    integer  :: nRight
    integer  :: leftElem
    real(wp) :: Phi(6, nCons)
    real(wp) :: st(2, 3, nCons)
    real(wp) :: w(2, 3, nCons)

    ! First, check that there are at least 6 faces (5 elements)
    if (nFaces < 6) then
        call printError("WENO35", &
                        "There must be at least 5 elements in the mesh.")
    end if

    ! Loop over all the faces
    do i = 1, size(faces)

        ! Skip undefined faces ('holes')
        if (faces(i)%ID < 0) then
            cycle
        end if

        ! Count how many elements are available for interpolation
        if (periodic) then
            nLeft  = 3
            nRight = 3
        else
            nLeft  = min(faces(i)%ID-1, 3)
            nRight = min(nFaces-faces(i)%ID, 3)
        end if

        ! Reference to the leftmost element of the stencil
        if (nLeft > 0) then
            leftElem = faces(i)%elemLeft
            do j = 1, nLeft-1
                leftElem = elems(leftElem)%elemLeft
            end do

        else if (nRight > 0) then
            leftElem = faces(i)%elemRight

        else
            call printError("WENO35", "There must be at least one element.")

        end if

        ! Fill in the 'Phi' values
        j = 3 - nLeft + 1
        do k = 1, nLeft + nRight
            Phi(j,:) = elems(leftElem)%Phi(1,:)   ! Only one point!
            leftElem = elems(leftElem)%elemRight
            j        = j + 1
        end do

        ! Calculate the values for each stencil and its weight
        call getStencils(Phi, nLeft, nRight, st, w)

        ! Reconstruct the values from the available stencils
        faces(i)%PhiL = w(eLeft,1,:) * st(eLeft,1,:) + &
                        w(eLeft,2,:) * st(eLeft,2,:) + &
                        w(eLeft,3,:) * st(eLeft,3,:)

        faces(i)%PhiR = w(eRight,1,:) * st(eRight,1,:) + &
                        w(eRight,2,:) * st(eRight,2,:) + &
                        w(eRight,3,:) * st(eRight,3,:)

    end do

end subroutine WENO35reconstruction

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> WENO stencils and sensor values according to Shu, Chi-Wang (2009). High
!  order weighted essentially nonoscillatory schemes for convection dominated
!  problems. SIAM Review, 51(1), 82–126.
!
!> @param[in]   Phi     average values at the elements of the stencil
!> @param[in]   nLeft   number of elements to the left of the face (max 3)
!> @param[in]   nRight  number of elements to the right of the face (max 3)
!> @param[out]  st      values at the interface of the stencils
!> @param[out]  weight  weights for the stencils
!···············································································
subroutine getStencils(Phi, nLeft, nRight, st, weight)
    !* Arguments *!
    integer,  intent(in)  :: nLeft
    integer,  intent(in)  :: nRight
    real(wp), intent(in)  :: Phi(:,:)
    real(wp), intent(out) :: st(:,:,:)
    real(wp), intent(out) :: weight(:,:,:)

    !* Local variables *!
    integer               :: i
    real(wp), allocatable :: weightSum(:)
    real(wp), allocatable :: beta(:,:,:)

    ! Initialise all the output values
    st     = 0.0_wp
    weight = 0.0_wp
    allocate(beta, source=weight)

    ! Positive wave propagation direction (coming from the left)
    if (nRight >= 2 .and. nLeft >= 1) then
        call stencil_3(Phi, eLeft, st, beta)
        weight(eLeft,3,:) = 1.0_wp

        if (nLeft >= 3) then
            call stencil_1(Phi, eLeft, st, beta)
            call stencil_2(Phi, eLeft, st, beta)

            ! WENO35 weights
            do i = 1, 3
                weight(eLeft,i,:) = gam(i) / (eps + beta(eLeft,i,:))**2
            end do

            weightSum = sum(weight(eLeft,:,:), 1)

            do i = 1, 3
                weight(eLeft,i,:) = weight(eLeft,i,:) / weightSum
            end do

        end if

    else if (nLeft >= 3) then
        call stencil_1(Phi, eLeft, st, beta)
        weight(eLeft,1,:) = 1.0_wp

    end if

    ! Negative wave propagation direction (coming from the right)
    if (nLeft >= 2 .and. nRight >= 1) then
        call stencil_3(Phi, eRight, st, beta)
        weight(eRight,3,:) = 1.0_wp

        if (nRight >= 3) then
            call stencil_1(Phi, eRight, st, beta)
            call stencil_2(Phi, eRight, st, beta)

            ! WENO35 weights
            do i = 1, 3
                weight(eRight,i,:) = gam(i) / (eps + beta(eRight,i,:))**2
            end do

            weightSum = sum(weight(eRight,:,:), 1)

            do i = 1, 3
                weight(eRight,i,:) = weight(eRight,i,:) / weightSum
            end do

        end if

    else if (nRight >= 3) then
        call stencil_1(Phi, eRight, st, beta)
        weight(eRight,1,:) = 1.0_wp

    end if

end subroutine getStencils

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> WENO stencil 1 and its sensor accoording to Shu, Chi-Wang (2009). High
!  order weighted essentially nonoscillatory schemes for convection dominated
!  problems. SIAM Review, 51(1), 82–126.
!
!> @param[in]   Phi     average values at the elements of the stencil
!> @param[in]   side    whether the stencil is for $\Phi_L$ or $\Phi_R$
!> @param[out]  st      values at the interface of the stencils
!> @param[out]  weight  weights for the stencils
!···············································································
subroutine stencil_1(Phi, side, st, beta)
    !* Arguments *!
    real(wp) ,intent(in)  :: Phi(:,:)
    integer,  intent(in)  :: side
    real(wp), intent(out) :: st(:,:,:)
    real(wp), intent(out) :: beta(:,:,:)

    if (side == eLeft) then

        st(eLeft,1,:) = 1.0_wp/3.0_wp*Phi(1,:) - 7.0_wp/6.0_wp*Phi(2,:) &
                      + 11.0_wp/6.0_wp*Phi(3,:)

        beta(eLeft,1,:) = 13.0_wp/12.0_wp*( Phi(1,:) - 2.0_wp*Phi(2,:)   &
                        + Phi(3,:) )**2 + 1.0_wp/4.0_wp * ( Phi(1,:)     &
                        - 4.0_wp*Phi(2,:) + 3.0_wp*Phi(3,:) )**2

    else if (side == eRight) then

        st(eRight,1,:) = 1.0_wp/3.0_wp*Phi(6,:) - 7.0_wp/6.0_wp*Phi(5,:) &
                       + 11.0_wp/6.0_wp*Phi(4,:)

        beta(eRight,1,:) = 13.0_wp/12.0_wp*( Phi(6,:) - 2.0_wp*Phi(5,:) &
                         + Phi(4,:) )**2 + 1.0_wp/4.0_wp*( Phi(6,:)     &
                         - 4.0_wp*Phi(5,:) + 3.0_wp*Phi(4,:) )**2
    end if

end subroutine stencil_1

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> WENO stencil 2 and its sensor accoording to Shu, Chi-Wang (2009). High
!  order weighted essentially nonoscillatory schemes for convection dominated
!  problems. SIAM Review, 51(1), 82–126.
!
!> @param[in]   Phi     average values at the elements of the stencil
!> @param[in]   side    whether the stencil is for $\Phi_L$ or $\Phi_R$
!> @param[out]  st      values at the interface of the stencils
!> @param[out]  weight  weights for the stencils
!···············································································
subroutine stencil_2(Phi, side, st, beta)
    !* Arguments *!
    real(wp), intent(in)  :: Phi(:,:)
    integer,  intent(in)  :: side
    real(wp), intent(out) :: st(:,:,:)
    real(wp), intent(out) :: beta(:,:,:)

    if (side == eLeft) then

        st(eLeft,2,:) = -1.0_wp/6.0_wp*Phi(2,:) + 5.0_wp/6.0_wp*Phi(3,:) &
                      + 1.0_wp/3.0_wp*Phi(4,:)

        beta(eLeft,2,:) = 13.0_wp/12.0_wp*( Phi(2,:) - 2.0_wp*Phi(3,:) &
                        + Phi(4,:) )**2 + 1.0_wp/4.0_wp*( Phi(2,:)-Phi(4,:) )**2

    else if (side == eRight) then

        st(eRight,2,:) = -1.0_wp/6.0_wp*Phi(5,:) + 5.0_wp/6.0_wp*Phi(4,:) &
                       + 1.0_wp/3.0_wp*Phi(3,:)

        beta(eRight,2,:) = 13.0_wp/12.0_wp*( Phi(5,:) - 2.0_wp*Phi(4,:) &
                         +Phi(3,:) )**2 + 1.0_wp/4.0_wp*( Phi(5,:)-Phi(3,:) )**2
    end if

end subroutine stencil_2

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> WENO stencil 3 and its sensor accoording to Shu, Chi-Wang (2009). High
!  order weighted essentially nonoscillatory schemes for convection dominated
!  problems. SIAM Review, 51(1), 82–126.
!
!> @param[in]   Phi     average values at the elements of the stencil
!> @param[in]   side    whether the stencil is for $\Phi_L$ or $\Phi_R$
!> @param[out]  st      values at the interface of the stencils
!> @param[out]  weight  weights for the stencils
!···············································································
subroutine stencil_3(Phi, side, st, beta)
    !* Arguments *!
    real(wp), intent(in)  :: Phi(:,:)
    integer,  intent(in)  :: side
    real(wp), intent(out) :: st(:,:,:)
    real(wp), intent(out) :: beta(:,:,:)

    if (side == eLeft) then

        st(eLeft,3,:) = 1.0_wp/3.0_wp*Phi(3,:) + 5.0_wp/6.0_wp*Phi(4,:) &
                      - 1.0_wp/6.0_wp*Phi(5,:)

        beta(eLeft,3,:) = 13.0_wp/12.0_wp*( Phi(3,:) - 2.0_wp*Phi(4,:)   &
                        +Phi(5,:) )**2 + 1.0_wp/4.0_wp*( 3.0_wp*Phi(3,:) &
                        - 4.0_wp*Phi(4,:) + Phi(5,:) )**2

    else if (side == eRight) then

        st(eRight,3,:) = 1.0_wp/3.0_wp*Phi(4,:) + 5.0_wp/6.0_wp*Phi(3,:) &
                       - 1.0_wp/6.0_wp*Phi(2,:)

        beta(eRight,3,:) = 13.0_wp/12.0_wp*( Phi(4,:) - 2.0_wp*Phi(3,:)    &
                         + Phi(2,:) )**2 + 1.0_wp/4.0_wp*( 3.0_wp*Phi(4,:) &
                         - 4.0_wp*Phi(3,:) + Phi(2,:) )**2
    end if

end subroutine stencil_3

end module WENO35
