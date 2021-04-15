module GradientOperators

    use Constants, only: wp, eLeft, eRight, eGaussLobatto
    use Physics,   only: NEQS, EulerFlux, TwoPointFlux
    use StdExpansionClass
    use StdExpansions
    use MeshClass
    use ElementClass, only: Elem_t

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: GradientDGSEMweak
    public :: GradientDGSEMstrong
    public :: DGSEMgradient_Int

    abstract interface
        function DGSEMgradient_Int(mesh, F, FstarLeft, FstarRight, &
                                   Fleft, Fright, elemID)
            import wp, Mesh_t
            type(Mesh_t), intent(inout) :: mesh
            real(wp),     intent(in)    :: F(:,:)
            real(wp),     intent(in)    :: FstarLeft(:)
            real(wp),     intent(in)    :: FstarRight(:)
            real(wp),     intent(in)    :: Fleft(:)
            real(wp),     intent(in)    :: Fright(:)
            integer,      intent(in)    :: elemID
            real(wp), allocatable       :: DGSEMgradient_Int(:,:)
        end function DGSEMgradient_Int
    end interface

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Gradient operator (weak) for the DGSEM.
!
!> @param[in]  mesh        mesh with all the elements of the domain
!> @param[in]  F           function to calculate the gradient
!> @param[in]  FstarLeft   numerical flux at the left boundary
!> @param[in]  FstarRight  numerical flux at the right boundary
!> @param[in]  Fleft       flux at the left boundary
!> @param[in]  Fright      flux at the right boundary
!> @param[in]  elemID      ID of the element where the operator acts
!
!> @return                 gradient (weak operator)
!···············································································
function GradientDGSEMweak(mesh, F, FstarLeft, FstarRight, &
                           Fleft, Fright, elemID) result(GradF)
    !* Arguments *!
    real(wp), intent(in) :: F(:,:)
    real(wp), intent(in) :: FstarLeft(:)
    real(wp), intent(in) :: FstarRight(:)
    real(wp), intent(in) :: FLeft(:)
    real(wp), intent(in) :: FRight(:)
    integer,  intent(in) :: elemID
    ! Derived types
    type(Mesh_t), intent(inout) :: mesh

    !* Return value *!
    real(wp), allocatable :: GradF(:,:)

    !* Local variabes *!
    integer :: i
    integer :: n
    type(Elem_t), pointer :: elem


    !---- Begin associate ----!
    elem => mesh%elems%at(elemID)
    associate(std => elem%std)

    n = std%n

    ! Interior values
    GradF = matmul(std%Dh, F)

    ! Add boundary flux
    if (std%hasBounds) then

        GradF(1,:) = GradF(1,:) - FstarLeft(:)  * std%iw(1)
        GradF(n,:) = GradF(n,:) + FstarRight(:) * std%iw(n)

    else

        do i = 1, NEQS
            GradF(:,i) = GradF(:,i) + (FstarRight(i) * std%bdMode(eRight,:) &
                       - FstarLeft(i) * std%bdMode(eLeft,:)) * std%iw
       end do

    end if

    end associate
    !----- End associate -----!

    nullify(elem)

end function GradientDGSEMweak

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Gradient operator (strong) for the DGSEM.
!
!> @param[in]  mesh        mesh with all the elements of the domain
!> @param[in]  F           function to calculate the gradient
!> @param[in]  FstarLeft   numerical flux at the left boundary
!> @param[in]  FstarRight  numerical flux at the right boundary
!> @param[in]  Fleft       flux at the left boundary
!> @param[in]  Fright      flux at the right boundary
!> @param[in]  elemID      ID of the element where the operator acts
!
!> @return                 gradient (weak operator)
!···············································································
function GradientDGSEMstrong(mesh, F, FstarLeft, FstarRight,&
                             Fleft, Fright, elemID) result(GradF)
    !* Arguments *!
    real(wp), intent(in) :: F(:,:)
    real(wp), intent(in) :: FstarLeft(:)
    real(wp), intent(in) :: FstarRight(:)
    real(wp), intent(in) :: Fleft(:)
    real(wp), intent(in) :: Fright(:)
    integer,  intent(in) :: elemID
    ! Derived types
    type(Mesh_t), intent(inout) :: mesh

    !* Return value *!
    real(wp), allocatable :: GradF(:,:)

    !* Local variabes *!
    integer :: i
    integer :: n
    type(Elem_t), pointer :: elem


    !---- Begin associate ----!
    elem => mesh%elems%at(elemID)
    associate(std => elem%std)

    n = std%n

    ! Interior values
    GradF = matmul(std%D, F)

    ! Add boundary flux
    if (std%hasBounds) then

        GradF(1,:) = GradF(1,:) - (FstarLeft-Fleft)   * std%iw(1)
        GradF(n,:) = GradF(n,:) + (FstarRight-Fright) * std%iw(n)

    else

        do i = 1, NEQS
            GradF(:,i) = GradF(:,i) + ( (FstarRight(i)-Fright(i))       &
                       * std%bdMode(eRight,:) - (FstarLeft(i)-Fleft(i)) &
                       * std%bdMode(eLeft,:) ) * std%iw
        end do

    end if

    end associate
    !----- End associate -----!

    nullify(elem)

end function GradientDGSEMstrong

end module GradientOperators
