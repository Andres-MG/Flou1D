!*******************************************************************************
!  MODULE: Gnuplot
!
!> @author
!> Andres Mateo
!
!> @brief
!> Routines to save data to a gnuplot (.gp) file.
!*******************************************************************************

module Gnuplot

    use Constants
    use Utilities
    use PDEclass
    use TruncationErrorClass
    use Physics
    use Setup_routines
    use ElementClass

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: Gnuplot_saveMesh
    public :: Gnuplot_saveMap

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Writes the variables 'vars' of all the elements of the mesh to the file
!> referred by 'meshUnit'.
!
!> @param[in]  meshUnit     unit of the file where the solution will be saved
!> @param[in]  baseName     name of the case
!> @param[in]  outCoords    nodes of the output mesh in std. space
!> @param[in]  vars         names of the variables to save
!> @param[in]  time         time of the solution
!> @param[in]  interpolate  .true. if the out nodes are not the mesh nodes
!···············································································
subroutine Gnuplot_saveMesh(meshUnit, baseName, outCoords, vars, time, &
                            interpolate)
    !* Arguments *!
    integer,          intent(in) :: meshUnit
    character(len=*), intent(in) :: baseName
    real(wp),         intent(in) :: outCoords(:)
    character(len=*), intent(in) :: vars(:)
    real(wp),         intent(in) :: time
    logical,          intent(in) :: interpolate

    !* Local variables *!
    integer                 :: nvars
    integer                 :: i
    integer                 :: j
    integer                 :: k
    logical                 :: available_exact
    real(wp), allocatable   :: coords(:)
    real(wp), allocatable   :: stdCoords(:)
    real(wp), allocatable   :: newPhi(:,:)
    real(wp), allocatable   :: newExact(:,:)
    real(wp), allocatable   :: newSrc(:,:)
    real(wp), allocatable   :: vals(:,:)
    character(len=CHAR_LEN) :: tfmt
    character(len=CHAR_LEN) :: fmat
    type(Elem_t), pointer   :: elem => null()

    ! Initialise
    nvars = size(vars)
    if (time < 1.0_wp) then
        tfmt = TIME0_FMT
    else
        tfmt = TIME1_FMT
    end if

    ! Write the header of the solution file
    write(meshUnit, '(a)')  "### CFD_SOLVER solution file ###"
    write(meshUnit, '(2a)')  "# Case: ", trim(baseName)
    write(meshUnit, '(a,'//tfmt//')')  "# Time: ", time
    write(meshUnit, '(2a)', advance="no") "# Variables: ", trim(vars(1))
    do i = 2, nvars
        write(meshUnit, '(2a)', advance="no") ", ", trim(vars(i))
    end do
    write(meshUnit, *) ''

    ! Choose saving format for reals
#ifdef SINGLE_PRECISION
    write(fmat, '(a, i0, a)') '(', nvars, 'e16.6)'
#else
    write(fmat, '(a, i0, a)') '(', nvars, 'e26.16)'
#endif

    ! Set the output coords if possible
    if (interpolate) then
        stdCoords = outCoords
    end if

    call PDE%mesh%elems%reset_last(0)
    do i = 1, PDE%mesh%elems%size()

        elem => PDE%mesh%elems%next()

        !---- Begin associate ----!
        associate(std  => elem%std)

        ! Get the output coordinates
        if (.not. interpolate) then
            stdCoords = std%x
        end if

        ! Get the values in the output mesh
        coords = elem%affineMap(stdCoords)
        newPhi = matmul(std%fineInp, elem%Phi)
        call reallocate(size(newPhi,1), size(newPhi,2), newExact)
        call reallocate(size(newPhi,1), size(newPhi,2), newSrc)

        available_exact = exactSolution(coords, time, newExact)
        call sourceTerm(coords, time, newSrc)

        call reallocate(size(newPhi,1), nvars-1, vals)
        vals(:,:NEQS)             = newPhi
        vals(:,NEQS+1)            = elem%aVis
        vals(:,NEQS+2)            = elem%sens
        vals(:,NEQS+3)            = std%n - 1
        vals(:,NEQS+4:2*NEQS+3)   = newExact
        vals(:,2*NEQS+4:3*NEQS+3) = newSrc
        vals(:,3*NEQS+4)          = getPhysicalEntropy(newPhi)

        ! Write header for the element
        write(meshUnit, '(a, i0)') "# ", elem%ID

        ! Save all the nodes (avoid warning for non-contiguous I/O)
        do j = 1, size(newPhi, dim=1)
            write(meshUnit, fmat) coords(j), (vals(j,k), k=1,size(vals,dim=2))
        end do

        ! Add spacing after each element
        write(meshUnit, *) ''
        write(meshUnit, *) ''

        end associate
        !----- End associate -----!

        nullify(elem)

    end do

end subroutine Gnuplot_saveMesh

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Writes the error estimation for different orders in each element.
!
!> @param[in]  mapUnit   unit of the file where the TE map will be saved
!> @param[in]  baseName  name of the case
!> @param[in]  vars      names of the variables to save
!> @param[in]  time      time of the solution
!···············································································
subroutine Gnuplot_saveMap(mapUnit, baseName, vars, time)
    !* Arguments *!
    integer,          intent(in) :: mapUnit
    character(len=*), intent(in) :: baseName
    character(len=*), intent(in) :: vars(:)
    real(wp),         intent(in) :: time

    !* Local variables *!
    integer                 :: i
    integer                 :: j
    integer                 :: k
    character(len=CHAR_LEN) :: tfmt
    character(len=CHAR_LEN) :: fmat
    type(Elem_t), pointer   :: elem => null()

    ! Initialise
    if (time < 1.0_wp) then
        tfmt = TIME0_FMT
    else
        tfmt = TIME1_FMT
    end if

    ! Header of the file
    write(mapUnit, '(a)')  "### CFD_SOLVER map file ###"
    write(mapUnit, '(2a)')  "# Case: ", trim(baseName)
    write(mapUnit, '(a,'//tfmt//')')  "# Time: ", time
    write(mapUnit, '(a)', advance="no") "# Variables: App. order"
    do i = 1, size(vars)
        write(mapUnit, '(3a)', advance="no") ", TE(", trim(vars(i)), ')'
    end do
    write(mapUnit, *) ''

    ! Set output format
#ifdef SINGLE_PRECISOIN
    write(fmat, '(a, i0, a)') '(i0, ', NEQS, 'e16.6)'
#else
    write(fmat, '(a, i0, a)') '(i0, ', NEQS, 'e26.16)'
#endif

    call PDE%mesh%elems%reset_last(0)
    do i = 1, PDE%mesh%elems%size()

        elem => PDE%mesh%elems%at(i)

        ! Write header for the element
        write(mapUnit, '(a, i0)') "# ", elem%ID

        ! Save all the truncation error estimates
        do j = TruncError%MapLimit+1, ubound(elem%Terr, dim=1)
            write(mapUnit, fmat) j-1, (elem%Terr(j,k), k=1,NEQS)
        end do

        ! Add spacing after each element
        write(mapUnit, *) ''
        write(mapUnit, *) ''

    end do

    nullify(elem)

end subroutine Gnuplot_saveMap

end module Gnuplot
