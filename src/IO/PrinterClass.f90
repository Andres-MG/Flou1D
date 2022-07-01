!*******************************************************************************
!  MODULE: PrinterClass
!
!> @author
!> Andres Mateo
!
!> @brief
!> Wrapper for the different input/output subroutines.
!*******************************************************************************

module PrinterClass

    use ExceptionsAndMessages
    use Constants
    use Utilities
    use Lagrange
    use Gnuplot
    use Tecplot
    use PDEclass
    use Physics
    use ElementClass, only: Elem_t
    use StdExpansionClass

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public types
    public :: Printer_t

    ! Explicitly define public variables
    public :: Printer

!···············································································
!> @class Printer_t
!
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Prints data to the screen and saves the results to the solution files.
!···············································································
    type :: Printer_t
        integer                 :: outResolution
        integer                 :: fileType
        real(wp), allocatable   :: fineNodes(:)
        character(len=CHAR_LEN) :: meshFile
        character(len=CHAR_LEN) :: mapFile
        character(len=CHAR_LEN) :: extension
        procedure(SaveMesh_Int), &
            pointer, nopass :: saveMesh_deferred
        procedure(SaveMap_Int), &
            pointer, nopass :: saveMap_deferred
    contains
        procedure :: construct  => Printer_constructor
        procedure :: saveMesh   => Printer_save_mesh
        procedure :: saveTEmap  => Printer_save_TE_map
        procedure :: endSummary => Printer_end_summary
    end type Printer_t

    abstract interface
        subroutine SaveMesh_Int(fileUnit, baseName, outCoords, vars, time, &
                                interpolate)
            import wp
            integer,          intent(in) :: fileUnit
            character(len=*), intent(in) :: baseName
            real(wp),         intent(in) :: outCoords(:)
            character(len=*), intent(in) :: vars(:)
            real(wp),         intent(in) :: time
            logical,          intent(in) :: interpolate
        end subroutine

        subroutine SaveMap_Int(fileUnit, baseName, vars, time)
            import wp
            integer,          intent(in) :: fileUnit
            character(len=*), intent(in) :: baseName
            character(len=*), intent(in) :: vars(:)
            real(wp),         intent(in) :: time
        end subroutine
    end interface

    type(Printer_t) :: Printer

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Saves the results stored in 'mesh' to a solution file with name
!> 'fileName_time.extension'.
!
!> @param[in]  meshName  name of the solution file
!> @param[in]  outRes    resolution of the output mesh in each element
!···············································································
subroutine Printer_constructor(this, meshName, outRes)
    !* Arguments *!
    character(len=*), intent(in) :: meshName
    integer,          intent(in) :: outRes
    ! Derived types
    class(Printer_t), intent(inout) :: this

    !* Local variables *!
    integer :: i
    integer :: extPosition
    character(len=:), allocatable :: baseName
    character(len=:), allocatable :: extension

    ! Extract extension
    extPosition = index(meshName, '.', back=.true.)
    baseName    = trim(adjustl(meshName(:extPosition-1)))
    extension   = trim(adjustl(meshName(extPosition:)))

    select case (extension)

    case (".gp")
        this%fileType = eGnuplot
        this%saveMesh_deferred => Gnuplot_saveMesh
        this%saveMap_deferred  => Gnuplot_saveMap

    case (".dat")
        this%fileType = eTecPlot
        this%saveMesh_deferred => TecPlot_saveMesh
        this%saveMap_deferred  => TecPlot_saveMap

    case default
        call printError("PrinterClass.f90", &
                        "Implemented export formats are '.gp' and '.dat'.")

    end select

    ! Store the extension and the files basename (time will be added later)
    this%extension = extension
    this%meshFile  = baseName
    this%mapFile   = baseName // "_TEmap_"

    ! Save output resolution
    this%outResolution = outRes

    if (outRes > 0) then
        this%fineNodes = [(-1.0_wp + 2.0_wp * (i-1)/(outRes-1), i = 1, outRes)]
    end if

end subroutine Printer_constructor

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Saves the results stored in 'mesh' to a solution file with name
!> 'meshFile_time.extension'.
!
!> @param[in]  time  time instant of the solution to be saved
!···············································································
subroutine Printer_save_mesh(this, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(Printer_t), intent(inout) :: this

    !* Local variables *!
    character(len=CHAR_LEN), allocatable :: keys(:)
    character(len=CHAR_LEN)              :: solutionFile
    character(len=CHAR_LEN)              :: infoMsg
    character(len=CHAR_LEN)              :: tfmt
    integer                              :: i
    integer                              :: unitMesh
    integer                              :: outRes
    logical                              :: interpolate
    real(wp),                allocatable :: outNodes(:)
    type(Elem_t),            pointer     :: elem
    class(StdExp_t),         pointer     :: std

    ! Initialization
    if (time < 1.0_wp) then
        tfmt = TIME0_FMT
    else
        tfmt = TIME1_FMT
    end if

    ! Allocate/update the interpolation matrices of the elements
    call PDE%mesh%elems%reset_last(0)
    do i = 1, PDE%mesh%elems%size()

        elem => PDE%mesh%elems%next()
        std  => elem%std

        ! Set the output resolution
        if (this%outResolution > 0) then
            outRes   = this%outResolution
            outNodes = this%fineNodes
            interpolate = .true.
        else
            outRes   = std%n
            outNodes = std%x
            interpolate = .false.
        end if

        ! Allocate if necessary
        if (.not. allocated(std%fineInp)) then

            allocate(std%fineInp(outRes, std%n))

            ! Compute interpolation matrix
            call polynomial_interpolation_matrix(std%n, std%x, std%wb, &
                                                 outRes, outNodes, std%fineInp)
        end if

    end do

    nullify(elem)
    nullify(std)

    ! Create 'keys' for the values to be saved (coords, vars)
    ! x, vars, visc, sensor, P, exact vars, S, entropy (s)
    allocate(keys(3*NEQS+5))
    keys(1)                 = 'x'
    keys(2:NEQS+1)          = VarNames
    keys(NEQS+2:NEQS+4)     = [character(len=CHAR_LEN) :: "artificial visc.", &
                                                          "sensor", "order" ]
    keys(NEQS+5:2*NEQS+4)   = "exact "//VarNames
    keys(2*NEQS+5:3*NEQS+4) = "source "//VarNames
    keys(3*NEQS+5)          = "entropy"

    ! Set name of the file
    write(solutionFile, '(a,'//tfmt//',a)') trim(this%meshFile)//'_', time, &
                                            trim(this%extension)

    ! And open the file
    open(newunit=unitMesh, file=solutionFile)

    ! Write to the file
    call this%saveMesh_deferred(unitMesh, this%meshFile, outNodes, &
                                keys, time, interpolate)

    ! And close the files
    close(unitMesh)

    ! Notify
    write(infoMsg, '(a,'//tfmt//',2a)') "Saved solution at t = ", time, &
                                        " in ", trim(solutionFile)
    call printInfo("PrinterClass.f90", infoMsg)

end subroutine Printer_save_mesh

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Saves the results stored in 'mesh' to a solution file with name
!> 'meshFile_time.extension'.
!
!> @param[in]  time  time instant of the solution to be saved
!···············································································
subroutine Printer_save_TE_map(this, time)
    !* Arguments *!
    real(wp), intent(in) :: time
    ! Derived types
    class(Printer_t), intent(inout) :: this

    !* Local variables *!
    character(len=CHAR_LEN) :: solutionFile
    character(len=CHAR_LEN) :: infoMsg
    character(len=CHAR_LEN) :: tfmt
    integer                 :: unitMap

    ! Initialization
    if (time < 1.0_wp) then
        tfmt = TIME0_FMT
    else
        tfmt = TIME1_FMT
    end if

    ! Set name of the file
    write(solutionFile, '(a,'//tfmt//',a)') trim(this%mapFile), time, &
                                            trim(this%extension)

    ! And open the file
    open(newunit=unitMap, file=solutionFile)

    ! Write to the file
    call this%saveMap_deferred(unitMap, this%mapFile, VarNames, time)

    ! And close the files
    close(unitMap)

    ! Notify
    write(infoMsg, '(a,'//tfmt//',2a)') "Saved TE map at t = ", time, &
                                        " in ", trim(solutionFile)
    call printInfo("PrinterClass.f90", infoMsg)

end subroutine Printer_save_TE_map

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Print a summary of the computation once it has been completed.
!
!> @param[in]  endTime        last simulated instant
!> @param[in]  executionTime  real cpu time the was running
!···············································································
subroutine Printer_end_summary(this, endTime, executionTime)
    !* Arguments *!
    real(wp), intent(in) :: endTime
    real(wp), intent(in) :: executionTime
    ! Derived types
    class(Printer_t), intent(in) :: this

    ! Just call the summary function of the PDE object
    call PDE%printEndSummary(endTime, executionTime)

end subroutine Printer_end_summary

end module PrinterClass
