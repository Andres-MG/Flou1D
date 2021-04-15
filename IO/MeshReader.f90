module MeshReader

    use Constants, only: CHAR_LEN
    use MeshClass, only: Mesh_t
    use ExceptionsAndMessages, only: printError

    implicit none

contains

function Read_mesh(filename, isHO) result(mesh)
    !* Arguments *!
    character(len=*), intent(in) :: filename
    logical,          intent(in) :: isHO

    !* Return value *!
    type(Mesh_t) :: mesh

    !* Local variables *!
    integer                 :: fd
    integer                 :: ios
    character(len=CHAR_LEN) :: line

    integer :: pos


    ! Open the file
    open(newunit=fd, file=trim(adjustl(filename)), &
         status="old", action="read", iostat=ios)
    if (ios /= 0) then
        call printError("MeshReader.f90", &
                        "The mesh file could not be opened.")
    end if

    ! Read the file
    do

        read(fd,'(a)') line
        call to_lower(line)

        ! First, find the main ZONE block

    end do

    ! Close the file
    close(fd)

end function Read_mesh

pure subroutine to_lower(string)
    !* Arguments *!
    character(len=*), intent(inout) :: string

    !* Local variables *!
    character(len=*), parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
    character(len=*), parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    integer :: pos
    integer :: i


    do i = 1, len_trim(string)
        pos = index(UPPER_CASE, string(i:i))
        if (pos > 0) string(i:i) = LOWER_CASE(pos:pos)
    end do

end subroutine to_lower

end module MeshReader
