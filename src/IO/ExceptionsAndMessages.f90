!*******************************************************************************
!  MODULE: ExceptionsAndMessages
!
!> @author
!> Andres Mateo
!
!> @brief
!> Useful routines to handle errors and warnings.
!*******************************************************************************

module ExceptionsAndMessages

    use Constants

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: printError
    public :: printWarning
    public :: printInfo
    public :: assertError
    public :: assertWarning
    public :: assertInfo

    ! Local variables
    character, parameter :: ESC = achar(27)

contains

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Prints an error message and stops the execution of the program.
!
!> @param[in]  failedMod  name of the module where the error took place
!> @param[in]  msg        error message to print
!···············································································
subroutine printError(failedMod, msg)
    !* Arguments *!
    character(len=*), intent(in) :: failedMod
    character(len=*), intent(in) :: msg

    write(stderr,'(a,a,a,a)') '['//ESC//"[31mError"//ESC//"[0m] ", &
                              failedMod, " -> ", msg
#ifndef NDEBUG
    error stop
#else
    stop
#endif

end subroutine printError

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Prints a warning message.
!
!> @param[in]  failedMod  name of the module where the warning took place
!> @param[in]  msg        warning message to print
!···············································································
subroutine printWarning(failedMod, msg)
    !* Arguments *!
    character(len=*), intent(in) :: failedMod
    character(len=*), intent(in) :: msg

    write(stderr,'(a,a,a,a)') '['//ESC//"[33mWarning"//ESC//"[0m] ", &
                              failedMod, ' -> ', msg

end subroutine printWarning

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Prints an informative message.
!
!> @param[in]  infoMod  name of the module where the error took place
!> @param[in]  msg      info message to print
!···············································································
subroutine printInfo(infoMod, msg)
    !* Arguments *!
    character(len=*), intent(in) :: infoMod
    character(len=*), intent(in) :: msg

    write(stdout,'(a,a,a,a)') '['//ESC//"[32mInfo"//ESC//"[0m] ", &
                              infoMod, " -> ", msg

end subroutine printInfo

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Prints an error message and stops the execution of the program (only debug).
!
!> @param[in]  failedMod  name of the module where the error took place
!> @param[in]  msg        error message to print
!···············································································
subroutine assertError(failedMod, msg)
    !* Arguments *!
    character(len=*), intent(in) :: failedMod
    character(len=*), intent(in) :: msg

#ifndef NDEBUG
    write(stderr,'(a,a,a,a)') '['//ESC//"[31mAssert Error"//ESC//"[0m] ", &
                              failedMod, " -> ", msg
    error stop
#endif

end subroutine assertError

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Prints a warning message (only debug).
!
!> @param[in]  failedMod  name of the module where the warning took place
!> @param[in]  msg        warning message to print
!···············································································
subroutine assertWarning(failedMod, msg)
    !* Arguments *!
    character(len=*), intent(in) :: failedMod
    character(len=*), intent(in) :: msg

#ifndef NDEBUG
    write(stderr,'(a,a,a,a)') '['//ESC//"[33mAssert Warning"//ESC//"[0m] ", &
                              failedMod, ' -> ', msg
#endif

end subroutine assertWarning

!···············································································
!> @author
!> Andres Mateo
!
!  DESCRIPTION
!> @brief
!> Prints an informative message (only debug).
!
!> @param[in]  infoMod  name of the module where the error took place
!> @param[in]  msg      info message to print
!···············································································
subroutine assertInfo(infoMod, msg)
    !* Arguments *!
    character(len=*), intent(in) :: infoMod
    character(len=*), intent(in) :: msg

#ifndef NDEBUG
    write(stdout,'(a,a,a,a)') '['//ESC//"[32mAssert Info"//ESC//"[0m] ", &
                              infoMod, " -> ", msg
#endif

end subroutine assertInfo

end module ExceptionsAndMessages
