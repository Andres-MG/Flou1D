module Utilities

    use Constants

    implicit none

    ! Defaults to private
    private

    ! Explicitly define public functions
    public :: almost_zero
    public :: logarithmicMean
    public :: linear_LSE_regression
    public :: semilogy_LSE_regression
    public :: loglog_LSE_regression
    public :: sigmoid
    public :: sinRamp
    public :: safeAllocate
    public :: reallocate
    public :: linear_mapping
    public :: linear_mapping_2D
    public :: undivided_differences

    interface sinRamp
        module procedure scalar_sinRamp
        module procedure vector_sinRamp
    end interface sinRamp

    interface safeAllocate
        module procedure safe_allocate_real_vect
        module procedure safe_allocate_real_array
    end interface safeAllocate

    interface reallocate
        module procedure reallocate_real_vect
        module procedure reallocate_real_array
        module procedure reallocate_mold_real_vect
        module procedure reallocate_mold_real_array
    end interface reallocate

    interface linear_mapping
        module procedure scalar_linear_mapping
        module procedure vector_linear_mapping
    end interface linear_mapping

    interface linear_mapping_2D
        module procedure scalar_linear_mapping_2D
        module procedure vector_linear_mapping_2D
    end interface linear_mapping_2D

contains

elemental function almost_zero(a, tol)
    !* Arguments *!
    real(wp),           intent(in) :: a
    real(wp), optional, intent(in) :: tol

    !* Return value *!
    logical :: almost_zero

    !* Local variables *!
    real(wp) :: tol_

    tol_ = merge(tol, 2.0_wp*FLT_EPS, present(tol))

    ! Absolute difference
    if (abs(a) <= tol_) then
        almost_zero = .true.
    else
        almost_zero = .false.
    end if
end function almost_zero

elemental function logarithmicMean(xL, xR) result(xM)
    !* Arguments *!
    real(wp), intent(in)  :: xL
    real(wp), intent(in)  :: xR

    !* Return values *!
    real(wp) :: xM

    !* Local variables *!
    real(wp), parameter :: eps = 0.01_wp  ! Magic constant!!
    real(wp) :: xi
    real(wp) :: f
    real(wp) :: f2
    real(wp) :: FF

    ! Set some values
    xi = xL / xR
    f  = (xi - 1.0_wp) / (xi + 1.0_wp)
    f2 = f * f

    ! Get the value of 'FF'
    if (f2 < eps) then

        ! Use a series expansion if f $xi \approx 1$
        FF = 1.0_wp + f2/3.0_wp + f2*f2/5.0_wp + f2*f2*f2/7.0_wp

    else

        ! Go with the natural logarithm otherwise
        FF = log(xi) / (2.0_wp * f)

    end if

    ! Final result
    xM = 0.5_wp * (xL + xR) / FF

end function logarithmicMean

subroutine linear_LSE_regression(x, y, m, n)
    !* Arguments *!
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: y(:)
    real(wp), intent(out) :: m
    real(wp), intent(out) :: n

    !* Local variables *!
    integer  :: nPts
    real(wp) :: sumX
    real(wp) :: sumY
    real(wp) :: sumXY
    real(wp) :: sumSquaresX
    real(wp) :: Ainv(2,2)
    real(wp) :: v(2)
    real(wp) :: b(2)

    ! Get size of the data set
    nPts = size(x)

    ! LSE matrix
    sumX        = sum(x)
    sumY        = sum(y)
    sumXY       = sum(x*y)
    sumSquaresX = sum(x*x)

    ! Setting up of the linear system
    Ainv(1,1) = nPts;  Ainv(1,2) = -sumX
    Ainv(2,1) = -sumX; Ainv(2,2) = sumSquaresX
    Ainv = Ainv / (nPts*sumSquaresX - sumX**2)

    b(1)   = sumXY
    b(2)   = sumY

    ! Linear system for the coefficients
    v = matmul(Ainv, b)
    m = v(1)
    n = v(2)

end subroutine linear_LSE_regression

subroutine semilogy_LSE_regression(x, y, m, n)
    !* Arguments *!
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: y(:)
    real(wp), intent(out) :: m
    real(wp), intent(out) :: n

    ! Pass the natural logarithm of 'y'
    call linear_LSE_regression(x, log(y), m, n)

    ! Undo the logarithm of the multiplicative constant
    n = exp(n)

end subroutine semilogy_LSE_regression

subroutine loglog_LSE_regression(x, y, m, n)
    !* Arguments *!
    real(wp), intent(in)  :: x(:)
    real(wp), intent(in)  :: y(:)
    real(wp), intent(out) :: m
    real(wp), intent(out) :: n

    ! Pass the natural logarithm of 'y'
    call linear_LSE_regression(log(x), log(y), m, n)

    ! Undo the logarithm of the multiplicative constant
    n = exp(n)

end subroutine loglog_LSE_regression

elemental function sigmoid(x)
    !* Arguments *!
    real(wp), intent(in) :: x

    !* Return values *!
    real(wp) :: sigmoid

    sigmoid = 1.0_wp / ( 1.0_wp + exp(-x) )

end function sigmoid

function scalar_sinRamp(x, dx)
    !* Arguments *!
    real(wp), intent(in) :: x
    real(wp), intent(in) :: dx

    !* Return values *!
    real(wp) :: scalar_sinRamp

    !* Local variables *!
    real(wp) :: halfDx

    halfDx = 0.5_wp * dx

    if (x < -halfDx) then
        scalar_sinRamp = 0.0_wp
    else if ( x > halfDx) then
        scalar_sinRamp = 1.0_wp
    else
        scalar_sinRamp = 0.5_wp * ( 1.0_wp + sin(PI*x/dx) )
    end if

end function scalar_sinRamp

function vector_sinRamp(x, dx)
    !* Arguments *!
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: dx

    !* Return values *!
    real(wp) :: vector_sinRamp(size(x))

    !* Local variables *!
    real(wp) :: halfDx

    halfDx = 0.5_wp * dx

    where (x < -halfDx)
        vector_sinRamp = 0.0_wp
    elsewhere ( x > halfDx)
        vector_sinRamp = 1.0_wp
    elsewhere
        vector_sinRamp = 0.5_wp * ( 1.0_wp + sin(PI*x/dx) )
    end where

end function vector_sinRamp

subroutine safe_allocate_real_vect(n, x)
    !* Arguments *!
    integer,               intent(in)    :: n
    real(wp), allocatable, intent(inout) :: x(:)

    if (.not. allocated(x)) then
        allocate(x(n))
    end if

end subroutine safe_allocate_real_vect

subroutine safe_allocate_real_array(m, n, A)
    !* Arguments *!
    integer,               intent(in)    :: m
    integer,               intent(in)    :: n
    real(wp), allocatable, intent(inout) :: A(:,:)

    if (.not. allocated(A)) then
        allocate(A(m, n))
    end if

end subroutine safe_allocate_real_array

subroutine reallocate_real_vect(n, x)
    !* Arguments *!
    integer,               intent(in)    :: n
    real(wp), allocatable, intent(inout) :: x(:)

    if (.not. allocated(x)) then
        allocate(x(n))

    else if (size(x) /= n) then
        deallocate(x)
        allocate(x(n))

    end if

end subroutine reallocate_real_vect

subroutine reallocate_real_array(m, n, A)
    !* Arguments *!
    integer,               intent(in)    :: m
    integer,               intent(in)    :: n
    real(wp), allocatable, intent(inout) :: A(:,:)

    if (.not. allocated(A)) then
        allocate(A(m,n))

    else if (size(A,1) /= m .or. size(A,2) /= n) then
        deallocate(A)
        allocate(A(m,n))

    end if

end subroutine reallocate_real_array

subroutine reallocate_mold_real_vect(x, mld)
    !* Arguments *!
    real(wp), allocatable, intent(inout) :: x(:)
    real(wp),              intent(in)    :: mld(:)

    if (.not. allocated(x)) then
        allocate(x, mold=mld)

    else if (size(x) /= size(mld)) then
        deallocate(x)
        allocate(x, mold=mld)

    end if

end subroutine reallocate_mold_real_vect

subroutine reallocate_mold_real_array(A, mld)
    !* Arguments *!
    real(wp), allocatable, intent(inout) :: A(:,:)
    real(wp),              intent(in)    :: mld(:,:)

    if (.not. allocated(A)) then
        allocate(A, mold=mld)

    else if (size(A,1) /= size(mld,1) .or. size(A,2) /= size(mld,2)) then
        deallocate(A)
        allocate(A, mold=mld)

    end if

end subroutine reallocate_mold_real_array

subroutine scalar_linear_mapping(a, b, z, x)
    !* Arguments *!
    real(wp), intent(in)  :: a
    real(wp), intent(in)  :: b
    real(wp), intent(in)  :: z
    real(wp), intent(out) :: x

    ! Mapped value
    x = 0.5_wp * ((b-a)*z + a+b)
end subroutine scalar_linear_mapping

subroutine vector_linear_mapping(a, b, z, x)
    !* Arguments *!
    real(wp), intent(in)  :: a
    real(wp), intent(in)  :: b
    real(wp), intent(in)  :: z(:)
    real(wp), intent(out) :: x(:)

    ! Mapped value
    x = 0.5_wp * ((b-a)*z + a+b)
end subroutine vector_linear_mapping

subroutine scalar_linear_mapping_2D(xlims, ylims, zx, zy, x)
    !* Arguments *!
    real(wp), intent(in)  :: xlims(4)
    real(wp), intent(in)  :: ylims(4)
    real(wp), intent(in)  :: zx
    real(wp), intent(in)  :: zy
    real(wp), intent(out) :: x(2)

    !* Local variables *!
    real(wp) :: x1(2)  !< Lower left
    real(wp) :: x2(2)  !< Lower right
    real(wp) :: x3(2)  !< Upper right
    real(wp) :: x4(2)  !< Upper left

    ! Initialisation
    x1 = [xlims(1), ylims(1)]
    x2 = [xlims(2), ylims(2)]
    x3 = [xlims(3), ylims(3)]
    x4 = [xlims(4), ylims(4)]

    x = 0.25_wp * (x1*(1.0_wp-zx)*(1.0_wp-zy) + x2*(1.0_wp+zx)*(1.0_wp-zy) &
      + x3*(1.0_wp+zx)*(1.0_wp+zy) + x4*(1.0_wp-zx)*(1.0_wp+zy))

end subroutine scalar_linear_mapping_2D

subroutine vector_linear_mapping_2D(xlims, ylims, zx, zy, x)
    !* Arguments *!
    real(wp), intent(in)  :: xlims(4)
    real(wp), intent(in)  :: ylims(4)
    real(wp), intent(in)  :: zx(:)
    real(wp), intent(in)  :: zy(:)
    real(wp), intent(out) :: x(:,:)  ! n rows, 2 columns

    !* Local variables *!
    integer  :: n
    integer  :: i
    real(wp) :: x1(2)  !< Lower left
    real(wp) :: x2(2)  !< Lower right
    real(wp) :: x3(2)  !< Upper right
    real(wp) :: x4(2)  !< Upper left

    ! Initialisation
    n  = size(zx)
    x1 = [xlims(1), ylims(1)]
    x2 = [xlims(2), ylims(2)]
    x3 = [xlims(3), ylims(3)]
    x4 = [xlims(4), ylims(4)]

    do i = 1, n
        x(i,:) = (x1*(1.0_wp-zx(i))*(1.0_wp-zy(i)) &
               + x2*(1.0_wp+zx(i))*(1.0_wp-zy(i))  &
               + x3*(1.0_wp+zx(i))*(1.0_wp+zy(i))  &
               + x4*(1.0_wp-zx(i))*(1.0_wp+zy(i))) / 4.0_wp
    end do

end subroutine vector_linear_mapping_2D

recursive function undivided_differences(x) result(dif)
    !* Arguments *!
    real(wp), intent(in) :: x(:)

    !* Return values *!
    real(wp) :: dif

    !* Local variables *!
    integer :: n

    n = size(x)

    if (n == 1) then
        dif = x(1)
    else
        dif = undivided_differences(x(2:n)) - undivided_differences(x(1:n-1))
    end if

end function undivided_differences

end module Utilities
