    implicit none

    private

    public :: __VECTORTYPE

    type, extends(__BASETYPE) :: __VECTORTYPE
        private
        integer :: length
        integer :: tailInd
        integer :: lastInd
#ifdef __ISEXTENDED
        type(__VALUETYPE), pointer :: vector(:)
#else
        __VALUETYPE,       pointer :: vector(:)
#endif
    contains
        procedure :: construct  => Vector_construct
        procedure :: reserve    => Vector_reserve
        procedure :: destruct   => Vector_destruct
        procedure :: size       => Vector_get_length
        procedure :: max_size   => Vector_max_size
        procedure :: push       => Vector_push
        procedure :: push_back  => Vector_push_back
        procedure :: pop        => Vector_pop
        procedure :: pop_back   => Vector_pop_back
        procedure :: insert     => Vector_insert
        procedure :: update     => Vector_update
        procedure :: erase      => Vector_erase
        procedure :: at         => Vector_at
        procedure :: val_at     => Vector_val_at
        procedure :: next       => Vector_next
        procedure :: prev       => Vector_prev
        procedure :: front      => Vector_front
        procedure :: back       => Vector_back
        procedure :: reset_last => Vector_reset_last_ptr
        final     :: Vector_destruct_final
    end type __VECTORTYPE

contains

subroutine Vector_construct(this, size)
    !* Arguments *!
    integer, intent(in) :: size
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this


    allocate(this%vector(size))

    this%tailInd = 0
    this%length  = size

end subroutine Vector_construct

subroutine Vector_reserve(this, size)
    !* Arguments *!
    integer, intent(in) :: size
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this

    !* Local variables *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), allocatable :: tmp(:)
#else
    __VALUETYPE,       allocatable :: tmp(:)
#endif


    ! Construct if not associated, but this should not be used
    if (.not. associated(this%vector)) then

        allocate(this%vector(size))
        this%tailInd = 0
        this%length  = size

    ! Reallocate the private vector
    else if (this%length < size) then

        allocate(tmp(size))
        tmp(:this%length) = this%vector

        deallocate(this%vector)   ! Cannot use move_alloc with pointers :(
        allocate(this%vector, source=tmp)

        deallocate(tmp)
        this%length = size

    end if

end subroutine Vector_reserve

function Vector_get_length(this) result(size)
    !* Arguments *!
    class(__VECTORTYPE), intent(in) :: this

    !* Return values *!
    integer :: size


    size = this%tailInd

end function Vector_get_length

function Vector_max_size(this) result(size)
    !* Arguments *!
    class(__VECTORTYPE), intent(in) :: this

    !* Return values *!
    integer :: size


    size = this%length

end function Vector_max_size

subroutine Vector_push(this, val)
    !* Arguments *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), intent(in) :: val
#else
    __VALUETYPE,       intent(in) :: val
#endif
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this


    ! Reallocate if needed
    if (this%tailInd+1 > this%length) then
        call this%reserve(this%length + __SIZE_INC)
    end if

    ! Shift all the vector to the 'right' and insert 'val' in the first position
    this%vector(2:this%tailInd+1) = this%vector(:this%tailInd)
    this%vector(1) = val
    this%tailInd = this%tailInd + 1

end subroutine Vector_push

subroutine Vector_push_back(this, val)
    !* Arguments *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), intent(in) :: val
#else
    __VALUETYPE,       intent(in) :: val
#endif
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this


    ! Reallocate if needed
    if (this%tailInd+1 > this%length) then
        call this%reserve(this%length + __SIZE_INC)
    end if

    ! Insert 'val' in the last position
    this%tailInd = this%tailInd + 1
    this%vector(this%tailInd) = val

end subroutine Vector_push_back

function Vector_pop(this) result(val)
    !* Arguments *!
    class(__VECTORTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE) :: val
#else
    __VALUETYPE       :: val
#endif


    val = this%vector(1)
    this%tailInd = this%tailInd - 1
    this%vector(:this%tailInd) = this%vector(2:this%tailInd+1)

end function Vector_pop

function Vector_pop_back(this) result(val)
    !* Arguments *!
    class(__VECTORTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE) :: val
#else
    __VALUETYPE       :: val
#endif


    val = this%vector(this%tailInd)
    this%tailInd = this%tailInd - 1

end function Vector_pop_back

subroutine Vector_insert(this, val, ind)
    !* Arguments *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), intent(in) :: val
#else
    __VALUETYPE,       intent(in) :: val
#endif
    integer,           intent(in) :: ind
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this


    ! Reallocate if necessary
    if (this%tailInd+1 > this%length) then
        call this%reserve(this%length + __SIZE_INC)
    end if

    ! Shift the vector and insert 'val'
    this%vector(ind+1:this%tailInd+1) = this%vector(ind:this%tailInd)
    this%tailInd = this%tailInd + 1
    this%vector(ind) = val

end subroutine Vector_insert

subroutine Vector_update(this, val, ind)
    !* Arguments *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), intent(in) :: val
#else
    __VALUETYPE,       intent(in) :: val
#endif
    integer,           intent(in) :: ind
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this


    this%vector(ind) = val

end subroutine Vector_update

subroutine Vector_erase(this, ind)
    !* Arguments *!
    integer, intent(in) :: ind
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this


    ! Shift the vector to the left
    this%vector(ind:this%tailInd-1) = this%vector(ind+1:this%tailInd)
    this%tailInd = this%tailInd - 1

end subroutine Vector_erase

function Vector_at(this, ind) result(val)
    !* Arguments *!
    integer, intent(in) :: ind
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), pointer :: val
#else
    __VALUETYPE,       pointer :: val
#endif


    val => this%vector(ind)
    this%lastInd = ind

end function Vector_at

function Vector_val_at(this, ind) result(val)
    !* Arguments *!
    integer, intent(in) :: ind
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE) :: val
#else
    __VALUETYPE       :: val
#endif


    val = this%vector(ind)
    this%lastInd = ind

end function Vector_val_at

function Vector_next(this) result(val)
    !* Arguments *!
    class(__VECTORTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), pointer :: val
#else
    __VALUETYPE,       pointer :: val
#endif


    ! Find the next element of the array
    this%lastInd = this%lastInd + 1
    if (this%lastInd > this%tailInd) then
        this%lastInd = 1
    end if

    val => this%vector(this%lastInd)

end function Vector_next

function Vector_prev(this) result(val)
    !* Arguments *!
    class(__VECTORTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), pointer :: val
#else
    __VALUETYPE,       pointer :: val
#endif


    ! Find the previous element of the array
    this%lastInd = this%lastInd - 1
    if (this%lastInd < 1) then
        this%lastInd = this%tailInd
    end if

    val => this%vector(this%lastInd)

end function Vector_prev

function Vector_front(this) result(val)
    !* Arguments *!
    class(__VECTORTYPE), intent(in) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), pointer :: val
#else
    __VALUETYPE,       pointer :: val
#endif


    val => this%vector(1)

end function Vector_front

function Vector_back(this) result(val)
    !* Arguments *!
    class(__VECTORTYPE), intent(in) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), pointer :: val
#else
    __VALUETYPE,       pointer :: val
#endif


    val => this%vector(this%tailInd)

end function Vector_back

subroutine Vector_reset_last_ptr(this, ind)
    !* Arguments *!
    integer, intent(in) :: ind
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this


    ! Ind = -1 might be useful for backward iterations
    if (ind == -1) then
        this%lastInd = 1

    ! Ind = 0 might be useful for forward iterations
    else if (ind == 0) then
        this%lastInd = this%tailInd

    else
        this%lastInd = ind

    end if

end subroutine Vector_reset_last_ptr

subroutine Vector_destruct_final(this)
    !* Arguments *!
    type(__VECTORTYPE), intent(inout) :: this


    if (associated(this%vector)) deallocate(this%vector)
    this%length  = 0
    this%tailInd = 0
    this%lastInd = 0

end subroutine Vector_destruct_final

subroutine Vector_destruct(this)
    !* Arguments *!
    class(__VECTORTYPE), intent(inout) :: this


    select type (this)
    type is (__VECTORTYPE)
        call Vector_destruct_final(this)
    end select

end subroutine Vector_destruct
