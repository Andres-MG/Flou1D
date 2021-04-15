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
    integer,           intent(in) :: size  ! Useless here
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this


    ! Check size
    if (size < 1) error stop "Size of __VECTORTYPE must be positive."

    if (associated(this%vector)) deallocate(this%vector)
    allocate(this%vector(size))

    this%tailInd   = 0
    this%length    = size

end subroutine Vector_construct

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


    error stop "__VECTORTYPE does not implement the push() method."

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


    ! Add one to the tail index
    this%tailInd = this%tailInd + 1

    if (this%tailInd > this%length) then
        this%tailInd = this%length
        error stop "The storage vector is already full."

    else
        this%vector(this%tailInd) = val

    end if

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


    error stop "__VECTORTYPE does not implement the pop() method."

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


    if (this%tailInd < 1) then
        this%tailInd = 0
        error stop "The storage vector is already empty."

    else
        val = this%vector(this%tailInd)
        this%tailInd = this%tailInd - 1

    end if

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


    if (ind > this%tailInd) this%tailInd = ind
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


    if (ind > this%tailInd) this%tailInd = ind
    this%vector(ind) = val

end subroutine Vector_update

subroutine Vector_erase(this, ind)
    !* Arguments *!
    integer, intent(in) :: ind
    ! Derived types
    class(__VECTORTYPE), intent(inout) :: this


    if (ind == this%tailInd) this%tailInd = this%tailInd - 1

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
