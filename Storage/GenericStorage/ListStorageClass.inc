    implicit none

    private

    public :: __LISTTYPE

    type, extends(__BASETYPE) :: __LISTTYPE
        private
        type(List_t), allocatable :: list
    contains
        procedure :: construct  => List_construct
        procedure :: destruct   => List_destruct
        procedure :: size       => List_get_length
        procedure :: max_size   => List_max_size
        procedure :: push       => List_push
        procedure :: push_back  => List_push_back
        procedure :: pop        => List_pop
        procedure :: pop_back   => List_pop_back
        procedure :: insert     => List_insert
        procedure :: update     => List_update
        procedure :: erase      => List_erase
        procedure :: at         => List_at
        procedure :: val_at     => List_val_at
        procedure :: next       => List_next
        procedure :: prev       => List_prev
        procedure :: front      => List_front
        procedure :: back       => List_back
        procedure :: reset_last => List_reset_last_ptr
        final     :: List_destruct_final
    end type __LISTTYPE

contains

subroutine List_construct(this, size)
    !* Arguments *!
    integer,           intent(in) :: size  ! Useless here
    ! Derived types
    class(__LISTTYPE), intent(inout) :: this


    if (allocated(this%list)) deallocate(this%list)
    allocate(this%list)
    call this%list%construct(size)

end subroutine List_construct

function List_get_length(this) result(size)
    !* Arguments *!
    class(__LISTTYPE), intent(in) :: this

    !* Return values *!
    integer :: size


    size = this%list%size()

end function List_get_length

function List_max_size(this) result(size)
    !* Arguments *!
    class(__LISTTYPE), intent(in) :: this

    !* Return values *!
    integer :: size


    ! Unknown limit for a dynamic list
    size = huge(1)

end function List_max_size

subroutine List_push(this, val)
    !* Arguments *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), intent(in) :: val
#else
    __VALUETYPE,       intent(in) :: val
#endif
    ! Derived types
    class(__LISTTYPE), intent(inout) :: this


    call this%list%push(val)

end subroutine List_push

subroutine List_push_back(this, val)
    !* Arguments *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), intent(in) :: val
#else
    __VALUETYPE,       intent(in) :: val
#endif
    ! Derived types
    class(__LISTTYPE), intent(inout) :: this


    call this%list%push_back(val)

end subroutine List_push_back

function List_pop(this) result(val)
    !* Arguments *!
    class(__LISTTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE) :: val
#else
    __VALUETYPE       :: val
#endif

    !* Local variables *!
    class(*), allocatable :: tmp


    tmp = this%list%pop()

    select type (tmp)
    type is (__VALUETYPE)
        val = tmp
    end select

end function List_pop

function List_pop_back(this) result(val)
    !* Arguments *!
    class(__LISTTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE) :: val
#else
    __VALUETYPE       :: val
#endif

    !* Local variables *!
    class(*), allocatable :: tmp


    tmp = this%list%pop_back()

    select type (tmp)
    type is (__VALUETYPE)
        val = tmp
    end select

end function List_pop_back

subroutine List_insert(this, val, ind)
    !* Arguments *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), intent(in) :: val
#else
    __VALUETYPE,       intent(in) :: val
#endif
    integer,           intent(in) :: ind
    ! Derived types
    class(__LISTTYPE), intent(inout) :: this


    call this%list%insert(val, ind)

end subroutine List_insert

subroutine List_update(this, val, ind)
    !* Arguments *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), intent(in) :: val
#else
    __VALUETYPE,       intent(in) :: val
#endif
    integer,           intent(in) :: ind
    ! Derived types
    class(__LISTTYPE), intent(inout) :: this


    call this%list%update(val, ind)

end subroutine List_update

subroutine List_erase(this, ind)
    !* Arguments *!
    integer, intent(in) :: ind
    ! Derived types
    class(__LISTTYPE), intent(inout) :: this


    call this%list%erase(ind)

end subroutine List_erase

function List_at(this, ind) result(val)
    !* Arguments *!
    integer, intent(in) :: ind
    ! Derived types
    class(__LISTTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), pointer :: val
#else
    __VALUETYPE,       pointer :: val
#endif

    !* Local variables *!
    class(*), pointer :: tmp


    tmp => this%list%at(ind)

    select type (tmp)
    type is (__VALUETYPE)
        val => tmp
    end select

    nullify(tmp)

end function List_at

function List_val_at(this, ind) result(val)
    !* Arguments *!
    integer, intent(in) :: ind
    ! Derived types
    class(__LISTTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE) :: val
#else
    __VALUETYPE       :: val
#endif

    !* Local variables *!
    class(*), pointer :: tmp


    tmp => this%list%at(ind)

    select type (tmp)
    type is (__VALUETYPE)
        val = tmp
    end select

    nullify(tmp)

end function List_val_at

function List_next(this) result(val)
    !* Arguments *!
    class(__LISTTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), pointer :: val
#else
    __VALUETYPE,       pointer :: val
#endif

    !* Local variables *!
    class(*), pointer :: tmp


    tmp => this%list%next()

    select type (tmp)
    type is (__VALUETYPE)
        val => tmp
    end select

    nullify(tmp)

end function List_next

function List_prev(this) result(val)
    !* Arguments *!
    class(__LISTTYPE), intent(inout) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), pointer :: val
#else
    __VALUETYPE,       pointer :: val
#endif

    !* Local variables *!
    class(*), pointer :: tmp


    tmp => this%list%prev()

    select type (tmp)
    type is (__VALUETYPE)
        val => tmp
    end select

    nullify(tmp)

end function List_prev

function List_front(this) result(val)
    !* Arguments *!
    class(__LISTTYPE), intent(in) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), pointer :: val
#else
    __VALUETYPE,       pointer :: val
#endif

    !* Local variables *!
    class(*), pointer :: tmp


    tmp => this%list%front()

    select type (tmp)
    type is (__VALUETYPE)
        val => tmp
    end select

    nullify(tmp)

end function List_front

function List_back(this) result(val)
    !* Arguments *!
    class(__LISTTYPE), intent(in) :: this

    !* Return values *!
#ifdef __ISEXTENDED
    type(__VALUETYPE), pointer :: val
#else
    __VALUETYPE,       pointer :: val
#endif

    !* Local variables *!
    class(*), pointer :: tmp


    tmp => this%list%back()

    select type (tmp)
    type is (__VALUETYPE)
        val => tmp
    end select

    nullify(tmp)

end function List_back

subroutine List_reset_last_ptr(this, ind)
    !* Arguments *!
    integer, intent(in) :: ind
    ! Derived types
    class(__LISTTYPE), intent(inout) :: this


    call this%list%reset_last(ind)

end subroutine List_reset_last_ptr

subroutine List_destruct_final(this)
    !* Arguments *!
    type(__LISTTYPE), intent(inout) :: this


    if (allocated(this%list)) deallocate(this%list)

end subroutine List_destruct_final

subroutine List_destruct(this)
    !* Arguments *!
    class(__LISTTYPE), intent(inout) :: this


    select type (this)
    type is (__LISTTYPE)
        call List_destruct_final(this)
    end select

end subroutine List_destruct
