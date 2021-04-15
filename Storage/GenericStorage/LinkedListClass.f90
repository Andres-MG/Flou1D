module LinkedListClass

    use Constants, only: wp

    implicit none

    private

    public :: Node_t
    public :: List_t

    type :: Node_t
        private
        class(*),      pointer :: data      => null()
        class(Node_t), pointer :: next_node => null()
        class(Node_t), pointer :: prev_node => null()
    contains
        procedure :: set        => Node_set_value
        procedure :: value      => Node_value
        procedure :: prev       => Node_prev_node
        procedure :: next       => Node_next_node
        procedure :: set_prev   => Node_set_prev_node
        procedure :: set_next   => Node_set_next_node
        procedure :: unset_next => Node_unset_next_node
        procedure :: unset_prev => Node_unset_prev_node
        final     :: Node_destruct
    end type Node_t

    type :: List_t
        private
        integer               :: length     = 0
        integer               :: halfLength = 0
        type(Node_t), pointer :: head => null()
        type(Node_t), pointer :: tail => null()
        type(Node_t), pointer :: last => null()
    contains
        procedure :: construct  => List_construct
        procedure :: destruct   => List_destruct
        procedure :: size       => List_get_length
        procedure :: push       => List_push
        procedure :: push_back  => List_push_back
        procedure :: pop        => List_pop
        procedure :: pop_back   => List_pop_back
        procedure :: insert     => List_insert
        procedure :: update     => List_update
        procedure :: erase      => List_erase
        procedure :: at         => List_at
        procedure :: next       => List_next
        procedure :: prev       => List_prev
        procedure :: front      => List_front
        procedure :: back       => List_back
        procedure :: reset_last => List_reset_last_ptr
        final     :: List_destruct_final
    end type List_t

contains

subroutine Node_set_value(this, val)
    !* Arguments *!
    class(*), intent(in) :: val
    ! Derived types
    class(Node_t), intent(inout) :: this


    if (associated(this%data)) deallocate(this%data)
    allocate(this%data, source=val)

end subroutine Node_set_value

function Node_value(this) result(val)
    !* Arguments *!
    class(Node_t), intent(in) :: this

    !* Return values *!
    class(*), pointer :: val


    val => this%data

end function Node_value

function Node_prev_node(this) result(node)
    !* Arguments *!
    class(Node_t), intent(in) :: this

    !* Return values *!
    class(Node_t), pointer :: node


    node => this%prev_node

end function Node_prev_node

function Node_next_node(this) result(node)
    !* Arguments *!
    class(Node_t), intent(in) :: this

    !* Return values *!
    class(Node_t), pointer :: node


    node => this%next_node

end function Node_next_node

subroutine Node_set_prev_node(this, node)
    !* Arguments *!
    class(Node_t),         intent(inout) :: this
    class(Node_t), target, intent(in)    :: node


    this%prev_node => node

end subroutine Node_set_prev_node

subroutine Node_set_next_node(this, node)
    !* Arguments *!
    class(Node_t),         intent(inout) :: this
    class(Node_t), target, intent(in)    :: node


    this%next_node => node

end subroutine Node_set_next_node

subroutine Node_unset_prev_node(this)
    !* Arguments *!
    class(Node_t), intent(inout) :: this


    nullify(this%prev_node)

end subroutine Node_unset_prev_node

subroutine Node_unset_next_node(this)
    !* Arguments *!
    class(Node_t), intent(inout) :: this


    nullify(this%next_node)

end subroutine Node_unset_next_node

subroutine Node_destruct(this)
    !* Arguments *!
    type(Node_t), intent(inout) :: this


    call this%unset_next()
    call this%unset_prev()
    if (associated(this%data)) deallocate(this%data)

end subroutine Node_destruct


subroutine List_construct(this, size)
    !* Arguments *!
    integer,  intent(in) :: size  ! Useless here
    ! Derived types
    class(List_t), intent(inout) :: this


    ! Check that everything is clear
    if (associated(this%head)) call this%destruct()

    allocate(this%head)
    this%tail => this%head
    this%length = 0
    this%halfLength = 0

end subroutine List_construct

function List_get_length(this) result(size)
    !* Arguments *!
    class(List_t), intent(in) :: this

    !* Return values *!
    integer :: size


    size = this%length

end function List_get_length

subroutine List_push(this, val)
    !* Arguments *!
    class(*), intent(in) :: val
    ! Derived types
    class(List_t), intent(inout) :: this

    !* Local variables *!
    type(Node_t), pointer :: tmp


    ! Handle the case of a constructed list without nodes
    if (this%length == 0 .and. associated(this%head)) then

        call this%head%set(val)
        this%length     = 1
        this%halfLength = 1

    else

        ! Create the node and connect it to the list
        allocate(tmp)
        call tmp%set(val)

        call tmp%set_next(this%head)
        call this%head%set_prev(tmp)
        this%length = this%length + 1
        this%halfLength = ceiling(0.5_wp*this%length)

        ! Reset the head
        this%head => tmp

        nullify(tmp)

    end if

end subroutine List_push

subroutine List_push_back(this, val)
    !* Arguments *!
    class(*), intent(in) :: val
    ! Derived types
    class(List_t), intent(inout) :: this

    !* Local variables *!
    type(Node_t), pointer :: tmp


    ! Handle the case of a constructed list without nodes
    if (this%length == 0 .and. associated(this%head)) then

        call this%head%set(val)
        this%length     = 1
        this%halfLength = 1

    else

        ! Create the node and connect it to the list
        allocate(tmp)
        call tmp%set(val)

        call tmp%set_prev(this%tail)
        call this%tail%set_next(tmp)
        this%length = this%length + 1
        this%halfLength = ceiling(0.5_wp*this%length)

        ! Reset the tail
        this%tail => tmp

        nullify(tmp)

    end if

end subroutine List_push_back

function List_pop(this) result(val)
    !* Arguments *!
    class(List_t), intent(inout) :: this

    !* Return values *!
    class(*), allocatable :: val

    !* Local variables *!
    type(Node_t), pointer :: tmp


    ! Get the value and delete the first element
    val = this%head%value()

    tmp => this%head%next()
    deallocate(this%head)
    this%head => tmp

    ! Handle possible empty list
    if (associated(tmp)) then
        call tmp%unset_prev()
    else
        nullify(this%tail)
        nullify(this%last)
    end if

    this%length = this%length - 1
    this%halfLength = ceiling(0.5_wp*this%length)

    ! Clean up
    nullify(tmp)

end function List_pop

function List_pop_back(this) result(val)
    !* Arguments *!
    class(List_t), intent(inout) :: this

    !* Return values *!
    class(*), allocatable :: val

    !* Local variables *!
    type(Node_t), pointer :: tmp


    ! Get the value and delete the last element
    val = this%tail%value()

    tmp => this%tail%prev()
    deallocate(this%tail)
    this%tail => tmp

    ! Handle possible empty list
    if (associated(tmp)) then
        call tmp%unset_next()
    else
        nullify(this%head)
        nullify(this%last)
    end if

    this%length = this%length - 1
    this%halfLength = ceiling(0.5_wp*this%length)

    ! Clean up
    nullify(tmp)

end function List_pop_back

subroutine List_insert(this, val, valInd)
    !* Arguments *!
    class(*), intent(in) :: val
    integer,  intent(in) :: valInd
    ! Derived types
    class(List_t), intent(inout) :: this

    !* Local variables *!
    integer :: i
    integer :: ind
    type(Node_t), pointer :: newNode
    type(Node_t), pointer :: tmp


    ! Check bounds
    if (valInd < 1) then
        ind = 1
    else if (valInd > this%length+1) then
        ind = this%length + 1
    else
        ind = valInd
    end if


    ! Handle the case of a constructed list without nodes
    if (ind == 1 .and. this%length == 0 .and. associated(this%head)) then

        call this%head%set(val)
        this%length     = this%length + 1
        this%halfLength = ceiling(0.5_wp*this%length)

    else

        ! Create the new node
        allocate(newNode)
        call newNode%set(val)
        this%length = this%length + 1
        this%halfLength = ceiling(0.5_wp*this%length)

        ! If 'ind=1' -> change the head
        if (ind == 1) then

            call newNode%set_next(this%head)
            call this%head%set_prev(newNode)
            this%head => newNode

        else

            ! Loop until the previous node otherwise
            tmp => this%head
            do i = 2, ind-1
                tmp => tmp%next()
            end do

            ! Link the node to the list
            call newNode%set_prev(tmp)
            call newNode%set_next(tmp%next())

            call tmp%set_next(newNode)

            ! Be careful with the tail node
            if (ind == this%length) then
                this%tail => newNode
            else
                tmp => newNode%next()
                call tmp%set_prev(newNode)
            end if

            nullify(tmp)

        end if

        ! Clean up
        nullify(newNode)

    end if

end subroutine List_insert

subroutine List_update(this, val, valInd)
    !* Arguments *!
    class(*), intent(in) :: val
    integer,  intent(in) :: valInd
    ! Derived types
    class(List_t), intent(inout) :: this

    !* Local variables *!
    integer :: i
    integer :: ind
    type(Node_t), pointer :: tmp


    ! Check bounds
    if (valInd < 1) then
        ind = 1
    else if (valInd > this%length+1) then
        ind = this%length + 1
    else
        ind = valInd
    end if

    ! Loop until the selected node
    tmp => this%head
    do i = 2, ind
        tmp => tmp%next()
    end do

    call tmp%set(val)

    ! Clean up
    nullify(tmp)

end subroutine List_update

subroutine List_erase(this, nodeInd)
    !* Arguments *!
    integer, intent(in) :: nodeInd
    ! Derived types
    class(List_t), intent(inout) :: this

    !* Local variables *!
    integer :: i
    integer :: ind
    type(Node_t), pointer :: tmp
    type(Node_t), pointer :: tmp2


    ! Check bounds
    if (nodeInd < 1) then
        ind = 1
    else if (nodeInd > this%length) then
        ind = this%length
    else
        ind = nodeInd
    end if

    ! Loop until the node to be erased
    tmp => this%head
    do i = 2, ind
        tmp => tmp%next()
    end do

    if (associated(tmp, this%head) .and. associated(tmp, this%tail)) then

        ! Clean everything if we erased the last node
        nullify(this%head)
        nullify(this%tail)
        nullify(this%last)

    else if (associated(tmp, this%head)) then

        this%head => this%head%next()
        call this%head%unset_prev()

    else if (associated(tmp, this%tail)) then

        this%tail => this%tail%prev()
        call this%tail%unset_next()

    else

        tmp2 => tmp%prev()
        call tmp2%set_next(tmp%next())

        tmp2 => tmp%next()
        call tmp2%set_prev(tmp%prev())

        nullify(tmp2)

    end if

    ! Remove the node
    deallocate(tmp)
    this%length = this%length - 1
    this%halfLength = ceiling(0.5_wp*this%length)

end subroutine List_erase

function List_at(this, nodeInd) result(val)
    !* Arguments *!
    integer, intent(in) :: nodeInd
    ! Derived types
    class(List_t), intent(inout) :: this

    !* Return values *!
    class(*), pointer :: val

    !* Local variables *!
    integer :: i
    integer :: ind
    type(Node_t), pointer :: tmp


    ! Check bounds
    if (nodeInd < 1) then
        ind = 1
    else if (nodeInd > this%length) then
        ind = this%length
    else
        ind = nodeInd
    end if

    ! Find the shortest way to get there
    if (ind - this%halfLength <= 0) then

        tmp => this%head
        do i = 2, ind
            tmp => tmp%next()
        end do

    else

        ind = this%length - (ind-1)
        tmp => this%tail
        do i = 2, ind
            tmp => tmp%prev()
        end do

    end if

    ! Set the 'last' pointer
    this%last => tmp

    ! Return value and clean up
    val => tmp%value()
    nullify(tmp)

end function List_at

function List_next(this) result(val)
    !* Arguments *!
    class(List_t), intent(inout) :: this

    !* Return values *!
    class(*), pointer :: val


    ! Point to the next node. rolling over if necessary
    this%last => this%last%next()
    if (.not. associated(this%last)) then
        this%last => this%head
    end if

    val => this%last%value()

end function List_next

function List_prev(this) result(val)
    !* Arguments *!
    class(List_t), intent(inout) :: this

    !* Return values *!
    class(*), pointer :: val


    ! Point to the previous node. rolling over if necessary
    this%last => this%last%prev()
    if (.not. associated(this%last)) then
        this%last => this%tail
    end if

    val => this%last%value()

end function List_prev

function List_front(this) result(val)
    !* Arguments *!
    class(List_t), intent(in) :: this

    !* Return values *!
    class(*), pointer :: val


    val => this%head%value()

end function List_front

function List_back(this) result(val)
    !* Arguments *!
    class(List_t), intent(in) :: this

    !* Return values *!
    class(*), pointer :: val


    val => this%tail%value()

end function List_back

subroutine List_reset_last_ptr(this, ind)
    !* Arguments *!
    integer, intent(in) :: ind
    ! Derived types
    class(List_t), intent(inout) :: this

    !* Local variables *!
    class(*), pointer :: tmp


    ! Ind = -1 might be useful for backward iterations
    if (ind == 1 .or. ind == -1) then
        this%last => this%head

    ! Ind = 0 might be useful for forward iterations
    else if (ind == this%length .or. ind == 0) then
        this%last => this%tail

    ! Not the best way..., but 'List_at' already sets the 'last' pointer
    else
        tmp => this%at(ind)
        nullify(tmp)

    end if

end subroutine List_reset_last_ptr

subroutine List_destruct_final(this)
    !* Arguments *!
    type(List_t), intent(inout) :: this

    !* Local variables *!
    type(Node_t), pointer :: tmp1
    type(Node_t), pointer :: tmp2


    ! Point to the tail
    tmp1 => this%tail

    ! And move backwards to the head
    do while (associated(tmp1))

        tmp2 => tmp1%prev()
        deallocate(tmp1)
        tmp1 => tmp2

    end do

    ! Final clean up
    this%length = 0
    this%halfLength = 0

    nullify(this%head)
    nullify(this%tail)
    nullify(this%last)

    nullify(tmp1)
    nullify(tmp2)

end subroutine List_destruct_final

subroutine List_destruct(this)
    !* Arguments *!
    class(List_t), intent(inout) :: this


    select type (list => this)
    type is (List_t)
        call List_destruct_final(list)
    end select

end subroutine List_destruct

end module LinkedListClass
