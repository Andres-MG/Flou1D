    implicit none

    private

    public :: __BASETYPE

    type, abstract :: __BASETYPE
    contains
        procedure(Construct_Int),  deferred :: construct
        procedure(Construct_Int),  deferred :: reserve
        procedure(Destruct_Int),   deferred :: destruct
        procedure(GetSize_Int),    deferred :: size
        procedure(GetSize_Int),    deferred :: max_size
        procedure(Push_Int),       deferred :: push
        procedure(Push_Int),       deferred :: push_back
        procedure(Pop_Int),        deferred :: pop
        procedure(Pop_Int),        deferred :: pop_back
        procedure(Insert_Int),     deferred :: insert
        procedure(Insert_Int),     deferred :: update
        procedure(Erase_Int),      deferred :: erase
        procedure(GetAt_Int),      deferred :: at
        procedure(GetValAt_Int),   deferred :: val_at
        procedure(GetPointer_Int), deferred :: next
        procedure(GetPointer_Int), deferred :: prev
        procedure(GetValue_Int),   deferred :: front
        procedure(GetValue_Int),   deferred :: back
        procedure(Reset_Int),      deferred :: reset_last
    end type __BASETYPE

    abstract interface
        subroutine Construct_Int(this, size)
            import __BASETYPE
            class(__BASETYPE), intent(inout) :: this
            integer,           intent(in)    :: size
        end subroutine Construct_Int

        function GetSize_Int(this)
            import __BASETYPE
            class(__BASETYPE), intent(in) :: this
            integer                       :: GetSize_Int
        end function GetSize_Int

        subroutine Push_Int(this, val)
            import __BASETYPE
#ifdef __ISEXTENDED
            import __VALUETYPE
#endif
            class(__BASETYPE), intent(inout) :: this
#ifdef __ISEXTENDED
            type(__VALUETYPE), intent(in)    :: val
#else
            __VALUETYPE,       intent(in)    :: val
#endif
        end subroutine Push_Int

        function Pop_Int(this)
            import __BASETYPE
#ifdef __ISEXTENDED
            import __VALUETYPE
#endif
            class(__BASETYPE), intent(inout) :: this
#ifdef __ISEXTENDED
            type(__VALUETYPE)                :: Pop_Int
#else
            __VALUETYPE                      :: Pop_Int
#endif
        end function Pop_Int

        function GetPointer_Int(this)
            import __BASETYPE
#ifdef __ISEXTENDED
            import __VALUETYPE
#endif
            class(__BASETYPE), intent(inout) :: this
#ifdef __ISEXTENDED
            type(__VALUETYPE), pointer       :: GetPointer_Int
#else
            __VALUETYPE,       pointer       :: GetPointer_Int
#endif
        end function GetPointer_Int

        function GetValue_Int(this)
            import __BASETYPE
#ifdef __ISEXTENDED
            import __VALUETYPE
#endif
            class(__BASETYPE), intent(in) :: this
#ifdef __ISEXTENDED
            type(__VALUETYPE), pointer    :: GetValue_Int
#else
            __VALUETYPE,       pointer    :: GetValue_Int
#endif
        end function GetValue_Int

        subroutine Insert_Int(this, val, ind)
            import __BASETYPE
#ifdef __ISEXTENDED
            import __VALUETYPE
#endif
            class(__BASETYPE), intent(inout) :: this
#ifdef __ISEXTENDED
            type(__VALUETYPE), intent(in)    :: val
#else
            __VALUETYPE,       intent(in)    :: val
#endif
            integer,           intent(in)    :: ind
        end subroutine Insert_Int

        subroutine Erase_Int(this, ind)
            import __BASETYPE
            class(__BASETYPE), intent(inout) :: this
            integer,           intent(in)    :: ind
        end subroutine Erase_Int

        function GetAt_Int(this, ind)
            import __BASETYPE
#ifdef __ISEXTENDED
            import __VALUETYPE
#endif
            class(__BASETYPE), intent(inout) :: this
            integer,           intent(in)    :: ind
#ifdef __ISEXTENDED
            type(__VALUETYPE), pointer       :: GetAt_Int
#else
            __VALUETYPE,       pointer       :: GetAt_Int
#endif
        end function GetAt_Int

        function GetValAt_Int(this, ind)
            import __BASETYPE
#ifdef __ISEXTENDED
            import __VALUETYPE
#endif
            class(__BASETYPE), intent(inout) :: this
            integer,           intent(in)    :: ind
#ifdef __ISEXTENDED
            type(__VALUETYPE)                :: GetValAt_Int
#else
            __VALUETYPE                      :: GetValAt_Int
#endif
        end function GetValAt_Int

        subroutine Reset_Int(this, ind)
            import __BASETYPE
            class(__BASETYPE), intent(inout) :: this
            integer,           intent(in)    :: ind
        end subroutine Reset_Int

        subroutine Destruct_Int(this)
            import __BASETYPE
            class(__BASETYPE), intent(inout) :: this
        end subroutine Destruct_Int
    end interface
