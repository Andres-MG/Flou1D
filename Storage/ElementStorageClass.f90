#define __ISEXTENDED
#define __BASETYPE   ElemSt_t
#define __VALUETYPE  Elem_t
#define __VECTORTYPE ElemVector_t
#define __LISTTYPE   ElemList_t

module ElementBaseStorageClass

    use ElementClass

#include "GenericStorage/BaseStorageClass.inc"

end module ElementBaseStorageClass

module ElementVectorClass

    use ElementBaseStorageClass
    use ElementClass

#include "GenericStorage/VectorStorageClass.inc"

end module ElementVectorClass

module ElementListClass

    use ElementBaseStorageClass
    use LinkedListClass
    use ElementClass

#include "GenericStorage/ListStorageClass.inc"

end module ElementListClass

module ElementStorageClass

    use ElementBaseStorageClass
    use ElementVectorClass
    use ElementListClass

end module ElementStorageClass
