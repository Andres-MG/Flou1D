#define __ISEXTENDED
#define __SIZE_INC   10
#define __BASETYPE   FaceSt_t
#define __VALUETYPE  Face_t
#define __VECTORTYPE FaceVector_t
#define __LISTTYPE   FaceList_t

module FaceBaseStorageClass

    use ElementClass

#include "GenericStorage/BaseStorageClass.inc"

end module FaceBaseStorageClass

module FaceVectorClass

    use FaceBaseStorageClass
    use ElementClass

#include "GenericStorage/VectorStorageClass.inc"

end module FaceVectorClass

module FaceListClass

    use FaceBaseStorageClass
    use LinkedListClass
    use ElementClass

#include "GenericStorage/ListStorageClass.inc"

end module FaceListClass

module FaceStorageClass

    use FaceBaseStorageClass
    use FaceListClass
    use FaceVectorClass

end module FaceStorageClass
