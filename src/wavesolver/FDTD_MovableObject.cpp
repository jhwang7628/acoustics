#include <wavesolver/FDTD_MovableObject.h> 

//##############################################################################
//##############################################################################
void FDTD_MovableObject::
PrintBoundingBox()
{
    std::cout << "minBound = " << _bboxWorld.minBound << std::endl; 
    std::cout << "maxBound = " << _bboxWorld.maxBound << std::endl; 
}

