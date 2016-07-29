#include <wavesolver/FDTD_MovableObject.h> 

//##############################################################################
//##############################################################################
void FDTD_MovableObject::
ApplyTranslation(const double &x, const double &y, const double &z)
{
    _modelingTransform.translate(Eigen::Vector3d(x,y,z)); 
    _modelingTransformInverse.translate(Eigen::Vector3d(-x,-y,-z)); 
    UpdateBoundingBox(); 
}

//##############################################################################
//##############################################################################
void FDTD_MovableObject::
ApplyRotation(const Quaternion<REAL> &quarternion)
{
}

//##############################################################################
//##############################################################################
void FDTD_MovableObject::
PrintBoundingBox()
{
    std::cout << "minBound = " << _bboxWorld.minBound << std::endl; 
    std::cout << "maxBound = " << _bboxWorld.maxBound << std::endl; 
}

//##############################################################################
//##############################################################################
void FDTD_MovableObject::
PrintTransformation()
{
    std::cout << "modeling matrix = \n"; 
    std::cout << _modelingTransform.matrix() << std::endl;
    std::cout << "inverse modeling matrix = \n"; 
    std::cout << _modelingTransformInverse.matrix() << std::endl;
}

