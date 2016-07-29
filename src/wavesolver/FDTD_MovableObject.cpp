#include <wavesolver/FDTD_MovableObject.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>

//##############################################################################
//##############################################################################
Eigen::Vector3d FDTD_MovableObject::
GetTranslation() const 
{
    return _modelingTransform.translation(); 
}

//##############################################################################
//##############################################################################
void FDTD_MovableObject::
GetRotationDegree(REAL &angle, Eigen::Vector3d &axis) const 
{
    const Eigen::AngleAxisd angleAxis(_modelingTransform.rotation()); 
    angle = RAD_2_DEG(angleAxis.angle()); 
    axis = angleAxis.axis(); 
}

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
    throw std::runtime_error("**ERROR** not implemented yet"); 
}

//##############################################################################
//##############################################################################
void FDTD_MovableObject::
SetTransform(const double &x, const double &y, const double &z, const double &angle, const double &rotationVectorx, const double &rotationVectory, const double &rotationVectorz)
{
    // reset
    _modelingTransform.setIdentity(); 
    _modelingTransformInverse.setIdentity(); 

    // apply translation
    _modelingTransform.translate(Eigen::Vector3d(x, y, z)); 

    // apply rotation
    Eigen::AngleAxisd rotation(angle, Eigen::Vector3d(rotationVectorx, rotationVectory, rotationVectorz));
    _modelingTransform = _modelingTransform * rotation; // first apply translation and then rotation

    // compute inverse and update bounding box
    _modelingTransformInverse = _modelingTransform.inverse(); 
    UpdateBoundingBox(); 
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

