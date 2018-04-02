#include <wavesolver/FDTD_MovableObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>
#include <utils/Conversions.h>

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
    _modelingTransform.pretranslate(Eigen::Vector3d(x,y,z));
    _modelingTransformInverse.translate(Eigen::Vector3d(-x,-y,-z));
    UpdateBoundingBox();
}

//##############################################################################
// scaling in the object frame
//##############################################################################
void FDTD_MovableObject::
ApplyScale(const REAL scale)
{
    _modelingTransform.prescale(scale);
    _modelingTransformInverse.scale(1./scale);
    UpdateBoundingBox();
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
ClearTransform()
{
    _modelingTransform.setIdentity();
    _modelingTransformInverse.setIdentity();
}

//##############################################################################
//##############################################################################
Vector3d FDTD_MovableObject::
WorldToObjectPoint(const Vector3d &worldPoint)
{
    Eigen::Vector3d position = Conversions::ToEigen(worldPoint);
    position = _modelingTransformInverse*position.eval();
    return Conversions::ToVector3(position);
}

//##############################################################################
//##############################################################################
Vector3d FDTD_MovableObject::
ObjectToWorldPoint(const Vector3d &objectPoint)
{
    Eigen::Vector3d position = Conversions::ToEigen(objectPoint);
    position = _modelingTransform*position.eval();
    return Conversions::ToVector3(position);
}

//##############################################################################
//##############################################################################
Vector3d FDTD_MovableObject::
WorldToObjectVector(const Vector3d &worldVector)
{
    Eigen::Vector3d vector_e = Conversions::ToEigen(worldVector);
    vector_e = _modelingTransformInverse.linear()*vector_e.eval();
    return Conversions::ToVector3(vector_e);
}

//##############################################################################
//##############################################################################
Vector3d FDTD_MovableObject::
ObjectToWorldVector(const Vector3d &objectVector)
{
    Eigen::Vector3d vector_e = Conversions::ToEigen(objectVector);
    vector_e = _modelingTransform.linear()*vector_e.eval();
    return Conversions::ToVector3(vector_e);
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

