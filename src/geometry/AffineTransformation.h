#ifndef AFFINETRANSFORMATION_H 
#define AFFINETRANSFORMATION_H 

#include <linearalgebra/Vector3.hpp> 

namespace Geometry
{

// reflect the input point with respect to a plane defined by point+normal
template <typename T> 
void Reflection(const Vector3<T> &originalPoint, const Vector3<T> &normal, const Vector3<T> &pointWithNormal, Vector3<T> &reflectedPoint)
{
    reflectedPoint = originalPoint - normal * (2.0*normal.dotProduct(originalPoint-pointWithNormal)); 
}


}

#endif
