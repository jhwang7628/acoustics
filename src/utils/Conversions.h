#ifndef CONVERSIONS_H 
#define CONVERSIONS_H 

#include <linearalgebra/Vector3.hpp> 
#include <Eigen/Dense> 

namespace Conversions
{
    template<typename T> 
    Eigen::Matrix<T,3,1> ToEigen(const Vector3<T> &vec)
    {
        return Eigen::Matrix<T,3,1>(vec.x,vec.y,vec.z); 
    }

    template<typename T> 
    Vector3<T> ToVector3(const Eigen::Matrix<T,3,1> &vec)
    {
        return Vector3<T>(vec[0],vec[1],vec[2]); 
    }
};

#endif
