#ifndef FDTD_PLANE_CONSTRAINT_H
#define FDTD_PLANE_CONSTRAINT_H
#include <memory>
#include "linearalgebra/Vector3.hpp"
#include "config.h"
//##############################################################################
// Class FDTD_PlaneContraint
//##############################################################################
class FDTD_PlaneConstraint
{
    private:
        int  _normalDirection = 1;
        int  _normalSign      = 1; 
        REAL _height; 
    public: 
        FDTD_PlaneConstraint(const REAL &h=(REAL)0.0)
            : _height(h)
        {}
        FDTD_PlaneConstraint(const int &dir, const int &sign, 
                             const REAL &h=(REAL)0.0)
            : _normalDirection(dir),
              _normalSign(sign),
              _height(h)
        {}
        inline Vector3d Normal() const
        {Vector3d nml(0,0,0); nml[_normalDirection]=_normalSign; return nml;}
        inline REAL SignedDistanceToPlane(const Vector3d &pos) const
        {return (pos[_normalDirection] - _height)*_normalSign;}
        inline REAL UnsignedDistanceToPlane(const Vector3d &pos) const
        {return fabs(SignedDistanceToPlane(pos));}
        inline bool Inside(const Vector3d &pos) const
        {return SignedDistanceToPlane(pos)<0;}
        int ReflectAgainstBoundary(const Vector3d &originalPoint, 
                                         Vector3d &reflectedPoint, 
                                         Vector3d &boundaryPoint, 
                                         Vector3d &erectedNormal, 
                                         REAL &distanceTravelled); 
};
using FDTD_PlaneConstraint_Ptr = std::shared_ptr<FDTD_PlaneConstraint>;
#endif
