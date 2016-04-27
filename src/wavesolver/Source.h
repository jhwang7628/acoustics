#ifndef WAVESOLVER_SOURCE_H 
#define WAVESOLVER_SOURCE_H 

#include <TYPES.h> 
#include <config.h> 
#include <wavesolver/FDTD_MovableObject.h>

//##############################################################################
// Source term for wave solver
//##############################################################################
class Source : public FDTD_MovableObject
{
    public: 
        Source(){}

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time)=0; 
        virtual void UpdateBoundingBox(); 
};

#endif 

