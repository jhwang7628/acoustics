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
        Source()
            : FDTD_MovableObject(SOURCE)
        {}

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle=-1)=0;
        virtual REAL Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)=0;
        virtual Vector3d Evaluate(const int &vertexID, const REAL &time)=0;
        virtual void UpdateBoundingBox();
        virtual void ResetUnionBox();
        // earliest nonzero source evaluation time in [startTime, infty)
        // default is no thresholding
        virtual REAL EarliestEventTime(const REAL startTime) const
        {return startTime;}
        // return whether the shader is zero within the bounding box.
        // this is used for the adaptive time-parallelization
        // default behavior is return false, meaning shader is dense in time
        virtual bool IsZero(const REAL t,
                            const bool checkBound = false,
                            const Vector3d &minBound = Vector3d(),
                            const Vector3d &maxBound = Vector3d()) const
        {return false;}

        virtual bool UpdateTime(const REAL time) {}; // default does nothing, bubbles object needs this to update the mesh
};

#endif

