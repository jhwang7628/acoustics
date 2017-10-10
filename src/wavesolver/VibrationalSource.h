#ifndef VIBRATIONAL_SOURCE_H 
#define VIBRATIONAL_SOURCE_H 

#include <TYPES.h> 
#include <linearalgebra/Vector3.hpp> 
#include <wavesolver/Source.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h> 

//##############################################################################
// Forward declaration
//##############################################################################
class FDTD_RigidObject; 

//##############################################################################
// Represents the surface acceleration that causes sound. 
//##############################################################################
class VibrationalSource : public Source
{
    protected:
        std::shared_ptr<FDTD_RigidObject> _owner; 

    public:
        VibrationalSource(){}
        VibrationalSource(RigidObjectPtr owner)
            : _owner(owner)
        {}

        // output the acceleration
        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle=-1)=0; 
        virtual REAL Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)=0; 
        virtual Vector3d Evaluate(const int &vertexID, const REAL &time)=0;
        virtual REAL EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)=0; 
        virtual REAL EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)=0; 
};

#endif
