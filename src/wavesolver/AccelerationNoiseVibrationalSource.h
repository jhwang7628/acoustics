#ifndef ACCELERATION_NOISE_VIBRATIONAL_SOURCE_H 
#define ACCELERATION_NOISE_VIBRATIONAL_SOURCE_H 

#include <TYPES.h> 
#include <wavesolver/Wavesolver_ConstantsAndTypes.h> 
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/FDTD_RigidObject.h> 
#include <wavesolver/FDTD_RigidSoundObject.h> 

//##############################################################################
// This class handles source evaluation for modal objects.
//##############################################################################
class AccelerationNoiseVibrationalSource : public VibrationalSource
{
    private: 
        RigidSoundObjectPtr _modalObjectOwner; 
    public:
        AccelerationNoiseVibrationalSource(RigidObjectPtr owner);

        static REAL ComputeS(const ImpulseSeriesObject::ImpactRecord &impulse, const REAL &time); 
        static REAL ComputeSDot(const ImpulseSeriesObject::ImpactRecord &impulse, const REAL &time); 

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle=-1); 
        virtual REAL Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time); 
        virtual Vector3d Evaluate(const int &vertexID, const REAL &time); 
        virtual REAL EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        virtual REAL EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time); 

        ///// debug method ///// 
        virtual REAL EvaluatePressureAnalytical(const Vector3d &position, const Vector3d &normal, const REAL &time, const REAL &density, const REAL &soundSpeed, const REAL &sphereRadius); 
};

#endif
