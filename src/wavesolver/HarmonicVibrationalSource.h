#ifndef HARMONIC_VIBRATIONAL_SOURCE_H 
#define HARMONIC_VIBRATIONAL_SOURCE_H 

#include <TYPES.h> 
#include <wavesolver/VibrationalSource.h> 
#include <wavesolver/FDTD_RigidObject.h> 

//##############################################################################
// Harmonic vibrations. 
//##############################################################################
class HarmonicVibrationalSource : public VibrationalSource
{
    private: 
        REAL    _omega; 
        REAL    _phase;

    public:
        HarmonicVibrationalSource(RigidObjectPtr owner, const REAL &omega, const REAL &phase) 
            : VibrationalSource(owner), _omega(omega), _phase(phase) 
        {}

        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        virtual REAL EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time); 
        virtual REAL EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time); 
};

#endif
