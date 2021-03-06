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
        REAL    _decayRate = -1.0;
        REAL    _decayCenter = 0.0;

        REAL    _t0 = std::numeric_limits<REAL>::lowest();
        REAL    _t1 = std::numeric_limits<REAL>::max();

    public:
        HarmonicVibrationalSource(RigidObjectPtr owner, const REAL &omega, const REAL &phase)
            : VibrationalSource(owner), _omega(omega), _phase(phase)
        {}

        HarmonicVibrationalSource(RigidObjectPtr owner, const REAL &omega, const REAL &phase, const REAL &decayRate, const REAL &decayCenter)
            : VibrationalSource(owner), _omega(omega), _phase(phase), _decayRate(decayRate), _decayCenter(decayCenter)
        {}

        inline void SetRange(const REAL t0, const REAL t1){_t0=t0; _t1=t1;}
        inline REAL decayFactor(const REAL &t){return exp(-_decayRate * (t-_decayCenter));}
        virtual REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle=-1);
        virtual REAL Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time);
        virtual Vector3d Evaluate(const int &vertexID, const REAL &time);
        virtual REAL EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time);
        virtual REAL EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time);
};

#endif
