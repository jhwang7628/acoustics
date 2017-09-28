#include <wavesolver/HarmonicVibrationalSource.h> 

//##############################################################################
//##############################################################################
REAL HarmonicVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle)
{
    if (time<_t0 || time>_t1)
        return 0.0;
    // monopole vibration consistent source
    //const REAL scale = _omega*_omega*0.00025;
    //const REAL val = sin(_omega*time + _phase) * scale;

    // dipole vibration consistent source
    const REAL scale = _omega*_omega*0.00025;
    const REAL val = sin(_omega*time + _phase) *scale*normal[2];
      
    // normal harmonic vibration
    //const REAL val = sin(_omega*time + _phase);
    if (_decayRate < 0)
        return val; 
    else
        return decayFactor(time) * val;
}

//##############################################################################
//##############################################################################
REAL HarmonicVibrationalSource::
Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)
{
    return Evaluate(Vector3d(), vertexNormal, time);
}

//##############################################################################
//##############################################################################
REAL HarmonicVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
}

//##############################################################################
//##############################################################################
REAL HarmonicVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
}
