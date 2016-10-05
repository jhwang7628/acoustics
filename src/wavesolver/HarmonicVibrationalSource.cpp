#include <wavesolver/HarmonicVibrationalSource.h> 

//##############################################################################
//##############################################################################
REAL HarmonicVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    const REAL val = cos(_omega*time + _phase); 
    if (_decayRate < 0)
        return val; 
    else
        return decayFactor(time) * val;
}



//##############################################################################
//##############################################################################
REAL HarmonicVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    const REAL val = 1./_omega * sin(_omega*time + _phase);
    if (_decayRate < 0)
        return val;
    else 
        return decayFactor(time) * val;
}

//##############################################################################
//##############################################################################
REAL HarmonicVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
}
