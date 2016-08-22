#include <wavesolver/HarmonicVibrationalSource.h> 

//##############################################################################
//##############################################################################
REAL HarmonicVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    return cos(_omega*time + _phase); 
}



//##############################################################################
//##############################################################################
REAL HarmonicVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    const REAL returnValue = 1./_omega * sin(_omega*time + _phase); 
    return returnValue; 
}

//##############################################################################
//##############################################################################
REAL HarmonicVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
}
