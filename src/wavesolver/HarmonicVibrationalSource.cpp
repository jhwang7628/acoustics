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
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time) // TODO! 
{
    throw std::runtime_error("**ERROR** not implemented"); 
}

//##############################################################################
//##############################################################################
REAL HarmonicVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time) // TODO!
{
    throw std::runtime_error("**ERROR** not implemented"); 
}
