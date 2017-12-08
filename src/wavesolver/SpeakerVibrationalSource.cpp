#include <wavesolver/SpeakerVibrationalSource.h> 

//##############################################################################
//##############################################################################
REAL SpeakerVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle)
{
    // Should not be called
    throw std::runtime_error("**ERROR** this is not implemented");
    return 0.0;
}

//##############################################################################
//##############################################################################
REAL SpeakerVibrationalSource::
Evaluate(const int &vertexID, const Vector3d &normal, const REAL &time)
{
    return 0.0;
}

//##############################################################################
//##############################################################################
Vector3d SpeakerVibrationalSource::
Evaluate(const int &vertexID, const REAL &time)
{
    return Vector3d(1.0, 1.0, 1.0);
}

//##############################################################################
//##############################################################################
REAL SpeakerVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
    return 0.0;
}

//##############################################################################
//##############################################################################
REAL SpeakerVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
    return 0.0;
}
