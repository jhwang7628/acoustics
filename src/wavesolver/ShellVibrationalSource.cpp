#include "wavesolver/ShellVibrationalSource.h"

//##############################################################################
// Function Evaluate
//##############################################################################
REAL ShellVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, 
         const int &hintTriangle)
{
}

//##############################################################################
// Function Evaluate
//##############################################################################
REAL ShellVibrationalSource::
Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)
{
    const Vector3d acc = Evaluate(vertexId, time); 
    return acc.dotProduct(vertexNormal); 
}

//##############################################################################
// Function Evaluate
//##############################################################################
Vector3d ShellVibrationalSource::
Evaluate(const int &vertexID, const REAL &time)
{
    _owner->UpdatePosAcc(time); 
    return _owner->GetVertexAcc(vertexID); 
}

//##############################################################################
// Function EvaluateVelocity
//##############################################################################
REAL ShellVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, 
                 const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
}

//##############################################################################
// Function EvaluateVelocity
//##############################################################################
REAL ShellVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, 
                     const REAL &time)
{
    throw std::runtime_error("**ERROR** not implemented"); 
}
