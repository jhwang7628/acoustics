#include <wavesolver/WaterVibrationalSource.h> 
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>

//##############################################################################
//##############################################################################
WaterVibrationalSource::
WaterVibrationalSource(RigidObjectPtr owner)
    : VibrationalSource(owner), _surfaceMesh(owner->GetMeshPtr())
{
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    return 0.0; 
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    return 0.0; 
}

//##############################################################################
//##############################################################################
REAL WaterVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    return 0.0; 
}

//##############################################################################
//##############################################################################
void WaterVibrationalSource::
Initialize(const std::string &wavFile)
{
    ReadOscillatorFromWav(wavFile); 
    ComputeVelocityAndAcceleration();
}

//##############################################################################
//##############################################################################
void WaterVibrationalSource::
ReadOscillatorFromWav(const std::string &wavFile)
{
}

//##############################################################################
//##############################################################################
void WaterVibrationalSource::
ComputeVelocityAndAcceleration()
{
}
