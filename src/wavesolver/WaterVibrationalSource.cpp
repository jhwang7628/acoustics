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
WaterVibrationalSource::
REAL Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
}

//##############################################################################
//##############################################################################
WaterVibrationalSource::
REAL EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
}

//##############################################################################
//##############################################################################
WaterVibrationalSource::
REAL EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
}
