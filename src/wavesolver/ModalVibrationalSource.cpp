#include <wavesolver/ModalVibrationalSource.h> 
#include <wavesolver/FDTD_RigidSoundObject.h>

ModalVibrationalSource::
ModalVibrationalSource(RigidObjectPtr owner)
    : VibrationalSource(owner) 
{
    _modalObjectOwner = std::static_pointer_cast<FDTD_RigidSoundObject>(owner);
}

//##############################################################################
//##############################################################################
REAL ModalVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    return _modalObjectOwner->SampleModalAcceleration(position, normal, time); 
}



//##############################################################################
//##############################################################################
REAL ModalVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    return _modalObjectOwner->SampleModalVelocity(position, normal, time); 
}

//##############################################################################
//##############################################################################
REAL ModalVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    return _modalObjectOwner->SampleModalDisplacement(position, normal, time); 
}
