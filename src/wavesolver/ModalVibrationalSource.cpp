#include <wavesolver/ModalVibrationalSource.h> 
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>

//##############################################################################
//##############################################################################
ModalVibrationalSource::
ModalVibrationalSource(RigidObjectPtr owner)
    : VibrationalSource(owner) 
{
    _modalObjectOwner = std::static_pointer_cast<FDTD_RigidSoundObject>(owner);
}

//##############################################################################
//##############################################################################
REAL ModalVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle)
{
    if (!_modalObjectOwner->IsModalObject())
        return 0.0;
#if DEBUG_PERFECT_MODAL_HARMONICS == 0
    return _modalObjectOwner->SampleModalAcceleration(position, normal, time, hintTriangle); 
#else
    return _modalObjectOwner->PerfectHarmonics_SampleModalAcceleration(DEBUG_PERFECT_MODAL_HARMONICS, position, normal, time); 
#endif
}

//##############################################################################
//##############################################################################
REAL ModalVibrationalSource::
Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)
{
    if (!_modalObjectOwner->IsModalObject())
        return 0.0;
#if DEBUG_PERFECT_MODAL_HARMONICS == 0
    return _modalObjectOwner->SampleModalAcceleration(vertexID, vertexNormal, time); 
#else
    return _modalObjectOwner->PerfectHarmonics_SampleModalAcceleration(DEBUG_PERFECT_MODAL_HARMONICS, vertexID, vertexNormal, time); 
#endif
}

//##############################################################################
//##############################################################################
Vector3d ModalVibrationalSource::
Evaluate(const int &vertexID, const REAL &time)
{
    if (!_modalObjectOwner->IsModalObject())
        return Vector3d();
#if DEBUG_PERFECT_MODAL_HARMONICS == 0
    return _modalObjectOwner->SampleModalAcceleration(vertexID, time); 
#else
    return _modalObjectOwner->PerfectHarmonics_SampleModalAcceleration(DEBUG_PERFECT_MODAL_HARMONICS, vertexID, time); 
#endif
}

//##############################################################################
//##############################################################################
REAL ModalVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    if (!_modalObjectOwner->IsModalObject())
        return 0.0;
#if DEBUG_PERFECT_MODAL_HARMONICS == 0
    return _modalObjectOwner->SampleModalVelocity(position, normal, time); 
#else
    return _modalObjectOwner->PerfectHarmonics_SampleModalVelocity(DEBUG_PERFECT_MODAL_HARMONICS, position, normal, time); 
#endif
}

//##############################################################################
//##############################################################################
REAL ModalVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    if (!_modalObjectOwner->IsModalObject())
        return 0.0;
    return _modalObjectOwner->SampleModalDisplacement(position, normal, time); 
}

//##############################################################################
//##############################################################################
REAL ModalVibrationalSource::
EarliestEventTime(const REAL startTime) const
{
    return _modalObjectOwner->GetFirstImpulseTime(startTime); 
}
