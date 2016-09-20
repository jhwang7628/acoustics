#include <wavesolver/AccelerationNoiseVibrationalSource.h> 
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>

//##############################################################################
//##############################################################################
AccelerationNoiseVibrationalSource::
AccelerationNoiseVibrationalSource(RigidObjectPtr owner)
    : VibrationalSource(owner) 
{
    _modalObjectOwner = std::static_pointer_cast<FDTD_RigidSoundObject>(owner);
}

//##############################################################################
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    std::cout << "Evaluate AN acceleration\n"; 
    return 0.0; 
}



//##############################################################################
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    std::cout << "Evaluate AN velocity\n"; 
    return 0.0; 
}

//##############################################################################
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    std::cout << "Evaluate AN displacement\n"; 
    return 0.0; 
}
