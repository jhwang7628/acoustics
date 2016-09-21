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
    std::vector<ImpulseSeriesObject::ImpactRecord> impactRecords; 
    _modalObjectOwner->GetImpulseWithinSupport(time, impactRecords); 

    Vector3d acceleration(0, 0, 0); 
    for (const auto &impulse : impactRecords) 
    {
        if (impulse.supportLength < SMALL_NUM)
            continue;

        const REAL S = sin(M_PI*(time - impulse.timestamp)/impulse.supportLength); 
        // translational acceleration
        acceleration += impulse.impactVector * (M_PI*S / (2.0*impulse.supportLength*_modalObjectOwner->Mass())); 
    }

    const REAL a_n = acceleration.dotProduct(normal); 
    return a_n;
}



//##############################################################################
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    std::cerr << "**WARNING** Evaluate AN velocity not implemented\n"; 
    return 0.0; 
}

//##############################################################################
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
EvaluateDisplacement(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
    std::cout << "**WARNING** Evaluate AN displacement\n"; 
    return 0.0; 
}
