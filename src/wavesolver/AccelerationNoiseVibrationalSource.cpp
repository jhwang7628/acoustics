#include <wavesolver/AccelerationNoiseVibrationalSource.h> 
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>

#define USE_GAUSSIAN_APPROXIMATION 1

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
ComputeS(const ImpulseSeriesObject::ImpactRecord &impulse, const REAL &time)
{
#if USE_GAUSSIAN_APPROXIMATION == 0
    return sin(M_PI*(time - impulse.timestamp)/impulse.supportLength); 
#else
    return exp(-6.0*pow((time - impulse.timestamp - impulse.supportLength/2.0)/impulse.supportLength, 2)); 
#endif
}

//##############################################################################
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
ComputeSDot(const ImpulseSeriesObject::ImpactRecord &impulse, const REAL &time)
{
#if USE_GAUSSIAN_APPROXIMATION == 0
    throw std::runtime_error("**ERROR** not implemented SDot computation for half sine"); 
    return 0.0; 
#else
    return -12.0 / pow(impulse.supportLength, 2)*(time - impulse.timestamp - impulse.supportLength/2.0)*ComputeS(impulse, time); 
#endif
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

        const REAL S = ComputeS(impulse, time); 
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

//##############################################################################
// This function evaluates the analytical expression for sphere acceleration
// noise. For more detail please see file
// "CS448Z_2015fa_HW1_RBD_AccelNoise.pdf" and Chadwick's PAN paper. 
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
EvaluatePressureAnalytical(const Vector3d &position, const Vector3d &normal, const REAL &time, const REAL &density, const REAL &soundSpeed, const REAL &sphereRadius)
{
    const Vector3d &x = position; // listening point 
    const Vector3d x0 = _modalObjectOwner->MeshCentroid(); 
    const REAL r = (x - x0).norm(); 
    const REAL rho_a3_over_2c_r = density * pow(sphereRadius, 3) / (2.0*soundSpeed*r); 

    std::vector<ImpulseSeriesObject::ImpactRecord> impactRecords; 
    _modalObjectOwner->GetImpulseWithinSupport(time, impactRecords); 

    REAL p = 0.0; 
    for (const auto &impulse : impactRecords) 
    {
        if (impulse.supportLength < SMALL_NUM)
            continue;

        const REAL delayedTime = time - r/soundSpeed; 
        const REAL S_dot = ComputeSDot(impulse, delayedTime); 
        const REAL cos_t = impulse.impactVector.dotProduct(x-x0) / impulse.impactVector.norm() / r; 
        // acceleration jerk
        const REAL jerk = M_PI*impulse.impactVector.norm() / (2.0*impulse.supportLength*_modalObjectOwner->Mass()) * S_dot; 
        p += rho_a3_over_2c_r * cos_t * jerk; 
    }

    return p; 
}
