#include <wavesolver/AccelerationNoiseVibrationalSource.h> 
#include <wavesolver/FDTD_RigidSoundObject.h>
#include <wavesolver/Wavesolver_ConstantsAndTypes.h>

#define USE_GAUSSIAN_APPROXIMATION 0

//##############################################################################
//##############################################################################
AccelerationNoiseVibrationalSource::
AccelerationNoiseVibrationalSource(RigidObjectPtr owner)
    : VibrationalSource(owner) 
{
    _modalObjectOwner = std::static_pointer_cast<FDTD_RigidSoundObject>(owner);
}

//##############################################################################
// S is the approximation of half-sine pulse. 
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
ComputeS(const ImpulseSeriesObject::ImpactRecord &impulse, const REAL &time)
{
#if USE_GAUSSIAN_APPROXIMATION == 0
    if (time <= impulse.timestamp+impulse.supportLength && time >= impulse.timestamp)
        return sin(M_PI*(time - impulse.timestamp)/impulse.supportLength); 
    else
        return 0.0;
#else
    return exp(-6.0*pow((time - impulse.timestamp - impulse.supportLength/2.0)/impulse.supportLength, 2)); 
#endif
}

//##############################################################################
// S_dot is the approximation of half-sine pulse time derivative. 
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
ComputeSDot(const ImpulseSeriesObject::ImpactRecord &impulse, const REAL &time)
{
#if USE_GAUSSIAN_APPROXIMATION == 0
    if (time <= impulse.timestamp+impulse.supportLength && time >= impulse.timestamp)
      return M_PI/impulse.supportLength * cos(M_PI*(time - impulse.timestamp)/impulse.supportLength); 
    else
        return 0.0;
#else
    return -12.0 / pow(impulse.supportLength, 2) * (time - impulse.timestamp - impulse.supportLength/2.0) * ComputeS(impulse, time); 
#endif
}

//##############################################################################
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
Evaluate(const Vector3d &position, const Vector3d &normal, const REAL &time, const int &hintTriangle)
{
    std::vector<ImpulseSeriesObject::ImpactRecord> impactRecords; 
    _modalObjectOwner->GetImpulseWithinSupport(time, impactRecords); 

    Vector3d acceleration(0, 0, 0); 
    for (const auto &impulse : impactRecords) 
    {
        if (impulse.supportLength < SMALL_NUM)
            continue;

        const Vector3d &J = impulse.impactVector; 
        const Vector3d r = impulse.impactPosition - _modalObjectOwner->CenterOfMass(); 
        const REAL S = ComputeS(impulse, time); 
        // translational acceleration
        acceleration += J * (M_PI*S / (2.0*impulse.supportLength*_modalObjectOwner->Mass())); 
        // rotational acceleration (Eq 13)
        const Vector3d alpha = _modalObjectOwner->PremultiplyInvInertiaTensor(r.crossProduct(J)) * (M_PI*S) / (2.0*impulse.supportLength); 
        acceleration += alpha.crossProduct(r); 
    }
    const REAL a_n = _modalObjectOwner->ObjectToWorldVector(acceleration).dotProduct(normal); 
    return a_n;
}

//##############################################################################
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
Evaluate(const int &vertexID, const Vector3d &vertexNormal, const REAL &time)
{
    return Evaluate(Vector3d(), vertexNormal, time);
}

//##############################################################################
// The velocity evaluation of AN source is ill-defined for the fdtd simulator.
//##############################################################################
REAL AccelerationNoiseVibrationalSource::
EvaluateVelocity(const Vector3d &position, const Vector3d &normal, const REAL &time)
{
#if 0 // using the integral of acceleration profile, but I don't think this is correct.
    std::vector<ImpulseSeriesObject::ImpactRecord> impactRecords; 
    _modalObjectOwner->GetImpulseWithinSupport(time, impactRecords); 

    Vector3d velocity(0, 0, 0); 
    for (const auto &impulse : impactRecords) 
    {
        if (impulse.supportLength < SMALL_NUM || time < impulse.timestamp)
            continue;

        if (time <= impulse.timestamp+impulse.supportLength)
            velocity += (impulse.impactVector / (2.0*_modalObjectOwner->Mass())) * (1.0 - cos(M_PI*(time - impulse.timestamp)/impulse.supportLength));
        else 
            velocity += impulse.impactVector / _modalObjectOwner->Mass(); 
    }

    const REAL v_n = _modalObjectOwner->ObjectToWorldVector(velocity).dotProduct(normal); 
    return v_n;
#endif 
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
    const REAL delayedTime = time - r/soundSpeed; 
    if (delayedTime < 0) // causality condition
        return 0.0; 
    const REAL rho_a3_over_2c_r = density * pow(sphereRadius, 3) / (2.0*soundSpeed*r); 

    std::vector<ImpulseSeriesObject::ImpactRecord> impactRecords; 
    _modalObjectOwner->GetImpulseWithinSupport(delayedTime, impactRecords); 

    REAL p = 0.0; 
    for (const auto &impulse : impactRecords) 
    {
        if (impulse.supportLength < SMALL_NUM)
            continue;

        const REAL S_dot = ComputeSDot(impulse, delayedTime); 
        const REAL cos_t = impulse.impactVector.dotProduct(x-x0) / impulse.impactVector.norm() / r; 
        // acceleration jerk
        const REAL jerk = M_PI*impulse.impactVector.norm() / (2.0*impulse.supportLength*_modalObjectOwner->Mass()) * S_dot; 
        p += rho_a3_over_2c_r * cos_t * jerk; 
    }

    return p; 
}
