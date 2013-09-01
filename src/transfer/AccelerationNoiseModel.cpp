//////////////////////////////////////////////////////////////////////
// AccelerationNoiseModel.cpp: Implementation of the
//                             AccelerationNoiseModel class
//
//////////////////////////////////////////////////////////////////////

#include "AccelerationNoiseModel.h"

#include <utils/IO.h>
#include <utils/trace.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
AccelerationNoiseModel::AccelerationNoiseModel()
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
AccelerationNoiseModel::~AccelerationNoiseModel()
{
}

//////////////////////////////////////////////////////////////////////
// For a single impact, adds the acceleration noise contribution
// to the given output signal.
//
// We require a GTS_TriMesh, ClosestPointMesh and RigidMesh as
// well as a pulse approximation model for all rigid body directions
//
// Also need the point position of the impact in each mesh, as well
// as the relative collision velocity, impact magnitude, contact normal
// and the current rotation state for each object.
//////////////////////////////////////////////////////////////////////
bool AccelerationNoiseModel::AddImpactNoise( const MeshSet &objectA,
                                             const MeshSet &objectB,
                                             const TwoObjectImpact &impactData,
                                             const Point3d &listeningPosition,
                                             /*
                                             RadialApproximation::AccelerationSet &pulseDataA,
                                             RadialApproximation::AccelerationSet &pulseDataB,
                                             */
                                             CompactRadialApproximation *pulseDataA,
                                             CompactRadialApproximation *pulseDataB,
                                             SampledFunction &outputSignal,
                                             ofstream *outputFile,
                                             REAL *impactLength,
                                             REAL *impactScale,
                                             ProxyManager *proxyManager,
                                             REAL collisionTimeScale )
{
    REAL                       timeScale;
    Vector3d                   impulseDirectionA;
    Vector3d                   impulseDirectionB;

    REAL                       inverseEffectiveMass;
    REAL                       forceScale;

    Vector3d                   translationalAccelerationA;
    Vector3d                   translationalAccelerationB;
    Vector3d                   rotationalAccelerationA;
    Vector3d                   rotationalAccelerationB;

    Point3d                    listeningPositionA;
    Point3d                    listeningPositionB;

    bool                       addingSound = false;

    REAL                       scaleA = 1.0;
    REAL                       scaleB = 1.0;

    impulseDirectionA.normalize();
    impulseDirectionB.normalize();

    // Get the Hertz timescale, and effective masses for both objects
    // in the given collision configuration
    impulseDirectionA
        = impactData._inverseRotationA.rotate( impactData._impulseDirection );
    impulseDirectionB
        = impactData._inverseRotationB.rotate( impactData._impulseDirection );

    // Figure out the acceleration magnitude in each direction
    // given a unit force
    objectA._rigidMesh.accelerationCoefficients( impactData._posA,
                                                 impulseDirectionA,
                                                 translationalAccelerationA,
                                                 rotationalAccelerationA );
    objectB._rigidMesh.accelerationCoefficients( impactData._posB,
                                                 impulseDirectionB,
                                                 translationalAccelerationB,
                                                 rotationalAccelerationB );

    // Get listening positions for each of the objects
    listeningPositionA = ObjectListeningPosition(
            listeningPosition,
            objectA._rigidMesh.centerOfMass(),
            impactData._posA,
            impactData._inverseRotationA );
    listeningPositionB = ObjectListeningPosition(
            listeningPosition,
            objectB._rigidMesh.centerOfMass(),
            impactData._posB,
            impactData._inverseRotationB );

    if ( proxyManager ) {
        // Use proxies for each object
        if ( objectA._useProxy ) {
            CompactRadialApproximation  *proxyA = NULL;

            proxyManager->getProxyData( objectA._objectID, proxyA, scaleA,
                    listeningPositionA,
                    translationalAccelerationA,
                    rotationalAccelerationA );

            pulseDataA = proxyA;
#if 0
            pulseDataA = pulseDataA != NULL ? proxyA : NULL;
#endif
        }

        if ( objectB._useProxy ) {
            CompactRadialApproximation  *proxyB = NULL;

            proxyManager->getProxyData( objectB._objectID, proxyB, scaleB,
                    listeningPositionB,
                    translationalAccelerationB,
                    rotationalAccelerationB );

            pulseDataB = proxyB;
#if 0
            pulseDataB = pulseDataB != NULL ? proxyB : NULL;
#endif
        }
    }

    // Figure out the impact time scale
    timeScale = HertzImpactTimeScale( objectA, objectB, impactData,
            impulseDirectionA, impulseDirectionB,
            inverseEffectiveMass );
    timeScale *= collisionTimeScale;

    // FIXME: clamp the time scale
    if ( pulseDataA && pulseDataB ) {
        timeScale = max( timeScale, max( pulseDataA->h() * scaleA,
                    pulseDataB->h() * scaleB ) );
    }
    else if ( pulseDataA ) {
        timeScale = max( timeScale, pulseDataA->h() * scaleA );
    }
    else if ( pulseDataB ) {
        timeScale = max( timeScale, pulseDataB->h() * scaleB );
    }
    else {
        return addingSound;
    }

    forceScale = ImpactForceScale( inverseEffectiveMass, timeScale,
            impactData._impulseMagnitude );


#if 0
    // FIXME:
    if ( forceScale < 0.1 )
    {
        timeScale = -1.0;
        forceScale = -1.0;
    }
#endif

    // Clamp the force scale
    REAL forceScaleA = forceScale;
#if 0
    forceScaleA = min( forceScale,
            impactData._impulseMagnitude / pulseDataA->h() );
#endif

    // Clamp the force scale
    REAL forceScaleB = forceScale;
#if 0
    forceScaleB = min( forceScale,
            impactData._impulseMagnitude / pulseDataB->h() );
#endif

    if ( outputFile )
    {
        ( *outputFile ) << timeScale << ' ' << forceScale << endl;
    }

    if ( impactLength && impactScale )
    {
        ( *impactLength ) = timeScale;
        ( *impactScale ) = forceScale;
    }

#if 0
    // FIXME:
    if ( forceScale < 0.0 )
    {
        return true;
    }
#endif

    addingSound = addingSound
        || AddObjectImpactNoise( impactData._impulseTime, timeScale, forceScaleA,
                translationalAccelerationA,
                rotationalAccelerationA,
                listeningPositionA, pulseDataA, outputSignal,
                scaleA );

    addingSound = addingSound
        || AddObjectImpactNoise( impactData._impulseTime, timeScale, forceScaleB,
                translationalAccelerationB,
                rotationalAccelerationB,
                listeningPositionB, pulseDataB, outputSignal,
                scaleB );

#if 0
    cout << SDUMP( timeScale ) << endl;
    cout << SDUMP( forceScale ) << endl;
    cout << SDUMP( objectA._rigidMesh.mass() ) << endl;
    cout << SDUMP( impactData._impulseMagnitude ) << endl;
    cout << SDUMP( translationalAccelerationA * forceScale ) << endl;
    cout << SDUMP( rotationalAccelerationA * forceScale ) << endl;
#endif

    return addingSound;
}

//////////////////////////////////////////////////////////////////////
// Adds impact noise for a single impact with a plane of infinite
// size and mass
//////////////////////////////////////////////////////////////////////
bool AccelerationNoiseModel::AddImpactNoise( const MeshSet &objectA,
                                             const PlaneImpact &impactData,
                                             const Point3d &listeningPosition,
                                             /*
                                             RadialApproximation::AccelerationSet &pulseDataA,
                                             */
                                             CompactRadialApproximation *pulseDataA,
                                             SampledFunction &outputSignal,
                                             ofstream *outputFile,
                                             REAL *impactLength,
                                             REAL *impactScale,
                                             ProxyManager *proxyManager,
                                             REAL collisionTimeScale )
{
    REAL                       timeScale;
    Vector3d                   impulseDirectionA;

    REAL                       inverseEffectiveMass;
    REAL                       forceScale;

    Vector3d                   translationalAccelerationA;
    Vector3d                   rotationalAccelerationA;

    Point3d                    listeningPositionA;

    REAL                       scaleA = 1.0;

    // Get the Hertz timescale, and effective masses for both objects
    // in the given collision configuration
    impulseDirectionA
        = impactData._inverseRotationA.rotate( impactData._impulseDirection );

    // Figure out the acceleration magnitude in each direction
    // given a unit force
    objectA._rigidMesh.accelerationCoefficients( impactData._posA,
            impulseDirectionA,
            translationalAccelerationA,
            rotationalAccelerationA );

    // Get listening positions for each of the objects
    listeningPositionA = ObjectListeningPosition(
            listeningPosition,
            objectA._rigidMesh.centerOfMass(),
            impactData._posA,
            impactData._inverseRotationA );

    if ( proxyManager && objectA._useProxy ) {
        // Use proxies for each object
        CompactRadialApproximation  *proxyA = NULL;

        proxyManager->getProxyData( objectA._objectID, proxyA, scaleA,
                listeningPositionA,
                translationalAccelerationA,
                rotationalAccelerationA );

        pulseDataA = proxyA;
#if 0
        pulseDataA = pulseDataA != NULL ? proxyA : NULL;
#endif
    }

    // Figure out the impact time scale
    timeScale = HertzImpactTimeScale( objectA, impactData,
            impulseDirectionA,
            inverseEffectiveMass );
    timeScale *= collisionTimeScale;

    // FIXME: clamp the time scale
    if ( pulseDataA ) {
        timeScale = max( timeScale, pulseDataA->h() * scaleA );
    }
    else {
        return false;
    }

    forceScale = ImpactForceScale( inverseEffectiveMass, timeScale,
            impactData._impulseMagnitude );

#if 0
    // Clamp the force scale
    forceScale = min( forceScale, impactData._impulseMagnitude / pulseDataA->h() );
#endif

#if 0
    // FIXME:
    if ( forceScale < 0.1 )
    {
        timeScale = -1.0;
        forceScale = -1.0;
    }
#endif

    if ( outputFile )
    {
        ( *outputFile ) << timeScale << ' ' << forceScale << endl;
    }

    if ( impactLength && impactScale )
    {
        ( *impactLength ) = timeScale;
        ( *impactScale ) = forceScale;
    }

#if 0
    // FIXME:
    if ( forceScale < 0.0 )
    {
        return true;
    }
#endif

    return AddObjectImpactNoise( impactData._impulseTime, timeScale, forceScale,
            translationalAccelerationA,
            rotationalAccelerationA,
            listeningPositionA, pulseDataA, outputSignal,
            scaleA );
}

//////////////////////////////////////////////////////////////////////
// Computes the force time scale for this collision and writes it to the
// given output stream
//
// Optionally rescale by the given factor
//////////////////////////////////////////////////////////////////////
AccelerationNoiseModel::ImpactDataPair
AccelerationNoiseModel::ImpactTimeScale( const MeshSet &objectA,
                                         const MeshSet &objectB,
                                         const TwoObjectImpact &impactData,
                                         REAL collisionTimeScale )
{
    REAL                       timeScale;
    Vector3d                   impulseDirectionA;
    Vector3d                   impulseDirectionB;

    REAL                       inverseEffectiveMass;
    REAL                       forceScale;

#if 0
    impulseDirectionA.normalize();
    impulseDirectionB.normalize();
#endif

    // Get the Hertz timescale, and effective masses for both objects
    // in the given collision configuration
    impulseDirectionA
        = impactData._inverseRotationA.rotate( impactData._impulseDirection );
    impulseDirectionB
        = impactData._inverseRotationB.rotate( impactData._impulseDirection );

    // Figure out the impact time scale
    timeScale = HertzImpactTimeScale( objectA, objectB, impactData,
            impulseDirectionA, impulseDirectionB,
            inverseEffectiveMass );
    timeScale *= collisionTimeScale;

    forceScale = ImpactForceScale( inverseEffectiveMass, timeScale, impactData._impulseMagnitude );

    return ImpactDataPair( timeScale, forceScale );
}

//////////////////////////////////////////////////////////////////////
// Uses Hertz impact theory to compute a time scale for the
// force pulse to be applied to the object.
//
// See [Johnson, 1985] page 353
//////////////////////////////////////////////////////////////////////
REAL AccelerationNoiseModel::HertzImpactTimeScale( const MeshSet &objectA,
                                                   const MeshSet &objectB,
                                                   const TwoObjectImpact &impactData,
                                                   const Vector3d &impulseDirectionA,
                                                   const Vector3d &impulseDirectionB,
                                                   REAL &inverseEffectiveMass )
{
    ClosestPointMesh::MeshPoint    closestPointA;
    ClosestPointMesh::MeshPoint    closestPointB;

    REAL                           curvatureA;
    REAL                           curvatureB;

    REAL                           inverseEffectiveMassA;
    REAL                           inverseEffectiveMassB;

    REAL                           inverseEffectiveRadius;

    REAL                           materialConstant;

    REAL                           timeScale;

    TRACE_ASSERT( impactData._relativeSpeed < 0.0 );

    // For the given impact positions (which may be interior to the meshes)
    // determine a nearby point on the object surface so that we can
    // use this for curvature calculations
    closestPointA = objectA._closestPointMesh.closestPoint( impactData._posA );
    closestPointB = objectB._closestPointMesh.closestPoint( impactData._posB );

    if ( closestPointA._triangleID >= 0 )
    {
        curvatureA = objectA._curvatureMesh.meanCurvature(
                closestPointA._triangleID,
                closestPointA._barycentricPosition );
    }
    else
    {
#if 0
        cerr << "Object impact: No closest point found" << endl;
#endif
        curvatureA = 0.0;
    }

    if ( closestPointB._triangleID >= 0 )
    {
        curvatureB = objectB._curvatureMesh.meanCurvature(
                closestPointB._triangleID,
                closestPointB._barycentricPosition );
    }
    else
    {
#if 0
        cerr << "Object impact: No closest point found" << endl;
#endif
        curvatureB = 0.0;
    }

    inverseEffectiveRadius = curvatureA + curvatureB;

    // Get effective masses for each of the two objects, given the impulse
    // direction on each object
    inverseEffectiveMassA = objectA._rigidMesh.inverseEffectiveMass(
            impactData._posA, impulseDirectionA );
    inverseEffectiveMassB = objectB._rigidMesh.inverseEffectiveMass(
            impactData._posB, impulseDirectionB );

    inverseEffectiveMass = inverseEffectiveMassA + inverseEffectiveMassB;

    materialConstant = HertzMaterialConstant(
            objectA._rigidMesh.youngsModulus(),
            objectB._rigidMesh.youngsModulus(),
            objectA._rigidMesh.poissonRatio(),
            objectB._rigidMesh.poissonRatio() );

    timeScale = 1.0 / inverseEffectiveMass;
    timeScale *= timeScale;
    timeScale *= inverseEffectiveRadius;
    timeScale *= materialConstant * materialConstant;
    timeScale /= -1.0 * impactData._relativeSpeed;
    timeScale = pow( timeScale, 0.2 );
    timeScale *= 2.87;

    return timeScale;
}

//////////////////////////////////////////////////////////////////////
// Hertz time scale for an impact with an ideally flat surface (eg.
// a ground plane)
//////////////////////////////////////////////////////////////////////
REAL AccelerationNoiseModel::HertzImpactTimeScale( const MeshSet &objectA,
                                                   const PlaneImpact &impactData,
                                                   const Vector3d &impulseDirectionA,
                                                   REAL &inverseEffectiveMass )
{
    ClosestPointMesh::MeshPoint    closestPointA;

    REAL                           curvatureA;
    REAL                           curvatureB = 0.0;

    REAL                           inverseEffectiveMassA;
    REAL                           inverseEffectiveMassB = 0.0;

    REAL                           inverseEffectiveRadius;

    REAL                           materialConstant;

    REAL                           timeScale;

    TRACE_ASSERT( impactData._relativeSpeed < 0.0 );

    // For the given impact positions (which may be interior to the meshes)
    // determine a nearby point on the object surface so that we can
    // use this for curvature calculations
    closestPointA = objectA._closestPointMesh.closestPoint( impactData._posA );

#if 0
    if ( closestPointA._triangleID < 0 ) {
        cerr << "Plane impact: No closest point found" << endl;
    }
#endif

    curvatureA = objectA._curvatureMesh.meanCurvature(
            closestPointA._triangleID,
            closestPointA._barycentricPosition );

    inverseEffectiveRadius = curvatureA + curvatureB;

    // Get effective masses for each of the two objects, given the impulse
    // direction on each object
    inverseEffectiveMassA = objectA._rigidMesh.inverseEffectiveMass(
            impactData._posA, impulseDirectionA );

    inverseEffectiveMass = inverseEffectiveMassA + inverseEffectiveMassB;

    materialConstant = HertzMaterialConstant(
            objectA._rigidMesh.youngsModulus(),
            impactData._planeYoungsModulus,
            objectA._rigidMesh.poissonRatio(),
            impactData._planePoissonRatio );

    timeScale = 1.0 / inverseEffectiveMass;
    timeScale *= timeScale;
    timeScale *= inverseEffectiveRadius;
    timeScale *= materialConstant * materialConstant;
    timeScale /= -1.0 * impactData._relativeSpeed;
    timeScale = pow( timeScale, 0.2 );
    timeScale *= 2.87;

    return timeScale;
}

//////////////////////////////////////////////////////////////////////
// Computes the 1 / E^{*} quantity from [Johnson, 1985] page 352
//////////////////////////////////////////////////////////////////////
REAL AccelerationNoiseModel::HertzMaterialConstant( REAL youngsModulusA,
                                                    REAL youngsModulusB,
                                                    REAL poissonRatioA,
                                                    REAL poissonRatioB )
{
    return ( 1.0 - poissonRatioA * poissonRatioA ) / youngsModulusA
        + ( 1.0 - poissonRatioB * poissonRatioB ) / youngsModulusB;
}

//////////////////////////////////////////////////////////////////////
// Adds impact noise for a single object
//////////////////////////////////////////////////////////////////////
bool AccelerationNoiseModel::AddObjectImpactNoise( REAL startTime, REAL timeScale, REAL forceScale,
                                                   const Vector3d &translationalAcceleration,
                                                   const Vector3d &rotationalAcceleration,
                                                   const Point3d &listeningPosition,
                                                   RadialApproximation::AccelerationSet &pulseData,
                                                   SampledFunction &outputSignal,
                                                   REAL scale )
{
    REAL                       accelerationCoefficients[] = {
                                                    forceScale * translationalAcceleration[ 0 ],
                                                    forceScale * translationalAcceleration[ 1 ],
                                                    forceScale * translationalAcceleration[ 2 ],
                                                    forceScale * rotationalAcceleration[ 0 ],
                                                    forceScale * rotationalAcceleration[ 1 ],
                                                    forceScale * rotationalAcceleration[ 2 ] };

    return RadialApproximation::AddSinePulse( pulseData, listeningPosition,
            // Pulse starting time and length
            startTime, timeScale,
            // Amplitude in each direction
            accelerationCoefficients,
            outputSignal );

}

//////////////////////////////////////////////////////////////////////
// Adds impact noise for a single object
//////////////////////////////////////////////////////////////////////
bool AccelerationNoiseModel::AddObjectImpactNoise( REAL startTime, REAL timeScale, REAL forceScale,
                                                   const Vector3d &translationalAcceleration,
                                                   const Vector3d &rotationalAcceleration,
                                                   const Point3d &listeningPosition,
                                                   CompactRadialApproximation *pulseData,
                                                   SampledFunction &outputSignal,
                                                   REAL scale )
{
    REAL                       accelerationCoefficients[] = {
                                                    forceScale * translationalAcceleration[ 0 ],
                                                    forceScale * translationalAcceleration[ 1 ],
                                                    forceScale * translationalAcceleration[ 2 ],
                                                    forceScale * rotationalAcceleration[ 0 ],
                                                    forceScale * rotationalAcceleration[ 1 ],
                                                    forceScale * rotationalAcceleration[ 2 ] };

    if ( pulseData )
    {
        return pulseData->addSinePulse( listeningPosition,
                // Pulse starting time and length
                startTime, timeScale,
                // Amplitude in each direction
                accelerationCoefficients,
                outputSignal,
                scale );
    }

    return false;
}
