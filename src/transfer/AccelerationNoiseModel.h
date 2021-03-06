//////////////////////////////////////////////////////////////////////
// AccelerationNoiseModel.h: Interface for the AccelerationNoiseModel
//                           class
//
//////////////////////////////////////////////////////////////////////

#ifndef ACCELERATION_NOISE_MODEL_H
#define ACCELERATION_NOISE_MODEL_H

#include <geometry/ClosestPointMesh.h>
#include <geometry/GTS_TriMesh.h>
#include <geometry/RigidMesh.h>
#include <geometry/TriangleMesh.hpp>

#include <io/ImpulseIO.h>

#include <linearalgebra/Quaternion.hpp>

#include <transfer/ProxyManager.h>
#include <transfer/PulseApproximation.h>

#include <TYPES.h>

#include <fstream>

//////////////////////////////////////////////////////////////////////
// AccelerationNoiseModel class
//
// A few routines for modelling acceleration noise from an object
//////////////////////////////////////////////////////////////////////
class AccelerationNoiseModel {
    public:
        struct MeshSet {
            MeshSet( const TriangleMesh<REAL> &mesh, const GTS_TriMesh &curvatureMesh,
                     const RigidMesh &rigidMesh,
                     const ClosestPointMesh &closestPointMesh,
                     int objectID, bool useProxy )
                : _mesh( mesh ),
                  _curvatureMesh( curvatureMesh ),
                  _rigidMesh( rigidMesh ),
                  _closestPointMesh( closestPointMesh ),
                  _objectID( objectID ),
                  _useProxy( useProxy )
            {
            }

            const TriangleMesh<REAL>    &_mesh;
            const GTS_TriMesh           &_curvatureMesh;
            const RigidMesh             &_rigidMesh;
            const ClosestPointMesh      &_closestPointMesh;

            int                          _objectID;
            bool                         _useProxy;
        };

        // For a single impact, adds the acceleration noise contribution
        // to the given output signal.
        //
        // We require a GTS_TriMesh, ClosestPointMesh and RigidMesh as
        // well as a pulse approximation model for all rigid body directions
        //
        // Also need the point position of the impact in each mesh, as well
        // as the relative collision velocity, impact magnitude, contact normal
        // and the current rotation state for each object.
        static bool AddImpactNoise( const MeshSet &objectA,
                                    const MeshSet &objectB,
                                    const ImpulseIO::TwoObjectImpact &impactData,
                                    const Point3d &listeningPosition,
                                    /*
                                    RadialApproximation::AccelerationSet &pulseDataA,
                                    RadialApproximation::AccelerationSet &pulseDataB,
                                    */
                                    CompactRadialApproximation *pulseDataA,
                                    CompactRadialApproximation *pulseDataB,
                                    SampledFunction &outputSignal,
                                    std::ofstream *outputFile = NULL,
                                    REAL *impactLength = NULL,
                                    REAL *impactScale = NULL,
                                    ProxyManager *proxyManager = NULL,
                                    REAL collisionTimeScale = 1.0 );

        // Adds impact noise for a single impact with a plane of infinite
        // size and mass
        static bool AddImpactNoise( const MeshSet &objectA,
                                    const ImpulseIO::PlaneImpact &impactData,
                                    const Point3d &listeningPosition,
                                    /*
                                    RadialApproximation::AccelerationSet &pulseDataA,
                                    */
                                    CompactRadialApproximation *pulseDataA,
                                    SampledFunction &outputSignal,
                                    std::ofstream *outputFile = NULL,
                                    REAL *impactLength = NULL,
                                    REAL *impactScale = NULL,
                                    ProxyManager *proxyManager = NULL,
                                    REAL collisionTimeScale = 1.0 );

        // Computes the force time scale and appropriate force scaling for this collision
        //
        // Optionally rescale by the given factor
        typedef std::pair<REAL, REAL> HertzImpactData;
        static HertzImpactData ImpactTimeScale( const MeshSet &objectA,
                                                const MeshSet &objectB,
                                                const ImpulseIO::TwoObjectImpact &impactData,
                                                REAL collisionTimeScale = 1.0 );

        // Computes the force time scale and appropriate force scaling for this collision
        //
        // Optionally rescale by the given factor
        static HertzImpactData ImpactTimeScale( const MeshSet &objectA,
                                                const ImpulseIO::PlaneImpact &impactData,
                                                REAL collisionTimeScale = 1.0 );


    protected:

    private:
        AccelerationNoiseModel();

        // Destructor
        virtual ~AccelerationNoiseModel();

        // Uses Hertz impact theory to compute a time scale for the
        // force pulse to be applied to the object.
        //
        // inverseEffectiveMass will be set by this function
        //
        // See [Johnson, 1985] page 353
        static REAL HertzImpactTimeScale( const MeshSet &objectA,
                                          const MeshSet &objectB,
                                          const ImpulseIO::TwoObjectImpact &impactData,
                                          const Vector3d &impulseDirectionA,
                                          const Vector3d &impulseDirectionB,
                                          REAL &inverseEffectiveMass );

        // Hertz time scale for an impact with an ideally flat surface (eg.
        // a ground plane)
        //
        // inverseEffectiveMass will be set by this function
        static REAL HertzImpactTimeScale( const MeshSet &objectA,
                                          const ImpulseIO::PlaneImpact &impactData,
                                          const Vector3d &impulseDirectionA,
                                          REAL &inverseEffectiveMass );


        // Computes the 1 / E^{*} quantity from [Johnson, 1985] page 352
        static REAL HertzMaterialConstant( REAL youngsModulusA,
                                           REAL youngsModulusB,
                                           REAL poissonRatioA,
                                           REAL poissonRatioB );

        // Adds impact noise for a single object
        static bool AddObjectImpactNoise( REAL startTime, REAL timeScale, REAL forceScale,
                                          const Vector3d &translationalAcceleration,
                                          const Vector3d &rotationalAcceleration,
                                          const Point3d &listeningPosition,
                                          RadialApproximation::AccelerationSet &pulseData,
                                          SampledFunction &outputSignal,
                                          REAL scale = 1.0 );

        // Adds impact noise for a single object
        static bool AddObjectImpactNoise( REAL startTime, REAL timeScale, REAL forceScale,
                                          const Vector3d &translationalAcceleration,
                                          const Vector3d &rotationalAcceleration,
                                          const Point3d &listeningPosition,
                                          CompactRadialApproximation *pulseData,
                                          SampledFunction &outputSignal,
                                          REAL scale = 1.0 );

        static REAL ImpactForceScale( REAL inverseEffectiveMass,
                                      REAL hertzTimeScale,
                                      REAL impulseScale )
        {
#if 0
            return ( M_PI * impulseScale )
                / ( 2.0 * hertzTimeScale * inverseEffectiveMass );
#endif
            // Match the time-averaged forces
            return ( M_PI * impulseScale ) / ( 2.0 * hertzTimeScale );
        }

        static Point3d ObjectListeningPosition( const Point3d &listeningPosition,
                                                const Point3d &restCenterOfMass,
                                                const Point3d &impactCenterOfMass,
                                                const Quat4d &impactInverseRotation )
        {
            Vector3d                centerOfMassOffset;

            // Offset from the current center of mass
            centerOfMassOffset = listeningPosition - impactCenterOfMass;

            // Rotate out of the impact frame
            centerOfMassOffset = impactInverseRotation.rotate( centerOfMassOffset );

            // Add back the original center of mass
            centerOfMassOffset += restCenterOfMass;

            return Point3d(centerOfMassOffset);
        }

};

#endif
