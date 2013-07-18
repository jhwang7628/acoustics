//////////////////////////////////////////////////////////////////////
// AccelerationNoiseModel.h: Interface for the AccelerationNoiseModel
//                           class
//
//////////////////////////////////////////////////////////////////////

#ifndef ACCELERATION_NOISE_MODEL_H
#define ACCELERATION_NOISE_MODEL_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/QUATERNION.h>
#include <linearalgebra/VECTOR.h>

#include <mesh/ClosestPointMesh.h>
#include <mesh/GTS_TriMesh.h>
#include <mesh/RigidMesh.h>
#include <mesh/TriMesh.h>

#include <transfer/ProxyManager.h>
#include <transfer/PulseApproximation.h>

#include <SETTINGS.h>
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
      MeshSet( const TriMesh &mesh, const GTS_TriMesh &curvatureMesh,
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

      const TriMesh           &_mesh;
      const GTS_TriMesh       &_curvatureMesh;
      const RigidMesh         &_rigidMesh;
      const ClosestPointMesh  &_closestPointMesh;

      int                      _objectID;
      bool                     _useProxy;
    };

    struct TwoObjectImpact {
      // Impact positions in each mesh, relative to the mesh's rest pose
      VEC3F                    _posA;
      VEC3F                    _posB;

      // Impulse information
      Real                     _impulseTime;
      VEC3F                    _impulseDirection;
      Real                     _impulseMagnitude;
      Real                     _relativeSpeed;

      // Rigid state information for the two meshes
      VEC3F                    _centerOfMassA;
      VEC3F                    _centerOfMassB;

      QUATERNION               _inverseRotationA;
      QUATERNION               _inverseRotationB;

    };

    struct PlaneImpact {
      // Impact position in the mesh
      VEC3F                    _posA;

      // Impulse information
      Real                     _impulseTime;
      VEC3F                    _impulseDirection;
      Real                     _impulseMagnitude;
      Real                     _relativeSpeed;

      // Rigid state information for the mesh
      VEC3F                    _centerOfMassA;
      QUATERNION               _inverseRotationA;

      // Also store material parameters for the ground plane here
      Real                     _planeYoungsModulus;
      Real                     _planePoissonRatio;

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
                                const TwoObjectImpact &impactData,
                                const VEC3F &listeningPosition,
                                /*
                                RadialApproximation::AccelerationSet &pulseDataA,
                                RadialApproximation::AccelerationSet &pulseDataB,
                                */
                                CompactRadialApproximation *pulseDataA,
                                CompactRadialApproximation *pulseDataB,
                                SampledFunction &outputSignal,
                                std::ofstream *outputFile = NULL,
                                Real *impactLength = NULL,
                                Real *impactScale = NULL,
                                ProxyManager *proxyManager = NULL,
                                Real collisionTimeScale = 1.0 );

    // Adds impact noise for a single impact with a plane of infinite
    // size and mass
    static bool AddImpactNoise( const MeshSet &objectA,
                                const PlaneImpact &impactData,
                                const VEC3F &listeningPosition,
                                /*
                                RadialApproximation::AccelerationSet &pulseDataA,
                                */
                                CompactRadialApproximation *pulseDataA,
                                SampledFunction &outputSignal,
                                std::ofstream *outputFile = NULL,
                                Real *impactLength = NULL,
                                Real *impactScale = NULL,
                                ProxyManager *proxyManager = NULL,
                                Real collisionTimeScale = 1.0 );

  protected:

  private:
    AccelerationNoiseModel();

    // Destructor
    virtual ~AccelerationNoiseModel();

    // Uses Hertz impact theory to compute a time scale for the
    // force pulse to be applied to the object.
    //
    // See [Johnson, 1985] page 353
    static Real HertzImpactTimeScale( const MeshSet &objectA,
                                      const MeshSet &objectB,
                                      const TwoObjectImpact &impactData,
                                      const VEC3F &impulseDirectionA,
                                      const VEC3F &impulseDirectionB,
                                      Real &inverseEffectiveMass );

    // Hertz time scale for an impact with an ideally flat surface (eg.
    // a ground plane)
    static Real HertzImpactTimeScale( const MeshSet &objectA,
                                      const PlaneImpact &impactData,
                                      const VEC3F &impulseDirectionA,
                                      Real &inverseEffectiveMass );
                                      

    // Computes the 1 / E^{*} quantity from [Johnson, 1985] page 352
    static Real HertzMaterialConstant( Real youngsModulusA,
                                       Real youngsModulusB,
                                       Real poissonRatioA,
                                       Real poissonRatioB );

    // Adds impact noise for a single object
    static bool AddObjectImpactNoise(
                          Real startTime, Real timeScale, Real forceScale,
                          const VEC3F &translationalAcceleration,
                          const VEC3F &rotationalAcceleration,
                          const VEC3F &listeningPosition,
                          RadialApproximation::AccelerationSet &pulseData,
                          SampledFunction &outputSignal,
                          Real scale = 1.0 );

    // Adds impact noise for a single object
    static bool AddObjectImpactNoise(
                          Real startTime, Real timeScale, Real forceScale,
                          const VEC3F &translationalAcceleration,
                          const VEC3F &rotationalAcceleration,
                          const VEC3F &listeningPosition,
                          CompactRadialApproximation *pulseData,
                          SampledFunction &outputSignal,
                          Real scale = 1.0 );

    static Real ImpactForceScale( Real inverseEffectiveMass,
                                  Real hertzTimeScale,
                                  Real impulseScale )
    {
#if 0
      return ( M_PI * impulseScale )
             / ( 2.0 * hertzTimeScale * inverseEffectiveMass );
#endif
      // Match the time-averaged forces
      return ( M_PI * impulseScale ) / ( 2.0 * hertzTimeScale );
    }

    static VEC3F ObjectListeningPosition(
                                const VEC3F &listeningPosition,
                                const VEC3F &restCenterOfMass,
                                const VEC3F &impactCenterOfMass,
                                const QUATERNION &impactInverseRotation )
    {
      VEC3F                  centerOfMassOffset;

      // Offset from the current center of mass
      centerOfMassOffset = listeningPosition - impactCenterOfMass;

      // Rotate out of the impact frame
      centerOfMassOffset = impactInverseRotation.rotate( centerOfMassOffset );

      // Add back the original center of mass
      centerOfMassOffset += restCenterOfMass;

      return centerOfMassOffset;
    }

};

#endif
