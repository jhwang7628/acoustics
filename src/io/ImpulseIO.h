//////////////////////////////////////////////////////////////////////
// ImpulseIO.h: Interface for the ImpulseIO
//
//////////////////////////////////////////////////////////////////////

#ifndef IMPULSE_IO_H
#define IMPULSE_IO_H

#include <geometry/Point3.hpp>

#include <linearalgebra/Quaternion.hpp>

#include <utils/MERSENNETWISTER.h>

#include <TYPES.h>

#include <iostream>
#include <vector>

//////////////////////////////////////////////////////////////////////
// ImpulseIO class
//
// Data structures and IO routines for rigid-body impulse data
//////////////////////////////////////////////////////////////////////
class ImpulseIO {
    public:
        // Impulse structures

        //////////////////////////////////////////////////////////////////////
        // Physical description for an impact between two objects
        //////////////////////////////////////////////////////////////////////
        struct TwoObjectImpact {
            // Impact positions in each mesh, relative to the mesh's rest pose
            Point3d                      _posA;
            Point3d                      _posB;

            // Impulse information
            REAL                         _impulseTime;
            Vector3d                     _impulseDirection;
            REAL                         _impulseMagnitude;
            REAL                         _relativeSpeed;

            // Rigid state information for the two meshes
            Point3d                      _centerOfMassA;
            Point3d                      _centerOfMassB;

            Quat4d                       _inverseRotationA;
            Quat4d                       _inverseRotationB;

        };

        //////////////////////////////////////////////////////////////////////
        // Physical description for an impact between an object and a static
        // ground plane
        //////////////////////////////////////////////////////////////////////
        struct PlaneImpact {
            // Impact position in the mesh
            Point3d                      _posA;

            // Impulse information
            REAL                         _impulseTime;
            Vector3d                     _impulseDirection;
            REAL                         _impulseMagnitude;
            REAL                         _relativeSpeed;

            // Rigid state information for the mesh
            Point3d                      _centerOfMassA;
            Quat4d                       _inverseRotationA;

            // Also store material parameters for the ground plane here
            REAL                         _planeYoungsModulus;
            REAL                         _planePoissonRatio;

        };

        //////////////////////////////////////////////////////////////////////
        // For representing an impulse between two objects
        //////////////////////////////////////////////////////////////////////
        struct ObjectImpulse {
          int                                      _objA_id;
          int                                      _objB_id;

          TwoObjectImpact                          _impactData;
        };

        //////////////////////////////////////////////////////////////////////
        // For representing an impulse between an object and the ground
        // plane
        //////////////////////////////////////////////////////////////////////
        struct PlaneImpulse {
          int                                      _objA_id;

          PlaneImpact                              _impactData;
        };

        //////////////////////////////////////////////////////////////////////
        // Represents a full set of impulses
        //////////////////////////////////////////////////////////////////////
        struct ImpulseSet {
            ImpulseSet()
                : _maxObjectImpulse( 0.0 ),
                  _maxPlaneImpulse( 0.0 )
            {
            }

            void addImpulse( const ObjectImpulse &impulse )
            {
                _maxObjectImpulse = max( _maxObjectImpulse,
                                         impulse._impactData._impulseMagnitude );
            }

            void addImpulse( const PlaneImpulse &impulse )
            {
                _maxPlaneImpulse = max( _maxPlaneImpulse,
                                         impulse._impactData._impulseMagnitude );
            }

            vector<ObjectImpulse>                  _objectImpulses;
            vector<PlaneImpulse>                   _planeImpulses;

            REAL                                   _maxObjectImpulse;
            REAL                                   _maxPlaneImpulse;
        };

    public:
        // Read impulses.  We assume that filename refers to a text file storing
        // impulse information.  If binary == true, then we will attempt to
        // read from a binary file named "filename".dat.  If the binary file
        // does not exist, we read from the text file, and write to a binary
        // format.
        static ImpulseSet ReadImpulses( const char *filename,
                                        bool randomizeImpulses, REAL timeStep,
                                        REAL groundYoungsModulus = PLANE_YOUNGS_MODULUS,
                                        REAL groundPoissonRatio = PLANE_POISSON_RATIO,
                                        bool binary = true );

    private:
        // Reads impulses from a text file format
        //
        // Optionally set randomizeImpulses == true to distribute impulse
        // times uniformly randomly within each time step
        static void ReadTextImpulses( const char *filename, ImpulseSet &impulses,
                                      REAL planeYoungsModulus, REAL planePoissonRatio,
                                      bool randomizeImpulses, REAL timeStep );

        // Reads an object-object impulse
        static void ReadImpulse( istream &input, ObjectImpulse &impulse,
                                 bool randomizeImpulses, REAL timeStep );

        // Reads an object-plane impulse
        static void ReadImpulse( istream &input, PlaneImpulse &impulse,
                                 bool randomizeImpulses, REAL timeStep );

        // Reads impulses from a binary file - returns true iff successful
        //
        // Uses "filename".dat as the binary file name
        static bool ReadBinaryImpulses( const char *filename, ImpulseSet &impulses );

        // Writes impulse data to a binary file
        static void WriteBinaryImpulses( const char *filename, const ImpulseSet &impulses );

    private:
        ImpulseIO() {}
        ~ImpulseIO() {}

        static const REAL            PLANE_YOUNGS_MODULUS;
        static const REAL            PLANE_POISSON_RATIO;

        // Random number generator
        static MERSENNETWISTER       GENERATOR;

};

#endif
