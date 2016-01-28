//////////////////////////////////////////////////////////////////////
// ImpulseIO.cpp: Implementation of the ImpulseIO class
//
//////////////////////////////////////////////////////////////////////

#include "ImpulseIO.h"

#include <utils/IO.h>
#include <utils/STLUtil.h>
#include <utils/trace.h>

#include <fstream>

using namespace std;

// Default material parameters for the ground
const REAL ImpulseIO::PLANE_YOUNGS_MODULUS  = 3.0e10;
const REAL ImpulseIO::PLANE_POISSON_RATIO   = 0.2;

MERSENNETWISTER ImpulseIO::GENERATOR;

//////////////////////////////////////////////////////////////////////
// Read impulses.  We assume that filename refers to a text file storing
// impulse information.  If binary == true, then we will attempt to
// read from a binary file named "filename".dat.  If the binary file
// does not exist, we read from the text file, and write to a binary
// format.
//////////////////////////////////////////////////////////////////////
ImpulseIO::ImpulseSet ImpulseIO::ReadImpulses( const char *filename,
                                               bool randomizeImpulses, REAL timeStep,
                                               REAL groundYoungsModulus, REAL groundPoissonRatio,
                                               bool binary )
{
    ImpulseSet               impulses;

    if ( binary ) {
        // Try to read/write impulses in a binary format
        if ( ReadBinaryImpulses( filename, impulses ) ) {
            return impulses;
        }
    }

    ReadTextImpulses( filename, impulses, groundYoungsModulus, groundPoissonRatio,
                      randomizeImpulses, timeStep );

    if ( binary ) {
        WriteBinaryImpulses( filename, impulses );
    }

    return impulses;
}

//////////////////////////////////////////////////////////////////////
// Reads impulses from a text file format
//////////////////////////////////////////////////////////////////////
void ImpulseIO::ReadTextImpulses( const char *filename, ImpulseSet &impulses,
                                  REAL planeYoungsModulus, REAL planePoissonRatio,
                                  bool randomizeImpulses, REAL timeStep )
{
    char                       impulseType;

    ifstream                   input( filename );

    ObjectImpulse              objectImpulse;
    PlaneImpulse               planeImpulse;

    planeImpulse._impactData._planeYoungsModulus = planeYoungsModulus;
    planeImpulse._impactData._planePoissonRatio = planePoissonRatio;

    int impulseNum = 0;

    TRACE_ASSERT( !input.fail() && input.good(), "Error reading text impulse file" );

    while ( !input.eof() ) {
        // Process the next impulse

        // Get the impulse type
        input >> impulseType;

        if ( impulseType == 'P' ) {
            ReadImpulse( input, objectImpulse, randomizeImpulses, timeStep );

            impulses.addImpulse( objectImpulse );
        } else if ( impulseType == 'C' ) {
            ReadImpulse( input, planeImpulse, randomizeImpulses, timeStep );

            impulses.addImpulse( planeImpulse );
        } else {
            cout << SDUMP( impulseType ) << endl;
            TRACE_ASSERT( NULL, "Should never get here" );
        }

        impulseNum++;
    }
}

//////////////////////////////////////////////////////////////////////
// Reads an object-object impulse
//////////////////////////////////////////////////////////////////////
void ImpulseIO::ReadImpulse( istream &input, ObjectImpulse &impulse,
                             bool randomizeImpulses, REAL timeStep )
{
    input >> impulse._impactData._impulseTime;

    if ( randomizeImpulses ) {
        impulse._impactData._impulseTime
            += GENERATOR.rand( timeStep ) - timeStep / 2.0;
    }

    // Object A information
    input >> impulse._objA_id;
    input >> impulse._impactData._posA[ 0 ];
    input >> impulse._impactData._posA[ 1 ];
    input >> impulse._impactData._posA[ 2 ];
    input >> impulse._impactData._centerOfMassA[ 0 ];
    input >> impulse._impactData._centerOfMassA[ 1 ];
    input >> impulse._impactData._centerOfMassA[ 2 ];
    input >> impulse._impactData._inverseRotationA.w;
    input >> impulse._impactData._inverseRotationA.v[ 0 ];
    input >> impulse._impactData._inverseRotationA.v[ 1 ];
    input >> impulse._impactData._inverseRotationA.v[ 2 ];

    // Object B information
    input >> impulse._objB_id;
    input >> impulse._impactData._posB[ 0 ];
    input >> impulse._impactData._posB[ 1 ];
    input >> impulse._impactData._posB[ 2 ];
    input >> impulse._impactData._centerOfMassB[ 0 ];
    input >> impulse._impactData._centerOfMassB[ 1 ];
    input >> impulse._impactData._centerOfMassB[ 2 ];
    input >> impulse._impactData._inverseRotationB.w;
    input >> impulse._impactData._inverseRotationB.v[ 0 ];
    input >> impulse._impactData._inverseRotationB.v[ 1 ];
    input >> impulse._impactData._inverseRotationB.v[ 2 ];

    // Other impulse information
    input >> impulse._impactData._relativeSpeed;
    input >> impulse._impactData._impulseMagnitude;
    input >> impulse._impactData._impulseDirection[ 0 ];
    input >> impulse._impactData._impulseDirection[ 1 ];
    input >> impulse._impactData._impulseDirection[ 2 ];
}

//////////////////////////////////////////////////////////////////////
// Reads an object-plane impulse
//////////////////////////////////////////////////////////////////////
void ImpulseIO::ReadImpulse( istream &input, PlaneImpulse &impulse,
                             bool randomizeImpulses, REAL timeStep )
{
    input >> impulse._impactData._impulseTime;

    if ( randomizeImpulses ) {
        impulse._impactData._impulseTime
            += GENERATOR.rand( timeStep ) - timeStep / 2.0;
    }

    // Object A information
    input >> impulse._objA_id;
    input >> impulse._impactData._posA[ 0 ];
    input >> impulse._impactData._posA[ 1 ];
    input >> impulse._impactData._posA[ 2 ];
    input >> impulse._impactData._centerOfMassA[ 0 ];
    input >> impulse._impactData._centerOfMassA[ 1 ];
    input >> impulse._impactData._centerOfMassA[ 2 ];
    input >> impulse._impactData._inverseRotationA.w;
    input >> impulse._impactData._inverseRotationA.v[ 0 ];
    input >> impulse._impactData._inverseRotationA.v[ 1 ];
    input >> impulse._impactData._inverseRotationA.v[ 2 ];

    // Other impulse information
    input >> impulse._impactData._relativeSpeed;
    input >> impulse._impactData._impulseMagnitude;
    input >> impulse._impactData._impulseDirection[ 0 ];
    input >> impulse._impactData._impulseDirection[ 1 ];
    input >> impulse._impactData._impulseDirection[ 2 ];
}

//////////////////////////////////////////////////////////////////////
// Reads impulses from a binary file - returns true iff successful
//
// Uses "filename".dat as the binary file name
//////////////////////////////////////////////////////////////////////
bool ImpulseIO::ReadBinaryImpulses( const char *filename, ImpulseSet &impulses )
{
    char                     binaryFilename[ 1024 ];
    FILE                    *f;

    //size_t                   bytes_read;

    sprintf( binaryFilename, "%s.dat", filename );

    f = fopen( binaryFilename, "rb" );
    if ( !f ) {
        cerr << "Could not load " << binaryFilename << endl;
        return false;
    }

    readVector( f, impulses._objectImpulses );
    readVector( f, impulses._planeImpulses );

    fread( (void*)&impulses._maxObjectImpulse, sizeof(REAL), 1, f );
    fread( (void*)&impulses._maxPlaneImpulse, sizeof(REAL), 1, f );

    fclose( f );

    return true;
}

//////////////////////////////////////////////////////////////////////
// Writes impulse data to a binary file
//////////////////////////////////////////////////////////////////////
void ImpulseIO::WriteBinaryImpulses( const char *filename, const ImpulseSet &impulses )
{
    char                     binaryFilename[ 1024 ];
    FILE                    *f;

    //size_t                   bytes_written;

    sprintf( binaryFilename, "%s.dat", filename );

    f = fopen( binaryFilename, "wb" );
    if ( !f ) {
        cerr << "Could not load " << binaryFilename << endl;
        return;
    }

    writeVector( f, impulses._objectImpulses );
    writeVector( f, impulses._planeImpulses );

    fwrite( (const void*)&impulses._maxObjectImpulse, sizeof(REAL), 1, f );
    fwrite( (const void*)&impulses._maxPlaneImpulse, sizeof(REAL), 1, f );

    fclose( f );
}
