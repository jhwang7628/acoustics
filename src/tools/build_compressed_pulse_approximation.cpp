//////////////////////////////////////////////////////////////////////
// build_compressed_approximation_system.cpp: Builds a multi term
//                                            approximation
//                                            represented with
//                                            wavelets
//
//////////////////////////////////////////////////////////////////////

#include <TYPES.h>

#include <distancefield/distanceField.h>

#include <math/InterpolationFunction.h>

#include <geometry/RigidMesh.h>
#include <geometry/TriangleMesh.hpp>

#include <parser/Parser.h>

#include <transfer/CompressedMultiTermApproximation.h>
#include <transfer/MultiTermApproximation.h>

#include <utils/IO.h>
#include <utils/MathUtil.h>

#include <iostream>
#include <string>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Build a multi-term approximation
//////////////////////////////////////////////////////////////////////
void buildPulseApproximationSystem( REAL fieldRadius,
                                    const Point3d &fieldCenter,
                                    const Parser::AccelerationPulseParms &parms,
                                    const char *suffix )
{
    RadialApproximation::AccelerationSet   pulseModel;
    FloatArray                             radii;
    InterpolationMitchellNetravali        *interp;

    string                                 dataFile;
    string                                 compressedOutputFile;

    int                                    listeningResolution;

    if ( parms._symmetrizeField ) {
        dataFile = parms._outputFile + suffix;
    }

    if ( parms._symmetrizeField ) {
        listeningResolution = parms._symmetricResolution;
        dataFile = parms._outputFile + suffix;
        compressedOutputFile = parms._compressedOutputFile + suffix;
    } else {
        listeningResolution = parms._listeningResolution;
        dataFile = parms._outputFile;
        compressedOutputFile = parms._compressedOutputFile;
    }

    // Build the listening radii set
    radii.push_back( fieldRadius );
    for ( int term_idx = 1; term_idx < parms._numShells; term_idx++ )
    {
        fieldRadius *= parms._radiusMultiplier;
        radii.push_back( fieldRadius );
    }

    interp = new InterpolationMitchellNetravali( parms._pulseTimeScale );

    for ( int accel_direction = TRANS_X; accel_direction < NUM_ACCEL_DIRECTIONS;
            accel_direction++ )
    {
        pulseModel._allFields[ accel_direction ] = new CompressedMultiTermApproximation(
                                                                    radii, parms._numTerms,
                                                                    dataFile,
                                                                    fieldCenter,
                                                                    listeningResolution,
                                                                    interp->h(),
                                                                    interp->supportLength(),
                                                                    parms._timeStepFrequency,
                                                                    accel_direction,
                                                                    parms._compressionTolerance );
    }

    MultiTermApproximation::WriteAccelerationSet( compressedOutputFile,
            pulseModel );

    delete interp;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main( int argc, char **argv )
{
    string                   fileName( "default.xml" );
    Parser                  *parser = NULL;
    TriangleMesh<REAL>      *mesh = NULL;
    RigidMesh               *rigidMesh = NULL;
    REAL                     listeningRadius;

    const char              *suffix = "_octant";

    Parser::AccelerationPulseParms parms;

    if ( argc >= 2 )
    {
        suffix = argv[ 1 ];
    }

    if ( argc >= 3 )
    {
        fileName = argv[ 2 ];
    }

    parser = Parser::buildParser( fileName );

    if ( !parser )
    {
        cerr << "ERROR: Could not build parser from " << fileName << endl;
        return 1;
    }

    mesh = parser->getMesh();

    if ( !mesh )
    {
        cerr << "ERROR: Could not build mesh" << endl;
        return 1;
    }

    rigidMesh = new RigidMesh( *mesh, parms._rigidPrefix, parms._rigidDensity );

    parms = parser->getAccelerationPulseParms();

    listeningRadius = parms._listeningRadius;
    listeningRadius *= mesh->boundingSphereRadius( rigidMesh->centerOfMass() );

    buildPulseApproximationSystem( listeningRadius, rigidMesh->centerOfMass(), parms, suffix );

    delete parser;
}
