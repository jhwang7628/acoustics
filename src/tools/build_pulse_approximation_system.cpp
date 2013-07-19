//////////////////////////////////////////////////////////////////////
// build_pulse_approximation_system.cpp: Builds the least squares
//                                       system for computing the
//                                       pulse approximation
//
//////////////////////////////////////////////////////////////////////

#include <TYPES.h>

#include <distancefield/distanceField.h>

#include <geometry/RigidMesh.h>
#include <geometry/TriangleMesh.hpp>

#include <math/InterpolationFunction.h>

#include <parser/Parser.h>

#include <transfer/MultiTermApproximation.h>
#include <transfer/PulseApproximation.h>

#include <utils/IO.h>
#include <utils/MathUtil.h>

#include <iostream>
#include <string>

using namespace std;

// FIXME
static REAL radiusMultiplier = 1.1892;
static int pointsPerShell;
static int numShells = 9;

//////////////////////////////////////////////////////////////////////
// Build a multi-term approximation
//////////////////////////////////////////////////////////////////////
void buildPulseApproximationSystem( REAL fieldRadius,
                                    const Point3d &fieldCenter,
                                    const Parser::AccelerationPulseParms &parms )
{
    RadialApproximation::AccelerationSet   pulseModel;
    FloatArray                             radii;
    InterpolationMitchellNetravali        *interp;
    REAL                                   svdTolerance = 0.005;
    //REAL                                   svdTolerance = 1e-1;

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
        pulseModel._allFields[ accel_direction ] = new MultiTermApproximation(
                                                                radii, parms._numTerms,
                                                                parms._outputFile,
                                                                fieldCenter,
                                                                parms._listeningResolution,
                                                                interp->h(),
                                                                interp->supportLength(),
                                                                parms._timeStepFrequency,
                                                                accel_direction,
                                                                svdTolerance );
    }

    MultiTermApproximation::WriteAccelerationSet( parms._multiTermOutputFile,
            pulseModel );

    delete interp;
}

#if 0
//////////////////////////////////////////////////////////////////////
// Build a multi-term approximation
//////////////////////////////////////////////////////////////////////
void buildPulseApproximationSystem( REAL fieldRadius,
                                    const Parser::AccelerationPulseParms &parms,
                                    const RadialApproximation::AccelerationSet &pulseModel )
{
    // FIXME: Let's try setting up the sparse system matrix here...
    REAL listeningRadius = fieldRadius;
    int startIndex_fixed = 0;
    int endIndex_fixed = 0;
    int startIndex, endIndex;
    MultiTermApproximation::GetRadiusIndexRange(
            listeningRadius,
            pulseModel._allFields[ 0 ]->signalLength(),
            parms._timeStepFrequency,
            startIndex, endIndex, 343.0 );
    startIndex_fixed = startIndex;
    endIndex_fixed = endIndex;

    for ( int i = 1; i < numShells; i++ )
    {
        listeningRadius *= radiusMultiplier;

        MultiTermApproximation::GetRadiusIndexRange(
                listeningRadius,
                pulseModel._allFields[ 0 ]->signalLength(),
                parms._timeStepFrequency,
                startIndex, endIndex, 343.0 );

        startIndex_fixed = min( startIndex, startIndex_fixed );
        endIndex_fixed = max( endIndex, endIndex_fixed );
    }

    SPARSE_MATRIX systemMatrix1(
            pulseModel._allFields[ 0 ]->signalLength() * numShells,
            endIndex_fixed - startIndex_fixed + 1 );
    SPARSE_MATRIX systemMatrix2(
            pulseModel._allFields[ 0 ]->signalLength() * numShells,
            endIndex_fixed - startIndex_fixed + 1 );

    listeningRadius = fieldRadius;
    FloatArray inputRadii;
    for ( int i = 0; i < numShells; i++ )
    {
        inputRadii.push_back( listeningRadius );
        MultiTermApproximation::AddRadiusTerms(
                listeningRadius,
                pulseModel._allFields[ 0 ]->signalLength(),
                1,
                parms._timeStepFrequency,
                systemMatrix1,
                startIndex_fixed,
                i * pulseModel._allFields[ 0 ]->signalLength(),
                343.0 );

        MultiTermApproximation::AddRadiusTerms(
                listeningRadius,
                pulseModel._allFields[ 0 ]->signalLength(),
                2,
                parms._timeStepFrequency,
                systemMatrix2,
                startIndex_fixed,
                i * pulseModel._allFields[ 0 ]->signalLength(),
                343.0 );

        listeningRadius *= radiusMultiplier;
    }

    cout << SDUMP( startIndex_fixed ) << endl;
    cout << SDUMP( endIndex_fixed ) << endl;

    MATRIX fullSystemMatrix;
    MultiTermApproximation::BuildMultiTermSystem(
            inputRadii,
            pulseModel._allFields[ 0 ]->signalLength(),
            2, /* 2 term model */
            parms._timeStepFrequency,
            fullSystemMatrix,
            343.0 );

    fullSystemMatrix.write( "testFullSystem.matrix" );

    // FIXME: test
    MATRIX tmp;
    tmp.read( "pulsedata_0.matrix" );
    MATRIX rhs;
    MultiTermApproximation::BuildMultiTermRHS( tmp, rhs, inputRadii.size() );

    rhs.write( "testRHS.matrix" );

    MATRIX Asub( 5 * pulseModel._allFields[ 0 ]->signalLength(),
            fullSystemMatrix.cols(), fullSystemMatrix.data() );
    MATRIX rhsSub( 5 * pulseModel._allFields[ 0 ]->signalLength(),
            rhs.cols(), rhs.data() );

    VECTOR singularValues;
    int rank;
    MATRIX::LeastSquares_TSVD( Asub, rhsSub, singularValues, rank, 0.002 );

    rhsSub.write( "testSolution.matrix" );

    systemMatrix1.writeToBinary( "test1.bcsm" );
    systemMatrix2.writeToBinary( "test2.bcsm" );
}
#endif

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main( int argc, char **argv )
{
    string                   fileName( "default.xml" );
    Parser                  *parser = NULL;
    TriangleMesh<REAL>      *mesh = NULL;
    RigidMesh               *rigidMesh = NULL;
    REAL                     listeningRadius;

    Parser::AccelerationPulseParms parms;

    if ( argc >= 2 )
    {
        fileName = argv[ 1 ];
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

    buildPulseApproximationSystem( listeningRadius, rigidMesh->centerOfMass(), parms );

#if 0
    // Read in pulse data from disk
    PulseApproximation::AccelerationSet  pulseModel;

    PulseApproximation::ReadAccelerationSet( parms._outputFile, pulseModel );

    buildPulseApproximationSystem( listeningRadius, parms, pulseModel );

    pulseModel.clear();
#endif

    delete parser;
}
