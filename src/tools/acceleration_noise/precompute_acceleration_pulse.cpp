//////////////////////////////////////////////////////////////////////
// precompute_acceleration_pulse.cpp: Precomputes the result of a
//                                    body accelerating over a short
//                                    time scale where the
//                                    acceleration pulse is modelled
//                                    using a simple interpolation
//                                    function
//
//////////////////////////////////////////////////////////////////////

#include <QtGui>

#include <config.h>
#include <TYPES.h>

#include <distancefield/distanceField.h>
#include <distancefield/FieldBuilder.h>

#include <geometry/RigidMesh.h>
#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/Vector3.hpp>

#include <math/InterpolationFunction.h>

#include <parser/Parser.h>

#include <transfer/PulseApproximation.h>

#include <ui/WaveViewer.h>

#include <utils/IO.h>
#include <utils/MathUtil.h>

#ifdef USE_CUDA
    #include <wavesolver/gpusolver/wrapper/cuda/CUDA_PAN_WaveSolver.h>
#else
    #include <wavesolver/PML_WaveSolver.h>
#endif

#include <wavesolver/WaveSolver.h>

#include <boost/bind.hpp>

#include <QGLViewer/qglviewer.h>

#include <GL/glut.h>

#include <iostream>
#include <string>

using namespace std;

// FIXME
static REAL radiusMultiplier = 1.1892;
static int pointsPerShell;
static int numShells = 9;

static Vector3d              ACCELERATION_DIRECTIONS[] = { Vector3d( 1.0, 0.0, 0.0 ),
                                                           Vector3d( 0.0, 1.0, 0.0 ),
                                                           Vector3d( 0.0, 0.0, 1.0 ),
                                                           Vector3d( 1.0, 0.0, 0.0 ),
                                                           Vector3d( 0.0, 1.0, 0.0 ),
                                                           Vector3d( 0.0, 0.0, 1.0 ) };

static Vector3d              CENTER_OF_MASS( 0.0, 0.0, 0.0 );

static REAL                  PULSE_HALF_PERIOD = 0.0;

bool                         ZERO_BC = true;

REAL boundaryEval( const Vector3d &x, const Vector3d &n, int obj_id, REAL t,
                   int field_id, InterpolationFunction *interp )
{
    REAL                       bcResult = 0.0;

    TRACE_ASSERT( obj_id == 0 );

    if ( t <= 2.0 * interp->supportLength() )
    {
        bcResult = interp->evaluate( t, interp->supportLength() );

        if ( field_id <= 2 )
        {
            bcResult *= n.dotProduct( ACCELERATION_DIRECTIONS[ field_id ] );
        }
        else
        {
#if 0
            bcResult *= ( n * ( cross( x - CENTER_OF_MASS,
                                ACCELERATION_DIRECTIONS[ field_id ] ) ) );
#endif
            bcResult *= n.dotProduct(
                            ( x - CENTER_OF_MASS ).crossProduct(
                                                        ACCELERATION_DIRECTIONS[ field_id ] ) );
        }

        if ( ZERO_BC )
        {
            ZERO_BC = false;
            cout << "Non-zero boundary condition!" << endl;
        }
    }

    return bcResult;
}

//////////////////////////////////////////////////////////////////////
// Callback for building the pulse approximation
//////////////////////////////////////////////////////////////////////
void buildPulseApproximation( const vector<vector<FloatArray> > &waveOutput,
        const Parser::AccelerationPulseParms &parms,
        const InterpolationFunction *interp,
        Vector3d fieldCenter, REAL fieldRadius )
{
#if 0
    PulseApproximation         pulseModel( waveOutput, fieldCenter, fieldRadius,
            parms._listeningResolution,
            interp->h(), interp->supportLength(),
            parms._timeStepFrequency,
            0 /* FIXME */ );

    pulseModel.write( parms._outputFile );
#endif
    RadialApproximation::AccelerationSet  pulseModel;

    for ( int accel_direction = 0; accel_direction < NUM_ACCEL_DIRECTIONS;
            accel_direction++ )
    {
        pulseModel._allFields[ accel_direction ] = new PulseApproximation(
                                                                    waveOutput, fieldCenter,
                                                                    fieldRadius,
                                                                    parms._listeningResolution,
                                                                    interp->h(),
                                                                    interp->supportLength(),
                                                                    parms._timeStepFrequency,
                                                                    accel_direction );
    }

    RadialApproximation::WriteAccelerationSet( parms._outputFile, pulseModel );

#if 0
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

    for ( int i = 0; i < 4; i++ )
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
    for ( int i = 0; i < numShells; i++ )
    {
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

    systemMatrix1.writeToBinary( "test1.bcsm" );
    systemMatrix2.writeToBinary( "test2.bcsm" );
#endif
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main( int argc, char **argv )
{
    QApplication             app( argc, argv );

    glutInit( &argc, argv );

    string                   fileName( "default.xml" );
    Parser                  *parser = NULL;
    TriangleMesh<REAL>      *mesh = NULL;
    RigidMesh               *rigidMesh = NULL;
    const DistanceField     *sdf = NULL;
    BoundingBox              fieldBBox;
    REAL                     cellSize;
    int                      cellDivisions;

    REAL                     endTime = -1.0;

    REAL                     listeningRadius;
    Vector3Array             listeningPositions;

    BoundaryEvaluator        boundaryCondition;

    REAL                     timeStep;

    WaveSolver::WriteCallback  callback;

    InterpolationFunction   *interp;

    Parser::AccelerationPulseParms parms;

    if ( argc >= 2 )
    {
        endTime = atof( argv[ 1 ] );
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

    parms = parser->getAccelerationPulseParms();

    rigidMesh = new RigidMesh( *mesh, parms._rigidPrefix, parms._rigidDensity );

    CENTER_OF_MASS = rigidMesh->centerOfMass();

#if 0
    sdf = mesh->buildDistanceField( parms._sdfResolution,
                                    parms._sdfFilePrefix.c_str() );
#endif
    sdf = DistanceFieldBuilder::BuildSignedDistanceField( parser->getMeshFileName().c_str(),
                                                          parms._sdfResolution,
                                                          parms._sdfFilePrefix.c_str() );

    fieldBBox = BoundingBox( sdf->bmin(), sdf->bmax() );

    // Scale this up to build a finite difference field
    fieldBBox *= parms._gridScale;

    cellDivisions = parms._gridResolution;

    cellSize = min( fieldBBox.axislength( 0 ),
            min( fieldBBox.axislength( 1 ), fieldBBox.axislength( 2 ) ) );
    cellSize /= (REAL)cellDivisions;

    cout << SDUMP( cellSize ) << endl;

    timeStep = 1.0 / (REAL)( parms._timeStepFrequency * parms._subSteps );

    cout << SDUMP( timeStep ) << endl;
    cout << SDUMP( parms._outputFile ) << endl;

    // Build a set of listening positions
    listeningRadius = parms._listeningRadius;
    listeningRadius *= mesh->boundingSphereRadius( rigidMesh->centerOfMass() );

    cout << SDUMP( rigidMesh->centerOfMass() );

    pointsPerShell = MathUtil::GenerateSpherePoints(
            rigidMesh->centerOfMass(), listeningRadius,
            parms._listeningResolution,
            listeningPositions );

    // FIXME: play around with something here..
    REAL newRadius = listeningRadius;
    cout << SDUMP( newRadius ) << endl;

    for ( int i = 1; i < parms._numShells; i++ )
    {
        newRadius *= parms._radiusMultiplier;
        cout << SDUMP( newRadius ) << endl;
        MathUtil::GenerateSpherePoints( rigidMesh->centerOfMass(), newRadius,
                parms._listeningResolution,
                listeningPositions );
    }

    // If we are using a fixed cell size, override it here
    if ( parms._fixedCellSize )
    {
        cout << "*** Using a fixed grid cell size ***" << endl;

        cellSize = parms._cellSize;

        // Figure out how big we need to make things to accomodate all
        // listening positions
        //
        // We add 30 here for the absorption layer
        cellDivisions = 2 * (int)ceil( newRadius / cellSize ) + 30;

        Vector3d minBound, maxBound;

        minBound = rigidMesh->centerOfMass();
        maxBound = rigidMesh->centerOfMass();
        for ( int i = 0; i < 3; i++ )
        {
            minBound[ i ] -= ( newRadius + 15.0 * (REAL)cellSize );
            maxBound[ i ] += ( newRadius + 15.0 * (REAL)cellSize );
        }

        fieldBBox.setMinBound( minBound );
        fieldBBox.setMaxBound( maxBound );
    }

    cout << SDUMP( pointsPerShell ) << endl;
    cout << SDUMP( listeningPositions.size() ) << endl;

#if 0
    // FIXME:
    for ( int i = 0; i < listeningPositions.size(); i++ )
    {
        cout << listeningPositions[ i ] << endl;
    }
#endif

    interp = new InterpolationMitchellNetravali( parms._pulseTimeScale );

    // Callback for writing out the pulse approximation
    callback = boost::bind( buildPulseApproximation, _1,
                            boost::ref( parms ), interp,
                            rigidMesh->centerOfMass(), listeningRadius );

#if 0
    WaveSolver                 solver( timeStep, fieldBBox, cellSize,
            *mesh, *sdf,
            true, /* Use leapfrog */
            0.0,
            true, /* Use boundary */
            true, /* Rasterize */
            &listeningPositions,
            //parms._outputFile.c_str(),
            NULL, /* No output file */
            &callback, /* Write callback*/
            parms._subSteps,
            6 /* accleration directions */ );
#endif
 #ifdef USE_CUDA
    CUDA_PAN_WaveSolver solver(timeStep,
                               fieldBBox, cellSize,
                               *mesh, CENTER_OF_MASS,
                               *sdf,
                               0.0,
                               &listeningPositions,
                               NULL,
                               parms._subSteps,
                               endTime,
                               parms._pulseTimeScale,
                               11,
                               1000000);
 #else
    PML_WaveSolver        solver( timeStep, fieldBBox, cellSize,
            *mesh, *sdf,
            0.0, /* distance tolerance */
            true, /* use boundary */
            &listeningPositions,
            NULL, /* No output file */
            &callback, /* Write callback */
            parms._subSteps,
            6, /* acceleration directions */
            endTime );

    solver.setPMLBoundaryWidth( 11.0, 1000000.0 );
 #endif

    boundaryCondition = boost::bind( boundaryEval, _1, _2, _3, _4, _5, interp );

    WaveWindow                 waveApp( solver );

    waveApp.setAccelerationFunction( &boundaryCondition );

    QWidget                   *window = waveApp.createWindow();

    window->show();

    return app.exec();
}
