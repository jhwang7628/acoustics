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

#include <distancefield/closestPointField.h>
#include <distancefield/FieldBuilder.h>

#include <geometry/RigidMesh.h>
#include <geometry/TriangleMesh.hpp>

#include <deformable/ModeData.h>

#include <linearalgebra/Vector3.hpp>

#include <math/InterpolationFunction.h>

#include <parser/Parser.h>

#include <transfer/PulseApproximation.h>

#include <ui/WaveViewer.h>

#include <utils/IO.h>
#include <utils/MathUtil.h>

#ifdef USE_CUDA
    #include <wavesolver/gpusolver/wrapper/cuda/CUDA_PAN_WaveSolver.h>
    #include <wavesolver/gpusolver/wrapper/cuda/CUDA_PAT_WaveSolver.h>
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
    ClosestPointField     *sdf = NULL;
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

    Parser::AcousticTransferParms parms;

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

    parms = parser->getAcousticTransferParms();

    rigidMesh = new RigidMesh( *mesh, parms._rigidPrefix, parms._rigidDensity );

    CENTER_OF_MASS = rigidMesh->centerOfMass();

    sdf = DistanceFieldBuilder::BuildSignedClosestPointField( parser->getMeshFileName().c_str(),
                                                              parms._sdfResolution,
                                                              parms._sdfFilePrefix.c_str() );

    fieldBBox = BoundingBox( sdf->bmin(), sdf->bmax() );

    // Scale this up to build a finite difference field
    fieldBBox *= parms._gridScale;

    cellDivisions = parms._gridResolution;

    cellSize = min( fieldBBox.axislength( 0 ),
            min( fieldBBox.axislength( 1 ), fieldBBox.axislength( 2 ) ) );
    cout << SDUMP(cellSize) << endl;
    cellSize /= (REAL)cellDivisions;

    cout << SDUMP( cellSize ) << endl;

    timeStep = 1.0 / (REAL)( parms._timeStepFrequency);

    cout << SDUMP( timeStep ) << endl;
    cout << SDUMP( parms._outputFile ) << endl;

    //listeningRadius *= mesh->boundingSphereRadius( rigidMesh->centerOfMass() );
    REAL radius = mesh->boundingSphereRadius( rigidMesh->centerOfMass() );


    cout << SDUMP( radius ) << endl;
    cout << SDUMP( CENTER_OF_MASS ) << endl;

    radius *= parms._radiusMultipole;
 #ifdef USE_CUDA
    // Modes are cooler
    // cellSize = 0.5/(REAL)cellDivisions;
    // cout << SDUMP( cellSize ) << endl;
    // CUDA_PAT_WaveSolver solver(timeStep,
    //                fieldBBox, cellSize,
    //                *mesh, CENTER_OF_MASS,
    //                *sdf,
    //                0.0,
    //                &listeningPositions, //listeningPositions
    //                NULL,
    //                parms._subSteps,
    //                endTime,
    //                2*acos(-1)*4000,
    //                parms._nbar,
    //                parms._radiusMultipole*0.05);
    ModeData modedata;
    modedata.read(parms._modeDataFile.c_str());
    REAL density = parms._rigidDensity;
    CUDA_PAT_WaveSolver solver(timeStep,
                   fieldBBox, cellSize,
                   *mesh, CENTER_OF_MASS,
                   *sdf,
                   0.0,
                   parms._mode, //Mode
                   modedata,
                   density,
                   &listeningPositions, //listeningPositions
                   NULL,
                   parms._subSteps,
                   endTime,
                   parms._nbar,
                   radius,
                   50,
                   100000);
 #else
    PML_WaveSolver        solver( timeStep, fieldBBox, cellSize,
            *mesh, *sdf),
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

    WaveWindow                 waveApp( solver );

    waveApp.setAccelerationFunction( NULL );

    QWidget                   *window = waveApp.createWindow();

    window->show();

    return app.exec();
}
