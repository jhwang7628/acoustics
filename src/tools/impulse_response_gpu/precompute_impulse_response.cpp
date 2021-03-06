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
#include <geometry/SimplePointSet.hpp> 

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

void readConfigFile(const char * filename, REAL * endTime, Vector3Array * listeningPositions, char * pattern, Vector3d & sound_source){
    FILE * fp;
    fp = fopen(filename, "r+");
    if(fp == NULL){
        printf("Config file %s does not exist\n", filename);
        exit(1);
    }
    fscanf(fp, "%s", pattern);
    fscanf(fp, "%lf", endTime);
    fscanf(fp, "%lf %lf %lf", &(sound_source[0]), &(sound_source[1]), &(sound_source[2]));
    int n;
    fscanf(fp, "%d", &n);
    for(int i = 0; i < n; i++){
        double x, y, z;
        fscanf(fp, "%lf %lf %lf", &x, &y, &z);
        listeningPositions->push_back(Vector3d(x*-1, y, z*-1));
    }
    // printf("SIZE_LSIT: %d\n", listeningPositions->size());
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
    ClosestPointField     *sdf = NULL;
    BoundingBox              fieldBBox;
    REAL                     cellSize;
    int                      cellDivisions;

    REAL                     endTime = -1.0;

    Vector3d                 sound_source;

    REAL                     listeningRadius;
    Vector3Array             listeningPositions;

    BoundaryEvaluator        boundaryCondition;

    REAL                     timeStep;

    WaveSolver::WriteCallback  callback;

    InterpolationFunction   *interp;

    Parser::AcousticTransferParms parms;


    // Extra GL data // 
    ExtraGLData extraGLData; 
    SimplePointSet<REAL> * fluidMshCentroids; 
    SimplePointSet<REAL> * fluidMshCentroids2;



    char pattern[100];

    if (argc != 5){
        printf("Usage: %s wavesetup listenpos fluidmshcentroid fluidmshcentroid2\n", argv[0]);
        exit(1);
    }

    fileName = argv[1];
    cout << "argv[2] = " << argv[2] << endl;
        
    readConfigFile(argv[2], &endTime, &listeningPositions, pattern, sound_source);

    printf("%s\n", pattern);

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

    // rigidMesh = new RigidMesh( *mesh, parms._rigidPrefix, parms._rigidDensity );

    CENTER_OF_MASS = Vector3d( 0.0, 0.0, 0.0 );

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

    REAL radius = mesh->boundingSphereRadius( CENTER_OF_MASS );


    cout << SDUMP( radius ) << endl;
    cout << SDUMP( CENTER_OF_MASS ) << endl;

    fluidMshCentroids  = new SimplePointSet<REAL>(argv[3]); 
    fluidMshCentroids2 = new SimplePointSet<REAL>(argv[4]); 
    extraGLData.pointset1 = fluidMshCentroids; 
    extraGLData.pointset2 = fluidMshCentroids2; 

 #ifdef USE_CUDA
    // Modes are cooler
    // cellSize = 0.5/(REAL)cellDivisions;
    // cout << SDUMP( cellSize ) << endl;
    CUDA_PAT_WaveSolver solver(timeStep,
                               sound_source,
                   fieldBBox, cellSize,
                   *mesh, CENTER_OF_MASS,
                   *sdf,
                   0.0,
                   &listeningPositions, //listeningPositions
                   NULL,
                   parms._subSteps,
                   endTime,
                   2*acos(-1)*4000,
                   parms._nbar,
                   radius,
                   20,
                   100000);
    // ModeData modedata;
    // modedata.read(parms._modeDataFile.c_str());
    // REAL density = parms._rigidDensity;
    // CUDA_PAT_WaveSolver solver(timeStep,
    //                fieldBBox, cellSize,
    //                *mesh, CENTER_OF_MASS,
    //                *sdf,
    //                0.0,
    //                parms._mode, //Mode
    //                modedata,
    //                density,
    //                &listeningPositions, //listeningPositions
    //                NULL,
    //                parms._subSteps,
    //                endTime,
    //                parms._nbar,
    //                radius,
    //                50,
    //                100000);
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

    WaveWindow                 waveApp( solver , &extraGLData);

    waveApp.setAccelerationFunction( NULL );

    QWidget                   *window = waveApp.createWindow();

    window->show();

    return app.exec();
}
