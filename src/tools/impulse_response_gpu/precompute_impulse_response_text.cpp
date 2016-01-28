/* 
 * CUDA port of file tools/acceleration_noise/precompute_acceleration_pulse.cpp
 * using GPU enabled wavesolver
 */

#include <config.h>
#include <TYPES.h>

#include <distancefield/closestPointField.h>
#include <distancefield/FieldBuilder.h>

//#include <geometry/RigidMesh.h>
#include <geometry/TriangleMesh.hpp>

//#include <deformable/ModeData.h>

#include <linearalgebra/Vector3.hpp>

#include <parser/Parser.h>

//#include <transfer/PulseApproximation.h>


//#include <utils/IO.h>
//#include <utils/MathUtil.h>
//#include <utils/Evaluator.h>

//#include <wavesolver/gpusolver/wrapper/cuda/CUDA_PAN_WaveSolver.h>
#include <wavesolver/gpusolver/wrapper/cuda/CUDA_PAT_WaveSolver.h>

#include <wavesolver/WaveSolver.h>

#include <boost/bind.hpp>

#include <iostream>
#include <string>



#if 1  
    #define _GNU_SOURCE 1  
    #include <fenv.h>
    static void __attribute__ ((constructor))
    trapfpe ()
      {
        /* Enable some exceptions.  At startup all exceptions are
         * masked.  */
        feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);  
      }
#endif  


void readConfigFile(const char * filename, REAL * endTime, 
                    Vector3Array * listeningPositions, 
                    char * pattern, 
                    Vector3d & sound_source); 

void writeData(const vector<vector<FloatArray> > & w, char * pattern, REAL endTime, int n); 

using namespace std;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main( int argc, char **argv )
{

#ifndef USE_CUDA
    printf("** ERROR: Not using CUDA\n");
    exit(1)
#endif


#if 0 //complete refactor configuration reader

    string                   fileName( "default.xml" );
    Parser                   *parser = NULL;
    Parser::ImpulseResponseParms parms;

    parser = Parser::buildParser( fileName );

    parms = parser->getImpulseResponseParms();

    REAL * tmp = (parms._sourcePositions)[0];
    cout << tmp[1] << endl;

#else

    if (argc != 3){
        printf("**Usage: %s solver_configuration_XML listening_configuration_file");
        exit(1);
    }

    string                   fileName( "default.xml" );
    Parser                  *parser = NULL;
    TriangleMesh<REAL>      *mesh = NULL;
    RigidMesh               *rigidMesh = NULL;
    ClosestPointField       *sdf = NULL;
    BoundingBox              fieldBBox;
    REAL                     cellSize;
    int                      cellDivisions;

    REAL                     endTime = -1.0;

    //REAL                     listeningRadius;
    Vector3Array             listeningPositions;
    Vector3d                 sound_source;
    //BoundaryEvaluator        boundaryCondition;

    REAL                     timeStep;

    Parser::AcousticTransferParms parms;

    fileName = argv[1];
    char pattern[100];
    readConfigFile(argv[2], &endTime, &listeningPositions, pattern, sound_source);
    printf("Queried end time for the simulation = %lf\n", endTime);

    /* Build parser from solver configuration. */
    parser = Parser::buildParser( fileName );
    if ( !parser )
    {
        cerr << "ERROR: Could not build parser from " << fileName << endl;
        return 1;
    }

    /* Build mesh. */
    mesh = parser->getMesh();
    if ( !mesh )
    {
        cerr << "ERROR: Could not build mesh" << endl;
        return 1;
    }

    parms = parser->getAcousticTransferParms();

    Vector3d CENTER_OF_MASS( 1.0, 0.0, 0.0 );

    sdf = DistanceFieldBuilder::BuildSignedClosestPointField( parser->getMeshFileName().c_str(),
                                                              parms._sdfResolution,
                                                              parms._sdfFilePrefix.c_str() );

    fieldBBox = BoundingBox( sdf->bmin(), sdf->bmax() );

    cout << "max, min = " << sdf->bmin() <<  ", " << sdf->bmax() << endl;

    // Scale this up to build a finite difference field
    fieldBBox *= parms._gridScale;

    cellDivisions = parms._gridResolution;

    cellSize = fieldBBox.minlength(); // use min length to check the CFL = c0*dt/dx
    cout << SDUMP(cellSize) << endl;
    cellSize /= (REAL)cellDivisions;

    cout << SDUMP( cellSize ) << endl;

    timeStep = 1.0 / (REAL)( parms._timeStepFrequency);


    REAL boundRadius = mesh->boundingSphereRadius( CENTER_OF_MASS );

    cout << SDUMP( timeStep ) << endl;
    cout << SDUMP( parms._outputFile ) << endl;
    cout << SDUMP( boundRadius ) << endl;
    cout << SDUMP( CENTER_OF_MASS ) << endl;
    cout << SDUMP( cellSize ) << endl; 

    boundRadius *= parms._radiusMultipole;

    /* Callback function for logging pressure at each time step. */
    WaveSolver::WriteCallback dacallback = boost::bind(writeData, _1, pattern, endTime, listeningPositions.size());

    //CUDA_PAT_WaveSolver solver(timeStep,
    //                           sound_source,
    //                           fieldBBox, cellSize,
    //                           *mesh, 
    //                           CENTER_OF_MASS,
    //                           *sdf,
    //                           0.0,
    //                           &listeningPositions,
    //                           &dacallback,
    //                           parms._subSteps,
    //                           endTime,
    //                           2*acos(-1)*4000,
    //                           parms._nbar,
    //                           boundRadius,
    //                           40,
    //                           100000);


    CUDA_PAT_WaveSolver solver(timeStep,
                               sound_source,
                               fieldBBox, cellSize,
                               *mesh, 
                               CENTER_OF_MASS,
                               *sdf,
                               0.0, // distance tolerance
                               &listeningPositions,
                               &dacallback,
                               parms._subSteps,
                               endTime,
                               QNAN_R, 
                               QNAN_I,
                               QNAN_R, 
                               40,
                               100000);



    bool ended = false;

    while(!ended){
        ended = !solver.stepSystem(NULL);
        REAL time = solver.currentSimTime();
        printf("%lf out of %lf (%lf %%)\n", time, endTime, time/endTime);
    }

    solver.writeWaveOutput();


#endif

    return 0;

}


/* 
 * Reading listening positions. 
 *
 * in: 
 *      filename: file name of read file. 
 *      endTime: end time of the simulation 
 *      listeningPosition: vector of R3 listening positions. 
 *      pattern: pattern of the written file path/name. 
 *      sound_source: sound source position
 */
void readConfigFile(const char * filename, REAL * endTime, Vector3Array * listeningPositions, char * pattern, Vector3d & sound_source)
{
    FILE * fp;
    fp = fopen(filename, "r+");
    if(fp == NULL){
        printf("Config file does not exist\n");
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
        listeningPositions->push_back(Vector3d(x, y, z));
    }
}

/* 
 * Write listened data. 
 */
void writeData(const vector<vector<FloatArray> > & w, char * pattern, REAL endTime, int n)
{
    int samples = w[0][0].size();
    char buffer[80];
    char fname[100];
    REAL timestep = 0;

    if(samples > 1){
        timestep = endTime/(samples-1);
    }

    REAL dtime = 0;

    double mabs = 0;

    for (int s = 0; s < samples; s++)
    {
        for(int i = 0; i < n; i++)
            mabs = max(mabs, abs(w[i][0][s]));
    }


    for (int s = 0; s < samples; s++)
    {
        sprintf(buffer, "%05d", s);
        sprintf(fname, pattern, buffer);
        FILE * fp = fopen(fname, "w+");
        fprintf(fp, "%d %.6lf\n", n, dtime);

        for(int i = 0; i < n; i++)
            fprintf(fp, "%.10lf\n", w[i][0][s]/mabs);

        fclose(fp);
        dtime = dtime + timestep;
    }
}
