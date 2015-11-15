/* 
 * CUDA port of file tools/acceleration_noise/precompute_acceleration_pulse.cpp
 * using GPU enabled wavesolver
 */


#include <config.h>
#include <TYPES.h>

#include "IO/IO.h"
#include <limits>

#include <math/InterpolationFunction.h>

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
//#include <wavesolver/gpusolver/wrapper/cuda/CUDA_PAT_WaveSolver.h>

#include <wavesolver/PML_WaveSolver.h>
#include <wavesolver/WaveSolver.h>
#include <wavesolver/WaveSolverPointData.h>

#include <boost/bind.hpp>

#include <iostream>
#include <string>


#include <unistd.h> 

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

/* Gaussian */
REAL Gaussian_3D( const Vector3d &evaluatePosition, const Vector3d &sourcePosition, const REAL stddev )
{
    //cout << "Initializing the field with Gaussian" << endl; 

    REAL value =       (evaluatePosition.x - sourcePosition.x) * (evaluatePosition.x - sourcePosition.x) 
                      +(evaluatePosition.y - sourcePosition.y) * (evaluatePosition.y - sourcePosition.y) 
                      +(evaluatePosition.z - sourcePosition.z) * (evaluatePosition.z - sourcePosition.z); 
    value = -value / ( stddev * stddev ); 
    value = exp( value ); 
    value *= 10000.0;
    return value; 
	
}


/* 
 * Harmonic Source with frequency w and velocity strength Sw, 
 * the pressure source strength will then be
 *
 * p = (-i*w) * rho * a^2/r * Sw * exp(-i*w*t)
 *
 * Modeled by a small ball radius given. (should be much smaller compared to
 * any length scale in the system)
 *
 */ 
REAL Harmonic_Source( const REAL &t, const Vector3d &sPosition, const Vector3d &tPosition, const REAL w, const REAL Sw, const REAL ballRadius, const REAL rho_0, const REAL phase)  
{

    const REAL r = (sPosition - tPosition).norm();  
    const REAL EPS = 1E-10;

    if ( r <= ballRadius && r > EPS)  // prevent singularity
    { 
        return -rho_0 *ballRadius *ballRadius / 4.0 / PI / r  * Sw * w * cos(w*t + phase);
    }
    else 
    {
        return numeric_limits<REAL>::quiet_NaN();
    }
}


void readConfigFile(const char * filename, REAL * endTime, 
                    Vector3Array * listeningPositions, 
                    char * pattern, 
                    Vector3d & sound_source); 

void writeData(const vector<REAL> & w, const REAL & timeStep, const int & timeStamp, const int & substep, char * pattern, REAL endTime, int n);

REAL boundaryEval( const Vector3d &x, const Vector3d &n, int obj_id, REAL t,
                   int field_id, InterpolationFunction *interp )
{
    REAL                       bcResult = 0.0;

    TRACE_ASSERT( obj_id == 0 );

    //if ( t <= 2.0 * interp->supportLength() )
    //{
    //    bcResult = interp->evaluate( t, interp->supportLength() );

    //    if ( field_id <= 2 )
    //    {
    //        bcResult *= n.dotProduct( ACCELERATION_DIRECTIONS[ field_id ] );
    //    }
    //    else
    //    {
    //        bcResult *= n.dotProduct(
    //                        ( x - CENTER_OF_MASS ).crossProduct(
    //                                                    ACCELERATION_DIRECTIONS[ field_id ] ) );
    //    }

    //    if ( ZERO_BC )
    //    {
    //        ZERO_BC = false;
    //        cout << "Non-zero boundary condition!" << endl;
    //    }
    //}

    return bcResult;
}

using namespace std;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main( int argc, char **argv )
{




    if (argc != 3){
        printf("**Usage: %s solver_configuration_XML listening_configuration_file", argv[0]);
        exit(1);
    }

    string                   fileName( "default.xml" );
    Parser                  *parser = NULL;
    TriangleMesh<REAL>      *mesh = NULL;
    //RigidMesh               *rigidMesh = NULL;
    ClosestPointField       *sdf = NULL;
    BoundingBox              fieldBBox;
    REAL                     cellSize;
    int                      cellDivisions;

    REAL                     endTime = -1.0;

    //REAL                     listeningRadius;
    Vector3Array             listeningPositions;
    Vector3d                 sourcePosition;
    BoundaryEvaluator        boundaryCondition;

    REAL                     timeStep;

    Parser::AcousticTransferParms parms;

    fileName = argv[1];
    char pattern[100];
    readConfigFile(argv[2], &endTime, &listeningPositions, pattern, sourcePosition);
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


    // Scale this up to build a finite difference field
    fieldBBox *= parms._gridScale;
    cout << "fieldBBox -> [min, max] = " << fieldBBox.minBound() <<  ", " << fieldBBox.maxBound() << endl;

    cellDivisions = parms._gridResolution;

    cellSize = fieldBBox.minlength(); // use min length to check the CFL = c0*dt/dx
    cellSize /= (REAL)cellDivisions;

    cout << SDUMP( cellSize ) << endl;

    timeStep = 1.0 / (REAL)( parms._timeStepFrequency);


    REAL boundRadius = mesh->boundingSphereRadius( CENTER_OF_MASS );

    cout << SDUMP( timeStep ) << endl;
    cout << SDUMP( parms._outputFile ) << endl;
    cout << SDUMP( boundRadius ) << endl;
    cout << SDUMP( CENTER_OF_MASS ) << endl;
    cout << SDUMP( cellSize ) << endl; 


    REAL CFL = parms._c * timeStep / cellSize; 
    cout << SDUMP( CFL ) << endl;

    boundRadius *= parms._radiusMultipole;

    /* Callback function for logging pressure at each time step. */
    PML_WaveSolver::WriteCallbackIndividual dacallback = boost::bind(writeData, _1, _2, _3, parms._subSteps, pattern, endTime, listeningPositions.size());

    // TODO bind ic_eval




   

    /// pass in listening position
    PML_WaveSolver        solver( timeStep, fieldBBox, cellSize,
                                  *mesh, *sdf,
                                  0.0, /* distance tolerance */
                                  true, /* use boundary */
                                  &listeningPositions,
                                  NULL, /* No output file */
                                  NULL, /* Write callback */
                                  NULL, //&dacallback, /* Write callback */
                                  parms._subSteps,
                                  1, /* acceleration directions */
                                  endTime );


    //WaveSolverPointData   rawData; 
    //PML_WaveSolver        solver( timeStep, fieldBBox, cellSize,
    //                              *mesh, *sdf,
    //                              &rawData,
    //                              0.0, /* distance tolerance */
    //                              true, /* use boundary */
    //                              &listeningPositions,
    //                              NULL, /* No output file */
    //                              NULL, /* Write callback */
    //                              parms._subSteps,
    //                              1, /* acceleration directions */
    //                              endTime );


    const REAL ic_stddev = 0.001; 
    //const Vector3d sourcePosition(sound_source); 
    cout << "Source position is set at " << sourcePosition << endl;
    //InitialConditionEvaluator initial = boost::bind(Gaussian_3D, _1, _2, ic_stddev); 
    InitialConditionEvaluator initial = boost::bind(Gaussian_3D, _1, sourcePosition, ic_stddev);
    solver.initSystemNontrivial( 0.0, &initial ); 
    //solver.setPMLBoundaryWidth( 11.0, 1000000.0 );
    solver.setPMLBoundaryWidth( 20.0, 100000.0 );

   
#if 0 // if harmonic sources 
    const REAL w = 2.0*3.1415926*2000; // 500 Hz test
    const REAL Sw = 1.0; // source strength
    const REAL ballRadius = cellSize*3; // threshold for discretizing point source position
    cout << SDUMP(ballRadius) << endl;
    const REAL phase = 0;
    HarmonicSourceEvaluator sourceFunction = boost::bind( Harmonic_Source, _1, sourcePosition, _2, w, Sw, ballRadius, parms._density, phase ); 
    solver.setHarmonicSource( &sourceFunction );  // WILL NEGATE ALL INITIALIZATION! 
#endif 


    InterpolationFunction * interp = new InterpolationMitchellNetravali( 0.1 );
    boundaryCondition = boost::bind( boundaryEval, _1, _2, _3, _4, _5, interp ); 


    


    /// write grid and pressure value
    Eigen::Vector3i minDrawBound; 
    Eigen::Vector3i maxDrawBound; 
    minDrawBound << 0, 0, 0; 
    maxDrawBound << solver.fieldDivisions()[0], solver.fieldDivisions()[1], solver.fieldDivisions()[2]; 
    printf( "domain size: (%u, %u, %u)\n", maxDrawBound[0], maxDrawBound[1], maxDrawBound[2] );

    int Ndivision = maxDrawBound[0] * maxDrawBound[1] * maxDrawBound[2];
    Eigen::MatrixXd vertexPosition( Ndivision, 3 ); // vertex position
    Eigen::MatrixXd vertexPressure( Ndivision, 1 ); // full pressure
    Eigen::MatrixXi vertexIndex( Ndivision, 3 ); // pressure vertex index


    int count = 0;
    for ( int kk=0; kk<maxDrawBound[2]; kk++ )
    {
        for ( int jj=0; jj<maxDrawBound[1]; jj++ ) 
        {
            for ( int ii=0; ii<maxDrawBound[0]; ii++ ) 
            {
                Tuple3i  vIndex( ii, jj, kk );
                Vector3d vPosition = solver.fieldPosition(  vIndex );

                vertexPosition( count, 0 ) = vPosition[0];
                vertexPosition( count, 1 ) = vPosition[1];
                vertexPosition( count, 2 ) = vPosition[2];

                vertexIndex( count, 0 ) = vIndex[0]; 
                vertexIndex( count, 1 ) = vIndex[1]; 
                vertexIndex( count, 2 ) = vIndex[2]; 

                count ++; 

            }
        }
    }

    printf("writing pressure vertex position\n");
    {

        char buffer[100]; 
        snprintf( buffer,100,pattern,"vertex_position.dat"); 
        IO::writeMatrixXd( vertexPosition, buffer, IO::BINARY );
        //ofstream of(buffer); 
        //of << std::setprecision(16) << std::fixed; 
        //of << vertexPosition << endl; 
        //of.close();
        //ofstream of(buffer); 
        //of << std::setprecision(16) << std::fixed; 
        //of << vertexPosition << endl; 
        //of.close();

    }

    printf("writing pressure vertex index\n");
    {

        char buffer[100]; 
        snprintf( buffer,100,pattern,"vertex_index.dat"); 
        IO::writeMatrixXi( vertexIndex, buffer, IO::BINARY );
        //ofstream of(buffer); 
        //of << std::setprecision(16) << std::fixed; 
        //of << vertexIndex << endl; 
        //of.close();

    }

    printf("writing solver settings\n"); 
    {

        char buffer[100]; 
        snprintf( buffer, 100, pattern, "solver_setting.txt"); 
        ofstream of(buffer); 

        of << "[fieldBBox.minBound] = " << fieldBBox.minBound() << endl;
        of << "[fieldBBox.maxBound] = " << fieldBBox.maxBound() << endl;
        of << "[solver.fieldDivisions] = " << maxDrawBound.transpose() << endl;
        of << SDUMP( cellSize ) << endl;
        of << SDUMP( CFL ) << endl;
        of << SDUMP( timeStep ) << endl;
        of << SDUMP( endTime ) << endl;
  
        of << "precision : double (hard-coded)" << endl;
        of << "binary write method : IO writeMatrixXd method (hard-coded)" << endl;
        of.close();

    }



    /// time step system
    bool continueStepping = true; 
    int nSteps = 0; 
    while ( continueStepping )
    {
        //continueStepping = solver.stepSystem( boundaryCondition, &sourceFunction ); 
        continueStepping = solver.stepSystem( boundaryCondition ); 

        if ( nSteps % parms._subSteps == 0 ) 
        {
            int count = 0;
            for ( int kk=0; kk<maxDrawBound[2]; kk++ )
            {
                for ( int jj=0; jj<maxDrawBound[1]; jj++ ) 
                {
                    for ( int ii=0; ii<maxDrawBound[0]; ii++ ) 
                    {
                        Tuple3i  vIndex( ii, jj, kk );

                        VECTOR vPressure;
                        solver.vertexPressure( vIndex, vPressure );

                        vertexPressure( count, 0 ) = vPressure[0];

                        count ++; 

                    }
                }
            }


            printf( "writing pressure %u\n", nSteps/parms._subSteps );
            {
                char buf[100]; 
                char pressure_buf[100]; 
                snprintf( pressure_buf, 100, "pressure_%u.dat", nSteps/parms._subSteps ); 
                snprintf( buf, 100, pattern, pressure_buf );

                IO::writeMatrixXd( vertexPressure, buf, IO::BINARY );
                //ofstream of(buf); 
                //of << std::setprecision(16) << std::fixed; 
                //of << vertexPressure << endl; 
                //of.close();
            }
        }


    nSteps ++; 


    }

    cout << "End of program. " << endl;



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
void writeData(const vector<REAL> & w, const REAL & timeStep, const int & timeStamp, const int & substep, char * pattern, REAL endTime, int n)
{
    char buffer[80];
    char fname[100];

    REAL nowTime = timeStep * (REAL) timeStamp; 
    sprintf(buffer, "%05d", timeStamp/substep);
    sprintf(fname, pattern, buffer);

    FILE * fp = fopen(fname, "w+");
    fprintf(fp, "%d %.6lf\n", n, nowTime);

    for (size_t ii=0; ii<w.size(); ii++)
        fprintf( fp, "%.16f\n", w[ii] );

    fclose(fp);
}



