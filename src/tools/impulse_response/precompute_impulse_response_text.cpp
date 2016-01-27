/* 
 * impulse-response solve using the finite-difference time-domain (FDTD) wave solver. 
 *
 * ref: 
 *  [1] Chadwick 2012, Precomputed Acceleration Noise for Improved Rigid-Body Sound
 *  [2] Liu 1997, The perfectly matched layer for acoustic waves in absorptive media
 *
 * Author: Jui-Hsien Wang
 *  
 *  
 */


#include <config.h>
#include <TYPES.h>

#include "IO/IO.h"
#include <limits>

#include <distancefield/closestPointField.h>
#include <distancefield/FieldBuilder.h>

#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/Vector3.hpp>

#include <parser/Parser.h>

#include <wavesolver/PML_WaveSolver.h>
#include <wavesolver/WaveSolver.h>

#include <grid/Grid.h>

#include <boost/bind.hpp>

#include <iostream>
#include <string>

#include <unistd.h> 

#include "vtkConverter/vtkConverter.h" 

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510582
#endif

//const double sqrt_pi_over_2 = sqrt(M_PI/2.0); 
const double sqrt_2_pi = sqrt(2.0*M_PI); 
const double sqrt_2 = sqrt(2.0); 

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

template <class T> 
Eigen::Matrix<T,3,1> ToEigenVector3(const Vector3<T> &vec)
{
    return Eigen::Matrix<T,3,1>(vec.x,vec.y,vec.z); 
}

/* Gaussian */
REAL Gaussian_3D( const Vector3d &evaluatePosition, const Vector3d &sourcePosition, const REAL stddev )
{
    //cout << "Initializing the field with Gaussian" << endl; 

    REAL value =       (evaluatePosition.x - sourcePosition.x) * (evaluatePosition.x - sourcePosition.x) 
                      +(evaluatePosition.y - sourcePosition.y) * (evaluatePosition.y - sourcePosition.y) 
                      +(evaluatePosition.z - sourcePosition.z) * (evaluatePosition.z - sourcePosition.z); 
    value = -value / ( 2.0 * stddev * stddev ); 
    value = exp( value ); 
    value *= 10000.0;
    return value; 
	
}

/* space-time gaussian with center at sourcePosition x offsetTime */
REAL Gaussian_4D( const Vector3d &evaluatePosition, const REAL &evaluateTime, const Vector3d &sourcePosition, const REAL &widthSpace, const REAL &widthTime, const REAL &offsetTime ) 
{
    REAL value =       (evaluatePosition.x - sourcePosition.x) * (evaluatePosition.x - sourcePosition.x) 
                      +(evaluatePosition.y - sourcePosition.y) * (evaluatePosition.y - sourcePosition.y) 
                      +(evaluatePosition.z - sourcePosition.z) * (evaluatePosition.z - sourcePosition.z); 
    value  = -value / ( 2.0 * widthSpace * widthSpace ); 
    value += -(evaluateTime - offsetTime)*(evaluateTime - offsetTime)/(2.0*widthTime*widthTime); 
    value  = exp( value ); 
    value *= 10000.0; // just to scale it up, the wave equation is linear

    return value; 
}

/* 
 * 3D gaussian in space x error function in time (integral of shifted gaussian) 
 *
 * */ 
REAL Gaussian_3D_erf_time( const Vector3d &evaluatePosition, const REAL &evaluateTime, const Vector3d &sourcePosition, const REAL &widthSpace, const REAL &widthTime, const REAL &offsetTime, const REAL &normalizeConstant ) 
{
    REAL value = (evaluatePosition.x - sourcePosition.x)*(evaluatePosition.x - sourcePosition.x) 
               + (evaluatePosition.y - sourcePosition.y)*(evaluatePosition.y - sourcePosition.y) 
               + (evaluatePosition.z - sourcePosition.z)*(evaluatePosition.z - sourcePosition.z); 
    value  = -value / (2.0*widthSpace*widthSpace); 
    value  = exp(value); 

    // integral of a gaussian in time: 
    // http://www.wolframalpha.com/input/?i=int%28+exp%28-%28x-a%29%5E2%2F2%2Fb%5E2+%29+%29+dx
    value *= erf((evaluateTime - offsetTime)/(sqrt_2*widthTime)) - erf(-offsetTime/(sqrt_2*widthTime)); 

    value *= normalizeConstant; // just to scale it up, the wave equation is linear

    return value; 
}

/* Evaluate distance between source and points in the grid for Point Gaussian initialization */
REAL PointGaussian_3D( const Vector3d &evaluatePosition, const Vector3d &sourcePosition )
{
    //cout << "Initializing the field with Gaussian" << endl; 

    REAL value =       (evaluatePosition.x - sourcePosition.x) * (evaluatePosition.x - sourcePosition.x) 
                      +(evaluatePosition.y - sourcePosition.y) * (evaluatePosition.y - sourcePosition.y) 
                      +(evaluatePosition.z - sourcePosition.z) * (evaluatePosition.z - sourcePosition.z); 

    return value; 
    
}

/* readjust the source position so that it lies at the center of some cell */ 
void AdjustSourcePosition( const Vector3d &minBound, const Vector3d &maxBound, const int &divisions, Vector3d &sourcePosition ) 
{
    Vector3d dx = (maxBound-minBound)/(double)divisions; 
    sourcePosition = Vector3d( min(maxBound.x, max(minBound.x, dx.x*(0.5 + floor(sourcePosition.x/dx.x)))),
                               min(maxBound.y, max(minBound.y, dx.y*(0.5 + floor(sourcePosition.y/dx.y)))),
                               min(maxBound.z, max(minBound.z, dx.z*(0.5 + floor(sourcePosition.z/dx.z)))));
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
        return -rho_0 *ballRadius *ballRadius / 4.0 / PI / r  * Sw * w * cos(w*t + phase);
    else 
        return numeric_limits<REAL>::quiet_NaN();
}

REAL boundaryEval( const Vector3d &x, const Vector3d &n, int obj_id, REAL t, int field_id, InterpolationFunction *interp )
{
    REAL bcResult = 0.0;

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


void writeData(const vector<REAL> & w, const REAL & timeStep, const int & timeStamp, const int & substep, const char * pattern, REAL endTime, int n);

void UnitTesting(); 

using namespace std;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main( int argc, char **argv )
{

    if (argc != 2){
        printf("**Usage: %s <solver_configuration_XML>\n", argv[0]);
        exit(1);
    }

    std::string fileName = std::string(argv[1]);

    Parser                  *parser = NULL;
    TriangleMesh<REAL>      *mesh   = NULL;
    ClosestPointField       *sdf    = NULL;
    BoundingBox              fieldBBox;
    REAL                     cellSize;
    int                      cellDivisions;


    REAL                     endTime = -1.0;

    Vector3Array             listeningPositions;
    Vector3d                 sourcePosition;
    BoundaryEvaluator        boundaryCondition;

    REAL                     timeStep;

    Parser::ImpulseResponseParms parms;


    /* Build parser from solver configuration. */
    parser = Parser::buildParser( fileName );
    if ( !parser ) throw std::runtime_error("**ERROR** Could not build parser from : "+fileName); 

    /* Build mesh. */
    mesh = parser->getMesh("impulse_response");
    if ( !mesh ) throw std::runtime_error("**ERROR** Could not build mesh from : " + fileName); 

    parms = parser->getImpulseResponseParms();

    endTime = parms._stopTime; 
    sourcePosition.x = parms._sourcePosition_x; 
    sourcePosition.y = parms._sourcePosition_y; 
    sourcePosition.z = parms._sourcePosition_z; 

    const char *pattern = parms._outputPattern.c_str();

    printf("Queried end time for the simulation = %lf\n", endTime);

    sdf = DistanceFieldBuilder::BuildSignedClosestPointField( parser->getMeshFileName().c_str(),
            parms._sdfResolution,
            parms._sdfFilePrefix.c_str() );

    // Scale this up to build a finite difference field
    fieldBBox  = BoundingBox(sdf->bmin(),sdf->bmax());
    fieldBBox *= parms._gridScale;

    cellDivisions = parms._gridResolution;

    cellSize = fieldBBox.minlength(); // use min length to check the CFL = c0*dt/dx
    cellSize /= (REAL)cellDivisions;

    timeStep = 1.0 / (REAL)( parms._timeStepFrequency);
    REAL CFL = parms._c * timeStep / cellSize; 

    cout << "fieldBBox -> [min, max] = " << fieldBBox.minBound() <<  ", " << fieldBBox.maxBound() << endl;
    cout << SDUMP(cellDivisions) << endl;
    cout << SDUMP( timeStep ) << endl;
    cout << SDUMP( cellSize ) << endl; 
    cout << SDUMP( CFL ) << endl;

    // for boundary interpolation : not used
    InterpolationFunction * interp = new InterpolationMitchellNetravali( 0.1 );
    boundaryCondition = boost::bind( boundaryEval, _1, _2, _3, _4, _5, interp );


    /* Callback function for logging pressure at each time step. */
    PML_WaveSolver::WriteCallbackIndividual dacallback = boost::bind(writeData, _1, _2, _3, parms._subSteps, pattern, endTime, listeningPositions.size());

    /// pass in listening position
    PML_WaveSolver        solver( timeStep, fieldBBox, cellSize,
                                  *mesh, *sdf,
                                  0.0, /* distance tolerance */
                                  true, /* use boundary */
                                  &listeningPositions,
                                  NULL, /* No output file */
                                  NULL, /* Write callback */
                                  NULL,// uncomment this if don't want to listen and comment the next one
                                  //&dacallback, /* Write callback */
                                  parms._subSteps,
                                  1, /* acceleration directions */
                                  endTime );

    AdjustSourcePosition(fieldBBox.minBound(), fieldBBox.maxBound(), cellDivisions, sourcePosition);
    ///// initialize system with initial condition /////
    //const REAL ic_stddev = 0.005; 
    //const Vector3d sourcePosition(sound_source); 
    //InitialConditionEvaluator initial = boost::bind(Gaussian_3D, _1, sourcePosition, ic_stddev);
    //InitialConditionEvaluator initial = boost::bind(PointGaussian_3D, _1, sourcePosition);
    //solver.initSystemNontrivial( 0.0, &initial ); 


    ///// initialize system with external source. ///// 
    cout << "Source position is set at " << sourcePosition << endl;
    solver.initSystem(0.0);
    const REAL widthSpace = cellSize; 
    //const REAL widthSpace = 0.2; 
    const REAL widthTime  = 1./10000/10.; 
    //const REAL widthTime  = timeStep*2.; 
    //const REAL offsetTime = 4.*widthTime; 
    const REAL offsetTime = 5.0*widthTime; // the center of the gaussian in time
    //const REAL offsetTime = 0.0; 
    const REAL normalizeConstant = 1.0 / pow(sqrt_2_pi*widthSpace,3); // for normalizing the gaussian 
    ExternalSourceEvaluator source = boost::bind(Gaussian_3D_erf_time, _1, _2, sourcePosition, widthSpace, widthTime, offsetTime, normalizeConstant); 

    const REAL PML_width=5.0; 
    const REAL PML_strength=1000000.0; 
    solver.SetExternalSource( &source ); 
    solver.setPMLBoundaryWidth( PML_width, PML_strength );
    //solver.setPMLBoundaryWidth( 10.0, 100000.0 );
    //solver.setPMLBoundaryWidth( 20.0, 100000.0 );


    // if harmonic sources 
    //const REAL w = 2.0*3.1415926*2000; // 500 Hz test
    //const REAL Sw = 1.0; // source strength
    //const REAL ballRadius = cellSize*3; // threshold for discretizing point source position
    //cout << SDUMP(ballRadius) << endl;
    //const REAL phase = 0;
    //HarmonicSourceEvaluator sourceFunction = boost::bind( Harmonic_Source, _1, sourcePosition, _2, w, Sw, ballRadius, parms._density, phase ); 
    //solver.setHarmonicSource( &sourceFunction );  // WILL NEGATE ALL INITIALIZATION! 
 

    /// write grid and pressure value
    Eigen::Vector3i minDrawBound(0,0,0); 
    Eigen::Vector3i maxDrawBound(solver.fieldDivisions()[0], solver.fieldDivisions()[1], solver.fieldDivisions()[2]); 
    printf( "domain size: (%u, %u, %u)\n", maxDrawBound[0], maxDrawBound[1], maxDrawBound[2] );

    const int Ndivision = maxDrawBound[0] * maxDrawBound[1] * maxDrawBound[2];
    Eigen::MatrixXd vertexPosition( Ndivision, 3 ); // vertex position
    std::shared_ptr<Eigen::MatrixXd> vertexPressure( new Eigen::MatrixXd(Ndivision, 1) ); // full pressure
    Eigen::MatrixXi vertexIndex( Ndivision, 3 ); // pressure vertex index

    
    // for vtk dump
    Eigen::Vector3d minBound = ToEigenVector3<double>(fieldBBox.minBound());
    Eigen::Vector3d maxBound = ToEigenVector3<double>(fieldBBox.maxBound());
    Eigen::Vector3i cellCount= ToEigenVector3<int>(solver.fieldDivisions()); 
    std::shared_ptr<Eigen::MatrixXd> cellCenteredPosition(new Eigen::MatrixXd()); 
    UniformGrid::GetAllCellCenterPosition(minBound, maxBound, cellCount,*cellCenteredPosition); 

    int count = 0;
    Vector3d vPosition;
    for ( int kk=0; kk<maxDrawBound[2]; kk++ )
        for ( int jj=0; jj<maxDrawBound[1]; jj++ ) 
            for ( int ii=0; ii<maxDrawBound[0]; ii++ ) 
            {
                Tuple3i  vIndex( ii, jj, kk );
                vPosition = solver.fieldPosition(  vIndex );

                vertexPosition( count, 0 ) = vPosition[0];
                vertexPosition( count, 1 ) = vPosition[1];
                vertexPosition( count, 2 ) = vPosition[2];

                vertexIndex( count, 0 ) = vIndex[0]; 
                vertexIndex( count, 1 ) = vIndex[1]; 
                vertexIndex( count, 2 ) = vIndex[2]; 


                VECTOR vPressure;
                solver.vertexPressure( vIndex, vPressure );

                (*vertexPressure)( count, 0 ) = vPressure[0];

                count ++; 

            }

    printf("writing solver settings\n"); 
    {

        char buffer[100]; 
        snprintf( buffer, 100, pattern, "solver_setting.txt"); 
        ofstream of(buffer); 

        of << std::setprecision(16) << std::fixed;
        of << "[sourceposition] = " << sourcePosition << endl;
        of << "[fieldBBox.minBound] = " << fieldBBox.minBound() << endl;
        of << "[fieldBBox.maxBound] = " << fieldBBox.maxBound() << endl;
        of << "[solver.fieldDivisions] = " << maxDrawBound.transpose() << endl;
        of << SDUMP( cellSize ) << endl;
        of << SDUMP( CFL ) << endl;
        of << SDUMP( timeStep ) << endl;
        of << SDUMP( endTime ) << endl;

        of << SDUMP( widthSpace ) << endl;
        of << SDUMP( widthTime ) << endl;
        of << SDUMP( offsetTime ) << endl;

        of << "precision : double (hard-coded)" << endl;
        of << "binary write method : IO writeMatrixXd method (hard-coded)" << endl;
        of.close();

    }


    printf("writing pressure vertex position\n");
    {

        char buffer[100]; 
        snprintf( buffer,100,pattern,"vertex_position.dat"); 
        IO::writeMatrixXd( vertexPosition, buffer, IO::BINARY );
    }

    printf("writing pressure vertex index\n");
    {

        char buffer[100]; 
        snprintf( buffer,100,pattern,"vertex_index.dat"); 
        IO::writeMatrixXi( vertexIndex, buffer, IO::BINARY );
    }

    printf("writing initial pressure\n");
    {
        char buf[100]; 
        char pressure_buf[100]; 
        snprintf( pressure_buf, 100, "pressure.dat" ); 
        snprintf( buf, 100, pattern, pressure_buf );
        IO::writeMatrixXd( *vertexPressure, buf, IO::BINARY );
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
                for ( int jj=0; jj<maxDrawBound[1]; jj++ ) 
                    for ( int ii=0; ii<maxDrawBound[0]; ii++ ) 
                    {
                        Tuple3i  vIndex( ii, jj, kk );

                        VECTOR vPressure;
                        solver.vertexPressure( vIndex, vPressure );

                        (*vertexPressure)( count, 0 ) = vPressure[0];

                        count ++; 

                    }


            printf( "writing pressure %u\n", nSteps/parms._subSteps );
            {
                char buf[100]; 
                char pressure_buf[100]; 
                snprintf( pressure_buf, 100, "pressure_%05u.dat", nSteps/parms._subSteps ); 
                snprintf( buf, 100, pattern, pressure_buf );

                IO::writeMatrixXd( *vertexPressure, buf, IO::BINARY );
                VTKConverter::VTKStructureGridWithScalarFromEigen(*cellCenteredPosition,*vertexPressure,std::string(buf)+".vtk","pressure",VTKConverter::BINARY,cellCount); 

                //UniformGrid::WriteVTKCellCenteredFromEigen(vertexPressure, minBound, maxBound, cellCount, std::string(buf), "test_data"); 
            }
        }


        nSteps ++; 


    }

    cout << "End of program. " << endl;



    return 0;

}

/* 
 * Write listened data. 
 */
void writeData(const vector<REAL> & w, const REAL & timeStep, const int & timeStamp, const int & substep, const char * pattern, REAL endTime, int n)
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



void UnitTesting(const char* filename)
{
    Parser  *parser = nullptr; 

    Parser::ImpulseResponseParms parms; 

    parser = Parser::buildParser( filename ); 

    parms = parser->getImpulseResponseParms(); 

    std::cout << SDUMP(parms._stopTime        ) << std::endl;
    std::cout << SDUMP(parms._outputPattern   ) << std::endl;
    std::cout << SDUMP(parms._listeningFile   ) << std::endl;
    std::cout << SDUMP(parms._sourcePosition_x) << std::endl;
    std::cout << SDUMP(parms._sourcePosition_y) << std::endl;
    std::cout << SDUMP(parms._sourcePosition_z) << std::endl;

}


