/* 
 * impulse-response solve using the finite-difference time-domain (FDTD) wave solver. 
 *
 * ref: 
 *  [1] Chadwick 2012, Precomputed Acceleration Noise for Improved Rigid-Body Sound
 *  [2] Liu 1997, The perfectly matched layer for acoustic waves in absorptive media
 *  
 */


#include <config.h>
#include <TYPES.h>

#include "utils/IO/IO.h"
#include <limits>

#include <distancefield/closestPointField.h>
#include <distancefield/FieldBuilder.h>

#include <geometry/TriangleMesh.hpp>

#include <linearalgebra/Vector3.hpp>

#include <parser/Parser.h>

#include <wavesolver/PML_WaveSolver.h>
#include <wavesolver/WaveSolver.h>
#include "wavesolver/ImpulseResponseTypesAndConstants.h" 

#include <grid/Grid.h>

#include <boost/bind.hpp>

#include <iostream>
#include <string>

#include <unistd.h> 

#include "vtkConverter/vtkConverter.h" 

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



void TestSignedDistanceField(ClosestPointField *sdf) 
{

    // sweep a z-plane
    const double z = 0.0; 
    const int N = 800;
    const Eigen::Vector3d cellCount(N,N,1); 
    Eigen::MatrixXd grid(N*N,3); 
    Eigen::MatrixXd normal(N*N,3); 
    Eigen::MatrixXd data(N*N,1); 

    Vector3d position; 
    Vector3d normalBuffer; 

    for (int jj=0; jj<N; jj++) 
        for (int ii=0; ii<N; ii++) 
        {

            const int index = jj*N+ii;
            grid(index,0) = -0.5 + (double)ii/(double)N; 
            grid(index,1) = -0.5 + (double)jj/(double)N; 
            grid(index,2) = z; 

            position.x = grid(index,0); 
            position.y = grid(index,1); 
            position.z = grid(index,2); 

            normalBuffer = sdf->gradient(position); 
            normalBuffer.normalize(); 

            normal(index,0) = normalBuffer.x; 
            normal(index,1) = normalBuffer.y; 
            normal(index,2) = normalBuffer.z; 

            //std::cout << grid.row(jj*N+ii).transpose() << std::endl;
            data(jj*N+ii,0) = sdf->distance(position); 
        }


    VTKConverter::VTKStructureGridWithScalarFromEigen(grid, data, "sdf.vtk", "sdf", VTKConverter::BINARY, cellCount);
    VTKConverter::VTKStructureGridWithScalarFromEigen(grid, normal.col(0), "normal_x.vtk", "normal_x", VTKConverter::BINARY, cellCount);
    VTKConverter::VTKStructureGridWithScalarFromEigen(grid, normal.col(1), "normal_y.vtk", "normal_y", VTKConverter::BINARY, cellCount);
    VTKConverter::VTKStructureGridWithScalarFromEigen(grid, normal.col(2), "normal_z.vtk", "normal_z", VTKConverter::BINARY, cellCount);


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

/* 
 * multiple sources 3D gaussian in space x error function in time (integral of shifted gaussian) 
 *
 * implemented using simple superposition principle
 *
 * */ 
REAL Gaussian_3D_erf_time_multiple( const Vector3d &evaluatePosition, const REAL &evaluateTime, const SourceVector &sources) 
{

    REAL returnValue=0; 
    for (size_t ii=0; ii<sources.size(); ii++)
    {
        if ((evaluatePosition - sources[ii].position).length() < sources[ii].widthSpace*3.0)
            returnValue += Gaussian_3D_erf_time(evaluatePosition, evaluateTime, sources[ii].position, sources[ii].widthSpace, sources[ii].widthTime, sources[ii].offsetTime, sources[ii].normalizeConstant); 
    }

    return returnValue; 
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
void AdjustSourcePosition( const Vector3d &minBound, const Vector3d &maxBound, const int &divisions, SourceVector &sources) 
{
    Vector3d dx = (maxBound-minBound)/(double)divisions; 
    for (size_t ii=0; ii<sources.size(); ii++) 
    {
        sources[ii].position.x = min(maxBound.x, max(minBound.x, dx.x*(0.5 + floor(sources[ii].position.x/dx.x)))); 
        sources[ii].position.y = min(maxBound.y, max(minBound.y, dx.y*(0.5 + floor(sources[ii].position.y/dx.y)))); 
        sources[ii].position.z = min(maxBound.z, max(minBound.z, dx.z*(0.5 + floor(sources[ii].position.z/dx.z)))); 
    } 
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

REAL boundaryEval_HarmonicPulsation(const Vector3d &x, const Vector3d &n, int obj_id, REAL t, int field_id, const REAL &w, const REAL &phase)
{
    //REAL bcResult = 0.0;

    //TRACE_ASSERT( obj_id == 0 );

    //return -1.0;
    return cos(w*t + phase); 
    //return cos(w*t); 
    //return 1.0;
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
    BoundingBox              fieldBBox;  // bounding box for distance field
    BoundingBox              solverBBox; // bounding box for wavesolver
    REAL                     cellSize;
    int                      cellDivisions;
    REAL                     endTime = -1.0;
    Vector3Array             listeningPositions;
    //Vector3d                 sourcePosition;
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
    //sourcePosition.x = parms._sourcePosition_x; 
    //sourcePosition.y = parms._sourcePosition_y; 
    //sourcePosition.z = parms._sourcePosition_z; 

    const char *pattern = parms._outputPattern.c_str();


    cellDivisions = parms._gridResolution;


    sdf = DistanceFieldBuilder::BuildSignedClosestPointField( parser->getMeshFileName().c_str(),
            parms._sdfResolution,
            parms._sdfFilePrefix.c_str()
            );


    //TestSignedDistanceField(sdf);


    fieldBBox = BoundingBox(sdf->bmin(), sdf->bmax()); 

    
    std::srand ( unsigned ( std::time(0) ) );


    if (parms._cellSize >= 1E-14)
    {

        cellSize = parms._cellSize; 
        const REAL fieldLength = cellSize*(REAL)cellDivisions; 
        const Vec3d fieldMin(-fieldLength/2.0, -fieldLength/2.0, -fieldLength/2.0); 
        const Vec3d fieldMax(+fieldLength/2.0, +fieldLength/2.0, +fieldLength/2.0); 
        solverBBox = BoundingBox(fieldMin, fieldMax); 
        //fieldBBox = BoundingBox(fieldMin, fieldMax); 

        //sdf = DistanceFieldBuilder::BuildSignedClosestPointField( parser->getMeshFileName().c_str(),
        //        parms._sdfResolution,
        //        parms._sdfFilePrefix.c_str(),
        //        fieldMin, 
        //        fieldMax
        //        );
    }
    else
    {
        //sdf = DistanceFieldBuilder::BuildSignedClosestPointField( parser->getMeshFileName().c_str(),
        //        parms._sdfResolution,
        //        parms._sdfFilePrefix.c_str() );

        //fieldBBox  = BoundingBox(sdf->bmin(),sdf->bmax());
        //fieldBBox *= parms._gridScale;

        // Scale this up to build a finite difference field
        solverBBox = BoundingBox(sdf->bmin(),sdf->bmax()); 
        solverBBox*= parms._gridScale; 

        cellSize = solverBBox.minlength(); // use min length to check the CFL = c0*dt/dx
        cellSize /= (REAL)cellDivisions;
    }


    timeStep = 1.0 / (REAL)( parms._timeStepFrequency);
    REAL CFL = parms._c * timeStep / cellSize; 

    // for boundary interpolation : not used
    InterpolationFunction * interp = new InterpolationMitchellNetravali( 0.1 );
    boundaryCondition = boost::bind( boundaryEval, _1, _2, _3, _4, _5, interp );
    //const REAL w = 2.0*M_PI*parms._f; // 1kHz source
    //const REAL phase = 0.0; // zero phase shift
    //COUT_SDUMP(parms._f); 
    //COUT_SDUMP(w); 
    //COUT_SDUMP(phase); 
    //boundaryCondition = boost::bind( boundaryEval_HarmonicPulsation, _1, _2, _3, _4, _5, w, phase ); 

    /// pass in listening position
    PML_WaveSolver        solver( timeStep, solverBBox, cellSize,
                                  *mesh, *sdf,
                                  parms._c,
                                  parms._density,
                                  0.0, /* distance tolerance */
                                  parms._useMesh, /* use boundary */
                                  &listeningPositions,
                                  NULL, /* No output file */
                                  NULL, /* Write callback */
                                  //&dacallback, /* Write callback */
                                  parms._subSteps,
                                  1, /* acceleration directions */
                                  endTime);

    solver.SetCornellBoxBoundaryCondition(parms._cornellBoxBoundaryCondition); 


    SourceVector &sources = parms._sources; 
    AdjustSourcePosition(solverBBox.minBound(), solverBBox.maxBound(), cellDivisions, sources);
    ExternalSourceEvaluator source = boost::bind(Gaussian_3D_erf_time_multiple, _1, _2, sources); 
    ///// initialize system with initial condition /////
    //const REAL ic_stddev = 0.005; 
    //const Vector3d sourcePosition(sound_source); 
    //InitialConditionEvaluator initial = boost::bind(Gaussian_3D, _1, sourcePosition, ic_stddev);
    //InitialConditionEvaluator initial = boost::bind(PointGaussian_3D, _1, sourcePosition);
    //solver.initSystemNontrivial( 0.0, &initial ); 


    ///// initialize system with external source. ///// 
    //cout << "Source position is set at " << sourcePosition << endl;
    solver.initSystem(0.0);
    //const REAL &widthTime = parms._sourceWidthTime; 
    ////const REAL widthTime  = 1./10000./32.; 
    ////const REAL widthTime  = timeStep; 
    ////const REAL offsetTime = 4.*widthTime; 
    //const REAL offsetTime = 1./20000.; 
    ////const REAL offsetTime = 16.0*widthTime; // the center of the gaussian in time
    ////const REAL offsetTime = 0.0; 
    //const REAL widthSpace = parms._c * widthTime; 
    ////const REAL widthSpace = 0.2; 
    //const REAL normalizeConstant = 1.0 / pow(sqrt_2_pi*widthSpace,3); // for normalizing the gaussian 
    //ExternalSourceEvaluator source = boost::bind(Gaussian_3D_erf_time, _1, _2, sourcePosition, widthSpace, widthTime, offsetTime, normalizeConstant); 

    const REAL PML_width=10.0; 
    const REAL PML_strength=100000.0; 
    solver.SetExternalSource( &source ); 
    solver.setPMLBoundaryWidth( PML_width, PML_strength );
    solver.SetGhostCellBoundary(parms._useGhostCellBoundary); 
    //solver.setPMLBoundaryWidth( 10.0, 100000.0 );
    //solver.setPMLBoundaryWidth( 20.0, 100000.0 );


    // if harmonic sources 
    //const REAL w = 2.0*3.1415926*2000; // 500 Hz test
    //const REAL Sw = 1.0; // source strength
    //const REAL ballRadius = cellSize*3; // threshold for discretizing point source position
    //cout << SDUMP(ballRadius) << endl;
    //const REAL phase = 0;
 

    /// write grid and pressure value
    Eigen::Vector3i maxDrawBound(solver.fieldDivisions()[0], solver.fieldDivisions()[1], solver.fieldDivisions()[2]); 
    printf( "domain size: (%u, %u, %u)\n", maxDrawBound[0], maxDrawBound[1], maxDrawBound[2] );

    const int Ndivision = maxDrawBound[0] * maxDrawBound[1] * maxDrawBound[2];
    std::shared_ptr<Eigen::MatrixXd> vertexPressure( new Eigen::MatrixXd(Ndivision, 1) ); // full pressure
    
    // for vtk dump
    Eigen::Vector3d minBound = ToEigenVector3<double>(solverBBox.minBound());
    Eigen::Vector3d maxBound = ToEigenVector3<double>(solverBBox.maxBound());
    Eigen::Vector3i cellCount= ToEigenVector3<int>(solver.fieldDivisions()); 
    std::shared_ptr<Eigen::MatrixXd> cellCenteredPosition(new Eigen::MatrixXd()); 
    UniformGrid::GetAllCellCenterPosition(minBound, maxBound, cellCount,*cellCenteredPosition); 

    printf("writing solver settings\n"); 
    {

        char buffer[500]; 
        snprintf( buffer, 500, pattern, "solver_setting.txt"); 
        ofstream of(buffer); 

        of << std::setprecision(16) << std::fixed;
        for (size_t ii=0; ii<sources.size(); ii++) 
        {
            of << "source " << ii << std::endl; 
            of << " position           = " << sources[ii].position          << std::endl; 
            of << " width_space        = " << sources[ii].widthSpace        << std::endl; 
            of << " width_time         = " << sources[ii].widthTime         << std::endl; 
            of << " offset_time        = " << sources[ii].offsetTime        << std::endl; 
            of << " normalize_constant = " << sources[ii].normalizeConstant << std::endl; 
        }
        //of << "[sourceposition] = " << sourcePosition << endl;
        of << "solverBBox.minBound = " << minBound << std::endl;
        of << "solverBBox.maxBound = " << maxBound << std::endl;
        of << "solver.fieldDivisions = " << cellCount.transpose() << std::endl;
        of << "cellSize = " << cellSize << std::endl; 
        of << "CFL = " << CFL << std::endl; 
        of << "timeStep = " << timeStep << std::endl; 
        of << "endTime = " << endTime << std::endl; 


        of << "precision : double (hard-coded)" << endl;
        of << "binary write method : IO writeMatrixXd method (hard-coded)" << endl;
        of.close();

    }

    Eigen::MatrixXd vertexPosition(Ndivision,3); 
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
                count ++;

            }

    printf("writing pressure vertex position\n");
    {
    
        char buffer[500];
        snprintf( buffer,500,pattern,"vertex_position.dat");
        IO::writeMatrixXd( vertexPosition, buffer, IO::BINARY );
    }

    /// time step system
    bool continueStepping = true; 
    int nSteps = 0; 
    //const int N_restart = 10; 
    while ( continueStepping )
    {
        continueStepping = solver.stepSystem( boundaryCondition ); 
        //continueStepping = solver.stepSystemWithRestart( boundaryCondition, N_restart ); 

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
                char buf[500]; 
                char pressure_buf[500]; 

                snprintf( pressure_buf, 500, "pressure_%05u.dat", nSteps/parms._subSteps ); 
                snprintf( buf, 500, pattern, pressure_buf );
                IO::writeMatrixXd( *vertexPressure, buf, IO::BINARY );

                //snprintf( pressure_buf, 100, "pressure_%05u.vtk", nSteps/parms._subSteps ); 
                //snprintf( buf, 100, pattern, pressure_buf ); 
                //VTKConverter::VTKStructureGridWithScalarFromEigen(*cellCenteredPosition,*vertexPressure,std::string(buf)+".vtk","pressure",VTKConverter::BINARY,cellCount); 
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

    for (size_t ii=0; ii<w.size(); ii++) fprintf( fp, "%.16f\n", w[ii] );

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
    //std::cout << SDUMP(parms._sourcePosition_x) << std::endl;
    //std::cout << SDUMP(parms._sourcePosition_y) << std::endl;
    //std::cout << SDUMP(parms._sourcePosition_z) << std::endl;

}


