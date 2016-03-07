//////////////////////////////////////////////////////////////////////
// PML_WaveSolver_T.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "PML_WaveSolver.h"
#include "utils/Evaluator.h"

#include "math.h"

#include <utils/IO.h>
#include <utils/STLUtil.h>


//////////////////////////////////////////////////////////////////////
// Constructor
//
// Provide the size of the domain (bbox), finite difference division
// size, and a signed distance function for the interior boundary.
//////////////////////////////////////////////////////////////////////
PML_WaveSolver::PML_WaveSolver( REAL timeStep,
                                const BoundingBox &bbox,
                                REAL cellSize,
                                const TriMesh &mesh,
                                const DistanceField &distanceField,
                                REAL distanceTolerance,
                                bool useBoundary,
                                const Vector3Array *listeningPositions,
                                const char *outputFile,
                                WriteCallback *callback,
                                WriteCallbackIndividual *callbacki,
                                int subSteps, int N,
                                REAL endTime )
: 
    _grid( bbox, cellSize, mesh, distanceField, distanceTolerance, N ),
    _timeStep( timeStep ),
    _timeIndex( 0 ),
    _subSteps( subSteps ),
    _endTime( endTime ),
    _N( N ),
    _cellSize(cellSize),
    _listeningPositions( listeningPositions ),
    _outputFile( outputFile ),
    _callback( callback ),
    _callbackInd( callbacki ),
    _zSlice( -1 )
{
    cout << "CPU PML wavesolver starting (config 1)... " << endl;
    _pFull.resizeAndWipe( _grid.numPressureCells(), N );
    _p[ 0 ].resizeAndWipe( _grid.numPressureCells(), N );
    _p[ 1 ].resizeAndWipe( _grid.numPressureCells(), N );
    _p[ 2 ].resizeAndWipe( _grid.numPressureCells(), N );
    _v[ 0 ].resizeAndWipe( _grid.numVelocityCellsX(), N );
    _v[ 1 ].resizeAndWipe( _grid.numVelocityCellsY(), N );
    _v[ 2 ].resizeAndWipe( _grid.numVelocityCellsZ(), N );

    _pLastTimestep.resizeAndWipe( _grid.numPressureCells(), N ); 
    _pThisTimestep.resizeAndWipe( _grid.numPressureCells(), N ); 
    _vThisTimestep[ 0 ].resizeAndWipe( _grid.numVelocityCellsX(), N );
    _vThisTimestep[ 1 ].resizeAndWipe( _grid.numVelocityCellsY(), N );
    _vThisTimestep[ 2 ].resizeAndWipe( _grid.numVelocityCellsZ(), N );

    _grid.initFieldRasterized( useBoundary );

    cout << SDUMP( cellSize ) << endl;

    if ( _listeningPositions )
    {
        _waveOutput.resize( _listeningPositions->size() );

        for ( size_t i = 0; i < _waveOutput.size(); i++ )
        {
            _waveOutput[ i ].resize( _N );
        }

        cout << "Setting " << _waveOutput.size() << " listening positions" << endl;
    }

    _listenerOutput.resizeAndWipe( _N );
}

PML_WaveSolver::PML_WaveSolver( REAL timeStep,
                                const BoundingBox &bbox,
                                REAL cellSize,
                                const TriMesh &mesh,
                                const DistanceField &distanceField,
                                REAL distanceTolerance,
                                bool useBoundary,
                                const Vector3Array *listeningPositions,
                                const char *outputFile,
                                WriteCallback *callback,
                                int subSteps, int N,
                                REAL endTime )
: 
    _grid( bbox, cellSize, mesh, distanceField, distanceTolerance, N ),
    _timeStep( timeStep ),
    _timeIndex( 0 ),
    _subSteps( subSteps ),
    _endTime( endTime ),
    _N( N ),
    _cellSize(cellSize),
    _listeningPositions( listeningPositions ),
    _outputFile( outputFile ),
    _callback( callback ),
    _zSlice( -1 )
{
    cout << "CPU PML wavesolver starting (config 1)... " << endl;
    _pFull.resizeAndWipe( _grid.numPressureCells(), N );
    _p[ 0 ].resizeAndWipe( _grid.numPressureCells(), N );
    _p[ 1 ].resizeAndWipe( _grid.numPressureCells(), N );
    _p[ 2 ].resizeAndWipe( _grid.numPressureCells(), N );
    _v[ 0 ].resizeAndWipe( _grid.numVelocityCellsX(), N );
    _v[ 1 ].resizeAndWipe( _grid.numVelocityCellsY(), N );
    _v[ 2 ].resizeAndWipe( _grid.numVelocityCellsZ(), N );

    _pLastTimestep.resizeAndWipe( _grid.numPressureCells(), N ); 
    _pThisTimestep.resizeAndWipe( _grid.numPressureCells(), N ); 
    _vThisTimestep[ 0 ].resizeAndWipe( _grid.numVelocityCellsX(), N );
    _vThisTimestep[ 1 ].resizeAndWipe( _grid.numVelocityCellsY(), N );
    _vThisTimestep[ 2 ].resizeAndWipe( _grid.numVelocityCellsZ(), N );

    _grid.initFieldRasterized( useBoundary );

    cout << SDUMP( cellSize ) << endl;

    if ( _listeningPositions )
    {
        _waveOutput.resize( _listeningPositions->size() );

        for ( size_t i = 0; i < _waveOutput.size(); i++ )
        {
            _waveOutput[ i ].resize( _N );
        }

        cout << "Setting " << _waveOutput.size() << " listening positions" << endl;
    }

    _listenerOutput.resizeAndWipe( _N );
}


//////////////////////////////////////////////////////////////////////
// Provide the size of the domain (bbox), finite difference division
// size, and a list of SDFs for the interior boundary
//////////////////////////////////////////////////////////////////////
PML_WaveSolver::PML_WaveSolver( REAL timeStep,
        const BoundingBox &bbox, REAL cellSize,
        vector<const TriMesh *> &meshes,
        vector<const DistanceField *> &boundaryFields,
        REAL distanceTolerance,
        bool useBoundary,
        const Vector3Array *listeningPositions,
        const char *outputFile,
        WriteCallback *callback,
        WriteCallbackIndividual *callbacki,
        int subSteps, int N,
        REAL endTime )
: 
    _grid( bbox, cellSize, meshes, boundaryFields, distanceTolerance, N ),
    _timeStep( timeStep ),
    _timeIndex( 0 ),
    _subSteps( subSteps ),
    _endTime( endTime ),
    _cellSize(cellSize),
    _listeningPositions( listeningPositions ),
    _outputFile( outputFile ),
    _callbackInd( callbacki ),
    _zSlice( -1 )
{
    cout << "CPU PML wavesolver starting (config 2)... " << endl;
    _pFull.resizeAndWipe( _grid.numPressureCells(), N );
    _p[ 0 ].resizeAndWipe( _grid.numPressureCells(), N );
    _p[ 1 ].resizeAndWipe( _grid.numPressureCells(), N );
    _p[ 2 ].resizeAndWipe( _grid.numPressureCells(), N );
    _v[ 0 ].resizeAndWipe( _grid.numVelocityCellsX(), N );
    _v[ 1 ].resizeAndWipe( _grid.numVelocityCellsY(), N );
    _v[ 2 ].resizeAndWipe( _grid.numVelocityCellsZ(), N );

    _pLastTimestep.resizeAndWipe( _grid.numPressureCells(), N ); 
    _pThisTimestep.resizeAndWipe( _grid.numPressureCells(), N ); 
    _vThisTimestep[ 0 ].resizeAndWipe( _grid.numVelocityCellsX(), N );
    _vThisTimestep[ 1 ].resizeAndWipe( _grid.numVelocityCellsY(), N );
    _vThisTimestep[ 2 ].resizeAndWipe( _grid.numVelocityCellsZ(), N );

    _grid.initFieldRasterized( useBoundary );

    cout << SDUMP( cellSize ) << endl;

    if ( _listeningPositions )
    {
        _waveOutput.resize( _listeningPositions->size() );

        for ( size_t i = 0; i < _waveOutput.size(); i++ )
        {
            _waveOutput[ i ].resize( _N );
        }

        cout << "Setting " << _waveOutput.size() << " listening positions" << endl;
    }

    _listenerOutput.resizeAndWipe( _N );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
PML_WaveSolver::~PML_WaveSolver()
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PML_WaveSolver::setPMLBoundaryWidth( REAL width, REAL strength )
{
    _grid.setPMLBoundaryWidth( width, strength );
}


//////////////////////////////////////////////////////////////////////
// set the harmonic source distribution for the wavesolver
//////////////////////////////////////////////////////////////////////
void PML_WaveSolver::setHarmonicSource( const HarmonicSourceEvaluator * hs_eval)
{

    // FIXME hack to initialize out the initial condition 
    initSystem( 0.0 );


    Timer<false> initHSTimer; 
    initHSTimer.start(); 


    if ( NULL == hs_eval ) 
    {
        cout << "** Warning ** harmonic source not set" << endl;
        return; 
    }

    const int N_pcell = _grid.numPressureCells(); 


    for (int ii=0; ii<N_pcell; ii++) 
    {
        const Vector3d fieldPosition = _grid.pressureFieldPosition( ii ); 

        const REAL testValue = (*hs_eval)(0.0, fieldPosition); 

        if ( ! std::isnan(testValue) )
        {
            hsIndex.push_back( ii );
        }

    }

    int N_hsIndex = hsIndex.size();

    if ( N_hsIndex < 1 ) 
    {
        cout << "** Warning ** no harmonic source has been set! Check the threshold " << endl; 
    }
    else 
    {
        printf("%u discrete harmonic sources have been identified.\n", N_hsIndex); 
    }

    initHSTimer.pause();
    printf("Initialize harmonic source distribution takes %f s.\n", initHSTimer.elapsed()); 

}


void PML_WaveSolver::SetWaveSolverPointDataPosition() 
{
    std::cout << "Set wave solver point data position." << std::endl;

    int NCell = _grid.numPressureCells(); 
    Eigen::MatrixXd pos; 
    pos.resize( NCell, 3 );
    pos.setZero();

    Eigen::MatrixXi cellIndex; 
    cellIndex.resize( NCell, 3 ); 
    cellIndex.setZero(); 

    ScalarField pressureField = _grid.pressureField(); 

    for (int ii=0; ii<NCell; ii++) 
    {
        Vector3d pressureFieldPosition = _grid.pressureFieldPosition( ii );

        pos( ii, 0 ) = pressureFieldPosition.x; 
        pos( ii, 1 ) = pressureFieldPosition.y; 
        pos( ii, 2 ) = pressureFieldPosition.z; 

        Tuple3i ind = pressureField.cellIndex( ii ); 

        cellIndex( ii, 0 ) = ind.x; 
        cellIndex( ii, 1 ) = ind.y; 
        cellIndex( ii, 2 ) = ind.z; 


    }
    Tuple3i div = pressureField.cellDivisions(); 
    Eigen::Vector3i cellDivisions; 
    cellDivisions << div.x, div.y, div.z; 

    //_rawData->SetPosition( pos, cellIndex, cellDivisions );

}


void PML_WaveSolver::AddCurrentPressureToWaveSolverPointData()
{
    std::cout << "write data to wave solver point data buffer" << std::endl;

    ScalarField pressureField = _grid.pressureField(); 

}


//////////////////////////////////////////////////////////////////////
// Time steps the system in the given interval
//////////////////////////////////////////////////////////////////////
void PML_WaveSolver::solveSystem( REAL startTime, REAL endTime,
                                  const BoundaryEvaluator &bcEvaluator )
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PML_WaveSolver::initSystem( REAL startTime )
{
    _currentTime = startTime;
    _timeIndex = 0;

    _pFull.clear();
    _p[ 0 ].clear();
    _p[ 1 ].clear();
    _p[ 2 ].clear();
    _v[ 0 ].clear();
    _v[ 1 ].clear();
    _v[ 2 ].clear();

    _sourceEvaluator = nullptr; 
}

/* 
 * Initialize the field data using non-trivial initial conditions. 
 */
void PML_WaveSolver::initSystemNontrivial( const REAL startTime, const InitialConditionEvaluator * ic_eval )
{

    if ( NULL == ic_eval ) 
    {
        cout << "** Warning ** Initial condition evaluator not set! " << endl; 
        return; 
    }

    Timer<false> initTimer; 
    initTimer.start();

    _currentTime = startTime;
    _timeIndex = 0;

    initSystem( startTime ); 

    const int N_pcell = _grid.numPressureCells(); 


    /// determine which is the closest vertex to put 1
    //int closestVertex = -1; 
    //REAL smallestDistance = 9999999.0; 
    //for (int ii=0; ii<N_pcell; ii++) 
    //{

    //    const Vector3d fieldPosition = _grid.pressureFieldPosition( ii ); 

    //    REAL distance = ( *ic_eval )( fieldPosition ); 

    //    if ( distance < smallestDistance ) 
    //    {
    //        closestVertex = ii; 
    //        smallestDistance = distance; 
    //    }

    //}
    //_pFull( closestVertex, 0 ) = 1.0;





    /// If want to use smooth Gaussian initialization
    for (int ii=0; ii<N_pcell; ii++) 
    {
        const Vector3d fieldPosition = _grid.pressureFieldPosition( ii ); 

        _pFull(ii, 0) = ( *ic_eval )( fieldPosition ); 

    }

    //SetWaveSolverPointDataPosition();

    initTimer.pause(); 
    printf("Initialize system with ICs takes %f s.\n", initTimer.elapsed()); 


}

//////////////////////////////////////////////////////////////////////
// Takes a single time step
//////////////////////////////////////////////////////////////////////
#include <omp.h>


bool PML_WaveSolver::stepSystemWithRestart(const BoundaryEvaluator &bcEvaluator, const int &N_restart)
{
    _stepTimer.start();
    double start = omp_get_wtime();
    bool restartStep = false; 

    if (_timeIndex % N_restart == 0 && _timeIndex != 0) 
        restartStep = true; 

    if (_timeIndex % N_restart == N_restart - 1)  // need to cache the result for the next step, which is a restart step
    {
        _pLastTimestep.parallelCopy(_pFull); 
    }
    else if (restartStep)  // cache the current step before stepping
    {
        std::cout << " restart step\n"; 
        _vThisTimestep[0].parallelCopy(_v[0]); 
        _vThisTimestep[1].parallelCopy(_v[1]); 
        _vThisTimestep[2].parallelCopy(_v[2]); 
        _pThisTimestep.parallelCopy(_pFull); 
    }

    stepLeapfrog(bcEvaluator);  // now _pFull stores i+1 step pressure

    if (restartStep) 
    {
        _grid.SmoothFieldInplace(_pLastTimestep, _pThisTimestep, _pFull, 0.25, 0.5, 0.25); 
        // now the pressure at step i is smoothed, recover v and continue stepping
        _v[0].parallelCopy(_vThisTimestep[0]); 
        _v[1].parallelCopy(_vThisTimestep[1]); 
        _v[2].parallelCopy(_vThisTimestep[2]); 
        //stepLeapfrog(bcEvaluator); 
    }

    _timeIndex += 1;
    _stepTimer.pause();

    //printf( "Average time step cost: %f ms\r", _stepTimer.getMsPerCycle() );
    printf( "Time step %d took %f s\n", _timeIndex, omp_get_wtime()-start);
    REAL t = (REAL)_timeIndex * _timeStep; 
    cout << SDUMP( t ) << endl; 
    _stepTimer.reset();

    if ( _endTime > 0.0 && (REAL)_timeIndex * _timeStep >= _endTime )
    {
        return false;
    }

    return true;

}

bool PML_WaveSolver::stepSystem( const BoundaryEvaluator &bcEvaluator )
{
    _stepTimer.start();
    double start = omp_get_wtime();
    stepLeapfrog( bcEvaluator );

    _timeIndex += 1;
    _stepTimer.pause();

    //printf( "Average time step cost: %f ms\r", _stepTimer.getMsPerCycle() );
    printf( "Time step %d took %f s\n", _timeIndex, omp_get_wtime()-start);
    REAL t = (REAL)_timeIndex * _timeStep; 
    cout << SDUMP( t ) << endl; 
    _stepTimer.reset();
#if 0
    printf( "Algebra took %f ms\n", _algebraTimer.getTotalSecs() );
    printf( "Memory operations took %f ms\n", _memoryTimer.getTotalSecs() );
    printf( "Writing took %f ms\n", _writeTimer.getTotalSecs() );
    printf( "Gradient computation took %f ms\n", _gradientTimer.getTotalSecs() );
    printf( "Divergence computation took %f ms\n",
            _divergenceTimer.getTotalSecs() );
#endif

    if ( _endTime > 0.0 && (REAL)_timeIndex * _timeStep >= _endTime )
    {
        return false;
    }

    return true;
}

// Overload the step system function to handle sources
bool PML_WaveSolver::stepSystem( const BoundaryEvaluator &bcEvaluator, const HarmonicSourceEvaluator *hsEval )
{
    _stepTimer.start();
    double start = omp_get_wtime();
    stepLeapfrog( bcEvaluator, hsEval );

    _timeIndex += 1;
    _stepTimer.pause();

    //printf( "Average time step cost: %f ms\r", _stepTimer.getMsPerCycle() );
    printf( "Time step %d took %f s with sources\n", _timeIndex, omp_get_wtime()-start);
    _stepTimer.reset();
#if 0
    printf( "Algebra took %f ms\n", _algebraTimer.getTotalSecs() );
    printf( "Memory operations took %f ms\n", _memoryTimer.getTotalSecs() );
    printf( "Writing took %f ms\n", _writeTimer.getTotalSecs() );
    printf( "Gradient computation took %f ms\n", _gradientTimer.getTotalSecs() );
    printf( "Divergence computation took %f ms\n",
            _divergenceTimer.getTotalSecs() );
#endif

    if ( _endTime > 0.0 && (REAL)_timeIndex * _timeStep >= _endTime )
    {
        return false;
    }

    return true;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PML_WaveSolver::vertexPressure( const Tuple3i &index,
        VECTOR &pressure )
{
    if ( pressure.size() != _N )
    {
        pressure.resizeAndWipe( _N );
    }

    MATRIX::copy( pressure.data(),
            // Move to the desired row
            _pFull.data() + _grid.pressureFieldVertexIndex( index ) * _N,
            1, _N );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PML_WaveSolver::writeWaveOutput() const
{
    if ( !_callback )
    {
        cerr << "** WARNING ** No write callback has been set" << endl;
        return;
    }

    ( *_callback )( _waveOutput );
}

REAL PML_WaveSolver::GetMaxCFL()
{

    const int N_vcellx = _grid.numVelocityCellsX(); 
    const int N_vcelly = _grid.numVelocityCellsY(); 
    const int N_vcellz = _grid.numVelocityCellsZ(); 

    if ( N_vcellx != N_vcelly || N_vcellx != N_vcellz ) 
    {
        cout << "** Warning ** number of velocity cells not equal" << endl; 
        return std::numeric_limits<REAL>::quiet_NaN(); 
    }

    const int N = N_vcellx; 

    REAL vmax = -1;  // should be positive

    for (int ii=0; ii<N; ii++) 
    {
        REAL vx, vy, vz, v; 

        vx = _v[0](ii, 0); 
        vy = _v[1](ii, 0); 
        vz = _v[2](ii, 0); 
        
        v = sqrt( vx*vx + vy*vy + vz*vz ); 

        vmax = max(v, vmax); 

    }

    //cout << "vmax = " << vmax << endl;

    REAL CFL = vmax * _timeStep / _cellSize; 
    //cout << SDUMP(CFL) << endl;


    return CFL; 
    

}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PML_WaveSolver::stepLeapfrog( const BoundaryEvaluator &bcEvaluator )
{
    //const MATRIX              *velocities[] = { &_v[ 0 ], &_v[ 1 ], &_v[ 2 ] };

#if 0
    printf( "p is %d x %d\n", _p.rows(), _p.cols() );
    printf( "v_x is %d x %d\n", _v[ 0 ].rows(), _v[ 0 ].cols() );
    printf( "v_y is %d x %d\n", _v[ 1 ].rows(), _v[ 1 ].cols() );
    printf( "v_z is %d x %d\n", _v[ 2 ].rows(), _v[ 2 ].cols() );
#endif

    // Update velocity in each direction
    double start = omp_get_wtime();
    _gradientTimer.start();
    _grid.PML_velocityUpdate( _pFull, bcEvaluator, _v[ 0 ], 0,
                              _currentTime, _timeStep );
    _grid.PML_velocityUpdate( _pFull, bcEvaluator, _v[ 1 ], 1,
                              _currentTime, _timeStep );
    _grid.PML_velocityUpdate( _pFull, bcEvaluator, _v[ 2 ], 2,
                              _currentTime, _timeStep );
    _gradientTimer.pause();
#if 0
    printf( "Max v_x: %f\n", _v[ 0 ].frobeniusNorm() );
    printf( "Max v_y: %f\n", _v[ 1 ].frobeniusNorm() );
    printf( "Max v_z: %f\n", _v[ 2 ].frobeniusNorm() );
#endif

    // Use the new velocity to update pressure
    start = omp_get_wtime();
    _divergenceTimer.start();
    _grid.PML_pressureUpdate( _v[ 0 ], _p[ 0 ], 0, _timeStep, WAVE_SPEED, _sourceEvaluator, _currentTime );
    _grid.PML_pressureUpdate( _v[ 1 ], _p[ 1 ], 1, _timeStep, WAVE_SPEED, _sourceEvaluator, _currentTime );
    _grid.PML_pressureUpdate( _v[ 2 ], _p[ 2 ], 2, _timeStep, WAVE_SPEED, _sourceEvaluator, _currentTime );
    _divergenceTimer.pause();
    printf( "pressure %d took %f s\n", _timeIndex, omp_get_wtime()-start);

    _algebraTimer.start();
#if 0
    _pFull.copyInplace( _p[ 0 ] );
    _pFull += _p[ 1 ];
    _pFull += _p[ 2 ];
#endif
#if 0
    _pFull.parallelCopy( _p[ 0 ] );
    _pFull.parallelAxpy( _p[ 1 ] );
    _pFull.parallelAxpy( _p[ 2 ] );
#endif
    _pFull.parallelCopyAdd( _p[ 0 ], _p[ 1 ], _p[ 2 ] );
    _algebraTimer.pause();

    _currentTime += _timeStep;


    //REAL MaxCFL = GetMaxCFL();
    //cout << SDUMP( MaxCFL ) << endl; 
    

#if 0
    printf( "Max pressure: %f\n", _pFull.frobeniusNorm() );
#endif

    _writeTimer.start();
    start = omp_get_wtime();
    if ( _listeningPositions && ( _timeIndex % _subSteps ) == 0 )
    {
        //REAL                     listenerOutput;
        const ScalarField       &field = _grid.pressureField();

        // MY IMPLEMENTATION //
        if ( _callbackInd )
        {
            vector<REAL> MyWaveOutput; 
            MyWaveOutput.resize( _listeningPositions->size() ); 
            for ( size_t ii=0; ii<_listeningPositions->size(); ii++) 
            {
                field.interpolateVectorField( _listeningPositions->at( ii ), 
                                              _pFull, _listenerOutput ); 
                  
                MyWaveOutput[ii] = _listenerOutput( 0 ); 
            }

            //REAL nowTime = (REAL)_timeIndex * _timeStep; 
            (*_callbackInd)( MyWaveOutput, _timeStep, _timeIndex );

        }
        // END MY IMPLEMENTATION //

        /* 
        for ( int i = 0; i < _listeningPositions->size(); i++ )
        {
            field.interpolateVectorField( _listeningPositions->at( i ),
                    _pFull, _listenerOutput );

            for ( int field_id = 0; field_id < _N; field_id++ )
            {
                listenerOutput = _listenerOutput( field_id );

                _waveOutput[ i ][ field_id ].push_back( listenerOutput );
            }

            if ( _outputFile && _timeIndex % 64 == 0 )
            {
                char                   buf[ 1024 ];

                for ( int field_id = 0; field_id < _N; field_id++ )
                {
                    sprintf( buf, "%s_field_%d_position_%03d.vector", _outputFile,
                            field_id, i );

                    writeVector( buf, _waveOutput[ i ][ field_id ] );
                }
            }
        }
        */
    }
    _writeTimer.pause();

    printf( "Writing %d took %f s\n", _timeIndex, omp_get_wtime()-start);
    if ( _zSlice >= 0 )
    {
        char                     buf[ 1024 ];

        sprintf( buf, "slicedata/slicedata_%04d.matrix", _timeIndex );

        _grid.sampleZSlice( _zSlice, _pFull, _sliceData );

        _sliceData.write( buf );
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PML_WaveSolver::stepLeapfrog( const BoundaryEvaluator &bcEvaluator, const HarmonicSourceEvaluator *hsEvaluator )
{

    Timer<false> evalSourceTimer; 
    evalSourceTimer.start(); 

    if ( NULL == hsEvaluator ) 
    {
        cout << "** Warning ** Harmonic source evaluator not set! " << endl; 
    }
    else 
    {

        for ( unsigned int ii=0; ii<hsIndex.size(); ii++ )
        {

            const Vector3d fieldPosition = _grid.pressureFieldPosition( hsIndex[ii] ); 
            const REAL pSource = (*hsEvaluator)( _currentTime, fieldPosition );
            _pFull( hsIndex[ii], 0 ) = pSource; 

            cout << "source " << ii << " : " << SDUMP( pSource ) << endl;

        }

    }

    evalSourceTimer.pause(); 
    printf("Evaluating harmonic source takes %f s.\n", evalSourceTimer.elapsed()); 

    stepLeapfrog( bcEvaluator ); 

}

