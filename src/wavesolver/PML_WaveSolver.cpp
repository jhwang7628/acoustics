//////////////////////////////////////////////////////////////////////
// PML_WaveSolver_T.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include <omp.h>
#include "math.h"
#include <utils/IO.h>
#include <utils/STLUtil.h>

#include <wavesolver/FDTD_Objects.h> 
#include "PML_WaveSolver.h"
#include "utils/Evaluator.h"
#include "utils/STL_Wrapper.h"

PML_WaveSolver::PML_WaveSolver( REAL timeStep,
                                const BoundingBox &bbox,
                                REAL cellSize,
                                const TriMesh &mesh,
                                const DistanceField &distanceField,
                                REAL waveSpeed, 
                                REAL density,
                                REAL distanceTolerance,
                                bool useBoundary,
                                const Vector3Array *listeningPositions,
                                const char *outputFile,
                                WriteCallback *callback,
                                int subSteps, int N,
                                REAL endTime) 
    : _waveSpeed(waveSpeed),
      _density(density),
      _grid( bbox, cellSize, mesh, distanceField, distanceTolerance, N ),
      _cellSize(cellSize),
      _subSteps( subSteps ),
      _endTime( endTime ),
      _timeStep( timeStep ),
      _timeIndex( 0 ),
      _N( N ),
      _zSlice( -1 ),
      _listeningPositions( listeningPositions ),
      _outputFile( outputFile ),
      _callback( callback )
{ 
    Reinitialize_PML_WaveSolver(useBoundary);
}

PML_WaveSolver::PML_WaveSolver( REAL timeStep,
                                const BoundingBox &bbox, REAL cellSize,
                                vector<const TriMesh *> &meshes,
                                vector<const DistanceField *> &boundaryFields,
                                REAL waveSpeed, 
                                REAL density, 
                                REAL distanceTolerance,
                                bool useBoundary,
                                const Vector3Array *listeningPositions,
                                const char *outputFile,
                                WriteCallback *callback,
                                int subSteps, int N,
                                REAL endTime) 
    : _waveSpeed(waveSpeed),
      _density(density),
      _grid( bbox, cellSize, meshes, boundaryFields, distanceTolerance, N ),
      _cellSize(cellSize),
      _subSteps( subSteps ),
      _endTime( endTime ),
      _timeStep( timeStep ),
      _timeIndex( 0 ),
      _N( N ),
      _zSlice( -1 ),
      _listeningPositions( listeningPositions ),
      _outputFile( outputFile ),
      _callback( callback )
{ 
    Reinitialize_PML_WaveSolver(useBoundary);
}

PML_WaveSolver::PML_WaveSolver(const PML_WaveSolver_Settings &settings, std::shared_ptr<FDTD_Objects> objects)
    : _waveSpeed(settings.soundSpeed), 
      _density(settings.airDensity), 
      _grid(BoundingBox(settings.cellSize,settings.cellDivisions),settings,objects),
      _cellSize(settings.cellSize), 
      _subSteps(settings.timeSavePerStep), 
      _endTime(settings.timeEnd), 
      _timeStep(settings.timeStepSize), 
      _timeIndex(0), 
      _N(1),
      _zSlice(-1),
      _listeningPositions(nullptr), 
      _outputFile(nullptr), 
      _callback(nullptr),
      _cornellBoxBoundaryCondition(settings.cornellBoxBoundaryCondition), 
      _useGhostCellBoundary(settings.useGhostCell),
      _objects(objects)
{
    Reinitialize_PML_WaveSolver(settings.useMesh); 
}

void PML_WaveSolver::Reinitialize_PML_WaveSolver(const bool &useBoundary)
{
    _pFull.resizeAndWipe( _grid.numPressureCells(),  _N );
    _p[ 0 ].resizeAndWipe( _grid.numPressureCells(), _N );
    _p[ 1 ].resizeAndWipe( _grid.numPressureCells(), _N );
    _p[ 2 ].resizeAndWipe( _grid.numPressureCells(), _N );
    _v[ 0 ].resizeAndWipe( _grid.numVelocityCellsX(),_N );
    _v[ 1 ].resizeAndWipe( _grid.numVelocityCellsY(),_N );
    _v[ 2 ].resizeAndWipe( _grid.numVelocityCellsZ(),_N );

    //_pLastTimestep.resizeAndWipe( _grid.numPressureCells(), _N ); 
    //_pThisTimestep.resizeAndWipe( _grid.numPressureCells(), _N ); 
    //_vThisTimestep[ 0 ].resizeAndWipe( _grid.numVelocityCellsX(), _N );
    //_vThisTimestep[ 1 ].resizeAndWipe( _grid.numVelocityCellsY(), _N );
    //_vThisTimestep[ 2 ].resizeAndWipe( _grid.numVelocityCellsZ(), _N );

    _grid.initFieldRasterized( useBoundary );

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
    _sourceEvaluator = nullptr; 
}

void PML_WaveSolver::SetCornellBoxBoundaryCondition(const bool &isOn)
{ 
    _cornellBoxBoundaryCondition = isOn; 
    _grid.SetCornellBoxBoundaryCondition(isOn); 
} 

void PML_WaveSolver::SetGhostCellBoundary(const bool &isOn)
{ 
    _useGhostCellBoundary = isOn; 
    _grid.SetGhostCellBoundary(isOn); 
} 

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
}

void PML_WaveSolver::initSystemNontrivial( const REAL startTime, const InitialConditionEvaluator * ic_eval )
{

    if ( NULL == ic_eval ) 
    {
        cout << "**WARNING** Initial condition evaluator not set! " << endl; 
        return; 
    }

    Timer<false> initTimer; 
    initTimer.start();

    _currentTime = startTime;
    _timeIndex = 0;

    initSystem( startTime ); 

    const int N_pcell = _grid.numPressureCells(); 

    // If want to use smooth Gaussian initialization
    for (int ii=0; ii<N_pcell; ii++) 
    {
        const Vector3d fieldPosition = _grid.pressureFieldPosition( ii ); 
        _pFull(ii, 0) = ( *ic_eval )( fieldPosition ); 
    }

    initTimer.pause(); 
    printf("Initialize system with ICs takes %f s.\n", initTimer.elapsed()); 
}

// TODO this currently produce bad results, maybe need to smooth velocity field as well
bool PML_WaveSolver::stepSystemWithRestart(const int &N_restart)
{
    _stepTimer.start();
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

    stepLeapfrog();  // now _pFull stores i+1 step pressure

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
    //printf( "Time step %d took %f s\n", _timeIndex, omp_get_wtime()-start);
    //_stepTimer.reset();

    if ( _endTime > 0.0 && (REAL)_timeIndex * _timeStep >= _endTime )
    {
        return false;
    }

    return true;

}

bool PML_WaveSolver::stepSystem(const BoundaryEvaluator &bcEvaluator)
{
    return false; 
}

bool PML_WaveSolver::stepSystem()
{
    _stepTimer.start();
    stepLeapfrog();
    _timeIndex += 1;
    _stepTimer.pause();

#ifdef DEBUG
    printf( "Average time step cost: %f sec\n", _stepTimer.getMsPerCycle() );
    printf( " - Algebra             took %f sec\n", _algebraTimer.getMsPerCycle() );
    printf( " - Memory operations   took %f sec\n", _memoryTimer.getMsPerCycle() );
    printf( " - Writing             took %f sec\n", _writeTimer.getMsPerCycle() );
    printf( " - Cell reclassify     took %f sec\n", _cellClassifyTimer.getMsPerCycle() );
    printf( " - Velocity update     took %f sec\n", _gradientTimer.getMsPerCycle() );
    printf( " - Pressure update     took %f sec\n", _divergenceTimer.getMsPerCycle() );

    if (_useGhostCellBoundary) 
        printf( " - Ghost-cell update took %f sec\n",  _ghostCellTimer.getMsPerCycle() );
#endif

    if ( _endTime > 0.0 && (REAL)_timeIndex * _timeStep >= _endTime )
        return false;

    return true;
}

void PML_WaveSolver::vertexPressure( const Tuple3i &index, VECTOR &pressure ) const 
{
    if ( pressure.size() != _N )
        pressure.resizeAndWipe( _N );

    MATRIX::copy( pressure.data(), _pFull.data() + _grid.pressureFieldVertexIndex( index ) * _N, 1, _N ); 
}

void PML_WaveSolver::writeWaveOutput() const
{
    if ( !_callback )
    {
        cerr << "**WARNING** No write callback has been set" << endl;
        return;
    }

    ( *_callback )( _waveOutput );
}

void PML_WaveSolver::stepLeapfrog()
{

    // reclassify cells occupied by objects
    _cellClassifyTimer.start(); 
    // FIXME DEBUG
//#ifdef DEBUG
//    _grid.classifyCellsDynamicAABB(true, _pFull, true);
//#else 
//    _grid.classifyCellsDynamicAABB(true, _pFull, false);
//#endif 
    _grid.classifyCellsDynamic(_pFull, _p, _v, true, false); 
    //_grid.initFieldRasterized(true);
    //_grid.classifyCellsDynamicAABB(true, _pFull, false);
    _cellClassifyTimer.pause(); 

    // deal with the fresh cell problem
    if (_useGhostCellBoundary)
    {
        _grid.FreshCellInterpolate(_pFull, _currentTime, _density); 
    }
    else 
    {
         _grid.InterpolatePressureCellHistory(_p[0], _currentTime, _density);  
         _grid.InterpolatePressureCellHistory(_p[1], _currentTime, _density);  
         _grid.InterpolatePressureCellHistory(_p[2], _currentTime, _density);  
         _grid.InterpolatePressureCellHistory(_pFull, _currentTime, _density);  
         //_grid.InterpolateVelocityCellHistory(_v[0], 0, _currentTime);
         //_grid.InterpolateVelocityCellHistory(_v[1], 1, _currentTime);
         //_grid.InterpolateVelocityCellHistory(_v[2], 2, _currentTime);
    }


    // Update velocity in each direction
    _gradientTimer.start();
    _grid.PML_velocityUpdate( _pFull, _v[ 0 ], 0, _currentTime, _timeStep, _density );
    _grid.PML_velocityUpdate( _pFull, _v[ 1 ], 1, _currentTime, _timeStep, _density );
    _grid.PML_velocityUpdate( _pFull, _v[ 2 ], 2, _currentTime, _timeStep, _density );
    _gradientTimer.pause();

    // Use the new velocity to update pressure
    _divergenceTimer.start();
    //_grid.PML_pressureUpdateFull( _v, _pFull, _timeStep, _waveSpeed, _sourceEvaluator, _currentTime, _density );
    _grid.PML_pressureUpdate( _v[ 0 ], _p[ 0 ], 0, _timeStep, _waveSpeed, _sourceEvaluator, _currentTime, _density );
    _grid.PML_pressureUpdate( _v[ 1 ], _p[ 1 ], 1, _timeStep, _waveSpeed, _sourceEvaluator, _currentTime, _density );
    _grid.PML_pressureUpdate( _v[ 2 ], _p[ 2 ], 2, _timeStep, _waveSpeed, _sourceEvaluator, _currentTime, _density );
    _divergenceTimer.pause();

    _algebraTimer.start();
    _pFull.parallelCopyAdd( _p[ 0 ], _p[ 1 ], _p[ 2 ] );
    _algebraTimer.pause();

    if (_useGhostCellBoundary)
    {
        _ghostCellTimer.start(); 
        _grid.PML_pressureUpdateGhostCells_Jacobi(_pFull, _timeStep, _waveSpeed, _currentTime, _density); 
        _ghostCellTimer.pause(); 
    }

    //_writeTimer.start();
    //if ( _listeningPositions && ( _timeIndex % _subSteps ) == 0 )
    //{
    //    REAL                     listenerOutput;
    //    const ScalarField       &field = _grid.pressureField();
    //    for ( size_t i = 0; i < _listeningPositions->size(); i++ )
    //    {
    //        field.interpolateVectorField( _listeningPositions->at( i ), _pFull, _listenerOutput );
    //        for ( int field_id = 0; field_id < _N; field_id++ )
    //        {
    //            listenerOutput = _listenerOutput( field_id );
    //            _waveOutput[ i ][ field_id ].push_back( listenerOutput );
    //        }
    //        if ( _outputFile && _timeIndex % 64 == 0 )
    //        {
    //            char buf[ 1024 ];
    //            for ( int field_id = 0; field_id < _N; field_id++ )
    //            {
    //                sprintf( buf, "%s_field_%d_position_%03lu.vector", _outputFile, field_id, i );
    //                writeVector( buf, _waveOutput[ i ][ field_id ] );
    //            }
    //        }
    //    }
    //}
    //_writeTimer.pause();
    //
    //if ( _zSlice >= 0 )
    //{
    //    char buf[ 1024 ];
    //    sprintf( buf, "slicedata/slicedata_%04d.matrix", _timeIndex );
    //    _grid.sampleZSlice( _zSlice, _pFull, _sliceData );
    //    _sliceData.write( buf );
    //}

    _currentTime += _timeStep;
}


REAL PML_WaveSolver::GetMaxCFL()
{
    const int N_vcellx = _grid.numVelocityCellsX(); 
    const int N_vcelly = _grid.numVelocityCellsY(); 
    const int N_vcellz = _grid.numVelocityCellsZ(); 

    assert(N_vcellx==N_vcelly&&N_vcellx==N_vcellz);

    const int &N = N_vcellx; 

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

    const REAL CFL = vmax * _timeStep / _cellSize; 

    return CFL; 
}

void PML_WaveSolver::PrintAllFieldExtremum()
{
    _grid.PrintFieldExtremum(_pFull,"_pFull"); 
    _grid.PrintFieldExtremum(_v[0],"_v[0]"); 
    _grid.PrintFieldExtremum(_v[1],"_v[1]"); 
    _grid.PrintFieldExtremum(_v[2],"_v[2]"); 
}

std::ostream &operator <<(std::ostream &os, const PML_WaveSolver &solver)
{
    const REAL CFL = solver._waveSpeed*solver._timeStep/solver._cellSize; 
    os << "--------------------------------------------------------------------------------\n" 
       << "Class PML_WaveSolver\n" 
       << "--------------------------------------------------------------------------------\n"
       << " wave speed      : " << solver._waveSpeed << "\n"
       << " density         : " << solver._density << "\n"
       << " cell size       : " << solver._cellSize << "\n"
       << " save per steps  : " << solver._subSteps << "\n"
       << " start time(TEMP): " << 0.0 << "\n"
       << " stop time       : " << solver._endTime << "\n" 
       << " time step size  : " << solver._timeStep << "\n"
       << " CFL             : " << CFL << "\n"
       << " cornell box boundary condition  : " << solver._cornellBoxBoundaryCondition << "\n"
       << " ghost cell boundary condition   : " << solver._useGhostCellBoundary << "\n"
       << "--------------------------------------------------------------------------------" 
       << std::flush; 
    return os; 
}



