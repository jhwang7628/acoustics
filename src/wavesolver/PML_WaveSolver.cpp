//////////////////////////////////////////////////////////////////////
// PML_WaveSolver_T.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "PML_WaveSolver.h"

#include <util/IO.h>
#include <util/STLUtil.h>

//////////////////////////////////////////////////////////////////////
// Constructor
//
// Provide the size of the domain (bbox), finite difference division
// size, and a signed distance function for the interior boundary.
//////////////////////////////////////////////////////////////////////
PML_WaveSolver::PML_WaveSolver( Real timeStep,
                                const BoundingBox &bbox,
                                Real cellSize,
                                const TriMesh &mesh,
                                const DistanceField &distanceField,
                                Real distanceTolerance,
                                bool useBoundary,
                                const Vector3Array *listeningPositions,
                                const char *outputFile,
                                WriteCallback *callback,
                                int subSteps, int N,
                                Real endTime )
  : _timeStep( timeStep ),
    _timeIndex( 0 ),
    _grid( bbox, cellSize, mesh, distanceField, distanceTolerance, N ),
    _listeningPositions( listeningPositions ),
    _outputFile( outputFile ),
    _callback( callback ),
    _subSteps( subSteps ),
    _N( N ),
    _endTime( endTime ),
    _zSlice( -1 )
{
  _pFull.resizeAndWipe( _grid.numPressureCells(), N );
  _p[ 0 ].resizeAndWipe( _grid.numPressureCells(), N );
  _p[ 1 ].resizeAndWipe( _grid.numPressureCells(), N );
  _p[ 2 ].resizeAndWipe( _grid.numPressureCells(), N );
  _v[ 0 ].resizeAndWipe( _grid.numVelocityCellsX(), N );
  _v[ 1 ].resizeAndWipe( _grid.numVelocityCellsY(), N );
  _v[ 2 ].resizeAndWipe( _grid.numVelocityCellsZ(), N );

  _grid.initFieldRasterized( useBoundary );

  cout << SDUMP( cellSize ) << endl;

  if ( _listeningPositions )
  {
    _waveOutput.resize( _listeningPositions->size() );

    for ( int i = 0; i < _waveOutput.size(); i++ )
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
PML_WaveSolver::PML_WaveSolver( Real timeStep,
                                const BoundingBox &bbox, Real cellSize,
                                vector<const TriMesh *> &meshes,
                                vector<const DistanceField *> &boundaryFields,
                                Real distanceTolerance,
                                bool useBoundary,
                                const Vector3Array *listeningPositions,
                                const char *outputFile,
                                WriteCallback *callback,
                                int subSteps, int N,
                                Real endTime )
  : _timeStep( timeStep ),
    _timeIndex( 0 ),
    _grid( bbox, cellSize, meshes, boundaryFields, distanceTolerance, N ),
    _listeningPositions( listeningPositions ),
    _outputFile( outputFile ),
    _callback( callback ),
    _subSteps( subSteps ),
    _endTime( endTime ),
    _zSlice( -1 )
{
  _pFull.resizeAndWipe( _grid.numPressureCells(), N );
  _p[ 0 ].resizeAndWipe( _grid.numPressureCells(), N );
  _p[ 1 ].resizeAndWipe( _grid.numPressureCells(), N );
  _p[ 2 ].resizeAndWipe( _grid.numPressureCells(), N );
  _v[ 0 ].resizeAndWipe( _grid.numVelocityCellsX(), N );
  _v[ 1 ].resizeAndWipe( _grid.numVelocityCellsY(), N );
  _v[ 2 ].resizeAndWipe( _grid.numVelocityCellsZ(), N );

  _grid.initFieldRasterized( useBoundary );

  cout << SDUMP( cellSize ) << endl;

  if ( _listeningPositions )
  {
    _waveOutput.resize( _listeningPositions->size() );

    for ( int i = 0; i < _waveOutput.size(); i++ )
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
void PML_WaveSolver::setPMLBoundaryWidth( Real width, Real strength )
{
  _grid.setPMLBoundaryWidth( width, strength );
}

//////////////////////////////////////////////////////////////////////
// Time steps the system in the given interval
//////////////////////////////////////////////////////////////////////
void PML_WaveSolver::solveSystem( Real startTime, Real endTime,
                                  const BoundaryEvaluator &bcEvaluator )
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PML_WaveSolver::initSystem( Real startTime )
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

//////////////////////////////////////////////////////////////////////
// Takes a single time step
//////////////////////////////////////////////////////////////////////
bool PML_WaveSolver::stepSystem( const BoundaryEvaluator &bcEvaluator )
{
  _stepTimer.tick();
  stepLeapfrog( bcEvaluator );

  _timeIndex += 1;
  _stepTimer.tock();

  //printf( "Average time step cost: %f ms\r", _stepTimer.getMsPerCycle() );
#if 0
  printf( "Time step %d took %f ms\n", _timeIndex, _stepTimer.getTotalSecs() );
  printf( "Algebra took %f ms\n", _algebraTimer.getTotalSecs() );
  printf( "Memory operations took %f ms\n", _memoryTimer.getTotalSecs() );
  printf( "Writing took %f ms\n", _writeTimer.getTotalSecs() );
  printf( "Gradient computation took %f ms\n", _gradientTimer.getTotalSecs() );
  printf( "Divergence computation took %f ms\n",
          _divergenceTimer.getTotalSecs() );
#endif

  if ( _endTime > 0.0 && (Real)_timeIndex * _timeStep >= _endTime )
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

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void PML_WaveSolver::stepLeapfrog( const BoundaryEvaluator &bcEvaluator )
{
  const MATRIX              *velocities[] = { &_v[ 0 ], &_v[ 1 ], &_v[ 2 ] };

#if 0
  printf( "p is %d x %d\n", _p.rows(), _p.cols() );
  printf( "v_x is %d x %d\n", _v[ 0 ].rows(), _v[ 0 ].cols() );
  printf( "v_y is %d x %d\n", _v[ 1 ].rows(), _v[ 1 ].cols() );
  printf( "v_z is %d x %d\n", _v[ 2 ].rows(), _v[ 2 ].cols() );
#endif

  // Update velocity in each direction
  _gradientTimer.tick();
  _grid.PML_velocityUpdate( _pFull, bcEvaluator, _v[ 0 ], 0,
                            _currentTime, _timeStep );
  _grid.PML_velocityUpdate( _pFull, bcEvaluator, _v[ 1 ], 1,
                            _currentTime, _timeStep );
  _grid.PML_velocityUpdate( _pFull, bcEvaluator, _v[ 2 ], 2,
                            _currentTime, _timeStep );
  _gradientTimer.tock();

#if 0
  printf( "Max v_x: %f\n", _v[ 0 ].frobeniusNorm() );
  printf( "Max v_y: %f\n", _v[ 1 ].frobeniusNorm() );
  printf( "Max v_z: %f\n", _v[ 2 ].frobeniusNorm() );
#endif

  // Use the new velocity to update pressure
  _divergenceTimer.tick();
  _grid.PML_pressureUpdate( _v[ 0 ], _p[ 0 ], 0, _timeStep, WAVE_SPEED );
  _grid.PML_pressureUpdate( _v[ 1 ], _p[ 1 ], 1, _timeStep, WAVE_SPEED );
  _grid.PML_pressureUpdate( _v[ 2 ], _p[ 2 ], 2, _timeStep, WAVE_SPEED );
  _divergenceTimer.tock();

  _algebraTimer.tick();
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
  _algebraTimer.tock();

  _currentTime += _timeStep;

#if 0
  printf( "Max pressure: %f\n", _pFull.frobeniusNorm() );
#endif

  _writeTimer.tick();
  if ( _listeningPositions && ( _timeIndex % _subSteps ) == 0 )
  {
    Real                     listenerOutput;
    const ScalarField       &field = _grid.pressureField();

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

          writeRealVector( buf, _waveOutput[ i ][ field_id ] );
        }
      }
    }
  }
  _writeTimer.tock();

  if ( _zSlice >= 0 )
  {
    char                     buf[ 1024 ];

    sprintf( buf, "slicedata/slicedata_%04d.matrix", _timeIndex );

    _grid.sampleZSlice( _zSlice, _pFull, _sliceData );

    _sliceData.write( buf );
  }
}

