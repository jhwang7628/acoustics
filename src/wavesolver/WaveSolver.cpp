//////////////////////////////////////////////////////////////////////
// WaveSolver_T.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "WaveSolver.h"

#include <utils/IO.h>
#include <utils/STLUtil.h>

//////////////////////////////////////////////////////////////////////
// Constructor
//
// Provide the size of the domain (bbox), finite difference division
// size, and a signed distance function for the interior boundary.
//////////////////////////////////////////////////////////////////////
WaveSolver::WaveSolver( REAL timeStep,
                        const BoundingBox &bbox, REAL cellSize,
                        const TriMesh &mesh,
                        const DistanceField &distanceField,
                        bool useLeapfrog,
                        REAL distanceTolerance,
                        bool useBoundary,
                        bool rasterize,
                        const Vector3Array *listeningPositions,
                        const char *outputFile,
                        WriteCallback *callback,
                        int subSteps, int N )
    : _timeStep( timeStep ),
      _timeIndex( 0 ),
      _laplacian( bbox, cellSize, mesh, distanceField, distanceTolerance, N ),
      _useLeapfrog( useLeapfrog ),
      _listeningPositions( listeningPositions ),
      _outputFile( outputFile ),
      _callback( callback ),
      _subSteps( subSteps ),
      _N( N ),
      _laplacianTimer( "Laplacian" ),
      _boundaryTimer( "Laplacian boundary" )
{
    if ( useLeapfrog )
    {
        _p.resizeAndWipe( _laplacian.numCells(), N );
        _v.resizeAndWipe( _laplacian.numCells(), N );
        _a.resizeAndWipe( _laplacian.numCells(), N );
    }
    else
    {
        // Zero everything out to start with (corresponds to zero initial
        // conditions)
        _p0.resizeAndWipe( _laplacian.numCells(), N );
        _p1.resizeAndWipe( _laplacian.numCells(), N );
        _p2.resizeAndWipe( _laplacian.numCells(), N );
    }

    _workspace.resizeAndWipe( _laplacian.numCells(), N );

    if ( rasterize )
    {
        _laplacian.initFieldRasterized( useBoundary );
    }
    else
    {
        _laplacian.initField( useBoundary );
    }

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
WaveSolver::WaveSolver( REAL timeStep,
        const BoundingBox &bbox, REAL cellSize,
        vector<const TriMesh *> &meshes,
        vector<const DistanceField *> &boundaryFields,
        bool useLeapfrog,
        REAL distanceTolerance,
        bool useBoundary,
        bool rasterize,
        const Vector3Array *listeningPositions,
        const char *outputFile,
        WriteCallback *callback,
        int subSteps, int N )
: _timeStep( timeStep ),
    _timeIndex( 0 ),
    _laplacian( bbox, cellSize, meshes, boundaryFields, distanceTolerance, N ),
    _useLeapfrog( useLeapfrog ),
    _listeningPositions( listeningPositions ),
    _outputFile( outputFile ),
    _callback( callback ),
    _subSteps( subSteps )
{
    if ( useLeapfrog )
    {
        _p.resizeAndWipe( _laplacian.numCells(), N );
        _v.resizeAndWipe( _laplacian.numCells(), N );
        _a.resizeAndWipe( _laplacian.numCells(), N );
    }
    else
    {
        // Zero everything out to start with (corresponds to zero initial
        // conditions)
        _p0.resizeAndWipe( _laplacian.numCells(), N );
        _p1.resizeAndWipe( _laplacian.numCells(), N );
        _p2.resizeAndWipe( _laplacian.numCells(), N );
    }

    _workspace.resizeAndWipe( _laplacian.numCells(), N );

    if ( rasterize )
    {
        _laplacian.initFieldRasterized( useBoundary );
    }
    else
    {
        _laplacian.initField( useBoundary );
    }

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
WaveSolver::~WaveSolver()
{
}

//////////////////////////////////////////////////////////////////////
// Time steps the system in the given interval
//////////////////////////////////////////////////////////////////////
void WaveSolver::solveSystem( REAL startTime, REAL endTime,
        const BoundaryEvaluator &bcEvaluator )
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveSolver::initSystem( REAL startTime,
        VECTOR *initialPressure,
        VECTOR *initialVelocity )
{
    _currentTime = startTime;
    _timeIndex = 0;

    if ( _useLeapfrog )
    {
        _p.clear();
        _v.clear();

        if ( initialPressure )
        {
            _p.copyInplace( *initialPressure );
        }

        if ( initialVelocity )
        {
            _v.copyInplace( *initialVelocity );
        }
    }
    else
    {
        _p0.clear();
        _p1.clear();
        _p2.clear();

        if ( initialPressure )
        {
            _p0.copyInplace( *initialPressure );
            _p1.copyInplace( _p0 );
            _p2.copyInplace( _p0 );
        }

        if ( initialVelocity )
        {
            _p1.copyInplace( *initialVelocity );
            _p1 *= _timeStep;
            _p1 += _p0;
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Takes a single time step
//////////////////////////////////////////////////////////////////////
bool WaveSolver::stepSystem( const BoundaryEvaluator &bcEvaluator )
{
    _stepTimer.start();
    if ( _useLeapfrog )
    {
        stepLeapfrog( bcEvaluator );
    }
    else
    {
        stepCenteredDifference( bcEvaluator );
    }

    _timeIndex += 1;
    _stepTimer.pause();

    printf( "Average time step cost: %f ms\r", _stepTimer.getMsPerCycle() );
#if 0
    printf( "Time step %d took %f ms\n", _timeIndex, _stepTimer.getTotalSecs() );
    printf( "Algebra took %f ms\n", _algebraTimer.getTotalSecs() );
    printf( "Memory operations took %f ms\n", _memoryTimer.getTotalSecs() );
    printf( "Writing took %f ms\n", _writeTimer.getTotalSecs() );
    printf( "Laplacian application took %f ms\n",
            _laplacianTimer.getTotalSecs() );
    printf( "Laplacian boundary application took %f ms\n",
            _boundaryTimer.getTotalSecs() );
#endif

    return true;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveSolver::vertexPressure( const Tuple3i &index, VECTOR &pressure )
{
    if ( pressure.size() != _N )
    {
        pressure.resizeAndWipe( _N );
    }

    if ( _useLeapfrog )
    {
        MATRIX::copy( pressure.data(),
                // Move to the desired row
                _p.data() + _laplacian.fieldVertexIndex( index ) * _N,
                1, _N );
    }
    else
    {
        MATRIX::copy( pressure.data(),
                // Move to the desired row
                _p2.data() + _laplacian.fieldVertexIndex( index ) * _N,
                1, _N );
    }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void WaveSolver::writeWaveOutput() const
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
void WaveSolver::stepLeapfrog( const BoundaryEvaluator &bcEvaluator )
{
    _algebraTimer.start();

    REAL                       dt2 = _timeStep * _timeStep;
    REAL                       c2 = WAVE_SPEED * WAVE_SPEED;
    REAL                       laplacianMultiplier = c2;

    _p.parallelAxpy( _timeStep, _v );
    _p.parallelAxpy( dt2 / 2.0, _a );

    _v.parallelAxpy( _timeStep / 2.0, _a );

    _currentTime += _timeStep;
    _algebraTimer.pause();

    // Get new acceleration at next time step
#if 0
    _memoryTimer.start();
    _a.clear();
    _memoryTimer.pause();
#endif

    _laplacianTimer.start();
    _laplacian.apply( _p, _a, laplacianMultiplier );
    _laplacianTimer.pause();

    _boundaryTimer.start();
    _laplacian.applyBoundary( _a, bcEvaluator, _currentTime, laplacianMultiplier );
    _boundaryTimer.pause();

    // Add to the velocity
    _algebraTimer.start();
    _v.parallelAxpy( _timeStep / 2.0, _a );
    _algebraTimer.pause();

    _writeTimer.start();
    if ( _listeningPositions && ( _timeIndex % _subSteps ) == 0 )
    {
        REAL                     listenerOutput;
        const ScalarField       &field = _laplacian.field();

        for ( int i = 0; i < _listeningPositions->size(); i++ )
        {
#if 0
            listenerOutput = field.interpolateScalarField(
                    _listeningPositions->at( i ), _p );
#endif
            field.interpolateVectorField( _listeningPositions->at( i ),
                    _p, _listenerOutput );

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
    }
    _writeTimer.pause();
}

//////////////////////////////////////////////////////////////////////
// Time stepping based on standard centered difference of 2nd
// time derivative
//////////////////////////////////////////////////////////////////////
void WaveSolver::stepCenteredDifference(
        const BoundaryEvaluator &bcEvaluator )
{
    // p2 = c^2 * dt^2 * L(p1) + ( 2 p1 - p0 )     NO! / dt^2
    REAL                       dt2 = _timeStep * _timeStep;
    REAL                       c2 = WAVE_SPEED * WAVE_SPEED;
    REAL                       laplacianMultiplier = c2 * dt2;

    _p2.clear();
    _p2.copyInplace( _p1 );
    _p2 *= 2.0;
    _p2 -= _p0;
    //_p2 /= dt2;

#if 0
    cout << SDUMP( _p2.absmax() ) << SDUMP( _p2.max_() );
    cout << SDUMP( _p2.min_() );
    cout << SDUMP( _currentTime ) << endl;
#endif

    _laplacian.apply( _p1, _p2, laplacianMultiplier );
    _laplacian.applyBoundary( _p2, bcEvaluator, _currentTime, laplacianMultiplier );

    _p0.copyInplace( _p1 );
    _p1.copyInplace( _p2 );

    _currentTime += _timeStep;
}
