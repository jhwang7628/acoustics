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
#include "utils/Conversions.h"

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
    Reinitialize_PML_WaveSolver(useBoundary, 0.0);
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
    Reinitialize_PML_WaveSolver(useBoundary, 0.0);
}

PML_WaveSolver::PML_WaveSolver(PML_WaveSolver_Settings_Ptr settings, std::shared_ptr<FDTD_Objects> objects)
    : _waveSpeed(settings->soundSpeed), 
      _density(settings->airDensity), 
      _grid(BoundingBox(settings->cellSize,settings->cellDivisions),settings,objects),
      _cellSize(settings->cellSize), 
      _subSteps(settings->timeSavePerStep), 
      _endTime(settings->timeEnd), 
      _timeStep(settings->timeStepSize), 
      _timeIndex(0), 
      _N(1),
      _zSlice(-1),
      _listeningPositions(nullptr ), 
      _outputFile(nullptr), 
      _callback(nullptr),
      _useGhostCellBoundary(settings->useGhostCell),
      _objects(objects), 
      _waveSolverSettings(settings)
{
    Reinitialize_PML_WaveSolver(settings->useMesh, 0.0); 
}

void PML_WaveSolver::Reinitialize_PML_WaveSolver(const bool &useBoundary, const REAL &startTime)
{
    _currentTime = startTime;
    _timeIndex = 0;

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

    _pLaplacian.resizeAndWipe(_grid.numPressureCells(), _N); 
    _pCollocated[0].resizeAndWipe(_grid.numPressureCells(), _N); 
    _pCollocated[1].resizeAndWipe(_grid.numPressureCells(), _N); 
    _pCollocated[2].resizeAndWipe(_grid.numPressureCells(), _N); 
    _pCollocatedInd = 1; 

    _grid.initFieldRasterized( useBoundary );
    _grid.classifyCellsDynamic(_pFull, _p, _pGhostCellsFull, _pGhostCells, _v, useBoundary, true); 
    _grid.ResetCellHistory(true); // set all history to valid

    //if ( _listeningPositions )
    //{
    //    _waveOutput.resize( _listeningPositions->size() );
    //    for ( size_t i = 0; i < _waveOutput.size(); i++ )
    //        _waveOutput[ i ].resize( _N );
    //    cout << "Setting " << _waveOutput.size() << " listening positions" << endl;
    //}

    _listenerOutput.resizeAndWipe( _N );
    _sourceEvaluator = nullptr; 
}

int PML_WaveSolver::numVelocityCells(const int &dim) const
{ 
    if (dim==0) 
        return _grid.numVelocityCellsX(); 
    else if (dim==1) 
        return _grid.numVelocityCellsY(); 
    else if (dim==2) 
        return _grid.numVelocityCellsZ(); 
    else 
        throw std::runtime_error("**ERROR** dimension for velocity cells out of bounds"); 
} 

void PML_WaveSolver::SetGhostCellBoundary(const bool &isOn)
{ 
    _useGhostCellBoundary = isOn; 
    _grid.SetGhostCellBoundary(isOn); 
} 

// not used anymore, see Reinitialize_PML_WaveSolver
void PML_WaveSolver::initSystem( REAL startTime )
{
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

void PML_WaveSolver::FetchScalarData(const MATRIX &scalar, const ScalarField &field, const Vector3Array &listeningPoints, Eigen::MatrixXd &data) 
{
    const int N = listeningPoints.size(); 
    if (N==0) return; 

    _writeTimer.start(); 
    data.resize(N, 1); 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int pt_idx=0; pt_idx<N; ++pt_idx) 
    {
        IntArray neighbours; 
        T_MLS mls; 
        std::vector<MLSPoint, P_ALLOCATOR> points; 
        std::vector<MLSVal, V_ALLOCATOR> attributes; 
        const Vector3d &point = listeningPoints.at(pt_idx);
        const MLSPoint evalPt = Conversions::ToEigen(point); 
        field.enclosingNeighbours(point, neighbours); 
        for (IntArray::iterator it=neighbours.begin(); it!=neighbours.end();)
        {
            if (std::isnan(scalar(*it, 0)))
                it = neighbours.erase(it); 
            else
            {
                const MLSPoint pt = Conversions::ToEigen(field.cellPosition(*it)); 
                MLSVal val; 
                val << scalar(*it, 0);
                points.push_back(pt); 
                attributes.push_back(val); 
                ++it; 
            }
        }
        if (neighbours.size() < 4)
            throw std::runtime_error("**ERROR** Interpolation error: cannot construct interpolant");
        else
        {
            const MLSVal mlsVal = mls.lookup(evalPt, points, attributes); 
            data(pt_idx, 0) = mlsVal(0, 0);
        }
    }
    _writeTimer.pause(); 
}

void PML_WaveSolver::FetchPressureData(const Vector3Array &listeningPoints, Eigen::MatrixXd &data, const int dim)
{
#if 0
    const int N = listeningPoints.size(); 
    if (N==0) return; 

    _writeTimer.start(); 
    const ScalarField &field = _grid.pressureField(); 
    const int N_resultDimension = 1; 
    data.resize(N, N_resultDimension); 
    VECTOR output(N_resultDimension); 
    BoundingBox pressureBoundingBox = _grid.PressureBoundingBox(); 
    for (int ii=0; ii<N; ++ii) 
    {
        if (!pressureBoundingBox.isInside(listeningPoints.at(ii)))
        {
            std::cout << "**WARNING** Listening position " << listeningPoints.at(ii) << " is out of bounds. Skipping\n"; 
            continue; 
        }

        if (dim == -1) 
#ifdef USE_COLLOCATED
            field.interpolateVectorField(listeningPoints.at(ii), _pCollocated[_pCollocatedInd], output); 
#else
            field.interpolateVectorField(listeningPoints.at(ii), _pFull, output); 
#endif
        else if (dim == 0)
            field.interpolateVectorField(listeningPoints.at(ii), _p[0], output); 
        else if (dim == 1)
            field.interpolateVectorField(listeningPoints.at(ii), _p[1], output); 
        else if (dim == 2)
            field.interpolateVectorField(listeningPoints.at(ii), _p[2], output); 

        for (size_t jj=0; jj<N_resultDimension; ++jj) 
            if (std::isnan(output(jj)))
                throw std::runtime_error("**ERROR** interpolation outputs nan"); 
            else
                data(ii, jj) = output(jj);
    }
    _writeTimer.pause(); 
#else
    if (dim == -1) 
#ifdef USE_COLLOCATED
        FetchScalarData(_pCollocated[_pCollocatedInd], _grid.pressureField(), listeningPoints, data); 
#else
        FetchScalarData(_pFull, _grid.pressureField(), listeningPoints, data); 
#endif
    else if (dim == 0)
        FetchScalarData(_p[0], _grid.pressureField(), listeningPoints, data); 
    else if (dim == 1)
        FetchScalarData(_p[1], _grid.pressureField(), listeningPoints, data); 
    else if (dim == 2)
        FetchScalarData(_p[2], _grid.pressureField(), listeningPoints, data); 
#endif 
}

void PML_WaveSolver::FetchVelocityData(const Vector3Array &listeningPoints, const int &dimension, Eigen::MatrixXd &data)
{
    const int N = listeningPoints.size(); 
    if (N==0) return; 

    _writeTimer.start(); 
    const ScalarField &field = _grid.velocityField(dimension); 
    MATRIX &v = _v[dimension];
    const int N_resultDimension = 1; 
    data.resize(N, N_resultDimension); 
    VECTOR output(N_resultDimension); 
    BoundingBox velocityBoundingBox = _grid.VelocityBoundingBox(dimension); 
    for (int ii=0; ii<N; ++ii) 
    {
        if (!velocityBoundingBox.isInside(listeningPoints.at(ii)))
            continue; 

        field.interpolateVectorField(listeningPoints.at(ii), v, output); 
        for (size_t jj=0; jj<N_resultDimension; ++jj) 
            data(ii, jj) = output(jj);
    }
    _writeTimer.pause(); 
}

// get neareset neighbour cell type
void PML_WaveSolver::FetchPressureCellType(const Vector3Array &listeningPoints, Eigen::MatrixXd &data)
{
    const int N = listeningPoints.size(); 
    if (N==0) return; 

    const ScalarField &field = _grid.pressureField(); 
    const int N_resultDimension = 1; 
    data.resize(N, N_resultDimension); 

    BoundingBox pressureBoundingBox = _grid.PressureBoundingBox(); 
    auto & grid = GetGrid(); 
    for (int ii=0; ii<N; ++ii) 
    {
        if (!pressureBoundingBox.isInside(listeningPoints.at(ii)))
        {
            std::cout << "**WARNING** Listening position " << listeningPoints.at(ii) << " is out of bounds. Skipping\n"; 
            continue; 
        }

        // grab the neighbour indices
        IntArray neighbours; 
        field.enclosingNeighbours(listeningPoints.at(ii), neighbours); 

        // find cell types
        REAL distance = std::numeric_limits<REAL>::max(); 
        int index = -1; 
        for (size_t nei_idx=0; nei_idx<neighbours.size(); ++nei_idx)
        {
            const REAL currentDistance = (field.cellPosition(neighbours.at(nei_idx)) - listeningPoints.at(ii)).lengthSqr();
            if ( currentDistance < distance)
            {
                distance = currentDistance; 
                index = neighbours.at(nei_idx); 
            }
        }
        data(ii, 0) = grid.PressureCellType(index); 
    }
}

void PML_WaveSolver::FetchCell(const int &cellIndex, MAC_Grid::Cell &cell) const 
{
#ifdef USE_COLLOCATED 
    _grid.GetCell(cellIndex, _p, _pCollocated[_pCollocatedInd], _pGhostCellsFull, _v, cell); 
    cell.laplacian = _pLaplacian(cell.index, 0); 
#else
    _grid.GetCell(cellIndex, _p, _pFull, _pGhostCellsFull, _v, cell); 
#endif
}

void PML_WaveSolver::SampleAxisAlignedSlice(const int &dim, const REAL &offset, std::vector<MAC_Grid::Cell> &sampledCells) const
{
#ifdef USE_COLLOCATED 
    _grid.SampleAxisAlignedSlice(dim, offset, _p, _pCollocated[_pCollocatedInd], _pGhostCellsFull, _v, sampledCells); 
#else
    _grid.SampleAxisAlignedSlice(dim, offset, _p, _pFull, _pGhostCellsFull, _v, sampledCells); 
#endif
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

    if ( _endTime > 0.0 && _currentTime >= _endTime )
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
#ifdef USE_COLLOCATED
    stepCollocated(); 
#else
    stepLeapfrog();
#endif
    _timeIndex += 1;
    _stepTimer.pause();

#ifdef DEBUG_PRINT
    printf( "Average time step cost: %f sec\n", _stepTimer.getMsPerCycle() );
    printf( " - Algebra             took %f sec\n", _algebraTimer.getMsPerCycle() );
    printf( " - Memory operations   took %f sec\n", _memoryTimer.getMsPerCycle() );
    printf( " - Writing             took %f sec\n", _writeTimer.getMsPerCycle() );
    printf( " - Cell reclassify     took %f sec\n", _cellClassifyTimer.getMsPerCycle() );
    printf( " - Velocity update     took %f sec\n", _gradientTimer.getMsPerCycle() );
    printf( " - Pressure update     took %f sec\n", _divergenceTimer.getMsPerCycle() );

    if (_useGhostCellBoundary) 
        printf( " - Ghost-cell update   took %f sec\n",  _ghostCellTimer.getMsPerCycle() );
#endif

    if ( _endTime > 0.0 && _currentTime >= _endTime )
        return false;

    return true;
}

bool PML_WaveSolver::stepSystemHalf(const int &flag)
{
    if (flag == 0) // step velocity
    {
        _grid.classifyCellsDynamic_FAST(_pFull, _p, _pGhostCellsFull, _pGhostCells, _v, _waveSolverSettings->useMesh, false);
        _grid.InterpolateFreshPressureCell(_p[0], _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshPressureCell(_p[1], _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshPressureCell(_p[2], _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshPressureCell(_pFull, _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshVelocityCell(_v[0], 0, _timeStep, _currentTime);
        _grid.InterpolateFreshVelocityCell(_v[1], 1, _timeStep, _currentTime);
        _grid.InterpolateFreshVelocityCell(_v[2], 2, _timeStep, _currentTime);
        _grid.PML_velocityUpdate( _pFull, _pGhostCellsFull, _v[ 0 ], 0, _currentTime, _timeStep, _density );
        _grid.PML_velocityUpdate( _pFull, _pGhostCellsFull, _v[ 1 ], 1, _currentTime, _timeStep, _density );
        _grid.PML_velocityUpdate( _pFull, _pGhostCellsFull, _v[ 2 ], 2, _currentTime, _timeStep, _density );
    }
    else // step pressure
    {
        _grid.PML_pressureUpdate( _v[ 0 ], _p[ 0 ], _pFull, 0, _timeStep, _waveSpeed, _sourceEvaluator, _currentTime, _density );
        _grid.PML_pressureUpdate( _v[ 1 ], _p[ 1 ], _pFull, 1, _timeStep, _waveSpeed, _sourceEvaluator, _currentTime, _density );
        _grid.PML_pressureUpdate( _v[ 2 ], _p[ 2 ], _pFull, 2, _timeStep, _waveSpeed, _sourceEvaluator, _currentTime, _density );
        _grid.UpdatePMLPressure(_p, _pFull); 
        _currentTime += _timeStep;
        _timeIndex += 1;
        if ( _endTime > 0.0 && _currentTime >= _endTime )
            return false;
    }
    return true;
}

void PML_WaveSolver::vertexPressure( const Tuple3i &index, VECTOR &pressure ) const 
{
    if ( pressure.size() != _N )
        pressure.resizeAndWipe( _N );

#ifdef USE_COLLOCATED
    MATRIX::copy( pressure.data(), _pCollocated[_pCollocatedInd].data() + _grid.pressureFieldVertexIndex( index ) * _N, 1, _N ); 
#else
    MATRIX::copy( pressure.data(), _pFull.data() + _grid.pressureFieldVertexIndex( index ) * _N, 1, _N ); 
#endif
}

void PML_WaveSolver::vertexVelocity( const Tuple3i &index, const int &dim, VECTOR &velocity ) const 
{
    if ( velocity.size() != _N )
        velocity.resizeAndWipe( _N );

    MATRIX::copy( velocity.data(), _v[dim].data() + _grid.velocityFieldVertexIndex( index, dim ) * _N, 1, _N ); 
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
    //_grid.classifyCellsDynamic(_pFull, _p, _pGhostCellsFull, _pGhostCells, _v, _waveSolverSettings->useMesh, true);
    _grid.classifyCellsDynamic_FAST(_pFull, _p, _pGhostCellsFull, _pGhostCells, _v, _waveSolverSettings->useMesh, false);
    _cellClassifyTimer.pause(); 

    if (_useGhostCellBoundary)
    {
        // interpolate fresh cells 
        _freshCellTimer.start(); 
        _grid.InterpolateFreshPressureCell(_p[0], _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshPressureCell(_p[1], _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshPressureCell(_p[2], _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshPressureCell(_pFull, _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshVelocityCell(_v[0], 0, _timeStep, _currentTime);
        _grid.InterpolateFreshVelocityCell(_v[1], 1, _timeStep, _currentTime);
        _grid.InterpolateFreshVelocityCell(_v[2], 2, _timeStep, _currentTime);
        _freshCellTimer.pause(); 
        
        // update ghost cells 
        _ghostCellTimer.start(); 
        //_grid.PML_pressureUpdateGhostCells_Jacobi(_p[0], _pGhostCells[0], _timeStep, _waveSpeed, _currentTime, _density); 
        //_grid.PML_pressureUpdateGhostCells_Jacobi(_p[1], _pGhostCells[1], _timeStep, _waveSpeed, _currentTime, _density); 
        //_grid.PML_pressureUpdateGhostCells_Jacobi(_p[2], _pGhostCells[2], _timeStep, _waveSpeed, _currentTime, _density); 
        _grid.PML_pressureUpdateGhostCells_Jacobi(_pFull, _pGhostCellsFull, _timeStep, _waveSpeed, _currentTime, _density); 
        _ghostCellTimer.pause(); 
    }
    else 
    {
        // interpolate fresh cells 
        _freshCellTimer.start(); 
        _grid.InterpolateFreshPressureCell(_p[0], _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshPressureCell(_p[1], _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshPressureCell(_p[2], _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshPressureCell(_pFull, _timeStep, _currentTime, _density);  
        _grid.InterpolateFreshVelocityCell(_v[0], 0, _timeStep, _currentTime);
        _grid.InterpolateFreshVelocityCell(_v[1], 1, _timeStep, _currentTime);
        _grid.InterpolateFreshVelocityCell(_v[2], 2, _timeStep, _currentTime);
        _freshCellTimer.pause(); 
    }

    // Update velocity in each direction
    _gradientTimer.start();
    _grid.PML_velocityUpdate( _pFull, _pGhostCellsFull, _v[ 0 ], 0, _currentTime, _timeStep, _density );
    _grid.PML_velocityUpdate( _pFull, _pGhostCellsFull, _v[ 1 ], 1, _currentTime, _timeStep, _density );
    _grid.PML_velocityUpdate( _pFull, _pGhostCellsFull, _v[ 2 ], 2, _currentTime, _timeStep, _density );
    _gradientTimer.pause();

    // Use the new velocity to update pressure
    _divergenceTimer.start();
    //_grid.PML_pressureUpdateFull( _v, _pFull, _timeStep, _waveSpeed, _sourceEvaluator, _currentTime, _density );
    _grid.PML_pressureUpdate( _v[ 0 ], _p[ 0 ], _pFull, 0, _timeStep, _waveSpeed, _sourceEvaluator, _currentTime, _density );
    _grid.PML_pressureUpdate( _v[ 1 ], _p[ 1 ], _pFull, 1, _timeStep, _waveSpeed, _sourceEvaluator, _currentTime, _density );
    _grid.PML_pressureUpdate( _v[ 2 ], _p[ 2 ], _pFull, 2, _timeStep, _waveSpeed, _sourceEvaluator, _currentTime, _density );
    _divergenceTimer.pause();

    _algebraTimer.start();
    _grid.UpdatePMLPressure(_p, _pFull); 
    _algebraTimer.pause();

    _currentTime += _timeStep;
}

void PML_WaveSolver::stepCollocated()
{
    MATRIX &pLast = _pCollocated[(_pCollocatedInd+2)%3]; 
    MATRIX &pCurr = _pCollocated[ _pCollocatedInd     ]; 
    MATRIX &pNext = _pCollocated[(_pCollocatedInd+1)%3]; 

    // reclassify cells occupied by objects
    _cellClassifyTimer.start(); 
    _grid.classifyCellsDynamic_FAST(_pFull, _pCollocated, _pGhostCellsFull, _pGhostCells, _v, _waveSolverSettings->useMesh, false);
    _cellClassifyTimer.pause(); 
    _freshCellTimer.start(); 
    _grid.InterpolateFreshPressureCell(pLast, _timeStep, _currentTime, _density);  
    _grid.InterpolateFreshPressureCell(pCurr, _timeStep, _currentTime, _density);  
    _freshCellTimer.pause(); 
    _ghostCellTimer.start(); 
    //_grid.PML_pressureUpdateGhostCells_Jacobi(pCurr, _pGhostCellsFull, _timeStep, _waveSpeed, _currentTime, _density); 
    _grid.PML_pressureUpdateGhostCells(pCurr, _pGhostCellsFull, _timeStep, _waveSpeed, _currentTime, _density); 
    _ghostCellTimer.pause(); 

    // Use the new velocity to update pressure
    _divergenceTimer.start();
    _grid.PML_velocityUpdateCollocated(_currentTime, _p, pCurr, _v); 
    _grid.pressureFieldLaplacianGhostCell(pCurr, _pGhostCellsFull, _pLaplacian); 
    _grid.PML_pressureUpdateCollocated(_currentTime, _v, _p, pLast, pCurr, pNext, _pLaplacian); 
    _pCollocatedInd = (_pCollocatedInd + 1)%3; 
    _divergenceTimer.pause();

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
#ifdef USE_COLLOCATED
    _grid.PrintFieldExtremum(_pCollocated[_pCollocatedInd],"_pFull"); 
#else
    _grid.PrintFieldExtremum(_pFull,"_pFull"); 
#endif
    _grid.PrintFieldExtremum(_p[0],"_p[0]"); 
    _grid.PrintFieldExtremum(_p[1],"_p[1]"); 
    _grid.PrintFieldExtremum(_p[2],"_p[2]"); 
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
       << " start time      : " << solver._currentTime << "\n"
       << " stop time       : " << solver._endTime << "\n" 
       << " time step size  : " << solver._timeStep << "\n"
       << " CFL             : " << CFL << "\n"
       << " ghost cell boundary condition   : " << solver._useGhostCellBoundary << "\n"
       << "--------------------------------------------------------------------------------" 
       << std::flush; 
    return os; 
}
