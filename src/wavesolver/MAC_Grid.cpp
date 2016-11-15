//////////////////////////////////////////////////////////////////////
// MAC_Grid.cpp: Implementation of the MAC_Grid class
//
//////////////////////////////////////////////////////////////////////

#include "MAC_Grid.h"

#include <utils/IO.h>
#include <utils/trace.h>
#include <utils/SimpleTimer.h>
#include <utils/STL_Wrapper.h> 
#include <utils/Conversions.h>

#include <math/LeastSquareSurface.h>
#include <math/MLSModeInterpolator.hpp>
#include <boost/timer/timer.hpp>

#include <distancefield/trilinearInterpolation.h> 
#include <linearalgebra/SparseLinearSystemSolver.h>
//#include <graph/UndirectedGraph.h> 

#include <unistd.h> 

#ifdef USE_OPENMP
#include <omp.h>
#endif

//##############################################################################
// Static variable initialize
//##############################################################################
int MAC_Grid::GhostCell::valuePointer = 0;

//##############################################################################
//##############################################################################
MAC_Grid::MAC_Grid( const BoundingBox &bbox, REAL cellSize,
                    const TriMesh &mesh,
                    const DistanceField &boundarySDF,
                    REAL distanceTolerance, int N ) 
    : _distanceTolerance( distanceTolerance ),
      _pressureField( bbox, cellSize ),
      _N( N ),
      _PML_absorptionWidth( 1.0 ),
      _PML_absorptionStrength( 0.0 )
{
    _boundaryFields.push_back( &boundarySDF );
    _boundaryMeshes.push_back( &mesh );

    Reinitialize_MAC_Grid(bbox, cellSize); 
}

MAC_Grid::MAC_Grid( const BoundingBox &bbox, REAL cellSize,
                    vector<const TriMesh *> &meshes,
                    vector<const DistanceField *> &boundaryFields,
                    REAL distanceTolerance, int N ) 
    : _boundaryFields( boundaryFields ),
      _boundaryMeshes( meshes ),
      _distanceTolerance(distanceTolerance),
      _pressureField( bbox, cellSize ),
      _N( N ),
      _PML_absorptionWidth( 1.0 ),
      _PML_absorptionStrength( 0.0 )
{
    Reinitialize_MAC_Grid(bbox, cellSize); 
}

MAC_Grid::MAC_Grid(const BoundingBox &bbox, PML_WaveSolver_Settings_Ptr settings, std::shared_ptr<FDTD_Objects> objects)
    : _pressureField(bbox,settings->cellSize), 
      _N(1), // no longer supports multiple fields, 
      _useGhostCellBoundary(settings->useGhostCell),
      _objects(objects), 
      _waveSolverSettings(settings)
{
    setPMLBoundaryWidth(settings->PML_width, settings->PML_strength);
    Reinitialize_MAC_Grid(bbox,settings->cellSize); 
}

void MAC_Grid::Reinitialize_MAC_Grid(const BoundingBox &bbox, const REAL &cellSize)
{
    Vector3d                   xMin, yMin, zMin;
    Tuple3i                    xDivs, yDivs, zDivs;

    for ( size_t mesh_idx = 0; mesh_idx < _boundaryMeshes.size(); mesh_idx++ )
    {
        printf( "Boundary mesh %lu has %lu vertices and %lu triangles\n",
                mesh_idx, _boundaryMeshes[ mesh_idx ]->vertices().size(),
                _boundaryMeshes[ mesh_idx ]->triangles().size() );
    }

    // Put velocity values on the exterior boundary
    xMin = bbox.minBound() - Vector3d( cellSize / 2.0, 0.0, 0.0 );
    yMin = bbox.minBound() - Vector3d( 0.0, cellSize / 2.0, 0.0 );
    zMin = bbox.minBound() - Vector3d( 0.0, 0.0, cellSize / 2.0 );

    xDivs = _pressureField.cellDivisions() + Tuple3i( 1, 0, 0 );
    yDivs = _pressureField.cellDivisions() + Tuple3i( 0, 1, 0 );
    zDivs = _pressureField.cellDivisions() + Tuple3i( 0, 0, 1 );

    _velocityField[ 0 ] = ScalarField( xMin, xDivs, cellSize );
    _velocityField[ 1 ] = ScalarField( yMin, yDivs, cellSize );
    _velocityField[ 2 ] = ScalarField( zMin, zDivs, cellSize );

    const int N_pressureCells = _pressureField.numCells(); 

    printf( "Initialized MAC_Grid with %d cells\n", N_pressureCells );

    // resize all necessary arrays defined everywhere
    _isBulkCell.resize(N_pressureCells, false);
    _pressureCellHasValidHistory.resize(N_pressureCells, true); 
    _containingObject.resize(N_pressureCells, -1);
    if (_useGhostCellBoundary)
        _isGhostCell.resize(N_pressureCells, false);
    for (int ii=0; ii<3; ++ii) 
    {
        const int N_velocityCells = _velocityField[ii].numCells(); 
        _isVelocityInterfacialCell[ii].resize(N_velocityCells, false); 
        _isVelocityBulkCell[ii].resize(N_velocityCells, false); 
        _velocityCellHasValidHistory[ii].resize(N_velocityCells, true); 
    }
    _cellSize = cellSize; 
}

MAC_Grid::~MAC_Grid()
{
}

void MAC_Grid::initFieldRasterized( bool useBoundary )
{
    cout << "MAC_Grid::initFieldRasterized: Classifying cells" << endl << endl;
    classifyCells(useBoundary); 
}

void MAC_Grid::setPMLBoundaryWidth( REAL width, REAL strength )
{
    _PML_absorptionWidth = width * _pressureField.cellSize();
    _PML_absorptionStrength = strength;
}

void MAC_Grid::pressureFieldLaplacian(const MATRIX &value, MATRIX &laplacian) const
{
    if (laplacian.rows()!=value.rows() || laplacian.cols()!=value.cols())
        laplacian.resizeAndWipe(value.rows(), value.cols()); 
    const auto &field = _pressureField; 
    const int N_cells = field.numCells(); 
    const REAL scale = 1.0/pow(_waveSolverSettings->cellSize,2); 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int cell_idx=0; cell_idx<N_cells; ++cell_idx)
    {
        if (!_isBulkCell.at(cell_idx))
            continue; 
        const Tuple3i cellIndices = field.cellIndex(cell_idx); 
        Tuple3i buf = cellIndices; 
        // skip if its boundary
        if (cellIndices[0]==0 || cellIndices[0]==_waveSolverSettings->cellDivisions-1 ||
            cellIndices[1]==0 || cellIndices[1]==_waveSolverSettings->cellDivisions-1 ||
            cellIndices[2]==0 || cellIndices[2]==_waveSolverSettings->cellDivisions-1)
            continue; 
        laplacian(cell_idx, 0) = -6.0*(value(cell_idx, 0)); // v_i
        for (int dim=0; dim<3; ++dim)
        {
            buf[dim] += 1; // v_i+1
            laplacian(cell_idx, 0) += value(field.cellIndex(buf), 0); 
            buf[dim] -= 2; // v_i-1
            laplacian(cell_idx, 0) += value(field.cellIndex(buf), 0); 
            buf[dim] += 1;
        }
        laplacian(cell_idx, 0) *= scale; 
    }
}

void MAC_Grid::pressureFieldLaplacianGhostCell(const MATRIX &value, const FloatArray &ghostCellValue, MATRIX &laplacian) const
{
    if (laplacian.rows()!=value.rows() || laplacian.cols()!=value.cols())
        laplacian.resizeAndWipe(value.rows(), value.cols()); 
    const auto &field = _pressureField; 
    const int N_cells = field.numCells(); 
    const REAL scale = 1.0/pow(_waveSolverSettings->cellSize,2); 
    laplacian.clear();

#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int cell_idx=0; cell_idx<N_cells; ++cell_idx)
    {
        if (!_isBulkCell.at(cell_idx))
            continue; 
        const Tuple3i cellIndices = field.cellIndex(cell_idx); 
        Tuple3i bufPos = cellIndices, bufNeg = cellIndices; 
        int buf_iPos = -1, buf_iNeg = -1; 
        // skip if its boundary
        if (cellIndices[0]==0 || cellIndices[0]==_waveSolverSettings->cellDivisions-1 ||
            cellIndices[1]==0 || cellIndices[1]==_waveSolverSettings->cellDivisions-1 ||
            cellIndices[2]==0 || cellIndices[2]==_waveSolverSettings->cellDivisions-1)
            continue; 
        for (int dim=0; dim<3; ++dim)
        {
            bufPos[dim] += 1; // v_i+1
            bufNeg[dim] -= 1; // v_i-1
            buf_iPos = field.cellIndex(bufPos); 
            buf_iNeg = field.cellIndex(bufNeg); 
#ifdef USE_FV
            const std::shared_ptr<GhostCell> &gc_pos = (_isGhostCell.at(buf_iPos) ? _ghostCellsCollection.at(buf_iPos) : nullptr); 
            const std::shared_ptr<GhostCell> &gc_neg = (_isGhostCell.at(buf_iNeg) ? _ghostCellsCollection.at(buf_iNeg) : nullptr); 
            REAL gc_pos_val=0.0, gc_neg_val=0.0; 
            if (gc_pos) gc_pos_val = (gc_pos->validSample ? gc_pos->values.at(gc_pos->valuePointer).at(dim*2  ) : value(cell_idx,0)); 
            if (gc_neg) gc_neg_val = (gc_neg->validSample ? gc_neg->values.at(gc_neg->valuePointer).at(dim*2+1) : value(cell_idx,0)); 
            if (_isGhostCell.at(buf_iPos) && _isGhostCell.at(buf_iNeg)) // both sides are ghost cell
                laplacian(cell_idx, 0) += (  gc_pos_val
                                           + gc_neg_val
                                           - 2.0*value(cell_idx, 0)) / (9.0/16.0); 
            else if (_isGhostCell.at(buf_iPos) && !_isGhostCell.at(buf_iNeg)) // only right side is ghost cell
                laplacian(cell_idx, 0) += (  gc_pos_val
                                           + 0.75*value(buf_iNeg, 0)
                                           - (7.0/4.0)*value(cell_idx, 0)) / (21.0/32.0); 
            else if (!_isGhostCell.at(buf_iPos) && _isGhostCell.at(buf_iNeg)) // only left side is ghost cell
                laplacian(cell_idx, 0) += (  0.75*value(buf_iPos, 0)
                                           + gc_neg_val
                                           - (7.0/4.0)*value(cell_idx, 0)) / (21.0/32.0); 
            else // both side is bulk
                laplacian(cell_idx, 0) += (  value(buf_iPos, 0)
                                           + value(buf_iNeg, 0)
                                           - 2.0*value(cell_idx, 0)); 
#else
            if (_isGhostCell.at(buf_iPos) && _isGhostCell.at(buf_iNeg)) // both sides are ghost cell
                laplacian(cell_idx, 0) += (  ghostCellValue.at(_ghostCellsChildren.at(_ghostCellsInverse.at(buf_iPos)).at(dim*2))
                                           + ghostCellValue.at(_ghostCellsChildren.at(_ghostCellsInverse.at(buf_iNeg)).at(dim*2+1))
                                           - 2.0*value(cell_idx, 0)) / (9.0/16.0); 
            else if (_isGhostCell.at(buf_iPos) && !_isGhostCell.at(buf_iNeg)) // only right side is ghost cell
                laplacian(cell_idx, 0) += (  ghostCellValue.at(_ghostCellsChildren.at(_ghostCellsInverse.at(buf_iPos)).at(dim*2))
                                           + 0.75*value(buf_iNeg, 0)
                                           - (7.0/4.0)*value(cell_idx, 0)) / (21.0/32.0); 
            else if (!_isGhostCell.at(buf_iPos) && _isGhostCell.at(buf_iNeg)) // only left side is ghost cell
                laplacian(cell_idx, 0) += (  0.75*value(buf_iPos, 0)
                                           + ghostCellValue.at(_ghostCellsChildren.at(_ghostCellsInverse.at(buf_iNeg)).at(dim*2+1))
                                           - (7.0/4.0)*value(cell_idx, 0)) / (21.0/32.0); 
            else // both side is bulk
                laplacian(cell_idx, 0) += (  value(buf_iPos, 0)
                                           + value(buf_iNeg, 0)
                                           - 2.0*value(cell_idx, 0)); 
#endif
            bufPos[dim] -= 1; 
            bufNeg[dim] += 1; 
        }
        laplacian(cell_idx, 0) *= scale; 
    }
}

void MAC_Grid::PML_velocityUpdate(const MATRIX &p, const FloatArray &pGC, MATRIX &v, int dimension, REAL t, REAL timeStep, REAL density )
{
    //const IntArray        &bulkCells                    = _velocityBulkCells[ dimension ];
    const IntArray        &interfacialCells             = _velocityInterfacialCells[ dimension ];
    const ScalarField     &field                        = _velocityField[ dimension ];
    const int              N_cells = field.numCells(); 
    const FloatArray      &interfaceBoundaryCoefficients= _interfacialBoundaryCoefficients[ dimension ];
    const REAL             minus_one_over_density       = -1.0 / density; 
    const REAL             one_over_three_quarters_cellsize      = 1.0 / (0.75 * _waveSolverSettings->cellSize); 
    const REAL             n_dt_over_dx_rho = -timeStep/(_pressureField.cellSize()*density);

    // Handle all bulk cells
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(static) default(shared)
    #endif
    for (int cell_idx = 0; cell_idx < N_cells; ++cell_idx)
    {
        if (!_isVelocityBulkCell[dimension].at(cell_idx))
            continue; 

        Tuple3i cell_indices = field.cellIndex( cell_idx );
        const Vector3d cell_position = field.cellPosition(cell_indices);
        const REAL      absorptionCoefficient   = PML_absorptionCoefficient(cell_position, _PML_absorptionWidth, dimension);
        const bool      inPML = (absorptionCoefficient > 1E-12) ? true : false; 

        // If we're on a boundary, then set the derivative to zero (since
        // we should have d_n p = 0 on the boundary
        if (cell_indices[dimension] == 0 || cell_indices[dimension] == field.cellDivisions()[dimension] - 1)
        {
            for ( int i = 0; i < _N; i++ )
                v( cell_idx, i ) = 0.0;
            continue;
        }

        if (inPML) 
        {
            const REAL updateCoefficient    = PML_velocityUpdateCoefficient(absorptionCoefficient, timeStep);
            const REAL gradientCoefficient  = PML_pressureGradientCoefficient(absorptionCoefficient, timeStep, _pressureField.cellSize(), density);
            for (int i=0; i<_N; ++i)
                v(cell_idx, i) *= updateCoefficient;

            // Sample the derivative of the pressure field in the necessary
            // direction
            int     neighbour_idx;
            neighbour_idx = _pressureField.cellIndex(cell_indices);

            for (int i=0; i<_N; ++i)
                v(cell_idx, i) += gradientCoefficient*p(neighbour_idx, i);

            cell_indices[ dimension ] -= 1;
            neighbour_idx = _pressureField.cellIndex(cell_indices);

            for (int i=0; i<_N; ++i)
                v(cell_idx, i) -= gradientCoefficient*p(neighbour_idx, i);
        }
        else // can collapse a lot of the terms
        {
            // updateCoefficient   = 1
            // gradientCoefficient = -dt/(dx*rho)
            int neighbour_idx;
            neighbour_idx = _pressureField.cellIndex(cell_indices);

            for (int i=0; i<_N; ++i)
                v(cell_idx, i) += n_dt_over_dx_rho*p(neighbour_idx, i);

            cell_indices[dimension] -= 1;
            neighbour_idx = _pressureField.cellIndex(cell_indices);

            for (int i=0; i<_N; ++i)
                v(cell_idx, i) -= n_dt_over_dx_rho*p(neighbour_idx, i);
        }
    }

    // Handle interfacial cells
    if (_useGhostCellBoundary)
    {
        //PML_velocityUpdateAux(p, v, dimension, t, timeStep, density, interfacialCells); 
        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(static) default(shared)
        #endif
        for ( size_t interfacial_cell_idx = 0; interfacial_cell_idx < interfacialCells.size(); interfacial_cell_idx++ )
        {
            const int cell_idx = interfacialCells[ interfacial_cell_idx ];
            const REAL boundaryDirection = _interfacialBoundaryDirections[dimension].at(interfacial_cell_idx); 

            // fetch regular pressure cell and subdivided ghost cell
            const Tuple3i velocityCellIndices = field.cellIndex(cell_idx); 
            Tuple3i leftPressureCellIndices = velocityCellIndices; 
            leftPressureCellIndices[dimension] -= 1;
            Tuple3i rightPressureCellIndices = velocityCellIndices; 

            // one index points to regular pressure cell, the other to
            // subdivided ghost cell
            int h_over_2_cell_index = 0, h_over_4_cell_index = 0; 
            int gcParent; 
            if (boundaryDirection > 0.5)  // it can take only 1 or -1 
            {
                h_over_2_cell_index = _pressureField.cellIndex(rightPressureCellIndices); 
                gcParent = _pressureField.cellIndex(leftPressureCellIndices); 
                // get child index in ghost cell pressure array. its position
                // in the tree depends on the boundaryDirection
                h_over_4_cell_index = _ghostCellsChildren.at(_ghostCellsInverse[gcParent]).at(dimension*2 + 1);
            }
            else if (boundaryDirection < -0.5)
            {
                h_over_2_cell_index = _pressureField.cellIndex(leftPressureCellIndices); 
                gcParent = _pressureField.cellIndex(rightPressureCellIndices); 
                // get child index in ghost cell pressure array. its position
                // in the tree depends on the boundaryDirection
                h_over_4_cell_index = _ghostCellsChildren.at(_ghostCellsInverse[gcParent]).at(dimension*2); 
            }
            else
            {
                throw std::runtime_error("**ERROR** boundary direction undefined."); 
            }

            if (h_over_4_cell_index == -1) 
                throw std::runtime_error("**ERROR** queried child node does not exist."); 

            // finite-difference estimate of acceleration using pressure at the
            // two locations.
            v(cell_idx, 0) += timeStep * minus_one_over_density * (p(h_over_2_cell_index, 0)-pGC.at(h_over_4_cell_index)) * one_over_three_quarters_cellsize * boundaryDirection; 
        }
    }
    else
    {
        //Handle boundary cells.
        //
        //Note: we assume here that all boundary cells have zero absorption
        //coefficient
        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(static) default(shared)
        #endif
        for ( size_t interfacial_cell_idx = 0; interfacial_cell_idx < interfacialCells.size(); interfacial_cell_idx++ )
        {
            int                  cell_idx = interfacialCells[ interfacial_cell_idx ];
            int                  objectID;
            REAL                 bcEval;
            REAL                 coefficient;

            Vector3d             normal( 0.0, 0.0, 0.0 ); // normal of the interfacial cell boundary
            Vector3d             x = field.cellPosition( cell_idx ); // position of the interfacial cell

            normal[dimension] = _interfacialBoundaryDirections[dimension].at(interfacial_cell_idx); 
            coefficient = interfaceBoundaryCoefficients[ interfacial_cell_idx ];
            objectID = _interfacialBoundaryIDs[dimension].at(interfacial_cell_idx);

            for ( int i = 0; i < _N; i++ )
            {
                bcEval = _objects->GetPtr(objectID)->EvaluateBoundaryAcceleration(x, normal, t);
                //bcEval = _objects->EvaluateNearestVibrationalSources(x,normal,t); 
                //bcEval = bc( x, normal, objectID, t, i );
                bcEval *= coefficient;
                v( cell_idx, i ) += timeStep * bcEval;
            }

        }
    }
}

void MAC_Grid::PML_velocityUpdateCollocated(const REAL &simulationTime, const MATRIX (&pDirectional)[3], const MATRIX &pFull, MATRIX (&v)[3])
{
    const int N_cells = _pmlVelocityCells.size(); 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int ii=0; ii<N_cells; ++ii)
    {
        const auto &cell = _pmlVelocityCells.at(ii);
        const int &cell_idx = cell.index; 
        // update velocity
        v[cell.dimension](cell_idx, 0) *= cell.updateCoefficient; 
        v[cell.dimension](cell_idx, 0) += cell.gradientCoefficient * pFull(cell.neighbour_p_right, 0); 
        v[cell.dimension](cell_idx, 0) -= cell.gradientCoefficient * pFull(cell.neighbour_p_left , 0); 
    }
}

void MAC_Grid::PML_pressureUpdateCollocated(const REAL &simulationTime, const MATRIX (&v)[3], MATRIX (&pDirectional)[3], MATRIX &pLast, MATRIX &pCurr, MATRIX &pNext, MATRIX &laplacian)
{
    const REAL &timeStep = _waveSolverSettings->timeStepSize; 
    const REAL one_over_dx = 1.0/_waveSolverSettings->cellSize; 
    const REAL c2_k2 = pow(_waveSolverSettings->soundSpeed*timeStep, 2); 
    const int N_cells = _pressureField.numCells(); 
    const int N_pmlCells = _pmlPressureCells.size(); 
    const bool evaluateExternalSource = _objects->HasExternalPressureSources();
      
    // first update PML region
    for (int ii=0; ii<N_pmlCells; ++ii)
    {
        const auto &cell = _pmlPressureCells.at(ii);
        const int &cell_idx = cell.index; 
        for (int dim=0; dim<3; ++dim)
        {
            pDirectional[dim](cell_idx, 0) *= cell.updateCoefficient[dim]; 
            pDirectional[dim](cell_idx, 0) += cell.divergenceCoefficient[dim] * v[dim](cell.neighbour_v_right[dim], 0) * one_over_dx;
            pDirectional[dim](cell_idx, 0) -= cell.divergenceCoefficient[dim] * v[dim](cell.neighbour_v_left[dim] , 0) * one_over_dx;
        }
        pNext(cell_idx, 0) = pDirectional[0](cell_idx, 0) + pDirectional[1](cell_idx, 0) + pDirectional[2](cell_idx, 0); 
    }

    // update laplacian
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(static) default(shared)
    #endif
    for (int cell_idx = 0; cell_idx < N_cells; ++cell_idx)
    {
        if (!_isBulkCell.at(cell_idx) || _isPMLCell.at(cell_idx))
            continue; 

        const Vector3d cell_position = _pressureField.cellPosition( cell_idx );
        pNext(cell_idx, 0) = pCurr(cell_idx, 0) * 2.0 - pLast(cell_idx, 0) + laplacian(cell_idx, 0) * c2_k2; 
        // evaluate external sources only happens not in PML
        // Liu Eq (16) f6x term
        if (evaluateExternalSource)
            pNext(cell_idx, 0) += _objects->EvaluatePressureSources(cell_position, cell_position, simulationTime+0.5*timeStep)*timeStep;
    }

}

void MAC_Grid::PML_pressureUpdate( const MATRIX &v, MATRIX &pDirectional, MATRIX &pFull, int dimension, REAL timeStep, REAL c, const ExternalSourceEvaluator *sourceEvaluator, const REAL simulationTime, REAL density )
{
    const int N_cells = numPressureCells(); 
    //const size_t bulkCellSize = _bulkCells.size();
    //const bool evaluateExternalSource = (sourceEvaluator != nullptr);
    const bool evaluateExternalSource = _objects->HasExternalPressureSources();
    const REAL n_rho_c_square_dt = -density*c*c*timeStep;
    const REAL v_cellSize_inv = 1./_velocityField[dimension].cellSize();

    // We only have to worry about bulk cells here
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(static) default(shared)
    #endif
    for (int cell_idx = 0; cell_idx < N_cells; ++cell_idx)
    {
        if (!_isBulkCell.at(cell_idx))
            continue; 

        const Vector3d  cell_position           = _pressureField.cellPosition( cell_idx );
        const REAL      absorptionCoefficient   = PML_absorptionCoefficient( cell_position, _PML_absorptionWidth, dimension );
        const bool      inPML = (absorptionCoefficient > 1E-12) ? true : false; 

        if (inPML) 
        {
            const REAL      directionalCoefficient  = PML_directionalCoefficient( absorptionCoefficient, timeStep );
            const REAL      updateCoefficient       = PML_pressureUpdateCoefficient( absorptionCoefficient, timeStep, directionalCoefficient );
            const REAL      divergenceCoefficient   = PML_divergenceCoefficient( density, c, directionalCoefficient );

            Tuple3i   cell_indices                  = _pressureField.cellIndex( cell_idx );
            int       neighbour_idx;

            // Scale the existing pressure value
            for ( int i = 0; i < _N; i++ )
                pDirectional(cell_idx, i) *= updateCoefficient;

            cell_indices[dimension] += 1;
            neighbour_idx = _velocityField[dimension].cellIndex(cell_indices);

            for ( int i = 0; i < _N; i++ )
                pDirectional(cell_idx, i) += divergenceCoefficient * v(neighbour_idx, i)/_velocityField[dimension].cellSize();

            cell_indices[dimension] -= 1;
            neighbour_idx = _velocityField[dimension].cellIndex(cell_indices);

            for ( int i = 0; i < _N; i++ )
                pDirectional(cell_idx, i) -= divergenceCoefficient * v(neighbour_idx, i)/_velocityField[dimension].cellSize();
        }
        else  // can collapse a lot of the terms 
        {
            // directionCoefficient = 1/dt
            // updateCoefficient    = 1
            // divergenceCoefficient= -rho c^2 *dt (cached)
            Tuple3i   cell_indices = _pressureField.cellIndex(cell_idx);
            int       neighbour_idx;

            cell_indices[dimension] += 1;
            neighbour_idx = _velocityField[dimension].cellIndex(cell_indices);

            for (int i=0; i<_N; i++)
                pDirectional(cell_idx, i) += n_rho_c_square_dt*v(neighbour_idx, i)*v_cellSize_inv;

            cell_indices[dimension] -= 1;
            neighbour_idx = _velocityField[dimension].cellIndex(cell_indices);

            for (int i=0; i<_N; i++)
                pDirectional(cell_idx, i) -= n_rho_c_square_dt*v(neighbour_idx, i)*v_cellSize_inv;

            // evaluate external sources only happens not in PML
            // Liu Eq (16) f6x term
            if (evaluateExternalSource)
                pDirectional(cell_idx, 0) += _objects->EvaluatePressureSources(cell_position, cell_position, simulationTime+0.5*timeStep)*timeStep;
                //for (int i = 0; i<_N; i++) 
                //    p(cell_idx, i) += (*sourceEvaluator)(cell_position, simulationTime+0.5*timeStep)*timeStep; 
        }
    }
}

void MAC_Grid::PML_pressureUpdateFull(const MATRIX *vArray, MATRIX &p, const REAL &timeStep, const REAL &c, const ExternalSourceEvaluator *sourceEvaluator, const REAL &simulationTime, const REAL &density )
{
    const int numPressureCells = _pressureField.numCells(); 
    //const int numPressureCells = _bulkCells.size(); 
    const bool evaluateExternalSource = _objects->HasExternalPressureSources();
    const REAL n_rho_c_square_dt = -density*c*c*timeStep;
    const REAL v_cellSize_inv[3] = {1./_velocityField[0].cellSize(),1./_velocityField[1].cellSize(),1./_velocityField[2].cellSize()};

    // We only have to worry about bulk cells here
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(static) default(shared)
    #endif
    for (int cell_idx = 0; cell_idx < numPressureCells; ++cell_idx)
    {
        if (!_isBulkCell[cell_idx])
            continue; 
        const Vector3d cell_position = _pressureField.cellPosition(cell_idx);
        //const Vector3d cell_position = _pressureField.cellPosition(_bulkCells[cell_idx]);

        for (int dimension=0; dimension<3; ++dimension) 
        {
            const REAL      absorptionCoefficient = PML_absorptionCoefficient(cell_position, _PML_absorptionWidth, dimension);
            const bool      inPML = (absorptionCoefficient > 1E-12) ? true : false; 
            const MATRIX   &v = vArray[dimension];

            if (inPML) 
            {
                const REAL      directionalCoefficient  = PML_directionalCoefficient( absorptionCoefficient, timeStep );
                const REAL      updateCoefficient       = PML_pressureUpdateCoefficient( absorptionCoefficient, timeStep, directionalCoefficient );
                const REAL      divergenceCoefficient   = PML_divergenceCoefficient( density, c, directionalCoefficient );

                Tuple3i   cell_indices                  = _pressureField.cellIndex(cell_idx);
                int       neighbour_idx;

                // Scale the existing pressure value
                for ( int i = 0; i < _N; i++ )
                    p(cell_idx, i) *= updateCoefficient; 

                cell_indices[dimension] += 1;
                neighbour_idx = _velocityField[dimension].cellIndex(cell_indices);

                for ( int i = 0; i < _N; i++ )
                    p(cell_idx, i) += divergenceCoefficient * v(neighbour_idx, i)/_velocityField[dimension].cellSize();

                cell_indices[dimension] -= 1;
                neighbour_idx = _velocityField[dimension].cellIndex(cell_indices);

                for ( int i = 0; i < _N; i++ )
                    p(cell_idx, i) -= divergenceCoefficient * v(neighbour_idx, i)/_velocityField[dimension].cellSize();
            }
            else  // can collapse a lot of the terms 
            {
                // directionCoefficient = 1/dt
                // updateCoefficient    = 1
                // divergenceCoefficient= -rho c^2 *dt (cached)
                Tuple3i   cell_indices = _pressureField.cellIndex(cell_idx);
                int       neighbour_idx;
                cell_indices[dimension] += 1;
                neighbour_idx = _velocityField[dimension].cellIndex(cell_indices);

                for (int i=0; i<_N; i++)
                    p(cell_idx, i) += n_rho_c_square_dt*v(neighbour_idx, i)*v_cellSize_inv[dimension];

                cell_indices[dimension] -= 1;
                neighbour_idx = _velocityField[dimension].cellIndex(cell_indices);

                for (int i=0; i<_N; i++)
                    p(cell_idx, i) -= n_rho_c_square_dt*v(neighbour_idx, i)*v_cellSize_inv[dimension];
                // evaluate external sources only happens not in PML
                // Liu Eq (16) f6x term
                if (evaluateExternalSource)
                    p(cell_idx,0) += _objects->EvaluatePressureSources(cell_position, cell_position, simulationTime+0.5*timeStep)*timeStep;
            }
        }
    }
}

void MAC_Grid::UpdatePMLPressure(MATRIX (&pDirectional)[3], MATRIX &pFull)
{
#if 0
    const int N_cells = numPressureCells(); 
    // update the full pressure in pml region
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(static) default(shared)
    #endif
    for (int cell_idx = 0; cell_idx < N_cells; ++cell_idx)
    {
        if (!_isBulkCell.at(cell_idx))
            continue; 

        const Vector3d  cell_position           = _pressureField.cellPosition( cell_idx );
        const REAL absorptionCoefficient_x   = PML_absorptionCoefficient( cell_position, _PML_absorptionWidth, 0);
        const REAL absorptionCoefficient_y   = PML_absorptionCoefficient( cell_position, _PML_absorptionWidth, 1);
        const REAL absorptionCoefficient_z   = PML_absorptionCoefficient( cell_position, _PML_absorptionWidth, 2);
        const bool inPML = (absorptionCoefficient_x > SMALL_NUM || 
                            absorptionCoefficient_y > SMALL_NUM || 
                            absorptionCoefficient_z > SMALL_NUM   ) ? true : false; 

        if (inPML) 
            pFull(cell_idx, 0) = pDirectional[0](cell_idx, 0) + pDirectional[1](cell_idx, 0) + pDirectional[2](cell_idx, 0); 
    }
#else
    pFull.parallelCopyAdd(pDirectional[0], pDirectional[1], pDirectional[2]);
#endif
}

void MAC_Grid::PML_pressureUpdateGhostCells(MATRIX &p, FloatArray &pGC, const REAL &timeStep, const REAL &c, const REAL &simulationTime, const REAL density)
{
    const int N_ghostCells = _ghostCellPositions.size(); 
    if (N_ghostCells == 0) return; 

#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int ghost_cell_idx=0; ghost_cell_idx<N_ghostCells; ++ghost_cell_idx)
    {
        const Vector3d &cellPosition = _ghostCellPositions.at(ghost_cell_idx); 
        int boundaryObject;
        REAL distance; 
        _objects->LowestObjectDistance(cellPosition, distance, boundaryObject); 

        // find reflection point, normal etc
        Vector3d boundaryPoint, imagePoint, erectedNormal; 
        auto object = _objects->GetPtr(boundaryObject); 
        object->ReflectAgainstBoundary(cellPosition, imagePoint, boundaryPoint, erectedNormal, distance); 
        const bool success = (_objects->LowestObjectDistance(imagePoint) >= DISTANCE_TOLERANCE); 
#pragma omp critical
        if (!success)
            std::cerr << "**WARNING** Reflection of ghost cell inside some objects: " << cellPosition << ". Proceed computation. \n"; 
        const REAL bcPressure = object->EvaluateBoundaryAcceleration(boundaryPoint, erectedNormal, simulationTime) * (-density); 
        const REAL weights = (object->DistanceToMesh(cellPosition) < DISTANCE_TOLERANCE ? -2.0*distance : -distance);  // finite-difference weight
        const REAL weightedPressure = bcPressure * weights; 

        // get the box enclosing the image point; 
        IntArray neighbours; 
        _pressureField.enclosingNeighbours(imagePoint, neighbours); 

        // figure out which dimension is this subdivied ghost cell TODO cache this
        const int cell_idx = _ghostCellParents.at(ghost_cell_idx); 
        int subdivideType = -1; 
        for (int ii=0; ii<6; ++ii)
        {
            if (_ghostCellsChildren.at(_ghostCellsInverse.at(cell_idx)).at(ii) == ghost_cell_idx)
            {
                subdivideType = ii;
                break;
            }
        }
        assert(subdivideType != -1); 
        const int subdivideDim = (int)floor(REAL(subdivideType)/2.0); 
        const int dir = (subdivideType % 2 == 0 ? -1 : 1); 

        REAL p_r; 
        if (object->GetOptionalAttributes().useRasterized) // use rasterized 
        {
            Tuple3i cellIndices = _pressureField.cellIndex(cell_idx); 
            cellIndices[subdivideDim] += dir; 
            p_r = p(_pressureField.cellIndex(cellIndices), 0); 
        }
        else // use MLS
        {
            T_MLS mls; 
            MLSPoint evalPt = Conversions::ToEigen(imagePoint); 
            std::vector<MLSPoint, P_ALLOCATOR> points; 
            std::vector<MLSVal, V_ALLOCATOR> attributes; 
            for (std::vector<int>::iterator it = neighbours.begin(); it != neighbours.end(); ) 
            {
                if (!_isBulkCell.at(*it))
                    it = neighbours.erase(it); 
                else
                {
                    const MLSPoint pt = Conversions::ToEigen(_pressureField.cellPosition(*it)); 
                    MLSVal val; 
                    val << p(*it, 0);
                    points.push_back(pt); 
                    attributes.push_back(val); 
                    ++it; 
                }
            }
            if (neighbours.size() >= 4)
            {
                const MLSVal mlsVal = mls.lookup(evalPt, points, attributes); 
                p_r = mlsVal(0, 0);
            }
            else
            {
                if (neighbours.size() != 0) // use closest neighbours
                {
                    std::vector<int>::iterator bestIt; 
                    REAL min_d2 = std::numeric_limits<REAL>::max();
                    for (std::vector<int>::iterator it=neighbours.begin(); it!=neighbours.end(); ++it)
                    {
                        const REAL d2 = (imagePoint - _pressureField.cellPosition(*it)).normSqr(); 
                        if (d2 < min_d2)
                        {
                            bestIt = it; 
                            min_d2 = d2; 
                        }
                    }
                    p_r = p(*bestIt, 0);
                }
                else // in this case, we don't have enough information, set this ghost cell to be its neighbour
                {
                    Tuple3i cellIndices = _pressureField.cellIndex(cell_idx); 
                    cellIndices[subdivideDim] += dir; 
                    p_r = p(_pressureField.cellIndex(cellIndices), 0); 
                }
            }
        }
        pGC.at(ghost_cell_idx) = p_r + weightedPressure; 
    }
}

//##############################################################################
// This function solves for the pressure ghost cells using Jacobi iteration. 
//
// NOTE: There are currently two cases this routine cannot properly resolve. Here
// I document the response of the system if these happen:
//  1) If reflection of a ghost cell ends up inside some other object, the solver
//     will simply print WARNING messages and proceed. The interpolant might be 
//     constructed using one or several solid cells, whose pressure values are 0.
//  2) If when constructing trilinear interpolant, there is one or more stencils 
//     are solid cells, 0 will be used to represent their values and WARNING
//     message will be printed. This treatment is consistent with problem 1).
//
// ..............................................................................
// The ghost-cell pressure update has several steps: 
//
//  1. Detect all ghost-cells, whose definition is cells that are solid but
//     has at least one fluid cell neighbours
//
//  2. For each ghost-cell GC, find the closest point on the boundary, BI, and
//     then compute an image point by extending BI, call this point IP. Note
//     that IP should be inside the fluid domain.
//
//     Ideally, GC-BI should be orthogonal to the level-set, in practice our
//     implementation does not guarantee that. Instead, for now I check the
//     orthogonality and throw an runtime error if angle is less than 80
//     degree. 
//
//  3. Compute interpolant to get pressure values at IP. Trilinear interpolant
//     is used if all eight neighbours are well-defined (uncoupled, not GC
//     itself). If one of the interpolation stencil is GC, replace a row of
//     Vandermonde matrix to incorporate the Neumann boundary condition sample.
//     Note that if the stencil is coupled to other ghost-cells, the inverse
//     of the Vandermonde can still be computed.
//
//     C = V^-1 p, C is the coefficient to the trilinear interpolant
//
//     p(x,y,z) = C^T b(x,y,z)
//
//     where 
//
//     b(x,y,z) = [xyz, xy, xz, yz, x, y, z, 1]^T
//
//  4. Finally we can compute beta's, such that 
//
//     p_IP = beta^T p_samples
//
//     p_samples are pressure values of eight neighbours. If i-th neighbour is
//     the GC cell itself, then the i-th component of p_samples is the Neumann
//     condition. 
//
//  5. We then enforce Neumann condition on the probe GC-BI-IP by
//     finite-difference. i.e., 
//
//     p_IP - p_GC = dl * dp/dn(BI)
//
//     dl is the length from GC-IP. 
//
//  6. Because the ghost-cells might be coupled, each evaluation of 5. is put
//     into a sparse matrix, and solved outside the loop. 
//
//  7. Update the pressure at ghost-cell indices. 
//
//##############################################################################
void MAC_Grid::PML_pressureUpdateGhostCells_Coupled( MATRIX &p, FloatArray &pGC, const REAL &timeStep, const REAL &c, const REAL &simulationTime, const REAL density)
{
    const int N_ghostCells = _ghostCellPositions.size(); 
    if (N_ghostCells == 0) return; 
    _ghostCellCoupledData.clear(); 
    _ghostCellCoupledData.resize(N_ghostCells); 
    int replacedCell = 0;

#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int ghost_cell_idx=0; ghost_cell_idx<N_ghostCells; ++ghost_cell_idx)
    {
        const int gcParentIndex = _ghostCellParents.at(ghost_cell_idx); 
        const Vector3d &cellPosition = _ghostCellPositions.at(ghost_cell_idx); 
        const int boundaryObject = _ghostCellBoundaryIDs.at(ghost_cell_idx); 
        //const int       cellIndex      = _ghostCells[ghost_cell_idx]; 
        //const Tuple3i   cellIndices    = _pressureField.cellIndex(cellIndex); 
        //const Vector3d  cellPosition   = _pressureField.cellPosition(cellIndices); 
        //const int       boundaryObject = _containingObject[cellIndex]; 

        // find point BI and IP in the formulation and then query boundary
        // condition
        Vector3d boundaryPoint, imagePoint, erectedNormal; 
        REAL distance; 
        auto object = _objects->GetPtr(boundaryObject); 
        object->ReflectAgainstBoundary(cellPosition, imagePoint, boundaryPoint, erectedNormal, distance); 
        const bool success = (_objects->LowestObjectDistance(imagePoint) >= DISTANCE_TOLERANCE); 

        if (!success)
            std::cerr << "**WARNING** Reflection of ghost cell inside some objects: " << cellPosition << ". Proceed computation. \n"; 
            //throw std::runtime_error("**ERROR** Ghost cell reflection still inside some objects. This is not handled."); 

        const REAL bcPressure = object->EvaluateBoundaryAcceleration(boundaryPoint, erectedNormal, simulationTime) * (-density); 
        const REAL weights = (distance >=0 ? -distance : 2.0*distance);  // finite-difference weight
        const REAL weightedPressure = bcPressure * weights; 

        // get the box enclosing the image point; 
        IntArray neighbours; 
        _pressureField.enclosingNeighbours(imagePoint, neighbours); 

        // hasGC  : has self as interpolation stencil
        int hasGC=-1; 

        // vandermonde matrix, see 2008 Mittals JCP paper Eq.18
        Eigen::MatrixXd V(8,8); 
        Tuple3i     indicesBuffer;
        Vector3d    positionBuffer; 
        Eigen::VectorXd pressureNeighbours(8);  // right hand side

        // some coupling information
        std::vector<int> coupledGhostCells; 
        std::vector<int> uncoupledGhostCellsNeighbours; // [0,neighbours.size)
        std::vector<int> coupledGhostCellsNeighbours; // [0,neighbours.size)

        int count_solid=0;
        for (size_t ii=0; ii<neighbours.size(); ii++) 
        {
            if (neighbours[ii] == gcParentIndex) 
            {
                hasGC = ii;
                FillVandermondeBoundary(ii, boundaryPoint, erectedNormal, V);
                pressureNeighbours(ii) = bcPressure; 
            }
            else  
            {
                // this cell contains a ghost cell.
                if (_isGhostCell[neighbours[ii]]) 
                {
                    // access its children and find the closest children index
                    const int gcArrayPosition = _ghostCellsInverse.at(neighbours[ii]); 
                    const IntArray &gcChildren = _ghostCellsChildren.at(gcArrayPosition); 

                    const int N_children = gcChildren.size(); 
                    int bestChild = -1; 
                    REAL bestChildDistance = std::numeric_limits<REAL>::max(); 
                    for (int c_idx=0; c_idx<N_children; ++c_idx) 
                    {
                        const int child = gcChildren.at(c_idx); 
                        if (child == -1) // child does not exist, skip
                            continue; 
                        const Vector3d &gcChildrenPosition = _ghostCellPositions.at(child); 
                        const REAL distanceSqr = (gcChildrenPosition - imagePoint).lengthSqr(); 
                        if (distanceSqr < bestChildDistance)
                        {
                            bestChild = child; 
                            bestChildDistance = distanceSqr; 
                        }
                    }
                    assert(bestChild != -1); 

                    coupledGhostCells.push_back(bestChild); 
                    //coupledGhostCells.push_back(_ghostCellsInverse.at(neighbours[ii])); 
                    coupledGhostCellsNeighbours.push_back(ii); 
                    FillVandermondeRegular(ii, _ghostCellPositions.at(bestChild), V); 
                    pressureNeighbours(ii) = pGC.at(bestChild); 
                }
                else
                {
                    uncoupledGhostCellsNeighbours.push_back(ii);

                    // use the regular pressure field position to fill
                    // vandermonde row
                    indicesBuffer = _pressureField.cellIndex(neighbours.at(ii)); 
                    positionBuffer= _pressureField.cellPosition(indicesBuffer); 
                    FillVandermondeRegular(ii, positionBuffer, V);
                    pressureNeighbours(ii) = p(neighbours.at(ii), 0); 
                    if (IsPressureCellSolid(neighbours[ii])) 
                        count_solid ++; 
                }
            }
        }

        if (count_solid > 0) // This case is currently handled by simply setting the contribution from solid cell to be zero
            std::cerr << "**WARNING** When constructing trilinear interpolant, " << count_solid << " solid cells have been used, whose pressure value = 0\n"; 

        // accumulate counter
#ifdef USE_OPENMP
#pragma omp critical
#endif
        if (hasGC != -1)
            replacedCell ++; 

        // want beta = V^-T b, equivalent to solving V^T beta = b
        Eigen::MatrixXd b(1,8); // row vector

        FillVandermondeRegular(0, imagePoint, b); 
        b.transposeInPlace();  // column vector
        V.transposeInPlace(); 

        // Linear solve: lu should be fastest, but in practice its marginally
        // faster in release mode. other options are householder qr and svd. 
        // condition number of the system is only available if svd is used.
        Eigen::VectorXd beta = V.partialPivLu().solve(b); 
        //Eigen::VectorXd beta = V.householderQr().solve(b); 
        //Eigen::JacobiSVD<Eigen::MatrixXd> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
        //Eigen::VectorXd beta = svd.solve(b); 
//        const double conditionNumber = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
//#ifdef USE_OPENMP
//#pragma omp critical
//#endif
//        if (conditionNumber > 1E10) 
//        {
//            std::cout << "**WARNING** condition number for the least square solve is = " << conditionNumber << "\n"
//                      << "            the solution can be inaccurate.\n"
//                      << "largest  singular value = " << svd.singularValues()(0) << "\n"
//                      << "smallest singular value = " << svd.singularValues()(svd.singularValues().size()-1) << "\n"
//                      << "problematic matrix: \n"
//                      << V << std::endl;
//        }

        // forms the sparse matrix
        JacobiIterationData &jacobiIterationData = _ghostCellCoupledData[ghost_cell_idx]; 
        jacobiIterationData.cellId = ghost_cell_idx; 
        for (size_t cc=0; cc<coupledGhostCells.size(); cc++) 
        {
            // only stores entries that is significant
            if (fabs(beta(coupledGhostCellsNeighbours[cc])) > GHOST_CELL_ENTRY_THRESHOLD)
            {
                jacobiIterationData.nnzIndex.push_back(coupledGhostCells[cc]); 
                jacobiIterationData.nnzValue.push_back(-beta(coupledGhostCellsNeighbours[cc])); 

                // log couple count to cell. this info is symmetric
                _ghostCellCoupledData.at(ghost_cell_idx).coupleCount ++; 
                _ghostCellCoupledData.at(coupledGhostCells[cc]).coupleCount ++; 
            }
        } 

        // RHS always exists even when no off-diagonal elements 
        double &RHS = jacobiIterationData.RHS; 
        RHS = weightedPressure;
        for (size_t uc=0; uc<uncoupledGhostCellsNeighbours.size(); uc++) 
            RHS += beta(uncoupledGhostCellsNeighbours[uc])*pressureNeighbours(uncoupledGhostCellsNeighbours[uc]); 

        // if GC itself is coupled, need to add this to the RHS
        if (hasGC != -1) 
            RHS += beta(hasGC)*pressureNeighbours(hasGC); 
    }

    // analyze sparsity pattern and deflate the linear system
    int maxCoupleCount = std::numeric_limits<int>::min(); 
    std::vector<int> coupledIndices; 
    for (int g_idx=0; g_idx<N_ghostCells; ++g_idx) 
    {
        const auto &data = _ghostCellCoupledData.at(g_idx); 
        maxCoupleCount = max(maxCoupleCount, data.coupleCount); 
        if (data.coupleCount == 0) 
            pGC.at(data.cellId) = data.RHS; 
        else
            coupledIndices.push_back(g_idx); 
    }
#ifdef DEBUG_PRINT
    std::cout << "\nGhost cell linear system analysis result: \n"; 
    std::cout << " Max cell couple count     : " << maxCoupleCount << std::endl;
    std::cout << " Decoupled cell percentage : " << (REAL)(N_ghostCells - coupledIndices.size()) / (REAL)N_ghostCells << "\n\n";
#endif

    ExamineJacobiMatrix();

    // solve the linear system
    const int maxIteration = GHOST_CELL_JACOBI_MAX_ITERATION;
    const int N_coupled = coupledIndices.size(); 
    FloatArray pGC_old; 
    REAL minResidual, maxResidual, maxOffDiagonal, oldResidual = -1; 
    int maxResidualEntry = -1; 
    for (int iteration=0; iteration<maxIteration; iteration++) 
    {
        // compute and print the residual
        ComputeGhostCellSolveResidual(pGC, minResidual, maxResidual, maxResidualEntry, maxOffDiagonal); 
        const REAL absResidual = max(fabs(minResidual), fabs(maxResidual)); 
        // early termination
        if (EQUAL_FLOATS(oldResidual, absResidual) || absResidual < 1E-8)
            break; 

        oldResidual = absResidual; 
        std::cout << "  iteration " << iteration << ": residual = " << std::setprecision(16) << absResidual << std::endl;

        // start iteration
        pGC_old = pGC; 
        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(static) default(shared)
        #endif
        for (int c_idx=0; c_idx<N_coupled; ++c_idx) 
        {
            const int ghost_cell_idx = coupledIndices.at(c_idx); 
            const JacobiIterationData &data = _ghostCellCoupledData.at(ghost_cell_idx); 
            // if there are nnz, update them 
            //p(_ghostCells[ghost_cell_idx],0) = data.RHS; 
            pGC.at(ghost_cell_idx) = data.RHS; 
            for (size_t cc=0; cc<data.nnzIndex.size(); cc++)
                pGC.at(ghost_cell_idx) -= data.nnzValue[cc]*pGC_old.at(data.nnzIndex[cc]); 
        }
    }
    std::cout << "  max residual = " << maxResidual << "; entry = " << maxResidualEntry << "; cell position = " << _ghostCellPositions.at(maxResidualEntry) << std::endl;
}

//##############################################################################
//##############################################################################
void MAC_Grid::
UpdateGhostCells_FV(MATRIX &p, const REAL &simulationTime)
{
    typedef std::unordered_map<int, std::shared_ptr<GhostCell> >::iterator Iterator_GC;
    typedef std::vector<GhostCell::BoundarySamples>::iterator Iterator_BS;
    typedef std::vector<GhostCell::VolumeSamples>::iterator Iterator_VS;
    typedef std::vector<TriangleIdentifier>::iterator         Iterator_TS;
    const int  &N_sub    = _waveSolverSettings->FV_boundarySubdivision; 
    const int   N_sub_2  = pow(N_sub, 2);
    const int   N_sub_3  = pow(N_sub, 3);
    const REAL &c        = _waveSolverSettings->soundSpeed; 
    const REAL &k        = _waveSolverSettings->timeStepSize; 
    const REAL &h        = _waveSolverSettings->cellSize; 
    const REAL &rho      = _waveSolverSettings->airDensity;
    const int  &cell_div = _waveSolverSettings->cellDivisions; 
    const REAL  h_sub    = h/(REAL)N_sub; 
    const REAL  h_over_2 = h/2.0;
    const REAL  h_sub_over_2 = h_sub/2.0; 
    const REAL  A_sub    = pow(h_sub, 2);
    const REAL  V_sub    = pow(h_sub, 3);
    const REAL  kc_2     = pow(k*c, 2); 
    const size_t N_totalBoundarySamples = 6*N_sub_2; 
    const size_t N_totalVolumeSamples = N_sub_3; 
    const int N_gc = _ghostCellsCollection.size();

#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int gc_i=0; gc_i<N_gc; ++gc_i)
    {
        auto gcm=_ghostCellsCollection.begin(); 
        std::advance(gcm, gc_i);
        auto &gc = gcm->second; 
        const Vector3d cellPosition = _pressureField.cellPosition(gc->parent_idx); 
        Tuple3i indicesBuffer = _pressureField.cellIndex(gc->parent_idx); 
        if (gc->boundarySamples.size() != N_totalBoundarySamples) // need to generate boundary samples
        {
            gc->boundarySamples.resize(N_totalBoundarySamples);
            for (int dim_0=0; dim_0<3; ++dim_0)
            {
                const int dim_1 = (dim_0+1)%3; 
                const int dim_2 = (dim_0+2)%3; 
                const REAL offset_dim_0[2] = {cellPosition[dim_0] - h_over_2, cellPosition[dim_0] + h_over_2}; 
                const int start_idx[2] = {(dim_0*2)*N_sub_2, (dim_0*2+1)*N_sub_2}; 
                const REAL lowCorner_dim_1 = cellPosition[dim_1] - h_over_2; 
                const REAL lowCorner_dim_2 = cellPosition[dim_2] - h_over_2; 
                for (int np=0; np<2; ++np) // negative/positive face
                {
                    Vector3d normal(0,0,0);
                    normal[dim_0] = (np==0 ? -1.0 : 1.0) * A_sub; 
                    indicesBuffer[dim_0] += (np==0 ? -1 : 1); 
                    assert(indicesBuffer[dim_0]>=0 && indicesBuffer[dim_0]<cell_div); 
                    const int neighbour_idx = _pressureField.cellIndex(indicesBuffer); 
                    const int gc_value_idx = (neighbour_idx > gc->parent_idx ? dim_0*2+1 : dim_0*2); 
                    indicesBuffer[dim_0] -= (np==0 ? -1 : 1); // reset buffer
                    for (int ii=0; ii<N_sub; ++ii) 
                    {
                        for (int jj=0; jj<N_sub; ++jj) 
                        {
                            Vector3d position;
                            position[dim_1] = lowCorner_dim_1 + h_sub_over_2 + (REAL)ii*h_sub; 
                            position[dim_2] = lowCorner_dim_2 + h_sub_over_2 + (REAL)jj*h_sub; 
                            position[dim_0] = offset_dim_0[np]; 
                            gc->boundarySamples.at(start_idx[np] + ii*N_sub + jj) = GhostCell::BoundarySamples(position, normal, neighbour_idx, gc_value_idx); 
                        }
                    }
                }
            }
        }
        if (gc->volumeSamples.size() != N_totalVolumeSamples) // need to generate volume samples
        {
            gc->volumeSamples.resize(N_totalVolumeSamples);
            const REAL lowCorner_x = cellPosition[0] - h_over_2; 
            const REAL lowCorner_y = cellPosition[1] - h_over_2; 
            const REAL lowCorner_z = cellPosition[2] - h_over_2; 
            for (int ii=0; ii<N_sub; ++ii)
            {
                for (int jj=0; jj<N_sub; ++jj)
                {
                    for (int kk=0; kk<N_sub; ++kk)
                    {
                        Vector3d position; 
                        position[0] = lowCorner_x + h_sub_over_2 + (REAL)ii*h_sub; 
                        position[1] = lowCorner_y + h_sub_over_2 + (REAL)jj*h_sub; 
                        position[2] = lowCorner_z + h_sub_over_2 + (REAL)kk*h_sub; 
                        gc->volumeSamples.at(ii*N_sub_2 + jj*N_sub + kk) = GhostCell::VolumeSamples(position);
                    }
                }
            }
        }
        // compute volume and pressure gradient
        REAL &volume = gc->volume; 
        REAL &dp_dn_dot_S = gc->dp_dn_dot_S; 
        volume=0.0; 
        dp_dn_dot_S=0.0;
        for (Iterator_VS sp=gc->volumeSamples.begin(); sp!=gc->volumeSamples.end(); ++sp)
        {
            sp->isBulk = (_objects->LowestObjectDistance(sp->position) < DISTANCE_TOLERANCE ? false : true);
            if (sp->isBulk)
                volume += V_sub;
        }
        for (Iterator_BS sp=gc->boundarySamples.begin(); sp!=gc->boundarySamples.end(); ++sp)
        {
            sp->isBulk = (_objects->LowestObjectDistance(sp->position) < DISTANCE_TOLERANCE ? false : true);
            if (sp->isBulk && !IsPressureCellSolid(sp->neighbour_idx))
            {
                REAL p_neighbour; 
                if (!_isGhostCell.at(sp->neighbour_idx))
                {
                    p_neighbour = p(sp->neighbour_idx, 0); 
                }
                else
                {
                    const std::shared_ptr<GhostCell> &neighbour_gc = _ghostCellsCollection.at(sp->neighbour_idx); 
                    const int neighbour_gc_value_idx = (sp->gc_value_idx%2==0 ? sp->gc_value_idx+1 : sp->gc_value_idx-1); 
                    if (neighbour_gc->validSample)
                        p_neighbour = neighbour_gc->values.at(neighbour_gc->valuePointer).at(neighbour_gc_value_idx); 
                    else
                        p_neighbour = gc->values.at(gc->valuePointer).at(sp->gc_value_idx); 
                }
                
                const REAL dp_dn_face = (sp->neighbour_idx > gc->parent_idx ? (p_neighbour - gc->values.at(gc->valuePointer).at(sp->gc_value_idx))/h
                                                                            : (p_neighbour - gc->values.at(gc->valuePointer).at(sp->gc_value_idx))/h);
                //volume += ((sp->position-cellPosition)/3.0).dotProduct(sp->normal);
                if (!IsPressureCellSolid(sp->neighbour_idx))
                    dp_dn_dot_S += dp_dn_face; 
            }
        }
        dp_dn_dot_S *= A_sub; // scale by boundary sample area
        for (Iterator_TS sp=gc->hashedTriangles->begin(); sp!=gc->hashedTriangles->end(); ++sp)
        {
            const auto &object = _objects->GetPtr(sp->objectID); 
            const auto &mesh   = object->GetMeshPtr();
            const Tuple3ui &triangle_vtx = mesh->triangle_ids(sp->triangleID);
            const REAL triangle_area = mesh->triangle_normal(sp->triangleID).norm(); 
            REAL dp_dn_e = 0.0; // sample three vertices 
            for (int tri_corner=0; tri_corner<3; ++tri_corner)
            {
                const int &tri_vtx = triangle_vtx[tri_corner]; 
                const Point3<REAL> vtx = object->ObjectToWorldPoint(mesh->vertex(tri_vtx)); 
                if (_objects->LowestObjectDistance(vtx) >= DISTANCE_TOLERANCE)
                {
                    const Vector3d normal = object->ObjectToWorldVector(mesh->normal(tri_vtx)); 
                    //volume -= ((vtx-cellPosition)/3.0).dotProduct(normal)/3.0*triangle_area; // flip normal; 3.0 for boundary integral; 3.0 for vertex avg
                    dp_dn_e += object->EvaluateBoundaryAcceleration(triangle_vtx[tri_corner], normal, simulationTime);
                }
            }
            dp_dn_e *= (-rho)/3.0; // average and convert acc to dp_dn
            dp_dn_dot_S -= dp_dn_e * triangle_area; // same normal flip as volume computation
        }
        if (volume < SMALL_NUM)
        {
            volume = SMALL_NUM; // clamp
            gc->validSample = false; 
            dp_dn_dot_S = 0;
        }
        else 
        {
            gc->validSample = true;
        }

        // update all ghost pressure samples; keep all sample values the same for now
        for (int sp=0; sp<6; ++sp)
        {
            const REAL &gc_p_last = gc->values.at((gc->valuePointer+2)%3).at(sp);
            const REAL &gc_p_curr = gc->values.at((gc->valuePointer  )%3).at(sp);
                  REAL &gc_p_next = gc->values.at((gc->valuePointer+1)%3).at(sp);
            if (volume > SMALL_NUM)
                gc_p_next = 2.0 * gc_p_curr - gc_p_last + kc_2*(dp_dn_dot_S/volume); 
            else
                gc_p_next = 0.0; 
        }
    }
    GhostCell::valuePointer = (GhostCell::valuePointer+1)%3;
}

void MAC_Grid::sampleZSlice( int slice, const MATRIX &p, MATRIX &sliceData )
{
    const Tuple3i             &divs = _pressureField.cellDivisions();

    if ( sliceData.rows() != divs[ 1 ] || sliceData.cols() != divs[ 0 ] )
    {
        sliceData.resizeAndWipe( divs[ 1 ], divs[ 0 ] );
    }

    for ( int i = 0; i < divs[ 0 ]; i++ )
        for ( int j = 0; j < divs[ 1 ]; j++ )
        {
            int                      cell_idx;

            cell_idx = _pressureField.cellIndex( Tuple3i( i, j, slice ) );

            // Sample just the first field
            sliceData( divs[ 1 ] - j - 1, i ) = p( cell_idx, 0 );
        }
}

void MAC_Grid::SampleAxisAlignedSlice(const int &dim, const REAL &offset, MATRIX const (&pDirectional)[3], const MATRIX &pFull, const FloatArray &pGC, const MATRIX (&v)[3], std::vector<Cell> &sampledCells) const
{
    const Tuple3i &divs = _pressureField.cellDivisions();
    const BoundingBox &bbox = _velocityField[dim].bbox(); 
    const REAL lowerCorner = bbox.minBound()[dim]; 
    const int dim1 = (dim+1)%3;
    const int dim2 = (dim+2)%3; 
    sampledCells.resize(divs[dim1] * divs[dim2]); 
    const int k = std::min<int>((int)((offset - lowerCorner)/_waveSolverSettings->cellSize), _waveSolverSettings->cellDivisions-1); 

    int count = 0;
    for (int i=0; i<divs[dim1]; ++i)
        for (int j = 0; j<divs[dim2]; ++j)
        {
            Tuple3i indices; 
            indices[dim ] = k; 
            indices[dim1] = i; 
            indices[dim2] = j; 
            const int cell_idx = _pressureField.cellIndex(indices); 
            GetCell(cell_idx, pDirectional, pFull, pGC, v, sampledCells.at(count)); 
            count++; 
        }
}

void MAC_Grid::SmoothFieldInplace(MATRIX &p1, MATRIX &p2, MATRIX &p3, REAL w1, REAL w2, REAL w3)
{
    p3 *= w3; 
    p3.parallelAxpy( w1, p1 ); 
    p3.parallelAxpy( w2, p2 );
}

void MAC_Grid::ComputeGhostCellInverseMap()
{
    _ghostCellsInverse.clear(); 
    for (size_t ii=0; ii<_ghostCells.size(); ii++)
        _ghostCellsInverse[_ghostCells[ii]] = ii; 
}

//##############################################################################
// This function interpolate the fresh pressure cell
//
// Note: 
//   1) fresh pressure cell will have to be outside the object, otherwise the last
//   step of the interpolation won't work: Need BI-CU-IP so CU is well-defined
//   for the interpolant. 
//
//   2) if the interpolation stencil does not have valid history and is the
//      cell itself, then throw exception because this should not happen, but if 
//      the stencil does not have valid history and is caused by other cells, 
//      zero will be used in the current implementation.
//
//##############################################################################
void MAC_Grid::InterpolateFreshPressureCell(MATRIX &p, const REAL &timeStep, const REAL &simulationTime, const REAL &density)
{
    const int N = _pressureField.numCells(); 

    int count = 0;
    for (int cell_idx = 0; cell_idx < N; ++cell_idx)
    {
        if (!_pressureCellHasValidHistory.at(cell_idx))
            count ++; 
    }

    for (int cell_idx = 0; cell_idx < N; ++cell_idx)
    {
        // if has valid history or if its not bulk, then no interpolation needed
        if (_pressureCellHasValidHistory.at(cell_idx) || !_isBulkCell.at(cell_idx))
            continue; 

        const Vector3d cellPosition = _pressureField.cellPosition(cell_idx);
        REAL distance; 
        int objectID; 
        _objects->LowestObjectDistance(cellPosition, distance, objectID); 
        
        if (distance < DISTANCE_TOLERANCE)
            throw std::runtime_error("**ERROR** Fresh cell inside some object. This shouldn't happen for pressure cells. Distance: " + std::to_string(distance)); 
        Vector3d imagePoint, boundaryPoint, erectedNormal; 
        REAL distanceTravelled; 
        _objects->Get(objectID).ReflectAgainstBoundary(cellPosition, imagePoint, boundaryPoint, erectedNormal, distanceTravelled); 

        // prepare interpolation stencils, this part is similar to
        // the vandermonde part in ghost cell pressure update
        IntArray neighbours; 
        _pressureField.enclosingNeighbours(imagePoint, neighbours); 

#if 0
        Eigen::MatrixXd V(8,8); 
        Tuple3i indicesBuffer;
        Vector3d positionBuffer; 
        Eigen::VectorXd RHS(8); 
        // should evaluate the pressure at the previous time-step
        const REAL boundaryPressureGradientNormal = _objects->Get(objectID).EvaluateBoundaryAcceleration(boundaryPoint, erectedNormal, simulationTime + 0.5*timeStep)*(-density); 
        //const REAL boundaryPressureGradientNormal = _objects->Get(objectID).EvaluateBoundaryAcceleration(boundaryPoint, erectedNormal, simulationTime)*(-density); 

        for (size_t row=0; row<neighbours.size(); row++) 
        {
            const int neighbour_idx = neighbours[row]; 
            const bool neighbourValid = _pressureCellHasValidHistory.at(neighbour_idx) ? true : false; 

            indicesBuffer = _pressureField.cellIndex(neighbour_idx); 
            positionBuffer= _pressureField.cellPosition(indicesBuffer); 
            if (neighbours[row] != cell_idx) // not self
            {
                FillVandermondeRegular(row, positionBuffer, V);
                if (neighbourValid)
                {
                    RHS(row) = p(neighbours[row], 0); 
                }
                else
                {
                    if (_containingObject.at(neighbour_idx) == objectID)
                    {
                        if (_waveSolverSettings->useGhostCell) // this shouldn't happen for ghost cells
                        {
                            std::cout << "invalid cell has all neighbours = \n";
                            STL_Wrapper::PrintVectorContent(std::cout, neighbours); 
                            throw std::runtime_error("**ERROR** Invalid self pressure cells encountered in interpolation."); 
                        }
                        else // using rasterized boundary, give it 0 for now FIXME: better way?
                        {
                            RHS(row) = 0.0;
                        }
                    }
                    else
                    {
                        RHS(row) = 0.0; 
                    }
                }
            }
            else 
            {
                FillVandermondeBoundary(row, boundaryPoint, erectedNormal, V);
                RHS(row) = boundaryPressureGradientNormal; 
            }
        }

        // coefficient for the interpolant
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::VectorXd C = svd.solve(RHS); 

        //// evaluate the interpolant at cell position
        Eigen::VectorXd coordinateVector(8); 
        FillVandermondeRegular(cellPosition, coordinateVector);
        const REAL pressureFreshCell = C.dot(coordinateVector); 
        p(cell_idx, 0) = pressureFreshCell; 

#else
        // remove invalid neighbours and ready the MLS interpolator
        T_MLS mls; 
        MLSPoint evalPt = Conversions::ToEigen(imagePoint); 
        std::vector<MLSPoint, P_ALLOCATOR> points; 
        std::vector<MLSVal, V_ALLOCATOR> attributes; 
        for (std::vector<int>::iterator it = neighbours.begin(); it != neighbours.end(); ) 
        {
            if (!_pressureCellHasValidHistory.at(*it))
                it = neighbours.erase(it); 
            else
            {
                const MLSPoint pt = Conversions::ToEigen(_pressureField.cellPosition(*it)); 
                MLSVal val; 
                val << p(*it, 0);
                points.push_back(pt); 
                attributes.push_back(val); 
                ++it; 
            }
        }
        if (neighbours.size() < 4)
        {
            if (neighbours.size() != 0) // worst case: use one of its neighbours
            {
                p(cell_idx, 0) = p(neighbours.at(0), 0); 
            }
            else // in this case, we just use nearest neighbours' pressure
            {
                Tuple3i indices = _pressureField.cellIndex(cell_idx); 
                int ind_i = -1, bestInd = -1; 
                for (int dim=0; dim<3; ++dim)
                {
                    indices[dim] += 1; 
                    ind_i = _pressureField.cellIndex(indices); 
                    if (_isBulkCell.at(ind_i))
                    {
                        bestInd = ind_i; 
                        break; 
                    }
                    indices[dim] -= 2; 
                    ind_i = _pressureField.cellIndex(indices); 
                    if (_isBulkCell.at(ind_i))
                    {
                        bestInd = ind_i; 
                        break;
                    }
                    indices[dim] += 1; 
                }
                if (bestInd != -1)
                    p(cell_idx, 0) = p(bestInd, 0); 
                else 
                    p(cell_idx, 0) = 0.0;
            }
        }
        else
        {
            const MLSVal mlsVal = mls.lookup(evalPt, points, attributes); 
            p(cell_idx, 0) = mlsVal(0, 0);
        }
#endif
    }
}

//##############################################################################
// Interpolate the fresh velocity cell
//
// Note: 
//   fresh velocity cell can be inside the object, so make sure the
//   interpolation procedure is robust
//
//   Only performed on interfacial cells, assuming solid cell cannot turn to
//   bulk directly without being interfacial first. 
//
//   If neighbour has invalid history, grab the nearest boundary condition to
//   that neighbour and use it as its velocity value.
//##############################################################################
void MAC_Grid::InterpolateFreshVelocityCell(MATRIX &v, const int &dim, const REAL &timeStep, const REAL &simulationTime)
{
    const int N = _velocityInterfacialCells[dim].size(); 
    // essentially same procedures as pressure cells
    for (int interfacial_cell_idx = 0; interfacial_cell_idx < N; ++interfacial_cell_idx) 
    {
        const int cell_idx = _velocityInterfacialCells[dim].at(interfacial_cell_idx);
        // if has valid history or if its solid, then no interpolation needed
        if (_velocityCellHasValidHistory[dim].at(cell_idx) || IsVelocityCellSolid(cell_idx, dim))
            continue; 
        const Vector3d cellPosition = _velocityField[dim].cellPosition(cell_idx);
        REAL distance; 
        int objectID; 
        _objects->LowestObjectDistance(cellPosition, distance, objectID); 
        Vector3d imagePoint, boundaryPoint, erectedNormal; 
        REAL dBuffer;
        _objects->Get(objectID).ReflectAgainstBoundary(cellPosition, imagePoint, boundaryPoint, erectedNormal, dBuffer); 

        // prepare interpolation stencils, this part is similar to
        // the vandermonde part in ghost cell pressure update
        IntArray neighbours; 
        Tuple3i indicesBuffer;
        Vector3d positionBuffer; 

        // interpolate using trilinear interpolant
        // boundary velocity will be scaled by the normal difference
        _velocityField[dim].enclosingNeighbours(imagePoint, neighbours); 

#if 0 // Mittal
        // evaluate the velocity source in the previous step for velocity
        Vector3d normal(0, 0, 0); 
        normal[dim] = 1.0;
        const REAL boundaryVelocity = _objects->Get(objectID).EvaluateBoundaryVelocity(boundaryPoint, normal, simulationTime - 0.5*timeStep); 
        Eigen::MatrixXd V(8,8); 
        Eigen::VectorXd RHS(8); 
        for (size_t row=0; row<neighbours.size(); row++) 
        {
            const int neighbour_idx = neighbours[row]; 
            const bool neighbourValid = _velocityCellHasValidHistory[dim].at(neighbour_idx) ? true : false; 

            //if (!_velocityCellHasValidHistory[dim].at(neighbour_idx) && neighbour_idx != cell_idx) 
            //    throw std::runtime_error("**ERROR** one of the interpolation stencil for velocity fresh cell has invalid history."); 

            indicesBuffer = _velocityField[dim].cellIndex(neighbours.at(row)); 
            positionBuffer= _velocityField[dim].cellPosition(indicesBuffer); 
            if (neighbours.at(row) != cell_idx) // not self
            {
                FillVandermondeRegular(row, positionBuffer, V);
                if (neighbourValid)
                {
                    RHS(row) = v(neighbours.at(row), 0); 
                }
                else // if not valid, try to get the closest boundary point and use the boundary condition there.
                {
                    const Vector3d cellPositionNeighbour = _velocityField[dim].cellPosition(neighbour_idx);
                    REAL distanceNeighbour; 
                    Vector3d imagePointNeighbour, boundaryPointNeighbour, erectedNormalNeighbour; 
                    _objects->Get(objectID).ReflectAgainstBoundary(cellPositionNeighbour, imagePointNeighbour, boundaryPointNeighbour, erectedNormalNeighbour, distanceNeighbour); 
                    Vector3d gradient; 
                    _objects->ObjectNormal(objectID, boundaryPointNeighbour, gradient);
                    gradient.normalize(); 
                    const REAL coefficient = gradient(dim); 
                    const REAL boundaryVelocityNeighbour = _objects->Get(objectID).EvaluateBoundaryVelocity(boundaryPointNeighbour, erectedNormalNeighbour, simulationTime - 0.5*timeStep) * coefficient;
                    RHS(row) = boundaryVelocityNeighbour; 
                }
            }
            else 
            {
                FillVandermondeRegular(row, boundaryPoint, V); 
                RHS(row) = boundaryVelocity; 
            }
        }
        // coefficient for the interpolant
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::VectorXd C = svd.solve(RHS); 

        Eigen::VectorXd coordinateVector(8); 
        FillVandermondeRegular(imagePoint, coordinateVector);
        const REAL v_IP = C.dot(coordinateVector); 
        const REAL v_BI = boundaryVelocity;
        const REAL distanceToMesh = _objects->GetPtr(objectID)->DistanceToMesh(imagePoint); 
        if (distanceToMesh > DISTANCE_TOLERANCE) // cell position is located outside the boundary: the erected normal is BI-CU-IP
            v(cell_idx, 0) = 0.5 * (v_IP + v_BI); 
        else // cell position is located inside the boundary: the erected normal is CU-BI-IP
            v(cell_idx, 0) = 2.0 * v_BI - v_IP;
#endif 
#if 0 // MLS
        // remove invalid neighbours and ready the MLS interpolator
        T_MLS mls; 
        MLSPoint evalPt = Conversions::ToEigen(imagePoint); 
        std::vector<MLSPoint, P_ALLOCATOR> points; 
        std::vector<MLSVal, V_ALLOCATOR> attributes; 
        for (std::vector<int>::iterator it = neighbours.begin(); it != neighbours.end(); ) 
        {
            if (!_velocityCellHasValidHistory[dim].at(*it))
                it = neighbours.erase(it); 
            else
            {
                MLSPoint pt = Conversions::ToEigen(_velocityField[dim].cellPosition(*it)); 
                MLSVal val = MLSVal(v(*it, 0)); 
                points.push_back(pt); 
                attributes.push_back(val); 
                ++it; 
            }
        }
        if (neighbours.size() < 4)
        {
            throw std::runtime_error("**ERROR** Velocity interpolation error: cannot construct interpolant");
        }
        else
        {
            MLSVal mlsVal = mls.lookup(evalPt, points, attributes); 
            v(cell_idx, 0) = mlsVal(0, 0);
        }
        std::cout << "v dim = " << dim << std::endl;
        COUT_SDUMP(cellPosition);
        COUT_SDUMP(v(cell_idx, 0));
#endif
#if 1
        // get nearest neighbour in the right direction
        const int offset = (erectedNormal[dim] > 0) ? 1 : -1; 
        indicesBuffer = _velocityField[dim].cellIndex(cell_idx); 
        indicesBuffer[dim] += offset; 
        v(cell_idx, 0) = v(_velocityField[dim].cellIndex(indicesBuffer), 0); 
#endif
    }
}

void MAC_Grid::classifyCells( bool useBoundary )
{
    const int numPressureCells = _pressureField.numCells();
    const int numVcells[] = { _velocityField[ 0 ].numCells(), _velocityField[ 1 ].numCells(), _velocityField[ 2 ].numCells() };
    const int maxNumCells = std::max<int>(std::max<int>(std::max<int>(numPressureCells, numVcells[0]), numVcells[1]), numVcells[2]);
    const REAL &dt = _waveSolverSettings->timeStepSize; 
    const REAL &dx = _waveSolverSettings->cellSize; 
    const REAL &c = _waveSolverSettings->soundSpeed; 
    const REAL &rho = _waveSolverSettings->airDensity; 
    Vector3d cellPos;
    IntArray neighbours;

    neighbours.reserve( ScalarField::NUM_NEIGHBOURS );

    _isBulkCell.clear();
    _isGhostCell.clear();
    _classified.clear(); 
    _isVelocityInterfacialCell[ 0 ].clear();
    _isVelocityInterfacialCell[ 1 ].clear();
    _isVelocityInterfacialCell[ 2 ].clear();
    _isVelocityBulkCell[ 0 ].clear();
    _isVelocityBulkCell[ 1 ].clear();
    _isVelocityBulkCell[ 2 ].clear();

    _isBulkCell.resize(numPressureCells, false);
    _isGhostCell.resize(numPressureCells, false);
    _isPMLCell.resize(numPressureCells, false); 
    _classified.resize(maxNumCells, false); 

    _isVelocityInterfacialCell[ 0 ].resize( numVcells[ 0 ], false );
    _isVelocityInterfacialCell[ 1 ].resize( numVcells[ 1 ], false );
    _isVelocityInterfacialCell[ 2 ].resize( numVcells[ 2 ], false );
    _isVelocityBulkCell[ 0 ].resize( numVcells[ 0 ], false );
    _isVelocityBulkCell[ 1 ].resize( numVcells[ 1 ], false );
    _isVelocityBulkCell[ 2 ].resize( numVcells[ 2 ], false );

    _ghostCells.clear();

    _velocityInterfacialCells[ 0 ].clear();
    _velocityInterfacialCells[ 1 ].clear();
    _velocityInterfacialCells[ 2 ].clear();

    _interfacialBoundaryIDs[ 0 ].clear();
    _interfacialBoundaryIDs[ 1 ].clear();
    _interfacialBoundaryIDs[ 2 ].clear();
    _interfacialBoundaryDirections[ 0 ].clear();
    _interfacialBoundaryDirections[ 1 ].clear();
    _interfacialBoundaryDirections[ 2 ].clear();
    _interfacialBoundaryCoefficients[ 0 ].clear();
    _interfacialBoundaryCoefficients[ 1 ].clear();
    _interfacialBoundaryCoefficients[ 2 ].clear();

    _containingObject.clear(); 
    _containingObject.resize(numPressureCells, -1);

    _pmlPressureCells.clear(); 
    _pmlVelocityCells.clear(); 

    if ( !useBoundary )
    {
        std::cout << "No classification (not using boundary)\n"; 
        for ( int cell_idx = 0; cell_idx < numPressureCells; cell_idx++ )
        {
            _isBulkCell[ cell_idx ] = true;
        }

        for ( int dimension = 0; dimension < 3; dimension++ )
        {
            for ( int cell_idx = 0; cell_idx < _velocityField[ dimension ].numCells(); cell_idx++ )
            {
                _isVelocityBulkCell[dimension][cell_idx] = true; 
                _isVelocityInterfacialCell[dimension][cell_idx] = false; 
            }
        }
        return;
    }

    // Start by classifying all pressure cells in the computational
    // domain as bulk cells
    for ( int cell_idx = 0; cell_idx < numPressureCells; cell_idx++ )
    {
        cellPos = _pressureField.cellPosition( cell_idx );

        // Check all boundary fields to see if this is a bulk cell
        _isBulkCell[ cell_idx ] = true;
        for ( size_t field_idx = 0; field_idx < _boundaryFields.size(); field_idx++ )
        {
            if ( _boundaryFields[ field_idx ]->distance( cellPos )
                    <= _distanceTolerance )
            {
                _isBulkCell[ cell_idx ] = false;
                _containingObject[ cell_idx ] = field_idx;
                break;
            }
        }

        // classify and cache pml parameters
        const int pmlDirection = InsidePML(cellPos, _PML_absorptionWidth); 
        if (pmlDirection >= 0)
        {
            PML_PressureCell cell; 
            cell.index = cell_idx; 
            REAL maxAbsorptionCoefficient = std::numeric_limits<REAL>::min(); 
            for (int dim=0; dim<3; ++dim)
            {
                //const int &dim = pmlDirection; 
                //cell.dimension = dim;
                const REAL absorptionCoefficient = PML_absorptionCoefficient(cellPos, _PML_absorptionWidth, dim); 
                maxAbsorptionCoefficient = std::max<REAL>(maxAbsorptionCoefficient, absorptionCoefficient); 
                const REAL directionalCoefficient = PML_directionalCoefficient(absorptionCoefficient, dt);
                cell.updateCoefficient[dim] = PML_pressureUpdateCoefficient(absorptionCoefficient, dt, directionalCoefficient);
                cell.divergenceCoefficient[dim] = PML_divergenceCoefficient(rho, c, directionalCoefficient); 

                // figure out neighbours
                Tuple3i cellIndices = _pressureField.cellIndex(cell_idx); 
                cell.neighbour_v_left[dim] = _velocityField[dim].cellIndex(cellIndices); 
                cellIndices[dim] += 1;
                cell.neighbour_v_right[dim] = _velocityField[dim].cellIndex(cellIndices); 
            }
            
            cell.absorptionCoefficient = maxAbsorptionCoefficient; 
            _pmlPressureCells.push_back(cell); 
            _isPMLCell.at(cell_idx) = true; 
        }
    }

    // Ghost cells are cells inside the boundary that have at least
    // one neighbour outside of the boundary
    for ( int cell_idx = 0; cell_idx < numPressureCells; cell_idx++ )
    {
        cellPos = _pressureField.cellPosition( cell_idx );

        if ( _isBulkCell[ cell_idx ] )
        {
            // Not inside the interior object
            continue;
        }

        _pressureField.cellNeighbours( cell_idx, neighbours );

        for ( size_t neighbour_idx = 0; neighbour_idx < neighbours.size();
                neighbour_idx++ )
        {
            if ( _isBulkCell[ neighbours[ neighbour_idx ] ] )
            {
                // We have a neighbour outside of the interior object, so
                // this is a ghost cell
                _isGhostCell[ cell_idx ] = true;
                _ghostCells.push_back(cell_idx); 

                break;
            }
        }
    }

    // Classify interfacial cells in each velocity field
    for ( int dimension = 0; dimension < 3; dimension++ )
    {
        for ( int cell_idx = 0; cell_idx < _velocityField[ dimension ].numCells();
                cell_idx++ )
        {
            Tuple3i                cell_coordinates;
            int                    pressure_cell_idx1;
            int                    pressure_cell_idx2;
            const Vector3d cellPosition = _velocityField[dimension].cellPosition(cell_idx);

            cell_coordinates = _velocityField[ dimension ].cellIndex( cell_idx );

            // We will assume that boundary cells are not interfacial
            if ( cell_coordinates[ dimension ] == 0
                    || cell_coordinates[ dimension ]
                    == _velocityField[ dimension ].cellDivisions()[ dimension ] - 1 )
            {
                _isVelocityBulkCell[ dimension ][ cell_idx ] = true;
                continue;
            }

            // Look at our neighbours in the pressure field
            cell_coordinates[ dimension ] -= 1;
            pressure_cell_idx1 = _pressureField.cellIndex( cell_coordinates );
            cell_coordinates[ dimension ] += 1;
            pressure_cell_idx2 = _pressureField.cellIndex( cell_coordinates );

            if ( _isBulkCell[ pressure_cell_idx1 ]
                    && _isBulkCell[ pressure_cell_idx2 ] )
            {
                // Both pressure cell neighbours are bulk cells, so this is
                // a bulk cell too
                _isVelocityBulkCell[ dimension ][ cell_idx ] = true;
            }
            else if ( _isBulkCell[ pressure_cell_idx1 ] 
                    && !_isBulkCell[ pressure_cell_idx2 ] )
            {
                // Only one neighbour is inside the domain, so this must
                // be an interfacial cell
                _isVelocityInterfacialCell[ dimension ][ cell_idx ] = true;

                // Get the object ID for the boundary that we are adjacent to
                int boundaryObject = _containingObject[ pressure_cell_idx2 ];

                TRACE_ASSERT( boundaryObject >= 0 );

                _velocityInterfacialCells[ dimension ].push_back( cell_idx );
                _interfacialBoundaryIDs[ dimension ].push_back( boundaryObject );
                //_interfacialBoundaryDirections[ dimension ].push_back( -1.0 );
                _interfacialBoundaryDirections[ dimension ].push_back( 1.0 );

                // Determine a scaling coefficient based on the angle between
                // the boundary normal and the rasterized boundary normal
                Vector3d normal( 0.0, 0.0, 0.0 );
                Vector3d gradient = _boundaryFields[ boundaryObject ]->gradient(cellPosition);

                normal[ dimension ] = 1.0;
                gradient.normalize();

                REAL coefficient = normal.dotProduct( gradient );

                //TRACE_ASSERT( coefficient >= 0.0 );

                _interfacialBoundaryCoefficients[ dimension ].push_back( coefficient );

            }
            else if ( !_isBulkCell[ pressure_cell_idx1 ]
                    && _isBulkCell[ pressure_cell_idx2 ] )
            {
                // Only one neighbour is inside the domain, so this must
                // be an interfacial cell
                _isVelocityInterfacialCell[ dimension ][ cell_idx ] = true;

                // Get the object ID for the boundary that we are adjacent to
                int boundaryObject = _containingObject[ pressure_cell_idx1 ];

                TRACE_ASSERT( boundaryObject >= 0 );

                _velocityInterfacialCells[ dimension ].push_back( cell_idx );
                _interfacialBoundaryIDs[ dimension ].push_back( boundaryObject );
                _interfacialBoundaryDirections[ dimension ].push_back( 1.0 );

                // Determine a scaling coefficient based on the angle between
                // the boundary normal and the rasterized boundary normal
                Vector3d normal( 0.0, 0.0, 0.0 );
                Vector3d gradient = _boundaryFields[ boundaryObject ]->gradient(cellPosition);

                normal[ dimension ] = 1.0;
                gradient.normalize();

                REAL coefficient = normal.dotProduct( gradient );

                //TRACE_ASSERT( coefficient >= 0.0 );

                _interfacialBoundaryCoefficients[ dimension ].push_back( coefficient );
            }

            if (InsidePML(cellPosition, _PML_absorptionWidth) >= 0)
            {
                const REAL absorptionCoefficient = PML_absorptionCoefficient(cellPosition, _PML_absorptionWidth, dimension); 
                const REAL updateCoefficient = PML_velocityUpdateCoefficient(absorptionCoefficient, dt);
                const REAL gradientCoefficient  = PML_pressureGradientCoefficient(absorptionCoefficient, dt, dx, rho);
                PML_VelocityCell cell; 
                cell.index = cell_idx; 
                cell.dimension = dimension; 
                cell.updateCoefficient = updateCoefficient; 
                cell.gradientCoefficient = gradientCoefficient; 

                // figure out neighbours
                Tuple3i cellIndices = _velocityField[dimension].cellIndex(cell_idx); 
                cell.neighbour_p_right = _pressureField.cellIndex(cellIndices); 
                cellIndices[dimension] -= 1;
                cell.neighbour_p_left  = _pressureField.cellIndex(cellIndices); 
                if      ( _isPMLCell.at(cell.neighbour_p_right) &&  _isPMLCell.at(cell.neighbour_p_left)) cell.neighbourInterior =  0; 
                else if ( _isPMLCell.at(cell.neighbour_p_right) && !_isPMLCell.at(cell.neighbour_p_left)) cell.neighbourInterior = -1; 
                else if (!_isPMLCell.at(cell.neighbour_p_right) &&  _isPMLCell.at(cell.neighbour_p_left)) cell.neighbourInterior =  1; 
                else throw std::runtime_error("**ERROR** both pressure neighbours are not pml, this should not happen"); 

                _pmlVelocityCells.push_back(cell);
            }
        }
    }

    printf( "MAC_Grid: classifyCells:\n" );
    printf( "\tFound %d ghost cells\n", (int)_ghostCells.size() );
    printf( "\tFound %d v_x interfacial cells\n", (int)_velocityInterfacialCells[ 0 ].size() );
    printf( "\tFound %d v_y interfacial cells\n", (int)_velocityInterfacialCells[ 1 ].size() );
    printf( "\tFound %d v_z interfacial cells\n", (int)_velocityInterfacialCells[ 2 ].size() );
    printf( "\tFound %d pml pressure cells\n", (int)_pmlPressureCells.size() );
    printf( "\tFound %d pml velocity cells\n", (int)_pmlVelocityCells.size() );

    if (_useGhostCellBoundary)
    {
        ComputeGhostCellInverseMap(); 
    }
}

//##############################################################################
// Classify types of cells at each step. need to perform the following steps: 
// 
//  1. clear all the arrays (cell types, boundary ids)
//  2. if useBoundary is off, set everying to bulk and return; otherwise, continues 
//  3. classify all pressure cells by querying all distance fields, note that 
//     by definition, if a cell centroid is claimed by two different objects, 
//     this means interpenetration has happened or it is a rounding error. In
//     either cases, we simply choose the first object that claims it and
//     continues. 
//  4. Depending on whether ghostCells is on, go to 4a (on) or 4b (off)
//      4a.(1) clear arrays related to ghost cells. 
//         (2) scan through the pressure cells and examine neighbours. if one of
//             the neigbhours is bulk, mark it as ghost cell and continues 
//         (3) note: in this case, all velocity cells are bulk cells 
//      4b.(1) clear arrays related to the velocity cells 
//         (2) classify velocity cells as either bulk or interfacial, depending on
//             its neighbouring pressure cells 
//
//  @param pFull full pressure array
//  @param p component pressure array
//  @param pGCFull ghost cell full pressure array, this array will be resized depending on classification results.
//  @param pGC same as above, for ghost cell components
//  @param v velocities
//  @param useBoundary
//  @param verbose If classification results are printed.
//
//##############################################################################
void MAC_Grid::classifyCellsDynamic(MATRIX &pFull, MATRIX (&p)[3], FloatArray &pGCFull, FloatArray (&pGC)[3], MATRIX (&v)[3], const bool &useBoundary, const bool &verbose)
{
    const int numPressureCells = _pressureField.numCells(); 
    Vector3d cellPos;
    IntArray neighbours;
    neighbours.reserve( ScalarField::NUM_NEIGHBOURS );

    // step 2 : not using boundary, classification is trivial
    if (!useBoundary)
    {
        if (verbose) 
            std::cout << "No classification (not using boundary)\n"; 
        for (int dimension = 0; dimension < 3; ++dimension)
        {
            for (int cell_idx = 0; cell_idx < _velocityField[dimension].numCells(); ++cell_idx)
            {
                _isVelocityBulkCell[dimension][cell_idx] = true; 
                _isVelocityInterfacialCell[dimension][cell_idx] = false; 
            }
        }
        return;
    }

    // step 3 : classify pressure bulk cells 

    for ( int cell_idx = 0; cell_idx < numPressureCells; cell_idx++ )
    {
        // before doing anything, first update history array, ghost cell is
        // considered having valid history
        if (IsPressureCellSolid(cell_idx))
            _pressureCellHasValidHistory.at(cell_idx) = false; 
        else 
            _pressureCellHasValidHistory.at(cell_idx) = true; 

        cellPos = _pressureField.cellPosition( cell_idx );
        // Check all boundary fields to see if this is a bulk cell
        int containObjectId = _objects->OccupyByObject(cellPos); 
        
        if (containObjectId < 0)
        {
            // check additional samples within the cell. use six subdivide
            // stencils for now.
            for (int dim=0; dim<3; ++dim)
            {
                Vector3d offset(0, 0, 0);
                offset[dim] += _waveSolverSettings->cellSize*0.25; 
                // check positive side
                containObjectId = _objects->OccupyByObject(cellPos + offset); 
                if (containObjectId >= 0)
                    break; 
                // check negative side
                containObjectId = _objects->OccupyByObject(cellPos - offset); 
                if (containObjectId >= 0)
                    break; 
            }
        }
        _containingObject.at(cell_idx) = containObjectId; 
        const bool newIsBulkCell = (_containingObject.at(cell_idx)>=0 ? false : true); 

        _isBulkCell.at(cell_idx) = newIsBulkCell; 

        // since we know for sure non-bulk is solid if ghost cell is off, can clear it now
        if (!newIsBulkCell && !_useGhostCellBoundary) 
        {
            pFull(cell_idx, 0) = 0;
            p[0](cell_idx, 0) = 0; 
            p[1](cell_idx, 0) = 0; 
            p[2](cell_idx, 0) = 0; 
        }
    }

    // step 4a : classify pressure ghost cells
    //if (_useGhostCellBoundary)  // FIXME debug 
    {
        _ghostCells.clear(); 
        _ghostCellsChildren.clear(); 

        // examine ghost cells 
        for ( int cell_idx = 0; cell_idx < numPressureCells; ++cell_idx )
        {
            // if it is bulk, it is not ghost
            if (_isBulkCell.at(cell_idx)) 
            {
                _isGhostCell.at(cell_idx) = false; 
                continue; 
            }

            bool newIsGhostCell = false; 
            _pressureField.cellNeighbours(cell_idx, neighbours);
            for (size_t nei = 0; nei<neighbours.size(); ++nei)
            {
                const int neighbour_idx = neighbours.at(nei); 
                if (_isBulkCell.at(neighbour_idx))
                {
                    // We have a bulk neighbour so this is a ghost cell
                    newIsGhostCell = true; 
                    _ghostCells.push_back(cell_idx); 
                    IntArray children(6, -1);  // always fully subdivided
                    _ghostCellsChildren.push_back(children); 
                    break;
                }
            }
            if (!newIsGhostCell)
            {
                pFull(cell_idx, 0) = 0;
                p[0](cell_idx, 0) = 0; 
                p[1](cell_idx, 0) = 0; 
                p[2](cell_idx, 0) = 0; 
            }
            _isGhostCell.at(cell_idx) = newIsGhostCell; 
        }
        ComputeGhostCellInverseMap(); 
    }

    // step 4b : classify velocity cells for rasterized boundary
    pGCFull.clear(); 
    _ghostCellParents.clear(); 
    _ghostCellPositions.clear(); 
    _ghostCellBoundaryIDs.clear(); 
    for (int dim=0; dim<3; ++dim)
    {
        _velocityInterfacialCells[dim].clear(); 
        _interfacialBoundaryIDs[dim].clear(); 
        _interfacialBoundaryDirections[dim].clear(); 
        _interfacialBoundaryCoefficients[dim].clear(); 
        _interfacialGhostCellID[dim].clear(); 
        pGC[dim].clear(); 
    }
    int numSolidCells[3] = {0, 0, 0}; 
    for (int dimension = 0; dimension < 3; dimension++)
    {
        const int numVelocityCells = _velocityField[dimension].numCells(); 
        for (int cell_idx = 0; cell_idx < numVelocityCells; ++cell_idx)
        {
            // first establish if it has valid history
            if (IsVelocityCellSolid(cell_idx, dimension)) 
            //if (!_isVelocityBulkCell[dimension].at(cell_idx)) 
                _velocityCellHasValidHistory[dimension].at(cell_idx) = false; 
            else
                _velocityCellHasValidHistory[dimension].at(cell_idx) = true; 

            Tuple3i cell_coordinates = _velocityField[ dimension ].cellIndex( cell_idx );
            // We will assume that boundary cells are not interfacial
            if ( cell_coordinates[ dimension ] == 0
              || cell_coordinates[ dimension ] == _velocityField[ dimension ].cellDivisions()[ dimension ] - 1 )
            {
                _isVelocityBulkCell[dimension].at(cell_idx) = true; 
                _isVelocityInterfacialCell[dimension][cell_idx] = false; 
                continue;
            }

            // Look at our neighbours in the pressure field
            cell_coordinates[ dimension ] -= 1;
            int pressure_cell_idx1 = _pressureField.cellIndex( cell_coordinates );
            cell_coordinates[ dimension ] += 1;
            int pressure_cell_idx2 = _pressureField.cellIndex( cell_coordinates );

            if (_isBulkCell[ pressure_cell_idx1 ] && _isBulkCell[ pressure_cell_idx2 ])
            {
                // Both pressure cell neighbours are bulk cells, so this is
                // a bulk cell too
                _isVelocityBulkCell[dimension].at(cell_idx) = true; 
                _isVelocityInterfacialCell[dimension].at(cell_idx) = false; 
            }
            else if ( _isBulkCell[ pressure_cell_idx1 ] 
                    && !_isBulkCell[ pressure_cell_idx2 ] )
            {
                // Only one neighbour is inside the domain, so this must
                // be an interfacial cell
                _isVelocityBulkCell[dimension].at(cell_idx) = false; 
                _isVelocityInterfacialCell[ dimension ][ cell_idx ] = true;

                // Get the object ID for the boundary that we are adjacent to
                const int boundaryObject = _containingObject[ pressure_cell_idx2 ];

                TRACE_ASSERT( boundaryObject >= 0 );

                _velocityInterfacialCells[ dimension ].push_back( cell_idx );
                _interfacialBoundaryIDs[ dimension ].push_back( boundaryObject );
                _interfacialBoundaryDirections[ dimension ].push_back( -1.0 );
                //_interfacialBoundaryDirections[ dimension ].push_back( 1.0 );

                // Determine a scaling coefficient based on the angle between
                // the boundary normal and the rasterized boundary normal
                Vector3d normal( 0.0, 0.0, 0.0 );
                Vector3d position = _velocityField[ dimension ].cellPosition( cell_idx );
                Vector3d gradient; 
                _objects->ObjectNormal(boundaryObject, position, gradient); 
                gradient.normalize();
                //Vector3d gradient = _boundaryFields[ boundaryObject ]->gradient( position );

                normal[ dimension ] = 1.0;
                REAL coefficient = normal.dotProduct( gradient );

                //TRACE_ASSERT( coefficient >= 0.0 );

                _interfacialBoundaryCoefficients[ dimension ].push_back( coefficient );

                // initialize a ghost cell dedicated to this interfacial cell
                const int ghostCellIndex = pGCFull.size(); 
                _interfacialGhostCellID[dimension].push_back(ghostCellIndex); 
                Vector3d ghostCellPosition = position;
                ghostCellPosition[dimension] += _waveSolverSettings->cellSize*0.25;
                _ghostCellParents.push_back(pressure_cell_idx2);
                _ghostCellPositions.push_back(ghostCellPosition); 
                _ghostCellBoundaryIDs.push_back(boundaryObject); 

                const int childArrayPosition = dimension*2; // this is the childArrayPosition-th child in the tree
                _ghostCellsChildren.at(_ghostCellsInverse[pressure_cell_idx2]).at(childArrayPosition) = ghostCellIndex; 

                pGC[0].push_back(0.0); 
                pGC[1].push_back(0.0); 
                pGC[2].push_back(0.0); 
                pGCFull.push_back(0.0);
            }
            else if ( !_isBulkCell[ pressure_cell_idx1 ]
                    && _isBulkCell[ pressure_cell_idx2 ] )
            {
                // Only one neighbour is inside the domain, so this must
                // be an interfacial cell
                _isVelocityBulkCell[dimension].at(cell_idx) = false; 
                _isVelocityInterfacialCell[dimension].at(cell_idx) = true;

                // Get the object ID for the boundary that we are adjacent to
                int boundaryObject = _containingObject[ pressure_cell_idx1 ];

                TRACE_ASSERT( boundaryObject >= 0 );

                _velocityInterfacialCells[ dimension ].push_back( cell_idx );
                _interfacialBoundaryIDs[ dimension ].push_back( boundaryObject );
                _interfacialBoundaryDirections[ dimension ].push_back( 1.0 );

                // Determine a scaling coefficient based on the angle between
                // the boundary normal and the rasterized boundary normal
                Vector3d normal( 0.0, 0.0, 0.0 );
                Vector3d position = _velocityField[ dimension ].cellPosition( cell_idx );
                Vector3d gradient; 
                _objects->ObjectNormal(boundaryObject, position, gradient); 
                gradient.normalize();

                normal[ dimension ] = 1.0;
                REAL coefficient = normal.dotProduct( gradient );

                //TRACE_ASSERT( coefficient >= 0.0 );

                _interfacialBoundaryCoefficients[ dimension ].push_back( coefficient );

                // initialize a ghost cell dedicated to this interfacial cell
                const int ghostCellIndex = pGCFull.size(); 
                _interfacialGhostCellID[dimension].push_back(ghostCellIndex); 
                Vector3d ghostCellPosition = position;
                ghostCellPosition[dimension] -= _waveSolverSettings->cellSize*0.25;
                _ghostCellParents.push_back(pressure_cell_idx1);
                _ghostCellPositions.push_back(ghostCellPosition); 
                _ghostCellBoundaryIDs.push_back(boundaryObject); 

                const int childArrayPosition = dimension*2 + 1; // this is the childArrayPosition-th child in the tree
                _ghostCellsChildren.at(_ghostCellsInverse[pressure_cell_idx1]).at(childArrayPosition) = ghostCellIndex; 

                pGC[0].push_back(0.0); 
                pGC[1].push_back(0.0); 
                pGC[2].push_back(0.0); 
                pGCFull.push_back(0.0);
            }
            else // both sides aren't bulk, this is a solid cell
            {
                _isVelocityBulkCell[dimension][cell_idx] = false; 
                _isVelocityInterfacialCell[dimension][cell_idx] = false;
                v[dimension](cell_idx, 0) = 0.0; // clear solid velocity cell
                numSolidCells[dimension] ++; 
            }
        }
    }

    if (verbose) 
    {
        printf( "MAC_Grid: classifyCellsDynamic:\n" );
        printf( "\tFound %d ghost cells\n", (int)_ghostCells.size() );
        printf( "\tFound %d v_x interfacial cells\n", (int)_velocityInterfacialCells[ 0 ].size() );
        printf( "\tFound %d v_y interfacial cells\n", (int)_velocityInterfacialCells[ 1 ].size() );
        printf( "\tFound %d v_z interfacial cells\n", (int)_velocityInterfacialCells[ 2 ].size() );
        printf( "\tFound %d v_x solid cells\n", numSolidCells[0] );
        printf( "\tFound %d v_y solid cells\n", numSolidCells[1] );
        printf( "\tFound %d v_z solid cells\n", numSolidCells[2] );
    }
}

//##############################################################################
//##############################################################################
void MAC_Grid::classifyCellsDynamic_FAST(MATRIX &pFull, MATRIX (&p)[3], FloatArray &pGCFull, FloatArray (&pGC)[3], MATRIX (&v)[3], const bool &useBoundary, const bool &verbose)
{
    if (!useBoundary) return;

    // get all bounding boxes and iteration range
    const int N = _objects->N(); 
    std::vector<ScalarField::RangeIndices> indices(N); 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int object_id=0; object_id<N; ++object_id)
    {
        const FDTD_MovableObject::BoundingBox &unionBBox = _objects->Get(object_id).GetUnionBBox();
        const Vector3d maxBound = unionBBox.maxBound + 2.0*_cellSize; 
        const Vector3d minBound = unionBBox.minBound - 2.0*_cellSize; 
        _pressureField.GetIterationBox(minBound, maxBound, indices[object_id]); 
    } 

    SetClassifiedSubset(_pressureField, N, indices, false); 
    for (int bbox_id=0; bbox_id<N; ++bbox_id)
    {
        const Vector3i &start = indices[bbox_id].startIndex; 
        const Vector3i &range = indices[bbox_id].dimensionIteration; 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
        for (int ii=start.x; ii<start.x+range.x; ++ii)
        for (int jj=start.y; jj<start.y+range.y; ++jj)
        for (int kk=start.z; kk<start.z+range.z; ++kk)
        {
            const Tuple3i cellIndices(ii,jj,kk);
            const int cell_idx = _pressureField.cellIndex(cellIndices); 
            if (_classified.at(cell_idx)) continue; 
            const Vector3d cellPos = _pressureField.cellPosition( cell_idx );

            // before doing anything, first update history array
            if (!_isBulkCell.at(cell_idx))
                _pressureCellHasValidHistory.at(cell_idx) = false; 
            else 
                _pressureCellHasValidHistory.at(cell_idx) = true; 

            // Check all the point samples in the cell. If any one of them are 
            // occupied, break and set the cell as occupied
            int containObjectId = _objects->OccupyByObject(cellPos); 
            if (containObjectId < 0)
            {
                for (int dim=0; dim<3; ++dim)
                {
                    Vector3d offset(0, 0, 0);
                    offset[dim] = _waveSolverSettings->cellSize*0.25; 
                    containObjectId = _objects->OccupyByObject(cellPos + offset); 
                    if (containObjectId >= 0) break; 
                    containObjectId = _objects->OccupyByObject(cellPos - offset); 
                    if (containObjectId >= 0) break; 
                }
            }
            _containingObject.at(cell_idx) = containObjectId; 

            const bool newIsBulkCell = (containObjectId>=0 ? false : true); 
            _isBulkCell.at(cell_idx) = newIsBulkCell; 
            if (!newIsBulkCell) 
            {
                pFull(cell_idx, 0) = 0;
                p[0](cell_idx, 0) = 0; 
                p[1](cell_idx, 0) = 0; 
                p[2](cell_idx, 0) = 0; 
            }
            _classified.at(cell_idx) = true; // toggle on to prevent reclassification
        }
    }

    // Classify pressure ghost cells
    SetClassifiedSubset(_pressureField, N, indices, false); 
    _ghostCells.clear(); 
    _ghostCellsChildren.clear(); 
    IntArray childrenIndex(6, -1);  // always fully subdivided
#ifdef USE_OPENMP
    const int N_max_threads = omp_get_max_threads();
#else
    const int N_max_threads = 1;
#endif
    std::vector<IntArray> ghostCellsThreads(N_max_threads);  // used to store ghost cell index for each thread
    for (int bbox_id=0; bbox_id<N; ++bbox_id)
    {
        const Vector3i &start = indices[bbox_id].startIndex; 
        const Vector3i &range = indices[bbox_id].dimensionIteration; 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
        for (int ii=start.x; ii<start.x+range.x; ++ii)
        for (int jj=start.y; jj<start.y+range.y; ++jj)
        for (int kk=start.z; kk<start.z+range.z; ++kk)
        {
            const Tuple3i cellIndices(ii,jj,kk);
            const int cell_idx = _pressureField.cellIndex(cellIndices); 
#ifdef USE_OPENMP
            const int thread_idx = omp_get_thread_num(); 
#else
            const int thread_idx = 0; 
#endif
            if (_classified.at(cell_idx)) continue; 
            // if it is bulk, it is not ghost
            if (_isBulkCell.at(cell_idx)) 
            {
                _isGhostCell.at(cell_idx) = false; 
                continue; 
            }

            bool newIsGhostCell = false; 
            IntArray neighbours;
            _pressureField.cellNeighbours(cell_idx, neighbours);
            for (size_t nei = 0; nei<neighbours.size(); ++nei)
            {
                const int neighbour_idx = neighbours.at(nei); 
                if (_isBulkCell.at(neighbour_idx))
                {
                    // We have a bulk neighbour so this is a ghost cell
                    newIsGhostCell = true; 
                    //IntArray children(6, -1);  // always fully subdivided
                    //_ghostCells.push_back(cell_idx); 
                    //_ghostCellsChildren.push_back(childrenIndex); 
                    ghostCellsThreads.at(thread_idx).push_back(cell_idx); 
                    break;
                }
            }
            _isGhostCell.at(cell_idx) = newIsGhostCell; 
            _classified.at(cell_idx) = true; // toggle on to prevent reclassification
        }
    }
    // push back ghost cells
    for (const IntArray &gcIndThread : ghostCellsThreads)
        for (const int &gcInd : gcIndThread)
        {
            _ghostCells.push_back(gcInd); 
            _ghostCellsChildren.push_back(childrenIndex); 
        }
    ComputeGhostCellInverseMap(); 

    // Classify velocity cells
    pGCFull.clear(); 
    _ghostCellParents.clear(); 
    _ghostCellPositions.clear(); 
    _ghostCellBoundaryIDs.clear(); 
    for (int dim=0; dim<3; ++dim)
    {
        _velocityInterfacialCells[dim].clear(); 
        _interfacialBoundaryIDs[dim].clear(); 
        _interfacialBoundaryDirections[dim].clear(); 
        _interfacialBoundaryCoefficients[dim].clear(); 
        _interfacialGhostCellID[dim].clear(); 
        pGC[dim].clear(); 
    }
    int numSolidCells[3] = {0, 0, 0}; 
    SetClassifiedSubset(_pressureField, N, indices, false); 
    std::vector<std::vector<GhostCellInfo> > ghostCellInfoThreads(N_max_threads);
    //int countGhostCells = 0;
    for (int dimension = 0; dimension < 3; dimension++)
    {
        for (int bbox_id=0; bbox_id<N; ++bbox_id)
        {
            const Vector3i &start = indices[bbox_id].startIndex; 
            const Vector3i &range = indices[bbox_id].dimensionIteration; 
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
            for (int ii=start.x; ii<start.x+range.x; ++ii)
            for (int jj=start.y; jj<start.y+range.y; ++jj)
            for (int kk=start.z; kk<start.z+range.z; ++kk)
            {
#ifdef USE_OPENMP
                const int thread_idx = omp_get_thread_num(); 
#else
                const int thread_idx = 0; 
#endif
                const Tuple3i cellIndices(ii,jj,kk);
                const int cell_idx = _velocityField[dimension].cellIndex(cellIndices); 
                if (_classified.at(cell_idx)) continue; 
                if (IsVelocityCellSolid(cell_idx, dimension)) 
                    _velocityCellHasValidHistory[dimension].at(cell_idx) = false; 
                else
                    _velocityCellHasValidHistory[dimension].at(cell_idx) = true; 

                Tuple3i cell_coordinates = _velocityField[ dimension ].cellIndex( cell_idx );
                // We will assume that boundary cells are not interfacial
                if ( cell_coordinates[ dimension ] == 0
                  || cell_coordinates[ dimension ] == _velocityField[ dimension ].cellDivisions()[ dimension ] - 1 )
                {
                    _isVelocityBulkCell[dimension].at(cell_idx) = true; 
                    _isVelocityInterfacialCell[dimension][cell_idx] = false; 
                    continue;
                }

                // Look at our neighbours in the pressure field
                cell_coordinates[ dimension ] -= 1;
                const int pressure_cell_idx1 = _pressureField.cellIndex( cell_coordinates );
                cell_coordinates[ dimension ] += 1;
                const int pressure_cell_idx2 = _pressureField.cellIndex( cell_coordinates );

                GhostCellInfo ghostCellInfo;
                bool infoUsed = false; 
                if (_isBulkCell[ pressure_cell_idx1 ] && _isBulkCell[ pressure_cell_idx2 ])
                {
                    // Both pressure cell neighbours are bulk cells, so this is
                    // a bulk cell too
                    _isVelocityBulkCell[dimension].at(cell_idx) = true; 
                    _isVelocityInterfacialCell[dimension].at(cell_idx) = false; 
                }
                else if ( _isBulkCell[ pressure_cell_idx1 ] 
                        && !_isBulkCell[ pressure_cell_idx2 ] )
                {
                    // Only one neighbour is inside the domain, so this must
                    // be an interfacial cell
                    _isVelocityBulkCell[dimension].at(cell_idx) = false; 
                    _isVelocityInterfacialCell[ dimension ][ cell_idx ] = true;

                    // Get the object ID for the boundary that we are adjacent to
                    const int boundaryObject = _containingObject[ pressure_cell_idx2 ];
                    TRACE_ASSERT( boundaryObject >= 0 );

                    // Determine a scaling coefficient based on the angle between
                    // the boundary normal and the rasterized boundary normal
                    Vector3d normal( 0.0, 0.0, 0.0 );
                    Vector3d position = _velocityField[ dimension ].cellPosition( cell_idx );
                    Vector3d gradient; 
                    _objects->ObjectNormal(boundaryObject, position, gradient); 
                    gradient.normalize();
                    const REAL coefficient = gradient[dimension];
                    // initialize a ghost cell dedicated to this interfacial cell
                    //const int ghostCellIndex = pGCFull.size(); 
                    //const int ghostCellIndex = countGhostCells; 
                    Vector3d ghostCellPosition = position;
                    ghostCellPosition[dimension] += _waveSolverSettings->cellSize*0.25;
                    const int childArrayPosition = dimension*2; // this is the childArrayPosition-th child in the tree
                    //_ghostCellsChildren.at(_ghostCellsInverse[pressure_cell_idx2]).at(childArrayPosition) = ghostCellIndex; 
                    ghostCellInfo.dim = dimension; 
                    ghostCellInfo.cellIndex = cell_idx; 
                    ghostCellInfo.boundaryObject = boundaryObject; 
                    ghostCellInfo.boundaryDirection = -1.0; 
                    ghostCellInfo.boundaryCoefficient = coefficient; 
                    ghostCellInfo.childArrayPosition = childArrayPosition; 
                    ghostCellInfo.ghostCellParent = pressure_cell_idx2; 
                    ghostCellInfo.ghostCellPosition = ghostCellPosition; 
                    ghostCellInfo.ghostCellBoundaryID = boundaryObject; 
                    infoUsed = true;
                    //_velocityInterfacialCells[ dimension ].push_back( cell_idx );
                    //_interfacialBoundaryIDs[ dimension ].push_back( boundaryObject );
                    //_interfacialBoundaryDirections[ dimension ].push_back( -1.0 );
                    //_interfacialBoundaryCoefficients[ dimension ].push_back( coefficient );
                    //_interfacialGhostCellID[dimension].push_back(ghostCellIndex); 
                    //_ghostCellParents.push_back(pressure_cell_idx2);
                    //_ghostCellPositions.push_back(ghostCellPosition); 
                    //_ghostCellBoundaryIDs.push_back(boundaryObject); 
                    //pGC[0].push_back(0.0); 
                    //pGC[1].push_back(0.0); 
                    //pGC[2].push_back(0.0); 
                    //pGCFull.push_back(0.0);
                }
                else if ( !_isBulkCell[ pressure_cell_idx1 ]
                        && _isBulkCell[ pressure_cell_idx2 ] )
                {
                    // Only one neighbour is inside the domain, so this must
                    // be an interfacial cell
                    _isVelocityBulkCell[dimension].at(cell_idx) = false; 
                    _isVelocityInterfacialCell[dimension].at(cell_idx) = true;

                    // Get the object ID for the boundary that we are adjacent to
                    int boundaryObject = _containingObject[ pressure_cell_idx1 ];
                    TRACE_ASSERT( boundaryObject >= 0 );

                    // Determine a scaling coefficient based on the angle between
                    // the boundary normal and the rasterized boundary normal
                    Vector3d normal( 0.0, 0.0, 0.0 );
                    Vector3d position = _velocityField[ dimension ].cellPosition( cell_idx );
                    Vector3d gradient; 
                    _objects->ObjectNormal(boundaryObject, position, gradient); 
                    gradient.normalize();
                    const REAL coefficient = gradient[dimension];
                    // initialize a ghost cell dedicated to this interfacial cell
                    //const int ghostCellIndex = pGCFull.size(); 
                    //const int ghostCellIndex = countGhostCells; 
                    Vector3d ghostCellPosition = position;
                    ghostCellPosition[dimension] -= _waveSolverSettings->cellSize*0.25;
                    const int childArrayPosition = dimension*2 + 1; // this is the childArrayPosition-th child in the tree
                    //_ghostCellsChildren.at(_ghostCellsInverse[pressure_cell_idx1]).at(childArrayPosition) = ghostCellIndex; 

                    ghostCellInfo.dim = dimension; 
                    ghostCellInfo.cellIndex = cell_idx; 
                    ghostCellInfo.boundaryObject = boundaryObject; 
                    ghostCellInfo.boundaryDirection = 1.0; 
                    ghostCellInfo.boundaryCoefficient = coefficient; 
                    ghostCellInfo.childArrayPosition = childArrayPosition; 
                    ghostCellInfo.ghostCellParent = pressure_cell_idx1; 
                    ghostCellInfo.ghostCellPosition = ghostCellPosition; 
                    ghostCellInfo.ghostCellBoundaryID = boundaryObject; 
                    infoUsed = true;
                    //_velocityInterfacialCells[ dimension ].push_back( cell_idx );
                    //_interfacialBoundaryIDs[ dimension ].push_back( boundaryObject );
                    //_interfacialBoundaryDirections[ dimension ].push_back( 1.0 );
                    //_interfacialBoundaryCoefficients[ dimension ].push_back( coefficient );
                    //_interfacialGhostCellID[dimension].push_back(ghostCellIndex); 
                    //_ghostCellParents.push_back(pressure_cell_idx1);
                    //_ghostCellPositions.push_back(ghostCellPosition); 
                    //_ghostCellBoundaryIDs.push_back(boundaryObject); 
                    //pGC[0].push_back(0.0); 
                    //pGC[1].push_back(0.0); 
                    //pGC[2].push_back(0.0); 
                    //pGCFull.push_back(0.0);
                }
                else // both sides aren't bulk, this is a solid cell
                {
                    _isVelocityBulkCell[dimension][cell_idx] = false; 
                    _isVelocityInterfacialCell[dimension][cell_idx] = false;
                    v[dimension](cell_idx, 0) = 0.0; // clear solid velocity cell
                    numSolidCells[dimension] ++; 
                }
                if (infoUsed) 
                    ghostCellInfoThreads.at(thread_idx).push_back(ghostCellInfo); 
                _classified.at(cell_idx) = true; // toggle on to prevent reclassification
            }
        }
        SetClassifiedSubset(_velocityField[dimension], N, indices, false); 
    }
    // push back all info gathered by each threads
    int gcCount = 0;
    for (const auto &gcInfoThread : ghostCellInfoThreads)
        for (const GhostCellInfo &gcInfo : gcInfoThread)
        {
            Push_Back_GhostCellInfo(gcCount, gcInfo, pGCFull, pGC);
            gcCount++;
        }

    if (verbose) 
    {
        printf( "MAC_Grid: classifyCellsDynamic:\n" );
        printf( "\tFound %d ghost cells\n", (int)_ghostCells.size() );
        printf( "\tFound %d v_x interfacial cells\n", (int)_velocityInterfacialCells[ 0 ].size() );
        printf( "\tFound %d v_y interfacial cells\n", (int)_velocityInterfacialCells[ 1 ].size() );
        printf( "\tFound %d v_z interfacial cells\n", (int)_velocityInterfacialCells[ 2 ].size() );
        printf( "\tFound %d v_x solid cells\n", numSolidCells[0] );
        printf( "\tFound %d v_y solid cells\n", numSolidCells[1] );
        printf( "\tFound %d v_z solid cells\n", numSolidCells[2] );
    }
}

//##############################################################################
// This function handles the cell identity updates for Finite-Volume formulation
// after objects new position has been established in the timestep. 
//##############################################################################
void MAC_Grid::classifyCellsFV(MATRIX &pFull, MATRIX (&p)[3], FloatArray &pGCFull, FloatArray (&pGC)[3], MATRIX (&v)[3], const bool &useBoundary, const bool &verbose)
{
    // reset all fields
    std::fill(_isGhostCell.begin(), _isGhostCell.end(), false);
    _ghostCellsCollection.clear();

    // hash triangles
    _fvMetaData.Clear(); 
    const int N_objects = _objects->N();
    for (int obj_idx=0; obj_idx<N_objects; ++obj_idx)
    {
        RigidSoundObjectPtr object = _objects->GetPtr(obj_idx); 
        std::shared_ptr<TriangleMesh<REAL> > mesh = object->GetMeshPtr(); 
        std::shared_ptr<TriangleMeshKDTree<REAL> > meshkd = std::dynamic_pointer_cast<TriangleMeshKDTree<REAL> >(mesh); 
        const int N_triangles = mesh->num_triangles(); 
        for (int t_idx=0; t_idx<N_triangles; ++t_idx)
        {
            const Vector3d centroid = object->ObjectToWorldPoint(meshkd->TriangleCentroid(t_idx)); 
            const Vector3d normal   = object->ObjectToWorldVector(meshkd->triangle_normal(t_idx)); 
            const int cell_idx = InPressureCell(centroid); 
            const TriangleIdentifier tri_id(obj_idx, t_idx, centroid, normal);
            auto &list = _fvMetaData.cellMap[cell_idx]; 
            if (!list) 
                list = std::make_shared<std::vector<TriangleIdentifier> >(1, tri_id); 
            else
                list->push_back(tri_id);
            _isGhostCell.at(cell_idx) = true; 
        }
    }

    // classify all non-ghost
    const int N_pcells = _pressureField.numCells(); 
    for (int cell_idx=0; cell_idx<N_pcells; ++cell_idx)
    {
        _pressureCellHasValidHistory.at(cell_idx) = (IsPressureCellSolid(cell_idx) ? false : true); 
        if (_isGhostCell.at(cell_idx)) 
        {
            auto &hashedTriangles = _fvMetaData.cellMap.at(cell_idx); 
            std::shared_ptr<GhostCell> gc = std::make_shared<GhostCell>(cell_idx, hashedTriangles); 
            _ghostCellsCollection[cell_idx] = gc; 
            _isBulkCell.at(cell_idx) = false; 
            continue;
        }
        const Vector3d cellPosition = _pressureField.cellPosition(cell_idx);
        const REAL sdfQuery = _objects->LowestObjectDistance(cellPosition);
        _isBulkCell.at(cell_idx) = (sdfQuery < DISTANCE_TOLERANCE ? false : true);
        if (!_isBulkCell.at(cell_idx))
        {
            pFull(cell_idx, 0) = 0;
            p[0](cell_idx, 0) = 0; 
            p[1](cell_idx, 0) = 0; 
            p[2](cell_idx, 0) = 0; 
        }
    }
}

void MAC_Grid::ComputeGhostCellSolveResidual(const FloatArray &p, REAL &minResidual, REAL &maxResidual, int &maxResidualEntry, REAL &maxOffDiagonalEntry)
{
    const int N_ghostCells = _ghostCellCoupledData.size(); 
    if (N_ghostCells == 0)
        return; 

    Eigen::VectorXd residual(N_ghostCells); 
    maxOffDiagonalEntry = std::numeric_limits<REAL>::min(); 
    for (int r_idx=0; r_idx<N_ghostCells; ++r_idx)
    {
        const JacobiIterationData &data = _ghostCellCoupledData.at(r_idx); 
        residual(r_idx) = data.RHS;
        const int N_entries = data.nnzIndex.size(); 
        residual(r_idx) -= p.at(r_idx); 
        for (int e_idx=0; e_idx<N_entries; ++e_idx)
        {
            residual(r_idx) -= data.nnzValue[e_idx]*p.at(data.nnzIndex[e_idx]); 
            maxOffDiagonalEntry = max(maxOffDiagonalEntry, fabs(data.nnzValue[e_idx])); 
        }
    }
    minResidual = residual.minCoeff(); 
    maxResidual = residual.maxCoeff(&maxResidualEntry); 
}

REAL MAC_Grid::PressureCellType(const int &idx) const
{
    if (_isPMLCell.at(idx))
        return -0.5; 

    if (_isBulkCell.at(idx)) // bulk
        return 0.5;  
    else if (!_isBulkCell.at(idx) && !_isGhostCell.at(idx)) // solid
        return 1; 
    else  // ghost
        return -1;
}

void MAC_Grid::visualizeClassifiedCells()
{
    std::ofstream of("ghost_cell.csv"); 
    Vector3d position; 

    for (size_t gg=0; gg<_ghostCells.size(); gg++)
    {
        position = _pressureField.cellPosition(_pressureField.cellIndex(_ghostCells[gg])); 
        of << position.x << ", " << position.y << ", " << position.z << std::endl; 
    }

    std::ofstream of2("bulk_cell.csv"); 
    for (int bb=0; bb<_pressureField.numCells(); bb++) 
    {
        if (!_isBulkCell[bb])
            continue; 

        position = _pressureField.cellPosition(_pressureField.cellIndex(bb)); 
        of2 << position.x << ", " << position.y << ", " << position.z << std::endl; 
    }

    of.close(); 
    of2.close();
}

//##############################################################################
// Return PML dimension
//##############################################################################
int MAC_Grid::InsidePML(const Vector3d &x, const REAL &absorptionWidth)
{
    const int &preset = _waveSolverSettings->boundaryConditionPreset; 
    const BoundingBox &bbox = _pressureField.bbox();
    const REAL pmlWidth = absorptionWidth + _waveSolverSettings->cellSize * 0.10; // add buffer for detection

    for (int dimension=0; dimension<3; ++dimension)
    {
        const REAL hMin = x[ dimension ] - bbox.minBound()[ dimension ];
        const REAL hMax = bbox.maxBound()[ dimension ] - x[ dimension ];
        switch (preset)
        {
            case 0: // no wall
                if (hMin <= pmlWidth)
                    return dimension; 
                else if (hMax <= pmlWidth)
                    return dimension; 
                break; 
            case 1: // wall on +x, +y, +z
                if (hMin <= pmlWidth)
                    return dimension; 
                break; 
            case 2: 
                if (dimension == 2 && hMax <= pmlWidth)
                    return dimension; 
                break; 
            default: 
                break; 
        }
    }
    return -1;
}

REAL MAC_Grid::PML_absorptionCoefficient( const Vector3d &x, REAL absorptionWidth, int dimension )
{
#if 1
    const int &preset = _waveSolverSettings->boundaryConditionPreset; 
    const BoundingBox &bbox = _pressureField.bbox();
    const REAL hMin = x[ dimension ] - bbox.minBound()[ dimension ];
    const REAL hMax = bbox.maxBound()[ dimension ] - x[ dimension ];
    const REAL a2 = absorptionWidth * absorptionWidth;
    const REAL distMin = absorptionWidth - hMin; 
    const REAL distMax = absorptionWidth - hMax; 

    switch (preset)
    {
        case 0: // no wall
            if (hMin <= absorptionWidth)
                return _PML_absorptionStrength * pow(distMin,2) / a2; 
            else if (hMax <= absorptionWidth)
                return _PML_absorptionStrength * pow(distMax,2) / a2; 
            else 
                return 0.0;
            break; 
        case 1: // wall on +x, +y, +z
            if (hMin <= absorptionWidth)
                return _PML_absorptionStrength * pow(distMin,2) / a2; 
            else 
                return 0.0;
            break; 
        case 2: 
            if (dimension != 2 || hMin <= absorptionWidth) // wall on all but +z
                return 0.0; 
            else if (dimension == 2 && hMax <= absorptionWidth)
                return _PML_absorptionStrength * pow(distMax,2) / a2; 
            break; 
        default: 
            break; 
    }
    return 0.0;
#else // linear profile
    const int &preset = _waveSolverSettings->boundaryConditionPreset; 
    const BoundingBox &bbox = _pressureField.bbox();
    const REAL hMin = x[ dimension ] - bbox.minBound()[ dimension ];
    const REAL hMax = bbox.maxBound()[ dimension ] - x[ dimension ];
    const REAL distMin = absorptionWidth - hMin; 
    const REAL distMax = absorptionWidth - hMax; 

    switch (preset)
    {
        case 0: // no wall
            if (hMin <= absorptionWidth)
                return _PML_absorptionStrength * distMin / absorptionWidth; 
            else if (hMax <= absorptionWidth)
                return _PML_absorptionStrength * distMax / absorptionWidth; 
            else 
                return 0.0;
            break; 
        case 1: // wall on +x, +y, +z
            if (hMin <= absorptionWidth)
                return _PML_absorptionStrength * distMin / absorptionWidth; 
            else 
                return 0.0;
            break; 
        case 2: 
            if (dimension != 2 || hMin <= absorptionWidth) // wall on all but +z
                return 0.0; 
            else if (dimension == 2 && hMax <= absorptionWidth)
                return _PML_absorptionStrength * distMax / absorptionWidth; 
            break; 
        default: 
            break; 
    }
    return 0.0;
#endif
}

void MAC_Grid::FindImagePoint(const Vector3d &cellPosition, const int &boundaryObjectID, Vector3d &closestPoint, Vector3d &imagePoint, Vector3d &erectedNormal)
{
    const DistanceField *csdf = _boundaryFields[boundaryObjectID]; 

    erectedNormal = csdf->gradient(cellPosition); 
    erectedNormal.normalize(); 

    const REAL distance = csdf->distance(cellPosition); 
    closestPoint = cellPosition - erectedNormal * (1.0*distance); 
    imagePoint   = cellPosition - erectedNormal * (2.0*distance); 

    // produce warning if the following scenario occurs 
    //  1. erected normal deviates from the boundary normal by too much. 
    //  2. image point is still inside the geometry. 
    if (true) 
    {
        Vector3d boundaryNormal = csdf->gradient(closestPoint); 
        boundaryNormal.normalize(); 
        if (erectedNormal.dotProduct(boundaryNormal) < 0.5)
        {
            std::cerr << "**WARNING** erected normal and true normal deviates. This might cause inaccuracy for the imposed Neumann boundary condition at cell position : " 
                      << cellPosition 
                      << "; the dot product is : " << erectedNormal.dotProduct(boundaryNormal) << std::endl; 
        }

        if (csdf->distance(imagePoint) < -1E-3)
            std::cerr << "**ERROR** image point very inside object" << std::endl; 
    }


}

void MAC_Grid::ResetCellHistory(const bool &valid)
{
    auto &p = _pressureCellHasValidHistory; 
    auto &v = _velocityCellHasValidHistory; 
    std::fill(p.begin(), p.end(), valid); 
    std::fill(v[0].begin(), v[0].end(), valid); 
    std::fill(v[1].begin(), v[1].end(), valid); 
    std::fill(v[2].begin(), v[2].end(), valid); 
}

void MAC_Grid::GetCell(const int &cellIndex, MATRIX const (&pDirectional)[3], const MATRIX &pFull, const FloatArray &pGC, const MATRIX (&v)[3], Cell &cell) const
{
    cell.index = cellIndex; 
    cell.r_identity = PressureCellType(cellIndex); 
    if (EQUAL_FLOATS(cell.r_identity, -0.5))
        cell.s_identity = std::string("Pressure PML Cell"); 
    else if (EQUAL_FLOATS(cell.r_identity, 0.5))
        cell.s_identity = std::string("Pressure Bulk Cell"); 
    else if (EQUAL_FLOATS(cell.r_identity, 1.0))
        cell.s_identity = std::string("Pressure Solid Cell"); 
    else if (EQUAL_FLOATS(cell.r_identity, -1.0))
    {
        cell.s_identity = std::string("Pressure Ghost Cell"); 
#ifdef USE_FV
        try
        {
            auto gc = _ghostCellsCollection.at(cellIndex); 
            cell.ghostCell = gc; 
            cell.gcValue.resize(6); 
            for (int a_idx=0; a_idx<6; ++a_idx)
                cell.gcValue.at(a_idx) = gc->values.at(gc->valuePointer).at(a_idx); 
        }
        catch (std::out_of_range &e)
        {
            cell.ghostCell = nullptr; 
            cell.gcValue.resize(6, 0.0);
        }
#else
        cell.gcValue.resize(6, -1); 
        for (int a_idx=0; a_idx<6; ++a_idx)
        {
            const int &gc_idx = _ghostCellsChildren.at(_ghostCellsInverse.at(cellIndex)).at(a_idx); 
            if (gc_idx != -1)
                cell.gcValue.at(a_idx) = pGC.at(gc_idx);
        }
#endif
    }

    cell.indices = _pressureField.cellIndex(cellIndex); 
    cell.centroidPosition = _pressureField.cellPosition(cellIndex); 
    cell.lowerCorner = cell.centroidPosition - _waveSolverSettings->cellSize/2.0; 
    cell.upperCorner = cell.centroidPosition + _waveSolverSettings->cellSize/2.0; 

    cell.pFull = pFull(cellIndex, 0); 
    cell.pDirectional[0] = pDirectional[0](cellIndex, 0); 
    cell.pDirectional[1] = pDirectional[1](cellIndex, 0); 
    cell.pDirectional[2] = pDirectional[2](cellIndex, 0); 

    // lower velocity cell
    Tuple3i &indicesBuffer = cell.indices; 

    // lower velocity cell
    cell.vx[0] = v[0](_velocityField[0].cellIndex(Tuple3i(indicesBuffer[0], indicesBuffer[1], indicesBuffer[2])), 0); 
    cell.vy[0] = v[1](_velocityField[1].cellIndex(Tuple3i(indicesBuffer[0], indicesBuffer[1], indicesBuffer[2])), 0); 
    cell.vz[0] = v[2](_velocityField[2].cellIndex(Tuple3i(indicesBuffer[0], indicesBuffer[1], indicesBuffer[2])), 0); 

    const int xUpper = std::min<int>(indicesBuffer.x + 1, _waveSolverSettings->cellDivisions); 
    const int yUpper = std::min<int>(indicesBuffer.y + 1, _waveSolverSettings->cellDivisions); 
    const int zUpper = std::min<int>(indicesBuffer.z + 1, _waveSolverSettings->cellDivisions); 

    // upper velocity cell
    cell.vx[1] = v[0](_velocityField[0].cellIndex(Tuple3i(xUpper          , indicesBuffer[1], indicesBuffer[2])), 0); 
    cell.vy[1] = v[1](_velocityField[1].cellIndex(Tuple3i(indicesBuffer[0], yUpper          , indicesBuffer[2])), 0); 
    cell.vz[1] = v[2](_velocityField[2].cellIndex(Tuple3i(indicesBuffer[0], indicesBuffer[1], zUpper          )), 0); 

    cell.h = _waveSolverSettings->cellSize; 
}

void MAC_Grid::SetClassifiedSubset(const ScalarField &field, const int &N, const std::vector<ScalarField::RangeIndices> &indices, const bool &state)
{
    for (int bbox_id=0; bbox_id<N; ++bbox_id){ 
        const Vector3i &start = indices[bbox_id].startIndex;
        const Vector3i &range = indices[bbox_id].dimensionIteration;
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
        for (int ii=start.x; ii<start.x+range.x; ++ii){
        for (int jj=start.y; jj<start.y+range.y; ++jj){
        for (int kk=start.z; kk<start.z+range.z; ++kk){
            const Tuple3i cellIndices(ii,jj,kk);
            const int cell_idx = field.cellIndex(cellIndices);
            _classified.at(cell_idx) = state;}}}}
}

void MAC_Grid::CheckClassified()
{
    int countTrue = 0; 
    for (int ii=0; ii<_pressureField.numCells(); ++ii){
        if (_classified.at(ii)) 
            countTrue ++; } 
    std::cout << "there are " << countTrue << " trues in _classified\n"; 
}

// not thread-safe!
void MAC_Grid::Push_Back_GhostCellInfo(const int &gcIndex, const GhostCellInfo &info, FloatArray &pGCFull, FloatArray (&pGC)[3])
{
    _ghostCellsChildren.at(_ghostCellsInverse[info.ghostCellParent]).at(info.childArrayPosition) = gcIndex; 
    _velocityInterfacialCells[info.dim].push_back(info.cellIndex);
    _interfacialBoundaryIDs[info.dim].push_back(info.boundaryObject);
    _interfacialBoundaryDirections[info.dim].push_back(info.boundaryDirection);
    _interfacialBoundaryCoefficients[info.dim].push_back(info.boundaryCoefficient);
    _interfacialGhostCellID[info.dim].push_back(gcIndex); 
    _ghostCellParents.push_back(info.ghostCellParent);
    _ghostCellPositions.push_back(info.ghostCellPosition); 
    _ghostCellBoundaryIDs.push_back(info.ghostCellBoundaryID); 
    pGC[0].push_back(0.0); 
    pGC[1].push_back(0.0); 
    pGC[2].push_back(0.0); 
    pGCFull.push_back(0.0);
}

//##############################################################################
// This function returns flat index of the pressure cell that contains the 
// queried position. 
//##############################################################################
int MAC_Grid::InPressureCell(const Vector3d &position)
{
    const Vector3d &pBBoxLower = _pressureField.bbox().minBound();
    const REAL &h = _waveSolverSettings->cellSize; 
    const int &div = _waveSolverSettings->cellDivisions; 
    const Tuple3i cellIndices = Tuple3i(std::max<int>(std::min<int>((int)((position.x - pBBoxLower.x)/h), div-1), 0),
                                        std::max<int>(std::min<int>((int)((position.y - pBBoxLower.y)/h), div-1), 0),
                                        std::max<int>(std::min<int>((int)((position.z - pBBoxLower.z)/h), div-1), 0)); 
    return _pressureField.cellIndex(cellIndices); 
}

void MAC_Grid::FillVandermondeRegular(const int &row, const Vector3d &cellPosition, Eigen::MatrixXd &V)
{
    assert(row<8); 

    const double &x = cellPosition.x; 
    const double &y = cellPosition.y; 
    const double &z = cellPosition.z; 

    V(row,0) = x*y*z; 
    V(row,1) = x*y  ; 
    V(row,2) = x  *z; 
    V(row,3) =   y*z; 
    V(row,4) = x    ; 
    V(row,5) =   y  ; 
    V(row,6) =     z; 
    V(row,7) =     1; 
}

void MAC_Grid::FillVandermondeRegular(const Vector3d &cellPosition, Eigen::VectorXd &V)
{
    Eigen::MatrixXd tmp(1,8); 
    FillVandermondeRegular(0, cellPosition, tmp); 
    for (int ii=0; ii<8; ii++)
        V(ii) = tmp(0,ii);
}

void MAC_Grid::FillVandermondeBoundary(const int &row, const Vector3d &boundaryPosition, const Vector3d &boundaryNormal, Eigen::MatrixXd &V)
{
    assert(row<8); 

    const double &x = boundaryPosition.x; 
    const double &y = boundaryPosition.y; 
    const double &z = boundaryPosition.z; 

    const double &nx= boundaryNormal.x; 
    const double &ny= boundaryNormal.y; 
    const double &nz= boundaryNormal.z; 

    V(row,0) = nx*y*z + ny*x*z + nz*x*y; 
    V(row,1) = nx*y + ny*x; 
    V(row,2) = nx*z + nz*x; 
    V(row,3) = ny*z + nz*y; 
    V(row,4) = nx; 
    V(row,5) = ny; 
    V(row,6) = nz; 
    V(row,7) = 0; 
}


//// debug methods //// 
void MAC_Grid::PrintFieldExtremum(const MATRIX &field, const std::string &fieldName) 
{
    REAL minField = std::numeric_limits<REAL>::max(); 
    REAL maxField = std::numeric_limits<REAL>::min(); 

    for (int col=0; col<field.cols(); col++) 
        for (int row=0; row<field.rows(); row++) 
        {
            minField = std::min<REAL>(minField, field(row,col)); 
            maxField = std::max<REAL>(maxField, field(row,col)); 
        }

    std::cout   << " ----------------------------- \n"
                << " INFORMATION FOR FIELD : " << fieldName << "\n" 
                << " ----------------------------- \n"
                << "   min: " << minField << "\n"
                << "   max: " << maxField << "\n" ; 

    std::cout << std::flush; 
}

void MAC_Grid::PrintGhostCellTreeInfo() 
{
    const int N_gc = _ghostCellsChildren.size(); 
    for (int g_idx=0; g_idx<N_gc; ++g_idx)
    {
        std::cout << g_idx << ": "; 
        STL_Wrapper::PrintVectorContent(std::cout, _ghostCellsChildren.at(g_idx)); 
    }
    std::cout << std::flush;
}

bool MAC_Grid::ExamineJacobiMatrix()
{
    // check diagonally dominance. 
    const int N_rows = _ghostCellCoupledData.size(); 
    REAL maxL1 = std::numeric_limits<REAL>::min(); 
    IntArray problemRows; 
    for (int r_idx=0; r_idx<N_rows; ++r_idx)
    {
        const int N_nz = _ghostCellCoupledData.at(r_idx).nnzValue.size(); 
        REAL rowL1 = 0;
        for (int c_idx=0; c_idx<N_nz; ++c_idx)
            rowL1 += abs(_ghostCellCoupledData.at(r_idx).nnzValue.at(c_idx)); 
        if (rowL1 > 1.0) 
            problemRows.push_back(r_idx);
        maxL1 = std::max<REAL>(maxL1, rowL1); 
    }
    std::cout << "--------------------------------------------------------------------------------\n" 
              << "MAC_Grid::ExamineJacobiMatrix\n" 
              << "--------------------------------------------------------------------------------\n"
              << " Max off-diagonal L1 norm: " << maxL1 << "\n"
              << "...............................\n";
    for (const int &row : problemRows) 
    {
        std::cout << " Problem rows: " << row << "\n"; 
        const auto &data = _ghostCellCoupledData.at(row); 
        std::cout << "   ";
        for (const REAL &val : data.nnzValue)
            std::cout << val << " "; 
        std::cout << "\n";
    }
    std::cout << "--------------------------------------------------------------------------------" << std::endl;
    return true;
}

//##############################################################################
// Estimate the energy by using
//  E(t) = 0.5 \int_V (p_t^2 + |grad(p)|^2) dV
// One can show analytically that E'(t) = 0 and thus energy is conserving for 
// the wave equation. We will ignore all the scaling 0.5, dV etc.
//##############################################################################
REAL MAC_Grid::EstimateEnergy(const MATRIX &pCurr, const MATRIX &pLast)
{
    const int N_pcells = numPressureCells(); 
    const REAL &h = _waveSolverSettings->cellSize;
    const REAL &k = _waveSolverSettings->timeStepSize;
    const int &N_div = _waveSolverSettings->cellDivisions;
    REAL E = 0.0;
#ifdef USE_OPENMP
    const int N_max_threads = omp_get_max_threads();
#else
    const int N_max_threads = 1;
#endif
    std::vector<REAL> E_threads(N_max_threads, 0.0);
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int cell_idx=0; cell_idx<N_pcells; ++cell_idx)
    {
        Tuple3i cellIndices = _pressureField.cellIndex(cell_idx); 
        int cell_p, cell_n;
        bool compute_energy = true; 
        for (int dim=0; dim<3; ++dim)
        {
            const int indBuf = cellIndices[dim]; 
            cellIndices[dim] += 1; // p_i+1
            cellIndices[dim] = std::min<int>(N_div-1, cellIndices[dim]);
            cell_p = _pressureField.cellIndex(cellIndices); 
            cellIndices[dim] -= 2; // p_i-1
            cellIndices[dim] = std::max<int>(0, cellIndices[dim]);
            cell_n = _pressureField.cellIndex(cellIndices); 
            cellIndices[dim] = indBuf; 
            if (!_isBulkCell.at(cell_p) || !_isBulkCell.at(cell_n))
                compute_energy = false;
        }
#ifdef USE_OPENMP
            const int thread_idx = omp_get_thread_num(); 
#else
            const int thread_idx = 0; 
#endif
        if (compute_energy)
        {
            for (int dim=0; dim<3; ++dim)
                E_threads.at(thread_idx) += pow((pCurr(cell_p, 0) - pCurr(cell_n, 0))/h, 2); 
            E_threads.at(thread_idx) += (pCurr(cell_idx, 0) - pLast(cell_idx, 0))/k; 
        }
    }
    E = std::accumulate(E_threads.begin(), E_threads.end(), 0.0);
    return E;
}

//############################################################################## 
// This function stores the ghost cell coupled matrix used in jacobi iteration
// to a file. 
//
// format: 
//
// N_ghostcells(int)
// START
// NNZ_row(int) 
// entry_0_location(int) entry_0_value(REAL)
// entry_1_location(int) entry_1_value(REAL)
// entry_2_location(int) entry_2_value(REAL)
//                    .
//                    .
//                    .
// NNZ_row(int)
// entry_0_location(int) entry_0_value(REAL)
// entry_1_location(int) entry_1_value(REAL)
// entry_2_location(int) entry_2_value(REAL)
//                    .
//                    .
//                    .
// END
//############################################################################## 
void MAC_Grid::ToFile_GhostCellCoupledMatrix(const std::string &filename)
{
    if (!_waveSolverSettings->useGhostCell)
        return; 

#ifdef USE_OPENMP
#pragma omp critical
#endif
    {
        std::ofstream of(filename.c_str()); 
        const int N_ghostCells = _ghostCellCoupledData.size(); 
        of << N_ghostCells << "\nSTART\n";

        for (int g_idx=0; g_idx<N_ghostCells; ++g_idx)
        {
            const auto &data = _ghostCellCoupledData.at(g_idx); 
            const int N_nnz = data.nnzIndex.size() + 1;  // in data diagonal entry not stored
            assert((int)data.nnzValue.size()+1 == N_nnz); 

            of << N_nnz << "\n";
            of << g_idx << " " << "1.0" << std::endl; // diagonal entry is always 1
            for (int e_idx=0; e_idx<N_nnz-1; ++e_idx)
                of << data.nnzIndex.at(e_idx) << " " << data.nnzValue.at(e_idx) << "\n";
        }
        of << "END" << std::endl; 
    }
}

//############################################################################## 
//############################################################################## 
void MAC_Grid::ToFile_GhostCellLinearSystem(const char *filename)
{
    std::ofstream of(filename); 
    for (const auto &data : _ghostCellCoupledData)
        of << data.RHS << std::endl;
}

std::ostream &operator <<(std::ostream &os, const MAC_Grid &grid)
{
    const Vector3d pressure_minBound = grid.pressureField().minBound(); 
    const Vector3d pressure_maxBound = grid.pressureField().maxBound(); 
    const Vector3d velocity_x_minBound = grid.velocityField(0).minBound(); 
    const Vector3d velocity_x_maxBound = grid.velocityField(0).maxBound(); 
    const Vector3d velocity_y_minBound = grid.velocityField(1).minBound(); 
    const Vector3d velocity_y_maxBound = grid.velocityField(1).maxBound(); 
    const Vector3d velocity_z_minBound = grid.velocityField(2).minBound(); 
    const Vector3d velocity_z_maxBound = grid.velocityField(2).maxBound(); 
    const Tuple3i &pressure_divisions = grid.pressureField().cellDivisions(); 
    const Tuple3i &velocity_x_divisions = grid.velocityField(0).cellDivisions(); 
    const Tuple3i &velocity_y_divisions = grid.velocityField(1).cellDivisions(); 
    const Tuple3i &velocity_z_divisions = grid.velocityField(2).cellDivisions(); 
    
    os << "--------------------------------------------------------------------------------\n" 
       << "Class MAC_Grid\n" 
       << "--------------------------------------------------------------------------------\n"
       << " pressure field min     : " << pressure_minBound.x << ", " << pressure_minBound.y << ", " << pressure_minBound.z << "\n"
       << " pressure field max     : " << pressure_maxBound.x << ", " << pressure_maxBound.y << ", " << pressure_maxBound.z << "\n"
       << " velocity field min x   : " << velocity_x_minBound.x << ", " << velocity_x_minBound.y << ", " << velocity_x_minBound.z << "\n"
       << " velocity field max x   : " << velocity_x_maxBound.x << ", " << velocity_x_maxBound.y << ", " << velocity_x_maxBound.z << "\n"
       << " velocity field min y   : " << velocity_y_minBound.x << ", " << velocity_y_minBound.y << ", " << velocity_y_minBound.z << "\n"
       << " velocity field max y   : " << velocity_y_maxBound.x << ", " << velocity_y_maxBound.y << ", " << velocity_y_maxBound.z << "\n"
       << " velocity field min z   : " << velocity_z_minBound.x << ", " << velocity_z_minBound.y << ", " << velocity_z_minBound.z << "\n"
       << " velocity field max z   : " << velocity_z_maxBound.x << ", " << velocity_z_maxBound.y << ", " << velocity_z_maxBound.z << "\n"
       << " number pressure cells  : " << pressure_divisions.x << "x" << pressure_divisions.y << "x" << pressure_divisions.z << " = " << grid.numPressureCells() << "\n"
       << " number velocity cells x: " << velocity_x_divisions.x << "x" << velocity_x_divisions.y << "x" << velocity_x_divisions.z << " = " << grid.numVelocityCellsX() << "\n"
       << " number velocity cells y: " << velocity_y_divisions.x << "x" << velocity_y_divisions.y << "x" << velocity_y_divisions.z << " = " << grid.numVelocityCellsY() << "\n"
       << " number velocity cells z: " << velocity_z_divisions.x << "x" << velocity_z_divisions.y << "x" << velocity_z_divisions.z << " = " << grid.numVelocityCellsZ() << "\n"
       << "--------------------------------------------------------------------------------" 
       << std::flush; 
    return os; 
}

std::ostream &operator <<(std::ostream &os, const MAC_Grid::Cell &cell)
{
    os << "--------------------------------------------------------------------------------\n" 
       << "Class MAC_Grid::Cell\n" 
       << "--------------------------------------------------------------------------------\n"
       << " index            : " << cell.index << "\n"
       << " indices          : " << cell.indices.x << ", " << cell.indices.y << ", " << cell.indices.z << "\n"
       << " identity         : " << cell.s_identity << "\n"
       << " centroid position: " << cell.centroidPosition.x << ", " << cell.centroidPosition.y << ", " << cell.centroidPosition.z << "\n"
       << " ......................" << "\n"
       << " pDirectional     : " << cell.pDirectional[0]<< ", " << cell.pDirectional[1] << ", " << cell.pDirectional[2] << "\n"
       << " pFull            : " << cell.pFull << "\n"
       << " vx               : " << cell.vx[0] << ", " << cell.vx[1] << "\n"
       << " vy               : " << cell.vy[0] << ", " << cell.vy[1] << "\n"
       << " vz               : " << cell.vz[0] << ", " << cell.vz[1] << "\n"
       << " divergence       : " << cell.Divergence() << "\n"
       << " laplacian        : " << cell.laplacian << "\n"
       << " ......................" << "\n";
    if (cell.s_identity.compare("Pressure Ghost Cell")==0)
    {
        os << " pGhostCell       : ["; 
        for (int ii=0; ii<5; ++ii) 
            os << cell.gcValue.at(ii) << ", "; 
        os << cell.gcValue.at(5) << "]\n"; 
#ifdef USE_FV
        if (cell.ghostCell)
        {
            os << " volume           : " << cell.ghostCell->volume      << "\n"; 
            os << " dp_dn_dot_S      : " << cell.ghostCell->dp_dn_dot_S << "\n"; 
        }
#endif
    }
    os << "--------------------------------------------------------------------------------" 
       << std::flush; 
    return os; 
}
