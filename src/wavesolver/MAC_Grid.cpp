//////////////////////////////////////////////////////////////////////
// MAC_Grid.cpp: Implementation of the MAC_Grid class
//
//////////////////////////////////////////////////////////////////////

#include "MAC_Grid.h"

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

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
#include <wavesolver/BoundaryInterface.h>

#include <queue>
#include <unistd.h> 

#ifdef USE_OPENMP
#include <omp.h>
#endif

//#define ENGQUIST_ORDER 1
#define ENGQUIST_ORDER 2

//##############################################################################
// Static variable initialize
//##############################################################################
//int MAC_Grid::GhostCell::valuePointer = 0;
//std::vector<Timer<false> > MAC_Grid::GhostCell::ghostCellTimers(20);

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

MAC_Grid::MAC_Grid(const BoundingBox &bbox, 
                   PML_WaveSolver_Settings_Ptr settings, 
                   std::shared_ptr<FDTD_Objects> objects)
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
    const Tuple3i &pFieldDivs = _pressureField.cellDivisions(); 
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
        if (cellIndices[0]==0 || cellIndices[0]==pFieldDivs[0]-1 ||
            cellIndices[1]==0 || cellIndices[1]==pFieldDivs[1]-1 ||
            cellIndices[2]==0 || cellIndices[2]==pFieldDivs[2]-1)
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
    const Tuple3i &pFieldDivs = _pressureField.cellDivisions(); 
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
        for (int dim=0; dim<3; ++dim)
        {
            int isBoundaryFace = 0; // detect boundary face
            if (cellIndices[dim]==0) 
                isBoundaryFace = -1;
            else if (cellIndices[dim]==pFieldDivs[dim]-1)
                isBoundaryFace =  1;

            if (isBoundaryFace!=+1) // v_i+1 if exists, otherwise v_i
                bufPos[dim] += 1; 
            if (isBoundaryFace!=-1) // v_i-1 if exists, otherwise v_i
                bufNeg[dim] -= 1; 
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
                laplacian(cell_idx, 0) += (  _ghostCells.at(GhostCell::MakeKey(cell_idx,buf_iPos))->pressure
                                           + _ghostCells.at(GhostCell::MakeKey(cell_idx,buf_iNeg))->pressure
                                           - 2.0*value(cell_idx, 0)) / (9.0/16.0); 
            else if (_isGhostCell.at(buf_iPos) && !_isGhostCell.at(buf_iNeg)) // only right side is ghost cell
                laplacian(cell_idx, 0) += (  _ghostCells.at(GhostCell::MakeKey(cell_idx,buf_iPos))->pressure
                                           + 0.75*value(buf_iNeg, 0)
                                           - (7.0/4.0)*value(cell_idx, 0)) / (21.0/32.0); 
            else if (!_isGhostCell.at(buf_iPos) && _isGhostCell.at(buf_iNeg)) // only left side is ghost cell
                laplacian(cell_idx, 0) += (  0.75*value(buf_iPos, 0)
                                           + _ghostCells.at(GhostCell::MakeKey(cell_idx,buf_iNeg))->pressure
                                           - (7.0/4.0)*value(cell_idx, 0)) / (21.0/32.0); 
            else // both side is bulk
                laplacian(cell_idx, 0) += (  value(buf_iPos, 0)
                                           + value(buf_iNeg, 0)
                                           - 2.0*value(cell_idx, 0)); 
#endif
            if (isBoundaryFace!=+1) // v_i+1 if exists, otherwise v_i
                bufPos[dim] -= 1; 
            if (isBoundaryFace!=-1) // v_i-1 if exists, otherwise v_i
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
        auto &cell = _pmlVelocityCells.at(ii);
        int &cell_idx = cell.index; 
        // update velocity
        v[cell.dimension](cell_idx, 0) *= cell.updateCoefficient; 
        v[cell.dimension](cell_idx, 0) += cell.gradientCoefficient * pFull(cell.neighbour_p_right, 0); 
        v[cell.dimension](cell_idx, 0) -= cell.gradientCoefficient * pFull(cell.neighbour_p_left , 0); 
        //cell.velocity *= cell.updateCoefficient; 
        //cell.velocity += cell.gradientCoefficient * pFull(cell.neighbour_p_right, 0); 
        //cell.velocity -= cell.gradientCoefficient * pFull(cell.neighbour_p_left , 0); 
    }
}

void MAC_Grid::PML_pressureUpdateCollocated(const REAL &simulationTime, const MATRIX (&v)[3], MATRIX (&pDirectional)[3], MATRIX &pLast, MATRIX &pCurr, MATRIX &pNext, MATRIX &currLaplacian, MATRIX &lastLaplacian)
{
    const REAL &timeStep = _waveSolverSettings->timeStepSize; 
    const REAL one_over_dx = 1.0/_waveSolverSettings->cellSize; 
    const REAL c2_k2 = pow(_waveSolverSettings->soundSpeed*timeStep, 2); 
    const REAL c_alp = _waveSolverSettings->soundSpeed*_waveSolverSettings->alpha; 
    const REAL kcalp = c_alp*timeStep; 
    const REAL lambda = _waveSolverSettings->soundSpeed*_waveSolverSettings->timeStepSize/_waveSolverSettings->cellSize; 
    const REAL lambda2 = pow(lambda,2); 
    const REAL lambda3 = pow(lambda,3); 
    const int N_cells = _pressureField.numCells(); 
    const int N_pmlCells = _pmlPressureCells.size(); 
    const bool evaluateExternalSource = _objects->HasExternalPressureSources();
    const Tuple3i &Ns = _pressureField.cellDivisions(); 
      
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
        const Tuple3i indices = _pressureField.cellIndex(cell_idx); 
#if   ENGQUIST_ORDER == 1
        const Tuple3i nml     = _pressureField.boundaryCellNormal(cell_idx); 
#elif ENGQUIST_ORDER == 2
        const Tuple3i nml     = _pressureField.boundaryCellNormal(cell_idx, 1); 
        const Tuple3i nml0    = _pressureField.boundaryCellNormal(cell_idx); 
        if (nml0.abs().sum()>0) continue; 
#endif

#define PCELL_IDX(di_,dj_,dk_,ii_,jj_,kk_) _pressureField.cellIndex(di_,dj_,dk_,ii_,jj_,kk_)
        const int btype = nml.abs().sum(); 
        const bool is_boundary = (btype > 0); 

        if (is_boundary)
        {
            // redefine dimension so that the normal is at x-direction (ii)
            // n: negative to boundary normal
            // p: positive to boundary normal
            // i.e., ii_p is the exterior
            int di=0, dj=0, dk=0; 
            int ii  , jj  , kk  ; 
            int ii_n, jj_n, kk_n; 
            int ii_p, jj_p, kk_p; 
            REAL p_ii_p, p_ii_n, p_jj_p, p_jj_n, p_kk_p, p_kk_n, p1; 
            if (btype == 1) // face
            {
                for (di=0; di<3; ++di)
                    if (nml[di]!=0) break; 
                dj = (di+1)%3; 
                dk = (di+2)%3; 
                ii = indices[di]; jj = indices[dj]; kk = indices[dk];
                ii_n = ii - nml[di]; 
                ii_p = ii + nml[di]; 
                jj_n = std::max(jj-1,0); jj_p = std::min(jj+1, Ns[dj]-1); 
                kk_n = std::max(kk-1,0); kk_p = std::min(kk+1, Ns[dk]-1); 
                // pressure of the neighbours
                int neighbour_idx; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii_n,jj,kk); 
                p_ii_n = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii,jj_p,kk); 
                p_jj_p = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii,jj_n,kk); 
                p_jj_n = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii,jj,kk_p); 
                p_kk_p = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii,jj,kk_n); 
                p_kk_n = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 

                // first, compute a regular pNext using ABC
#if   ENGQUIST_ORDER == 1
                p1 = (2./lambda2 - 6.)     *pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0)
                   + (lambda - 1.0)/lambda2*pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0)
                   + 2.*p_ii_n + p_jj_p + p_jj_n + p_kk_p + p_kk_n;
                p1 *= (lambda2/(1. + lambda)); 
                p_ii_p = (pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0) - p1)/lambda + p_ii_n;
#elif ENGQUIST_ORDER == 2
                p1 = (2./lambda2 + 4./lambda - 6. - 4.*lambda) * pCurr(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk  ),0)
                   - (1./lambda2 + 2./lambda                 ) * pLast(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk  ),0)
                   +                                             pLast(PCELL_IDX(di,dj,dk,ii_p,jj  ,kk  ),0)
                   -                                             pLast(PCELL_IDX(di,dj,dk,ii_n,jj  ,kk  ),0)
                   + 2.                                        * pCurr(PCELL_IDX(di,dj,dk,ii_n,jj  ,kk  ),0)
                   + (lambda + 1.                            ) *(pCurr(PCELL_IDX(di,dj,dk,ii  ,jj_p,kk  ),0)
                                                               + pCurr(PCELL_IDX(di,dj,dk,ii  ,jj_n,kk  ),0)
                                                               + pCurr(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk_p),0)
                                                               + pCurr(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk_n),0)); 
                p1 *= (lambda2/(1. + 2.*lambda));
                p_ii_p = pLast(PCELL_IDX(di,dj,dk,ii_p,jj,kk),0)
                       + pCurr(PCELL_IDX(di,dj,dk,ii_n,jj,kk),0) 
                       - pLast(PCELL_IDX(di,dj,dk,ii_n,jj,kk),0) 
                       - 2./lambda*(p1 + pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0) - 2.*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0))
                       + lambda *  (pCurr(PCELL_IDX(di,dj,dk,ii  ,jj_p,kk  ),0)
                                +   pCurr(PCELL_IDX(di,dj,dk,ii  ,jj_n,kk  ),0)
                                +   pCurr(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk_p),0)
                                +   pCurr(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk_n),0)
                                -4.*pCurr(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk  ),0)); 
#endif
                pCurr(PCELL_IDX(di,dj,dk,ii_p,jj,kk),0) = p_ii_p;

                // fetch for possible contribution through interface (with other sim units)
                REAL p1_src = 0.0;
                REAL dbuf; 
                REAL alpha=-1.0; 
                for (auto &interface : _boundaryInterfaces)
                {
                    if (interface->GetDirection()!=di)
                        continue; 
                    const bool hasNeighbour = interface->GetOtherCellPressure(*grid_id, cell_idx, dbuf); 
                    if (hasNeighbour)
                    {
                        alpha = interface->GetBlendCoeff(simulationTime); 
                        p1_src += dbuf; 
                    }
                }
                // alpha-blend the neighbour pressure with ABC to smooth out discont.
                // then do a simple six-point Laplacian computation
                REAL p_blend; 
                if      (alpha<0)                  p_blend = p_ii_p; 
                else if (alpha>=0.0 && alpha<=1.0) p_blend = (1.0 - alpha)*p_ii_p + alpha*p1_src; 
                else                               p_blend = p1_src; 
                pNext(cell_idx,0) = 
                    + 2.0*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0) 
                    -     pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0)
                    + lambda2*(  
                        + p_ii_n + p_jj_p + p_jj_n + p_kk_p + p_kk_n
                        + p_blend
                        - 6.*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0));
            }
            else if (btype == 2) // edge
            {
                for (di=0; di<3; ++di) if (nml[di]!=0) break;
                for (dj=di+1; dj<3; ++dj) if (nml[dj]!=0) break; 
                if (dj!=di+1) dk = (di+1)%3; 
                else          dk = (di+2)%3; 
                ii = indices[di]; jj = indices[dj]; kk = indices[dk];
                ii_n = ii - nml[di]; 
                jj_n = jj - nml[dj]; 
                ii_p = ii + nml[di]; 
                jj_p = jj + nml[dj]; 
                kk_n = std::max(kk-1,0); kk_p = std::min(kk+1, Ns[dk]-1); 
                // pressure of the neighbours
                int neighbour_idx; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii_n,jj,kk); 
                p_ii_n = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii,jj_n,kk); 
                p_jj_n = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii,jj,kk_p); 
                p_kk_p = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii,jj,kk_n); 
                p_kk_n = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 

                // first, compute a regular pNext using ABC
//#if   ENGQUIST_ORDER == 1
                p1 = (2./lambda2 - 6. )        *pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0)
                   + (2.0*lambda - 1.0)/lambda2*pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0)
                   + 2.*p_ii_n + 2.*p_jj_n + p_kk_p + p_kk_n;
                p1 *= (lambda2/(1. + 2.*lambda)); 
                p_ii_p = (pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0) - p1)/lambda + p_ii_n;
                p_jj_p = (pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0) - p1)/lambda + p_jj_n;
//#elif ENGQUIST_ORDER == 2
//                p1 = (2. + 6.*lambda - 6.*lambda2 - 2.*lambda3)* pCurr(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk  ),0)
//                   + lambda2                                   *(pLast(PCELL_IDX(di,dj,dk,ii_p,jj  ,kk  ),0)
//                                                                +pLast(PCELL_IDX(di,dj,dk,ii  ,jj_p,kk  ),0))
//                   + (lambda3 + lambda2 - lambda + 1.         )*(pCurr(PCELL_IDX(di,dj,dk,ii_n,jj  ,kk  ),0)
//                                                                +pCurr(PCELL_IDX(di,dj,dk,ii  ,jj_n,kk  ),0))
//                   + ((1.+lambda)*lambda2                     )*(pCurr(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk_p),0)
//                                                                +pCurr(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk_n),0));
//                p1 *= (1./(1.+3.*lambda)); 
//                p1 += -pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0); 
//                if (EQUAL_FLOATS(lambda,1.))
//                    throw std::runtime_error("**WARNING** lambda is too close to 1, this formulation is singular");
//                REAL p_sp = -(p1 + pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0) - 2.*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0))
//                          - (lambda /4.)*(-   (pLast(PCELL_IDX(di,dj,dk,ii_p,jj,kk  ),0) + pLast(PCELL_IDX(di,dj,dk,ii,jj_p,kk  ),0))
//                                          -   (pCurr(PCELL_IDX(di,dj,dk,ii_n,jj,kk  ),0) + pCurr(PCELL_IDX(di,dj,dk,ii,jj_n,kk  ),0))
//                                          +   (pLast(PCELL_IDX(di,dj,dk,ii_n,jj,kk  ),0) + pLast(PCELL_IDX(di,dj,dk,ii,jj_n,kk  ),0)))
//                          + (lambda2/4.)*(    (pCurr(PCELL_IDX(di,dj,dk,ii_n,jj,kk  ),0) + pCurr(PCELL_IDX(di,dj,dk,ii,jj_n,kk  ),0))
//                                          +2.*(pCurr(PCELL_IDX(di,dj,dk,ii  ,jj,kk_p),0) + pLast(PCELL_IDX(di,dj,dk,ii,jj  ,kk_n),0))
//                                          -8.* pCurr(PCELL_IDX(di,dj,dk,ii  ,jj,kk  ),0));
//                p_sp *= (4./lambda/(1.-lambda));
//                p_ii_p = p_sp / 2.0; 
//                p_jj_p = p_sp / 2.0;
//#endif
                pCurr(PCELL_IDX(di,dj,dk,ii_p,jj  ,kk),0) = p_ii_p;
                pCurr(PCELL_IDX(di,dj,dk,ii  ,jj_p,kk),0) = p_jj_p;
                pCurr(PCELL_IDX(di,dj,dk,ii_p,jj_p,kk),0) = 0.0;

                // fetch for possible contribution through interface (with other sim units)
                REAL p1_src = 0.0, p2_src = 0.0; 
                REAL dbuf; 
                REAL alpha1=-1.0, alpha2=-1.0; 
                for (auto &interface : _boundaryInterfaces)
                {
                    if (interface->GetDirection()==di)
                    {
                        const bool hasNeighbour = interface->GetOtherCellPressure(*grid_id, cell_idx, dbuf); 
                        if (hasNeighbour)
                        {
                            alpha1 = interface->GetBlendCoeff(simulationTime); 
                            p1_src += dbuf; 
                        }
                    }
                    else if (interface->GetDirection()==dj)
                    {
                        const bool hasNeighbour = interface->GetOtherCellPressure(*grid_id, cell_idx, dbuf); 
                        if (hasNeighbour)
                        {
                            alpha2 = interface->GetBlendCoeff(simulationTime); 
                            p2_src += dbuf; 
                        }
                    }
                }
                // alpha-blend the neighbour pressure with ABC to smooth out discont.
                // then do a simple six-point Laplacian computation
                REAL p1_blend, p2_blend; 
                if      (alpha1<0)                   p1_blend = p_ii_p; 
                else if (alpha1>=0.0 && alpha1<=1.0) p1_blend = (1.0 - alpha1)*p_ii_p + alpha1*p1_src; 
                else                                 p1_blend = p1_src; 
                if      (alpha2<0)                   p2_blend = p_jj_p; 
                else if (alpha2>=0.0 && alpha2<=1.0) p2_blend = (1.0 - alpha2)*p_jj_p + alpha2*p2_src; 
                else                                 p2_blend = p2_src; 
                pNext(cell_idx,0) = 
                    + 2.0*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0) 
                    -     pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0)
                    + lambda2*(  
                        + p_ii_n  + p_jj_n + p_kk_p + p_kk_n
                        + p1_blend + p2_blend
                        - 6.*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0));
            }
            else if (btype == 3) // corner
            {
                di=0; dj=1; dk=2;
                ii = indices[di]; jj = indices[dj]; kk = indices[dk];
                ii_n = ii - nml[di]; 
                jj_n = jj - nml[dj]; 
                kk_n = kk - nml[dk]; 
                ii_p = ii + nml[di]; 
                jj_p = jj + nml[dj]; 
                kk_p = kk + nml[dk]; 
                // pressure of the neighbours
                int neighbour_idx; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii_n,jj,kk); 
                p_ii_n = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii,jj_n,kk); 
                p_jj_n = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 
                neighbour_idx = PCELL_IDX(di,dj,dk,ii,jj,kk_n); 
                p_kk_n = (!_isGhostCell.at(neighbour_idx)) ? pCurr(neighbour_idx,0)
                                                           : _ghostCells.at(GhostCell::MakeKey(cell_idx,neighbour_idx))->pressure; 

                // first, compute a regular pNext using ABC
//#if   ENGQUIST_ORDER == 1
                p1 = (2./lambda2 - 6. )        *pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0)
                   + (3.0*lambda - 1.0)/lambda2*pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0)
                   + 2.*p_ii_n + 2.*p_jj_n + 2.*p_kk_n;
                p1 *= (lambda2/(1. + 3.*lambda)); 
                p_ii_p = (pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0) - p1)/lambda + p_ii_n;
                p_jj_p = (pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0) - p1)/lambda + p_jj_n;
                p_kk_p = (pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0) - p1)/lambda + p_kk_n;
//#elif ENGQUIST_ORDER == 2
//                const REAL p_spnm = pLast(PCELL_IDX(di,dj,dk,ii_p,jj  ,kk  ),0)
//                                  + pLast(PCELL_IDX(di,dj,dk,ii  ,jj_p,kk  ),0)
//                                  + pLast(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk_p),0);
//                const REAL p_smn0 = pCurr(PCELL_IDX(di,dj,dk,ii_n,jj  ,kk  ),0)
//                                  + pCurr(PCELL_IDX(di,dj,dk,ii  ,jj_n,kk  ),0)
//                                  + pCurr(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk_n),0);
//                const REAL p_smnm = pLast(PCELL_IDX(di,dj,dk,ii_n,jj  ,kk  ),0)
//                                  + pLast(PCELL_IDX(di,dj,dk,ii  ,jj_n,kk  ),0)
//                                  + pLast(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk_n),0);
//                p1 = -(pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0) - 2.*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0))
//                     -(lambda/6.)*(-p_spnm - p_smn0 + p_smnm)
//                     +(lambda/6.)*( p_smn0 - 6.*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0)); 
//                p1 *= 6.*lambda/(1.+4.*lambda); 
//                p1 += (1.-2.*lambda)/(1.+4.*lambda)*(2.*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0)
//                                                    -   pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0));
//                if (EQUAL_FLOATS(lambda,0.5))
//                    throw std::runtime_error("**WARNING** lambda is too close to 0.5, this formulation is singular");
//                REAL p_sp = -(p1 + pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0) - 2.*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0))
//                          - (lambda /6.)*(-p_spnm - p_smn0 + p_smnm)
//                          + (lambda2/3.)*( p_smn0 - 6.*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0));
//                p_sp *= (6./lambda/(1.-2.*lambda));
//                p_ii_p = 1./3.*p_sp; 
//                p_jj_p = 1./3.*p_sp; 
//                p_kk_p = 1./3.*p_sp; 
//#endif
                pCurr(PCELL_IDX(di,dj,dk,ii_p,jj  ,kk  ),0) = p_ii_p;
                pCurr(PCELL_IDX(di,dj,dk,ii  ,jj_p,kk  ),0) = p_jj_p;
                pCurr(PCELL_IDX(di,dj,dk,ii  ,jj  ,kk_p),0) = p_kk_p;
                pCurr(PCELL_IDX(di,dj,dk,ii_p,jj_p,kk_p),0) = 0.0;

                // fetch for possible contribution through interface (with other sim units)
                REAL p1_src = 0.0, p2_src = 0.0, p3_src = 0.0; 
                REAL dbuf; 
                REAL alpha1=-1.0, alpha2=-1.0, alpha3=-1.0; 
                for (auto &interface : _boundaryInterfaces)
                {
                    if (interface->GetDirection()==di)
                    {
                        const bool hasNeighbour = interface->GetOtherCellPressure(*grid_id, cell_idx, dbuf); 
                        if (hasNeighbour)
                        {
                            alpha1 = interface->GetBlendCoeff(simulationTime); 
                            p1_src += dbuf; 
                        }
                    }
                    else if (interface->GetDirection()==dj)
                    {
                        const bool hasNeighbour = interface->GetOtherCellPressure(*grid_id, cell_idx, dbuf); 
                        if (hasNeighbour)
                        {
                            alpha2 = interface->GetBlendCoeff(simulationTime); 
                            p2_src += dbuf; 
                        }
                    }
                    else if (interface->GetDirection()==dk)
                    {
                        const bool hasNeighbour = interface->GetOtherCellPressure(*grid_id, cell_idx, dbuf); 
                        if (hasNeighbour)
                        {
                            alpha3 = interface->GetBlendCoeff(simulationTime); 
                            p3_src += dbuf; 
                        }
                    }
                }
                // alpha-blend the neighbour pressure with ABC to smooth out discont.
                // then do a simple six-point Laplacian computation
                REAL p1_blend, p2_blend, p3_blend; 
                if      (alpha1<0)                   p1_blend = p_ii_p; 
                else if (alpha1>=0.0 && alpha1<=1.0) p1_blend = (1.0 - alpha1)*p_ii_p + alpha1*p1_src; 
                else                                 p1_blend = p1_src; 
                if      (alpha2<0)                   p2_blend = p_jj_p; 
                else if (alpha2>=0.0 && alpha2<=1.0) p2_blend = (1.0 - alpha2)*p_jj_p + alpha2*p2_src; 
                else                                 p2_blend = p2_src; 
                if      (alpha3<0)                   p3_blend = p_kk_p; 
                else if (alpha3>=0.0 && alpha3<=1.0) p3_blend = (1.0 - alpha3)*p_kk_p + alpha3*p3_src; 
                else                                 p3_blend = p3_src; 
                pNext(cell_idx,0) = 
                    + 2.0*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0) 
                    -     pLast(PCELL_IDX(di,dj,dk,ii,jj,kk),0)
                    + lambda2*(  
                        + p_ii_n  + p_jj_n + p_kk_n
                        + p1_blend + p2_blend + p3_blend
                        - 6.*pCurr(PCELL_IDX(di,dj,dk,ii,jj,kk),0));
            }
            else 
            {
                throw std::runtime_error("**ERROR** Cell " + std::to_string(cell_idx) 
                        + " has more than 3 boundary faces"); 
            }
        }
        else  // not boundary cells
        {
            // update normal cell
            if (_waveSolverSettings->useAirViscosity && _pressureCellHasValidHistory.at(cell_idx))
                pNext(cell_idx, 0) = 2.0*pCurr(cell_idx, 0) 
                                   -     pLast(cell_idx, 0) 
                                   - lastLaplacian(cell_idx, 0) *  kcalp
                                   + currLaplacian(cell_idx, 0) * (c2_k2+kcalp); 
            else
                pNext(cell_idx, 0) = 2.0*pCurr(cell_idx, 0) 
                                   -     pLast(cell_idx, 0) 
                                   + currLaplacian(cell_idx, 0) * c2_k2; 
        }
#undef PCELL_IDX

        // evaluate external sources only happens not in PML
        // Liu Eq (16) f6x term
        if (evaluateExternalSource)
            pNext(cell_idx, 0) += 
                _objects->EvaluatePressureSources(cell_position, 
                                                  cell_position, 
                                                  simulationTime+0.5*timeStep)*timeStep;
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
    const int N_ghostCells = _ghostCells.size() + _boundaryGhostCells.size(); 
    const auto &boundaryType = _waveSolverSettings->boundaryHandlingType; 
    if (N_ghostCells == 0) return; 

#ifdef USE_OPENMP
    const int N_max_threads = omp_get_max_threads();
#else
    const int N_max_threads = 1;
#endif

    // copy the indices for openmp access
    std::vector<GhostCell_Key>          gcIndices(_ghostCells.size()); 
    std::vector<BoundaryGhostCell_Key> bgcIndices(_boundaryGhostCells.size()); 
    int count=0; 
    for (auto &m : _ghostCells)
    {
        gcIndices[count] = m.first; 
        ++count; 
    }
    count=0; 
    for (auto &m : _boundaryGhostCells)
    {
        bgcIndices[count] = m.first; 
        ++count; 
    }

    std::unordered_map<GhostCell_Key, int> mapGC; // to index into matrices
    for (const auto &m : _ghostCells)
    {
        const int entry = mapGC.size(); 
        mapGC[m.first] = entry; 
    }
    Eigen::SparseMatrix<REAL, Eigen::ColMajor> matGC(N_ghostCells, N_ghostCells); // matrix
    Eigen::Matrix<REAL, Eigen::Dynamic, 1>     rhsGC(N_ghostCells);               // right hand side
    typedef Eigen::Triplet<double> Triplet; 
    std::vector<Triplet>                       entGC;                             // entries
    std::vector<std::vector<Triplet>>          threadGCEntries;                   // for openmp
    threadGCEntries.resize(N_max_threads); 
    for (auto &v : threadGCEntries)
        v.reserve(N_ghostCells*8);

#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int ii=0; ii<N_ghostCells; ++ii)
    {
        GhostCell_Key key_gc; 
        BoundaryGhostCell_Key key_bgc; 
#ifdef USE_OPENMP
            const int thread_idx = omp_get_thread_num(); 
#else
            const int thread_idx = 0; 
#endif

        if (ii<_ghostCells.size()) key_gc  =  gcIndices.at(ii); 
        else                       key_bgc = bgcIndices.at(ii - _ghostCells.size());
        auto &gc = (ii<_ghostCells.size() ? _ghostCells.at(key_gc) : _boundaryGhostCells.at(key_bgc)); 
        const Vector3d &cellPosition = gc->position; 

        // find neighbour pressure as reference
        REAL p_neig; 
        if (!gc->boundary) p_neig = p(gc->neighbourCell, 0); 
        else               gc->interface->GetOtherCellPressure(*grid_id, gc->ownerCell, p_neig);

        // if its plane constraint, its always perfectly reflecting
        if (gc->type == 1) 
        {
            const int row = mapGC.at(key_gc); 
            threadGCEntries.at(thread_idx).push_back(Triplet(row, row, -1.0)); // A(r,r) = 1.0
#pragma omp critical
            rhsGC(row) = -p_neig; 
            continue; 
        }
            
        // find the sample object, triangle and reflection point
        int boundaryObject, closestTriangle;
        Vector3d boundaryPoint, imagePoint, erectedNormal; 
        REAL distance; 
        _objects->LowestObjectDistance(cellPosition, distance, boundaryObject); 
        auto object = _objects->GetPtr(boundaryObject); 
        if (object->Type() == SHELL_OBJ)
        {
            // shell is deformable, first we collect a list of triangles in the neighbourhood of the 
            // cell, then use this list to find reflection point.
            // this effectively skips the kd-tree, since its only built for static object
            auto shell_object = std::dynamic_pointer_cast<FDTD_ShellObject>(object); 
            const int shell_object_id = std::stoi(shell_object->GetMeshName()); 
            IntArray cell_list;
            _pressureField.cell26Neighbours(gc->ownerCell, cell_list);
            cell_list.push_back(gc->ownerCell); 
            std::set<int> triangles;
            for (const auto &cell : cell_list)
                for (const auto tri_id : _cellTriangles.at(cell))
                    if (tri_id.objectID == shell_object_id)
                        triangles.insert(tri_id.triangleID); 
            closestTriangle = shell_object->ReflectAgainstBoundary(cellPosition,
                                                                   triangles,
                                                                   imagePoint,
                                                                   boundaryPoint, 
                                                                   erectedNormal, 
                                                                   distance); 
        }
        else if (   gc->cache 
                 && gc->cache->tid.objectID == boundaryObject)
        {
            closestTriangle = object->ReflectAgainstBoundary(cellPosition, imagePoint, boundaryPoint, erectedNormal, 
                                                             distance, gc->cache->tid.triangleID); 
            gc->cache->tid.triangleID = closestTriangle; 
        } 
        else 
        {
            closestTriangle = object->ReflectAgainstBoundary(cellPosition, imagePoint, boundaryPoint, erectedNormal, 
                                                             distance); 
            if (!gc->cache) 
                gc->cache = std::make_unique<GhostCell_Cache>(); 
            gc->cache->tid.objectID = boundaryObject;
            gc->cache->tid.triangleID = closestTriangle; 
        }
        const Tuple3ui triangle = object->GetMeshPtr()->triangle_ids(closestTriangle); 

        // sample the acceleration and compute the pressure at ghost cell position if
        // rasterization is used
        const REAL &h = _waveSolverSettings->cellSize; 
        Vector3d grad_p;
        for (int ii=0; ii<3; ++ii)
            grad_p += object->EvaluateBoundaryAcceleration(triangle[ii], simulationTime); 
        grad_p *= (-density/3.0); 
        const Vector3d cell_n = gc->CellNormal(); 
        const REAL grad_p_nr = grad_p.dotProduct(cell_n); 
        const REAL grad_p_ne = grad_p.dotProduct(erectedNormal); 
        const REAL pg = p_neig - h*grad_p_nr; 

        const REAL weights = (object->DistanceToMesh(cellPosition) < DISTANCE_TOLERANCE ? 
                -2.0*fabs(distance) : -fabs(distance));  // finite-difference weight
        REAL weightedPressure = grad_p_ne*weights;  // p_r - p_g = grad_n(p) * h

        // get the box enclosing the image point; 
        IntArray neighbours; 
        neighbours.reserve(8); 
        _pressureField.enclosingNeighbours(imagePoint, neighbours); 

        // figure out whether high-order boundary is needed
        // if not, then pressure(this ghost cell) = pg
        bool use_immerse_interface = 
            (neighbours.size() == 8 ) && 
            (object->Type() != SHELL_OBJ) &&
            (boundaryType == PML_WaveSolver_Settings::BoundaryHandling::FULLY_COUPLED);

        if (!use_immerse_interface)
        {
            const int row = mapGC.at(key_gc); 
            threadGCEntries.at(thread_idx).push_back(Triplet(row, row, 1.0)); 
            rhsGC(row) = pg; 
            continue; 
        }

        // high-order needed, do complex solve
        Eigen::MatrixXd V(8,8); // vandermonde matrix
        Eigen::MatrixXd v(1,8); // eval point
        Eigen::MatrixXd r(8,1); // rhs
        Eigen::MatrixXd w(8,1); // weights
        Eigen::MatrixXd c(8,1); // polynomial coefficients

        int                         s_i = -1; // indicate whether self is involved
        std::set<int>               s_k;      // set of known   cells in terms of local idx
        std::map<int,GhostCell_Key> s_u;      // map of unknown cells in terms of local idx

        // prepare for the coordinate transformation to
        // make the interpolation happens in [-1,1]^3
        Vector3d origin(0.,0.,0.); 
        for (int nn=0; nn<8; ++nn)
            origin += _pressureField.cellPosition(neighbours.at(nn));
        origin /= 8.0;

        const int &c_idx = gc->ownerCell; 
        FillVandermondeRegularS(0, imagePoint, v, origin, h);
        for (int nn=0; nn<8; ++nn)
        {
            const int n_idx = neighbours.at(nn); 
            if (n_idx == c_idx) 
            {
                FillVandermondeBoundaryS(nn, boundaryPoint, erectedNormal, V, origin, h); 
                r(nn) = grad_p_ne*h; // scaling due to coordinate transformation
                s_i = nn;
            }
            else if (!_isGhostCell.at(n_idx))
            {
                FillVandermondeRegularS( nn, pressureFieldPosition(n_idx), V, origin, h); 
                r(nn) = p(n_idx, 0); 
                s_k.insert(nn); 
            }
            else 
            {
                if (IsPressureCellSolid(n_idx))
                    throw std::runtime_error("**ERROR** solid!");
                // get nearest ghost cells to form interpolator
                if (_ghostCellsTree.at(n_idx).size() == 0) 
                    throw std::runtime_error("**ERROR** ghost cell pressure update failed due to insufficient children");
                REAL          min_d2  = std::numeric_limits<REAL>::max(); 
                GhostCell_Key min_key; 
                for (auto &key : _ghostCellsTree.at(n_idx))
                {
                    auto &gc_n = _ghostCells.at(key); // neighbour ghost cells
                    const REAL d2 = (gc_n->position - imagePoint).normSqr();
                    if (d2 < min_d2)
                    {
                        min_d2  = d2; 
                        min_key = key; 
                    }
                }
                FillVandermondeRegularS( nn, _ghostCells.at(min_key)->position, V, origin, h); 
                s_u[nn] = min_key; 
            }
        }
        v.transposeInPlace(); 
        V.transposeInPlace(); 

        Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::FullPivHouseholderQRPreconditioner> svd(V,
                Eigen::ComputeFullU | Eigen::ComputeFullV);
        const double cond = svd.singularValues()(0)
            / svd.singularValues()(svd.singularValues().size()-1);

        w = svd.solve(v); 
        // if the solution is not good, fall back 
        if ((cond>CONDN_LIM) && s_i != -1)
        {
            use_immerse_interface = false; 
            assert(s_i != -1); 
            w.setZero(); 
            r.setZero(); 
            s_k.clear(); 
            s_u.clear(); 
            weightedPressure = pg;
        }

        // either way, w should be well-posed now
        for (int nn=0; nn<8; ++nn)
            assert(std::isfinite(w(nn)));

        // form linear system 
        REAL rhs = -weightedPressure; 
        const int row = mapGC.at(key_gc); 
        threadGCEntries.at(thread_idx).push_back(Triplet(row, row, -1.0)); // A(r,r) = -1.0
        for (const int  &i_k : s_k)
            rhs -= w(i_k)*r(i_k); 
        if (s_i != -1)
            rhs -= w(s_i)*r(s_i); 
        for (const auto &m_u : s_u)
        {
            const int i_k = m_u.first; 
            const int col = mapGC.at(m_u.second); 
            threadGCEntries.at(thread_idx).push_back(Triplet(row, col, w(i_k))); 
        }
        rhsGC(row) = rhs; 
    }
    // actually fill the sparse matrix
    entGC.reserve(threadGCEntries.at(0).size()*N_max_threads); 
    for (const auto &v : threadGCEntries) 
        entGC.insert(entGC.end(), v.begin(), v.end()); 
    matGC.setFromTriplets(entGC.begin(), entGC.end()); 
    
    // linear solve
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> slnGC(N_ghostCells); 
    Eigen::BiCGSTAB<Eigen::SparseMatrix<REAL>> solver; 
    solver.compute(matGC); 
    slnGC = solver.solve(rhsGC); 
    std::cout << "norm(rhs)        = " << rhsGC.norm()        << std::endl; 
    std::cout << "#iterations      = " << solver.iterations() << std::endl; 
    std::cout << "#estimated error = " << solver.error()      << std::endl; 
    std::cout << "sol_max ="   << slnGC.maxCoeff() 
              << "; sol_min =" << slnGC.minCoeff() << std::endl;
    //Eigen::SparseLU<Eigen::SparseMatrix<REAL, Eigen::ColMajor> > solver; 
    //solver.analyzePattern(matGC); 
    //solver.factorize(matGC); 
    //slnGC = solver.solve(rhsGC); 
    for (auto &m : _ghostCells)
    {
        auto &gc = m.second; 
        gc->pressure = slnGC(mapGC.at(gc->MakeKey()));
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
        const int closestTriangleIndex = object->ReflectAgainstBoundary(cellPosition, imagePoint, boundaryPoint, erectedNormal, distance); 
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
#ifdef USE_FV
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
        //auto gcm = loop_array[gc_i];
        auto gcm = _ghostCellsCollection.begin(); 
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
                    assert(indicesBuffer[dim_0]>=0 && indicesBuffer[dim_0]<_pressureField.cellDivisions()[dim_0]); 
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
                
                const REAL dp_dn_face = (p_neighbour - gc->values.at(gc->valuePointer).at(sp->gc_value_idx))/h;
                //volume += ((sp->position-cellPosition)/3.0).dotProduct(sp->normal);
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
#endif

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
    const int k = std::min<int>((int)((offset - lowerCorner)/_waveSolverSettings->cellSize), divs[dim]-1); 

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
    //_ghostCellsInverse.clear(); 
    //for (size_t ii=0; ii<_ghostCells.size(); ii++)
    //    _ghostCellsInverse[_ghostCells[ii]] = ii; 
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
            throw std::runtime_error("**ERROR** Fresh cell inside some object. This shouldn't happen for pressure cells. objectID: " + std::to_string(objectID)); 
        Vector3d imagePoint, boundaryPoint, erectedNormal; 
        REAL distanceTravelled; 
        auto object = _objects->GetPtr(objectID); 
        int closestTriangle;
        if (object->Type() == SHELL_OBJ)
        {
            // this procedure is the same as in pressureUpdateGhostCells
            auto shell_object = std::dynamic_pointer_cast<FDTD_ShellObject>(object); 
            const int shell_object_id = std::stoi(shell_object->GetMeshName()); 
            IntArray cell_list;
            _pressureField.cell26Neighbours(cell_idx, cell_list);
            cell_list.push_back(cell_idx); 
            std::set<int> triangles;
            for (const auto &cell : cell_list)
                for (const auto tri_id : _cellTriangles.at(cell))
                    if (tri_id.objectID == shell_object_id)
                        triangles.insert(tri_id.triangleID); 
            closestTriangle = shell_object->ReflectAgainstBoundary(cellPosition,
                                                                   triangles,
                                                                   imagePoint,
                                                                   boundaryPoint, 
                                                                   erectedNormal, 
                                                                   distanceTravelled); 
        }
        else
        {
            closestTriangle = object->ReflectAgainstBoundary(cellPosition, imagePoint, boundaryPoint, erectedNormal, distanceTravelled); 
        }
        // prepare interpolation stencils, this part is similar to
        // the vandermonde part in ghost cell pressure update
        IntArray neighbours; 
        _pressureField.enclosingNeighbours(imagePoint, neighbours); 

#if 0 // Mittal
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

#endif
#if 1 // MLS
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
                REAL min_d2 = std::numeric_limits<REAL>::max(); 
                int  min_id = -1; 
                for (const int & n : neighbours)
                {
                    const REAL d2 = (imagePoint - _pressureField.cellPosition(n)).normSqr();
                    if (d2 < min_d2)
                    {
                        min_d2 = d2; 
                        min_id = n; 
                    }
                }
                assert(min_id != -1);
                p(cell_idx, 0) = p(min_id, 0); 
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
            const Tuple3ui triangle = object->GetMeshPtr()->triangle_ids(closestTriangle); 
            Vector3d acc; 
            for (int ii=0; ii<3; ++ii)
                acc += object->EvaluateBoundaryAcceleration(triangle[ii], simulationTime); 
            acc /= 3.0; 
            const REAL an = acc.dotProduct(erectedNormal); 
            p(cell_idx, 0) = mlsVal(0, 0) + density*an*fabs(distanceTravelled);
        }
#endif
#if 0
        REAL p1 = 0.0; 
        REAL p2 = 0.0;
        REAL n  = 0.0;
        Tuple3i indices = _pressureField.cellIndex(cell_idx); 
        int n_idx; 
        for (int dd=0; dd<3; ++dd)
        {
            indices[dd] += 2; 
            n_idx = _pressureField.cellIndex(indices); 
            p2 += p(n_idx, 0); 
            indices[dd] -= 4; 
            n_idx = _pressureField.cellIndex(indices); 
            p2 += p(n_idx, 0); 
            indices[dd] += 2; 
        }
        for (int nn=0; nn<neighbours.size(); ++nn)
        {
            if (_isBulkCell.at(neighbours.at(nn)))
            {
                p1 += p(neighbours.at(nn), 0); 
                n += 1.0; 
            }
        }
        p(cell_idx, 0) = (1./pow(n,2) + 2./n)*p1 - 1./n*p2; 
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
        int pmlDirection; 
        if (InsidePML(cellPos, _PML_absorptionWidth, pmlDirection))
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
           
            cell.position = cellPos; 
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
                //_ghostCells.push_back(cell_idx);  // FIXME debug

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

            int flag;
            if (InsidePML(cellPosition, _PML_absorptionWidth, flag))
            {
                const REAL absorptionCoefficient = PML_absorptionCoefficient(cellPosition, _PML_absorptionWidth, dimension); 
                const REAL updateCoefficient = PML_velocityUpdateCoefficient(absorptionCoefficient, dt);
                const REAL gradientCoefficient  = PML_pressureGradientCoefficient(absorptionCoefficient, dt, dx, rho);
                PML_VelocityCell cell; 
                cell.index = cell_idx; 
                cell.dimension = dimension; 
                cell.updateCoefficient = updateCoefficient; 
                cell.gradientCoefficient = gradientCoefficient; 
                cell.position = cellPosition;

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
                    //_ghostCells.push_back(cell_idx);  // FIXME debug
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
    int N = _objects->N(); 
    std::vector<ScalarField::RangeIndices> indices(N); 
    auto &objects = _objects->GetRigidObjects(); 
    int bbox_id = 0;
    for (auto &m : objects)
    {
        const FDTD_MovableObject::BoundingBox &unionBBox = m.second->GetUnionBBox();
        const Vector3d maxBound = unionBBox.maxBound + 2.0*_cellSize; 
        const Vector3d minBound = unionBBox.minBound - 2.0*_cellSize; 
        _pressureField.GetIterationBox(minBound, maxBound, indices[bbox_id]); 
        ++ bbox_id; 
    } 


    // FIXME debug START
    N = 1; 
    indices.clear(); 
    ScalarField::RangeIndices tmp; 
    tmp.startIndex.set(0,0,0); 
    tmp.dimensionIteration.set(_pressureField.cellDivisions()[0],
                               _pressureField.cellDivisions()[1],
                               _pressureField.cellDivisions()[2]); 
    indices.push_back(std::move(tmp)); 
    // FIXME debug END


    std::cout << "test 7-1\n"; // FIXME debug
    std::fill(_isBulkCell.begin(), _isBulkCell.end(), true); // FIXME debug
    std::fill(_isGhostCell.begin(), _isGhostCell.end(), false); // FIXME debug 
    std::fill(_containingObject.begin(), _containingObject.end(), -1); // FIXME debug

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

    std::cout << "test 7-2\n"; // FIXME debug

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

    std::cout << "test 7-3\n"; // FIXME debug

    // push back ghost cells
    for (const IntArray &gcIndThread : ghostCellsThreads)
        for (const int &gcInd : gcIndThread)
        {
            //_ghostCells.push_back(gcInd); // FIXME 
            _ghostCellsChildren.push_back(childrenIndex); 
        }
    std::cout << "test 7-4\n"; // FIXME debug
    ComputeGhostCellInverseMap(); 

    std::cout << "test 7-5\n"; // FIXME debug
    // Classify velocity cells
    pGCFull.clear(); 
    _ghostCellParents.clear(); 
    _ghostCellChildArrayPositions.clear(); 
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
    std::cerr << "test 7-6\n"; // FIXME debug
    for (int dimension = 0; dimension < 3; dimension++)
    {
        // FIXME debug START
        N = 1; 
        indices.clear(); 
        ScalarField::RangeIndices tmp; 
        tmp.startIndex.set(0,0,0); 
        tmp.dimensionIteration.set(_velocityField[dimension].cellDivisions()[0],
                                   _velocityField[dimension].cellDivisions()[1],
                                   _velocityField[dimension].cellDivisions()[2]); 
        indices.push_back(std::move(tmp)); 
        // FIXME debug END
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
try{
                if (_classified.at(cell_idx)) continue;
                if (IsVelocityCellSolid(cell_idx, dimension)) 
                    _velocityCellHasValidHistory[dimension].at(cell_idx) = false; 
                else
                    _velocityCellHasValidHistory[dimension].at(cell_idx) = true; 
} catch (std::out_of_range &e) {
    std::cerr << "seg 1\n"; 
    exit(1); 
}

                Tuple3i cell_coordinates = _velocityField[ dimension ].cellIndex( cell_idx );
try{
                // We will assume that boundary cells are not interfacial
                if ( cell_coordinates[ dimension ] == 0
                  || cell_coordinates[ dimension ] == _velocityField[ dimension ].cellDivisions()[ dimension ] - 1 )
                {
                    _isVelocityBulkCell[dimension].at(cell_idx) = true; 
                    _isVelocityInterfacialCell[dimension].at(cell_idx) = false; 
                    continue;
                }
} catch (std::out_of_range &e) {
    std::cerr << "seg 2\n"; 
    exit(1); 
}


                // Look at our neighbours in the pressure field
                cell_coordinates[ dimension ] -= 1;
                const int pressure_cell_idx1 = _pressureField.cellIndex( cell_coordinates );
                cell_coordinates[ dimension ] += 1;
                const int pressure_cell_idx2 = _pressureField.cellIndex( cell_coordinates );
try{
} catch (std::out_of_range &e) {
    std::cerr << "seg 3\n"; 
    exit(1); 
}


                GhostCellInfo ghostCellInfo;
                bool infoUsed = false; 
try{
                if (_isBulkCell.at(pressure_cell_idx1) && _isBulkCell.at(pressure_cell_idx2))
                {
                    // Both pressure cell neighbours are bulk cells, so this is
                    // a bulk cell too
                    _isVelocityBulkCell[dimension].at(cell_idx) = true; 
                    _isVelocityInterfacialCell[dimension].at(cell_idx) = false; 
                }
                else if (   _isBulkCell.at(pressure_cell_idx1)
                        && !_isBulkCell.at(pressure_cell_idx2) )
                {
                    // Only one neighbour is inside the domain, so this must
                    // be an interfacial cell
                    _isVelocityBulkCell[dimension].at(cell_idx) = false; 
                    _isVelocityInterfacialCell[ dimension ].at(cell_idx) = true;

                    // Get the object ID for the boundary that we are adjacent to
                    const int boundaryObject = _containingObject.at(pressure_cell_idx2);
                    TRACE_ASSERT( boundaryObject >= 0 );

                    // Determine a scaling coefficient based on the angle between
                    // the boundary normal and the rasterized boundary normal
                    Vector3d normal( 0.0, 0.0, 0.0 );
                    Vector3d position = _velocityField[ dimension ].cellPosition( cell_idx );
                    Vector3d gradient; 
try {
                    _objects->ObjectNormal(boundaryObject, position, gradient); 
} catch (std::out_of_range &e) {
    std::cerr << "seg 4-1\n"; 
    std::cerr << e.what() << std::endl;
    COUT_SDUMP(boundaryObject);
    exit(1); 
}
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
                else if ( !_isBulkCell.at(pressure_cell_idx1)
                        && _isBulkCell.at(pressure_cell_idx2) )
                {
                    // Only one neighbour is inside the domain, so this must
                    // be an interfacial cell
                    _isVelocityBulkCell[dimension].at(cell_idx) = false; 
                    _isVelocityInterfacialCell[dimension].at(cell_idx) = true;

                    // Get the object ID for the boundary that we are adjacent to
                    int boundaryObject = _containingObject.at(pressure_cell_idx1);
                    TRACE_ASSERT( boundaryObject >= 0 );

                    // Determine a scaling coefficient based on the angle between
                    // the boundary normal and the rasterized boundary normal
                    Vector3d normal( 0.0, 0.0, 0.0 );
                    Vector3d position = _velocityField[ dimension ].cellPosition( cell_idx );
                    Vector3d gradient; 
try{
                    _objects->ObjectNormal(boundaryObject, position, gradient); 
} catch (std::out_of_range &e) {
    std::cerr << "seg 4-2\n"; 
    std::cerr << e.what() << std::endl;
    COUT_SDUMP(boundaryObject);
    exit(1); 
}
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
                    _isVelocityBulkCell[dimension].at(cell_idx) = false; 
                    _isVelocityInterfacialCell[dimension].at(cell_idx) = false;
                    v[dimension](cell_idx, 0) = 0.0; // clear solid velocity cell
                    numSolidCells[dimension] ++; 
                }
} catch (std::out_of_range &e) {
    std::cerr << "seg 4\n"; 
    std::cerr << e.what() << std::endl;
    exit(1); 
}

try{
                if (infoUsed) 
                    ghostCellInfoThreads.at(thread_idx).push_back(ghostCellInfo); 
                _classified.at(cell_idx) = true; // toggle on to prevent reclassification
} catch (std::out_of_range &e) {
    std::cerr << "seg 5\n"; 
    exit(1); 
}
            }
        }
try{
        SetClassifiedSubset(_velocityField[dimension], N, indices, false); 

} catch (std::out_of_range &e) {
    std::cerr << "seg 6\n"; 
    exit(1); 
}
    }
    std::cout << "test 7-7\n"; // FIXME debug
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
        printf( "MAC_Grid: classifyCellsDynamic_FAST:\n" );
        printf( "\tFound %d ghost cells\n", (int)_ghostCells.size() );
        printf( "\tFound %d v_x interfacial cells\n", (int)_velocityInterfacialCells[ 0 ].size() );
        printf( "\tFound %d v_y interfacial cells\n", (int)_velocityInterfacialCells[ 1 ].size() );
        printf( "\tFound %d v_z interfacial cells\n", (int)_velocityInterfacialCells[ 2 ].size() );
        printf( "\tFound %d v_x solid cells\n", numSolidCells[0] );
        printf( "\tFound %d v_y solid cells\n", numSolidCells[1] );
        printf( "\tFound %d v_z solid cells\n", numSolidCells[2] );
    }

    std::cout << "test 7-7\n"; // FIXME debug
    ClearUnusedCache();
}

//##############################################################################
// Function classifyCells_FAST
//   This function classify cells. 
//   NOTE: Only works for collocated scheme + non-pml impl. 
//#############################################################################
void MAC_Grid::classifyCells_FAST(MATRIX (&pCollocated)[3], const bool &verbose)
{
    // set history valid for interpolation
    const int numCells = _pressureField.numCells(); 
    std::copy(_isBulkCell.begin(), _isBulkCell.end(), 
              _pressureCellHasValidHistory.begin()); 

    // set default field values
    std::fill(_isBulkCell.begin(), _isBulkCell.end(), true); 
    std::fill(_isGhostCell.begin(), _isGhostCell.end(), false); 
    std::fill(_classified.begin(), _classified.end(), false); 
    if (_cellTriangles.size() != numPressureCells())
        _cellTriangles.resize(numPressureCells()); 

    // copy the cache to map and clear ghost cells
    for (auto &m : _ghostCells)
        _ghostCellsCached.emplace(
                std::make_pair(m.second->MakeKey(), std::move(m.second->cache))); 
    _ghostCells.clear(); 
    _boundaryGhostCells.clear(); // boundary ghost cells utilize no cache

    // we only have to check bounding box around each objects
    // first get the bounding box
    const int N = _objects->N(); 
    auto &objects = _objects->GetRigidObjects(); 
    int count = 0;
#if 1 // faster 
    std::vector<ScalarField::RangeIndices> bbox_rast(N); 
    const auto &constraints = _objects->GetConstraints(); 
    for (auto &m : objects)
    {
        const FDTD_MovableObject::BoundingBox &unionBBox = m.second->GetUnionBBox();
        Vector3d maxBound = unionBBox.maxBound; 
        Vector3d minBound = unionBBox.minBound; 
        // exclude volumes covered by all the constraints (where sdf is negative)
        for (const auto &c : constraints) 
        {
            const auto &plane = c.second; 
            const int nmlDir = plane->NormalDirection(); 
            const int nmlSgn = plane->NormalSign(); 
            if (nmlSgn > 0) minBound[nmlDir] = std::max(minBound[nmlDir], plane->Height()); 
            else            maxBound[nmlDir] = std::min(maxBound[nmlDir], plane->Height()); 
        }
        _pressureField.GetIterationBox(minBound, maxBound, bbox_rast[count]); 
        ++ count; 
    } 
    // added a layer around constraint so that ghost cells can be detected
    const REAL NEG_INF = std::numeric_limits<REAL>::lowest(); 
    const REAL POS_INF = std::numeric_limits<REAL>::max(); 
    for (const auto &c : constraints)
    {
        const auto &plane = c.second; 
        const int nmlDir = plane->NormalDirection(); 
        const int nmlSgn = plane->NormalSign(); 
        Vector3d minBound(NEG_INF, NEG_INF, NEG_INF);
        Vector3d maxBound(POS_INF, POS_INF, POS_INF);
        if (nmlSgn > 0)
            maxBound[nmlDir] = plane->Height(); 
        else
            minBound[nmlDir] = plane->Height(); 
        ScalarField::RangeIndices plane_box; 
        _pressureField.GetIterationBox(minBound, maxBound, plane_box); 
        bbox_rast.push_back(plane_box); 
    }
#else 
    std::vector<ScalarField::RangeIndices> bbox_rast(1); 
    _pressureField.GetIterationBox(PressureBoundingBox().minBound(), 
                                   PressureBoundingBox().maxBound(), 
                                   bbox_rast[0]); 
#endif
    // helper lambdas
    auto isFluid = [&,this](const Vector3d &pos, int &type){
        const bool isfluid = !(_objects->OccupyByConstraint(pos))
            && _objects->OccupyByObject(pos)<0
            && _objects->OccupyByObject(pos+Vector3d(1.,0.,0.)*0.25*_waveSolverSettings->cellSize)<0 
            && _objects->OccupyByObject(pos-Vector3d(1.,0.,0.)*0.25*_waveSolverSettings->cellSize)<0 
            && _objects->OccupyByObject(pos+Vector3d(0.,1.,0.)*0.25*_waveSolverSettings->cellSize)<0 
            && _objects->OccupyByObject(pos-Vector3d(0.,1.,0.)*0.25*_waveSolverSettings->cellSize)<0 
            && _objects->OccupyByObject(pos+Vector3d(0.,0.,1.)*0.25*_waveSolverSettings->cellSize)<0 
            && _objects->OccupyByObject(pos-Vector3d(0.,0.,1.)*0.25*_waveSolverSettings->cellSize)<0;
        if      (!isfluid &&   _objects->OccupyByConstraint(pos))  type = 1;
        else if (!isfluid && !(_objects->OccupyByConstraint(pos))) type = 0;
        return isfluid;};

    // initialize threading containers
#ifdef USE_OPENMP
    const int N_max_threads = omp_get_max_threads();
#else
    const int N_max_threads = 1;
#endif
    std::vector<GhostCellType>           thread_GCMap(N_max_threads); 
    std::vector<BoundaryGhostCellType>   thread_BGCMap(N_max_threads); 
    std::vector<unordered_map<int, int>> thread_GCType(N_max_threads);
    std::vector<IntArray> thread_candidate_cells(N_max_threads); 
    for (auto &tcc : thread_candidate_cells)
        tcc.reserve(numCells/4); 

    // loop through each of the box and set bulk cells
    for (const auto &bbox : bbox_rast)
    {
        const Vector3i &start = bbox.startIndex; 
        const Vector3i &range = bbox.dimensionIteration; 
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
            const int cell_idx = _pressureField.cellIndex(cellIndices); 
            if (_classified.at(cell_idx)) continue;
            const Vector3d pos = _pressureField.cellPosition(cell_idx);
            int type; 
            bool isSolid = false; 
            std::set<TriangleIdentifier, TIComp> out_tris_shell; 
            if (!isFluid(pos, type)) // for rigid objects
            {
                isSolid = true; 
                thread_GCType.at(thread_idx)[cell_idx] = type;
            }
            else if (_objects->TriangleCubeIntersection( // for shells
                       pos, Vector3d(_waveSolverSettings->cellSize/2.0,
                                     _waveSolverSettings->cellSize/2.0,
                                     _waveSolverSettings->cellSize/2.0),
                       out_tris_shell))
            {
                isSolid = true; 
                _cellTriangles.at(cell_idx) = std::move(out_tris_shell);
                thread_GCType.at(thread_idx)[cell_idx] = 2;
            }

            if (isSolid)
            {
                _isBulkCell.at(cell_idx) = false; 
                thread_candidate_cells.at(thread_idx).push_back(cell_idx); 
            }
            _classified.at(cell_idx) = true; 
        }
    }

    // concatenate results from different threads
    for (int tid=1; tid<N_max_threads; ++tid)
        thread_candidate_cells.at(0).insert(thread_candidate_cells.at(0).end(),
                                            thread_candidate_cells.at(tid).begin(),
                                            thread_candidate_cells.at(tid).end()); 
    for (int tid=1; tid<N_max_threads; ++tid)
        thread_GCType.at(0).insert(thread_GCType.at(tid).begin(),
                                   thread_GCType.at(tid).end()); 

    // loop through all non-bulk and find ghost cells
    const IntArray &candidate_cells = thread_candidate_cells.at(0); 
    const auto     &gcTypes         = thread_GCType.at(0);
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int ii=0; ii<candidate_cells.size(); ++ii)
    {
        const int cell_idx = candidate_cells.at(ii);
        pCollocated[0](cell_idx,0) = 0.0; // set all pressure data to zero
        pCollocated[1](cell_idx,0) = 0.0;
        pCollocated[2](cell_idx,0) = 0.0;
        IntArray neighbours, topology; 
        neighbours.reserve(6); topology.reserve(6); 
        _pressureField.cellNeighbours(cell_idx, neighbours, topology); 
        int is_ghost = false; 
        for (int ii=0; ii<neighbours.size(); ++ii)
        {
            const int ns = neighbours.at(ii); 
            const int to = topology.at(ii); 
            if (_isBulkCell.at(ns))
            {
                is_ghost = true; 
                //int gctype; isFluid(_pressureField.cellPosition(cell_idx), gctype); 
                //auto gc = MakeGhostCell(cell_idx, ns, to, gctype); 
                auto gc = MakeGhostCell(cell_idx, ns, to, gcTypes.at(cell_idx)); 
#ifdef USE_OPENMP
                const int thread_idx = omp_get_thread_num(); 
#else
                const int thread_idx = 0; 
#endif
                thread_GCMap.at(thread_idx)[gc->MakeKey()] = std::move(gc); 
            } 
        }
        _isGhostCell.at(cell_idx) = is_ghost; 
        // if this cell is boundary, then need to check across interface if existed
        if (_pressureField.isBoundaryCell(cell_idx))
        {
            const Tuple3i nml = _pressureField.boundaryCellNormal(cell_idx); 
            int otherCell, topology = 0; 
            for (int dd=0; dd<3; ++dd)
            {
                for (auto &interface : _boundaryInterfaces)
                {
                    if (interface->GetDirection()!=dd)
                        continue; 
                    assert(grid_id); 
                    const bool hasNeighbour = interface->GetOtherCell(*grid_id, cell_idx, otherCell); 
                    const auto &otherGrid = interface->GetSimUnit(
                            interface->GetOtherSolver(*grid_id))->simulator->GetGrid(); 
                    if (hasNeighbour && otherGrid.IsPressureCellBulk(otherCell))
                    {
                        _isGhostCell.at(cell_idx) = true; 

                        if      (nml[0] == -1) topology = -1; // boundary at x=0
                        else if (nml[0] ==  1) topology =  1; // boundary at x=l
                        if      (nml[1] == -1) topology = -2; // boundary at y=0
                        else if (nml[1] ==  1) topology =  2; // boundary at y=l
                        if      (nml[2] == -1) topology = -3; // boundary at z=0
                        else if (nml[2] ==  1) topology =  3; // boundary at z=l
                        int gctype; isFluid(_pressureField.cellPosition(cell_idx), gctype); 
                        auto gc = MakeGhostCell(cell_idx, otherCell, topology, gctype); 
                        gc->boundary = true;
                        gc->ownerSolverId = *grid_id; 
                        gc->neighbourSolverId = interface->GetOtherSolver(*grid_id); 
                        gc->interface = interface; 
#ifdef USE_OPENMP
                        const int thread_idx = omp_get_thread_num(); 
#else
                        const int thread_idx = 0; 
#endif
                        thread_BGCMap.at(thread_idx)[gc->MakeKey_Boundary()] = std::move(gc); 
                    }
                }
            }
        }
    }
    // combine containers for each thread
    using IT_GC = GhostCellType::iterator; 
    using MV_IT = std::move_iterator<IT_GC>;
    for (auto &thread_map : thread_GCMap)
        _ghostCells.insert(MV_IT(thread_map.begin()), MV_IT(thread_map.end())); 
    _ghostCellsTree.clear(); 
    for (auto &m : _ghostCells)
        _ghostCellsTree[m.second->ownerCell].insert(m.first); 
    using IT_BGC = BoundaryGhostCellType::iterator; 
    using MV_BIT = std::move_iterator<IT_BGC>;
    for (auto &thread_map : thread_BGCMap)
        _boundaryGhostCells.insert(MV_BIT(thread_map.begin()), MV_BIT(thread_map.end())); 
    // any member in the map is not being identified at the current step, 
    // therefore we should remove all of them (including null ptrs)
    _ghostCellsCached.clear();
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
    auto &objects = _objects->GetRigidObjects(); 
    for (auto &m : objects)
    {
        auto &object = m.second; 
        std::shared_ptr<TriangleMesh<REAL> > mesh = object->GetMeshPtr(); 
        std::shared_ptr<TriangleMeshKDTree<REAL> > meshkd = std::dynamic_pointer_cast<TriangleMeshKDTree<REAL> >(mesh); 
        const int N_triangles = mesh->num_triangles(); 
        for (int t_idx=0; t_idx<N_triangles; ++t_idx)
        {
            const Vector3d centroid = object->ObjectToWorldPoint(meshkd->TriangleCentroid(t_idx)); 
            const Vector3d normal   = object->ObjectToWorldVector(meshkd->triangle_normal(t_idx)); 
            const int cell_idx = InPressureCell(centroid); 
            const TriangleIdentifier tri_id(m.first, t_idx, centroid, normal);
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
            auto gc = std::make_shared<GhostCell_Deprecated>(cell_idx, hashedTriangles); 
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

REAL MAC_Grid::PressureCellType(const int &idx, const BoundingBox *sceneBox) const
{
    if (_isPMLCell.at(idx))
    {
        const Vector3d cellPos = _pressureField.cellPosition(idx); 
        REAL maxAbsorptionCoefficient = std::numeric_limits<REAL>::min(); 
        for (int dim=0; dim<3; ++dim)
        {
            const REAL coeff = PML_absorptionCoefficient(cellPos, _PML_absorptionWidth, dim, sceneBox); 
            maxAbsorptionCoefficient = std::max<REAL>(maxAbsorptionCoefficient, coeff); 
        }
        return -0.5*maxAbsorptionCoefficient/_waveSolverSettings->PML_strength; 
    }

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

    for (auto &gg : _ghostCells)
    {
        position = gg.second->position; 
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
bool MAC_Grid::InsidePML(const Vector3d &x, const REAL &absorptionWidth, int &flag, const BoundingBox *sceneBox)
{
    const int &preset = _waveSolverSettings->boundaryConditionPreset; 
    const BoundingBox &bbox = (sceneBox ? *sceneBox : _pressureField.bbox());
    const REAL pmlWidth = absorptionWidth + 1E-12; // add buffer for detection

    for (int dimension=0; dimension<3; ++dimension)
    {
        const REAL hMin = x[ dimension ] - bbox.minBound()[ dimension ];
        const REAL hMax = bbox.maxBound()[ dimension ] - x[ dimension ];
        switch (preset)
        {
            case 0: // no wall
                if (hMin <= pmlWidth) {
                    flag = dimension; return true; 
                }
                else if (hMax <= pmlWidth) {
                    flag = dimension; return true;
                }
                break; 
            case 1: // wall on +x, +z
                if (hMin <= pmlWidth || (dimension == 1 && hMax <= pmlWidth)) {
                    flag = dimension; return true; 
                }
                break; 
            case 2: 
                if (dimension == 2 && hMax <= pmlWidth) {
                    flag = dimension; return true;  
                }
                break; 
            default: 
                break; 
        }
    }
    return false;
}

REAL MAC_Grid::PML_absorptionCoefficient(const Vector3d &x, REAL absorptionWidth, int dimension, const BoundingBox *sceneBox) const
{
#if 1
    const int &preset = _waveSolverSettings->boundaryConditionPreset; 
    const BoundingBox &bbox = (sceneBox ? *sceneBox : _pressureField.bbox());
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
        case 1: // wall on +x, +z
            if (hMin <= absorptionWidth)
                return _PML_absorptionStrength * pow(distMin,2) / a2; 
            else if (dimension == 1 && hMax <= absorptionWidth)
                return _PML_absorptionStrength * pow(distMax,2) / a2; 
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
    cell.indices = _pressureField.cellIndex(cellIndex); 
    cell.centroidPosition = _pressureField.cellPosition(cellIndex); 
    cell.lowerCorner = cell.centroidPosition - _waveSolverSettings->cellSize/2.0; 
    cell.upperCorner = cell.centroidPosition + _waveSolverSettings->cellSize/2.0; 

    if (EQUAL_FLOATS(cell.r_identity, 0.5))
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
        int count=0;
        for (int dim=0; dim<3; ++dim)
        {
            for (int dir=-1; dir<=1; dir+=2)
            {
                Tuple3i neighbour = cell.indices; 
                neighbour[dim] += dir; 
                const int neighbour_idx = _pressureField.cellIndex(neighbour); 
                auto key = GhostCell::MakeKey(cellIndex, neighbour_idx); 
                if (_ghostCells.find(key) != _ghostCells.end())
                    cell.gcValue[count] = _ghostCells.at(key)->pressure; 
            }
        }
#endif
    }
    else 
    {
        cell.s_identity = std::string("Pressure PML Cell"); 
    }

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

    const int xUpper = std::min<int>(indicesBuffer.x + 1, _pressureField.cellDivisions()[0]); 
    const int yUpper = std::min<int>(indicesBuffer.y + 1, _pressureField.cellDivisions()[1]); 
    const int zUpper = std::min<int>(indicesBuffer.z + 1, _pressureField.cellDivisions()[2]); 

    // upper velocity cell
    cell.vx[1] = v[0](_velocityField[0].cellIndex(Tuple3i(xUpper          , indicesBuffer[1], indicesBuffer[2])), 0); 
    cell.vy[1] = v[1](_velocityField[1].cellIndex(Tuple3i(indicesBuffer[0], yUpper          , indicesBuffer[2])), 0); 
    cell.vz[1] = v[2](_velocityField[2].cellIndex(Tuple3i(indicesBuffer[0], indicesBuffer[1], zUpper          )), 0); 

    cell.h = _waveSolverSettings->cellSize; 
}

void MAC_Grid::SetClassifiedSubset(const ScalarField &field, const int &N, const std::vector<ScalarField::RangeIndices> &indices, const bool &state)
{
    std::fill(_classified.begin(), _classified.end(), state); 
    //for (int bbox_id=0; bbox_id<N; ++bbox_id){ 
    //    const Vector3i &start = indices[bbox_id].startIndex;
    //    const Vector3i &range = indices[bbox_id].dimensionIteration;
    //    for (int ii=start.x; ii<start.x+range.x; ++ii){
    //    for (int jj=start.y; jj<start.y+range.y; ++jj){
    //    for (int kk=start.z; kk<start.z+range.z; ++kk){
    //        const Tuple3i cellIndices(ii,jj,kk);
    //        const int cell_idx = field.cellIndex(cellIndices); 
    //        _classified.at(cell_idx) = state;}}}}
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
    _ghostCellsChildren.at(_ghostCellsInverse.at(info.ghostCellParent)).at(info.childArrayPosition) = gcIndex; 
    _ghostCellChildArrayPositions.push_back(info.childArrayPosition); 
    _velocityInterfacialCells[info.dim].push_back(info.cellIndex); // not used in collocated
    _interfacialBoundaryIDs[info.dim].push_back(info.boundaryObject); // not used in collocated
    _interfacialBoundaryDirections[info.dim].push_back(info.boundaryDirection); // not used in collocated
    _interfacialBoundaryCoefficients[info.dim].push_back(info.boundaryCoefficient); // not used in collocated
    _interfacialGhostCellID[info.dim].push_back(gcIndex);  // not used in collocated
    _ghostCellParents.push_back(info.ghostCellParent);
    _ghostCellPositions.push_back(info.ghostCellPosition); 
    _ghostCellBoundaryIDs.push_back(info.ghostCellBoundaryID);  // not used in collocated
    pGC[0].push_back(0.0); // can be optimized
    pGC[1].push_back(0.0); // can be optimized
    pGC[2].push_back(0.0); // can be optimized
    pGCFull.push_back(0.0);// can be optimized
}

//##############################################################################
// Function Push_Back_GhostCellInfo
//  @param cell cell index
//  @param neighbour neighbour index
//  @param topology topology of the neighbour to this cell, returned value by
//                  ScalarField::cellNeighbours(int, IntArray, IntArray)
//  @param gctype 0: ghost cell for objects, 1: for constraint
//##############################################################################
std::unique_ptr<MAC_Grid::GhostCell> MAC_Grid::
MakeGhostCell(const int cell, const int neighbour, 
              const int topology, const int gctype)
{
    auto gc = std::make_unique<GhostCell>(); 
    gc->ownerCell = cell; 
    gc->neighbourCell = neighbour; 
    gc->topology = topology; 
    gc->type = gctype; 
    gc->pressure = 0.0; 
    gc->position = _pressureField.cellPosition(cell); 
    const int dim = abs(topology)-1; 
    gc->position[dim] += (topology > 0 ? 1.0 : -1.0)
                      *  0.25*_waveSolverSettings->cellSize; 
    auto key = gc->MakeKey(); 
    auto it = _ghostCellsCached.find(key); 
    if (it != _ghostCellsCached.end())
        gc->cache = std::move(it->second); 
    else 
        gc->cache = std::make_unique<GhostCell_Cache>(); 
    return gc; 
}


//##############################################################################
// This function returns flat index of the pressure cell that contains the 
// queried position. 
//##############################################################################
int MAC_Grid::InPressureCell(const Vector3d &position)
{
    const Vector3d &pBBoxLower = _pressureField.bbox().minBound();
    const REAL &h = _waveSolverSettings->cellSize; 
    const Tuple3i &divs = _pressureField.cellDivisions(); 
    const Tuple3i cellIndices = Tuple3i(std::max<int>(std::min<int>((int)((position.x - pBBoxLower.x)/h), divs[0]-1), 0),
                                        std::max<int>(std::min<int>((int)((position.y - pBBoxLower.y)/h), divs[1]-1), 0),
                                        std::max<int>(std::min<int>((int)((position.z - pBBoxLower.z)/h), divs[2]-1), 0)); 
    return _pressureField.cellIndex(cellIndices); 
}

//##############################################################################
// Function UpdatePMLAbsorptionCoeffs
//   This function updates all the pml cells coefficients when the scene center
//   changes.
//##############################################################################
void MAC_Grid::UpdatePMLAbsorptionCoeffs(const BoundingBox &sceneBox) 
{
    std::cout << "Update PML absorption coeffs using box: " 
              << sceneBox.minBound() << "; " 
              << sceneBox.maxBound() << std::endl;
    // constants alias
    const REAL &dt = _waveSolverSettings->timeStepSize; 
    const REAL &dx = _waveSolverSettings->cellSize; 
    const REAL &c = _waveSolverSettings->soundSpeed; 
    const REAL &rho = _waveSolverSettings->airDensity; 
    // update pressure PML cells coeff
    for (auto &p_cell : _pmlPressureCells)
    {
        REAL maxAbsorptionCoefficient = std::numeric_limits<REAL>::min(); 
        for (int dim=0; dim<3; ++dim)
        {
            const REAL absorptionCoefficient = PML_absorptionCoefficient(p_cell.position, _PML_absorptionWidth, dim, &sceneBox); 
            maxAbsorptionCoefficient = std::max<REAL>(maxAbsorptionCoefficient, absorptionCoefficient); 
            const REAL directionalCoefficient = PML_directionalCoefficient(absorptionCoefficient, dt);
            p_cell.updateCoefficient[dim] = PML_pressureUpdateCoefficient(absorptionCoefficient, dt, directionalCoefficient);
            p_cell.divergenceCoefficient[dim] = PML_divergenceCoefficient(rho, c, directionalCoefficient); 
        }
        p_cell.absorptionCoefficient = maxAbsorptionCoefficient; 
    }
    // update velocity PML cells coeff
    for (auto &v_cell : _pmlVelocityCells)
    {
        const REAL absorptionCoefficient = PML_absorptionCoefficient(v_cell.position, _PML_absorptionWidth, v_cell.dimension); 
        const REAL updateCoefficient = PML_velocityUpdateCoefficient(absorptionCoefficient, dt);
        const REAL gradientCoefficient  = PML_pressureGradientCoefficient(absorptionCoefficient, dt, dx, rho);
        v_cell.updateCoefficient = updateCoefficient; 
        v_cell.gradientCoefficient = gradientCoefficient; 
    }
}

//##############################################################################
// Function RemoveOldPML
//##############################################################################
void MAC_Grid::RemoveOldPML(const BoundingBox &sceneBox)
{
    std::cout << "Remove Old PML\n"; 
    // flag all cells need removal
    int flag; 
    int removeP=0, removeV=0;
    for (auto &p_cell : _pmlPressureCells) 
    {
        p_cell.needRemove = 
            (!InsidePML(p_cell.position, _PML_absorptionWidth, flag, &sceneBox)); 
        _isPMLCell.at(p_cell.index)  = (!p_cell.needRemove); 
        _isBulkCell.at(p_cell.index) = ( p_cell.needRemove); 
        if (p_cell.needRemove) ++removeP;
    }
    for (auto &v_cell : _pmlVelocityCells)
    {
        v_cell.needRemove = 
            (!InsidePML(v_cell.position, _PML_absorptionWidth, flag, &sceneBox)); 
        if (v_cell.needRemove) ++removeV;
    }
    _pmlPressureCells.erase(std::remove_if(_pmlPressureCells.begin(), 
                                           _pmlPressureCells.end(), 
                                           [](const auto &a){return a.needRemove;}),
                            _pmlPressureCells.end()); 
    _pmlVelocityCells.erase(std::remove_if(_pmlVelocityCells.begin(), 
                                           _pmlVelocityCells.end(), 
                                           [](const auto &a){return a.needRemove;}),
                            _pmlVelocityCells.end()); 
    std::cout << " remove pressure: " << removeP << "; now = " << _pmlPressureCells.size() << std::endl; 
    std::cout << " remove velocity: " << removeV << "; now = " << _pmlVelocityCells.size() << std::endl; 
}

//##############################################################################
// Function UpdatePML
//   This function updates all the pml cells 
//##############################################################################
void MAC_Grid::UpdatePML(const BoundingBox &sceneBox) 
{
    // flag all the cells that need removal
    RemoveOldPML(sceneBox); 
    UpdatePMLAbsorptionCoeffs(sceneBox); 
    //TODO 
}

void MAC_Grid::FillVandermondeRegularS(const int &row, const Vector3d &cellPosition, Eigen::MatrixXd &V,
                                       const Vector3d &origin, const REAL &h)
{
    Vector3d np = cellPosition; 
    np -= origin; 
    np /= (h); 
    np += 0.5;
    FillVandermondeRegular(row, np, V); 
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

void MAC_Grid::FillVandermondeBoundaryS(const int &row, const Vector3d &boundaryPosition, const Vector3d &boundaryNormal, Eigen::MatrixXd &V, const Vector3d &origin, const REAL &h)
{
    Vector3d np = boundaryPosition; 
    np -= origin; 
    np /= (h); 
    np += 0.5; 

    FillVandermondeBoundary(row, np, boundaryNormal, V); 
}

void MAC_Grid::ClearUnusedCache()
{
    typedef std::unordered_map<GhostCellId, std::unordered_map<int, TriangleIdentifier>>::iterator It_Outer; 
    typedef std::unordered_map<int, TriangleIdentifier>::iterator It_Inner; 
    It_Outer it_o = _ghostCellPreviousTriangles.begin(); 
    while (it_o != _ghostCellPreviousTriangles.end())
    {
        It_Inner it_i = it_o->second.begin(); 
        while (it_i != it_o->second.end())
        {
            // clear unused
            if (!it_i->second.active) it_i = it_o->second.erase(it_i); 
            else                      ++it_i; 
        }
        if (it_o->second.size()==0) 
        {
            it_o = _ghostCellPreviousTriangles.erase(it_o); 
        } 
        else 
        {
            ++it_o; 
        }
    }
}

void MAC_Grid::FillBoundaryFreshCellGrid(const int &dim, const int &ind, MATRIX &pCurr, const MATRIX &pLast)
{
    if (!(ind==0 || ind==(_pressureField.cellDivisions()[dim]-1)))
        throw std::runtime_error("**ERROR** Wrong MAC_GRid::FillBoundaryFreshCellGrid input: ind");
    const auto &field = _pressureField; 
    const Tuple3i &div = field.cellDivisions(); 
    const int &d0 =  dim; 
    const int  d1 = (dim+1)%3; 
    const int  d2 = (dim+2)%3; 
    const int interior_dir = (ind==0 ? 1 : -1); 
    const double &h = _waveSolverSettings->cellSize; 
    const double &k = _waveSolverSettings->timeStepSize; 
    const double &c = _waveSolverSettings->soundSpeed; 
    const double lambda = k*c/h; 
    const double lambda2 = pow(lambda,2); 
    Tuple3i indices, indices_int; 
    // fill data in i0
    for (int ii=0; ii<div[d1]; ++ii)
    for (int jj=0; jj<div[d2]; ++jj)
    {
        indices[d0] = ind; 
        indices[d1] = ii;
        indices[d2] = jj; 
        const int cell_idx_off0 = field.cellIndex(indices); 
        indices[d0] = ind + interior_dir;
        indices[d1] = ii; 
        indices[d2] = jj; 
        const int cell_idx_off1 = field.cellIndex(indices);
        indices[d0] = ind + interior_dir*2;
        indices[d1] = ii; 
        indices[d2] = jj; 
        const int cell_idx_off2 = field.cellIndex(indices);
        const int  ii_p = std::min(ii+1, div[d1]-1); 
        const int  ii_n = std::max(ii-1, 0      ); 
        const int  jj_p = std::min(jj+1, div[d2]-1); 
        const int  jj_n = std::max(jj-1, 0      ); 
        const int &kk   = ind; 
        const int  kk_p = kk+interior_dir; 
#define PCELL_IDX(_i,_j,_k) field.cellIndex(d1,d2,d0,_i,_j,_k)
#if ENGQUIST_ORDER == 1
        const double pNext = lambda2/(1.+ lambda)*(
                (2./lambda2 - 6.)*pCurr(PCELL_IDX(ii,jj,kk),0)
                + pCurr(PCELL_IDX(ii  ,jj_p,kk  ),0) + pCurr(PCELL_IDX(ii  ,jj_n,kk  ),0)
                + pCurr(PCELL_IDX(ii_p,jj  ,kk  ),0) + pCurr(PCELL_IDX(ii_n,jj  ,kk  ),0)
                +2.*pCurr(PCELL_IDX(ii,jj,kk_p),0) +(lambda - 1.0)/lambda2*pLast(PCELL_IDX(ii,jj,kk),0));
        pCurr(cell_idx_off0,0) = pCurr(cell_idx_off2,0) + 1.0/lambda*(pLast(cell_idx_off1,0)-pNext); 
#else
        //pCurr(cell_idx_off0,0) = 0.0;
        pCurr(cell_idx_off0,0) = pCurr(cell_idx_off1,0); 
#endif
#undef PCELL_IDX
    }
}

void MAC_Grid::GetAllBoundaryCells(const int &dimension, const int &sign, std::vector<int> &indices, std::vector<Vector3d> &positions)
{
    const Tuple3i divs = pressureFieldDivisions(); 
    const int &d = dimension; // type less..
    const int dx = (d+1)%3; 
    const int dy = (d+2)%3; 
    const int kk = (sign>0 ? divs[d]-1 : 0); 
    const int Nx = divs[dx];
    const int Ny = divs[dy];
    const int NxNy = Nx*Ny; 
    Tuple3i indicesBuf; 
    indicesBuf[d] = kk; 
    if (indices.size() != NxNy)   indices.resize(NxNy); 
    if (positions.size() != NxNy) positions.resize(NxNy); 
    for (int ii=0; ii<Nx; ++ii)
        for (int jj=0; jj<Ny; ++jj)
        {
            const int ind = ii*Nx + jj; 
            indicesBuf[dx] = ii; 
            indicesBuf[dy] = jj;
            indices.at(ind) = pressureFieldVertexIndex(indicesBuf); 
            positions.at(ind) = pressureFieldPosition(indicesBuf); 
        }
}

bool MAC_Grid::BoundaryGhostCellPressure(const std::string &solver_a, const int &cell_a, const int &cell_b, REAL &pressure_b) const 
{
    assert(grid_id); 
    const auto key = GhostCell::MakeKey_Boundary(solver_a, *grid_id, cell_a, cell_b); 
    auto it = _boundaryGhostCells.find(key); 
    if (it == _boundaryGhostCells.end())
        return false; 
    const auto &gc = it->second; 
    assert(gc->boundary); 
    pressure_b = gc->pressure; 

    return true; 
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
    const Tuple3i &pFieldDivs = _pressureField.cellDivisions();
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
            cellIndices[dim] = std::min<int>(pFieldDivs[dim]-1, cellIndices[dim]);
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
    else if (cell.s_identity.compare("Pressure PML Cell")==0)
    {
        os << " absorptionCoefficient: " << cell.r_identity/(-0.5) << "\n"; 
    }
    os << "--------------------------------------------------------------------------------" 
       << std::flush; 
    return os; 
}
