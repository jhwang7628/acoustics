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
#include <boost/timer/timer.hpp>

#include <distancefield/trilinearInterpolation.h> 
#include <linearalgebra/SparseLinearSystemSolver.h>
//#include <graph/UndirectedGraph.h> 

#include <unistd.h> 

#ifdef USE_OPENMP
#include <omp.h>
#endif

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

MAC_Grid::MAC_Grid(const BoundingBox &bbox, const PML_WaveSolver_Settings &settings, std::shared_ptr<FDTD_Objects> objects)
    : _pressureField(bbox,settings.cellSize), 
      _N(1), // no longer supports multiple fields, 
      _cornellBoxBoundaryCondition(settings.cornellBoxBoundaryCondition), 
      _useGhostCellBoundary(settings.useGhostCell),
      _objects(objects)
{
    setPMLBoundaryWidth(settings.PML_width, settings.PML_strength);
    Reinitialize_MAC_Grid(bbox,settings.cellSize); 
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
    _toggledBulkCells.resize(N_pressureCells, 0); 
    _toggledGhostCells.resize(N_pressureCells, 0); 
    _pressureCellHasValidHistory.resize(N_pressureCells, true); 
    _containingObject.resize(N_pressureCells, -1);
    if (_useGhostCellBoundary)
    {
        _isGhostCell.resize(N_pressureCells, false);
    }
    else 
    {
        for (int ii=0; ii<3; ++ii) 
        {
            const int N_velocityCells = _velocityField[ii].numCells(); 
            _isVelocityInterfacialCell[ii].resize(N_velocityCells, false); 
            _isVelocityBulkCell[ii].resize(N_velocityCells, false); 
            _toggledVelocityInterfacialCells[ii].resize(N_velocityCells, 0); 
            _velocityCellHasValidHistory[ii].resize(N_velocityCells, true); 
        }
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

void MAC_Grid::velocityDerivative( const MATRIX &p,
                                   const BoundaryEvaluator &bc,
                                   MATRIX &dV_dt, int dimension,
                                   REAL t, REAL alpha ) const
{
    const IntArray        &bulkCells = _velocityBulkCells[ dimension ];
    const IntArray        &interfacialCells
        = _velocityInterfacialCells[ dimension ];
    const ScalarField     &field = _velocityField[ dimension ];
    const IntArray        &interfaceBoundaryIDs
        = _interfacialBoundaryIDs[ dimension ];
    const FloatArray      &interfaceBoundaryDirections
        = _interfacialBoundaryDirections[ dimension ];

    // Handle all bulk cells
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for ( size_t bulk_cell_idx = 0; bulk_cell_idx < bulkCells.size();
            bulk_cell_idx++ )
    {
        int                  cell_idx = bulkCells[ bulk_cell_idx ];
        int                  neighbour_idx;
        Tuple3i              cell_indices = field.cellIndex( cell_idx );

        // If we're on a boundary, then set the derivative to zero (since
        // we should have d_n p = 0 on the boundary
        if ( cell_indices[ dimension ] == 0
                || cell_indices[ dimension ] == field.cellDivisions()[ dimension ] - 1 )
        {
            for ( int i = 0; i < _N; i++ )
            {
                dV_dt( cell_idx, i ) = 0.0;
            }
        }

        // Sample the derivative of the pressure field in the necessary
        // direction
        neighbour_idx = _pressureField.cellIndex( cell_indices );

        for ( int i = 0; i < _N; i++ )
        {
            dV_dt( cell_idx, i )
                += alpha * p( neighbour_idx, i ) / _pressureField.cellSize();
        }

        cell_indices[ dimension ] -= 1;
        neighbour_idx = _pressureField.cellIndex( cell_indices );

        for ( int i = 0; i < _N; i++ )
        {
            dV_dt( cell_idx, i )
                -= alpha * p( neighbour_idx, i ) / _pressureField.cellSize();
        }
    }

    // Handle boundary cells
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for ( size_t interfacial_cell_idx = 0;
            interfacial_cell_idx < interfacialCells.size();
            interfacial_cell_idx++ )
    {
        int                  cell_idx = interfacialCells[ interfacial_cell_idx ];
        int                  objectID;
        REAL                 bcEval;

        Vector3d             normal( 0.0, 0.0, 0.0 );
        Vector3d             x = field.cellPosition( cell_idx );

        normal[ dimension ] = interfaceBoundaryDirections[ interfacial_cell_idx ];

        objectID = interfaceBoundaryIDs[ interfacial_cell_idx ];

        for ( int i = 0; i < _N; i++ )
        {
            bcEval = bc( x, normal, objectID, t, i );

            dV_dt( cell_idx, i ) += alpha * bcEval;
        }
    }
}

void MAC_Grid::pressureDerivative( const MATRIX *v[ 3 ], MATRIX &p,
                                   REAL alpha ) const
{
    // We only have to worry about bulk cells here
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for ( size_t bulk_cell_idx = 0; bulk_cell_idx < _bulkCells.size();
            bulk_cell_idx++ )
    {
        int                      cell_idx = _bulkCells[ bulk_cell_idx ];
        int                      neighbour_idx;
        Tuple3i                  cell_indices;

        cell_indices = _pressureField.cellIndex( cell_idx );

        // Add the divergence term in each direction
        for ( int dimension = 0; dimension < 3; dimension++ )
        {
            const MATRIX          &vd = *( v[ dimension ] );

            cell_indices[ dimension ] += 1;
            neighbour_idx = _velocityField[ dimension ].cellIndex( cell_indices );

            for ( int i = 0; i < _N; i++ )
            {
                REAL addedValue = alpha * vd( neighbour_idx, i );
                addedValue /= _velocityField[ dimension ].cellSize();

                p( cell_idx, i ) += alpha * vd( neighbour_idx, i )
                    / _velocityField[ dimension ].cellSize();
            }

            cell_indices[ dimension ] -= 1;
            neighbour_idx = _velocityField[ dimension ].cellIndex( cell_indices );

            for ( int i = 0; i < _N; i++ )
            {
                REAL addedValue = -1.0 * alpha * vd( neighbour_idx, i );
                addedValue /= _velocityField[ dimension ].cellSize();

                p( cell_idx, i ) -= alpha * vd( neighbour_idx, i )
                    / _velocityField[ dimension ].cellSize();
            }
        }
    }
}

void MAC_Grid::PML_velocityUpdate( const MATRIX &p, MATRIX &v, int dimension, REAL t, REAL timeStep, REAL density )
{
    const IntArray        &bulkCells                    = _velocityBulkCells[ dimension ];
    const IntArray        &interfacialCells             = _velocityInterfacialCells[ dimension ];
    const ScalarField     &field                        = _velocityField[ dimension ];
    const FloatArray      &interfaceBoundaryCoefficients= _interfacialBoundaryCoefficients[ dimension ];

    // Handle all bulk cells
    PML_velocityUpdateAux(p, v, dimension, t, timeStep, density, bulkCells); 

    // Handle interfacial cells
    if (_useGhostCellBoundary)
    {
        PML_velocityUpdateAux(p, v, dimension, t, timeStep, density, interfacialCells); 
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
            //int                  objectID;
            REAL                 bcEval;
            REAL                 coefficient;

            Vector3d             normal( 0.0, 0.0, 0.0 ); // normal of the interfacial cell boundary
            Vector3d             x = field.cellPosition( cell_idx ); // position of the interfacial cell

            normal[ dimension ] = 1.0;
            coefficient = interfaceBoundaryCoefficients[ interfacial_cell_idx ];
            //objectID = interfaceBoundaryIDs[ interfacial_cell_idx ];

            for ( int i = 0; i < _N; i++ )
            {
                bcEval = _objects->EvaluateNearestVibrationalSources(x,normal,t); 
                //bcEval = bc( x, normal, objectID, t, i );
                bcEval *= coefficient;
                v( cell_idx, i ) += timeStep * bcEval;
            }
        }
    }
}

void MAC_Grid::PML_velocityUpdateAux( const MATRIX &p, MATRIX &v, int dimension, REAL t, REAL timeStep, REAL density, const IntArray &bulkCells )
{
    const ScalarField &field = _velocityField[ dimension ];
    const REAL n_dt_over_dx_rho = -timeStep/(_pressureField.cellSize()*density);
    const int bulkCellSize = bulkCells.size(); 

    // Handle all bulk cells
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(static) default(shared)
    #endif
    for (int bulk_cell_idx = 0; bulk_cell_idx < bulkCellSize; ++bulk_cell_idx)
    {
        const int       cell_idx                = bulkCells[bulk_cell_idx];
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
}

void MAC_Grid::PML_pressureUpdate( const MATRIX &v, MATRIX &p, int dimension, REAL timeStep, REAL c, const ExternalSourceEvaluator *sourceEvaluator, const REAL simulationTime, REAL density )
{
    const size_t bulkCellSize = _bulkCells.size();
    //const bool evaluateExternalSource = (sourceEvaluator != nullptr);
    const bool evaluateExternalSource = _objects->HasExternalPressureSources();
    const REAL n_rho_c_square_dt = -density*c*c*timeStep;
    const REAL v_cellSize_inv = 1./_velocityField[dimension].cellSize();

    // We only have to worry about bulk cells here
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(static) default(shared)
    #endif
    for (size_t bulk_cell_idx = 0; bulk_cell_idx < bulkCellSize; ++bulk_cell_idx)
    {
        const int       cell_idx                = _bulkCells[ bulk_cell_idx ];
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
                p(cell_idx, i) += n_rho_c_square_dt*v(neighbour_idx, i)*v_cellSize_inv;

            cell_indices[dimension] -= 1;
            neighbour_idx = _velocityField[dimension].cellIndex(cell_indices);

            for (int i=0; i<_N; i++)
                p(cell_idx, i) -= n_rho_c_square_dt*v(neighbour_idx, i)*v_cellSize_inv;


            // evaluate external sources only happens not in PML
            // Liu Eq (16) f6x term
            if (evaluateExternalSource)
                p(cell_idx,0) += _objects->EvaluatePressureSources(cell_position, cell_position, simulationTime+0.5*timeStep)*timeStep;
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

// TODO can optimize the sparse linear system setup
//
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
void MAC_Grid::PML_pressureUpdateGhostCells( MATRIX &p, const REAL &timeStep, const REAL &c, const REAL &simulationTime, const REAL density)
{
//    // for the ghost-cell coupling
//    SparseLinearSystemSolver solver(_ghostCells.size()); 
//
//#ifdef USE_OPENMP
//#pragma omp parallel for schedule(static) default(shared)
//#endif
//    for (size_t ghost_cell_idx=0; ghost_cell_idx<_ghostCells.size(); ghost_cell_idx++) 
//    {
//
//        const int       cellIndex      = _ghostCells[ghost_cell_idx]; 
//        const Tuple3i   cellIndices    = _pressureField.cellIndex(cellIndex); 
//        const Vector3d  cellPosition   = _pressureField.cellPosition(cellIndices); 
//        const int       boundaryObject = _containingObject[cellIndex]; 
//
//        // find point BI and IP in the formulation
//        Vector3d boundaryPoint, imagePoint, erectedNormal; 
//        FindImagePoint(cellPosition, boundaryObject, boundaryPoint, imagePoint, erectedNormal); 
//
//        // get the box enclosing the image point; 
//        IntArray neighbours; 
//        _pressureField.enclosingNeighbours(imagePoint, neighbours); 
//
//        assert(neighbours.size()==8);
//
//        // hasGC  : has self as interpolation stencil
//        int hasGC=-1; 
//
//        std::vector<int> coupledGhostCells; 
//        std::vector<int> uncoupledGhostCellsNeighbours; // [0,neighbours.size)
//        std::vector<int> coupledGhostCellsNeighbours; // [0,neighbours.size)
//
//        for (size_t ii=0; ii<neighbours.size(); ii++) 
//        {
//            if (neighbours[ii] == cellIndex) 
//            {
//                hasGC = ii;
//            }
//            else  
//            {
//                // if not itself but a ghost cell then add an edge to graph
//                // adjacency matrix
//                if (_isGhostCell[neighbours[ii]]) 
//                {
//                    coupledGhostCells.push_back(_ghostCellsInverse.at(neighbours[ii])); 
//                    coupledGhostCellsNeighbours.push_back(ii); 
//                }
//                else
//                    uncoupledGhostCellsNeighbours.push_back(ii);
//            }
//        }
//
//        // vandermonde matrix, see 2008 Mittals JCP paper Eq.18
//        Eigen::MatrixXd V(8,8); 
//        Tuple3i     indicesBuffer;
//        Vector3d    positionBuffer; 
//        Eigen::VectorXd pressureNeighbours(8);  // right hand side
//
//        // evaluate at boundarPoint, the boundary is prescribing normal acceleration so scaling is needed for pressure Neumann
//        const double bcEval = _objects->EvaluateVibrationalSources(boundaryPoint, erectedNormal, simulationTime)*(-density);
//        //const double bcEval = bc(boundaryPoint, erectedNormal, boundaryObject, simulationTime, 0) * (-density);  
//
//        for (size_t row=0; row<neighbours.size(); row++) 
//        {
//            indicesBuffer = _pressureField.cellIndex(neighbours[row]); 
//            positionBuffer= _pressureField.cellPosition(indicesBuffer); 
//            
//            if ((int)row!=hasGC)
//            {
//                FillVandermondeRegular(row,positionBuffer, V);
//                pressureNeighbours(row) = p(neighbours[row],0); 
//            }
//            else 
//            {
//                FillVandermondeBoundary(row, boundaryPoint, erectedNormal, V);
//                pressureNeighbours(row) = bcEval; 
//            }
//        }
//
//        // want beta = V^-T b, equivalent to solving V^T x = b, and then x^T
//        Eigen::MatrixXd b(1,8); // row vector
//
//        FillVandermondeRegular(0,imagePoint,b); 
//        b.transposeInPlace();  // column vector
//        V.transposeInPlace(); 
//
//        // TODO svd solve is needed? can I use QR? 
//        Eigen::JacobiSVD<Eigen::MatrixXd> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
//        Eigen::VectorXd beta = svd.solve(b);
//        //Eigen::VectorXd beta = V.householderQr().solve(b); 
//        //const double conditionNumber = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
//        //if (conditionNumber > 1E5 && replaced) 
//        //{
//        //    std::cout << "**WARNING** condition number for the least square solve is = " << conditionNumber << "\n"
//        //              << "            the solution can be inaccurate.\n"
//        //              << "largest  singular value = " << svd.singularValues()(0) << "\n"
//        //              << "smallest singular value = " << svd.singularValues()(svd.singularValues().size()-1) << "\n"
//        //              << "problematic matrix: \n"
//        //              << V << std::endl;
//        //}
//        //
//#pragma omp critical
//        {
//
//            // fill in the sparse matrix
//            
//            // start by filling the diagonal entry
//            solver.StageEntry(ghost_cell_idx, ghost_cell_idx, 1.0);
//
//            // next fill the off-diagonal terms by iterating through the coupled
//            // ghost cells for this row. its coefficient is defined by the i-th
//            // entry of beta. i=1,...,8 are the index of the neighbour
//            for (size_t cc=0; cc<coupledGhostCells.size(); cc++) 
//                solver.StageEntry(ghost_cell_idx, coupledGhostCells[cc], -beta(coupledGhostCellsNeighbours[cc])); 
//
//
//            // FIXME debug
//            //double offdiagonalSum=0.0; 
//            //for (size_t cc=0; cc<coupledGhostCells.size(); cc++) 
//            //    offdiagonalSum += fabs(beta(coupledGhostCellsNeighbours[cc])); 
//
//            //if (offdiagonalSum>1.0)
//            //    std::cerr << "**WARNING** not diagonally dominant\n"; 
//
//            // next we compute the right-hand-side (RHS) of the equation. looking
//            // at the equation, it will be the neumann condition + any uncoupled
//            // betas (see google doc for eqution)
//            double RHS = - (imagePoint - cellPosition).length() * bcEval; 
//            for (size_t uc=0; uc<uncoupledGhostCellsNeighbours.size(); uc++) 
//                RHS += beta(uncoupledGhostCellsNeighbours[uc])*pressureNeighbours(uncoupledGhostCellsNeighbours[uc]); 
//
//            // if GC itself is coupled, need to add this to the RHS
//            if (hasGC != -1) 
//            {
//                RHS += beta(hasGC)*pressureNeighbours(hasGC); 
//            }
//
//            // finally we fill the rhs of the sparse linear system
//            solver.FillInRHS(ghost_cell_idx, RHS); 
//        }
//
//        // actually fill-in all the staged entry of the sparse matrix. after this
//        // step the matrix is ready for solve. 
//        {
//            boost::timer::auto_cpu_timer t(" sparse system fill-in takes %w sec\n"); 
//            solver.FillIn(); 
//        }
//
//        // the actual solve. unfortunately eigen supports only dense RHS solve. 
//        Eigen::VectorXd x; 
//        {
//            boost::timer::auto_cpu_timer t(" linear system solve takes %w sec\n");
//            solver.Solve(x); 
//            //solver.DenseSolve(x); 
//        }
//        // put the solution into matrix p
//        {
//            boost::timer::auto_cpu_timer t(" recover solution to dense vector takes %w sec\n");
//            for (size_t ghost_cell_idx=0; ghost_cell_idx<_ghostCells.size(); ghost_cell_idx++) 
//                p(_ghostCells[ghost_cell_idx],0) = x(ghost_cell_idx);
//        }
//    }
}

void MAC_Grid::PML_pressureUpdateGhostCells_Jacobi( MATRIX &p, const REAL &timeStep, const REAL &c, const REAL &simulationTime, const REAL density)
{

    _ghostCellCoupledData.clear(); 
    _ghostCellCoupledData.resize(_ghostCells.size()); 

    //debug FIXME 
    std::ofstream of, of2, of3; 
    bool write = false; 

    if ((int)(simulationTime/timeStep)%4==0)
    {
        write = true; 
        of.open("work/near-field/data_test/boundaryPoint_"+std::to_string((int)(simulationTime/timeStep/4.0))+".csv"); 
        of2.open("work/near-field/data_test/ghostPoint_"+std::to_string((int)(simulationTime/timeStep/4.0))+".csv"); 
    }
      
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (size_t ghost_cell_idx=0; ghost_cell_idx<_ghostCells.size(); ghost_cell_idx++) 
    {

        const int       cellIndex      = _ghostCells[ghost_cell_idx]; 
        const Tuple3i   cellIndices    = _pressureField.cellIndex(cellIndex); 
        const Vector3d  cellPosition   = _pressureField.cellPosition(cellIndices); 
        //const int       boundaryObject = _containingObject[cellIndex]; 

        // find point BI and IP in the formulation
        Vector3d boundaryPoint, imagePoint, erectedNormal; 
        REAL accumulatedBoundaryConditionValue; 
        //FindImagePoint(cellPosition, boundaryObject, boundaryPoint, imagePoint, erectedNormal); 
        _objects->ReflectAgainstAllBoundaries(cellPosition, simulationTime, imagePoint, boundaryPoint, erectedNormal, accumulatedBoundaryConditionValue, density, 1);

        // debug FIXME
#pragma omp critical
        {
        if (write)
            of << boundaryPoint.x << "," << boundaryPoint.y << "," << boundaryPoint.z << std::endl;
            of2 << cellPosition.x << "," << cellPosition.y << "," << cellPosition.z << std::endl;
        }

        // get the box enclosing the image point; 
        IntArray neighbours; 
        _pressureField.enclosingNeighbours(imagePoint, neighbours); 

        // hasGC  : has self as interpolation stencil
        // coupled: if at least one interpolation stencil is another ghost-cell 
        int hasGC=-1; 
        //int coupled=0; 

        std::vector<int> coupledGhostCells; 
        std::vector<int> uncoupledGhostCellsNeighbours; // [0,neighbours.size)
        std::vector<int> coupledGhostCellsNeighbours; // [0,neighbours.size)

        for (size_t ii=0; ii<neighbours.size(); ii++) 
        {
            if (neighbours[ii] == cellIndex) 
            {
                hasGC = ii;
            }
            else  
            {
                // if not itself but a ghost cell then add an edge to graph
                // adjacency matrix
                if (_isGhostCell[neighbours[ii]]) 
                {
                    coupledGhostCells.push_back(_ghostCellsInverse.at(neighbours[ii])); 
                    coupledGhostCellsNeighbours.push_back(ii); 
                }
                else
                    uncoupledGhostCellsNeighbours.push_back(ii);
            }
        }

        // vandermonde matrix, see 2008 Mittals JCP paper Eq.18
        Eigen::MatrixXd V(8,8); 
        Tuple3i     indicesBuffer;
        Vector3d    positionBuffer; 
        Eigen::VectorXd pressureNeighbours(8);  // right hand side
 
        // evaluate at boundarPoint, the boundary is prescribing normal acceleration so scaling is needed for pressure Neumann
        const REAL bcEval = _objects->EvaluateNearestVibrationalSources(boundaryPoint, erectedNormal, simulationTime)*(-density);
        //const double bcEval = bc(boundaryPoint, erectedNormal, boundaryObject, simulationTime, 0) * (-density);  

        for (size_t row=0; row<neighbours.size(); row++) 
        {
            indicesBuffer = _pressureField.cellIndex(neighbours[row]); 
            positionBuffer= _pressureField.cellPosition(indicesBuffer); 
            
            if ((int)row!=hasGC)
            {
                FillVandermondeRegular(row,positionBuffer, V);
                pressureNeighbours(row) = p(neighbours[row],0); 
            }
            else 
            {
                FillVandermondeBoundary(row, boundaryPoint, erectedNormal, V);
                pressureNeighbours(row) = bcEval; 
            }
        }

        // want beta = V^-T b, equivalent to solving V^T beta = b
        Eigen::MatrixXd b(1,8); // row vector

        FillVandermondeRegular(0,imagePoint,b); 
        b.transposeInPlace();  // column vector
        V.transposeInPlace(); 

        // TODO svd solve is needed? can I use QR? 
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::VectorXd beta = svd.solve(b); 
        //Eigen::VectorXd beta = V.householderQr().solve(b); 
        //const double conditionNumber = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
        //if (conditionNumber > 1E5 && replaced) 
        //{
        //    std::cout << "**WARNING** condition number for the least square solve is = " << conditionNumber << "\n"
        //              << "            the solution can be inaccurate.\n"
        //              << "largest  singular value = " << svd.singularValues()(0) << "\n"
        //              << "smallest singular value = " << svd.singularValues()(svd.singularValues().size()-1) << "\n"
        //              << "problematic matrix: \n"
        //              << V << std::endl;
        //}
        
        JacobiIterationData &jacobiIterationData = _ghostCellCoupledData[ghost_cell_idx]; 
        for (size_t cc=0; cc<coupledGhostCells.size(); cc++) 
        {
            if (fabs(beta(coupledGhostCellsNeighbours[cc])) > 1E-14)
            {
              jacobiIterationData.nnzIndex.push_back(coupledGhostCells[cc]); 
              jacobiIterationData.nnzValue.push_back(-beta(coupledGhostCellsNeighbours[cc])); 
            }
        } 

        // RHS always exists even when no off-diagonal elements 
        double &RHS = jacobiIterationData.RHS; 
        RHS = accumulatedBoundaryConditionValue;
        //RHS = - (imagePoint - cellPosition).length() * bcEval; 
        for (size_t uc=0; uc<uncoupledGhostCellsNeighbours.size(); uc++) 
            RHS += beta(uncoupledGhostCellsNeighbours[uc])*pressureNeighbours(uncoupledGhostCellsNeighbours[uc]); 

        // if GC itself is coupled, need to add this to the RHS
        if (hasGC != -1) 
            RHS += beta(hasGC)*pressureNeighbours(hasGC); 
    }

    // debug FIXME
    if (write)
    {
        of.close();
        of2.close();
    }


    const int maxIteration = 50;
    for (int iteration=0; iteration<maxIteration; iteration++) 
    {
        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(static) default(shared)
        #endif
        for (size_t ghost_cell_idx=0; ghost_cell_idx<_ghostCells.size(); ghost_cell_idx++) 
        {
            const JacobiIterationData &data = _ghostCellCoupledData[ghost_cell_idx]; 
            // if there are nnz, update them 
            p(_ghostCells[ghost_cell_idx],0) = data.RHS; 
            for (size_t cc=0; cc<data.nnzIndex.size(); cc++)
                p(_ghostCells[ghost_cell_idx],0) -= data.nnzValue[cc]*p(_ghostCells[data.nnzIndex[cc]],0); 
        }
    }
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

void MAC_Grid::FreshCellInterpolate(MATRIX &p, const REAL &simulationTime, const REAL &density)
{
//    const int N_cells = _toggledBulkCells.size(); 
//    const int N = _objects->N(); 
//    std::vector<ScalarField::RangeIndices> indices(N); 
//    for (int bbox_id = 0; bbox_id<N; ++bbox_id) 
//    {
//        const FDTD_MovableObject::BoundingBox &unionBox = _objects->Get(bbox_id).GetUnionBBox();
//        // enlarge the box by a little bit in order to include the fresh cells
//        const Vector3d maxBound = unionBox.maxBound + 2.0*_cellSize; 
//        const Vector3d minBound = unionBox.minBound - 2.0*_cellSize; 
//        _pressureField.GetIterationBox(minBound, maxBound, indices[bbox_id]); 
//    }
//
//#ifdef USE_OPENMP
//#pragma omp parallel for schedule(static) default(shared)
//#endif
//    for (int bbox_id=0; bbox_id<N; ++bbox_id)
//    {
//        Vector3d cellPosition;
//        int ii,jj,kk; 
//        FOR_ALL_3D_GRID_VECTOR3(indices[bbox_id].startIndex, indices[bbox_id].dimensionIteration, ii, jj, kk)
//        {
//            const Tuple3i cellIndices(ii,jj,kk);
//            const int cell_idx = _pressureField.cellIndex(cellIndices); 
//            cellPosition = _pressureField.cellPosition(cell_idx);
//
//            if (_toggledBulkCells[cell_idx] == 1) // this is the cell which turns into bulk cells
//            {
//                REAL distance; 
//                int objectID; 
//                _objects->LowestObjectDistance(cellPosition, distance, objectID); 
//                if (distance < DISTANCE_TOLERANCE)
//                    throw std::runtime_error("**ERROR** fresh cell inside some object. may be classification rounding error. distance: " + std::to_string(distance)); 
//                Vector3d imagePoint, boundaryPoint, erectedNormal; 
//                REAL distanceTravelled; 
//                _objects->Get(objectID).FindImageFreshCell(cellPosition, imagePoint, boundaryPoint, erectedNormal, distanceTravelled); 
//
//                // prepare interpolation stencils, this part is similar to
//                // the vandermonde part in ghost cell pressure update
//                IntArray neighbours; 
//                _pressureField.enclosingNeighbours(imagePoint, neighbours); 
//                Eigen::MatrixXd V(8,8); 
//                Tuple3i indicesBuffer;
//                Vector3d positionBuffer; 
//                Eigen::VectorXd RHS(8); 
//                const REAL bcEval = _objects->EvaluateVibrationalSources(boundaryPoint, erectedNormal, simulationTime)*(-density);
//
//                for (size_t row=0; row<neighbours.size(); row++) 
//                {
//                    // TODO 
//                    if (!_isBulkCell[neighbours[row]]) 
//                        throw std::runtime_error("**ERROR** one of the interpolation stencil for fresh cell is not bulk, this is not handled"); 
//                    indicesBuffer = _pressureField.cellIndex(neighbours[row]); 
//                    positionBuffer= _pressureField.cellPosition(indicesBuffer); 
//                    if (neighbours[row] != cell_idx) // not self
//                    {
//                        FillVandermondeRegular(row, positionBuffer, V);
//                        RHS(row) = p(neighbours[row], 0); 
//                    }
//                    else 
//                    {
//                        FillVandermondeBoundary(row, boundaryPoint, erectedNormal, V);
//                        RHS(row) = bcEval; 
//                    }
//                }
//
//                // coefficient for the interpolant
//                Eigen::JacobiSVD<Eigen::MatrixXd> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
//                //Eigen::VectorXd C = V.lu().solve(RHS); 
//                Eigen::VectorXd C = svd.solve(RHS); 
//
//                // evaluate the interpolant at cell position
//                Eigen::VectorXd coordinateVector(8); 
//                FillVandermondeRegular(cellPosition, coordinateVector);
//                const REAL pressureFreshCell = C.dot(coordinateVector); 
//                p(cell_idx, 0) = pressureFreshCell; 
//
//                // set the neighbouring new ghost cells into the same value as
//                // the fresh cell
//                for (size_t nei=0; nei<neighbours.size(); nei++) 
//                {
//                    std::cout << _toggledGhostCells[neighbours[nei]] << " "; 
//                    if (_toggledGhostCells[neighbours[nei]] == 1) 
//                    {
//                        exit(1);
//                    }
//                }
//                std::cout << std::endl;
//            }
//        } 
//    }

    //// FIXME DEBUG ////
      
    //const int N_cells = _toggledBulkCells.size(); 
    //const int N = _objects->N(); 
    //std::vector<ScalarField::RangeIndices> indices(N); 
    //for (int bbox_id = 0; bbox_id<N; ++bbox_id) 
    //{
    //    const FDTD_MovableObject::BoundingBox &unionBox = _objects->Get(bbox_id).GetUnionBBox();
    //    // enlarge the box by a little bit in order to include the fresh cells
    //    const Vector3d maxBound = unionBox.maxBound + 2.0*_cellSize; 
    //    const Vector3d minBound = unionBox.minBound - 2.0*_cellSize; 
    //    _pressureField.GetIterationBox(minBound, maxBound, indices[bbox_id]); 
    //}

    const int numPressureCells = _pressureField.numCells();
    for (int cell_idx = 0; cell_idx<numPressureCells; ++cell_idx) 
    {
        Vector3d cellPosition = _pressureField.cellPosition(cell_idx);
        if (_toggledBulkCells[cell_idx] == 1) // this is the fresh pressure cell
        {
            REAL distance; 
            int objectID; 
            _objects->LowestObjectDistance(cellPosition, distance, objectID); 
            if (distance < DISTANCE_TOLERANCE)
                throw std::runtime_error("**ERROR** fresh cell inside some object. may be classification rounding error. distance: " + std::to_string(distance)); 
            Vector3d imagePoint, boundaryPoint, erectedNormal; 
            REAL distanceTravelled; 
            _objects->Get(objectID).FindImageFreshCell(cellPosition, imagePoint, boundaryPoint, erectedNormal, distanceTravelled); 

            // prepare interpolation stencils, this part is similar to
            // the vandermonde part in ghost cell pressure update
            IntArray neighbours; 
            _pressureField.enclosingNeighbours(imagePoint, neighbours); 
            Eigen::MatrixXd V(8,8); 
            Tuple3i indicesBuffer;
            Vector3d positionBuffer; 
            Eigen::VectorXd RHS(8); 
            const REAL bcEval = _objects->Get(objectID).EvaluateBoundaryAcceleration(boundaryPoint, erectedNormal, simulationTime)*(-density);

            for (size_t row=0; row<neighbours.size(); row++) 
            {
                // TODO 
                if (!_isBulkCell[neighbours[row]]) 
                    throw std::runtime_error("**ERROR** one of the interpolation stencil for fresh cell is not bulk, this is not handled"); 
                indicesBuffer = _pressureField.cellIndex(neighbours[row]); 
                positionBuffer= _pressureField.cellPosition(indicesBuffer); 
                if (neighbours[row] != cell_idx) // not self
                {
                    FillVandermondeRegular(row, positionBuffer, V);
                    RHS(row) = p(neighbours[row], 0); 
                }
                else 
                {
                    FillVandermondeBoundary(row, boundaryPoint, erectedNormal, V);
                    RHS(row) = bcEval; 
                }
            }

            // coefficient for the interpolant
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
            //Eigen::VectorXd C = V.lu().solve(RHS); 
            Eigen::VectorXd C = svd.solve(RHS); 

            // evaluate the interpolant at cell position
            Eigen::VectorXd coordinateVector(8); 
            FillVandermondeRegular(cellPosition, coordinateVector);
            const REAL pressureFreshCell = C.dot(coordinateVector); 
            p(cell_idx, 0) = pressureFreshCell; 

            // reuse the neighbour int array to get neighbours of the current
            // fresh cell
            // set the neighbouring new ghost cells into the same value as
            // the fresh cell
            _pressureField.cellNeighbours(cell_idx, neighbours); 
            for (size_t nei=0; nei<neighbours.size(); nei++) 
            {
                std::cout << neighbours[nei] << " : " << _toggledGhostCells[neighbours[nei]] << " "; 
                if (_toggledGhostCells[neighbours[nei]] == 1) 
                {
                    p(neighbours[nei], 0) = p(cell_idx, 0); 
                }
            }
            std::cout << std::endl;
        }
    }
    //// DEBUG END FIXME ////
} 

void MAC_Grid::InterpolateFreshPressureCell_Rasterized(MATRIX &p, const REAL &timeStep, const REAL &simulationTime, const REAL &density)
{
    const int N = _pressureField.numCells(); 
    for (int cell_idx = 0; cell_idx < N; ++cell_idx)
    {
        // if has valid history or if its solid, then no interpolation needed
        if (_pressureCellHasValidHistory.at(cell_idx) || IsPressureCellSolid(cell_idx))
            continue; 

        const Vector3d cellPosition = _pressureField.cellPosition(cell_idx);
        REAL distance; 
        int objectID; 
        _objects->LowestObjectDistance(cellPosition, distance, objectID); 
        if (distance < DISTANCE_TOLERANCE)
            throw std::runtime_error("**ERROR** Fresh cell inside some object. This shouldn't happen for pressure cells. Distance: " + std::to_string(distance)); 
        Vector3d imagePoint, boundaryPoint, erectedNormal; 
        REAL distanceTravelled; 
        _objects->Get(objectID).FindImageFreshCell(cellPosition, imagePoint, boundaryPoint, erectedNormal, distanceTravelled); 

        // prepare interpolation stencils, this part is similar to
        // the vandermonde part in ghost cell pressure update
        IntArray neighbours; 
        _pressureField.enclosingNeighbours(imagePoint, neighbours); 
        Eigen::MatrixXd V(8,8); 
        Tuple3i indicesBuffer;
        Vector3d positionBuffer; 
        Eigen::VectorXd RHS(8); 
        const REAL boundaryPressureGradientNormal = _objects->Get(objectID).EvaluateBoundaryAcceleration(boundaryPoint, erectedNormal, simulationTime - 0.5*timeStep)*(-density); 

        for (size_t row=0; row<neighbours.size(); row++) 
        {
            const int neighbour_idx = neighbours[row]; 
            if (!_pressureCellHasValidHistory.at(neighbour_idx) && neighbour_idx != cell_idx)
                throw std::runtime_error("**ERROR** one of the interpolation stencil for pressure fresh cell has invalid history."); 

            indicesBuffer = _pressureField.cellIndex(neighbours[row]); 
            positionBuffer= _pressureField.cellPosition(indicesBuffer); 
            if (neighbours[row] != cell_idx) // not self
            {
                FillVandermondeRegular(row, positionBuffer, V);
                RHS(row) = p(neighbours[row], 0); 
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

        // evaluate the interpolant at cell position
        Eigen::VectorXd coordinateVector(8); 
        FillVandermondeRegular(cellPosition, coordinateVector);
        const REAL pressureFreshCell = C.dot(coordinateVector); 
        p(cell_idx, 0) = pressureFreshCell; 
    }
}

//##############################################################################
// Interpolate the fresh velocity cell
//
// Note: 
//   fresh velocity cell can be inside the object, so make sure the
//   interpolation procedure is robust
//##############################################################################
void MAC_Grid::InterpolateFreshVelocityCell_Rasterized(MATRIX &v, const int &dim, const REAL &timeStep, const REAL &simulationTime)
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
        REAL distanceToMesh;
        _objects->Get(objectID).FindImageFreshCell(cellPosition, imagePoint, boundaryPoint, erectedNormal, distanceToMesh); 

        // prepare interpolation stencils, this part is similar to
        // the vandermonde part in ghost cell pressure update
        IntArray neighbours; 
        Tuple3i indicesBuffer;
        Vector3d positionBuffer; 

        // interpolate using trilinear interpolant
        // boundary velocity will be scaled by the normal difference
        _velocityField[dim].enclosingNeighbours(imagePoint, neighbours); 
        const REAL boundaryVelocity = _objects->Get(objectID).EvaluateBoundaryVelocity(boundaryPoint, erectedNormal, simulationTime - timeStep) * _interfacialBoundaryCoefficients[dim].at(interfacial_cell_idx); 
        Eigen::MatrixXd V(8,8); 
        Eigen::VectorXd RHS(8); 
        int row_cell_idx = -1; 
        for (size_t row=0; row<neighbours.size(); row++) 
        {
            const int neighbour_idx = neighbours[row]; 

            if (!_velocityCellHasValidHistory[dim].at(neighbour_idx) && neighbour_idx != cell_idx) 
                throw std::runtime_error("**ERROR** one of the interpolation stencil for velocity fresh cell has invalid history."); 

            indicesBuffer = _velocityField[dim].cellIndex(neighbours.at(row)); 
            positionBuffer= _velocityField[dim].cellPosition(indicesBuffer); 
            if (neighbours.at(row) != cell_idx) // not self
            {
                FillVandermondeRegular(row, positionBuffer, V);
                RHS(row) = v(neighbours.at(row), 0); 
            }
            else 
            {
                row_cell_idx = row; 
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
        if (distanceToMesh > DISTANCE_TOLERANCE) // cell position is located outside the boundary: the erected normal is BI-CU-IP
            v(cell_idx, 0) = 0.5 * (v_IP + v_BI); 
        else // cell position is located inside the boundary: the erected normal is CU-BI-IP
            v(cell_idx, 0) = 2.0 * v_BI - v_IP;
    }
}

void MAC_Grid::InterpolateFreshPressureCell_GhostCell(MATRIX &p, const REAL &timeStep, const REAL &simulationTime, const REAL &density)
{
    const int N = _pressureField.numCells(); 
    for (int cell_idx = 0; cell_idx < N; ++cell_idx)
    {
        if (_toggledBulkCells[cell_idx] != 1) // no history problem, skip
            continue; 

        const Vector3d cellPosition = _pressureField.cellPosition(cell_idx);
        REAL distance; 
        int objectID; 
        _objects->LowestObjectDistance(cellPosition, distance, objectID); 
        if (distance < DISTANCE_TOLERANCE)
            throw std::runtime_error("**ERROR** fresh cell inside some object. may be classification rounding error. distance: " + std::to_string(distance)); 
        Vector3d imagePoint, boundaryPoint, erectedNormal; 
        REAL distanceTravelled; 
        _objects->Get(objectID).FindImageFreshCell(cellPosition, imagePoint, boundaryPoint, erectedNormal, distanceTravelled); 

        // prepare interpolation stencils, this part is similar to
        // the vandermonde part in ghost cell pressure update
        IntArray neighbours; 
        _pressureField.enclosingNeighbours(imagePoint, neighbours); 
        Eigen::MatrixXd V(8,8); 
        Tuple3i indicesBuffer;
        Vector3d positionBuffer; 
        Eigen::VectorXd RHS(8); 
        const REAL boundaryPressureGradientNormal = _objects->Get(objectID).EvaluateBoundaryAcceleration(boundaryPoint, erectedNormal, simulationTime - 0.5*timeStep)*(-density); 

        for (size_t row=0; row<neighbours.size(); row++) 
        {
            // TODO 
            if (!_isBulkCell[neighbours[row]]) 
                throw std::runtime_error("**ERROR** one of the interpolation stencil for fresh cell is not bulk, this is not handled"); 
            indicesBuffer = _pressureField.cellIndex(neighbours[row]); 
            positionBuffer= _pressureField.cellPosition(indicesBuffer); 
            if (neighbours[row] != cell_idx) // not self
            {
                FillVandermondeRegular(row, positionBuffer, V);
                RHS(row) = p(neighbours[row], 0); 
            }
            else 
            {
                FillVandermondeBoundary(row, boundaryPoint, erectedNormal, V);
                RHS(row) = boundaryPressureGradientNormal; 
            }
        }

        // coefficient for the interpolant
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
        //Eigen::VectorXd C = V.lu().solve(RHS); 
        Eigen::VectorXd C = svd.solve(RHS); 

        // evaluate the interpolant at cell position
        Eigen::VectorXd coordinateVector(8); 
        FillVandermondeRegular(cellPosition, coordinateVector);
        const REAL pressureFreshCell = C.dot(coordinateVector); 
        p(cell_idx, 0) = pressureFreshCell; 
    }
}

//##############################################################################
// Interpolate the fresh velocity cell
//
// Note: 
//   fresh velocity cell can be inside the object, so make sure the
//   interpolation procedure is robust
//##############################################################################
void MAC_Grid::InterpolateFreshVelocityCell_GhostCell(MATRIX &v, const int &dim, const REAL &timeStep, const REAL &simulationTime)
{
    const int N = _velocityInterfacialCells[dim].size(); 
    // essentially same procedures as pressure cells
    for (int interfacial_cell_idx = 0; interfacial_cell_idx < N; ++interfacial_cell_idx) 
    {
        const int cell_idx = _velocityInterfacialCells[dim].at(interfacial_cell_idx);
        if (_toggledVelocityInterfacialCells[dim].at(cell_idx) != 1) // no history problem, skip
            continue; 

        const Vector3d cellPosition = _velocityField[dim].cellPosition(cell_idx);
        REAL distance; 
        int objectID; 
        _objects->LowestObjectDistance(cellPosition, distance, objectID); 
        Vector3d imagePoint, boundaryPoint, erectedNormal; 
        REAL distanceToMesh;
        _objects->Get(objectID).FindImageFreshCell(cellPosition, imagePoint, boundaryPoint, erectedNormal, distanceToMesh); 

        // prepare interpolation stencils, this part is similar to
        // the vandermonde part in ghost cell pressure update
        IntArray neighbours; 
        Tuple3i indicesBuffer;
        Vector3d positionBuffer; 

        // interpolate using trilinear interpolant
        // boundary velocity will be scaled by the normal difference
        _velocityField[dim].enclosingNeighbours(imagePoint, neighbours); 
        const REAL boundaryVelocity = _objects->Get(objectID).EvaluateBoundaryVelocity(boundaryPoint, erectedNormal, simulationTime - timeStep) * _interfacialBoundaryCoefficients[dim].at(interfacial_cell_idx); 
        Eigen::MatrixXd V(8,8); 
        Eigen::VectorXd RHS(8); 
        int row_cell_idx = -1; 
        for (size_t row=0; row<neighbours.size(); row++) 
        {
            int neighbourHasHistory = 0;
            if (!_isVelocityBulkCell[dim].at(neighbours[row]))
            {
                if (!_isVelocityInterfacialCell[dim].at(neighbours[row]))
                {
                    neighbourHasHistory = 1; 
                }
                else  // is interfacial, but still need to check whether it also becomes interfacial at the same time step
                {
                    // neighbouring interfacial cell without valid history
                    if (_toggledVelocityInterfacialCells[dim].at(neighbours[row]) == 1 && neighbours.at(row) != cell_idx) 
                        neighbourHasHistory = 2; 
                }
            }
            if (neighbourHasHistory != 0) 
                throw std::runtime_error("**ERROR** one of the interpolation stencil for velocity fresh cell has invalid history. Error code: " + std::to_string(neighbourHasHistory)); 

            indicesBuffer = _velocityField[dim].cellIndex(neighbours.at(row)); 
            positionBuffer= _velocityField[dim].cellPosition(indicesBuffer); 
            if (neighbours.at(row) != cell_idx) // not self
            {
                FillVandermondeRegular(row, positionBuffer, V);
                RHS(row) = v(neighbours.at(row), 0); 
            }
            else 
            {
                row_cell_idx = row; 
                FillVandermondeRegular(row, boundaryPoint, V); 
                RHS(row) = boundaryVelocity; 
            }
        }
        // coefficient for the interpolant
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
        //Eigen::VectorXd C = V.lu().solve(RHS); 
        Eigen::VectorXd C = svd.solve(RHS); 

        Eigen::VectorXd coordinateVector(8); 
        FillVandermondeRegular(imagePoint, coordinateVector);
        const REAL v_IP = C.dot(coordinateVector); 
        const REAL v_BI = boundaryVelocity;
        if (distanceToMesh > DISTANCE_TOLERANCE) // cell position is located outside the boundary: the erected normal is BI-CU-IP
            v(cell_idx, 0) = 0.5 * (v_IP + v_BI); 
        else // cell position is located inside the boundary: the erected normal is CU-BI-IP
            v(cell_idx, 0) = 2.0 * v_BI - v_IP;
    }
}

void MAC_Grid::classifyCells( bool useBoundary )
{
    int      numPressureCells = _pressureField.numCells();
    int      numVcells[] = { _velocityField[ 0 ].numCells(), _velocityField[ 1 ].numCells(), _velocityField[ 2 ].numCells() };
    Vector3d cellPos;
    IntArray neighbours;

    neighbours.reserve( ScalarField::NUM_NEIGHBOURS );

    _isBulkCell.clear();
    _isGhostCell.clear();
    _isVelocityInterfacialCell[ 0 ].clear();
    _isVelocityInterfacialCell[ 1 ].clear();
    _isVelocityInterfacialCell[ 2 ].clear();
    _isVelocityBulkCell[ 0 ].clear();
    _isVelocityBulkCell[ 1 ].clear();
    _isVelocityBulkCell[ 2 ].clear();

    _isBulkCell.resize( numPressureCells, false );
    _isGhostCell.resize( numPressureCells, false );

    _isVelocityInterfacialCell[ 0 ].resize( numVcells[ 0 ], false );
    _isVelocityInterfacialCell[ 1 ].resize( numVcells[ 1 ], false );
    _isVelocityInterfacialCell[ 2 ].resize( numVcells[ 2 ], false );
    _isVelocityBulkCell[ 0 ].resize( numVcells[ 0 ], false );
    _isVelocityBulkCell[ 1 ].resize( numVcells[ 1 ], false );
    _isVelocityBulkCell[ 2 ].resize( numVcells[ 2 ], false );

    _bulkCells.clear();
    _ghostCells.clear();

    _velocityInterfacialCells[ 0 ].clear();
    _velocityInterfacialCells[ 1 ].clear();
    _velocityInterfacialCells[ 2 ].clear();
    _velocityBulkCells[ 0 ].clear();
    _velocityBulkCells[ 1 ].clear();
    _velocityBulkCells[ 2 ].clear();

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

    if ( !useBoundary )
    {
        std::cout << "No classification (not using boundary)\n"; 
        for ( int cell_idx = 0; cell_idx < numPressureCells; cell_idx++ )
        {
            _isBulkCell[ cell_idx ] = true;
            _bulkCells.push_back( cell_idx );
        }

        for ( int dimension = 0; dimension < 3; dimension++ )
        {
            for ( int cell_idx = 0; cell_idx < _velocityField[ dimension ].numCells(); cell_idx++ )
            {
                _isVelocityBulkCell[dimension][cell_idx] = true; 
                _isVelocityInterfacialCell[dimension][cell_idx] = false; 
                _velocityBulkCells[dimension].push_back(cell_idx); 
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

        if ( _isBulkCell[ cell_idx ] )
        {
            _bulkCells.push_back( cell_idx );
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

            cell_coordinates = _velocityField[ dimension ].cellIndex( cell_idx );

            // We will assume that boundary cells are not interfacial
            if ( cell_coordinates[ dimension ] == 0
                    || cell_coordinates[ dimension ]
                    == _velocityField[ dimension ].cellDivisions()[ dimension ] - 1 )
            {
                _isVelocityBulkCell[ dimension ][ cell_idx ] = true;
                _velocityBulkCells[ dimension ].push_back( cell_idx );

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

                _velocityBulkCells[ dimension ].push_back( cell_idx );
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
                Vector3d position = _velocityField[ dimension ].cellPosition( cell_idx );
                Vector3d gradient = _boundaryFields[ boundaryObject ]->gradient( position );

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
                Vector3d position = _velocityField[ dimension ].cellPosition( cell_idx );
                Vector3d gradient = _boundaryFields[ boundaryObject ]->gradient( position );

                normal[ dimension ] = 1.0;
                gradient.normalize();

                REAL coefficient = normal.dotProduct( gradient );

                //TRACE_ASSERT( coefficient >= 0.0 );

                _interfacialBoundaryCoefficients[ dimension ].push_back( coefficient );

            }
        }
    }


    printf( "MAC_Grid: classifyCells:\n" );
    printf( "\tFound %d bulk cells\n", (int)_bulkCells.size() );
    printf( "\tFound %d ghost cells\n", (int)_ghostCells.size() );
    printf( "\tFound %d v_x interfacial cells\n",
            (int)_velocityInterfacialCells[ 0 ].size() );
    printf( "\tFound %d v_y interfacial cells\n",
            (int)_velocityInterfacialCells[ 1 ].size() );
    printf( "\tFound %d v_z interfacial cells\n",
            (int)_velocityInterfacialCells[ 2 ].size() );
    printf( "\tFound %d v_x bulk cells\n", (int)_velocityBulkCells[ 0 ].size() );
    printf( "\tFound %d v_y bulk cells\n", (int)_velocityBulkCells[ 1 ].size() );
    printf( "\tFound %d v_z bulk cells\n", (int)_velocityBulkCells[ 2 ].size() );

    if (_useGhostCellBoundary)
    {
        visualizeClassifiedCells(); 
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
//##############################################################################
void MAC_Grid::classifyCellsDynamic(MATRIX &pFull, MATRIX (&p)[3], MATRIX (&v)[3], const bool &useBoundary, const bool &verbose)
{
    const int numPressureCells = _pressureField.numCells(); 
    Vector3d cellPos;
    IntArray neighbours;
    neighbours.reserve( ScalarField::NUM_NEIGHBOURS );

    // step 1 : clear relevant arrays
    _bulkCells.clear(); 
    for (int dim=0; dim<3; ++dim) 
        _velocityBulkCells[dim].clear(); 
    //std::fill(_isBulkCell.begin(), _isBulkCell.end(), true); 
    //std::fill(_containingObject.begin(), _containingObject.end(), -1); 
    // need to know whether at last time step history is valid

    // step 2 : not using boundary, classification is trivial
    if (!useBoundary)
    {
        if (verbose) 
            std::cout << "No classification (not using boundary)\n"; 
        for (int cell_idx = 0; cell_idx < numPressureCells; ++cell_idx)
            _bulkCells.push_back( cell_idx );
        for (int dimension = 0; dimension < 3; ++dimension)
        {
            for (int cell_idx = 0; cell_idx < _velocityField[dimension].numCells(); ++cell_idx)
            {
                _isVelocityBulkCell[dimension][cell_idx] = true; 
                _isVelocityInterfacialCell[dimension][cell_idx] = false; 
                _velocityBulkCells[dimension].push_back(cell_idx); 
            }
        }
        return;
    }

    // step 3 : classify pressure bulk cells 
    for ( int cell_idx = 0; cell_idx < numPressureCells; cell_idx++ )
    {
        // first establish if it has valid history
        if (_isBulkCell.at(cell_idx)) // was bulk -> valid history
        {
            _pressureCellHasValidHistory.at(cell_idx) = true; 
        }
        else 
        {
            if (_useGhostCellBoundary && _isGhostCell.at(cell_idx)) // was ghost -> valid history
                _pressureCellHasValidHistory.at(cell_idx) = true;
            else 
                _pressureCellHasValidHistory.at(cell_idx) = false;
        }

        cellPos = _pressureField.cellPosition( cell_idx );
        // Check all boundary fields to see if this is a bulk cell
        _containingObject.at(cell_idx) = _objects->OccupyByObject(cellPos); 
        const bool newIsBulkCell = (_containingObject.at(cell_idx)>=0 ? false : true); 

        _isBulkCell.at(cell_idx) = newIsBulkCell; 
        if (newIsBulkCell)
            _bulkCells.push_back(cell_idx); 

        // clear history for solid pressure
        if (!newIsBulkCell && !_useGhostCellBoundary) 
        {
            pFull(cell_idx, 0) = 0;
            p[0](cell_idx, 0) = 0; 
            p[1](cell_idx, 0) = 0; 
            p[2](cell_idx, 0) = 0; 
        }

        //if (_isBulkCell[cell_idx] && !newIsBulkCell)  // turning to solid cell
        //{
        //    _toggledBulkCells[cell_idx] = -1; 
        //    pFull(cell_idx, 0) = 0; // clear pressure solid cell
        //    p[0](cell_idx, 0) = 0; 
        //    p[1](cell_idx, 0) = 0; 
        //    p[2](cell_idx, 0) = 0; 
        //} 
        //else if (!_isBulkCell[cell_idx] && newIsBulkCell) // turning to bulk cell
        //    _toggledBulkCells[cell_idx] = 1; 
        //else
        //    _toggledBulkCells[cell_idx] = 0; 

        //if (!newIsBulkCell) 
        //{
        //    _isBulkCell[cell_idx] = false; 
        //    // clear the pressure cells that are bulk
        //    _containingObject[cell_idx] = indexOccupyObject; 
        //}
        //else 
        //{
        //    _isBulkCell[cell_idx] = true; 
        //    _containingObject[cell_idx] = -1; 
        //    _bulkCells.push_back( cell_idx );
        //}
    }

    // step 4a : classify pressure ghost cells if flag is on
    if (_useGhostCellBoundary) 
    {
        std::fill(_toggledGhostCells.begin(), _toggledGhostCells.end(), 0); 
        _ghostCells.clear(); 
        // examine ghost cells 
        for ( int cell_idx = 0; cell_idx < numPressureCells; ++cell_idx )
        {
            bool newIsGhostCell = false; 
            if (_isBulkCell[cell_idx]) // transition to bulk, if this is from solid then we need to worry, if from ghost cell ignore
            {
                if (_isGhostCell[cell_idx] == false && _toggledBulkCells[cell_idx] != 0) // it changes and not changed from ghost, so must change from solid
                {
                    throw std::runtime_error("**ERROR** turning to bulk directly from solid. not handled."); 
                }
                else
                {
                    _toggledGhostCells[cell_idx] = 0; 
                    _isGhostCell[cell_idx] = false; 
                }
                continue;
            }
            _pressureField.cellNeighbours( cell_idx, neighbours );
            for (size_t neighbour_idx = 0; neighbour_idx<neighbours.size(); ++neighbour_idx)
            {
                if (_isBulkCell[neighbours[neighbour_idx]])
                {
                    // We have a neighbour outside of the interior object, so
                    // this is a ghost cell
                    newIsGhostCell = true; 
                    _ghostCells.push_back(cell_idx); 
                    break;
                }
            }
            if (_isGhostCell[cell_idx] && !newIsGhostCell)
                _toggledGhostCells[cell_idx] = -1;
            else if (!_isGhostCell[cell_idx] && newIsGhostCell) 
                _toggledGhostCells[cell_idx] = 1; 
            else 
                _toggledGhostCells[cell_idx] = 0;
            _isGhostCell[ cell_idx ] = newIsGhostCell;
        }
        ComputeGhostCellInverseMap(); 
        // TODO can make smarter 
        for (int dim=0; dim<3; ++dim)
            for (int ii=0; ii<_velocityField[dim].numCells(); ++ii)
                _velocityBulkCells[dim].push_back(ii);
    }
    // step 4b : classify velocity cells for rasterized boundary
    else
    {
        for (int dim=0; dim<3; ++dim)
        {
            //std::fill(_isVelocityBulkCell[dim].begin(), _isVelocityBulkCell[dim].end(), true); 
            //std::fill(_isVelocityInterfacialCell[dim].begin(), _isVelocityInterfacialCell[dim].end(), false); 
            std::fill(_toggledVelocityInterfacialCells[dim].begin(), _toggledVelocityInterfacialCells[dim].end(), 0); 
            _velocityInterfacialCells[dim].clear(); 
            _interfacialBoundaryIDs[dim].clear(); 
            _interfacialBoundaryDirections[dim].clear(); 
            _interfacialBoundaryCoefficients[dim].clear(); 
            _velocityBulkCells[dim].clear(); 
        }
        for (int dimension = 0; dimension < 3; dimension++)
        {
            for (int cell_idx = 0; cell_idx < _velocityField[dimension].numCells(); cell_idx++)
            {
                Tuple3i cell_coordinates;
                int     pressure_cell_idx1;
                int     pressure_cell_idx2;

                cell_coordinates = _velocityField[ dimension ].cellIndex( cell_idx );

                // first establish if it has valid history
                if (_isVelocityBulkCell[dimension].at(cell_idx) || _isVelocityInterfacialCell[dimension].at(cell_idx)) 
                    _velocityCellHasValidHistory[dimension].at(cell_idx) = true; 
                else
                    _velocityCellHasValidHistory[dimension].at(cell_idx) = false; 

                // We will assume that boundary cells are not interfacial
                if ( cell_coordinates[ dimension ] == 0
                  || cell_coordinates[ dimension ] == _velocityField[ dimension ].cellDivisions()[ dimension ] - 1 )
                {
                    _velocityBulkCells[ dimension ].push_back( cell_idx );
                    _isVelocityInterfacialCell[dimension][cell_idx] = false; 
                    continue;
                }

                // Look at our neighbours in the pressure field
                cell_coordinates[ dimension ] -= 1;
                pressure_cell_idx1 = _pressureField.cellIndex( cell_coordinates );
                cell_coordinates[ dimension ] += 1;
                pressure_cell_idx2 = _pressureField.cellIndex( cell_coordinates );

                if (_isBulkCell[ pressure_cell_idx1 ] && _isBulkCell[ pressure_cell_idx2 ])
                {
                    // Both pressure cell neighbours are bulk cells, so this is
                    // a bulk cell too
                    _isVelocityInterfacialCell[dimension][cell_idx] = false; 
                    _isVelocityBulkCell[dimension][cell_idx] = true; 
                    _velocityBulkCells[ dimension ].push_back( cell_idx );
                }
                else if ( _isBulkCell[ pressure_cell_idx1 ] 
                        && !_isBulkCell[ pressure_cell_idx2 ] )
                {
                    // Only one neighbour is inside the domain, so this must
                    // be an interfacial cell
                    if (_isVelocityBulkCell[dimension][cell_idx] == false && _isVelocityInterfacialCell[dimension][cell_idx] == false)
                        _toggledVelocityInterfacialCells[dimension][cell_idx] = 1; // was solid cell
                    _isVelocityInterfacialCell[ dimension ][ cell_idx ] = true;

                    // Get the object ID for the boundary that we are adjacent to
                    const int boundaryObject = _containingObject[ pressure_cell_idx2 ];

                    TRACE_ASSERT( boundaryObject >= 0 );

                    _velocityInterfacialCells[ dimension ].push_back( cell_idx );
                    _interfacialBoundaryIDs[ dimension ].push_back( boundaryObject );
                    //_interfacialBoundaryDirections[ dimension ].push_back( -1.0 );
                    _interfacialBoundaryDirections[ dimension ].push_back( 1.0 );

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

                    _isVelocityBulkCell[dimension][cell_idx] = false; 
                    _interfacialBoundaryCoefficients[ dimension ].push_back( coefficient );
                }
                else if ( !_isBulkCell[ pressure_cell_idx1 ]
                        && _isBulkCell[ pressure_cell_idx2 ] )
                {
                    // Only one neighbour is inside the domain, so this must
                    // be an interfacial cell
                    if (_isVelocityBulkCell[dimension][cell_idx] == false && _isVelocityInterfacialCell[dimension][cell_idx] == false) 
                        _toggledVelocityInterfacialCells[dimension][cell_idx] = 1; // was solid cell
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
                    Vector3d position = _velocityField[ dimension ].cellPosition( cell_idx );
                    Vector3d gradient; 
                    _objects->ObjectNormal(boundaryObject, position, gradient); 
                    gradient.normalize();

                    normal[ dimension ] = 1.0;
                    REAL coefficient = normal.dotProduct( gradient );

                    //TRACE_ASSERT( coefficient >= 0.0 );

                    _isVelocityBulkCell[dimension][cell_idx] = false; 
                    _interfacialBoundaryCoefficients[ dimension ].push_back( coefficient );
                }
                else // both sides aren't bulk, this is a solid cell
                {
                    if (_isVelocityInterfacialCell[dimension][cell_idx] = true)
                        _toggledVelocityInterfacialCells[dimension][cell_idx] = -1; 
                    _isVelocityBulkCell[dimension][cell_idx] = false; 
                    _isVelocityInterfacialCell[dimension][cell_idx] = false;

                    v[dimension](cell_idx, 0) = 0.0; // clear solid velocity cell
                }
            }
        }
    }

    if (verbose) 
    {
        printf( "MAC_Grid: classifyCellsDynamic:\n" );
        printf( "\tFound %d bulk cells\n", (int)_bulkCells.size() );
        printf( "\tFound %d ghost cells\n", (int)_ghostCells.size() );
        printf( "\tFound %d v_x interfacial cells\n",
                (int)_velocityInterfacialCells[ 0 ].size() );
        printf( "\tFound %d v_y interfacial cells\n",
                (int)_velocityInterfacialCells[ 1 ].size() );
        printf( "\tFound %d v_z interfacial cells\n",
                (int)_velocityInterfacialCells[ 2 ].size() );
        printf( "\tFound %d v_x bulk cells\n", (int)_velocityBulkCells[ 0 ].size() );
        printf( "\tFound %d v_y bulk cells\n", (int)_velocityBulkCells[ 1 ].size() );
        printf( "\tFound %d v_z bulk cells\n", (int)_velocityBulkCells[ 2 ].size() );
    }
}

//##############################################################################
// Classify types of cells at each step using axis-aligned bounding box 
//
// Note: assumed cells has been classified properly and just need to be updated 
// base on the new object positions. Therefore, not all cells will be checked.
//
///##############################################################################
void MAC_Grid::classifyCellsDynamicAABB(const bool &useBoundary, MATRIX &p, const bool &verbose)
{
    _ghostCells.clear();

    // if not using boundaries, return immediately
    if (!useBoundary)
        return;

    // get all bounding boxes and iteration range
    const int N = _objects->N(); 
    std::vector<ScalarField::RangeIndices> indices(N); 
    for (int object_id=0; object_id<N; ++object_id)
    {
        const FDTD_MovableObject::BoundingBox &unionBBox = _objects->Get(object_id).GetUnionBBox();

        const Vector3d maxBound = unionBBox.maxBound + 2.0*_cellSize; 
        const Vector3d minBound = unionBBox.minBound - 2.0*_cellSize; 
        _pressureField.GetIterationBox(minBound, maxBound, indices[object_id]); 
    } 

    // classify only the subset indicated by bounding box
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (int bbox_id=0; bbox_id<N; ++bbox_id)
    {
        Vector3d cellPos;
        int ii,jj,kk; 
        FOR_ALL_3D_GRID_VECTOR3(indices[bbox_id].startIndex, indices[bbox_id].dimensionIteration, ii, jj, kk)
        {
            const Tuple3i cellIndices(ii,jj,kk);
            const int cell_idx = _pressureField.cellIndex(cellIndices); 
            cellPos = _pressureField.cellPosition(cell_idx);
            _containingObject[cell_idx] = _objects->OccupyByObject(cellPos); 
            // need to update the id of the cell 
            const bool newIsBulkCell = (_containingObject[cell_idx]>=0 ? false : true); 
            if (_isBulkCell[cell_idx] && !newIsBulkCell)  // turning to solid cell
                _toggledBulkCells[cell_idx] = -1; 
            else if (!_isBulkCell[cell_idx] && newIsBulkCell) // turning to bulk cell
                _toggledBulkCells[cell_idx] = 1; 
            else // identity unchanged
                _toggledBulkCells[cell_idx] = 0; 
            _isBulkCell[cell_idx] = newIsBulkCell; 

            if (!newIsBulkCell)  // reset the pressure for all solid cells; 
                p(cell_idx, 0) = 0.0;
        }
    }

    // step 4a 
    if (_useGhostCellBoundary) 
    {
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
        for (int bbox_id=0; bbox_id<N; ++bbox_id)
        {
            int ii,jj,kk; 
            IntArray neighbours;
            neighbours.reserve( ScalarField::NUM_NEIGHBOURS );
            FOR_ALL_3D_GRID_VECTOR3(indices[bbox_id].startIndex, indices[bbox_id].dimensionIteration, ii, jj, kk)
            {
                const Tuple3i cellIndices(ii,jj,kk);
                const int cell_idx = _pressureField.cellIndex(cellIndices); 
                bool newIsGhostCell = false; 
                // classify ghost cells
                if (!_isBulkCell[cell_idx]) 
                {
                    _pressureField.cellNeighbours(cell_idx, neighbours);
                    const int neighbourSize = neighbours.size();
                    for (int neighbour_idx = 0; neighbour_idx<neighbourSize; ++neighbour_idx)
                    {
                        if (_isBulkCell[neighbours[neighbour_idx]])
                        {
                            // We have a neighbour outside of the interior object, so
                            // this is a ghost cell
                            newIsGhostCell = true; 
#pragma omp critical
                            _ghostCells.push_back(cell_idx); 
                            break;
                        }
                    }
                }
                if (_isGhostCell[cell_idx] || newIsGhostCell) 
                    std::cout << "toggle ghost cell = " << _isGhostCell[cell_idx] << " " << newIsGhostCell << std::endl;
                if (_isGhostCell[cell_idx] && !newIsGhostCell)
                    _toggledGhostCells[cell_idx] = -1;
                else if (!_isGhostCell[cell_idx] && newIsGhostCell) 
                    _toggledGhostCells[cell_idx] = 1; 
                else 
                    _toggledGhostCells[cell_idx] = 0;
                _isGhostCell[cell_idx] = newIsGhostCell; 
            }
        }
        STL_Wrapper::VectorSortAndTrimInPlace(_ghostCells);
        ComputeGhostCellInverseMap(); 
    }
    // step 4b 
    else 
    {
        throw std::runtime_error("**ERROR** not supporting rasterized cell reclassification. you can use the non-AABB accelerated version"); 
    }

    if (verbose) 
    {
        printf( "MAC_Grid: classifyCellsDynamicAABB:\n" );
        printf( "\tFound %d bulk cells\n", (int)_bulkCells.size() );
        printf( "\tFound %d ghost cells\n", (int)_ghostCells.size() );
        printf( "\tFound %d v_x interfacial cells\n",
                (int)_velocityInterfacialCells[ 0 ].size() );
        printf( "\tFound %d v_y interfacial cells\n",
                (int)_velocityInterfacialCells[ 1 ].size() );
        printf( "\tFound %d v_z interfacial cells\n",
                (int)_velocityInterfacialCells[ 2 ].size() );
        printf( "\tFound %d v_x bulk cells\n", (int)_velocityBulkCells[ 0 ].size() );
        printf( "\tFound %d v_y bulk cells\n", (int)_velocityBulkCells[ 1 ].size() );
        printf( "\tFound %d v_z bulk cells\n", (int)_velocityBulkCells[ 2 ].size() );
    }
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

REAL MAC_Grid::PML_absorptionCoefficient( const Vector3d &x, REAL absorptionWidth, int dimension )
{
    //return 0.0;

    if (dimension!=2 && _cornellBoxBoundaryCondition) // skip all but z-direction face 
        return 0.0;

    const BoundingBox         &bbox = _pressureField.bbox();
    REAL                       h;
    REAL                       a2 = absorptionWidth * absorptionWidth;

    h = x[ dimension ] - bbox.minBound()[ dimension ];
    if ( h <= absorptionWidth )
    {
        REAL                     dist = absorptionWidth - h;

        //return _PML_absorptionStrength * ( absorptionWidth - h ) / absorptionWidth;
        return _PML_absorptionStrength * dist * dist / a2;
    }

    h = bbox.maxBound()[ dimension ] - x[ dimension ];
    if ( h <= absorptionWidth )
    {
        REAL                     dist = absorptionWidth - h;
        if (_cornellBoxBoundaryCondition)
            return 0.0; // skip face at +z face
        else 
            return _PML_absorptionStrength * dist * dist / a2;
    }

    return 0.0;
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
       << " number pressure cells  : " << grid.numPressureCells() << "\n"
       << " number velocity cells x: " << grid.numVelocityCellsX() << "\n"
       << " number velocity cells y: " << grid.numVelocityCellsY() << "\n"
       << " number velocity cells z: " << grid.numVelocityCellsZ() << "\n"
       << "--------------------------------------------------------------------------------" 
       << std::flush; 
    return os; 
}








