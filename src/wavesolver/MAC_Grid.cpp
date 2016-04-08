//////////////////////////////////////////////////////////////////////
// MAC_Grid.cpp: Implementation of the MAC_Grid class
//
//////////////////////////////////////////////////////////////////////

#include "MAC_Grid.h"

#include <utils/IO.h>
#include <utils/trace.h>

#include <math/LeastSquareSurface.h>
#include <utils/STL_Wrapper.h> 
#include <boost/timer/timer.hpp>

#include <distancefield/trilinearInterpolation.h> 
#include <linearalgebra/SparseLinearSystemSolver.h>
//#include <graph/UndirectedGraph.h> 

#include <unistd.h> 

#ifdef USE_OPENMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor
//
// Provide the size of the domain, finite difference division
// size, and a signed distance function for the interior boundary
//////////////////////////////////////////////////////////////////////
MAC_Grid::MAC_Grid( const BoundingBox &bbox, REAL cellSize,
                    const TriMesh &mesh,
                    const DistanceField &boundarySDF,
                    REAL distanceTolerance, int N )
    : 
      _distanceTolerance( distanceTolerance ),
      _pressureField( bbox, cellSize ),
      _ghostCellsInverseComputed(false),
      _N( N ),
      _PML_absorptionWidth( 1.0 ),
      _PML_absorptionStrength( 0.0 )
{
    Vector3d                   xMin, yMin, zMin;
    Tuple3i                    xDivs, yDivs, zDivs;

    _boundaryFields.push_back( &boundarySDF );
    _boundaryMeshes.push_back( &mesh );

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

    printf( "Initialized MAC_Grid with %d cells\n", _pressureField.numCells() );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
MAC_Grid::MAC_Grid( const BoundingBox &bbox, REAL cellSize,
                    vector<const TriMesh *> &meshes,
                    vector<const DistanceField *> &boundaryFields,
                    REAL distanceTolerance, int N )
    : 
      _boundaryFields( boundaryFields ),
      _boundaryMeshes( meshes ),
      _distanceTolerance(distanceTolerance),
      _pressureField( bbox, cellSize ),
      _ghostCellsInverseComputed(false),
      _N( N ),
      _PML_absorptionWidth( 1.0 ),
      _PML_absorptionStrength( 0.0 )
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

    printf( "Initialized MAC_Grid with %d cells\n", _pressureField.numCells() );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
MAC_Grid::~MAC_Grid()
{
}

//////////////////////////////////////////////////////////////////////
// Initializes field using a rasterized version of the boundary
//////////////////////////////////////////////////////////////////////
void MAC_Grid::initFieldRasterized( bool useBoundary )
{
    cout << "MAC_Grid::initFieldRasterized: Classifying cells" << endl << endl;
    classifyCells( useBoundary );
}

//////////////////////////////////////////////////////////////////////
// Width is the number of cells we wish to absorb in
//////////////////////////////////////////////////////////////////////
void MAC_Grid::setPMLBoundaryWidth( REAL width, REAL strength )
{
    _PML_absorptionWidth = width * _pressureField.cellSize();
    _PML_absorptionStrength = strength;
}

//////////////////////////////////////////////////////////////////////
// Gets the derivative for the given velocity field, based on the
// input pressure value
//////////////////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////////////////////////
// Gets the derivative of the pressure field using the current
// velocity field values in each direction
//////////////////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////////////////////////
// Performs a velocity update in the given direction, as detailed
// by Liu et al. (equation (14))
//////////////////////////////////////////////////////////////////////
void MAC_Grid::PML_velocityUpdate( const MATRIX &p, const BoundaryEvaluator &bc,
                                   MATRIX &v, int dimension,
                                   REAL t, REAL timeStep, REAL density )
{
    const IntArray        &bulkCells                    = _velocityBulkCells[ dimension ];
    const IntArray        &interfacialCells             = _velocityInterfacialCells[ dimension ];
    const ScalarField     &field                        = _velocityField[ dimension ];
    const IntArray        &interfaceBoundaryIDs         = _interfacialBoundaryIDs[ dimension ];
    const FloatArray      &interfaceBoundaryDirections  = _interfacialBoundaryDirections[ dimension ];
    const FloatArray      &interfaceBoundaryCoefficients= _interfacialBoundaryCoefficients[ dimension ];

    //// TODO debug quick fix for getting the acceleration
    //MATRIX v_old = v; 

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
        Vector3d             cell_position = field.cellPosition( cell_indices );
        REAL                 absorptionCoefficient;
        REAL                 updateCoefficient;
        REAL                 gradientCoefficient;

        absorptionCoefficient = PML_absorptionCoefficient( cell_position,
                _PML_absorptionWidth,
                dimension );
        updateCoefficient = PML_velocityUpdateCoefficient( absorptionCoefficient,
                timeStep );
        gradientCoefficient = PML_pressureGradientCoefficient(
                absorptionCoefficient,
                timeStep,
                _pressureField.cellSize(),
                density );

        for ( int i = 0; i < _N; i++ )
        {
            v( cell_idx, i ) *= updateCoefficient;
        }

        // If we're on a boundary, then set the derivative to zero (since
        // we should have d_n p = 0 on the boundary
        if ( cell_indices[ dimension ] == 0
                || cell_indices[ dimension ] == field.cellDivisions()[ dimension ] - 1 )
        {
            for ( int i = 0; i < _N; i++ )
            {
                v( cell_idx, i ) = 0.0;
            }

            continue;
        }

        // Sample the derivative of the pressure field in the necessary
        // direction
        neighbour_idx = _pressureField.cellIndex( cell_indices );

        for ( int i = 0; i < _N; i++ )
        {
            v( cell_idx, i )
                += gradientCoefficient * p( neighbour_idx, i );
            /// _pressureField.cellSize();
        }

        cell_indices[ dimension ] -= 1;
        neighbour_idx = _pressureField.cellIndex( cell_indices );

        for ( int i = 0; i < _N; i++ )
        {
            v( cell_idx, i )
                -= gradientCoefficient * p( neighbour_idx, i );
            /// _pressureField.cellSize();
        }
    }



    if (_useGhostCellBoundary)
    {

        // with ghost-cell approach the interfacial cell can be treated as regular
        // cell 
        #ifdef USE_OPENMP
        #pragma omp parallel for schedule(static) default(shared)
        #endif
        for ( size_t interfacial_cell_idx = 0; interfacial_cell_idx < interfacialCells.size();
                interfacial_cell_idx++ )
        {
            int                  cell_idx = interfacialCells[ interfacial_cell_idx ];
            int                  neighbour_idx;
            Tuple3i              cell_indices = field.cellIndex( cell_idx );
            Vector3d             cell_position = field.cellPosition( cell_indices );
            REAL                 absorptionCoefficient;
            REAL                 updateCoefficient;
            REAL                 gradientCoefficient;

            absorptionCoefficient = PML_absorptionCoefficient( cell_position,
                    _PML_absorptionWidth,
                    dimension );
            updateCoefficient = PML_velocityUpdateCoefficient( absorptionCoefficient,
                    timeStep );
            gradientCoefficient = PML_pressureGradientCoefficient(
                    absorptionCoefficient,
                    timeStep,
                    _pressureField.cellSize(),
                    density );

            for ( int i = 0; i < _N; i++ )
            {
                v( cell_idx, i ) *= updateCoefficient;
            }

            // If we're on a boundary, then set the derivative to zero (since
            // we should have d_n p = 0 on the boundary
            if ( cell_indices[ dimension ] == 0
                    || cell_indices[ dimension ] == field.cellDivisions()[ dimension ] - 1 )
            {
                for ( int i = 0; i < _N; i++ )
                {
                    v( cell_idx, i ) = 0.0;
                }

                continue;
            }


            // Sample the derivative of the pressure field in the necessary
            // direction
            neighbour_idx = _pressureField.cellIndex( cell_indices );

            for ( int i = 0; i < _N; i++ )
            {
                v( cell_idx, i )
                    += gradientCoefficient * p( neighbour_idx, i );
                /// _pressureField.cellSize();
            }


            cell_indices[ dimension ] -= 1;
            neighbour_idx = _pressureField.cellIndex( cell_indices );

            for ( int i = 0; i < _N; i++ )
            {
                v( cell_idx, i )
                    -= gradientCoefficient * p( neighbour_idx, i );
                /// _pressureField.cellSize();
            }

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
        for ( size_t interfacial_cell_idx = 0;
                interfacial_cell_idx < interfacialCells.size();
                interfacial_cell_idx++ )
        {
            int                  cell_idx = interfacialCells[ interfacial_cell_idx ];
            int                  objectID;
            REAL                 bcEval;
            REAL                 coefficient;

            Vector3d             normal( 0.0, 0.0, 0.0 ); // normal of the interfacial cell boundary
            Vector3d             x = field.cellPosition( cell_idx ); // position of the interfacial cell

            normal[ dimension ] = interfaceBoundaryDirections[ interfacial_cell_idx ];

            coefficient = interfaceBoundaryCoefficients[ interfacial_cell_idx ];

            objectID = interfaceBoundaryIDs[ interfacial_cell_idx ];

            for ( int i = 0; i < _N; i++ )
            {
                bcEval = bc( x, normal, objectID, t, i );
                bcEval *= coefficient;

                v( cell_idx, i ) += timeStep * bcEval;
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Performs a pressure update for the given pressure direction,
void MAC_Grid::PML_pressureUpdate( const MATRIX &v, MATRIX &p, int dimension,
                                   REAL timeStep, REAL c, const ExternalSourceEvaluator *sourceEvaluator, const REAL simulationTime, REAL density )
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
        Vector3d                 cell_position;
        REAL                     absorptionCoefficient;
        REAL                     directionalCoefficient;
        REAL                     updateCoefficient;
        REAL                     divergenceCoefficient;

        cell_indices = _pressureField.cellIndex( cell_idx );
        cell_position = _pressureField.cellPosition( cell_idx );
        absorptionCoefficient = PML_absorptionCoefficient( cell_position,
                _PML_absorptionWidth,
                dimension );
        directionalCoefficient = PML_directionalCoefficient( absorptionCoefficient,
                timeStep );
        updateCoefficient = PML_pressureUpdateCoefficient( absorptionCoefficient,
                timeStep,
                directionalCoefficient );
        divergenceCoefficient = PML_divergenceCoefficient( density, c,
                directionalCoefficient );

        // Scale the existing pressure value
        for ( int i = 0; i < _N; i++ )
        {
            p( cell_idx, i ) *= updateCoefficient;
        }

        cell_indices[ dimension ] += 1;
        neighbour_idx = _velocityField[ dimension ].cellIndex( cell_indices );

        for ( int i = 0; i < _N; i++ )
        {
            p( cell_idx, i ) += divergenceCoefficient * v( neighbour_idx, i )
                / _velocityField[ dimension ].cellSize();
        }

        cell_indices[ dimension ] -= 1;
        neighbour_idx = _velocityField[ dimension ].cellIndex( cell_indices );

        for ( int i = 0; i < _N; i++ )
        {
            p( cell_idx, i ) -= divergenceCoefficient * v( neighbour_idx, i )
                / _velocityField[ dimension ].cellSize();
        }


        // evaluate external sources 
        // Liu Eq (16) f6x term
        if ( sourceEvaluator != nullptr )
        {
            for ( int i = 0; i < _N; i++ ) 
            {
                p( cell_idx, i ) += (*sourceEvaluator)( cell_position, simulationTime+0.5*timeStep ) / directionalCoefficient; 
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
void MAC_Grid::PML_pressureUpdateGhostCells( MATRIX &p, const REAL &timeStep, const REAL &c, const BoundaryEvaluator &bc, const REAL &simulationTime, const REAL density)
{

    //std::ostream *of_image = new std::ofstream("image_points.csv"); 
    //std::ostream *of_boundary = new std::ofstream("boundary_points.csv"); 

    // for the ghost-cell coupling
    SparseLinearSystemSolver solver(_ghostCells.size()); 

#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for (size_t ghost_cell_idx=0; ghost_cell_idx<_ghostCells.size(); ghost_cell_idx++) 
    {

        const int       cellIndex      = _ghostCells[ghost_cell_idx]; 
        const Tuple3i   cellIndices    = _pressureField.cellIndex(cellIndex); 
        const Vector3d  cellPosition   = _pressureField.cellPosition(cellIndices); 
        const int       boundaryObject = _containingObject[cellIndex]; 

        // find point BI and IP in the formulation
        Vector3d boundaryPoint, imagePoint, erectedNormal; 

        FindImagePoint(cellPosition, boundaryObject, boundaryPoint, imagePoint, erectedNormal); 
        //*of_boundary << boundaryPoint.x << "," << boundaryPoint.y << "," << boundaryPoint.z << std::endl;
        //*of_image << imagePoint[0] << "," << imagePoint[1] << "," << imagePoint[2] << std::endl;


        // get the box enclosing the image point; 
        IntArray neighbours; 
        _pressureField.enclosingNeighbours(imagePoint, neighbours); 

        assert(neighbours.size()==8);

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

        const double bcEval = bc(boundaryPoint, erectedNormal, boundaryObject, simulationTime, 0) * (-density);  // evaluate at boundarPoint, the boundary is prescribing normal acceleration so scaling is needed for pressure Neumann

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

        // want beta = V^-T b, equivalent to solving V^T x = b, and then x^T
        Eigen::MatrixXd b(1,8); // row vector

        FillVandermondeRegular(0,imagePoint,b); 
        b.transposeInPlace();  // column vector
        V.transposeInPlace(); 

        // TODO svd solve is needed? can I use QR? 
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(V, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::VectorXd beta = svd.solve(b); 
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
        //
#pragma omp critical
        {

            // fill in the sparse matrix
            
            // start by filling the diagonal entry
            solver.StageEntry(ghost_cell_idx, ghost_cell_idx, 1.0);

            // next fill the off-diagonal terms by iterating through the coupled
            // ghost cells for this row. its coefficient is defined by the i-th
            // entry of beta. i=1,...,8 are the index of the neighbour
            for (size_t cc=0; cc<coupledGhostCells.size(); cc++) 
                solver.StageEntry(ghost_cell_idx, coupledGhostCells[cc], -beta(coupledGhostCellsNeighbours[cc])); 

            // next we compute the right-hand-side (RHS) of the equation. looking
            // at the equation, it will be the neumann condition + any uncoupled
            // betas (see google doc for eqution)
            double RHS = - (imagePoint - cellPosition).length() * bcEval; 
            for (size_t uc=0; uc<uncoupledGhostCellsNeighbours.size(); uc++) 
                RHS += beta(uncoupledGhostCellsNeighbours[uc])*pressureNeighbours(uncoupledGhostCellsNeighbours[uc]); 

            // if GC itself is coupled, need to add this to the RHS
            if (hasGC != -1) 
            {
                RHS += beta(hasGC)*pressureNeighbours(hasGC); 
            }

            // finally we fill the rhs of the sparse linear system
            solver.FillInRHS(ghost_cell_idx, RHS); 
        }


    }


    //reinterpret_cast<std::ofstream*>(of_image)->close();
    //reinterpret_cast<std::ofstream*>(of_boundary)->close();
    //exit(1);

    // actually fill-in all the staged entry of the sparse matrix. after this
    // step the matrix is ready for solve. 
    {
        boost::timer::auto_cpu_timer t(" sparse system fill-in takes %w sec\n"); 
        solver.FillIn(); 
    }

    // the actual solve. unfortunately eigen supports only dense RHS solve. 
    Eigen::VectorXd x; 
    {
        boost::timer::auto_cpu_timer t(" linear system solve takes %w sec\n");
        solver.Solve(x); 
        //solver.DenseSolve(x); 
    }


    //COUT_SDUMP(x.maxCoeff());
    //COUT_SDUMP(x.minCoeff());
    //solver.PrintMatrixDense(std::cout); 
    //solver.PrintVectorDense(std::cout);

    // put the solution into matrix p
    {
        boost::timer::auto_cpu_timer t(" recover solution to dense vector takes %w sec\n");
        for (size_t ghost_cell_idx=0; ghost_cell_idx<_ghostCells.size(); ghost_cell_idx++) 
            p(_ghostCells[ghost_cell_idx],0) = x(ghost_cell_idx);
    }
}



//////////////////////////////////////////////////////////////////////
// Samples data from a z slice of the finite difference grid and
// puts it in to a matrix
//////////////////////////////////////////////////////////////////////
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


//////////////////////////////////////////////////////////////////////
// Smooth field given weights
//////////////////////////////////////////////////////////////////////
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

    _ghostCellsInverseComputed = true; 
}

//////////////////////////////////////////////////////////////////////
// Classifies cells as either a bulk cell, ghost cell, or
// interfacial cell
//////////////////////////////////////////////////////////////////////
void MAC_Grid::classifyCells( bool useBoundary )
{
    int                        numPressureCells = _pressureField.numCells();
    int                        numVcells[] = { _velocityField[ 0 ].numCells(),
        _velocityField[ 1 ].numCells(),
        _velocityField[ 2 ].numCells() };
    Vector3d                   cellPos;
    IntArray                   neighbours;

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

    //IntArray                   containingObject( numPressureCells, -1 );

    if ( !useBoundary )
    {
        std::cout << "classifying cells: not using boundary\n"; 
        for ( int cell_idx = 0; cell_idx < numPressureCells; cell_idx++ )
        {
            _isBulkCell[ cell_idx ] = true;
            _bulkCells.push_back( cell_idx );
        }

        // Classify interfacial cells in each velocity field
        for ( int dimension = 0; dimension < 3; dimension++ )
        {
            for ( int cell_idx = 0; cell_idx < _velocityField[ dimension ].numCells(); cell_idx++ )
            {
                Tuple3i                cell_coordinates;
                int                    pressure_cell_idx1;
                int                    pressure_cell_idx2;

                cell_coordinates = _velocityField[ dimension ].cellIndex( cell_idx );

                // We will assume that boundary cells are not interfacial
                if (  cell_coordinates[ dimension ] == 0 || 
                        cell_coordinates[ dimension ] == _velocityField[ dimension ].cellDivisions()[ dimension ] - 1 )
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

                if ( _isBulkCell[ pressure_cell_idx1 ] && 
                        _isBulkCell[ pressure_cell_idx2 ] )
                {
                    // Both pressure cell neighbours are bulk cells, so this is
                    // a bulk cell too
                    _isVelocityBulkCell[ dimension ][ cell_idx ] = true;

                    _velocityBulkCells[ dimension ].push_back( cell_idx );
                }
                else if (  _isBulkCell[ pressure_cell_idx1 ] && 
                        !_isBulkCell[ pressure_cell_idx2 ] )
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

                    normal[ dimension ] = 1.0;
                    gradient.normalize();

                    REAL coefficient = normal.dotProduct( gradient );

                    //TRACE_ASSERT( coefficient >= 0.0 );

                    _interfacialBoundaryCoefficients[ dimension ].push_back( coefficient );
                }
                else if ( !_isBulkCell[ pressure_cell_idx1 ] && 
                        _isBulkCell[ pressure_cell_idx2 ] )
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

                //Vector3d closestPoint;

                //const ClosestPointField *sdf = reinterpret_cast<const ClosestPointField*>(_boundaryFields[boundaryObject]); 
                //sdf->closestPoint( position, closestPoint ); 
                //COUT_SDUMP(position-closestPoint); 
                //COUT_SDUMP((position-closestPoint).length()); 

                //normal[ dimension ] = -1.0;
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


    // debug
    //for (size_t ii=0; ii<_ghostCells.size(); ii++) 
    //{
    //    const Vector3d cellPosition = _pressureField.cellPosition(_pressureField.cellIndex(_ghostCells[ii])); 
    //    COUT_SDUMP(ii); 
    //    COUT_SDUMP(cellPosition);
    //}

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


    visualizeClassifiedCells(); 
    ComputeGhostCellInverseMap(); 

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
        position = _pressureField.cellPosition(_pressureField.cellIndex(_bulkCells[bb])); 
        of2 << position.x << ", " << position.y << ", " << position.z << std::endl; 
    }


    of.close(); 
    of2.close();
}

//////////////////////////////////////////////////////////////////////
// Returns the absorption coefficient along a certain
// dimension for a point in space.
//
// FIXME: For now, we will use a quadratic profile here, though
// we may need to try something more complex later on.
//
// if cornellBoxBoundary is on, the PML has been disabled except for the face at -z
// , because we want to make a acoustic cornell box (3/1, 2016)
//////////////////////////////////////////////////////////////////////
REAL MAC_Grid::PML_absorptionCoefficient( const Vector3d &x, REAL absorptionWidth,
        int dimension )
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

#if 0
        cout << "returning "
            << _PML_absorptionStrength * ( absorptionWidth - h )
            / absorptionWidth << endl;
#endif
        //return _PML_absorptionStrength * ( absorptionWidth - h ) / absorptionWidth;
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









