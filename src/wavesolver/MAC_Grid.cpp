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
#include <graph/UndirectedGraph.h> 

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
    const IntArray        &bulkCells = _velocityBulkCells[ dimension ];
    const IntArray        &interfacialCells = _velocityInterfacialCells[ dimension ];
    const ScalarField     &field = _velocityField[ dimension ];
    const IntArray        &interfaceBoundaryIDs = _interfacialBoundaryIDs[ dimension ];
    const FloatArray      &interfaceBoundaryDirections
                                        = _interfacialBoundaryDirections[ dimension ];
    const FloatArray      &interfaceBoundaryCoefficients
                                        = _interfacialBoundaryCoefficients[ dimension ];

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

    // Handle boundary cells.
    //
    // Note: we assume here that all boundary cells have zero absorption
    // coefficient
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




        // attempt to do interpolation/extrapolation on acceleration // 
          


          
        //// TODO debug
        //REAL                 bcEvalOld; 

        //// compute the closest point to the interfacial face center. 
        //// TODO this part can be optimized easily debug
        //Vector3d             closestPoint;
        //const ClosestPointField *csdf = reinterpret_cast<const ClosestPointField*>(_boundaryFields[objectID]); 
        //csdf->closestPoint(x, closestPoint); 

        //IntArray             neighbourIndices; 
        //IntArray             neighbourIndicesTrimmed; 


        //for ( int i = 0; i < _N; i++ )
        //{
        //    bcEvalOld = bc( x, normal, objectID, t, i );
        //    bcEvalOld *= coefficient;

        //      
        //    // evaluate the boundary condition on the closest point -> exact 
        //    bcEval = bc( closestPoint, normal, objectID, t, i );
        //    bcEval *= coefficient; // TODO this part still cast away the tangential component

        //    // search for local neighbors
        //    //field.cell26Neighbours(cell_idx, neighbourIndices); 
        //    field.cellNeighbours(cell_idx, neighbourIndices); 

        //    // first filter out all interfacial ones, only keep the bulk cells
        //    for (size_t jj=0; jj<neighbourIndices.size(); jj++) 
        //    {
        //        if (_isVelocityInterfacialCell[dimension][neighbourIndices[jj]])
        //            continue; 
        //        neighbourIndicesTrimmed.push_back(neighbourIndices[jj]); 
        //    }

        //    // currently cannot handling underdetermined system
        //    if (neighbourIndicesTrimmed.size() < 3) 
        //    {
        //        field.AddCornerNeighbours(cell_idx, neighbourIndices); 
        //        neighbourIndicesTrimmed.clear(); 
        //        for (size_t jj=0; jj<neighbourIndices.size(); jj++) 
        //        {
        //            if (_isVelocityInterfacialCell[dimension][neighbourIndices[jj]])
        //                continue; 
        //            neighbourIndicesTrimmed.push_back(neighbourIndices[jj]); 
        //        }
        //        if (neighbourIndicesTrimmed.size() < 3) 
        //            throw std::runtime_error("**ERROR** have less than 3 bulk neighbours in interfacial interpolation. this will cause underdetermined system"); 
        //    }

        //    const int N_neighbours = neighbourIndicesTrimmed.size(); 
        //    //const int N_neighbours = 5; 
        //    // forming the least-square matrix and sample values. the last row
        //    // will store the boundary velocity samples
        //    Eigen::MatrixXd samplePoints(N_neighbours+1, 3); 
        //    Eigen::VectorXd sampleValues(N_neighbours+1); 

        //    Vector3d neighbourPositionBuffer; 
        //    int      neighbourIndexBuffer; 

        //    // TODO debug
        //    //std::random_shuffle(neighbourIndicesTrimmed.begin(), neighbourIndicesTrimmed.end()); 
        //    for (int jj=0; jj<N_neighbours; jj++) 
        //    //for (int jj=0; jj<3; jj++) 
        //    {

        //        neighbourIndexBuffer = neighbourIndicesTrimmed[jj]; 

        //        neighbourPositionBuffer = field.cellPosition(neighbourIndexBuffer); 

        //        samplePoints(jj,0) = neighbourPositionBuffer.x; 
        //        samplePoints(jj,1) = neighbourPositionBuffer.y; 
        //        samplePoints(jj,2) = neighbourPositionBuffer.z; 

        //        // TODO debug
        //        sampleValues(jj)    = v(neighbourIndexBuffer, i) - v_old(neighbourIndexBuffer, i) / timeStep; 
        //    }

        //    samplePoints(N_neighbours,0) = closestPoint.x; 
        //    samplePoints(N_neighbours,1) = closestPoint.y; 
        //    samplePoints(N_neighbours,2) = closestPoint.z; 
        //    sampleValues(N_neighbours)   = bcEval; 


        //    WeightedLeastSquareSurfaceLinear3D surface; 
        //    //LeastSquareSurfaceLinear3D surface; 
        //    //surface.ComputeCoefficients(samplePoints, sampleValues); 

        //    // TODO god damn debug
        //    Eigen::Vector3d xEigen; 
        //    xEigen(0) = x.x; 
        //    xEigen(1) = x.y; 
        //    xEigen(2) = x.z; 


        //    //bcEval = surface.Evaluate(xEigen); 
        //    
        //    // manually set weights
        //    Eigen::VectorXd weights = Eigen::VectorXd::Ones(samplePoints.rows()); 
        //    weights(N_neighbours) = 5;
        //    surface.SetGaussianSpacing(0.02); 
        //    surface.SetWeights(weights);
        //    //surface.SetInverseQuadraticEps(0.1);
        //    bcEval = surface.Evaluate(samplePoints, sampleValues, xEigen); 



        //    // TODO debug
        //    Eigen::VectorXd difference(samplePoints.rows()); 
        //    for (int jj=0; jj<samplePoints.rows(); jj++)
        //    {
        //        Eigen::Vector3d diff;
        //        diff[0] = samplePoints(jj,0) - xEigen(0); 
        //        diff[1] = samplePoints(jj,1) - xEigen(1); 
        //        diff[2] = samplePoints(jj,2) - xEigen(2); 
        //        
        //        difference(jj) = diff.norm(); 
        //    }



        //    // TODO debug
        //    if (interfacial_cell_idx == 111) 
        //    {

        //        COUT_SDUMP(coefficient);
        //        COUT_SDUMP(x);
        //        COUT_SDUMP(difference.transpose()); 
        //        COUT_SDUMP(sampleValues.transpose()); 
        //        // debug
        //        std::cout << "weights = " << surface.GetWeights().transpose() << std::endl;; 
        //        std::cout << "bcEval: " << bcEval << std::endl; 
        //        std::cout << "bcEvalOld: " << bcEvalOld << std::endl;
        //        std::cout << "difference least square makes: " << bcEval - bcEvalOld << std::endl; 
        //    }



        //    v( cell_idx, i ) += timeStep * bcEval;
        //}
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

void MAC_Grid::PML_pressureUpdateGhostCells( MATRIX &p, const REAL &timeStep, const REAL &c, const BoundaryEvaluator &bc, const REAL &simulationTime, const REAL density)
{

    boost::timer::auto_cpu_timer timer("GC update takes %w sec\n"); 

    std::cout << "pressure update for the ghost cells" << std::endl;



    // use a graph Laplacian to store which ghost cells are coupled
    UndirectedGraph ghostCellCoupling(numPressureCells()); 

    // loop through all ghost cells to get the uncoupled ghost cell values. 
    // construct the graph for the coupled ghost cells 
    for (size_t ghost_cell_idx=0; ghost_cell_idx<_ghostCells.size(); ghost_cell_idx++) 
    {

        const int       cellIndex      = _ghostCells[ghost_cell_idx]; 
        const Tuple3i   cellIndices    = _pressureField.cellIndex(cellIndex); 
        const Vector3d  cellPosition   = _pressureField.cellPosition(cellIndices); 
        const int       boundaryObject = _containingObject[cellIndex]; 

        // find point BI and IP in the formulation
        Vector3d boundaryPoint, imagePoint, erectedNormal; 
        FindImagePoint(cellPosition, boundaryObject, boundaryPoint, imagePoint, erectedNormal); 

        // get the box enclosing the image point; 
        IntArray neighbours; 
        _pressureField.enclosingNeighbours(imagePoint, neighbours); 

        assert(neighbours.size()==8);

        // hasGC  : has self as interpolation stencil
        // coupled: if at least one interpolation stencil is another ghost-cell 
        int hasGC=-1; 
        int coupled=0; 


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
                    coupled += 1;

                    // assume symmetry
                    ghostCellCoupling.AddEdgeToStage(cellIndex, neighbours[ii]); 
                    //tripletList.push_back(Triplet(cellIndex,neighbours[ii],1)); 
                    //tripletList.push_back(Triplet(neighbours[ii],cellIndex,1));
                }
            }
        }

        
        //std::cout << "for ghost cell " << cellIndex << " at position " << cellPosition << ": \n"; 
        //std::cout << "hasGC=" << hasGC << std::endl; 
        //if (coupled != 0)
        //{
        //    std::cout << "for ghost cell " << cellIndex << " at position " << cellPosition << ": \n"; 
        //    std::cout << "coupled=" << coupled << std::endl << std::endl; 
        //}

        // vandermonde matrix, see 2008 Mittals JCP paper Eq.18
        Eigen::MatrixXd V(8,8); 
        Tuple3i     indicesBuffer;
        Vector3d    positionBuffer; 
        Eigen::VectorXd pressureNeighbours(8); 

        // FIXME need to pass in boundary evaluator
        const double bcEval = bc(boundaryPoint, erectedNormal, boundaryObject, simulationTime, 0);  // evaluate at boundarPoint

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
        Eigen::VectorXd beta = V.householderQr().solve(b);  // this gives beta=V^-T b

        if (coupled == 0) // not coupled, can solve now and we are done
        {
            p(cellIndex, 0) = beta.dot(pressureNeighbours) - (cellPosition - imagePoint).length() * bcEval; 
        } 

    }

    // set the off-diagonal term of the laplacian 
    ghostCellCoupling.ConstructEdges(); 

    std::vector<std::vector<int> > allComponents; 
    ghostCellCoupling.GetAllComponentsVerbose(allComponents); 


    Vector3d positionBuffer; 
    const int numberComponents = allComponents.size(); 
    for (int component=0; component<numberComponents; component++) 
    {
        std::ofstream of(("component_"+to_string(component)+".csv").c_str()); 
        for (size_t ii=0; ii<allComponents[component].size(); ii++) 
        {
            positionBuffer = _pressureField.cellPosition(_pressureField.cellIndex(allComponents[component][ii])); 
            of << positionBuffer.x << ", " << positionBuffer.y << ", " << positionBuffer.z << std::endl; 
        }
        of.close();
    }

    // FIXME debug
    std::vector<int> decoupledIndex; 
    for (size_t jj=0; jj<_ghostCells.size(); jj++)
    {

        bool coupled = false; 

        for (int component=0; component<numberComponents; component++) 
        {
            for (size_t ii=0; ii<allComponents[component].size(); ii++) 
            {
                if (_ghostCells[jj] == allComponents[component][ii]) 
                    coupled = true; 
            }
        }

        if (!coupled)
            decoupledIndex.push_back(_ghostCells[jj]); 
    }

    std::ofstream of("decoupled_component.csv"); 
    for (int &p : decoupledIndex)
    {
        positionBuffer = _pressureField.cellPosition(_pressureField.cellIndex(p)); 
        of << positionBuffer.x << ", " << positionBuffer.y << ", " << positionBuffer.z << std::endl; 
    }
    of.close();


    usleep(1000);



    //ghostCellLaplacian.setFromTriplets(tripletList.begin(), tripletList.end()); 

    //COUT_SDUMP(ghostCellLaplacian.cols()); 
    //COUT_SDUMP(ghostCellLaplacian.outerSize()); 

    //// checking the sparse laplacian 
    //{
    //    boost::timer::auto_cpu_timer t("iterate through the matrix columns %w sec\n"); 
    //    for (int col=0; col<ghostCellLaplacian.outerSize(); col++) 
    //    {
    //        //const SparseMatrix colVector = ghostCellLaplacian.block(0,col,numPressureCells(),1); 
    //        //const int nnz = colVector.nonZeros(); 
    //        //if (nnz != 0) 
    //        //{
    //        //    std::cout << "column " << col << " for the sparse matrix has " << nnz << " non-zero entries" << std::endl;
    //        //}

    //        for (SparseMatrix::InnerIterator it(ghostCellLaplacian,col); it; ++it)
    //        {
    //            std::cout << " [" << it.row() << "," << it.col() << "] = " << it.value() << std::endl;
    //        }


    //    }
    //}








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
    for (size_t bb=0; bb<_ghostCells.size(); bb++) 
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
    // FIXME reinterpret
    const ClosestPointField *csdf = reinterpret_cast<const ClosestPointField*>(_boundaryFields[boundaryObjectID]); 
    csdf->closestPoint(cellPosition, closestPoint); 

    Vector3d sdfNormal = csdf->gradient(closestPoint); 
    sdfNormal.normalize();

    erectedNormal = (closestPoint - cellPosition); 
    const double angle = acos(sdfNormal.dotProduct(erectedNormal)/erectedNormal.length())*180./M_PI; 

    if (angle > 10)
        std::cerr << "**WARNING** angle between the erected normal and true normal is greater than 10 degrees, might be inaccurate. angle = " << angle << " for ghost cell at : " << cellPosition.x << "," << cellPosition.y << "," << cellPosition.z << std::endl;


    imagePoint = cellPosition + erectedNormal*2.0; 

    if (csdf->distance(imagePoint) < 0)
        std::cerr << "**WARNING** image point is still inside the object, could be that object is severely underresolved." << std::endl;


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









