//////////////////////////////////////////////////////////////////////
// Laplacian.cpp: Implementation of the Laplacian class
//
//////////////////////////////////////////////////////////////////////

#include "Laplacian.h"

#include <utils/utils_IO.h>
#include <utils/trace.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor
//
// Provide the size of the domain, finite difference division
// size, and a signed distance function for the interior boundary
//////////////////////////////////////////////////////////////////////
Laplacian::Laplacian( const BoundingBox &bbox, REAL cellSize,
        const TriMesh &mesh,
        const DistanceField &boundarySDF,
        REAL distanceTolerance, int N )
: _field( bbox, cellSize ),
    _distanceTolerance( distanceTolerance ),
    _N( N )
{
    _boundaryFields.push_back( &boundarySDF );
    _boundaryMeshes.push_back( &mesh );

    printf( "Initialized Laplacian with %d cells\n", _field.numCells() );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
Laplacian::Laplacian( const BoundingBox &bbox, REAL cellSize,
        vector<const TriMesh *> &meshes,
        vector<const DistanceField *> &boundaryFields,
        REAL distanceTolerance, int N )
: _field( bbox, cellSize ),
    _boundaryMeshes( meshes ),
    _boundaryFields( boundaryFields ),
    _N( N )
{
    printf( "Initialized Laplacian with %d cells\n", _field.numCells() );

    for ( int mesh_idx = 0; mesh_idx < _boundaryMeshes.size(); mesh_idx++ )
    {
        printf( "Boundary mesh %d has %d vertices and %d triangles\n",
                mesh_idx, (int)_boundaryMeshes[ mesh_idx ]->vertices().size(),
                (int)_boundaryMeshes[ mesh_idx ]->triangles().size() );
    }
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
Laplacian::~Laplacian()
{
}

//////////////////////////////////////////////////////////////////////
// Initializes the field boundary condition, and produces
// a Laplacian discretization
//////////////////////////////////////////////////////////////////////
void Laplacian::initField( bool useBoundary )
{
    cout << "Laplacian::initField: Classifying cells" << endl << endl;
    classifyCells( useBoundary );

    cout << "Laplacian::initField: Building ghost cell coefficients" << endl;
    buildGhostCellCoefficients();

    cout << "Laplacian::initField: Building interfacial coefficients" << endl;
    buildInterfacialCoefficients();

#if 0
    SPARSE_MATRIX              L;

    constructLaplacianMatrix( L );

    L.writeToBinary( "laplacian.bcsm" );
#endif
}

//////////////////////////////////////////////////////////////////////
// Initializes field using a rasterized version of the boundary
//////////////////////////////////////////////////////////////////////
void Laplacian::initFieldRasterized( bool useBoundary )
{
    cout << "Laplacian::initFieldRasterized: Classifying cells" << endl << endl;
    classifyCells( useBoundary );

    cout << "Laplacian::initFieldRasterized: "
        "Building interfacial coefficients" << endl;
    buildInterfacialCoefficientsRasterized();
}

//////////////////////////////////////////////////////////////////////
// Applies the laplacian to field p, and places the result in Lp
//////////////////////////////////////////////////////////////////////
void Laplacian::apply( const MATRIX &p, MATRIX &Lp, REAL alpha ) const
{
    const Tuple3i             &cellDivisions = _field.cellDivisions();
    REAL                       cellSize2;

    cellSize2 = _field.cellSize();
    cellSize2 *= cellSize2;

    if ( p.cols() != _N || Lp.cols() != _N )
    {
        cout << SDUMP( p.cols() ) << endl;
        cout << SDUMP( Lp.cols() ) << endl;
        cout << SDUMP( _N ) << endl;
        cout << SDUMP( _N ) << endl;
        TRACE_ASSERT( p.cols() == _N && Lp.cols() == _N );
    }

    // Handle all bulk cells
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for ( int bulk_cell_idx = 0; bulk_cell_idx < _bulkCells.size();
            bulk_cell_idx++ )
    {
        int                      cell_idx = _bulkCells[ bulk_cell_idx ];
        int                      neighbour_idx;
        Tuple3i                  cell_indices = _field.cellIndex( cell_idx );

        REAL                     coefficient = 6.0;

        if ( cell_indices[ 0 ] == 0 || cell_indices[ 0 ] == cellDivisions[ 0 ] - 1 )
        {
            coefficient -= 1.0;
        }
        if ( cell_indices[ 1 ] == 0 || cell_indices[ 1 ] == cellDivisions[ 1 ] - 1 )
        {
            coefficient -= 1.0;
        }
        if ( cell_indices[ 2 ] == 0 || cell_indices[ 2 ] == cellDivisions[ 2 ] - 1 )
        {
            coefficient -= 1.0;
        }

        //Lp( cell_idx ) -= alpha * 6.0 * p( cell_idx ) / cellSize2;
        //Lp( cell_idx ) -= alpha * coefficient * p( cell_idx ) / cellSize2;
#if 0
        MATRIX::axpy( Lp.data() + cell_idx * _N, /* Desired row of Lp */
                p.data() + cell_idx * _N, /* Desired row of p */
                1, _N, /* Copy size */
                -1.0 * alpha * coefficient / cellSize2 /* multiplier */ );
#endif
        for ( int i = 0; i < _N; i++ )
        {
            Lp( cell_idx, i ) = -alpha * coefficient * p( cell_idx, i ) / cellSize2;
        }

        for ( int direction = GRID_X; direction < NUM_DIRECTIONS; direction++ )
        {
            cell_indices[ direction ] += 1;

            if ( cell_indices[ direction ] < cellDivisions[ direction ] )
            {
                neighbour_idx = _field.cellIndex( cell_indices );

                //Lp( cell_idx ) += alpha * p( neighbour_idx ) / cellSize2;
#if 0
                MATRIX::axpy( Lp.data() + cell_idx * _N, /* Row of Lp */
                        p.data() + neighbour_idx * _N, /* Row of p */
                        1, _N, /* Copy size */
                        alpha / cellSize2 /* multiplier */ );
#endif
                for ( int i = 0; i < _N; i++ )
                {
                    Lp( cell_idx, i ) += alpha * p( neighbour_idx, i ) / cellSize2;
                }
            }

            cell_indices[ direction ] -= 2;

            if ( cell_indices[ direction ] >= 0 )
            {
                neighbour_idx = _field.cellIndex( cell_indices );

                //Lp( cell_idx ) += alpha * p( neighbour_idx ) / cellSize2;
#if 0
                MATRIX::axpy( Lp.data() + cell_idx * _N, /* Row of Lp */
                        p.data() + neighbour_idx * _N, /* Row of p */
                        1, _N, /* Copy size */
                        alpha / cellSize2 );
#endif
                for ( int i = 0; i < _N; i++ )
                {
                    Lp( cell_idx, i ) += alpha * p( neighbour_idx, i ) / cellSize2;
                }
            }

            cell_indices[ direction ] += 1;
        }
    }

#if 0
    // FIXME
    REAL averageCoef = 0.0;
    REAL minCoef = 0.0;
    REAL maxCoef = 0.0;

    REAL averageOffDiagCoef = 0.0;
    long int numOffDiagCoefs = 0;
#endif

    // Handle all interfacial cells
#ifdef USE_OPENMP
#pragma omp parallel for schedule(static) default(shared)
#endif
    for ( int interfacial_cell_idx = 0;
            interfacial_cell_idx < _interfacialCells.size();
            interfacial_cell_idx++ )
    {
        int                      cell_idx;
        int                      neighbour_idx;

        cell_idx = _interfacialCells[ interfacial_cell_idx ];

        const KeyValueArray     &coefficients
            = _interfacialCoefficients[ interfacial_cell_idx ];

        for ( int i = 0; i < _N; i++ )
        {
            Lp( cell_idx, i ) = 0.0;
        }

        for ( int coef_idx = 0; coef_idx < coefficients.size(); coef_idx++ )
        {
            const KeyValuePair    &coef = coefficients[ coef_idx ];

            TRACE_ASSERT(
                    _isBulkCell[ coef.first ] || _isInterfacialCell[ coef.first ],
                    "Invalid application cell" );

            //Lp( cell_idx ) += alpha * coef.second * p( coef.first );
#if 0
            MATRIX::axpy( Lp.data() + cell_idx * _N, /* Row of Lp */
                    p.data() + coef.first * _N, /* Row of p */
                    1, _N, /* Copy size */
                    alpha * coef.second /* multiplier */ );
#endif
            for ( int i = 0; i < _N; i++ )
            {
                Lp( cell_idx, i ) += alpha * coef.second * p( coef.first, i );
            }

#if 0
            // FIXME
            if ( coef.first == cell_idx )
            {
                averageCoef += coef.second;
                minCoef = min( coef.second, minCoef );
                maxCoef = -max( -coef.second, maxCoef );
            }
            else
            {
                averageOffDiagCoef += coef.second;
                numOffDiagCoefs += 1;
            }
#endif
        }
    }

#if 0
    averageCoef /= (REAL)( _interfacialCells.size() );
    averageOffDiagCoef /= (REAL)numOffDiagCoefs;
#endif

#if 0
    printf( "Bulk cell coefficient (%f, %f)\n"
            "Interfacial cell coefficient (%f, %f, %f, %f)\n",
            -6.0 / cellSize2, 1.0 / cellSize2,
            averageCoef, minCoef, maxCoef, averageOffDiagCoef );
#endif
}

//////////////////////////////////////////////////////////////////////
// Forms the boundary condition contribution to Lp resulting from
// the application of the Laplacian
//////////////////////////////////////////////////////////////////////
void Laplacian::applyBoundary( MATRIX &Lp, const BoundaryEvaluator &bc,
        REAL t, REAL alpha )
{
#if 0
    TRACE_ASSERT( _interfacialCells.size() > 0, "No interfacial cells!" );
#endif

    bool anyNonZero = false;

    // This only applies to interfacial cells; bulk cells have no boundary
    // contribution
    for ( int interfacial_cell_idx = 0;
            interfacial_cell_idx < _interfacialCells.size();
            interfacial_cell_idx++ )
    {
        int                      cell_idx;

        cell_idx = _interfacialCells[ interfacial_cell_idx ];

        const BoundaryPointArray &coefficients
            = _interfacialBoundaryCoefficients[ interfacial_cell_idx ];

        TRACE_ASSERT( coefficients.size() > 0, "Empty coefficient list" );

        for ( int coef_idx = 0; coef_idx < coefficients.size(); coef_idx++ )
        {
            const BoundaryPoint   &point = coefficients[ coef_idx ];

            for ( int field_id = 0; field_id < _N; field_id++ )
            {
                // Get a boundary condition for each field
                REAL bcEval = bc( point._x, point._n,
                        _interfacialBoundaryIDs[ interfacial_cell_idx ],
                        t, field_id );

#if 0
                // FIXME
                if ( interfacial_cell_idx == 10 )
                {
                    cout << "   Boundary point " << SDUMP( point._x );
                    cout << SDUMP( point._n ) << SDUMP( bcEval ) << endl;
                }

                if ( bcEval != 0.0 )
                {
                    anyNonZero = true;
                }
#endif

#if 0
                if ( alpha * point._boundaryCoefficient * bcEval > 0 )
                {
                    cout << "A non-zero coefficient" << endl;
                    cout << "    " << alpha * point._boundaryCoefficient * bcEval << endl;
                }
#endif

                Lp( cell_idx, field_id ) += alpha * point._boundaryCoefficient * bcEval;
            }
        }
    }

#if 0
    if ( !anyNonZero )
    {
        cout << "Zero boundary condition" << endl;
    }
#endif
}

#if 0
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void Laplacian::constructLaplacianMatrix( SPARSE_MATRIX &L )
{
    const Tuple3i             &cellDivisions = _field.cellDivisions();
    REAL                       cellSize2;

    int                        nCells = _field.numCells();

    cellSize2 = _field.cellSize();
    cellSize2 *= cellSize2;

    L.resize( _field.numCells(), _field.numCells() );
    L.clear();

    // Handle all bulk cells
    for ( int bulk_cell_idx = 0; bulk_cell_idx < _bulkCells.size();
            bulk_cell_idx++ )
    {
        int                      cell_idx = _bulkCells[ bulk_cell_idx ];
        int                      neighbour_idx;
        Tuple3i                  cell_indices = _field.cellIndex( cell_idx );

        L( cell_idx, cell_idx ) -= 6.0 / cellSize2;

        for ( int direction = GRID_X; direction < NUM_DIRECTIONS; direction++ )
        {
            cell_indices[ direction ] += 1;

            if ( cell_indices[ direction ] < cellDivisions[ direction ] )
            {
                neighbour_idx = _field.cellIndex( cell_indices );

                TRACE_ASSERT( neighbour_idx >= 0 && neighbour_idx < nCells,
                        "Index out of range" );

                L( cell_idx, neighbour_idx ) += 1.0 / cellSize2;
            }

            cell_indices[ direction ] -= 2;

            if ( cell_indices[ direction ] >= 0 )
            {
                neighbour_idx = _field.cellIndex( cell_indices );

                TRACE_ASSERT( neighbour_idx >= 0 && neighbour_idx < nCells,
                        "Index out of range" );

                L( cell_idx, neighbour_idx ) += 1.0 / cellSize2;
            }

            cell_indices[ direction ] += 1;
        }
    }

    // Handle all interfacial cells
    for ( int interfacial_cell_idx = 0;
            interfacial_cell_idx < _interfacialCells.size();
            interfacial_cell_idx++ )
    {
        int                      cell_idx;
        int                      neighbour_idx;

        cell_idx = _interfacialCells[ interfacial_cell_idx ];

        const KeyValueArray     &coefficients
            = _interfacialCoefficients[ interfacial_cell_idx ];

        for ( int coef_idx = 0; coef_idx < coefficients.size(); coef_idx++ )
        {
            const KeyValuePair    &coef = coefficients[ coef_idx ];

            TRACE_ASSERT( coef.first >= 0 && coef.first < nCells,
                    "Index out of range" );

            L( cell_idx, coef.first ) += coef.second;
        }
    }
}
#endif

//////////////////////////////////////////////////////////////////////
// Classifies cells as either a bulk cell, ghost cell, or
// interfacial cell
//////////////////////////////////////////////////////////////////////
void Laplacian::classifyCells( bool useBoundary )
{
    int                        numCells = _field.numCells();
    Vector3d                   cellPos;
    IntArray                   neighbours;

    neighbours.reserve( ScalarField::NUM_NEIGHBOURS );

    _isBulkCell.clear();
    _isInterfacialCell.clear();
    _isGhostCell.clear();

    _isBulkCell.resize( numCells, false );
    _isInterfacialCell.resize( numCells, false );
    _isGhostCell.resize( numCells, false );

    _bulkCells.clear();
    _interfacialCells.clear();
    _ghostCells.clear();

    _ghostCoefficientPointers.clear();
    _ghostCoefficientPointers.resize( numCells, -1 );

    IntArray                   containingObject( numCells, -1 );

    if ( !useBoundary )
    {
        for ( int cell_idx = 0; cell_idx < numCells; cell_idx++ )
        {
            _isBulkCell[ cell_idx ] = true;
            _bulkCells.push_back( cell_idx );
        }

        return;
    }

    // Start by classifying all cells in the computational domain
    // as bulk cells.  This will change once we figure out which
    // cells are ghost cells.
    for ( int cell_idx = 0; cell_idx < numCells; cell_idx++ )
    {
        cellPos = _field.cellPosition( cell_idx );

        // Check all boundary fields to see if this is a bulk cell
        _isBulkCell[ cell_idx ] = true;
        for ( int field_idx = 0; field_idx < _boundaryFields.size(); field_idx++ )
        {
            if ( _boundaryFields[ field_idx ]->distance( cellPos )
                    <= _distanceTolerance )
            {
                _isBulkCell[ cell_idx ] = false;
                containingObject[ cell_idx ] = field_idx;
                break;
            }
        }
#if 0
        if ( _boundaryFields[ 0 ]->distance( cellPos ) > _distanceTolerance )
        {
            _isBulkCell[ cell_idx ] = true;
        }
#endif
    }

    // Ghost cells are cells inside the boundary that have at least
    // one neighbour outside of the boundary
    for ( int cell_idx = 0; cell_idx < numCells; cell_idx++ )
    {
        cellPos = _field.cellPosition( cell_idx );

#if 0
        if ( _boundaryFields[ 0 ]->distance( cellPos ) > _distanceTolerance )
        {
            // Not inside the interior object
            continue;
        }
#endif
        if ( _isBulkCell[ cell_idx ] )
        {
            // Not inside the interior object
            continue;
        }

        _field.cellNeighbours( cell_idx, neighbours );

        for ( int neighbour_idx = 0; neighbour_idx < neighbours.size();
                neighbour_idx++ )
        {
            if ( _isBulkCell[ neighbours[ neighbour_idx ] ] )
            {
                // We have a neighbour outside of the interior object, so
                // this is a ghost cell
                _isGhostCell[ cell_idx ] = true;

                _ghostCells.push_back( cell_idx );

                break;
            }
        }
    }

    // Finally, any cell outside of the interior object, but adjacent
    // to a ghost cell is classified as an interfacial cell
    for ( int cell_idx = 0; cell_idx < numCells; cell_idx++ )
    {
        if ( !_isBulkCell[ cell_idx ] )
        {
            continue;
        }

        _field.cellNeighbours( cell_idx, neighbours );

        for ( int neighbour_idx = 0; neighbour_idx < neighbours.size();
                neighbour_idx++ )
        {
            if ( _isGhostCell[ neighbours[ neighbour_idx ] ] )
            {
                // This cell is adjacent to a ghost cell, so it is actually
                // an interfacial cell and not a bulk cell
                _isInterfacialCell[ cell_idx ] = true;
                _isBulkCell[ cell_idx ] = false;

                _interfacialCells.push_back( cell_idx );

                int boundaryObject = containingObject[ neighbours[ neighbour_idx ] ];
                TRACE_ASSERT( boundaryObject >= 0 );

                _interfacialBoundaryIDs.push_back( boundaryObject );

#if 0
                cout << "Classifying cell at position "
                    << _field.cellPosition( cell_idx )
                    << " as an interfacial cell with boundary ID "
                    << boundaryObject << endl;
#endif

                break;
            }
        }

        // If this is still marked as a bulk cell, add it to the list
        if ( _isBulkCell[ cell_idx ] )
        {
            _bulkCells.push_back( cell_idx );
        }
    }

    // FIXME: Removed this for now, since it is not necessary
    // for the rasterized computation
#if 0
    checkGhostCellValidity();
#endif

    printf( "Laplacian: classifyCells:\n" );
    printf( "\tFound %d bulk cells\n", (int)_bulkCells.size() );
    printf( "\tFound %d interfacial cells\n", (int)_interfacialCells.size() );
    printf( "\tFound %d ghost cells\n\n", (int)_ghostCells.size() );
}

//////////////////////////////////////////////////////////////////////
// Checks the validity of ghost cells.  A ghost cell is only
// valid if we can extrapolate along its gradient to the boundary
// and then out a distance of 2 * h, and lie in a region governed
// only by ghost, bulk or interfacial cells
//
// If we find any invalid cells, we will add new ghost cells
//////////////////////////////////////////////////////////////////////
void Laplacian::checkGhostCellValidity()
{
    Vector3d                   cellPos;
    REAL                       sdfDistance;
    Vector3d                   normal;
    Vector3d                   x;
    KeyValueArray              interpolationCoefficients;

    for ( int ghost_cell_idx = 0; ghost_cell_idx < _ghostCells.size();
            ghost_cell_idx++ )
    {
        cellPos = _field.cellPosition( _ghostCells[ ghost_cell_idx ] );

        sdfDistance = abs( _boundaryFields[ 0 ]->distance( cellPos ) );

        normal = _boundaryFields[ 0 ]->gradient( cellPos );
        normal.normalize();

        x = cellPos + normal * sdfDistance;

        for ( int i = 0; i < 2; i++ )
        {
            x += normal * _field.cellSize();

            _field.interpolationCoefficients( x, interpolationCoefficients );

            for ( int j = 0; j < interpolationCoefficients.size(); j++ )
            {
                if ( !insideDomain( interpolationCoefficients[ j ].first ) )
                {
                    // Add this as a ghost cell
                    _isGhostCell[ interpolationCoefficients[ j ].first ] = true;
                    _ghostCells.push_back( interpolationCoefficients[ j ].first );
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Builds all ghost cell coefficients
//////////////////////////////////////////////////////////////////////
void Laplacian::buildGhostCellCoefficients()
{
    map<int,int>               ghostIndices;
    MATRIX                     ghostCellSystem;
    vector<KeyValueArray>      rhsCoefficients( _ghostCells.size() );
    vector<BoundaryPointArray> boundaryCoefficients( _ghostCells.size() );
    vector<map<int,REAL> >     finalRHSCoefficients( _ghostCells.size() );

    vector<BoundaryPointArray> finalBoundaryCoefficients( _ghostCells.size() );

    ghostCellSystem.resizeAndWipe( _ghostCells.size(), _ghostCells.size() );

    // Start by building the system rows and RHS contributions
    for ( int cell_idx = 0; cell_idx < _field.numCells(); cell_idx++ )
    {
        if ( _isGhostCell[ cell_idx ] )
        {
            buildGhostCellCoefficients( cell_idx );
        }
    }

    TRACE_ASSERT( _ghostCoefficients.size() == _ghostCells.size(),
            "Size mismatch" );

    // Map ghost cell indices to system indices
    for ( int cell_idx = 0; cell_idx < _ghostCells.size(); cell_idx++ )
    {
        ghostIndices[ _ghostCells[ cell_idx ] ] = cell_idx;
    }

    // Go through each ghost cell coefficient list, add contributions
    // from ghost cells to the system, and contributions from non-ghost
    // cells to the right-hand side
    for ( int ghost_cell_idx = 0; ghost_cell_idx < _ghostCells.size();
            ghost_cell_idx++ )
    {
        const KeyValueArray     &coefficients = _ghostCoefficients[ ghost_cell_idx ];

        for ( int coefficient_idx = 0; coefficient_idx < coefficients.size();
                coefficient_idx++ )
        {
            int                    cell_idx = coefficients[ coefficient_idx ].first;
            REAL                   entry = coefficients[ coefficient_idx ].second;

            if ( _isGhostCell[ cell_idx ] )
            {
                int                  entry_idx = ghostIndices[ cell_idx ];

                ghostCellSystem( ghost_cell_idx, entry_idx ) += entry;
            }
            else
            {
                // Subtract this from the right hand side
                rhsCoefficients[ ghost_cell_idx ].push_back(
                        KeyValuePair( cell_idx, -1.0 * entry ) );
            }
        }

        const BoundaryPoint &ghostBoundary
            = _ghostBoundaryCoefficients[ ghost_cell_idx ].back();

        boundaryCoefficients[ ghost_cell_idx ].push_back(
                BoundaryPoint( ghostBoundary._x, ghostBoundary._n,
                    -1.0 * ghostBoundary._boundaryCoefficient ) );
    }

    // FIXME
    ghostCellSystem.write( "ghostSystem.matrix" );

    // Solve for coefficients
    ghostCellSystem.invert();

    // FIXME
    ghostCellSystem.write( "ghostInverse.matrix" );

    int nnz = 0;

    // Build a total list of coefficients
    for ( int row_idx = 0; row_idx < _ghostCells.size(); row_idx++ )
    {
        map<int,REAL>           &cellTerms = finalRHSCoefficients[ row_idx ];
        BoundaryPointArray      &bcTerms = finalBoundaryCoefficients[ row_idx ];

        for ( int col_idx = 0; col_idx < _ghostCells.size(); col_idx++ )
        {
            if ( ghostCellSystem( row_idx, col_idx ) == 0.0 )
            {
                continue;
            }

            const KeyValueArray       &rhsTerms = rhsCoefficients[ col_idx ];
            const BoundaryPointArray  &rhsBCterms = boundaryCoefficients[ col_idx ];

            REAL                       entry = ghostCellSystem( row_idx, col_idx );

            nnz += 1;

            // Add all terms from this entry on the right hand side of the
            // system, scaled by the inverse system entry
            for ( int term_idx = 0; term_idx < rhsTerms.size(); term_idx++ )
            {
                const KeyValuePair    &coef = rhsTerms[ term_idx ];

                cellTerms[ coef.first ] += entry * coef.second;
            }

            // Do the same for boundary terms.  These should be unique though,
            // so we shouldn't need a map.
            for ( int term_idx = 0; term_idx < rhsBCterms.size(); term_idx++ )
            {
                const BoundaryPoint   &bcTerm = rhsBCterms[ term_idx ];

                bcTerms.push_back( BoundaryPoint(
                            bcTerm._x, bcTerm._n, entry * bcTerm._boundaryCoefficient ) );
            }
        }
    }

    printf( "Found %d ghost cell relationships\n", nnz );

    // Assemble the final cell/boundary coefficients for ghost cells
    for ( int ghost_cell_idx = 0; ghost_cell_idx < _ghostCells.size();
            ghost_cell_idx++ )
    {
        KeyValueArray       &coefs = _ghostCoefficients[ ghost_cell_idx ];
        BoundaryPointArray  &bcCoefs = _ghostBoundaryCoefficients[ ghost_cell_idx ];

        map<int,REAL>       &coefMap = finalRHSCoefficients[ ghost_cell_idx ];

        coefs.clear();
        bcCoefs.clear();

        bcCoefs = finalBoundaryCoefficients[ ghost_cell_idx ];

        for ( map<int,REAL>::const_iterator i = coefMap.begin();
                i != coefMap.end(); i++ )
        {
            coefs.push_back( KeyValuePair( i->first, i->second ) );
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Builds the coefficient list for a single ghost cell.
//
// NOTE: Currently this just uses the distance field and
// not the mesh itself.
//////////////////////////////////////////////////////////////////////
void Laplacian::buildGhostCellCoefficients( int cell_idx )
{
    Vector3d                   cellPos = _field.cellPosition( cell_idx );
    REAL                       sdfDistance
        = _boundaryFields[ 0 ]->distance( cellPos );
    Vector3d                   normal = _boundaryFields[ 0 ]->gradient( cellPos );
    map<int,REAL>              coefficients;
    Vector3d                   x1, x2, x3;
    KeyValueArray              interpolationCoefficients_x2;
    KeyValueArray              interpolationCoefficients_x3;
    REAL                       boundaryCoefficient;
    REAL                       x_coefficient;
    REAL                       x2_coefficient;
    REAL                       x3_coefficient;

    // Distance values from [Marelli et al.]
    REAL                       d1 = abs( sdfDistance );
    REAL                       d2 = d1 + _field.cellSize();
    REAL                       d3 = d2 + _field.cellSize();

    TRACE_ASSERT( _isGhostCell[ cell_idx ], "Not a ghost cell" );

    normal.normalize();

    coefficients.clear();

    // x1 is just the estimated surface position.  x2 and x3
    // are offset from this by h and 2h, respectively
    x1 = cellPos + normal * d1;
    x2 = x1 + normal * _field.cellSize();
    x3 = x2 + normal * _field.cellSize();

    // Evaluate x2 and x3 in terms of grid cells
    _field.interpolationCoefficients( x2, interpolationCoefficients_x2 );
    _field.interpolationCoefficients( x3, interpolationCoefficients_x3 );

    // Evaluate coefficients for the cell position itself, x2, x3 and
    // the boundary contribution
    x_coefficient = ghost_pG_Coefficient( d1, d2 );
    x2_coefficient = ghost_pG2_Coefficient( d1, d2 );
    x3_coefficient = ghost_pG3_Coefficient( d1, d2 );
    boundaryCoefficient = ghostBoundaryCoefficient( d1, d2, _field.cellSize() );

    coefficients[ cell_idx ] = x_coefficient;

    TRACE_ASSERT( interpolationCoefficients_x2.size() == 8
            && interpolationCoefficients_x3.size() == 8,
            "Wrong size" );

    // Iterate over the coefficients for the interpolated values of
    // pG2 and pG3 to build a final coefficient list for the ghost cell
    // position.
    //
    // This builds a row/rhs contribution in the ghost cell system
    for ( int i = 0; i < interpolationCoefficients_x2.size(); i++ )
    {
        KeyValuePair            &coef = interpolationCoefficients_x2[ i ];
        int                      interp_idx = coef.first;

        // FIXME
        TRACE_ASSERT( coef.second >= -1e-14 && coef.second <= 1.00000000000001,
                "Invalid interpolation coefficient" );

        if ( interp_idx == cell_idx )
        {
            cout << "Found current ghost cell" << endl;
        }
        else if ( _isGhostCell[ interp_idx ] )
        {
            cout << "Found other ghost cell" << endl;
        }

        if ( !insideDomain( interp_idx ) )
        {
            printf( "Cell %d is not valid.  Its"
                    " distance value is %f\n", interp_idx,
                    _boundaryFields[ 0 ]->distance(
                        _field.cellPosition( interp_idx ) ) );
            printf( "Cell spacing is %f\n", _field.cellSize() );
            printf( "Ghost cell distance is %f\n",
                    _boundaryFields[ 0 ]->distance(
                        _field.cellPosition( cell_idx ) ) );

            Vector3d badCellPos = _field.cellPosition( interp_idx );
            Vector3d direction = _boundaryFields[ 0 ]->gradient( badCellPos );
            direction.normalize();
            REAL badCellDistance = _boundaryFields[ 0 ]->distance( badCellPos );

            badCellPos += abs( badCellDistance ) * direction;
            cout << SDUMP( badCellPos ) << "   ";
            cout << SDUMP( _boundaryFields[ 0 ]->distance( badCellPos ) ) << endl;

            badCellPos += _field.cellSize() * direction;
            cout << SDUMP( badCellPos ) << "   ";
            cout << SDUMP( _boundaryFields[ 0 ]->distance( badCellPos ) ) << endl;
        }

        coefficients[ interp_idx ] += x2_coefficient * coef.second;
    }

    for ( int i = 0; i < interpolationCoefficients_x3.size(); i++ )
    {
        KeyValuePair            &coef = interpolationCoefficients_x3[ i ];
        int                      interp_idx = coef.first;

        coefficients[ interp_idx ] += x3_coefficient * coef.second;
    }

    // Build the final coefficient list
#if 0
    _ghostCoefficientPointers.push_back( (int)_ghostCoefficients.size() );
#endif
    _ghostCoefficientPointers.at( cell_idx ) = (int)_ghostCoefficients.size();
    _ghostCoefficients.push_back( KeyValueArray() );

    _ghostBoundaryCoefficients.push_back( BoundaryPointArray() );

    KeyValueArray       &ghostCoefficients = _ghostCoefficients.back();
    BoundaryPointArray  &boundaryCoefficients = _ghostBoundaryCoefficients.back();

    for ( map<int,REAL>::const_iterator i = coefficients.begin();
            i != coefficients.end(); i++ )
    {
        ghostCoefficients.push_back( KeyValuePair( i->first, i->second ) );
    }

    boundaryCoefficients.push_back( BoundaryPoint( x1, normal,
                boundaryCoefficient ) );
}

//////////////////////////////////////////////////////////////////////
// Given an interfacial point, builds coefficients for the second
// derivative in one direction
//////////////////////////////////////////////////////////////////////
void Laplacian::buildDirectionalCoefficients(
        int cell_idx, GridDirection direction,
        KeyValueArray &coefficients,
        BoundaryPointArray &boundaryCoefficients )
{
    Tuple3i                    neighbourIdx = _field.cellIndex( cell_idx );
    int                        leftNeighbour;
    int                        rightNeighbour;

    KeyValueArray              leftCoefficients;
    KeyValueArray              rightCoefficients;

    BoundaryPointArray         leftBoundaryCoefficients;
    BoundaryPointArray         rightBoundaryCoefficients;

    // Coefficients from [Marelli et al.]
    REAL                       leftRatio;             // xi_{-}
REAL                       rightRatio;            // xi_{+}
REAL                       leftMultiplier;        // alpha_{-}
REAL                       rightMultiplier;       // alpha_{+}
REAL                       totalWeight = 0.0;     // gamma

REAL                       centralCoefficient = 0.0;

const REAL                 MIN_DISTANCE_RATIO = 0.01;

neighbourIdx[ direction ] -= 1;
leftNeighbour = _field.cellIndex( neighbourIdx );

TRACE_ASSERT( neighbourIdx[ direction ] >= 0
        && neighbourIdx[ direction ] < _field.cellDivisions()[ direction ],
        "Neighbour index does not lie inside of the grid" );

neighbourIdx[ direction ] += 2;
rightNeighbour = _field.cellIndex( neighbourIdx );

TRACE_ASSERT( neighbourIdx[ direction ] >= 0
        && neighbourIdx[ direction ] < _field.cellDivisions()[ direction ],
        "Neighbour index does not lie inside of the grid" );

TRACE_ASSERT( _isInterfacialCell[ cell_idx ], "Not an interfacial cell" );

if ( _isBulkCell[ leftNeighbour ] || _isInterfacialCell[ leftNeighbour ] )
{
    // Not much to do here
    leftRatio = 1.0;

    leftCoefficients.push_back( KeyValuePair( leftNeighbour, 1.0 ) );
}
else
{
    leftRatio = buildIntersectionCoefficients( cell_idx, leftNeighbour,
            leftCoefficients,
            leftBoundaryCoefficients );

    // Clamping as described in [Marella et al.]
    leftRatio = max( leftRatio, MIN_DISTANCE_RATIO );
}

leftMultiplier = 1.0 / leftRatio;
totalWeight += leftRatio / 2.0;

if ( _isBulkCell[ rightNeighbour ] || _isInterfacialCell[ rightNeighbour ] )
{
    // Not much to do here
    rightRatio = 1.0;

    rightCoefficients.push_back( KeyValuePair( rightNeighbour, 1.0 ) );
}
else
{
    rightRatio = buildIntersectionCoefficients( cell_idx, rightNeighbour,
            rightCoefficients,
            rightBoundaryCoefficients );

    // Clamping as described in [Marella et al.]
    rightRatio = max( rightRatio, MIN_DISTANCE_RATIO );
}

rightMultiplier = 1.0 / rightRatio;
totalWeight += rightRatio / 2.0;

coefficients.clear();
boundaryCoefficients.clear();

// Add scaled coefficients to the list.  We won't worry about
// guaranteeing uniqueness here, since we can just take care
// of that at the end.
for ( int coef_idx = 0; coef_idx < rightCoefficients.size(); coef_idx++ )
{
    KeyValuePair            &coef = rightCoefficients[ coef_idx ];

    coef.second *= rightMultiplier;
    coef.second /= totalWeight;
    coef.second /= _field.cellSize() * _field.cellSize();

    coefficients.push_back( coef );
}
for ( int coef_idx = 0; coef_idx < leftCoefficients.size(); coef_idx++ )
{
    KeyValuePair            &coef = leftCoefficients[ coef_idx ];

    coef.second *= leftMultiplier;
    coef.second /= totalWeight;
    coef.second /= _field.cellSize() * _field.cellSize();

    coefficients.push_back( coef );
}

// FIXME: Check the sign on the boundary coefficients (we may need
// to scale everything here by -1)
for ( int coef_idx = 0; coef_idx < rightBoundaryCoefficients.size();
        coef_idx++ )
{
    BoundaryPoint           &coef = rightBoundaryCoefficients[ coef_idx ];

    coef._boundaryCoefficient *= rightMultiplier;
    coef._boundaryCoefficient /= totalWeight;
    coef._boundaryCoefficient /= _field.cellSize() * _field.cellSize();

    boundaryCoefficients.push_back( coef );
}
for ( int coef_idx = 0; coef_idx < leftBoundaryCoefficients.size();
        coef_idx++ )
{
    BoundaryPoint           &coef = leftBoundaryCoefficients[ coef_idx ];

    coef._boundaryCoefficient *= leftMultiplier;
    coef._boundaryCoefficient /= totalWeight;
    coef._boundaryCoefficient /= _field.cellSize() * _field.cellSize();

    boundaryCoefficients.push_back( leftBoundaryCoefficients[ coef_idx ] );
}

// Add the contribution for this cell
centralCoefficient = -1.0 * ( leftMultiplier + rightMultiplier );
centralCoefficient /= totalWeight;
centralCoefficient /= _field.cellSize() * _field.cellSize();

coefficients.push_back( KeyValuePair( cell_idx, centralCoefficient ) );
}

//////////////////////////////////////////////////////////////////////
// Builds the coefficients resulting from the line joining cell_idx
// and neighbour_idx and the boundary sdf
//////////////////////////////////////////////////////////////////////
REAL Laplacian::buildIntersectionCoefficients(
        int cell_idx, int neighbour_idx,
        KeyValueArray &coefficients,
        BoundaryPointArray &boundaryCoefficients )
{
    Vector3d                   x;
    Vector3d                   d;

    REAL                       boundaryDistance;
    Vector3d                   intersection;

    REAL                       distanceRatio;

    // Figure out where the line joining these two cells crosses
    // the boundary
    x = _field.cellPosition( cell_idx );
    d = _field.cellPosition( neighbour_idx );
    d -= x;

    boundaryDistance = boundaryIntersection( x, d );

    intersection = x + boundaryDistance * d;

    buildIntersectionCoefficients( intersection, coefficients,
            boundaryCoefficients );

    distanceRatio = ( intersection - x ).length();
    distanceRatio /= _field.cellSize();

    return distanceRatio;
}

//////////////////////////////////////////////////////////////////////
// Given an intersection point with the boundary, build coefficients
// by extrapolating a normal line out from the boundary and constraining
// the point at the boundary to obey a Neumann BC
//////////////////////////////////////////////////////////////////////
void Laplacian::buildIntersectionCoefficients(
        const Vector3d &intersection,
        KeyValueArray &coefficients,
        BoundaryPointArray &boundaryCoefficients )
{
    KeyValueArray              interpolationCoefficients;
    map<int,REAL>              finalCoefficients;
    int                        cell_idx, cell_idx2;
    REAL                       interp_coef, interp_coef2;

    // Points x1 and x2, following equations (38b) and (38c) from
    // [Marelli et al.]
    Vector3d                   x[ 2 ];
    Vector3d                   normal;

    // From equation (39)
    REAL                       coefs[] = { 4.0 / 3.0, -1.0 / 3.0 };
    REAL                       boundaryCoef = -2.0 / 3.0;

    normal = _boundaryFields[ 0 ]->gradient( intersection );
    normal.normalize();

    x[ 0 ] = intersection + _field.cellSize() * normal;
    x[ 1 ] = x[ 0 ] + _field.cellSize() * normal;

    // Add the boundary point at the intersection
    boundaryCoefficients.push_back( BoundaryPoint(
                intersection, normal,
                boundaryCoef * _field.cellSize() ) );

    // Add the coefficients for these nodes
    for ( int i = 0; i < 2; i++ )
    {
        _field.interpolationCoefficients( x[ i ], interpolationCoefficients );

        for ( int interp_idx = 0; interp_idx < interpolationCoefficients.size();
                interp_idx++ )
        {
            cell_idx = interpolationCoefficients[ interp_idx ].first;
            interp_coef = interpolationCoefficients[ interp_idx ].second;

            if ( _isBulkCell[ cell_idx ] || _isInterfacialCell[ cell_idx ] )
            {
                finalCoefficients[ cell_idx ] += interp_coef * coefs[ i ];
            }
            else if ( _isGhostCell[ cell_idx ] )
            {
                const KeyValueArray &ghostCoefficients
                    = _ghostCoefficients.at( _ghostCoefficientPointers.at( cell_idx ) );

                // Add each cell in the ghost coefficient list, scaling their
                // coefficients by the inerpolation coefficient
                for ( int ghost_idx = 0; ghost_idx < ghostCoefficients.size();
                        ghost_idx++ )
                {
                    cell_idx2 = ghostCoefficients[ ghost_idx ].first;
                    interp_coef2 = ghostCoefficients[ ghost_idx ].second;

                    TRACE_ASSERT(
                            _isBulkCell[ cell_idx2 ] || _isInterfacialCell[ cell_idx2 ],
                            "Not a valid cell" );

                    finalCoefficients[ cell_idx2 ]
                        += interp_coef * interp_coef2 * coefs[ i ];
                }

                const BoundaryPointArray &ghostBoundaryCoefficients
                    = _ghostBoundaryCoefficients.at(
                            _ghostCoefficientPointers[ cell_idx ] );

                // Add ghost boundary points to the mix
                for ( int ghost_idx = 0; ghost_idx < ghostBoundaryCoefficients.size();
                        ghost_idx++ )
                {
                    boundaryCoefficients.push_back(
                            ghostBoundaryCoefficients[ ghost_idx ] );

                    // Rescale the coefficient
                    boundaryCoefficients.back()._boundaryCoefficient *=
                        interp_coef * coefs[ i ];
                }
            }
            else
            {
                TRACE_ASSERT( NULL, "Unclassified cell" );
            }
        }
    }

    // Add to the coefficient list
    for ( map<int,REAL>::const_iterator i = finalCoefficients.begin();
            i != finalCoefficients.end(); i++ )
    {
        coefficients.push_back( KeyValuePair( i->first, i->second ) );
    }
}

//////////////////////////////////////////////////////////////////////
// Returns the distance along vector d, originating at x at which
// an intersection with the boundary field occurs.  The
// result is assumed to lie within (0,1).
//
// x is assumed to lie outside of the boundary field, and x + d is
// assumed to lie inside.
//
// Intersection is determined via bisection.
//////////////////////////////////////////////////////////////////////
REAL Laplacian::boundaryIntersection( Vector3d x, Vector3d d, REAL tolerance )
{
    REAL                       distance = d.length();
    REAL                       startDistance = distance;
    Vector3d                   startPoint = x;
    Vector3d                   midPoint;

    TRACE_ASSERT( _boundaryFields[ 0 ]->distance( x ) > 0.0
            && _boundaryFields[ 0 ]->distance( x + d ) < 0.0,
            "Initial bisection condition violated" );

    while ( distance > tolerance * _field.cellSize() )
    {
        // FIXME
        REAL xd = _boundaryFields[ 0 ]->distance( x );
        REAL md = _boundaryFields[ 0 ]->distance( x + d );
        TRACE_ASSERT( _boundaryFields[ 0 ]->distance( x ) > 0.0
                && _boundaryFields[ 0 ]->distance( x + d ) <= 0.0,
                "Bisection condition violated" );

        d /= 2.0;
        midPoint = x + d;

        if ( _boundaryFields[ 0 ]->distance( midPoint ) > 0.0 )
        {
            //x += d / 2.0;
            x = midPoint;
        }

        distance = d.length();
    }

    return ( ( x + d / 2.0 - startPoint ).length() / startDistance );
}

//////////////////////////////////////////////////////////////////////
// Builds coefficients for all interface points
//////////////////////////////////////////////////////////////////////
void Laplacian::buildInterfacialCoefficients()
{
    KeyValueArray              directionCoefficients;
    BoundaryPointArray         directionBoundaryCoefficients;

    _interfacialCoefficients.resize( _interfacialCells.size() );
    _interfacialBoundaryCoefficients.resize( _interfacialCells.size() );

    for ( int cell_idx = 0; cell_idx < _interfacialCells.size(); cell_idx++ )
    {
        // Note: this assumes that no interfacial cells lie on the
        //       grid boundary.
        map<int,REAL>            fullCoefficientSet;

        if ( _field.isBoundaryCell( _interfacialCells[ cell_idx ] ) )
        {
            cout << "Interfacial cell " << _interfacialCells[ cell_idx ]
                << " at position "
                << _field.cellPosition( _interfacialCells[ cell_idx ] )
                << " is on the boundary" << endl;

            TRACE_ASSERT( !_field.isBoundaryCell( _interfacialCells[ cell_idx ] ),
                    "Interfacial cell cannot lie on the domain boundary" );
        }

        KeyValueArray           &coefficients
            = _interfacialCoefficients[ cell_idx ];
        BoundaryPointArray      &boundaryCoefficients
            = _interfacialBoundaryCoefficients[ cell_idx ];

        coefficients.clear();
        boundaryCoefficients.clear();

        // Get coefficients in each direction
        for ( int direction = GRID_X; direction < NUM_DIRECTIONS; direction += 1 )
        {
            buildDirectionalCoefficients( _interfacialCells[ cell_idx ],
                    (GridDirection)direction,
                    directionCoefficients,
                    directionBoundaryCoefficients );

            // Put the cell coefficients in a map so that we can guarantee
            // uniqueness in the final list
            for ( int coef_idx = 0; coef_idx < directionCoefficients.size();
                    coef_idx++ )
            {
                const KeyValuePair  &coef = directionCoefficients[ coef_idx ];

                fullCoefficientSet[ coef.first ] += coef.second;
            }

            // Add the boundary conditions to the list for this cell
            for ( int coef_idx = 0; coef_idx < directionBoundaryCoefficients.size();
                    coef_idx++ )
            {
                const BoundaryPoint &coef = directionBoundaryCoefficients[ coef_idx ];

                boundaryCoefficients.push_back( coef );
            }
        }

        TRACE_ASSERT( boundaryCoefficients.size() > 0,
                "Empty boundary coefficients list " );

        // Now that we have unique coefficients for each point involved,
        // build the full list of coefficients
        for ( map<int,REAL>::const_iterator coef_iter = fullCoefficientSet.begin();
                coef_iter != fullCoefficientSet.end(); coef_iter++ )
        {
            coefficients.push_back( KeyValuePair( coef_iter->first,
                        coef_iter->second ) );
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Builds coefficients for interface points assuming a rasterized
// version of the input shape
//////////////////////////////////////////////////////////////////////
void Laplacian::buildInterfacialCoefficientsRasterized()
{
    REAL                       cellSize, cellSize2;
    Tuple3i                    neighbourIndex;
    int                        neighbourIndexFlat;

    Vector3d                   normal;
    Vector3d                   cellPosition, neighbourPosition;

    cellSize = _field.cellSize();
    cellSize2 = cellSize * cellSize;

    _interfacialCoefficients.resize( _interfacialCells.size() );
    _interfacialBoundaryCoefficients.resize( _interfacialCells.size() );

    for ( int cell_idx = 0; cell_idx < _interfacialCells.size(); cell_idx++ )
    {
        // Note: this assumes that no interfacial cells lie on the
        //       grid boundary.
        map<int,REAL>            fullCoefficientSet;
        Tuple3i                  cellIndex;

        cellPosition = _field.cellPosition( _interfacialCells[ cell_idx ] );
        neighbourIndex = _field.cellIndex( _interfacialCells[ cell_idx ] );

#if 0
        if ( _field.isBoundaryCell( _interfacialCells[ cell_idx ] ) )
        {
            cout << "Interfacial cell " << _interfacialCells[ cell_idx ]
                << " at position "
                << _field.cellPosition( _interfacialCells[ cell_idx ] )
                << " is on the boundary" << endl;

            TRACE_ASSERT( !_field.isBoundaryCell( _interfacialCells[ cell_idx ] ),
                    "Interfacial cell cannot lie on the domain boundary" );
        }
#endif

        KeyValueArray           &coefficients
            = _interfacialCoefficients[ cell_idx ];
        BoundaryPointArray      &boundaryCoefficients
            = _interfacialBoundaryCoefficients[ cell_idx ];

        coefficients.clear();
        boundaryCoefficients.clear();

        // Push the diagonal coefficient
        fullCoefficientSet[ _interfacialCells[ cell_idx ] ] = -6.0 / cellSize2;

        // Get coefficients in each direction
        for ( int direction = GRID_X; direction < NUM_DIRECTIONS; direction += 1 )
        {
            neighbourIndex[ direction ] += 1;

            if ( neighbourIndex[ direction ]
                    >= _field.cellDivisions()[ direction ] - 1 )
            {
                fullCoefficientSet[ _interfacialCells[ cell_idx ] ] += 1.0 / cellSize2;
            }
            else
            {
                neighbourIndexFlat = _field.cellIndex( neighbourIndex );

                if ( _isGhostCell[ neighbourIndexFlat ] )
                {
                    normal = Vector3d( 0.0, 0.0, 0.0 );
                    normal[ direction ] = -1.0;

                    // Get the mid point and assume a grid-aligned boundary
                    // condition exists here
                    neighbourPosition = _field.cellPosition( neighbourIndex );
                    neighbourPosition += cellPosition;
                    neighbourPosition /= 2.0;

                    // Evaluate neighbour cell using boundary condition
                    fullCoefficientSet[ _interfacialCells[ cell_idx ] ] += 1.0 / cellSize2;

                    boundaryCoefficients.push_back( BoundaryPoint(
                                neighbourPosition, normal,
                                -1.0 / cellSize ) );
                }
                else
                {
                    fullCoefficientSet[ neighbourIndexFlat ] += 1.0 / cellSize2;
                }
            }

            neighbourIndex[ direction ] -= 2;

            if ( neighbourIndex[ direction ] < 1 )
            {
                fullCoefficientSet[ _interfacialCells[ cell_idx ] ] += 1.0 / cellSize2;
            }
            else
            {
                neighbourIndexFlat = _field.cellIndex( neighbourIndex );

                if ( _isGhostCell[ neighbourIndexFlat ] )
                {
                    normal = Vector3d( 0.0, 0.0, 0.0 );
                    normal[ direction ] = 1.0;

                    // Get the mid point and assume a grid-aligned boundary
                    // condition exists here
                    neighbourPosition = _field.cellPosition( neighbourIndex );
                    neighbourPosition += cellPosition;
                    neighbourPosition /= 2.0;

                    // Evaluate neighbour cell using boundary condition
                    fullCoefficientSet[ _interfacialCells[ cell_idx ] ] += 1.0 / cellSize2;

                    boundaryCoefficients.push_back( BoundaryPoint(
                                neighbourPosition, normal,
                                -1.0 / cellSize ) );
                }
                else
                {
                    fullCoefficientSet[ neighbourIndexFlat ] += 1.0 / cellSize2;
                }
            }

            neighbourIndex[ direction ] += 1;
        }

        // Now that we have unique coefficients for each point involved,
        // build the full list of coefficients
        for ( map<int,REAL>::const_iterator coef_iter = fullCoefficientSet.begin();
                coef_iter != fullCoefficientSet.end(); coef_iter++ )
        {
            coefficients.push_back( KeyValuePair( coef_iter->first,
                        coef_iter->second ) );
        }
    }
}

//////////////////////////////////////////////////////////////////////
// Ghost cell coefficient calculations.
//
// Based on (41) and (39) (with a non-zero boundary condition) in
// [Marelli et al.]
//////////////////////////////////////////////////////////////////////
REAL Laplacian::ghostBoundaryCoefficient( REAL d1, REAL d2, REAL h )
{
    REAL                       coef;

    coef = d2 * d2 - 2.0 * d1 * d2;
    coef /= ( d2 - d1 ) * ( d2 - d1 );
    coef *= h;
    coef *= 2.0 / 3.0;

    return coef;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
REAL Laplacian::ghost_pG_Coefficient( REAL d1, REAL d2 )
{
    return 1.0;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
REAL Laplacian::ghost_pG2_Coefficient( REAL d1, REAL d2 )
{
    REAL                       coef;

    coef = d2 * d2 - 2.0 * d1 * d2;
    coef *= 4.0 / 3.0;
    coef += d1 * d1;
    coef *= -1.0;
    coef /= ( d2 - d1 ) * ( d2 - d1 );

    return coef;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
REAL Laplacian::ghost_pG3_Coefficient( REAL d1, REAL d2 )
{
    REAL                       coef;

    coef = d2 * d2 - 2.0 * d1 * d2;
    coef /= ( d2 - d1 ) * ( d2 - d1 );
    coef /= 3.0;

    return coef;
}
