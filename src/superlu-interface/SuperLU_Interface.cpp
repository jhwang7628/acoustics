//////////////////////////////////////////////////////////////////////
// SuperLU_Interface.cpp: Implementation of the SuperLU_Interface
//                        class
//
//////////////////////////////////////////////////////////////////////

#include "SuperLU_Interface.h"

#include <utils/trace.h>

#include <iostream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
SuperLU_Solver::SuperLU_Solver( const SPARSE_MATRIX &A )
    : _columnPermutation( NULL ),
      _rowPermutation( NULL ),
      _nzval( NULL ),
      _rowindices( NULL ),
      _colptr( NULL )
{
    SPARSE_MATRIX::SparseColumnMatrix AsparseCol;

    A.constructSparseColumnCopy( AsparseCol );

    init( AsparseCol );
}

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
SuperLU_Solver::SuperLU_Solver( const SPARSE_MATRIX::SparseColumnMatrix &A )
    : _columnPermutation( NULL ),
      _rowPermutation( NULL ),
      _nzval( NULL ),
      _rowindices( NULL ),
      _colptr( NULL )
{
    init( A );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
SuperLU_Solver::~SuperLU_Solver()
{
    SUPERLU_FREE( _rowPermutation );
    SUPERLU_FREE( _columnPermutation );

    // This should also take care of _nzval, etc.
    Destroy_CompCol_Matrix( &_A );
    Destroy_SuperNode_Matrix( &_L );
    Destroy_CompCol_Matrix( &_U );
    StatFree( &_stat );
}

//////////////////////////////////////////////////////////////////////
// Solves the system with the given right hand side
//
// This function assumes that rhs and solution have the correct size
//////////////////////////////////////////////////////////////////////
void SuperLU_Solver::solve( const REAL *rhs, REAL *solution )
{
    SuperMatrix              B;
    int                      info;
    REAL                    *rhsStorage;

    if ( !(rhsStorage = doubleMalloc( _A.nrow )) ) {
        TRACE_ASSERT( NULL, "Error allocating RHS storage" );
    }

    // Copy
    for ( int i = 0; i < _A.nrow; i++ ) {
        rhsStorage[ i ] = rhs[ i ];
    }

    dCreate_Dense_Matrix( &B, _A.nrow, 1 /* Single RHS */, rhsStorage,
                          _A.nrow, /* Screwed up leading dimension due to column-major */
                          SLU_DN, SLU_D, SLU_GE /* Various SuperLU definitions */ );

    if ( _options.Fact == DOFACT ) {
        // The actual driver routine to factor (if necessary) then solve the system
        dgssv( &_options, &_A, _columnPermutation, _rowPermutation, &_L, &_U, &B, &_stat, &info );
    } else {
        dgstrs( NOTRANS, &_L, &_U, _columnPermutation, _rowPermutation, &B, &_stat, &info );
    }

    // Copy back to the solution
    for ( int i = 0; i < _A.nrow; i++ ) {
        solution[ i ] = rhsStorage[ i ];
    }

#if 0
    dPrint_CompCol_Matrix( "A", &_A );
    dPrint_CompCol_Matrix( "U", &_U );
    dPrint_SuperNode_Matrix( "L", &_L );

    printf( "Row permutation: \n" );
    for ( int i = 0; i < _A.nrow; i++ ) {
        printf( "  %d --> %d\n", i, _rowPermutation[ i ] );
    }
    printf( "\n" );

    printf( "Column permutation: \n" );
    for ( int i = 0; i < _A.ncol; i++ ) {
        printf( "  %d --> %d\n", i, _columnPermutation[ i ] );
    }
    printf( "\n" );
#endif

    SUPERLU_FREE( rhsStorage );
    Destroy_SuperMatrix_Store( &B );
}

//////////////////////////////////////////////////////////////////////
// Solves the system with the given right hand side
//////////////////////////////////////////////////////////////////////
void SuperLU_Solver::solve( const VECTOR &rhs, VECTOR &solution )
{
    // Extra error checking
    TRACE_ASSERT( rhs.size() == _A.nrow, "RHS size mismatch" );

    if ( solution.size() != rhs.size() ) {
        solution.resizeAndWipe( rhs.size() );
    }

    solve( rhs.data(), solution.data() );
}

//////////////////////////////////////////////////////////////////////
// Initialize from a matrix in compressed column format
//////////////////////////////////////////////////////////////////////
void SuperLU_Solver::init( const SPARSE_MATRIX::SparseColumnMatrix &A )
{
    TRACE_ASSERT( A._nrow == A._ncol, "Cannot factor a non-square matrix" );

    // Copy the data fields from the sparse column matrix, so that we
    // have a copy of them
    if ( !(_nzval = doubleMalloc(A._nzmax)) ) {
        TRACE_ASSERT( "Error building nzval array" );
    } else {
        memcpy( (void *)_nzval, (void *)A._x, A._nzmax * sizeof( REAL ) );
    }

    if ( !(_rowindices = intMalloc(A._nzmax)) ) {
        TRACE_ASSERT( "Error building rowindices array" );
    } else {
        memcpy( (void *)_rowindices, (void *)A._i, A._nzmax * sizeof( int ) );
    }

    if ( !(_colptr = intMalloc(A._ncol + 1)) ) {
        TRACE_ASSERT( "Error building colptr array" );
    } else {
        memcpy( (void *)_colptr, (void *)A._p, (A._ncol + 1) * sizeof( int ) );
    }

    // Initialize a SuperLU matrix in compressed column format.  Fortunately,
    // this exactly matches the format of our SparseColumnMatrix
    dCreate_CompCol_Matrix( &_A, A._nrow, A._ncol, A._nzmax,
                            _nzval, // Non-zero values
                            _rowindices, // Non-zero row indices
                            _colptr, // Column pointers
                            SLU_NC, SLU_D, SLU_GE /* Various SuperLU definitions */ );

    // Initialize SuperLU options and stats structs.  Just use defaults for now
    set_default_options( &_options );
    StatInit( &_stat );

    // Initialize permutation vectors
    if ( !(_rowPermutation = intMalloc( A._nrow )) ) {
        TRACE_ASSERT( NULL, "Malloc failed for _rowPermutation");
    }

    if ( !(_columnPermutation = intMalloc( A._ncol )) ) {
        TRACE_ASSERT( NULL, "Malloc failed for _columnPermutation" );
    }

    // Solve the system with a dummy right-hand side.  This will generate the LU factorization
    VECTOR                   rhs( A._nrow, 1.0 );
    VECTOR                   solution;
    REAL                     rhsNorm = rhs.norm2();

    solve( rhs, solution );

    // Check the solution
    SPARSE_MATRIX::matrixMultiply( A, solution.data(), rhs.data(),
                                   // Subtract product from rhs without clearing rhs
                                   -1.0, false );

    printf( "SuperLU_Solver::init: Linear system solution yields relative residual "
            "error of %e (this should be close to zero)\n",
            rhs.norm2() / rhsNorm );

    // Finally, set our options to reflect the fact that the matrix has already been factored
    _options.Fact = FACTORED;

#if 0
    // FIXME: remove this.  Only here for testing purposes
    rhs = VECTOR( A._nrow, 3.0 );
    solution = VECTOR( A._nrow );
    rhsNorm = rhs.norm2();

    solve( rhs, solution );

    SPARSE_MATRIX::matrixMultiply( A, solution.data(), rhs.data(),
                                   // Subtract product from rhs without clearing rhs
                                   -1.0, false );

    printf( "SuperLU LU linear system solution yields relative residual error of %f\n",
            rhs.norm2() / rhsNorm );
#endif
}
