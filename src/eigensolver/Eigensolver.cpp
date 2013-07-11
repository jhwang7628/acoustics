//////////////////////////////////////////////////////////////////////
// Eigensolver.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "Eigensolver.h"

#include <utils/trace.h>

#include <iostream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
Eigensolver::Eigensolver( const SPARSE_MATRIX &K, const SPARSE_MATRIX &M, double shift )
    : _K( K ),
      _linearSolver( NULL )
{
    TRACE_ASSERT( K._nrow == K._ncol && M._nrow == M._ncol && K._nrow == M._nrow,
                  "Invalid matrix dimensions in eigenvalue problem" );

    set_shift( shift );

    set_n( K.rows() );

    M.constructSparseColumnCopy( _M );

    // This might be kind of slow, due to the logarithmic access time for
    // each matrix element
    _inverseMatrix = M;
    _inverseMatrix *= -1.0 * shift;
    _inverseMatrix += K;

    _linearSolver = new SuperLU_Solver( _inverseMatrix );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
Eigensolver::~Eigensolver()
{
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void Eigensolver::times_OP1( double* x, double* opx )
{
    _linearSolver->solve( x, opx );
}

//////////////////////////////////////////////////////////////////////
// Applies the matrix (K - sigma_ * M)^{-1} * M to the input vector
// and places the result in opx.  Also puts M * x in Mx
//////////////////////////////////////////////////////////////////////
void Eigensolver::times_OP2( double* x, double* opx, double* Mx )
{
    //cout << "Calling times_OP2" << endl;

    // Start by multiplying with M
    times_M( x, Mx );

    // Next, apply the inverse of (K - sigma_ * M)
    times_OP1( Mx, opx );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void Eigensolver::times_M( double* x, double* Mx )
{
    SPARSE_MATRIX::matrixMultiply( _M, x, Mx, 1.0 );
}
